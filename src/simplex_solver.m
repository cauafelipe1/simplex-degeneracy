function [x, obj_val, y, metrics] = simplex_solver(A, b, c, pivot_rule_handle, tol, max_iter)
% simplex_solver
%  solves a linear programming problem in Vanderbei canonical form using the Two-Phase Revised Simplex Method.
%   [x, obj_val, y, metrics] = simplex_solver(A, b, c, pivot_rule_handle, tol, max_iter)
%   maximizes zeta = c^T x subject to A * x <= b and x >= 0, where A is m x n (sparse or dense),
%   b is m x 1, and c is n x 1.
%
%   hint: In the Revised Simplex Method (Vanderbei, Chapter 6), the extended system is
%    [A  I_m] [x; w] = b,   [x; w] >= 0
%   partitioned into basic indices B_idx (|B_idx| = m) and non-basic indices N_idx (|N_idx| = n).
%   Instead of explicitly inverting B or forming B^(-1)*N, each iteration solves three m x m systems
%   via Octave's native backslash operator (\):
%    1. Primal basic values:   B * x_B = b               => x_B = B \ b
%    2. Dual multipliers:      B^T * y = c_B             => y   = B' \ c_B   (BTRAN)
%    3. Dual reduced costs:    z_N = N^T * y - c_N                           (Pricing)
%    4. Primal step direction: B * d_B = N(:, j_in)      => d_B = B \ N(:, j_in) (FTRAN)
%   When min(b) < -tol (infeasible origin), Phase I solves the auxiliary problem max -x_0
%   subject to A * x + I_m * w - 1_m * x_0 = b via a single special initial pivot that swaps x_0
%   with the most negative slack variable w_{i_min} (i_min = argmin b_i), restoring primal
%   feasibility in one step before transitioning to Phase II.

    t_start = tic();

    % 1. Defensive Programming and Input Validation
    if nargin < 3 || isempty(A) || isempty(b) || isempty(c)
        error("simplex_solver: arguments A, b, and c are required and cannot be empty");
    endif

    [m, n] = size(A);
    b = b(:);
    c = c(:);

    if numel(b) != m
        error(sprintf("simplex_solver: dimension mismatch between A (%d rows) and b (%d elements)", ...
                      m, numel(b)));
    endif

    if numel(c) != n
        error(sprintf("simplex_solver: dimension mismatch between A (%d cols) and c (%d elements)", ...
                      n, numel(c)));
    endif

    if nargin < 4 || isempty(pivot_rule_handle)
        pivot_rule_handle = @standard_rule;
    elseif ischar(pivot_rule_handle)
        rule_name = lower(strtrim(pivot_rule_handle));
        if strcmp(rule_name, "standard") || strcmp(rule_name, "standard_rule") || strcmp(rule_name, "dantzig")
            pivot_rule_handle = @standard_rule;
        elseif strcmp(rule_name, "bland") || strcmp(rule_name, "bland_rule")
            pivot_rule_handle = @bland_rule;
        elseif strcmp(rule_name, "lexico") || strcmp(rule_name, "lexico_rule") || strcmp(rule_name, "lexicographic")
            pivot_rule_handle = @lexico_rule;
        else
            error(sprintf("simplex_solver: unknown pivot rule identifier '%s'", pivot_rule_handle));
        endif
    elseif ~isa(pivot_rule_handle, "function_handle")
        error("simplex_solver: pivot_rule_handle must be a function handle or rule name string");
    endif

    if nargin < 5 || isempty(tol)
        tol = 1e-9;
    endif

    if nargin < 6 || isempty(max_iter)
        max_iter = 10000;
    endif

    % 2. Extended Matrix System Assembly [A  I_m]
    if issparse(A)
        I_m = speye(m);
    else
        I_m = eye(m);
    endif

    A_bar = [A, I_m];
    c_bar = [c; zeros(m, 1)];

    N_idx = (1:n)';
    B_idx = ((n + 1):(n + m))';

    phase1_iters = 0;
    phase1_degen = 0;

    % 3. Phase I (Auxiliary Problem Initialization when min(b) < -tol)
    if min(b) < -tol
        art_idx = n + m + 1;
        if issparse(A_bar)
            a_0 = sparse(1:m, ones(m, 1), -ones(m, 1), m, 1);
        else
            a_0 = -ones(m, 1);
        endif

        A_aux = [A_bar, a_0];
        c_aux = [zeros(n + m, 1); -1.0];

        N_aux = [N_idx; art_idx];
        B_aux = B_idx;

        % Special initial pivot: x_0 (at N_aux(n+1)) enters, most negative slack w_{i_min} leaves
        [~, i_min] = min(b);
        leave_var = B_aux(i_min);
        B_aux(i_min) = art_idx;
        N_aux(n + 1) = leave_var;
        phase1_iters = 1;

        % Run Revised Simplex on Auxiliary Problem (max -x_0)
        [B_aux, N_aux, x_B_aux, ~, p1_status, p1_loop_iters, p1_degen] = ...
            run_revised_simplex_phase(A_aux, b, c_aux, B_aux, N_aux, ...
                                      pivot_rule_handle, tol, max_iter - 1, art_idx);

        phase1_iters = phase1_iters + p1_loop_iters;
        phase1_degen = p1_degen;

        if strcmp(p1_status, "MAX_ITER")
            x = nan(n, 1);
            obj_val = nan;
            y = nan(m, 1);
            metrics = build_metrics("MAX_ITER", phase1_iters, phase1_iters, 0, ...
                                    phase1_degen, toc(t_start), nan, B_aux, N_aux);
            return;
        endif

        % Determine value of artificial variable x_0* at Phase I termination
        pos_x0 = find(B_aux == art_idx, 1);
        if ~isempty(pos_x0)
            x0_val = x_B_aux(pos_x0);
        else
            x0_val = 0.0;
        endif

        % Infeasibility Check: if x_0* > tol, the original problem is infeasible
        if x0_val > tol
            x = nan(n, 1);
            obj_val = -inf;
            y = nan(m, 1);
            metrics = build_metrics("INFEASIBLE", phase1_iters, phase1_iters, 0, ...
                                    phase1_degen, toc(t_start), nan, B_aux, N_aux);
            return;
        endif

        % If x_0 is still in the basis at degenerate value x_0* == 0, pivot it out cleanly
        if ~isempty(pos_x0)
            B_mat_aux = A_aux(:, B_aux);
            orig_nonbasic_pos = find(N_aux <= (n + m));
            N_orig_cols = A_aux(:, N_aux(orig_nonbasic_pos));

            if issparse(B_mat_aux)
                e_i0 = sparse(pos_x0, 1, 1.0, m, 1);
            else
                e_i0 = zeros(m, 1);
                e_i0(pos_x0) = 1.0;
            endif

            % Compute row pos_x0 of B^(-1) * N_orig via BTRAN: row_i0 = (B' \ e_i0)' * N_orig
            u_i0 = B_mat_aux' \ e_i0;
            piv_row = (u_i0' * N_orig_cols)';
            [max_piv, best_rel] = max(abs(piv_row));

            if max_piv > tol
                j_clean = orig_nonbasic_pos(best_rel);
                B_aux(pos_x0) = N_aux(j_clean);
                N_aux(j_clean) = art_idx;
            else
                % Redundant constraint row: replace x_0 with a slack variable not currently in B_aux
                unused_slacks = setdiff(((n + 1):(n + m))', B_aux);
                j_clean = find(N_aux == unused_slacks(1), 1);
                B_aux(pos_x0) = N_aux(j_clean);
                N_aux(j_clean) = art_idx;
            endif
        endif

        % Drop artificial variable art_idx from N_aux to restore (n+m)-variable partition
        B_idx = B_aux;
        N_idx = N_aux(N_aux != art_idx);
    endif

    % 4. Phase II (Main Primal Revised Simplex Loop with Original Objective c_bar)
    remaining_iters = max(1, max_iter - phase1_iters);
    [B_idx, N_idx, x_B, y, p2_status, phase2_iters, phase2_degen] = ...
        run_revised_simplex_phase(A_bar, b, c_bar, B_idx, N_idx, ...
                                  pivot_rule_handle, tol, remaining_iters, 0);

    total_iters = phase1_iters + phase2_iters;
    total_degen = phase1_degen + phase2_degen;

    % 5. Post-Processing and Solution Reconstruction
    if strcmp(p2_status, "OPTIMAL")
        x_ext = zeros(n + m, 1);
        x_ext(B_idx) = x_B;
        x = x_ext(1:n);
        obj_val = c' * x;
        duality_gap = abs(b' * y - obj_val);
    elseif strcmp(p2_status, "UNBOUNDED")
        x = nan(n, 1);
        obj_val = inf;
        duality_gap = nan;
    else
        % MAX_ITER reached (e.g., cycling under standard_rule)
        x_ext = zeros(n + m, 1);
        x_ext(B_idx) = x_B;
        x = x_ext(1:n);
        obj_val = c' * x;
        duality_gap = abs(b' * y - obj_val);
    endif

    cpu_elapsed = toc(t_start);
    metrics = build_metrics(p2_status, total_iters, phase1_iters, phase2_iters, ...
                            total_degen, cpu_elapsed, duality_gap, B_idx, N_idx);
endfunction


function [B_idx, N_idx, x_B, y, status, iter_count, degen_count] = ...
    run_revised_simplex_phase(A_bar, b, c_bar, B_idx, N_idx, pivot_rule_handle, tol, max_iter, art_var_idx)
% run_revised_simplex_phase
%  executes the core Revised Simplex iterative loop (Steps 1 to 9 of Vanderbei Chapter 6).

    iter_count = 0;
    degen_count = 0;
    status = "MAX_ITER";

    % Disable singular/near-singular matrix console spam during intense degeneracy benchmarks
    warn_state = warning("off", "Octave:singular-matrix");
    unwind_protect
        while true
            % Extract current basic and non-basic submatrices and cost vectors
            B = A_bar(:, B_idx);
            N = A_bar(:, N_idx);
            c_B = c_bar(B_idx);
            c_N = c_bar(N_idx);

            % Factor basis matrix B once per iteration (sparse: P*B*Q = L*U, dense: P*B = L*U)
            lu_f = factor_basis(B);

            % Step 1: Solve primal basic system B * x_B = b
            x_B = solve_primal_lu(lu_f, b);
            % Clamp tiny floating-point negative residuals on degenerate basic variables
            x_B(x_B < 0 & x_B >= -tol) = 0.0;

            % Early exit for Phase I: if artificial variable x_0 is non-basic or zero, Phase I is optimal
            if art_var_idx > 0
                pos_art = find(B_idx == art_var_idx, 1);
                if isempty(pos_art) || x_B(pos_art) <= tol
                    y = solve_dual_lu(lu_f, c_B);
                    status = "OPTIMAL";
                    break;
                endif
            endif

            % Step 2 (BTRAN): Solve dual system B^T * y = c_B
            y = solve_dual_lu(lu_f, c_B);

            % Step 3 (Pricing): Compute dual reduced costs z_N = N^T * y - c_N
            z_N = N' * y - c_N;

            % Steps 4 & 5: Optimality test and entering variable selection
            [j_in, k_in] = pivot_rule_handle("entering", z_N, N_idx, tol);
            if j_in == 0
                status = "OPTIMAL";
                break;
            endif

            % Check iteration limit before performing pivot
            if iter_count >= max_iter
                status = "MAX_ITER";
                break;
            endif

            % Step 6 (FTRAN): Solve primal step direction B * d_B = N(:, j_in)
            d_B = solve_primal_lu(lu_f, N(:, j_in));

            % Steps 7 & 8: Ratio test and leaving variable selection (reusing lu_f)
            [i_out, l_out, t_star, is_degen] = ...
                pivot_rule_handle("leaving", x_B, d_B, B_idx, B, tol, lu_f);

            if i_out == 0
                status = "UNBOUNDED";
                break;
            endif

            % Phase I stabilization: if x_0 is basic and tied for the minimum ratio, prioritize x_0 leaving
            if art_var_idx > 0
                pos_art = find(B_idx == art_var_idx, 1);
                if ~isempty(pos_art) && d_B(pos_art) > tol
                    r_art = x_B(pos_art) / d_B(pos_art);
                    if abs(r_art - t_star) <= tol
                        i_out = pos_art;
                        l_out = art_var_idx;
                    endif
                endif
            endif

            % Step 9: Update basis partition and telemetry counters
            B_idx(i_out) = k_in;
            N_idx(j_in) = l_out;

            iter_count = iter_count + 1;
            if is_degen
                degen_count = degen_count + 1;
            endif
        endwhile
    unwind_protect_cleanup
        warning(warn_state);
    end_unwind_protect
endfunction


function lu_f = factor_basis(B)
% factor_basis
%  computes a single LU factorization of the basis matrix B for reuse across
%  primal (B*x=b), BTRAN (B'*y=c_B), FTRAN (B*d_B=N_j), and lexicographic tie-breaking.

    lu_f = struct();
    lu_f.is_sparse = issparse(B);
    if lu_f.is_sparse
        [lu_f.L, lu_f.U, lu_f.P, lu_f.Q] = lu(B);
    else
        [lu_f.L, lu_f.U, lu_f.P] = lu(B);
        lu_f.Q = [];
    endif
endfunction


function sol = solve_primal_lu(lu_f, rhs)
% solve_primal_lu
%  solves B * sol = rhs using precomputed LU factors.

    if lu_f.is_sparse
        sol = lu_f.Q * (lu_f.U \ (lu_f.L \ (lu_f.P * rhs)));
    else
        sol = lu_f.U \ (lu_f.L \ (lu_f.P * rhs));
    endif
endfunction


function sol = solve_dual_lu(lu_f, rhs)
% solve_dual_lu
%  solves B^T * sol = rhs using precomputed LU factors.

    if lu_f.is_sparse
        sol = lu_f.P' * (lu_f.L' \ (lu_f.U' \ (lu_f.Q' * rhs)));
    else
        sol = lu_f.P' * (lu_f.L' \ (lu_f.U' \ rhs));
    endif
endfunction


function metrics = build_metrics(status, total_iters, p1_iters, p2_iters, degen_pivots, cpu_time, gap, B_idx, N_idx)
% build_metrics
%  assembles the standardized telemetry struct returned by simplex_solver.

    metrics = struct();
    metrics.status = status;
    metrics.iterations = total_iters;
    metrics.phase1_iterations = p1_iters;
    metrics.phase2_iterations = p2_iters;
    metrics.degenerate_pivots = degen_pivots;
    metrics.cpu_time = cpu_time;
    metrics.duality_gap = gap;
    metrics.B_idx = B_idx;
    metrics.N_idx = N_idx;
endfunction
