function [idx_pos, var_id, t_star, is_degenerate] = lexico_rule(mode, arg1, arg2, arg3, arg4, tol, lu_factors)
% lexico_rule
%  implements the Matrix Lexicographic (Perturbation) anti-cycling pivot rule.
%   [j_in, k_in] = lexico_rule("entering", z_N, N_idx, tol)
%   selects the non-basic column position j_in in {1..n} with the most negative dual
%   reduced cost z_N(j_in) < -tol (identical to Dantzig's largest-coefficient rule).
%
%   [i_out, l_out, t_star, is_degenerate] = lexico_rule("leaving", x_B, d_B, B_idx, B, tol, lu_factors)
%   performs the primal ratio test t* = min { x_B(i) / d_B(i) | d_B(i) > tol } and,
%   whenever a tie occurs (|T| > 1), breaks ties lexicographically by solving the
%   perturbed linear systems B * v_q = e_q column by column (q = 1..m) using cached LU
%   factors of B (or factoring B once via lu(B)) until a unique leaving variable remains.
%
%   hint: Under Vanderbei's symbolic perturbation b(eps) = b + I_m * eps with
%   0 < eps_m << ... << eps_1 << 1 (Chapter 3, Section 4 & Chapter 6), the perturbed
%   basic primal solution is:
%    x_B*(eps) = B^(-1) b + B^(-1) eps = x_B* + sum_{q=1}^m v_q * eps_q
%   where v_q = B \ e_q is the q-th column of B^(-1). For tied rows i in T, the
%   q-th order perturbed ratio is r_{i,q} = v_q(i) / d_B(i). Because rank(B^(-1)) = m,
%   no two rows of B^(-1) can be positively proportional across all columns q = 1..m
%   (Theorem 3.2), guaranteeing that |T| shrinks to 1 and cycling is impossible.

    % 1. Defensive validation of mode argument
    if nargin < 1 || isempty(mode) || ~ischar(mode)
        error("lexico_rule: first argument 'mode' must be 'entering' or 'leaving'");
    endif

    mode_clean = lower(strtrim(mode));

    % 2. Dispatch according to stage ("entering" vs "leaving")
    if strcmp(mode_clean, "entering")
        if nargin < 3 || isempty(arg1) || isempty(arg2)
            error("lexico_rule ('entering'): requires z_N and N_idx arguments");
        endif

        z_N = arg1(:);
        N_idx = arg2(:);
        n = numel(z_N);

        if numel(N_idx) != n
            error(sprintf("lexico_rule ('entering'): dimension mismatch between z_N (%d) and N_idx (%d)", ...
                          n, numel(N_idx)));
        endif

        if nargin < 4 || isempty(arg3)
            tol = 1e-9;
        else
            tol = arg3;
        endif

        % Identify candidate non-basic variables with negative reduced cost (z_N < -tol)
        cand_pos = find(z_N < -tol);
        if isempty(cand_pos)
            idx_pos = 0;
            var_id = 0;
            t_star = 0.0;
            is_degenerate = false;
            return;
        endif

        % Choose candidate with the most negative reduced cost z_N(j)
        [~, best_rel] = min(z_N(cand_pos));
        idx_pos = cand_pos(best_rel);
        var_id = N_idx(idx_pos);
        t_star = 0.0;
        is_degenerate = false;

    elseif strcmp(mode_clean, "leaving")
        if nargin < 5 || isempty(arg1) || isempty(arg2) || isempty(arg3) || isempty(arg4)
            error("lexico_rule ('leaving'): requires x_B, d_B, B_idx, and basis matrix B arguments");
        endif

        x_B = arg1(:);
        d_B = arg2(:);
        B_idx = arg3(:);
        B = arg4;
        m = numel(x_B);

        if numel(d_B) != m || numel(B_idx) != m
            error(sprintf("lexico_rule ('leaving'): dimension mismatch among x_B (%d), d_B (%d), and B_idx (%d)", ...
                          m, numel(d_B), numel(B_idx)));
        endif

        [m_B, n_B] = size(B);
        if m_B != m || n_B != m
            error(sprintf("lexico_rule ('leaving'): basis matrix B must be %d x %d (got %d x %d)", ...
                          m, m, m_B, n_B));
        endif

        if nargin < 6 || isempty(tol)
            tol = 1e-9;
        endif

        % Identify eligible rows with strictly positive primal direction (d_B > tol)
        eligible = find(d_B > tol);
        if isempty(eligible)
            idx_pos = 0;
            var_id = 0;
            t_star = inf;
            is_degenerate = false;
            return;
        endif

        % 0-th order comparison: standard primal ratio test
        ratios = x_B(eligible) ./ d_B(eligible);
        min_ratio = min(ratios);
        t_star = max(0.0, min_ratio);
        is_degenerate = (t_star <= tol);

        tie_mask = abs(ratios - min_ratio) <= tol;
        tie_set = eligible(tie_mask);

        % Lazy column-by-column lexicographic tie-breaking using single LU factorization
        if numel(tie_set) > 1
            if nargin >= 7 && ~isempty(lu_factors) && isstruct(lu_factors)
                lu_f = lu_factors;
            else
                lu_f = struct();
                lu_f.is_sparse = issparse(B);
                if lu_f.is_sparse
                    [lu_f.L, lu_f.U, lu_f.P, lu_f.Q] = lu(B);
                else
                    [lu_f.L, lu_f.U, lu_f.P] = lu(B);
                    lu_f.Q = [];
                endif
            endif

            for q = 1:m
                if lu_f.is_sparse
                    e_q = sparse(q, 1, 1.0, m, 1);
                    v_q = lu_f.Q * (lu_f.U \ (lu_f.L \ (lu_f.P * e_q)));
                else
                    e_q = zeros(m, 1);
                    e_q(q) = 1.0;
                    v_q = lu_f.U \ (lu_f.L \ (lu_f.P * e_q));
                endif

                r_q = v_q(tie_set) ./ d_B(tie_set);
                min_rq = min(r_q);
                tie_set = tie_set(abs(r_q - min_rq) <= tol);

                % Early exit as soon as a unique leaving row is isolated
                if numel(tie_set) == 1
                    break;
                endif
            endfor
        endif

        idx_pos = tie_set(1);
        var_id = B_idx(idx_pos);

    else
        error(sprintf("lexico_rule: unknown mode '%s' (expected 'entering' or 'leaving')", mode));
    endif
endfunction
