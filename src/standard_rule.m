function [idx_pos, var_id, t_star, is_degenerate] = standard_rule(mode, arg1, arg2, arg3, arg4, tol, lu_factors)
% standard_rule
%  implements Dantzig's standard (largest-coefficient) pivot rule for the Revised Simplex Method.
%   [j_in, k_in] = standard_rule("entering", z_N, N_idx, tol)
%   selects the non-basic column position j_in in {1..n} with the most negative dual
%   reduced cost z_N(j_in) < -tol (equivalent to the largest positive objective
%   coefficient c_bar_j = -z_N(j) > tol in dictionary notation).
%
%   [i_out, l_out, t_star, is_degenerate] = standard_rule("leaving", x_B, d_B, B_idx, B, tol, lu_factors)
%   performs the primal ratio test t* = min { x_B(i) / d_B(i) | d_B(i) > tol } and
%   selects the first basic position i_out in {1..m} achieving the minimum ratio.
%   The optional lu_factors argument is accepted for polymorphism with lexico_rule.
%
%   hint: Dantzig's rule maximizes the instantaneous rate of objective improvement per
%   unit step of the entering variable. However, at degenerate vertices (where x_B(i) = 0
%   and t* = 0), breaking ties arbitrarily by selecting the first tied position in the
%   basis can lead to infinite cycling, as demonstrated in Vanderbei's 6-dictionary
%   counterexample (Chapter 3, Section 3.3).

    % 1. Defensive validation of mode argument
    if nargin < 1 || isempty(mode) || ~ischar(mode)
        error("standard_rule: first argument 'mode' must be 'entering' or 'leaving'");
    endif

    mode_clean = lower(strtrim(mode));

    % 2. Dispatch according to stage ("entering" vs "leaving")
    if strcmp(mode_clean, "entering")
        if nargin < 3 || isempty(arg1) || isempty(arg2)
            error("standard_rule ('entering'): requires z_N and N_idx arguments");
        endif

        z_N = arg1(:);
        N_idx = arg2(:);
        n = numel(z_N);

        if numel(N_idx) != n
            error(sprintf("standard_rule ('entering'): dimension mismatch between z_N (%d) and N_idx (%d)", ...
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
        if nargin < 4 || isempty(arg1) || isempty(arg2) || isempty(arg3)
            error("standard_rule ('leaving'): requires x_B, d_B, and B_idx arguments");
        endif

        x_B = arg1(:);
        d_B = arg2(:);
        B_idx = arg3(:);
        m = numel(x_B);

        if numel(d_B) != m || numel(B_idx) != m
            error(sprintf("standard_rule ('leaving'): dimension mismatch among x_B (%d), d_B (%d), and B_idx (%d)", ...
                          m, numel(d_B), numel(B_idx)));
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

        % Compute primal ratios and find minimum step size t*
        ratios = x_B(eligible) ./ d_B(eligible);
        min_ratio = min(ratios);
        t_star = max(0.0, min_ratio);
        is_degenerate = (t_star <= tol);

        % Identify tie set T = { i in eligible | |r_i - min_ratio| <= tol }
        tie_mask = abs(ratios - min_ratio) <= tol;
        tie_set = eligible(tie_mask);

        % Standard rule breaks ties by picking the first position in basis order
        idx_pos = tie_set(1);
        var_id = B_idx(idx_pos);

    else
        error(sprintf("standard_rule: unknown mode '%s' (expected 'entering' or 'leaving')", mode));
    endif
endfunction
