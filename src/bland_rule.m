function [idx_pos, var_id, t_star, is_degenerate] = bland_rule(mode, arg1, arg2, arg3, arg4, tol, lu_factors)
% bland_rule
%  implements Bland's smallest-subscript anti-cycling pivot rule (1977).
%   [j_in, k_in] = bland_rule("entering", z_N, N_idx, tol)
%   among all non-basic candidates with negative dual reduced cost (z_N(j) < -tol),
%   selects the position j_in whose global variable index k_in = N_idx(j_in) is minimal:
%    k_in = min { N_idx(j) | z_N(j) < -tol }.
%
%   [i_out, l_out, t_star, is_degenerate] = bland_rule("leaving", x_B, d_B, B_idx, B, tol, lu_factors)
%   performs the primal ratio test t* = min { x_B(i) / d_B(i) | d_B(i) > tol } and,
%   among all basic rows tied in T = { i | |r_i - t*| <= tol }, selects the row i_out
%   whose global variable index l_out = B_idx(i_out) is minimal:
%    l_out = min { B_idx(i) | i in T }.
%   The optional lu_factors argument is accepted for polymorphism with lexico_rule.
%
%   hint: Theorem 3.3 (Vanderbei, Chapter 3, Section 5) proves via a contradiction
%   argument on fickle variables that the Simplex Method always terminates in a finite
%   number of steps when both entering and leaving variables follow Bland's smallest-index
%   rule, completely preventing cycling without any linear system overhead.

    % 1. Defensive validation of mode argument
    if nargin < 1 || isempty(mode) || ~ischar(mode)
        error("bland_rule: first argument 'mode' must be 'entering' or 'leaving'");
    endif

    mode_clean = lower(strtrim(mode));

    % 2. Dispatch according to stage ("entering" vs "leaving")
    if strcmp(mode_clean, "entering")
        if nargin < 3 || isempty(arg1) || isempty(arg2)
            error("bland_rule ('entering'): requires z_N and N_idx arguments");
        endif

        z_N = arg1(:);
        N_idx = arg2(:);
        n = numel(z_N);

        if numel(N_idx) != n
            error(sprintf("bland_rule ('entering'): dimension mismatch between z_N (%d) and N_idx (%d)", ...
                          n, numel(N_idx)));
        endif

        if nargin < 4 || isempty(arg3)
            tol = 1e-9;
        else
            tol = arg3;
        endif

        % Identify all candidate positions with strictly negative reduced cost
        cand_pos = find(z_N < -tol);
        if isempty(cand_pos)
            idx_pos = 0;
            var_id = 0;
            t_star = 0.0;
            is_degenerate = false;
            return;
        endif

        % Select the candidate with the smallest global variable index N_idx(j)
        [var_id, min_rel] = min(N_idx(cand_pos));
        idx_pos = cand_pos(min_rel);
        t_star = 0.0;
        is_degenerate = false;

    elseif strcmp(mode_clean, "leaving")
        if nargin < 4 || isempty(arg1) || isempty(arg2) || isempty(arg3)
            error("bland_rule ('leaving'): requires x_B, d_B, and B_idx arguments");
        endif

        x_B = arg1(:);
        d_B = arg2(:);
        B_idx = arg3(:);
        m = numel(x_B);

        if numel(d_B) != m || numel(B_idx) != m
            error(sprintf("bland_rule ('leaving'): dimension mismatch among x_B (%d), d_B (%d), and B_idx (%d)", ...
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

        % Break ties by choosing the basic variable with the smallest global index B_idx(i)
        [var_id, min_rel] = min(B_idx(tie_set));
        idx_pos = tie_set(min_rel);

    else
        error(sprintf("bland_rule: unknown mode '%s' (expected 'entering' or 'leaving')", mode));
    endif
endfunction
