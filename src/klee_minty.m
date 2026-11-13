function [A, b, c, x_opt_exact, obj_opt_exact] = klee_minty(n)
% klee_minty
%  generates the n-dimensional Klee-Minty deformed hypercube linear program (1972).
%   [A, b, c, x_opt_exact, obj_opt_exact] = klee_minty(n) constructs the matrices
%   A (n x n lower triangular), b (n x 1), and c (n x 1) for Vanderbei's canonical form:
%    max  sum_{j=1}^n 10^(n-j) x_j
%    s.t. 2 * sum_{j=1}^{i-1} 10^(i-j) x_j + x_i <= 100^(i-1),  i = 1..n
%         x_j >= 0,                                             j = 1..n
%
%   hint: Geometrically, the feasible region is a slightly perturbed hypercube
%   [0, 1] x [0, 100] x ... x [0, 100^(n-1)] with 2^n vertices (Vanderbei, Chapter 4,
%   Section 5, Eq. 4.1). Under Dantzig's largest-coefficient rule (@standard_rule),
%   the Simplex Method is forced to traverse a Hamiltonian path (Gray code) visiting
%   every single vertex of the polytope, executing exactly 2^n - 1 non-degenerate
%   iterations before reaching the unique optimum x* = [0; ...; 0; 100^(n-1)] with
%   objective value zeta* = 100^(n-1).

    % 1. Defensive Programming and Input Validation
    if nargin < 1 || isempty(n)
        error("klee_minty: dimension argument 'n' is required and cannot be empty");
    endif

    if ~isscalar(n) || ~isnumeric(n) || n < 1 || floor(n) != n
        error("klee_minty: dimension 'n' must be a positive integer scalar (n >= 1)");
    endif

    % 2. Vectorized Construction of Technological Matrix A (n x n)
    [I_mat, J_mat] = ndgrid(1:n, 1:n);
    lower_mask = (I_mat > J_mat);
    A = eye(n) + 2.0 * lower_mask .* (10.0 .^ (I_mat - J_mat));

    % 3. Right-Hand Side Vector b and Objective Cost Vector c
    exponents_b = (0:(n - 1))';
    b = 100.0 .^ exponents_b;

    exponents_c = ((n - 1):-1:0)';
    c = 10.0 .^ exponents_c;

    % 4. Exact Analytical Optimal Solution
    x_opt_exact = zeros(n, 1);
    x_opt_exact(n) = 100.0 ^ (n - 1);
    obj_opt_exact = 100.0 ^ (n - 1);
endfunction
