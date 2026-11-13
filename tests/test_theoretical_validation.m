function test_theoretical_validation()
% test_theoretical_validation
%  executes the Phase 3 theoretical validation and duality certification suite.
%   test_theoretical_validation() runs automated assertions verifying:
%    1. Chapter 3 (Degeneracy & Cycling): 6-dictionary cycling under @standard_rule
%       versus finite convergence to zeta* = 1 under @bland_rule and @lexico_rule.
%    2. Chapter 4 (Exponential Complexity): exact 2^n - 1 iteration count on the
%       Klee-Minty hypercube for n in {3, 4, 5, 6} under @standard_rule.
%    3. Chapter 5 (Duality Theory): Primal/Dual feasibility, Strong Duality (zero gap),
%       and Complementary Slackness conditions (x_j * z_j = 0, y_i * w_i = 0).
%
%   hint: Theorem 5.2 (Strong Duality) and Theorem 5.3 (Complementary Slackness)
%   provide an independent algebraic certificate of global optimality for any candidate
%   primal-dual pair (x*, y*) without trusting the solver's internal trajectory.

    repo_root = fileparts(fileparts(mfilename("fullpath")));
    addpath(fullfile(repo_root, "src"));
    tol = 1e-9;

    printf("================================================================================\n");
    printf("FASE 3: SUITE DE VALIDACAO TEORICA E CERTIFICACAO DUAL (VANDERBEI CAPS. 3, 4, 5)\n");
    printf("================================================================================\n\n");

    % -------------------------------------------------------------------------
    % BLOCO 1: Contraexemplo Classico de Ciclagem de 6 Dicionarios (Cap. 3)
    % -------------------------------------------------------------------------
    printf("[Bloco 1] Validando Contraexemplo de Ciclagem de 6 Dicionarios (Cap. 3, Sec. 3.3)...\n");
    A_cyc = [0.5, -5.5, -2.5,  9.0; ...
             0.5, -1.5, -0.5,  1.0; ...
             1.0,  0.0,  0.0,  0.0];
    b_cyc = [0.0; 0.0; 1.0];
    c_cyc = [10.0; -57.0; -9.0; -24.0];
    x_cyc_exact = [1.0; 0.0; 1.0; 0.0];
    z_cyc_exact = 1.0;

    [~, z_std, ~, m_std] = simplex_solver(A_cyc, b_cyc, c_cyc, @standard_rule, tol, 36);
    [x_bld, z_bld, y_bld, m_bld] = simplex_solver(A_cyc, b_cyc, c_cyc, @bland_rule, tol, 36);
    [x_lex, z_lex, y_lex, m_lex] = simplex_solver(A_cyc, b_cyc, c_cyc, @lexico_rule, tol, 36);

    assert(strcmp(m_std.status, "MAX_ITER"), "Standard rule deveria ciclar e atingir MAX_ITER");
    assert(m_std.degenerate_pivots == 36, "Todos os pivos no ciclo de Standard devem ser degenerados");
    assert(abs(z_std - 0.0) <= tol, "Objetivo durante o ciclo na origem deve permanecer 0");

    assert(strcmp(m_bld.status, "OPTIMAL"), "Regra de Bland deve convergir para OPTIMAL");
    assert(abs(z_bld - z_cyc_exact) <= tol, "Regra de Bland deve atingir zeta* = 1");
    assert(norm(x_bld - x_cyc_exact, inf) <= tol, "Solucao primal de Bland incorreta");
    verify_optimality_certificate(A_cyc, b_cyc, c_cyc, x_bld, y_bld, z_bld, tol, "Ciclagem (Bland)");

    assert(strcmp(m_lex.status, "OPTIMAL"), "Metodo Lexicografico deve convergir para OPTIMAL");
    assert(abs(z_lex - z_cyc_exact) <= tol, "Metodo Lexicografico deve atingir zeta* = 1");
    assert(norm(x_lex - x_cyc_exact, inf) <= tol, "Solucao primal Lexicografica incorreta");
    verify_optimality_certificate(A_cyc, b_cyc, c_cyc, x_lex, y_lex, z_lex, tol, "Ciclagem (Lexico)");

    printf("  -> Standard: status=%s (%d/%d pivos degenerados, zeta=%.2f) [CICLO CONFIRMADO]\n", ...
           m_std.status, m_std.degenerate_pivots, m_std.iterations, z_std);
    printf("  -> Bland:    status=%s (%d iters, %d degen, zeta*=%.4f, gap=%.2e) [OK]\n", ...
           m_bld.status, m_bld.iterations, m_bld.degenerate_pivots, z_bld, m_bld.duality_gap);
    printf("  -> Lexico:   status=%s (%d iters, %d degen, zeta*=%.4f, gap=%.2e) [OK]\n\n", ...
           m_lex.status, m_lex.iterations, m_lex.degenerate_pivots, z_lex, m_lex.duality_gap);

    % -------------------------------------------------------------------------
    % BLOCO 2: Complexidade Exponencial Exata 2^n - 1 em Klee-Minty (Cap. 4)
    % -------------------------------------------------------------------------
    printf("[Bloco 2] Validando Complexidade Exponencial 2^n - 1 no Cubo de Klee-Minty (Cap. 4)...\n");
    for n_dim = 3:6
        [A_km, b_km, c_km, x_km_ex, z_km_ex] = klee_minty(n_dim);
        expected_iters = (2 ^ n_dim) - 1;

        [x_km, z_km, y_km, m_km] = simplex_solver(A_km, b_km, c_km, @standard_rule, tol, 1000);

        assert(strcmp(m_km.status, "OPTIMAL"), sprintf("Klee-Minty n=%d deve atingir OPTIMAL", n_dim));
        assert(m_km.iterations == expected_iters, ...
               sprintf("Klee-Minty n=%d executou %d iters (esperado 2^%d - 1 = %d)", ...
                       n_dim, m_km.iterations, n_dim, expected_iters));
        assert(m_km.degenerate_pivots == 0, "Klee-Minty nao possui pivos degenerados");
        assert(abs(z_km - z_km_ex) / max(1.0, abs(z_km_ex)) <= tol, "Valor otimo Klee-Minty incorreto");
        assert(norm(x_km - x_km_ex, inf) / max(1.0, norm(x_km_ex, inf)) <= tol, "Vetor x* Klee-Minty incorreto");

        verify_optimality_certificate(A_km, b_km, c_km, x_km, y_km, z_km, tol, ...
                                      sprintf("Klee-Minty n=%d", n_dim));

        printf("  -> n=%d: iters=%2d (2^%d - 1 = %2d), degen=%d, zeta*=%.1e, gap=%.2e [OK]\n", ...
               n_dim, m_km.iterations, n_dim, expected_iters, m_km.degenerate_pivots, z_km, m_km.duality_gap);
    endfor
    printf("\n");

    % -------------------------------------------------------------------------
    % BLOCO 3: Certificacao Dual em Problemas de Duas Fases (Cap. 2 + Cap. 5)
    % -------------------------------------------------------------------------
    printf("[Bloco 3] Validando Certificados Duais e Folgas Complementares em Duas Fases (Cap. 5)...\n");
    A_p2 = [-1.0,  1.0; ...
            -1.0, -2.0; ...
             0.0,  1.0];
    b_p2 = [-1.0; -2.0; 1.0];
    c_p2 = [-2.0; -1.0];

    [x_p2, z_p2, y_p2, m_p2] = simplex_solver(A_p2, b_p2, c_p2, @lexico_rule, tol, 100);
    assert(strcmp(m_p2.status, "OPTIMAL"), "Problema de 2 fases deve ser OPTIMAL");
    verify_optimality_certificate(A_p2, b_p2, c_p2, x_p2, y_p2, z_p2, tol, "Cap. 2 Duas Fases");
    printf("  -> Duas Fases (Cap. 2 Sec. 4.6): zeta*=%.4f, p1=%d, p2=%d, folgas complementares [OK]\n\n", ...
           z_p2, m_p2.phase1_iterations, m_p2.phase2_iterations);

    printf("================================================================================\n");
    printf("SUCESSO: TODAS AS ASSERCOES TEORICAS E CERTIFICADOS DUAIS FORAM VALIDADOS!\n");
    printf("================================================================================\n");
endfunction


function verify_optimality_certificate(A, b, c, x, y, obj_val, tol, label_str)
% verify_optimality_certificate
%  verifies Primal Feasibility, Dual Feasibility, Strong Duality, and Complementary Slackness.

    scale_factor = max([1.0, norm(b, inf), norm(c, inf), abs(obj_val)]);
    scaled_tol = tol * scale_factor;

    % 1. Primal Feasibility: x >= 0 and w = b - A*x >= 0
    w = b - A * x;
    assert(all(x >= -scaled_tol), ...
           sprintf("[%s] Violacao de nao-negatividade primal: min(x) = %.3e", label_str, min(x)));
    assert(all(w >= -scaled_tol), ...
           sprintf("[%s] Violacao de restricao primal (w = b - Ax): min(w) = %.3e", label_str, min(w)));

    % 2. Dual Feasibility: y >= 0 and z = A^T*y - c >= 0
    z = A' * y - c;
    assert(all(y >= -scaled_tol), ...
           sprintf("[%s] Violacao de nao-negatividade dual: min(y) = %.3e", label_str, min(y)));
    assert(all(z >= -scaled_tol), ...
           sprintf("[%s] Violacao de restricao dual (z = A'y - c): min(z) = %.3e", label_str, min(z)));

    % 3. Strong Duality (Theorem 5.2): c^T x* == b^T y*
    primal_obj = c' * x;
    dual_obj = b' * y;
    duality_gap = abs(dual_obj - primal_obj);
    assert(duality_gap <= scaled_tol, ...
           sprintf("[%s] Gap de Dualidade Forte excessivo: %.3e", label_str, duality_gap));

    % 4. Complementary Slackness (Theorem 5.3): x_j * z_j == 0 and y_i * w_i == 0
    comp_primal = max(abs(x .* z));
    comp_dual = max(abs(y .* w));
    assert(comp_primal <= scaled_tol, ...
           sprintf("[%s] Violacao de folga complementar (x .* z): %.3e", label_str, comp_primal));
    assert(comp_dual <= scaled_tol, ...
           sprintf("[%s] Violacao de folga complementar (y .* w): %.3e", label_str, comp_dual));
endfunction
