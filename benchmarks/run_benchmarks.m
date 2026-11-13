function results = run_benchmarks(mode, max_iter_netlib, tol)
% run_benchmarks
%  executes the empirical degeneracy benchmark suite across synthetic and NETLIB instances.
%   results = run_benchmarks(mode, max_iter_netlib, tol) evaluates @standard_rule,
%   @bland_rule, and @lexico_rule across benchmark problems, prints a formatted
%   comparative table, and exports results to data/processed/.
%
%   Input arguments (all optional):
%    mode            - "fast" (default): Cycling_6Dict, KleeMinty_n5, NETLIB_degen2 (~32s)
%                      "full": also includes large-scale NETLIB_degen3 (2220 x 1818)
%    max_iter_netlib - iteration limit for NETLIB instances (default: 15000 in fast, 25000 in full)
%    tol             - numerical tolerance (default: 1e-9)
%
%   hint: This benchmark empirically contrasts combinatorial iteration efficiency against
%   per-iteration linear algebra overhead under severe primal degeneracy:
%    - @standard_rule: fast per pivot, but susceptible to cycling/stalling on degenerate facets.
%    - @bland_rule:    zero tie-breaking overhead, guaranteed finite termination (Theorem 3.3),
%                      but sacrifices steep objective ascent by picking the minimum index.
%    - @lexico_rule:   preserves steepest-coefficient ascent and escapes degenerate plateaus
%                      rapidly (Theorem 3.2), trading fewer iterations for extra LU solves.

    repo_root = fileparts(fileparts(mfilename("fullpath")));
    addpath(fullfile(repo_root, "src"));

    % 1. Defensive Programming and Default Arguments
    if nargin < 1 || isempty(mode)
        mode = "fast";
    endif
    mode_clean = lower(strtrim(mode));
    if ~strcmp(mode_clean, "fast") && ~strcmp(mode_clean, "full")
        error(sprintf("run_benchmarks: unknown mode '%s' (expected 'fast' or 'full')", mode));
    endif

    if nargin < 2 || isempty(max_iter_netlib)
        if strcmp(mode_clean, "full")
            max_iter_netlib = 25000;
        else
            max_iter_netlib = 15000;
        endif
    endif

    if nargin < 3 || isempty(tol)
        tol = 1e-9;
    endif

    data_mps_dir = fullfile(repo_root, "data", "mps");
    data_proc_dir = fullfile(repo_root, "data", "processed");
    if exist(data_proc_dir, "dir") != 7
        mkdir(data_proc_dir);
    endif

    % 2. Assemble Benchmark Datasets
    datasets = struct([]);
    ds_count = 0;

    % Dataset 1: Classical 6-Dictionary Cycling Counterexample (Cap. 3, Sec. 3.3)
    ds_count = ds_count + 1;
    datasets(ds_count).name = "Cycling_6Dict";
    datasets(ds_count).A = [0.5, -5.5, -2.5,  9.0; ...
                            0.5, -1.5, -0.5,  1.0; ...
                            1.0,  0.0,  0.0,  0.0];
    datasets(ds_count).b = [0.0; 0.0; 1.0];
    datasets(ds_count).c = [10.0; -57.0; -9.0; -24.0];
    datasets(ds_count).max_iter = 100;

    % Dataset 2: Klee-Minty Deformed Hypercube n = 5 (Cap. 4, Sec. 5)
    [A_km5, b_km5, c_km5] = klee_minty(5);
    ds_count = ds_count + 1;
    datasets(ds_count).name = "KleeMinty_n5";
    datasets(ds_count).A = A_km5;
    datasets(ds_count).b = b_km5;
    datasets(ds_count).c = c_km5;
    datasets(ds_count).max_iter = 500;

    % Dataset 3: NETLIB DEGEN2 (auto-parse .mps if .mat is missing)
    [has_d2, d2_data] = load_netlib_dataset("degen2", data_mps_dir, data_proc_dir);
    if has_d2
        ds_count = ds_count + 1;
        datasets(ds_count).name = "NETLIB_degen2";
        datasets(ds_count).A = d2_data.A;
        datasets(ds_count).b = d2_data.b;
        datasets(ds_count).c = d2_data.c;
        datasets(ds_count).max_iter = max_iter_netlib;
    endif

    % Dataset 4 (full mode only): NETLIB DEGEN3
    if strcmp(mode_clean, "full")
        [has_d3, d3_data] = load_netlib_dataset("degen3", data_mps_dir, data_proc_dir);
        if has_d3
            ds_count = ds_count + 1;
            datasets(ds_count).name = "NETLIB_degen3";
            datasets(ds_count).A = d3_data.A;
            datasets(ds_count).b = d3_data.b;
            datasets(ds_count).c = d3_data.c;
            datasets(ds_count).max_iter = max_iter_netlib;
        endif
    endif

    % 3. Define Pivot Rules to Compare
    rules = {
        "Standard", @standard_rule;
        "Bland",    @bland_rule;
        "Lexico",   @lexico_rule
    };

    % 4. Print Benchmark Header
    sep_line = repmat("=", 1, 122);
    sub_line = repmat("-", 1, 122);
    printf("\n%s\n", sep_line);
    printf("REVISED SIMPLEX BENCHMARK [%s mode] - DEGENERACY & ANTI-CYCLING PIVOT RULES (GNU OCTAVE)\n", ...
           upper(mode_clean));
    printf("%s\n", sep_line);
    printf("%-15s %-13s %-10s %-10s %16s %18s %11s %11s %9s\n", ...
           "Dataset", "Dim (m x n)", "Rule", "Status", "Obj (zeta*)", ...
           "Iters (F1+F2)", "Degen (%)", "Dual Gap", "CPU (s)");
    printf("%s\n", sub_line);

    results = struct([]);
    rec_idx = 0;

    % 5. Execute Benchmark Matrix (Datasets x Rules)
    for d_idx = 1:numel(datasets)
        ds = datasets(d_idx);
        [m_ds, n_ds] = size(ds.A);
        dim_str = sprintf("%d x %d", m_ds, n_ds);

        for r_idx = 1:size(rules, 1)
            rule_label = rules{r_idx, 1};
            rule_handle = rules{r_idx, 2};

            [~, obj_val, ~, met] = simplex_solver(ds.A, ds.b, ds.c, rule_handle, tol, ds.max_iter);

            degen_pct = 100.0 * met.degenerate_pivots / max(1, met.iterations);
            iter_str = sprintf("%d (%d+%d)", met.iterations, met.phase1_iterations, met.phase2_iterations);

            if isnan(obj_val)
                obj_str = "N/A (Phase I)";
            else
                obj_str = sprintf("%.4f", obj_val);
            endif

            if isnan(met.duality_gap)
                gap_str = "-";
            else
                gap_str = sprintf("%.2e", met.duality_gap);
            endif

            printf("%-15s %-13s %-10s %-10s %16s %18s %9.1f%% %11s %9.3f\n", ...
                   ds.name, dim_str, rule_label, met.status, obj_str, ...
                   iter_str, degen_pct, gap_str, met.cpu_time);
            fflush(stdout);

            rec_idx = rec_idx + 1;
            results(rec_idx).dataset = ds.name;
            results(rec_idx).m = m_ds;
            results(rec_idx).n = n_ds;
            results(rec_idx).rule = rule_label;
            results(rec_idx).status = met.status;
            results(rec_idx).obj_val = obj_val;
            results(rec_idx).total_iters = met.iterations;
            results(rec_idx).phase1_iters = met.phase1_iterations;
            results(rec_idx).phase2_iters = met.phase2_iterations;
            results(rec_idx).degen_pivots = met.degenerate_pivots;
            results(rec_idx).degen_pct = degen_pct;
            results(rec_idx).duality_gap = met.duality_gap;
            results(rec_idx).cpu_time = met.cpu_time;
        endfor
        printf("%s\n", sub_line);
    endfor

    % 6. Export Consolidated Results (.mat and .csv) in data/processed/
    out_mat = fullfile(data_proc_dir, "benchmark_results.mat");
    out_csv = fullfile(data_proc_dir, "benchmark_summary.csv");

    save("-mat7-binary", out_mat, "results");

    fid_csv = fopen(out_csv, "w");
    if fid_csv != -1
        fprintf(fid_csv, "dataset,m,n,rule,status,obj_val,total_iters,phase1_iters,phase2_iters,degen_pivots,degen_pct,duality_gap,cpu_time_sec\n");
        for k = 1:numel(results)
            r = results(k);
            fprintf(fid_csv, "%s,%d,%d,%s,%s,%.6f,%d,%d,%d,%d,%.2f,%.4e,%.4f\n", ...
                    r.dataset, r.m, r.n, r.rule, r.status, r.obj_val, ...
                    r.total_iters, r.phase1_iters, r.phase2_iters, ...
                    r.degen_pivots, r.degen_pct, r.duality_gap, r.cpu_time);
        endfor
        fclose(fid_csv);
    endif

    printf("Results exported to:\n");
    printf("  -> %s\n", out_mat);
    printf("  -> %s\n", out_csv);
    printf("%s\n", sep_line);
endfunction


function [found, data] = load_netlib_dataset(name, data_mps_dir, data_proc_dir)
% load_netlib_dataset
%  loads a preprocessed .mat NETLIB instance or converts it from .mps if available.

    mat_path = fullfile(data_proc_dir, [name, ".mat"]);
    mps_path = fullfile(data_mps_dir, [name, ".mps"]);

    if exist(mat_path, "file") == 2
        data = load(mat_path);
        found = true;
        return;
    endif

    if exist(mps_path, "file") == 2
        printf("Preprocessing %s.mps into sparse .mat format...\n", name);
        [A, b, c, raw_mps] = mps_parser(mps_path, mat_path, true);
        data = struct("A", A, "b", b, "c", c, "raw_mps", raw_mps);
        found = true;
        return;
    endif

    printf("Warning: NETLIB instance '%s' not found. Run 'bash scripts/download_netlib.sh' to include it.\n", name);
    data = struct();
    found = false;
endfunction
