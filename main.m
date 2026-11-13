% main.m
%  Main entrypoint for the Simplex Degeneracy validation and benchmark suite.
%
%  Usage from terminal:
%   octave --no-gui main.m           % Run theoretical validation + fast benchmarks
%   octave --no-gui main.m --test    % Run only theoretical validation tests
%   octave --no-gui main.m --bench   % Run only fast benchmark suite
%   octave --no-gui main.m --full    % Run validation + full benchmark (incl. degen3)
%
%  Usage from Octave interactive prompt:
%   run main.m

clear; clc;

repo_root = fileparts(mfilename("fullpath"));
if isempty(repo_root)
    repo_root = pwd();
endif

addpath(fullfile(repo_root, "src"));
addpath(fullfile(repo_root, "tests"));
addpath(fullfile(repo_root, "benchmarks"));

args = argv();
run_test = true;
run_bench = true;
bench_mode = "fast";

if ~isempty(args)
    if any(strcmp(args, "--help")) || any(strcmp(args, "-h"))
        printf("Simplex Degeneracy Benchmark Suite (GNU Octave)\n");
        printf("Usage:\n");
        printf("  octave --no-gui main.m          Run theoretical validation + fast benchmarks\n");
        printf("  octave --no-gui main.m --test   Run theoretical validation tests only\n");
        printf("  octave --no-gui main.m --bench  Run fast benchmark suite only\n");
        printf("  octave --no-gui main.m --full   Run validation + full benchmarks (with degen3)\n");
        return;
    elseif any(strcmp(args, "--test"))
        run_test = true;
        run_bench = false;
    elseif any(strcmp(args, "--bench"))
        run_test = false;
        run_bench = true;
    elseif any(strcmp(args, "--full"))
        run_test = true;
        run_bench = true;
        bench_mode = "full";
    endif
endif

if run_test
    test_theoretical_validation();
endif

if run_bench
    run_benchmarks(bench_mode);
endif
