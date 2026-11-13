OCTAVE ?= octave --no-gui

.PHONY: all test bench bench-full data clean help

all:
	$(OCTAVE) main.m

test:
	$(OCTAVE) main.m --test

bench:
	$(OCTAVE) main.m --bench

bench-full:
	$(OCTAVE) main.m --full

data:
	bash scripts/download_netlib.sh

clean:
	rm -f octave-workspace data/processed/*.mat data/processed/*.csv

help:
	@echo "Available targets:"
	@echo "  make all         - Run theoretical validation + fast benchmarks"
	@echo "  make test        - Run theoretical validation suite only (< 1s)"
	@echo "  make bench       - Run fast benchmarks (Cycling, Klee-Minty, NETLIB degen2)"
	@echo "  make bench-full  - Run full benchmarks including NETLIB degen3"
	@echo "  make data        - Download and preprocess NETLIB degen2 & degen3 instances"
	@echo "  make clean       - Remove generated .mat/.csv files and octave-workspace"
