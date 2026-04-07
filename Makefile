SHELL := /bin/bash

CXX ?= g++
PYTHON ?= python3
BUILD_DIR ?= build
PARAMS_FILE ?= params.yaml
PREDICTOR_MK := $(BUILD_DIR)/predictor.mk

COMMON_FLAGS := -std=c++20 -O3
CBP_WARN_FLAGS := -Wall -Wextra -pedantic -Wold-style-cast -Werror -Wno-deprecated-declarations -Wno-mismatched-tags
# Keep reference build less strict because upstream reference code is warning-heavy.
REFERENCE_WARN_FLAGS := -Wall -Wextra -pedantic -Wno-deprecated-declarations -Wno-mismatched-tags
# Profile build is less strict since it's instrumentation code
PROFILE_WARN_FLAGS := -Wall -Wextra -pedantic -Wno-deprecated-declarations -Wno-mismatched-tags
EXTRA_COMMON_FLAGS ?=
EXTRA_CBP_FLAGS ?=

# for trace analysis
TRACE_DIR := traces
TRACE_READER := trace_files/trace_reader.hpp

-include $(PREDICTOR_MK)
PREDICTOR_TYPE ?= tage<>
TRACE ?= ./gcc_test_trace.gz
TRACE_NAME ?= test
WARMUP ?= 1000
MEASURE ?= 40000
REGION_SHIFT ?= 12


.PHONY: all help cbp reference predictor-config run-cbp run-reference trace-analyze run-trace-analyze ahead-block-analyze run-ahead-block-analyze cbp-profile-acc cbp-profile-acc-regions cbp-profile-analyze cbp-profile-analyze-regions cbp-monitor monitor-vis compare compare-all quick-eval quick-eval-all test-tage-compile gsweep gsweep-report sweep sweep-report save list-saved clean clean-cbp clean-monitor clean-profile clean-reference clean-out

all: cbp reference

help:
	@echo "Targets:"
	@echo "  make cbp                    Build CBP simulator (uses $(PARAMS_FILE) if present)"
	@echo "  make reference              Build CBP2025 reference predictor"
	@echo "  make run-cbp                Run cbp on TRACE with direct args"
	@echo "  make run-reference          Run reference on TRACE with direct args"
	@echo "  make cbp-profile                Build profiling harness"
	@echo "  make cbp-profile-acc            Run profiler in accuracy mode (fast iteration)"
	@echo "  make cbp-profile-acc-regions    Accuracy mode with per-region breakdown"
	@echo "  make cbp-profile-analyze        Run profiler in analyze mode (per-function breakdown)"
	@echo "  make cbp-profile-analyze-regions Analyze mode with per-region breakdown"
	@echo "  make cbp-monitor                Build+run with TAGE monitor (output to MONITOR_OUT)"
	@echo "  make quick-eval                 Run current predictor on 20 representative traces"
	@echo "  make quick-eval-all             Compare PRED_A vs PRED_B on representative traces"
	@echo "  make trace-analyze              Build trace analysis tool"
	@echo "  make run-trace-analyze          Run trace analyzer on TRACE"
	@echo "  make block-analyze              Build block analyzer for ahead banking"
	@echo "  make run-block-analyze          Run block analyzer (BLOCK_FW, BLOCK_BANKS, BLOCK_INSTR)"
	@echo "  make ahead-block-analyze        Build ahead block analyzer (N-capped blocks, secondary tag)"
	@echo "  make run-ahead-block-analyze    Run ahead analyzer (AHEAD_LINEINST, AHEAD_N, AHEAD_BANKS, AHEAD_INSTR)"
	@echo "  make predictor-config           Generate $(PREDICTOR_MK) from $(PARAMS_FILE)"
	@echo "  make save                   Bookmark current cbp binary (SAVE_NAME=label)"
	@echo "  make list-saved             List bookmarked binaries"
	@echo "  make clean                  Remove generated build artifacts (preserves saved)"
	@echo
	@echo "Variables you can override:"
	@echo "  TRACE=...            Trace file (default: $(TRACE))"
	@echo "  REGION_SHIFT=...     PC region bit shift for analysis (default: $(REGION_SHIFT))"
	@echo "                       12=4KB pages, 21=2MB, 30=1GB"
	@echo "  TRACE_NAME=... WARMUP=... MEASURE=..."
	@echo "  (all default from $(PARAMS_FILE); CLI values still override)"
	@echo "  EXTRA_COMMON_FLAGS=... EXTRA_CBP_FLAGS=..."
	@echo "  PARAMS_FILE=... PREDICTOR_TYPE=... CXX=... PYTHON=..."

$(BUILD_DIR):
	mkdir -p $@

out:
	mkdir -p $@

predictor-config: $(PREDICTOR_MK)

$(PREDICTOR_MK): scripts/gen_predictor_config.py $(PARAMS_FILE) | $(BUILD_DIR)
	$(PYTHON) scripts/gen_predictor_config.py --input $(PARAMS_FILE) --output $@

$(BUILD_DIR)/cbp: cbp.cpp cbp.hpp branch_predictor.hpp $(TRACE_READER) harcom.hpp $(wildcard predictors/*.hpp) $(PREDICTOR_MK) | $(BUILD_DIR)
	$(CXX) $(COMMON_FLAGS) $(EXTRA_COMMON_FLAGS) $(CBP_WARN_FLAGS) $(EXTRA_CBP_FLAGS) -Itrace_files -o $@ cbp.cpp -lz -DPREDICTOR='$(PREDICTOR_TYPE)'

cbp: $(BUILD_DIR)/cbp

$(BUILD_DIR)/reference: reference.cpp $(TRACE_READER) seznec_cbp2025.h | $(BUILD_DIR)
	$(CXX) $(COMMON_FLAGS) $(EXTRA_COMMON_FLAGS) $(REFERENCE_WARN_FLAGS) -Itrace_files -o $@ reference.cpp -lz

reference: $(BUILD_DIR)/reference

$(BUILD_DIR)/trace-analyze: trace_files/trace_analyzer.cpp trace_files/trace_reader.hpp | $(BUILD_DIR)
	$(CXX) $(COMMON_FLAGS) $(EXTRA_COMMON_FLAGS) $(CBP_WARN_FLAGS) -Itrace_files -o $@ trace_files/trace_analyzer.cpp -lz

trace-analyze: $(BUILD_DIR)/trace-analyze

$(BUILD_DIR)/block-analyze: trace_files/block_analyzer.cpp trace_files/trace_reader.hpp | $(BUILD_DIR)
	$(CXX) $(COMMON_FLAGS) $(EXTRA_COMMON_FLAGS) -Wall -Wextra -Itrace_files -o $@ trace_files/block_analyzer.cpp -lz

block-analyze: $(BUILD_DIR)/block-analyze

BLOCK_FW ?= 16
BLOCK_BANKS ?= 8
BLOCK_INSTR ?= 0
run-block-analyze: $(BUILD_DIR)/block-analyze | out
	./$(BUILD_DIR)/block-analyze $(TRACE) $(BLOCK_FW) $(BLOCK_BANKS) $(BLOCK_INSTR) | tee out/block_analysis.txt

$(BUILD_DIR)/ahead-block-analyze: trace_files/ahead_block_analyzer.cpp trace_files/trace_reader.hpp | $(BUILD_DIR)
	$(CXX) $(COMMON_FLAGS) $(EXTRA_COMMON_FLAGS) -Wall -Wextra -Itrace_files -o $@ trace_files/ahead_block_analyzer.cpp -lz

ahead-block-analyze: $(BUILD_DIR)/ahead-block-analyze

AHEAD_LINEINST ?= 64
AHEAD_N ?= 7
AHEAD_BANKS ?= 8
AHEAD_INSTR ?= 0
run-ahead-block-analyze: $(BUILD_DIR)/ahead-block-analyze | out
	./$(BUILD_DIR)/ahead-block-analyze $(TRACE) $(AHEAD_LINEINST) $(AHEAD_N) $(AHEAD_BANKS) $(AHEAD_INSTR) | tee out/ahead_block_analysis.txt

run-cbp: $(BUILD_DIR)/cbp
	./$(BUILD_DIR)/cbp $(TRACE) $(TRACE_NAME) $(WARMUP) $(MEASURE)

run-reference: $(BUILD_DIR)/reference
	./$(BUILD_DIR)/reference $(TRACE) $(TRACE_NAME) $(WARMUP) $(MEASURE)

run-trace-analyze: $(BUILD_DIR)/trace-analyze | out
	./$(BUILD_DIR)/trace-analyze $(TRACE) $(REGION_SHIFT) | tee out/trace_analysis.txt

$(BUILD_DIR)/cbp-profile: cbp_profile.cpp cbp.hpp branch_predictor.hpp $(TRACE_READER) harcom.hpp $(wildcard predictors/*.hpp) $(PREDICTOR_MK) | $(BUILD_DIR)
	$(CXX) $(COMMON_FLAGS) $(EXTRA_COMMON_FLAGS) $(PROFILE_WARN_FLAGS) -Itrace_files -o $@ cbp_profile.cpp -lz -DPREDICTOR='$(PREDICTOR_TYPE)' -DENABLE_REGION_PROFILING

cbp-profile-acc: $(BUILD_DIR)/cbp-profile | out
	./$(BUILD_DIR)/cbp-profile --format csv --mode acc --no-score $(TRACE) $(TRACE_NAME) $(WARMUP) $(MEASURE) 1> out/profile.csv 2> out/profile_acc.txt
	@echo "=== Accuracy Analysis ===" && tail -20 out/profile_acc.txt

cbp-profile-acc-regions: $(BUILD_DIR)/cbp-profile | out
	./$(BUILD_DIR)/cbp-profile --format csv --mode acc --no-score --regions $(TRACE) $(TRACE_NAME) $(WARMUP) $(MEASURE) 1> out/profile.csv 2> out/profile_acc.txt
	@echo "=== Accuracy Analysis with Regions ===" && tail -40 out/profile_acc.txt

cbp-profile-analyze: $(BUILD_DIR)/cbp-profile | out
	./$(BUILD_DIR)/cbp-profile --format csv --mode analyze --profile $(TRACE) $(TRACE_NAME) $(WARMUP) $(MEASURE) 1> out/profile.csv 2> out/profile_analyze.txt
	@echo "=== Per-Function Analysis ===" && tail -30 out/profile_analyze.txt

cbp-profile-analyze-regions: $(BUILD_DIR)/cbp-profile | out
	./$(BUILD_DIR)/cbp-profile --format csv --mode analyze --profile --regions $(TRACE) $(TRACE_NAME) $(WARMUP) $(MEASURE) 1> out/profile.csv 2> out/profile_analyze.txt
	@echo "=== Per-Function Analysis with Regions ===" && tail -40 out/profile_analyze.txt

# TAGE Monitor: build with SW instrumentation and run
# Requires CHEATING_MODE for val<> data extraction. Zero cost when not enabled.
MONITOR_FLAGS := -DTAGE_MONITOR -DCHEATING_MODE -DREAD_WRITE_RAM
MONITOR_OUT ?= out/monitor.txt

$(BUILD_DIR)/cbp-monitor: cbp.cpp cbp.hpp branch_predictor.hpp $(TRACE_READER) harcom.hpp $(wildcard predictors/*.hpp predictors/custom/*.hpp) $(PREDICTOR_MK) | $(BUILD_DIR)
	$(CXX) $(COMMON_FLAGS) $(EXTRA_COMMON_FLAGS) $(PROFILE_WARN_FLAGS) $(MONITOR_FLAGS) -Itrace_files -o $@ cbp.cpp -lz -DPREDICTOR='$(PREDICTOR_TYPE)'

MONITOR_MEASURE ?= $(MEASURE)
PROFILE_MEASURE ?= 100000

cbp-monitor: $(BUILD_DIR)/cbp-monitor $(BUILD_DIR)/cbp-profile | out
	@echo "=== Running Monitor ($(MONITOR_MEASURE) instr) ==="
	./$(BUILD_DIR)/cbp-monitor $(TRACE) $(TRACE_NAME) $(WARMUP) $(MONITOR_MEASURE) 1> out/monitor_csv.txt 2> $(MONITOR_OUT)
	# @echo "=== Running Profile Analyze ($(PROFILE_MEASURE) instr) ==="
	# ./$(BUILD_DIR)/cbp-profile --format csv --mode analyze --profile $(TRACE) $(TRACE_NAME) $(WARMUP) $(PROFILE_MEASURE) 1> out/profile.csv 2> out/profile_analyze.txt
	@echo ""
	# @echo "=== Profile ===" && tail -15 out/profile_analyze.txt
	@echo ""
	@echo "=== Monitor ===" && cat $(MONITOR_OUT)
	@echo ""
	@echo "Files: $(MONITOR_OUT), out/monitor_csv.txt"

monitor-vis:
	python3 scripts/monitor_vis.py out/monitor_csv.txt

# Compare two predictors side-by-side.
# Single trace:  make compare TRACE=path/to/trace.gz
# All traces:    make compare-all TRACE_DIR=path/to/traces/
PRED_A ?= Tage<>
PRED_B ?= tage<>
TRACE_DIR ?= ./traces
compare:
	scripts/compare_predictors.sh '$(PRED_A)' '$(PRED_B)' '$(TRACE)' $(WARMUP) $(MEASURE) '$(EXTRA_COMMON_FLAGS) $(EXTRA_CBP_FLAGS)'

compare-all:
	scripts/compare_predictors.sh --dir '$(PRED_A)' '$(PRED_B)' '$(TRACE_DIR)' $(WARMUP) $(MEASURE) '$(EXTRA_COMMON_FLAGS) $(EXTRA_CBP_FLAGS)'

$(BUILD_DIR)/test-tage-compile: tests/test_tage_compile.cpp predictors/custom/Tage.hpp predictors/custom/TageTable.hpp harcom.hpp | $(BUILD_DIR)
	$(CXX) $(COMMON_FLAGS) $(EXTRA_COMMON_FLAGS) $(CBP_WARN_FLAGS) -Itrace_files -o $@ $< -lz

test-tage-compile: $(BUILD_DIR)/test-tage-compile

# Quick evaluation on representative trace subset (20 traces)
QUICK_JOBS ?= $(shell nproc)
QUICK_OUT ?= out/quick

quick-eval: $(BUILD_DIR)/cbp
	scripts/quick_eval.sh ./$(BUILD_DIR)/cbp $(TRACE_DIR) $(QUICK_OUT) $(QUICK_JOBS)

# Compare two predictors on the representative subset
quick-eval-all: | $(BUILD_DIR)
	@echo "=== Building and evaluating: $(PRED_A) ==="
	$(CXX) $(COMMON_FLAGS) $(EXTRA_COMMON_FLAGS) $(CBP_WARN_FLAGS) $(EXTRA_CBP_FLAGS) -Itrace_files -o $(BUILD_DIR)/quick_a cbp.cpp -lz -DPREDICTOR='$(PRED_A)'
	scripts/quick_eval.sh ./$(BUILD_DIR)/quick_a $(TRACE_DIR) out/quick_a_results $(QUICK_JOBS)
	@echo ""
	@echo "=== Building and evaluating: $(PRED_B) ==="
	$(CXX) $(COMMON_FLAGS) $(EXTRA_COMMON_FLAGS) $(CBP_WARN_FLAGS) $(EXTRA_CBP_FLAGS) -Itrace_files -o $(BUILD_DIR)/quick_b cbp.cpp -lz -DPREDICTOR='$(PRED_B)'
	scripts/quick_eval.sh ./$(BUILD_DIR)/quick_b $(TRACE_DIR) out/quick_b_results $(QUICK_JOBS)

GSWEEP_CONFIGS ?= 60
GSWEEP_JOBS ?= $(shell nproc)
GSWEEP_SEED ?= 42

gsweep:
	$(PYTHON) scripts/gaussian_sweep.py \
		-n $(GSWEEP_CONFIGS) --seed $(GSWEEP_SEED) -j $(GSWEEP_JOBS) \
		--trace-dir $(TRACE_DIR) \
		--extra-flags '$(EXTRA_COMMON_FLAGS) $(EXTRA_CBP_FLAGS)' --resume

gsweep-report:
	$(PYTHON) scripts/gaussian_sweep.py --report out/gsweep/results.csv --top $(GSWEEP_CONFIGS)

SWEEP_TIER ?= 1
SWEEP_JOBS ?= $(shell nproc)
SWEEP_WARMUP ?= 1000
SWEEP_MEASURE ?= 10000

sweep:
	mkdir -p out/sweep
	$(PYTHON) scripts/tage_sweep.py \
		--tier $(SWEEP_TIER) --jobs $(SWEEP_JOBS) \
		--traces $(wildcard traces/*_trace.gz) \
		--warmup $(SWEEP_WARMUP) --measure $(SWEEP_MEASURE) \
		--build-dir out/sweep/bin --output out/sweep/results.csv \
		--extra-flags '$(EXTRA_COMMON_FLAGS) $(EXTRA_CBP_FLAGS)' --resume

sweep-report:
	$(PYTHON) scripts/tage_sweep.py --report out/sweep/results.csv --top 20

SAVE_DIR := $(BUILD_DIR)/saved
SAVE_NAME ?= $(shell date +%Y%m%d_%H%M%S)

$(SAVE_DIR):
	mkdir -p $@

save: $(BUILD_DIR)/cbp | $(SAVE_DIR)
	cp $(BUILD_DIR)/cbp $(SAVE_DIR)/$(SAVE_NAME)
	@echo "Saved to $(SAVE_DIR)/$(SAVE_NAME)"

list-saved:
	@ls -lt $(SAVE_DIR)/ 2>/dev/null || echo "No saved binaries"

clean:
	rm -f $(BUILD_DIR)/cbp $(BUILD_DIR)/reference $(BUILD_DIR)/trace-analyze
	rm -f $(BUILD_DIR)/block-analyze $(BUILD_DIR)/ahead-block-analyze
	rm -f $(BUILD_DIR)/cbp-profile $(BUILD_DIR)/cbp-monitor
	rm -f $(BUILD_DIR)/test-tage-compile $(BUILD_DIR)/quick_a $(BUILD_DIR)/quick_b
	rm -f $(PREDICTOR_MK)
	rm -rf out/*

# Clean individual targets: make clean-cbp, clean-monitor, clean-profile, clean-out
clean-cbp:
	rm -f $(BUILD_DIR)/cbp $(PREDICTOR_MK)
clean-monitor:
	rm -f $(BUILD_DIR)/cbp-monitor
clean-profile:
	rm -f $(BUILD_DIR)/cbp-profile
clean-reference:
	rm -f $(BUILD_DIR)/reference
clean-out:
	rm -rf out/*
