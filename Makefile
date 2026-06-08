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

# Hash of PREDICTOR_TYPE + flags so parallel runs with different configs get separate binaries.
# Same config reuses its cached binary.
PRED_HASH := $(shell echo '$(PREDICTOR_TYPE)|$(EXTRA_COMMON_FLAGS)|$(EXTRA_CBP_FLAGS)' | md5sum | cut -c1-8)


.PHONY: all help cbp reference predictor-config run-cbp run-reference trace-analyze run-trace-analyze ahead-block-analyze run-ahead-block-analyze cbp-profile-acc cbp-profile-acc-regions cbp-profile-analyze cbp-profile-analyze-regions cbp-monitor monitor-vis eval-monitor eval-monitor-vis compare compare-all group-frag-analyze run-group-frag-analyze quick-eval quick-eval-all test-tage-compile gsweep gsweep-report sweep sweep-report gradient gradient-report debug-print run-debug-print cbp-1c cbp-2c cbp-both quick-eval-1c quick-eval-2c quick-eval-both debug-print-1c debug-print-2c run-debug-print-1c run-debug-print-2c psweep-init psweep-generate psweep-timing psweep-eval psweep-run psweep-resume psweep-report save list-saved clean clean-cbp clean-monitor clean-profile clean-reference clean-out

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
	@echo "  make eval-monitor               Run monitor on 20 representative traces"
	@echo "  make eval-monitor-vis           Generate visualizations from eval-monitor output"
	@echo "  make cbp-1c / cbp-2c            Build 1-cycle / 2-cycle competition configs"
	@echo "  make cbp-both                   Build both competition configs"
	@echo "  make quick-eval-1c / -2c        Evaluate 1-cycle / 2-cycle on representative traces"
	@echo "  make quick-eval-both            Build+evaluate both configs side by side"
	@echo "  make debug-print                Build with DEBUG_PRINT (signal timing dumps)"
	@echo "  make run-debug-print            Build+run debug print (warmup=2, measure=2)"
	@echo "  make run-debug-print-1c / -2c   Debug print for 1-cycle / 2-cycle configs"
	@echo "  make quick-eval                 Run current predictor on 20 representative traces"
	@echo "  make quick-eval-all             Compare PRED_A vs PRED_B on representative traces"
	@echo "  make trace-analyze              Build trace analysis tool"
	@echo "  make run-trace-analyze          Run trace analyzer on TRACE"
	@echo "  make block-analyze              Build block analyzer for ahead banking"
	@echo "  make run-block-analyze          Run block analyzer (BLOCK_FW, BLOCK_BANKS, BLOCK_INSTR)"
	@echo "  make ahead-block-analyze        Build ahead block analyzer (N-capped blocks, secondary tag)"
	@echo "  make run-ahead-block-analyze    Run ahead analyzer (AHEAD_LINEINST, AHEAD_N, AHEAD_BANKS, AHEAD_INSTR)"
	@echo "  make run-group-frag-analyze     Run group fragmentation analyzer (FRAG_LINEINST, FRAG_N, FRAG_INSTR)"
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

$(BUILD_DIR)/cbp-$(PRED_HASH): cbp.cpp cbp.hpp branch_predictor.hpp $(TRACE_READER) harcom.hpp $(wildcard predictors/*.hpp) $(PREDICTOR_MK) | $(BUILD_DIR)
	$(CXX) $(COMMON_FLAGS) $(EXTRA_COMMON_FLAGS) $(CBP_WARN_FLAGS) $(EXTRA_CBP_FLAGS) -Itrace_files -o $@ cbp.cpp -lz -DPREDICTOR='$(PREDICTOR_TYPE)'

cbp: $(BUILD_DIR)/cbp-$(PRED_HASH)

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

$(BUILD_DIR)/group-frag-analyze: trace_files/group_frag_analyzer.cpp trace_files/trace_reader.hpp | $(BUILD_DIR)
	$(CXX) $(COMMON_FLAGS) $(EXTRA_COMMON_FLAGS) -Wall -Wextra -Itrace_files -o $@ trace_files/group_frag_analyzer.cpp -lz

group-frag-analyze: $(BUILD_DIR)/group-frag-analyze

FRAG_LINEINST ?= 256
FRAG_N ?= 7
FRAG_INSTR ?= 0
run-group-frag-analyze: $(BUILD_DIR)/group-frag-analyze | out
	./$(BUILD_DIR)/group-frag-analyze $(TRACE) $(FRAG_LINEINST) $(FRAG_N) $(FRAG_INSTR) | tee out/group_frag_analysis.txt

run-cbp: $(BUILD_DIR)/cbp-$(PRED_HASH)
	./$(BUILD_DIR)/cbp-$(PRED_HASH) $(TRACE) $(TRACE_NAME) $(WARMUP) $(MEASURE)

run-reference: $(BUILD_DIR)/reference
	./$(BUILD_DIR)/reference $(TRACE) $(TRACE_NAME) $(WARMUP) $(MEASURE)

run-trace-analyze: $(BUILD_DIR)/trace-analyze | out
	./$(BUILD_DIR)/trace-analyze $(TRACE) $(REGION_SHIFT) | tee out/trace_analysis.txt

$(BUILD_DIR)/cbp-profile-$(PRED_HASH): cbp_profile.cpp cbp.hpp branch_predictor.hpp $(TRACE_READER) harcom.hpp $(wildcard predictors/*.hpp) $(PREDICTOR_MK) | $(BUILD_DIR)
	$(CXX) $(COMMON_FLAGS) $(EXTRA_COMMON_FLAGS) $(PROFILE_WARN_FLAGS) -Itrace_files -o $@ cbp_profile.cpp -lz -DPREDICTOR='$(PREDICTOR_TYPE)' -DENABLE_REGION_PROFILING

cbp-profile-acc: $(BUILD_DIR)/cbp-profile-$(PRED_HASH) | out
	./$(BUILD_DIR)/cbp-profile-$(PRED_HASH) --format csv --mode acc --no-score $(TRACE) $(TRACE_NAME) $(WARMUP) $(MEASURE) 1> out/profile.csv 2> out/profile_acc.txt
	@echo "=== Accuracy Analysis ===" && tail -20 out/profile_acc.txt

cbp-profile-acc-regions: $(BUILD_DIR)/cbp-profile-$(PRED_HASH) | out
	./$(BUILD_DIR)/cbp-profile-$(PRED_HASH) --format csv --mode acc --no-score --regions $(TRACE) $(TRACE_NAME) $(WARMUP) $(MEASURE) 1> out/profile.csv 2> out/profile_acc.txt
	@echo "=== Accuracy Analysis with Regions ===" && tail -40 out/profile_acc.txt

cbp-profile-analyze: $(BUILD_DIR)/cbp-profile-$(PRED_HASH) | out
	./$(BUILD_DIR)/cbp-profile-$(PRED_HASH) --format csv --mode analyze --profile $(TRACE) $(TRACE_NAME) $(WARMUP) $(MEASURE) 1> out/profile.csv 2> out/profile_analyze.txt
	@echo "=== Per-Function Analysis ===" && tail -30 out/profile_analyze.txt

cbp-profile-analyze-regions: $(BUILD_DIR)/cbp-profile-$(PRED_HASH) | out
	./$(BUILD_DIR)/cbp-profile-$(PRED_HASH) --format csv --mode analyze --profile --regions $(TRACE) $(TRACE_NAME) $(WARMUP) $(MEASURE) 1> out/profile.csv 2> out/profile_analyze.txt
	@echo "=== Per-Function Analysis with Regions ===" && tail -40 out/profile_analyze.txt

# TAGE Monitor: build with SW instrumentation and run
# Requires CHEATING_MODE for val<> data extraction. Zero cost when not enabled.
MONITOR_FLAGS := -DTAGE_MONITOR -DCHEATING_MODE -DFREE_FANOUT
MONITOR_OUT ?= out/monitor.txt

$(BUILD_DIR)/cbp-monitor-$(PRED_HASH): cbp.cpp cbp.hpp branch_predictor.hpp $(TRACE_READER) harcom.hpp $(wildcard predictors/*.hpp predictors/custom/*.hpp) $(PREDICTOR_MK) | $(BUILD_DIR)
	$(CXX) $(COMMON_FLAGS) $(EXTRA_COMMON_FLAGS) $(PROFILE_WARN_FLAGS) $(MONITOR_FLAGS) -Itrace_files -o $@ cbp.cpp -lz -DPREDICTOR='$(PREDICTOR_TYPE)'

MONITOR_MEASURE ?= $(MEASURE)
PROFILE_MEASURE ?= 100000

cbp-monitor: $(BUILD_DIR)/cbp-monitor-$(PRED_HASH) $(if $(WITH_PROFILE),$(BUILD_DIR)/cbp-profile-$(PRED_HASH)) | out
	@echo "=== Running Monitor ($(MONITOR_MEASURE) instr) ==="
ifdef WITH_PROFILE
	./$(BUILD_DIR)/cbp-profile-$(PRED_HASH) --format csv --mode analyze --profile $(TRACE) $(TRACE_NAME) $(WARMUP) $(PROFILE_MEASURE) 1> out/profile.csv 2> out/profile_analyze.txt &
endif
	./$(BUILD_DIR)/cbp-monitor-$(PRED_HASH) $(TRACE) $(TRACE_NAME) $(WARMUP) $(MONITOR_MEASURE) 1> out/monitor_csv.txt 2> $(MONITOR_OUT)
ifdef WITH_PROFILE
	wait
endif
	@echo ""
	@echo "=== Monitor ===" && cat $(MONITOR_OUT)
ifdef WITH_PROFILE
	@echo ""
	@echo "=== Profile ===" && tail -15 out/profile_analyze.txt
endif
	@echo ""
	@echo "Files: $(MONITOR_OUT), out/monitor_csv.txt"

monitor-vis:
	python3 scripts/monitor_vis.py $(MONITOR_OUT)

# Evaluate monitor on the representative trace set (same 20 traces as quick-eval)
EVAL_MONITOR_PRED ?= TageAhead1C
EVAL_MONITOR_JOBS ?= $(shell nproc)
EVAL_MONITOR_OUT ?= out/eval_monitor

$(BUILD_DIR)/cbp-monitor-eval: cbp.cpp cbp.hpp branch_predictor.hpp $(TRACE_READER) harcom.hpp $(wildcard predictors/*.hpp predictors/custom/*.hpp) | $(BUILD_DIR)
	$(CXX) $(COMMON_FLAGS) $(EXTRA_COMMON_FLAGS) $(PROFILE_WARN_FLAGS) $(MONITOR_FLAGS) -Itrace_files -o $@ cbp.cpp -lz -DPREDICTOR='$(EVAL_MONITOR_PRED)'

eval-monitor: $(BUILD_DIR)/cbp-monitor-eval | out
	scripts/eval_monitor.sh ./$(BUILD_DIR)/cbp-monitor-eval $(TRACE_DIR) $(EVAL_MONITOR_OUT) $(EVAL_MONITOR_JOBS)

eval-monitor-vis:
	python3 scripts/monitor_vis.py --dir $(EVAL_MONITOR_OUT)

# Debug print: build+run with DEBUG_PRINT to dump timing of all main signals
DEBUG_PRINT_MEASURE ?= 20

$(BUILD_DIR)/debug-print-$(PRED_HASH): cbp.cpp cbp.hpp branch_predictor.hpp $(TRACE_READER) harcom.hpp $(wildcard predictors/*.hpp predictors/custom/*.hpp) $(PREDICTOR_MK) | $(BUILD_DIR)
	$(CXX) $(COMMON_FLAGS) $(EXTRA_COMMON_FLAGS) $(CBP_WARN_FLAGS) -DDEBUG_PRINT $(EXTRA_CBP_FLAGS) -Itrace_files -o $@ cbp.cpp -lz -DPREDICTOR='$(PREDICTOR_TYPE)'

debug-print: $(BUILD_DIR)/debug-print-$(PRED_HASH)

run-debug-print: $(BUILD_DIR)/debug-print-$(PRED_HASH)
	./$(BUILD_DIR)/debug-print-$(PRED_HASH) $(TRACE) $(TRACE_NAME) 2 $(DEBUG_PRINT_MEASURE) 2>&1

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

# ---- Competition track builds: 1-cycle and 2-cycle ----
cbp-1c: | $(BUILD_DIR)
	$(CXX) $(COMMON_FLAGS) $(EXTRA_COMMON_FLAGS) $(CBP_WARN_FLAGS) $(EXTRA_CBP_FLAGS) -Itrace_files -o $(BUILD_DIR)/cbp-1c cbp.cpp -lz -DPREDICTOR='TageAhead1C'

cbp-2c: | $(BUILD_DIR)
	$(CXX) $(COMMON_FLAGS) $(EXTRA_COMMON_FLAGS) $(CBP_WARN_FLAGS) $(EXTRA_CBP_FLAGS) -Itrace_files -o $(BUILD_DIR)/cbp-2c cbp.cpp -lz -DPREDICTOR='TageAhead2C'

cbp-both: cbp-1c cbp-2c

quick-eval-1c: cbp-1c
	scripts/quick_eval.sh ./$(BUILD_DIR)/cbp-1c $(TRACE_DIR) out/quick_1c $(QUICK_JOBS)

quick-eval-2c: cbp-2c
	scripts/quick_eval.sh ./$(BUILD_DIR)/cbp-2c $(TRACE_DIR) out/quick_2c $(QUICK_JOBS)

quick-eval-both: | $(BUILD_DIR)
	@echo "=== Building 1-cycle ==="
	$(CXX) $(COMMON_FLAGS) $(EXTRA_COMMON_FLAGS) $(CBP_WARN_FLAGS) $(EXTRA_CBP_FLAGS) -Itrace_files -o $(BUILD_DIR)/cbp-1c cbp.cpp -lz -DPREDICTOR='TageAhead1C'
	@echo "=== Building 2-cycle ==="
	$(CXX) $(COMMON_FLAGS) $(EXTRA_COMMON_FLAGS) $(CBP_WARN_FLAGS) $(EXTRA_CBP_FLAGS) -Itrace_files -o $(BUILD_DIR)/cbp-2c cbp.cpp -lz -DPREDICTOR='TageAhead2C'
	@echo "=== Evaluating 1-cycle ==="
	scripts/quick_eval.sh ./$(BUILD_DIR)/cbp-1c $(TRACE_DIR) out/quick_1c $(QUICK_JOBS)
	@echo ""
	@echo "=== Evaluating 2-cycle ==="
	scripts/quick_eval.sh ./$(BUILD_DIR)/cbp-2c $(TRACE_DIR) out/quick_2c $(QUICK_JOBS)

debug-print-1c: | $(BUILD_DIR)
	$(CXX) $(COMMON_FLAGS) $(EXTRA_COMMON_FLAGS) $(CBP_WARN_FLAGS) -DDEBUG_PRINT $(EXTRA_CBP_FLAGS) -Itrace_files -o $(BUILD_DIR)/debug-print-1c cbp.cpp -lz -DPREDICTOR='TageAhead1C'

debug-print-2c: | $(BUILD_DIR)
	$(CXX) $(COMMON_FLAGS) $(EXTRA_COMMON_FLAGS) $(CBP_WARN_FLAGS) -DDEBUG_PRINT $(EXTRA_CBP_FLAGS) -Itrace_files -o $(BUILD_DIR)/debug-print-2c cbp.cpp -lz -DPREDICTOR='TageAhead2C'

run-debug-print-1c: debug-print-1c
	./$(BUILD_DIR)/debug-print-1c $(TRACE) $(TRACE_NAME) 2 $(DEBUG_PRINT_MEASURE) 2>&1

run-debug-print-2c: debug-print-2c
	./$(BUILD_DIR)/debug-print-2c $(TRACE) $(TRACE_NAME) 2 $(DEBUG_PRINT_MEASURE) 2>&1

# Quick evaluation on representative trace subset (20 traces)
QUICK_JOBS ?= $(shell nproc)
QUICK_OUT ?= out/quick

quick-eval: $(BUILD_DIR)/cbp-$(PRED_HASH)
	scripts/quick_eval.sh ./$(BUILD_DIR)/cbp-$(PRED_HASH) $(TRACE_DIR) $(QUICK_OUT) $(QUICK_JOBS)

# Energy breakdown: logic vs fanout vs wiring vs RAM
# Usage: make diff-energy PRED_A=TageAheadHC_IR PRED_B=example_tage
# Single predictor: make diff-energy PRED_A=TageAheadHC_IR
DIFF_ENERGY_JOBS ?= $(QUICK_JOBS)
diff-energy:
	$(PYTHON) scripts/energy_breakdown.py '$(PRED_A)' $(if $(PRED_B),'$(PRED_B)') \
		--traces $(TRACE_DIR) --jobs $(DIFF_ENERGY_JOBS) \
		--warmup $(WARMUP) --measure $(MEASURE) \
		--build-dir $(BUILD_DIR)

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

GRADIENT_CONFIG ?= configs/gradient_default.yaml
GRADIENT_JOBS ?= $(shell nproc)
GRADIENT_PERTURBATION ?= 1

gradient:
	mkdir -p out/gradient
	$(PYTHON) scripts/gradient_ascent.py \
		--config $(GRADIENT_CONFIG) \
		--trace-dir $(TRACE_DIR) \
		-j $(GRADIENT_JOBS) \
		--perturbation $(GRADIENT_PERTURBATION) \
		--extra-flags '$(EXTRA_COMMON_FLAGS) $(EXTRA_CBP_FLAGS)' --resume

gradient-report:
	$(PYTHON) scripts/gradient_ascent.py --report out/gradient/results.csv --top 20

# TageAhead iterative parameter sweep (ALWAYS x PERTURBED, two-phase timing+eval)
PSWEEP_CONFIG ?= sweep_config.json
PSWEEP_JOBS ?= $(shell nproc)
PSWEEP_TOP ?= 20

psweep-init:
	$(PYTHON) scripts/param_sweep.py --config $(PSWEEP_CONFIG) --init

psweep-generate:
	$(PYTHON) scripts/param_sweep.py --config $(PSWEEP_CONFIG) --generate -j $(PSWEEP_JOBS)

psweep-timing:
	$(PYTHON) scripts/param_sweep.py --config $(PSWEEP_CONFIG) --timing -j $(PSWEEP_JOBS)

psweep-eval:
	$(PYTHON) scripts/param_sweep.py --config $(PSWEEP_CONFIG) --eval -j $(PSWEEP_JOBS)

psweep-run:
	$(PYTHON) scripts/param_sweep.py --config $(PSWEEP_CONFIG) --run -j $(PSWEEP_JOBS)

psweep-resume:
	$(PYTHON) scripts/param_sweep.py --config $(PSWEEP_CONFIG) --resume -j $(PSWEEP_JOBS)

psweep-report:
	$(PYTHON) scripts/param_sweep.py --config $(PSWEEP_CONFIG) --report --top $(PSWEEP_TOP)

SAVE_DIR := $(BUILD_DIR)/saved
SAVE_NAME ?= $(shell date +%Y%m%d_%H%M%S)

$(SAVE_DIR):
	mkdir -p $@

save: $(BUILD_DIR)/cbp-$(PRED_HASH) | $(SAVE_DIR)
	cp $(BUILD_DIR)/cbp-$(PRED_HASH) $(SAVE_DIR)/$(SAVE_NAME)
	@echo "Saved to $(SAVE_DIR)/$(SAVE_NAME)"

list-saved:
	@ls -lt $(SAVE_DIR)/ 2>/dev/null || echo "No saved binaries"

clean:
	rm -f $(BUILD_DIR)/cbp-* $(BUILD_DIR)/reference $(BUILD_DIR)/trace-analyze
	rm -f $(BUILD_DIR)/block-analyze $(BUILD_DIR)/ahead-block-analyze
	rm -f $(BUILD_DIR)/test-tage-compile $(BUILD_DIR)/quick_a $(BUILD_DIR)/quick_b
	rm -f $(BUILD_DIR)/timing-probe-*
	rm -f $(BUILD_DIR)/debug-print-*
	rm -f $(PREDICTOR_MK)

# Clean individual targets: make clean-cbp, clean-monitor, clean-profile, clean-out
clean-cbp:
	rm -f $(BUILD_DIR)/cbp-[0-9a-f]* $(PREDICTOR_MK)
clean-monitor:
	rm -f $(BUILD_DIR)/cbp-monitor-*
clean-profile:
	rm -f $(BUILD_DIR)/cbp-profile-*
clean-reference:
	rm -f $(BUILD_DIR)/reference
clean-out:
