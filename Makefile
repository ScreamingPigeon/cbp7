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
EXTRA_COMMON_FLAGS ?=
EXTRA_CBP_FLAGS ?=

-include $(PREDICTOR_MK)
PREDICTOR_TYPE ?= tage<>
TRACE ?= ./gcc_test_trace.gz
TRACE_NAME ?= test
WARMUP ?= 1000000
MEASURE ?= 40000000

.PHONY: all help cbp reference predictor-config run-cbp run-reference harness clean

all: cbp reference

help:
	@echo "Targets:"
	@echo "  make cbp              Build CBP simulator (uses $(PARAMS_FILE) if present)"
	@echo "  make reference        Build CBP2025 reference predictor"
	@echo "  make run-cbp          Run cbp on TRACE with direct args"
	@echo "  make run-reference    Run reference on TRACE with direct args"
	@echo "  make harness          Analyze hardware geometry of predictor"
	@echo "  make predictor-config Generate $(PREDICTOR_MK) from $(PARAMS_FILE)"
	@echo "  make clean            Remove generated build artifacts"
	@echo
	@echo "Variables you can override:"
	@echo "  TRACE=... TRACE_NAME=... WARMUP=... MEASURE=..."
	@echo "  (all default from $(PARAMS_FILE); CLI values still override)"
	@echo "  EXTRA_COMMON_FLAGS=... EXTRA_CBP_FLAGS=..."
	@echo "  PARAMS_FILE=... PREDICTOR_TYPE=... CXX=... PYTHON=..."
	@echo
	@echo "Harness targets:"
	@echo "  make harness PREDICTOR_TYPE='Custom<>'  Analyze Custom predictor"
	@echo "  make harness PREDICTOR_TYPE='tage<>'    Analyze TAGE predictor"
	@echo "  Optional: add --detailed flag for verbose output"

$(BUILD_DIR):
	mkdir -p $@

predictor-config: $(PREDICTOR_MK)

$(PREDICTOR_MK): scripts/gen_predictor_config.py $(PARAMS_FILE) | $(BUILD_DIR)
	$(PYTHON) scripts/gen_predictor_config.py --input $(PARAMS_FILE) --output $@

cbp: cbp.cpp cbp.hpp branch_predictor.hpp trace_reader.hpp harcom.hpp $(wildcard predictors/*.hpp) $(PREDICTOR_MK)
	$(CXX) $(COMMON_FLAGS) $(EXTRA_COMMON_FLAGS) $(CBP_WARN_FLAGS) $(EXTRA_CBP_FLAGS) -o $@ cbp.cpp -lz -DPREDICTOR='$(PREDICTOR_TYPE)'

reference: reference.cpp trace_reader.hpp seznec_cbp2025.h
	$(CXX) $(COMMON_FLAGS) $(EXTRA_COMMON_FLAGS) $(REFERENCE_WARN_FLAGS) -o $@ reference.cpp -lz

run-cbp: cbp
	./cbp $(TRACE) $(TRACE_NAME) $(WARMUP) $(MEASURE) > out/$(TRACE_NAME).out

run-reference: reference
	./reference $(TRACE) $(TRACE_NAME) $(WARMUP) $(MEASURE)

harness: build/harness_executable
	@echo "=== HARCOM Hardware Geometry Analysis ==="
	@echo "Predictor: $(PREDICTOR_TYPE)"
	@echo
	./build/harness_executable

build/harness_executable: agent/hw_lib/harness_instantiate.cpp agent/hw_lib/harness_main.cpp cbp.hpp branch_predictor.hpp harcom.hpp $(wildcard predictors/*.hpp) $(wildcard agent/hw_lib/*.hpp) | $(BUILD_DIR)
	$(CXX) $(COMMON_FLAGS) $(EXTRA_COMMON_FLAGS) $(CBP_WARN_FLAGS) $(EXTRA_CBP_FLAGS) \
	  -o $@ agent/hw_lib/harness_instantiate.cpp agent/hw_lib/harness_main.cpp \
	  -lz -DPREDICTOR='$(PREDICTOR_TYPE)'

clean:
	rm -f cbp reference build/harness_executable
	rm -f $(PREDICTOR_MK)
