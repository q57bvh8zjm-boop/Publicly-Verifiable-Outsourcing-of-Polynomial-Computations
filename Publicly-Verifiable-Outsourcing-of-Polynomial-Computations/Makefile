# Publicly Verifiable Outsourcing of Polynomial Computations
# Core protocols and fixed public polynomial classifiers.
CC ?= gcc
PYTHON ?= python3
CLASSIFICATION_OUTPUT ?= results/measured/local-classification
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    ifneq (,$(wildcard /opt/homebrew/include/flint/flint.h))
        FLINT_INCLUDE = -I/opt/homebrew/include -I/opt/homebrew/include/flint
        FLINT_LIB_PATH = -L/opt/homebrew/lib
    else ifneq (,$(wildcard /usr/local/include/flint/flint.h))
        FLINT_INCLUDE = -I/usr/local/include -I/usr/local/include/flint
        FLINT_LIB_PATH = -L/usr/local/lib
    else
        FLINT_INCLUDE = -I/usr/include/flint
        FLINT_LIB_PATH =
    endif
else
    FLINT_INCLUDE = -I/usr/include/flint
    FLINT_LIB_PATH =
endif
CFLAGS ?= -Wall -O2 -std=c99
CFLAGS += $(FLINT_INCLUDE) -Ilib/common -Iprotocols
LIBS = $(FLINT_LIB_PATH) -lflint -lgmp -lsodium -lm
LIB_DIR = lib/common
PROTOCOLS_DIR = protocols
BUILD_DIR = build
COMMON_OBJ = $(BUILD_DIR)/mpoly.o $(BUILD_DIR)/crypto_util.o $(BUILD_DIR)/extended_field.o $(BUILD_DIR)/hss_paillier.o
CORE_OBJ = $(BUILD_DIR)/gamma_fe_sh.o $(BUILD_DIR)/gamma_ov_sh.o $(BUILD_DIR)/gamma_ov_hss.o $(BUILD_DIR)/protocol2.o
all: core $(BUILD_DIR)/benchmark

core: $(COMMON_OBJ) $(CORE_OBJ)

test: $(BUILD_DIR)/test_mpoly $(BUILD_DIR)/test_protocols
	$(BUILD_DIR)/test_mpoly
	$(BUILD_DIR)/test_protocols

test-hss: $(BUILD_DIR)/test_protocols
	$(BUILD_DIR)/test_protocols --hss

test-hss-high: $(BUILD_DIR)/test_protocols
	$(BUILD_DIR)/test_protocols --hss-high

validate:
	$(PYTHON) scripts/validate_artifact.py

check: validate all test

benchmark: $(BUILD_DIR)/benchmark

reported-plots:
	$(PYTHON) scripts/plot_reported_results.py

figures figure2 figure3 figure4 figure5: reported-plots

$(BUILD_DIR):
	mkdir -p $@

$(BUILD_DIR)/%.o: $(LIB_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

$(BUILD_DIR)/%.o: $(PROTOCOLS_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

$(BUILD_DIR)/test_mpoly: tests/test_mpoly.c $(COMMON_OBJ) $(CORE_OBJ)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

$(BUILD_DIR)/test_protocols: tests/test_protocols.c $(COMMON_OBJ) $(CORE_OBJ)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

$(BUILD_DIR)/benchmark: benchmarks/benchmark.c $(COMMON_OBJ) $(CORE_OBJ)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
	rm -rf $(BUILD_DIR)

-include $(COMMON_OBJ:.o=.d) $(CORE_OBJ:.o=.d)
.PHONY: all core test test-hss test-hss-high validate check benchmark reported-plots figures figure2 figure3 figure4 figure5 clean

# Compile the classification adapter with the -O3 setting used for Table 7.
CLASSIFICATION_SRC = applications/private-classification/src/classification_bench.c $(LIB_DIR)/crypto_util.c $(LIB_DIR)/mpoly.c $(LIB_DIR)/extended_field.c $(PROTOCOLS_DIR)/gamma_fe_sh.c $(PROTOCOLS_DIR)/gamma_ov_sh.c
CLASSIFICATION_HEADERS = $(wildcard $(LIB_DIR)/*.h) $(PROTOCOLS_DIR)/gamma_fe_sh.h $(PROTOCOLS_DIR)/gamma_ov_sh.h

$(BUILD_DIR)/classification_bench: $(CLASSIFICATION_SRC) $(CLASSIFICATION_HEADERS) | $(BUILD_DIR)
	$(CC) $(CFLAGS) -O3 $(CLASSIFICATION_SRC) -o $@ $(LIBS)

classification-check classification-benchmark: $(BUILD_DIR)/classification_bench
	$(PYTHON) applications/private-classification/run.py --output $(CLASSIFICATION_OUTPUT)

.PHONY: classification-check classification-benchmark
