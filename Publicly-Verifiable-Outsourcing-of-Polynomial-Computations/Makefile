# Publicly Verifiable Outsourcing of Polynomial Computations
# Build only source files that are actually present in this artifact.
CC ?= gcc
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
ifneq ($(PV_PREFIX),)
    FLINT_INCLUDE = -I$(PV_PREFIX)/include -I$(PV_PREFIX)/include/flint
    FLINT_LIB_PATH = -L$(PV_PREFIX)/lib -Wl,-rpath,$(PV_PREFIX)/lib
endif
CFLAGS ?= -Wall -O2 -std=c99
CFLAGS += $(FLINT_INCLUDE) -Ilib/common -Iprotocols
LIBS = $(FLINT_LIB_PATH) -lflint -lgmp -lsodium -lm
LIB_DIR = lib/common
PROTOCOLS_DIR = protocols
BUILD_DIR = build
COMMON_OBJ = $(BUILD_DIR)/mpoly.o $(BUILD_DIR)/crypto_util.o $(BUILD_DIR)/extended_field.o $(BUILD_DIR)/hss_paillier.o
CORE_OBJ = $(BUILD_DIR)/gamma_fe_sh.o $(BUILD_DIR)/gamma_ov_sh.o $(BUILD_DIR)/gamma_ov_hss.o $(BUILD_DIR)/protocol2.o $(BUILD_DIR)/pi3.o
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
	python3 scripts/validate_artifact.py

check: validate all test

benchmark: $(BUILD_DIR)/benchmark

client-benchmark: $(BUILD_DIR)/client_benchmark

var_benchmark: $(BUILD_DIR)/benchmark
	$(BUILD_DIR)/benchmark Gamma_OV_Sh 4 2 1 3 var

reported-plots:
	python3 scripts/plot_reported_results.py

VM_RESULTS ?= results/measured/vm-paper
VM_PLOTS ?= results/plots/vm-paper

figures:
	python3 scripts/plot_vm_results.py $(VM_RESULTS) --output $(VM_PLOTS)

figure2 figure3 figure4 figure5:
	python3 scripts/plot_vm_results.py $(VM_RESULTS) --output $(VM_PLOTS) --figures $(patsubst figure%,%,$@)

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

$(BUILD_DIR)/client_benchmark: benchmarks/client_benchmark.c $(COMMON_OBJ) $(CORE_OBJ)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

$(BUILD_DIR)/vm_wire.o: benchmarks/vm/wire.c benchmarks/vm/wire.h | $(BUILD_DIR)
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

$(BUILD_DIR)/vm_client: benchmarks/vm/client.c $(BUILD_DIR)/vm_wire.o $(COMMON_OBJ) $(CORE_OBJ)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

$(BUILD_DIR)/vm_worker: benchmarks/vm/worker.c $(BUILD_DIR)/vm_wire.o $(COMMON_OBJ) $(CORE_OBJ)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

$(BUILD_DIR)/microbenchmark: benchmarks/microbenchmark.c $(BUILD_DIR)/vm_wire.o $(COMMON_OBJ) $(CORE_OBJ)
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

experiments: $(BUILD_DIR)/vm_client $(BUILD_DIR)/vm_worker $(BUILD_DIR)/microbenchmark
	python3 scripts/write_build_manifest.py --compiler '$(CC)' --cflags '$(CFLAGS)' --libs '$(LIBS)'

clean:
	rm -rf $(BUILD_DIR)

-include $(COMMON_OBJ:.o=.d) $(CORE_OBJ:.o=.d) $(BUILD_DIR)/vm_wire.d
.PHONY: experiments all core test test-hss test-hss-high validate check benchmark client-benchmark var_benchmark reported-plots figures figure2 figure3 figure4 figure5 clean
