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
CFLAGS ?= -Wall -O2 -std=c99
CFLAGS += $(FLINT_INCLUDE) -Ilib/common -Iprotocols
LIBS = $(FLINT_LIB_PATH) -lflint -lgmp -lsodium -lm
LIB_DIR = lib/common
PROTOCOLS_DIR = protocols
BUILD_DIR = build
COMMON_OBJ = $(BUILD_DIR)/mpoly.o $(BUILD_DIR)/extended_field.o
CORE_OBJ = $(BUILD_DIR)/scheme1.o $(BUILD_DIR)/scheme2.o $(BUILD_DIR)/scheme3.o $(BUILD_DIR)/protocol2.o
all: dirs $(COMMON_OBJ) $(CORE_OBJ) scheme3_pir_optimize protocol2_pir_optimize
	@echo "Core schemes and exploratory PIR programs compiled successfully."

test: dirs $(BUILD_DIR)/test_mpoly
	$(BUILD_DIR)/test_mpoly

validate:
	python3 scripts/validate_artifact.py

check: validate test all
	@echo "Artifact validation, tests, and build completed successfully."
core: dirs $(COMMON_OBJ) $(CORE_OBJ)
	@echo "Core scheme objects compiled successfully."
scheme3_pir_optimize: dirs $(COMMON_OBJ) $(BUILD_DIR)/scheme3.o $(PROTOCOLS_DIR)/scheme3_pir_optimize.c
	$(CC) $(CFLAGS) $(PROTOCOLS_DIR)/scheme3_pir_optimize.c $(BUILD_DIR)/scheme3.o $(COMMON_OBJ) -o $(BUILD_DIR)/scheme3_pir_optimize $(LIBS)
protocol2_pir_optimize: dirs $(COMMON_OBJ) $(BUILD_DIR)/protocol2.o $(PROTOCOLS_DIR)/protocol2_pir_optimize.c
	$(CC) $(CFLAGS) $(PROTOCOLS_DIR)/protocol2_pir_optimize.c $(BUILD_DIR)/protocol2.o $(COMMON_OBJ) -o $(BUILD_DIR)/protocol2_pir_optimize $(LIBS)
var_benchmark: dirs $(COMMON_OBJ) $(BUILD_DIR)/scheme1.o $(BUILD_DIR)/scheme2.o $(BUILD_DIR)/scheme3.o var_benchmark.c
	$(CC) $(CFLAGS) var_benchmark.c $(BUILD_DIR)/scheme1.o $(BUILD_DIR)/scheme2.o $(BUILD_DIR)/scheme3.o $(COMMON_OBJ) -o $(BUILD_DIR)/var_benchmark $(LIBS)
reported-plots:
	python3 scripts/plot_reported_results.py

figures: reported-plots
	@echo "Generated Figures 2, 3, 4, and 5 under results/plots/."

figure2: reported-plots
	@echo "Generated results/plots/figure2_low_degree.png"

figure3: reported-plots
	@echo "Generated results/plots/figure3_thresholds.png"

figure4: reported-plots
	@echo "Generated results/plots/figure4_d4.png, figure4_d8.png, and figure4_d10.png"

figure5: reported-plots
	@echo "Generated results/plots/figure5_gamma1_vs_pi3.png"
dirs:
	@mkdir -p $(BUILD_DIR)
$(BUILD_DIR)/mpoly.o: $(LIB_DIR)/mpoly.c $(LIB_DIR)/mpoly.h
	$(CC) $(CFLAGS) -c $< -o $@
$(BUILD_DIR)/extended_field.o: $(LIB_DIR)/extended_field.c $(LIB_DIR)/extended_field.h
	$(CC) $(CFLAGS) -c $< -o $@
$(BUILD_DIR)/scheme1.o: $(PROTOCOLS_DIR)/scheme1.c $(PROTOCOLS_DIR)/scheme1.h $(LIB_DIR)/mpoly.h $(LIB_DIR)/extended_field.h
	$(CC) $(CFLAGS) -c $< -o $@
$(BUILD_DIR)/scheme2.o: $(PROTOCOLS_DIR)/scheme2.c $(PROTOCOLS_DIR)/scheme2.h $(LIB_DIR)/mpoly.h
	$(CC) $(CFLAGS) -c $< -o $@
$(BUILD_DIR)/scheme3.o: $(PROTOCOLS_DIR)/scheme3.c $(PROTOCOLS_DIR)/scheme3.h $(LIB_DIR)/mpoly.h
	$(CC) $(CFLAGS) -c $< -o $@
$(BUILD_DIR)/protocol2.o: $(PROTOCOLS_DIR)/protocol2.c $(PROTOCOLS_DIR)/protocol2.h
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/test_mpoly: tests/test_mpoly.c $(COMMON_OBJ) $(BUILD_DIR)/scheme1.o $(BUILD_DIR)/protocol2.o
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)
clean:
	rm -rf $(BUILD_DIR)
	@mkdir -p $(BUILD_DIR)
	@touch $(BUILD_DIR)/.gitkeep
	@echo "Cleanup complete."
.PHONY: all check core dirs clean figure2 figure3 figure4 figure5 figures reported-plots test validate var_benchmark
