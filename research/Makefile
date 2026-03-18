# Multi-server Low-degree Polynomial Verifiable Computation Project
# Simplified Makefile for Schemes 1, 2, 3, and 4

CC = gcc

# Detect system and set FLINT path
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

CFLAGS = -Wall -O2 -std=c99 $(FLINT_INCLUDE) -Ilib/common -Iprotocols
LIBS = $(FLINT_LIB_PATH) -lflint -lgmp -lsodium -lm

# Directories
LIB_DIR = lib/common
PROTOCOLS_DIR = protocols
BUILD_DIR = build

# Source files
COMMON_SRC = $(LIB_DIR)/mpoly.c $(LIB_DIR)/extended_field.c
SCHEME1_SRC = $(PROTOCOLS_DIR)/scheme1.c
SCHEME2_SRC = $(PROTOCOLS_DIR)/scheme2.c
SCHEME3_SRC = $(PROTOCOLS_DIR)/scheme3.c
SCHEME4_SRC = $(PROTOCOLS_DIR)/scheme4.c
OPTIMIZE_SRC = $(PROTOCOLS_DIR)/scheme3_pir_optimize.c
PROTOCOL2_SRC = $(PROTOCOLS_DIR)/protocol2.c
P2_OPTIMIZE_SRC = $(PROTOCOLS_DIR)/protocol2_pir_optimize.c

# Object files
COMMON_OBJ = $(BUILD_DIR)/mpoly.o $(BUILD_DIR)/extended_field.o
SCHEME_OBJS = $(BUILD_DIR)/scheme1.o $(BUILD_DIR)/scheme2.o $(BUILD_DIR)/scheme3.o $(BUILD_DIR)/scheme4.o $(BUILD_DIR)/protocol2.o

# Targets
all: dirs $(COMMON_OBJ) $(SCHEME_OBJS) scheme3_pir_optimize protocol2_pir_optimize
	@echo "All schemes compiled successfully as object files."

scheme3_pir_optimize: dirs $(COMMON_OBJ) $(BUILD_DIR)/scheme3.o $(OPTIMIZE_SRC)
	$(CC) $(CFLAGS) $(OPTIMIZE_SRC) $(BUILD_DIR)/scheme3.o $(COMMON_OBJ) -o $(BUILD_DIR)/scheme3_pir_optimize $(LIBS)

protocol2_pir_optimize: dirs $(COMMON_OBJ) $(BUILD_DIR)/protocol2.o $(P2_OPTIMIZE_SRC)
	$(CC) $(CFLAGS) $(P2_OPTIMIZE_SRC) $(BUILD_DIR)/protocol2.o $(COMMON_OBJ) -o $(BUILD_DIR)/protocol2_pir_optimize $(LIBS)

$(BUILD_DIR)/protocol2.o: $(PROTOCOLS_DIR)/protocol2.c $(PROTOCOLS_DIR)/protocol2.h
	$(CC) $(CFLAGS) -c $< -o $@

dirs:
	@mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%.o: $(LIB_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/scheme1.o: $(PROTOCOLS_DIR)/scheme1.c $(PROTOCOLS_DIR)/scheme1.h
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/scheme2.o: $(PROTOCOLS_DIR)/scheme2.c $(PROTOCOLS_DIR)/scheme2.h
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/scheme3.o: $(PROTOCOLS_DIR)/scheme3.c $(PROTOCOLS_DIR)/scheme3.h
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/scheme4.o: $(PROTOCOLS_DIR)/scheme4.c $(PROTOCOLS_DIR)/scheme4.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR)
	@echo "Cleanup complete."

.PHONY: all dirs clean