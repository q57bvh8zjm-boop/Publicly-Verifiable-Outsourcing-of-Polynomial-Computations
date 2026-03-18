# Artifacts for Paper Evaluation

This document provides a map of the artifacts associated with the paper's contribution and instructions on how to access and evaluate them.

## 1. Core Contribution: Scheme Implementations
The primary artifacts are the C implementations of the four schemes discussed in the paper. These are located in the `protocols/` directory:

- **Scheme 1 (Variable Input Dimensions, Private Verification)**  
  - Path: [`protocols/scheme1.c`](./protocols/scheme1.c), [`protocols/scheme1.h`](./protocols/scheme1.h)
  - Contribution: Efficient private verification using extended field techniques.

- **Scheme 2 (Variable Input Dimensions, Public Verification)**  
  - Path: [`protocols/scheme2.c`](./protocols/scheme2.c), [`protocols/scheme2.h`](./protocols/scheme2.h)
  - Contribution: Support for public verification using cryptographic group operations.

- **Scheme 3 (Low-Degree Polynomials via HSS)**  
  - Path: [`protocols/scheme3.c`](./protocols/scheme3.c), [`protocols/scheme3.h`](./protocols/scheme3.h)
  - Contribution: Leveraging Paillier-based Homomorphic Secret Sharing for polynomial evaluation.

- **Scheme 4 (Combinatorial Grouping)**  
  - Path: [`protocols/scheme4.c`](./protocols/scheme4.c), [`protocols/scheme4.h`](./protocols/scheme4.h)
  - Contribution: Optimization of server count and computation via combinatorial covering designs.

## 2. Shared Libraries and Dependencies
The schemes rely on a common set of utility functions located in `lib/common/`:

- **Multivariate Polynomials**: [`lib/common/mpoly.c`](./lib/common/mpoly.c) - Core logic for polynomial representation and operations.
- **Extended Field Arithmetic**: [`lib/common/extended_field.c`](./lib/common/extended_field.c) - Optimized arithmetic for $F_{q^2}$ used in Scheme 1.

## 3. Build and Evaluation System
The project includes a unified build system to verify the integrity and compilability of the code:

- **Build Script**: [`Makefile`](./Makefile) - Contains automated instructions to compile all schemes into verified object files.
- **Instructions**:
  ```bash
  # Compile all artifacts
  make all
  ```

## 4. Documentation
- **High-level Overview**: [`README.md`](./README.md) - Provides an overview of the project, environment requirements, and mapping between filenames and paper sections.

## How to Access
Reviewers can access each artifact by exploring the directory structure listed above. All paths are relative to the root of the research repository.
