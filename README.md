# Segmented Sieve of Eratosthenes with Wheel-30 Optimization

This repository provides an implementation of the segmented Sieve of Eratosthenes enhanced with Wheel-30 optimization.  
The algorithm is designed for large-scale prime enumeration and prime-counting experiments, verified up to 1e11 and  
potentially scalable to 1e12.  

It is intended for reproducibility in computational number theory research, offering a clear reference implementation  
that balances correctness, memory efficiency, and scalability.

---

## Features

- **Segmented sieve**: Enables sieving in fixed-size blocks, reducing memory usage and allowing scalability to very large bounds.
- **Wheel-30 optimization**: Skips multiples of 2, 3, and 5, reducing sieve operations by a factor of 8.
- **Checkpoint system**: Supports prime-counting validation against known values of π(x).
- **Progress logging**: Records segment-by-segment progress with elapsed time and ETA.
- **Meta file output**: Saves checkpoints and progress logs in `meta.txt` for reproducibility.

---

## Usage

### Build
```bash
cargo build --release
```

### Run
```bash
cargo run --release -- <n> <step>
```

- `<n>`: Upper bound of sieve (default: `1e9`).
- `<step>`: Interval for writing π(x) values to `pi_index.bin` (default: `1e6`).

### Example
```bash
cargo run --release -- 1000000000 1000000
```

This computes all primes up to 1e9, writing:
- `primes.bin`: Binary list of primes (little-endian `u32`).
- `pi_index.bin`: Prime-count function sampled every `<step>`.
- `sieve_w30.bitset`: Bitset representation of the sieve.
- `meta.txt`: Human-readable log of progress and checkpoints.

---

## Checkpoints

The implementation includes known prime-counting values for verification:

- π(1,000,000)   = 78,498  
- π(10,000,000)  = 664,579  
- π(100,000,000) = 5,761,455  
- π(1,000,000,000) = 50,847,534  
- π(10,000,000,000) = 455,052,511  

These are automatically checked during execution.

---

## License

This project is licensed under the MIT License.  
Copyright (c) 2025 Hannah Cheng
