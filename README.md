# Segmented Eratosthenes sieve with a wheel of 30

This repository contains a Rust implementation of a segmented Sieve of Eratosthenes. It keeps only candidates coprime to `2*3*5=30`, using the eight wheel residues `1, 7, 11, 13, 17, 19, 23, 29`. It is intended for reproducible prime-enumeration experiments, not as a replacement for a specialized prime-counting library.

## Algorithm

1. A conventional sieve generates base primes up to `floor(sqrt(n))`.
2. `[0,n]` is processed in fixed-size segments; the default segment contains 9,000,000 integers.
3. For each base prime `p >= 7`, multiples starting at `max(p*p, ceil(low/p)*p)` are marked. Only wheel residues are touched.
4. Unmarked candidates are emitted; 2, 3, and 5 are emitted separately.

The wheel removes seven of every thirty positions before sieving. The segmented sieve still has `O(n log log n)` marking work and `O(sqrt(n) + segment_size)` working memory.

## Build and run

```sh
cargo build --release
cargo run --release -- <n> <step>
```

`n` is the inclusive upper bound (default `1_000_000_000`); `step` is the sampling interval for `pi_index.bin` (default `1_000_000`). Example:

```sh
cargo run --release -- 10000000 1000000
```

Outputs are written to the current directory:

| File | Format |
| --- | --- |
| `primes.bin` | increasing primes as little-endian `u32` values |
| `pi_index.bin` | running prime-count samples as little-endian `u64` values |
| `sieve_w30.bitset` | `W30BIT\\0` header followed by one wheel byte per 30-block; bit 1 means composite |
| `meta.txt` | parameters, checkpoint results, progress, and final `pi(n)` |

## Correctness checks

Known values checked at runtime when the selected bound reaches them:

| x | pi(x) |
| ---: | ---: |
| 1,000,000 | 78,498 |
| 10,000,000 | 664,579 |
| 100,000,000 | 5,761,455 |
| 1,000,000,000 | 50,847,534 |
| 10,000,000,000 | 455,052,511 |

These checks are not a proof for every input. The repository currently has no unit tests; `cargo test --release` checks compilation and reports zero test cases.

## MacBook test

Measured on the repository's MacBook environment on 2026-09-12:

| Command | Result | Wall time |
| --- | --- | ---: |
| `cargo run --release -- 1000000 100000` | `pi(1,000,000)=78,498` | 0.85 s (including first release build) |
| `cargo run --release -- 10000000 1000000` | `pi(10,000,000)=664,579`; 1M checkpoint OK | 0.19 s |

The second run reused the release binary. Times vary with hardware, storage, and background load.

## Numeric limit

The bound is parsed as `u64`, but `primes.bin` stores primes as `u32`. The current output format is therefore safe only while every emitted prime fits in `u32` (`n <= u32::MAX`). This repository does not claim a verified `1e11` run. Larger bounds require changing the output format and auditing every `u32` conversion.

## License

MIT. Copyright (c) 2025 Hannah Cheng.
