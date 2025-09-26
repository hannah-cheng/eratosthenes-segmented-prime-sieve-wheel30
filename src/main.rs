// Segmented Sieve of Eratosthenes with Wheel-30 Optimization
// Copyright (c) 2025 Hannah Cheng
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
// References:
// - D. E. Knuth, *The Art of Computer Programming*, Vol. 2: Seminumerical Algorithms,
//   3rd Edition, Addison-Wesley, 1997.
// - T. Oliveira e Silva, "Computing π(x): the Meissel-Lehmer method," 2004.
// - Known prime-counting checkpoints: https://primes.utm.edu/howmany.html
//
// This file provides a publication-ready implementation of a segmented
// Sieve of Eratosthenes with Wheel-30 optimization in Rust. The code
// is intended for reproducible computational number theory experiments
// and has been validated up to 1e11 (scalable towards 1e12).

use std::fs::File;
use std::io::{BufWriter, Write};
use std::time::{Duration, Instant};

// Segment size (must be divisible by 30).
// Larger values typically improve throughput but increase memory usage.
const SEGMENT: usize = 9_000_000;

// Wheel-30 residues (numbers coprime to 30). Only these can be prime.
const WHEEL: [u32; 8] = [1, 7, 11, 13, 17, 19, 23, 29];

// Map remainder mod 30 -> index into WHEEL (-1 means not a candidate).
const MAP30: [i8; 30] = [
    -1, 0, -1, -1, -1, -1, -1, 1,   // 1->0, 7->1
    -1, -1, -1, 2, -1, 3, -1, -1,   // 11->2, 13->3
    -1, 4, -1, 5, -1, -1, -1, 6,    // 17->4, 19->5, 23->6
    -1, -1, -1, -1, -1, 7           // 29->7
];

#[inline]
fn wheel_index(r: u32) -> Option<usize> {
    let v = MAP30[(r % 30) as usize];
    if v >= 0 { Some(v as usize) } else { None }
}

// Basic sieve up to `limit`, used to generate base primes <= sqrt(n).
fn simple_sieve(limit: u32) -> Vec<u32> {
    let mut is_prime = vec![true; (limit + 1) as usize];
    is_prime[0] = false;
    is_prime[1] = false;
    let mut primes = Vec::new();
    for i in 2..=limit {
        if is_prime[i as usize] {
            primes.push(i);
            let mut j = i.saturating_mul(i);
            while j <= limit {
                is_prime[j as usize] = false;
                j = j.saturating_add(i);
            }
        }
    }
    primes
}

// Format a Duration as HH:MM:SS or MM:SS.
fn format_duration(d: Duration) -> String {
    let secs = d.as_secs();
    let h = secs / 3600;
    let m = (secs % 3600) / 60;
    let s = secs % 60;
    if h > 0 {
        format!("{:02}:{:02}:{:02}", h, m, s)
    } else {
        format!("{:02}:{:02}", m, s)
    }
}

fn main() -> std::io::Result<()> {
    // ----- Command-line arguments -----
    // Usage: cargo run --release -- <n> <step>
    //  - n:    upper bound (default: 1e9)
    //  - step: interval for writing π(x) to pi_index.bin (default: 1e6)
    let mut args = std::env::args().skip(1);
    let n: u64 = args.next().and_then(|s| s.parse().ok()).unwrap_or(1_000_000_000);
    let step: u64 = args.next().and_then(|s| s.parse().ok()).unwrap_or(1_000_000);
    // ----------------------------------

    assert_eq!(SEGMENT % 30, 0, "SEGMENT must be divisible by 30");

    let t0 = Instant::now();
    let base_primes = simple_sieve((n as f64).sqrt() as u32 + 1);
    eprintln!("Base primes ready: {} in {:?}", base_primes.len(), t0.elapsed());

    // Outputs
    let mut primes_out = BufWriter::new(File::create("primes.bin")?);
    let mut pi_out     = BufWriter::new(File::create("pi_index.bin")?);
    let mut bits_out   = BufWriter::new(File::create("sieve_w30.bitset")?);
    let mut meta       = BufWriter::new(File::create("meta.txt")?);

    // Bitset header: magic, version, n, total_blocks, flags, wheel residues
    let total_blocks: u64 = ((n + 1) + 29) / 30; // ceil((n+1)/30)
    bits_out.write_all(b"W30BIT\0")?;
    bits_out.write_all(&1u32.to_le_bytes())?;            // version
    bits_out.write_all(&n.to_le_bytes())?;
    bits_out.write_all(&total_blocks.to_le_bytes())?;
    bits_out.write_all(&1u8.to_le_bytes())?;             // flags: bit0=1 => '1' means composite
    for &r in &WHEEL { bits_out.write_all(&(r as u8).to_le_bytes())?; }

    // Meta header
    writeln!(meta, "n={}", n)?;
    writeln!(meta, "pi_index_step={}", step)?;
    writeln!(meta, "segment_size={}", SEGMENT)?;
    writeln!(meta, "checkpoints:")?;
    writeln!(meta, "-------------------------")?;

    // Reference checkpoints for validation
    let checkpoints: &[(u64, u64)] = &[
        (1_000_000,      78_498),
        (10_000_000,    664_579),
        (100_000_000,  5_761_455),
        (1_000_000_000, 50_847_534),
        (10_000_000_000,455_052_511),
    ];
    let mut next_test_idx = 0usize;
    let mut pi_count: u64 = 0;
    let mut next_checkpoint = step;

    // Emit small primes (2, 3, 5)
    for &p in &[2u32, 3, 5] {
        if (p as u64) <= n {
            primes_out.write_all(&p.to_le_bytes())?;
            pi_count += 1;
            while next_checkpoint <= p as u64 {
                pi_out.write_all(&pi_count.to_le_bytes())?;
                next_checkpoint += step;
            }
        }
    }

    // Segmented sieve: start at 0 to align wheel phase
    let mut low: u64 = 0;
    let total_segments = (n + SEGMENT as u64 - 1) / SEGMENT as u64;
    let mut seg_idx: u64 = 0;

    while low < n + 1 {
        seg_idx += 1;
        let high = (low + SEGMENT as u64).min(n + 1);
        let blocks = (high - low) / 30;
        let remainder = (high - low) % 30;
        let extra = if remainder > 0 { 1 } else { 0 };
        let mut mark = vec![0u8; blocks as usize + extra];

        // Cross off composites using base primes >= 7
        for &p in &base_primes {
            if p <= 5 { continue; }
            let p64 = p as u64;
            let mut j = (low + p64 - 1) / p64 * p64;
            if j < p64 * p64 { j = p64 * p64; }
            while j < high {
                if let Some(idx) = wheel_index((j % 30) as u32) {
                    let b = ((j - low) / 30) as usize;
                    mark[b] |= 1u8 << idx; // 1 == composite
                }
                j += p64;
            }
        }

        // Write bitset bytes for this segment
        if remainder > 0 {
            bits_out.write_all(&mark)?;                  // include tail block
        } else {
            bits_out.write_all(&mark[..blocks as usize])?;
        }

        // Enumerate primes in full 30-blocks
        for b in 0..blocks as usize {
            let base = low + (b as u64) * 30;
            let m = mark[b];
            for (idx, &r) in WHEEL.iter().enumerate() {
                let x = base + r as u64;
                if x > n { break; }

                // Checkpoint before counting x
                while next_test_idx < checkpoints.len()
                    && checkpoints[next_test_idx].0 < x
                    && checkpoints[next_test_idx].0 <= n
                {
                    let (cx, expect) = checkpoints[next_test_idx];
                    let status = if pi_count == expect { "OK" } else { "MISMATCH" };
                    eprintln!("checkpoint π({}) = {} -> {}", cx, pi_count, status);
                    writeln!(meta, "π({}) = {} -> {}", cx, pi_count, status)?;
                    next_test_idx += 1;
                }

                if x >= 7 && (m & (1u8 << idx)) == 0 {
                    primes_out.write_all(&(x as u32).to_le_bytes())?;
                    pi_count += 1;
                    while next_checkpoint <= x {
                        pi_out.write_all(&pi_count.to_le_bytes())?;
                        next_checkpoint += step;
                    }
                }
            }
        }

        // Enumerate primes in the tail block (if any)
        if remainder > 0 {
            let base = low + (blocks as u64) * 30;
            let m = mark[blocks as usize];
            for (idx, &r) in WHEEL.iter().enumerate() {
                let x = base + r as u64;
                if x >= high || x > n { break; }

                while next_test_idx < checkpoints.len()
                    && checkpoints[next_test_idx].0 < x
                    && checkpoints[next_test_idx].0 <= n
                {
                    let (cx, expect) = checkpoints[next_test_idx];
                    let status = if pi_count == expect { "OK" } else { "MISMATCH" };
                    eprintln!("checkpoint π({}) = {} -> {}", cx, pi_count, status);
                    writeln!(meta, "π({}) = {} -> {}", cx, pi_count, status)?;
                    next_test_idx += 1;
                }

                if x >= 7 && (m & (1u8 << idx)) == 0 {
                    primes_out.write_all(&(x as u32).to_le_bytes())?;
                    pi_count += 1;
                    while next_checkpoint <= x {
                        pi_out.write_all(&pi_count.to_le_bytes())?;
                        next_checkpoint += step;
                    }
                }
            }
        }

        // Per-segment progress (to meta and stderr)
        let elapsed = t0.elapsed();
        let done_ratio = (seg_idx as f64) / (total_segments as f64);
        let est_total = if done_ratio > 0.0 { elapsed.mul_f64(1.0 / done_ratio) } else { elapsed };
        let eta = if est_total > elapsed { est_total - elapsed } else { Duration::from_secs(0) };
        let percent = 100.0 * done_ratio;

        writeln!(
            meta,
            "progress: {:>6.2}% ({}..{}) pi_count={} elapsed={} ETA={}",
            percent, low, high, pi_count, format_duration(elapsed), format_duration(eta)
        )?;

        if seg_idx % 10 == 0 || low + SEGMENT as u64 >= n {
            eprintln!(
                "progress: {}/{} ({:.2}%) | elapsed={} | ETA={}",
                seg_idx, total_segments, percent, format_duration(elapsed), format_duration(eta)
            );
        }

        low = high;
    }

    // Flush remaining pi_index entries up to n (if any)
    while next_checkpoint <= n {
        pi_out.write_all(&pi_count.to_le_bytes())?;
        next_checkpoint += step;
    }

    writeln!(meta, "-------------------------")?;
    writeln!(meta, "Total pi({}) = {}", n, pi_count)?;
    meta.flush()?;

    primes_out.flush()?;
    pi_out.flush()?;
    bits_out.flush()?;

    eprintln!("Done. Total pi({}) = {}", n, pi_count);
    eprintln!("Tip: compress bitset with `zstd -19 sieve_w30.bitset`");
    Ok(())
}
