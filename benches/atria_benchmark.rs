use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId};
use std::hint::black_box;

// The crate exports are under the library name from Cargo.toml
extern crate ATriaPlugin as atria;
use atria::{ATriaPlugin, cpu_floyd};
use pluma_plugin_trait::PluMAPlugin;

/// Benchmark the core Floyd-Warshall algorithm with different matrix sizes
fn bench_cpu_floyd(c: &mut Criterion) {
    let mut group = c.benchmark_group("cpu_floyd");

    // Test various matrix sizes (2N x 2N where N is the number of bacteria)
    for size in [10, 20, 50, 100, 252].iter() {
        let n = *size;
        let mut graph = vec![0.0f32; n * n];

        // Initialize with some test data - diagonal = 1, some edges with weights
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    graph[i * n + j] = 1.0;
                } else if (i + j) % 3 == 0 {
                    graph[i * n + j] = 0.5;
                } else if (i + j) % 5 == 0 {
                    graph[i * n + j] = -0.3;
                }
            }
        }

        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, &n| {
            b.iter(|| {
                let mut g = graph.clone();
                cpu_floyd(black_box(&mut g), black_box(n));
            });
        });
    }

    group.finish();
}

/// Benchmark the full ATria run with the test CSV file
fn bench_atria_run(c: &mut Criterion) {
    let mut group = c.benchmark_group("atria_full_run");

    // Benchmark with the test file (126 bacteria)
    group.sample_size(10); // Reduce sample size for long-running benchmark

    group.bench_function("corrP.never.csv (126 bacteria)", |b| {
        b.iter_with_setup(
            || {
                let mut plugin = ATriaPlugin::default();
                plugin.input("./tests/corrP.never.csv".to_string()).unwrap();
                plugin
            },
            |mut plugin| {
                plugin.run().unwrap();
            },
        );
    });

    group.finish();
}

/// Benchmark input parsing
fn bench_input_parsing(c: &mut Criterion) {
    c.bench_function("input_parsing_126_bacteria", |b| {
        b.iter(|| {
            let mut plugin = ATriaPlugin::default();
            plugin.input(black_box("./tests/corrP.never.csv".to_string())).unwrap();
        });
    });
}

/// Benchmark the pay calculation loop (isolated)
fn bench_pay_calculation(c: &mut Criterion) {
    let mut group = c.benchmark_group("pay_calculation");

    for gsize in [10, 50, 126].iter() {
        let n = gsize * 2;
        let graph = vec![0.5f32; n * n];

        group.bench_with_input(BenchmarkId::from_parameter(gsize), gsize, |b, &gsize| {
            b.iter(|| {
                let mut h_pay = vec![0.0f32; gsize];
                for i in 0..gsize {
                    let mut pay = 0.0f32;
                    for j in 0..n {
                        pay += graph[(i * 2) * n + j];
                    }
                    pay -= 1.0;
                    h_pay[i] = pay;
                }
                black_box(h_pay)
            });
        });
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_cpu_floyd,
    bench_atria_run,
    bench_input_parsing,
    bench_pay_calculation
);
criterion_main!(benches);
