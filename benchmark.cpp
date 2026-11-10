//go:build ignore

#include "tabulated-function.hpp"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <format>
#include <functional>
#include <iostream>
#include <new>
#include <string_view>

struct MemoryStats {
    std::size_t alloc_count{0};
    std::size_t bytes_allocated{0};
};

inline thread_local bool g_track_memory{false};
inline thread_local MemoryStats g_mem_stats{};

void* operator new(std::size_t size) {
    if (g_track_memory) {
        g_mem_stats.alloc_count++;
        g_mem_stats.bytes_allocated += size;
    }
    void* ptr = std::malloc(size);
    if (!ptr) throw std::bad_alloc();
    return ptr;
}

void operator delete(void* ptr) noexcept {
    std::free(ptr);
}

void operator delete(void* ptr, std::size_t) noexcept {
    std::free(ptr);
}

void* operator new[](std::size_t size) {
    if (g_track_memory) {
        g_mem_stats.alloc_count++;
        g_mem_stats.bytes_allocated += size;
    }
    void* ptr = std::malloc(size);
    if (!ptr) throw std::bad_alloc();
    return ptr;
}

void operator delete[](void* ptr) noexcept {
    std::free(ptr);
}

void operator delete[](void* ptr, std::size_t) noexcept {
    std::free(ptr);
}

namespace {

template <typename T>
inline void DoNotOptimize(T const& val) {
#if defined(__GNUC__) || defined(__clang__)
    asm volatile("" : : "g"(val) : "memory");
#else
    volatile char const dummy = *reinterpret_cast<char const*>(&val);
    (void)dummy;
#endif
}

void RunBenchmark(std::string_view name, std::function<void(std::size_t)> fn, std::size_t warmups = 5, std::chrono::milliseconds target_duration = std::chrono::milliseconds(1000)) {
    // Warmup
    for (std::size_t i = 0; i < warmups; ++i) {
        fn(100);
    }

    // Iteration estimate
    std::size_t iters = 1000;
    auto start = std::chrono::steady_clock::now();
    fn(iters);
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start);

    if (elapsed.count() > 0) {
        double ns_per_it = static_cast<double>(elapsed.count()) / static_cast<double>(iters);
        double desired_ns = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(target_duration).count());
        iters = std::max<std::size_t>(1000, static_cast<std::size_t>(desired_ns / std::max(ns_per_it, 1.0)));
    }

    g_mem_stats = {};
    g_track_memory = true;
    start = std::chrono::steady_clock::now();
    fn(iters);
    elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start);
    g_track_memory = false;

    double ns_per_op = static_cast<double>(elapsed.count()) / static_cast<double>(iters);
    double bytes_per_op = static_cast<double>(g_mem_stats.bytes_allocated) / static_cast<double>(iters);
    double allocs_per_op = static_cast<double>(g_mem_stats.alloc_count) / static_cast<double>(iters);

    std::cout << std::format("{:<32} {:>10} {:>12.2f} ns/op {:>10.0f} B/op {:>10.0f} allocs/op\n",
                             name, iters, ns_per_op, bytes_per_op, allocs_per_op);
}

} // namespace

int main() {
    std::cout << std::format("{:<32} {:>10} {:>18} {:>15} {:>17}\n", "Benchmark", "Iterations", "Time", "Memory", "Allocs");
    std::cout << std::string(96, '-') << "\n";

    // 1. Sequential point insertion
    RunBenchmark("BenchmarkAddPointSequential", [](std::size_t iters) {
        for (std::size_t i = 0; i < iters; ++i) {
            tabulatedfunction::TabulatedFunction f;
            for (int j = 0; j < 1000; ++j) {
                f.AddPoint(static_cast<double>(j), static_cast<double>(j * 10));
            }
            DoNotOptimize(f);
        }
    });

    // 2. Random/binary-search point insertion
    {
        tabulatedfunction::TabulatedFunction base;
        for (int j = 0; j < 1000; ++j) {
            base.AddPoint(static_cast<double>(j) * 2.0, static_cast<double>(j));
        }

        RunBenchmark("BenchmarkAddPointRandom", [&base](std::size_t iters) {
            tabulatedfunction::TabulatedFunction f = base;
            for (std::size_t i = 0; i < iters; ++i) {
                double x = static_cast<double>(i % 999) * 2.0 + 1.0;
                f.AddPoint(x, x * 1.5);
            }
            DoNotOptimize(f);
        });
    }

    // 3. F(x) Evaluation (Binary Search + Linear Interpolation)
    {
        tabulatedfunction::TabulatedFunction f;
        for (int j = 0; j < 1000; ++j) {
            f.AddPoint(static_cast<double>(j), std::sin(static_cast<double>(j) * 0.1));
        }

        RunBenchmark("BenchmarkFEvaluation", [&f](std::size_t iters) {
            double sum = 0.0;
            for (std::size_t i = 0; i < iters; ++i) {
                double x = static_cast<double>(i % 1000) + 0.5;
                sum += f.F(x);
            }
            DoNotOptimize(sum);
        });
    }

    // 4. Trapolate with Cosine mode
    {
        tabulatedfunction::TabulatedFunction f;
        for (int j = 0; j < 1000; ++j) {
            f.AddPoint(static_cast<double>(j), std::sin(static_cast<double>(j) * 0.1));
        }

        RunBenchmark("BenchmarkTrapolateCosine", [&f](std::size_t iters) {
            double sum = 0.0;
            for (std::size_t i = 0; i < iters; ++i) {
                double x = static_cast<double>(i % 1000) + 0.5;
                sum += f.Trapolate(x, tabulatedfunction::Trapolation::Cosine);
            }
            DoNotOptimize(sum);
        });
    }

    // 5. Definite Integration (Simpson's 3-point rule)
    {
        tabulatedfunction::TabulatedFunction f;
        for (int j = 0; j < 1000; ++j) {
            f.AddPoint(static_cast<double>(j), static_cast<double>(j * j));
        }

        RunBenchmark("BenchmarkIntegrate", [&f](std::size_t iters) {
            double total = 0.0;
            for (std::size_t i = 0; i < iters; ++i) {
                total += f.Integrate();
            }
            DoNotOptimize(total);
        });
    }

    // 6. Numerical Differentiation
    {
        tabulatedfunction::TabulatedFunction f_orig;
        for (int j = 0; j < 1000; ++j) {
            f_orig.AddPoint(static_cast<double>(j), static_cast<double>(j * j));
        }

        RunBenchmark("BenchmarkDerivative", [&f_orig](std::size_t iters) {
            for (std::size_t i = 0; i < iters; ++i) {
                tabulatedfunction::TabulatedFunction f = f_orig;
                f.Derivative();
                DoNotOptimize(f);
            }
        });
    }

    // 7. 3-Point Smoothing
    {
        tabulatedfunction::TabulatedFunction f_orig;
        for (int j = 0; j < 1000; ++j) {
            f_orig.AddPoint(static_cast<double>(j), std::sin(static_cast<double>(j)));
        }

        RunBenchmark("BenchmarkSmooth", [&f_orig](std::size_t iters) {
            for (std::size_t i = 0; i < iters; ++i) {
                tabulatedfunction::TabulatedFunction f = f_orig;
                f.Smooth();
                DoNotOptimize(f);
            }
        });
    }

    // 8. Midpoint Subdivision (MorePoints)
    {
        tabulatedfunction::TabulatedFunction f_orig;
        for (int j = 0; j < 500; ++j) {
            f_orig.AddPoint(static_cast<double>(j), static_cast<double>(j * 2));
        }

        RunBenchmark("BenchmarkMorePoints", [&f_orig](std::size_t iters) {
            for (std::size_t i = 0; i < iters; ++i) {
                tabulatedfunction::TabulatedFunction f = f_orig;
                f.MorePoints();
                DoNotOptimize(f);
            }
        });
    }

    // 9. Pointwise Multiplication
    {
        tabulatedfunction::TabulatedFunction f1_orig, f2_orig;
        for (int j = 0; j < 500; ++j) {
            f1_orig.AddPoint(static_cast<double>(j), static_cast<double>(j));
            f2_orig.AddPoint(static_cast<double>(j) + 0.5, static_cast<double>(j * 2));
        }

        RunBenchmark("BenchmarkMultiply", [&f1_orig, &f2_orig](std::size_t iters) {
            for (std::size_t i = 0; i < iters; ++i) {
                tabulatedfunction::TabulatedFunction f1 = f1_orig;
                tabulatedfunction::TabulatedFunction f2 = f2_orig;
                f1.Multiply(f2);
                DoNotOptimize(f1);
            }
        });
    }

    return 0;
}