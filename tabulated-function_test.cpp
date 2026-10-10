#include "tabulated-function.hpp"

#include <cmath>
#include <filesystem>
#include <format>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace tabulatedfunction;

constexpr double float64EqualityThreshold = 1e-9;

bool almostEqual(double a, double b) {
    if (std::isnan(a) && std::isnan(b)) {
        return true;
    }
    return std::abs(a - b) <= float64EqualityThreshold;
}

#define TEST_CHECK(cond) \
    do { \
        if (!(cond)) { \
            throw std::runtime_error(std::format("Condition failed: ({}) at {}:{}", #cond, __FILE__, __LINE__)); \
        } \
    } while (false)

void test_new() {
    TabulatedFunction f;
    TEST_CHECK(f.Order == 1);
    TEST_CHECK(f.trapolation == Trapolation::Linear);
    TEST_CHECK(f.GetNdots() == 0);
    TEST_CHECK(!f.changed);
}

void test_add_point_and_f() {
    TabulatedFunction f;
    f.SetOrder(1);

    TEST_CHECK(std::isnan(f.F(0.0)));

    f.AddPoint(1, 10, 0);
    f.AddPoint(3, 30, 0);
    f.AddPoint(2, 20, 0);

    TEST_CHECK(f.GetNdots() == 3);
    TEST_CHECK(f.X[0] == 1 && f.X[1] == 2 && f.X[2] == 3);

    struct TestCase {
        double x;
        double expected;
    };
    const std::vector<TestCase> testCases = {
        {1.0, 10.0},
        {2.0, 20.0},
        {3.0, 30.0},
        {1.5, 15.0},
        {2.5, 25.0},
        {0.0, 10.0},
        {4.0, 30.0},
    };

    for (const auto& tc : testCases) {
        TEST_CHECK(almostEqual(f.F(tc.x), tc.expected));
    }

    const double oldY = f.AddPoint(2, 22, 1);
    TEST_CHECK(almostEqual(oldY, 20.0));
    TEST_CHECK(almostEqual(f.F(2), 22.0));
}

void test_getters() {
    TabulatedFunction f;
    f.AddPoint(0, 5, 0);
    f.AddPoint(10, -5, 0);
    f.AddPoint(5, 15, 0);

    (void)f.F(1.0);

    TEST_CHECK(almostEqual(f.GetXmin(), 0.0));
    TEST_CHECK(almostEqual(f.GetXmax(), 10.0));
    TEST_CHECK(almostEqual(f.GetYmin(), -5.0));
    TEST_CHECK(almostEqual(f.GetYmax(), 15.0));
    TEST_CHECK(f.GetNdots() == 3);
    TEST_CHECK(almostEqual(f.GetStep(), 5.0));
}

void test_orders() {
    TabulatedFunction f;
    f.AddPoint(0, 0, 0);
    f.AddPoint(1, 1, 0);
    f.AddPoint(2, 0, 0);

    f.SetOrder(0);
    TEST_CHECK(almostEqual(f.F(0.5), 0.0));
    TEST_CHECK(almostEqual(f.F(1.5), 1.0));
    TEST_CHECK(almostEqual(f.F(-1.0), 0.0));
    TEST_CHECK(almostEqual(f.F(3.0), 0.0));

    f.SetOrder(1);
    TEST_CHECK(almostEqual(f.F(0.5), 0.5));
    TEST_CHECK(almostEqual(f.F(1.5), 0.5));
}

void test_trapolation_linear() {
    TabulatedFunction f;
    f.SetTrapolation(Trapolation::Linear);
    f.AddPoint(0, 0, 0);
    f.AddPoint(2, 2, 0);
    TEST_CHECK(almostEqual(f.F(1), 1.0));
}

void test_trapolate_direct() {
    TabulatedFunction f;
    TEST_CHECK(std::isnan(f.Trapolate(10, Trapolation::Linear)));

    f.AddPoint(0, 10, 0);
    f.AddPoint(5, 20, 0);
    f.AddPoint(10, 30, 0);

    TEST_CHECK(almostEqual(f.Trapolate(0, Trapolation::Shift), 15.0));
    TEST_CHECK(almostEqual(f.Trapolate(10, Trapolation::Shift), 25.0));
    TEST_CHECK(almostEqual(f.Trapolate(5, Trapolation::Shift), 20.0));
}

void test_interpolate_edge_cases_and_panic() {
    TabulatedFunction f;
    f.X = {1.0, 1.0};
    f.Y = {10.0, 20.0};

    TEST_CHECK(almostEqual(f._interpolate(1.0, 0, 1, Trapolation::Linear), 10.0));
    TEST_CHECK(almostEqual(f._interpolate(1.0, 0, 1, Trapolation::Cosine), 10.0));

    TabulatedFunction f2;
    f2.AddPoint(0, 0, 0);
    f2.AddPoint(10, 10, 0);
    (void)f2.Trapolate(5, Trapolation::MinMax);

    const double yMinMax = f2._interpolate(5, -1, 5, Trapolation::MinMax);
    TEST_CHECK(!std::isnan(yMinMax));

    bool threw = false;
    try {
        (void)f2._interpolate(5, 0, 1, static_cast<Trapolation>(999));
    } catch (const std::runtime_error&) {
        threw = true;
    }
    TEST_CHECK(threw);
}

void test_calculus_border_cases() {
    TabulatedFunction f0;
    f0.Derivative();
    TEST_CHECK(f0.GetNdots() == 0);

    TabulatedFunction f1;
    f1.AddPoint(1, 50, 0);
    f1.Derivative();
    TEST_CHECK(f1.Y[0] == 0.0);

    TabulatedFunction f2;
    f2.AddPoint(0, 10, 0);
    f2.AddPoint(2, 20, 0);
    f2.Derivative();
    TEST_CHECK(almostEqual(f2.Y[0], 5.0) && almostEqual(f2.Y[1], 5.0));

    f0.Integral();
    TEST_CHECK(f0.GetNdots() == 0);
    f1.Integral();
    TEST_CHECK(f1.Y[0] == 0.0);

    TabulatedFunction f2Int;
    f2Int.AddPoint(0, 4, 0);
    f2Int.AddPoint(3, 4, 0);
    f2Int.Integral();
    TEST_CHECK(almostEqual(f2Int.Y[0], 0.0) && almostEqual(f2Int.Y[1], 12.0));

    TEST_CHECK(almostEqual(f0.Integrate(), 0.0));
    TEST_CHECK(almostEqual(f1.Integrate(), 0.0));

    TabulatedFunction fTrap;
    fTrap.AddPoint(0, 2, 0);
    fTrap.AddPoint(4, 6, 0);
    TEST_CHECK(almostEqual(fTrap.Integrate(), 16.0));

    TabulatedFunction fNonUniform;
    fNonUniform.AddPoint(0, 0, 0);
    fNonUniform.AddPoint(1, 1, 0);
    fNonUniform.AddPoint(3, 9, 0);
    fNonUniform.AddPoint(4, 16, 0);
    const double val = fNonUniform.Integrate();
    TEST_CHECK(!std::isnan(val) && val > 0.0);
}

void test_normalise_border_cases() {
    TabulatedFunction fZeros;
    fZeros.AddPoint(0, 0, 0);
    fZeros.AddPoint(1, 0, 0);
    fZeros.Normalise();
    TEST_CHECK(fZeros.Y[0] == 0.0 && fZeros.Y[1] == 0.0);

    TabulatedFunction fNeg;
    fNeg.AddPoint(0, -2, 0);
    fNeg.AddPoint(1, -10, 0);
    fNeg.Normalise();
    TEST_CHECK(almostEqual(fNeg.Y[0], -0.2) && almostEqual(fNeg.Y[1], -1.0));
}

void test_multiply_border_cases() {
    TabulatedFunction empty;
    TabulatedFunction populated;
    populated.AddPoint(0, 5, 0);
    populated.AddPoint(1, 10, 0);

    empty.Multiply(populated);
    TEST_CHECK(empty.GetNdots() == 0);

    populated.Multiply(TabulatedFunction{});
    TEST_CHECK(populated.GetNdots() == 0);
}

void test_can_insert_point_and_expand_border_cases() {
    TabulatedFunction f;
    TEST_CHECK(f.canInsertPoint(5.0));

    f.AddPoint(0, 0, 0);
    f.AddPoint(10, 10, 0);
    (void)f.GetStep();

    TEST_CHECK(!f.canInsertPoint(0.0) && !f.canInsertPoint(10.0));
    TEST_CHECK(!f.canInsertPoint(2.0));

    const size_t dotsBefore = f.GetNdots();
    f.Expand(0);
    TEST_CHECK(f.GetNdots() == dotsBefore);
}

void test_draw_ps_border_cases() {
    TabulatedFunction fSingle;
    fSingle.AddPoint(5, 5, 0);
    const auto tmpPath = std::filesystem::temp_directory_path() / "test_single_ps.ps";
    TEST_CHECK(fSingle.DrawPS(tmpPath.string()));
    std::filesystem::remove(tmpPath);

    TEST_CHECK(!fSingle.DrawPS("/invalid_nonexistent_directory/file.ps"));
}

void test_dump_and_string() {
    TabulatedFunction f;
    const std::string s = f.to_string();
    TEST_CHECK(!s.empty());

    f.AddPoint(0, 0, 1);
    f.AddPoint(1, 1, 2);
    const Dump d = f.ToDump();
    TabulatedFunction f2;
    f2.FromDump(d);

    TEST_CHECK(f.GetNdots() == f2.GetNdots());
    for (size_t i = 0; i < f.GetNdots(); ++i) {
        TEST_CHECK(almostEqual(f.X[i], f2.X[i]));
        TEST_CHECK(almostEqual(f.Y[i], f2.Y[i]));
        TEST_CHECK(f.epoch[i] == f2.epoch[i]);
    }
}

void test_trapolation_opposite() {
    TabulatedFunction f;
    f.SetTrapolation(Trapolation::Opposite);
    f.AddPoint(0, 10, 0);
    f.AddPoint(10, 20, 0);

    TEST_CHECK(almostEqual(f.F(5), 20.0));
    TEST_CHECK(almostEqual(f.F(6), 10.0));
    TEST_CHECK(almostEqual(f.F(4), 20.0));
    (void)f.GetYmax();
    TEST_CHECK(almostEqual(f.F(15), 10.0));
}

void test_trapolation_nearest() {
    TabulatedFunction f;
    f.SetTrapolation(Trapolation::Nearest);
    f.AddPoint(0, 10, 0);
    f.AddPoint(10, 20, 0);

    TEST_CHECK(almostEqual(f.F(2), 10.0));
    TEST_CHECK(almostEqual(f.F(8), 20.0));
    TEST_CHECK(almostEqual(f.F(5), 20.0));
    TEST_CHECK(almostEqual(f.F(-5), 10.0));
    TEST_CHECK(almostEqual(f.F(15), 20.0));
}

void test_trapolation_cosine() {
    TabulatedFunction f;
    f.SetTrapolation(Trapolation::Cosine);
    f.AddPoint(0, 0, 0);
    f.AddPoint(10, 100, 0);
    TEST_CHECK(almostEqual(f.F(5), 50.0));
}

void test_trapolation_shift() {
    TabulatedFunction f;
    f.SetTrapolation(Trapolation::Shift);
    f.AddPoint(0, 10, 0);
    f.AddPoint(2, 20, 0);
    TEST_CHECK(almostEqual(f.F(1), 15.0));
    TEST_CHECK(almostEqual(f.F(-1), 10.0));
    TEST_CHECK(almostEqual(f.F(3), 20.0));
}

void test_trapolation_min_max() {
    TabulatedFunction f;
    f.SetTrapolation(Trapolation::MinMax);
    f.AddPoint(0, 0, 0);
    f.AddPoint(10, 10, 0);
    (void)f.F(5);
    const double y = f.F(5);
    TEST_CHECK(!std::isnan(y) && y >= 0.0 && y <= 10.0);
}

void test_trapolation_cosine_extrapolation() {
    TabulatedFunction f;
    f.SetTrapolation(Trapolation::Cosine);
    f.AddPoint(0, 10, 0);
    f.AddPoint(10, 20, 0);
    TEST_CHECK(almostEqual(f.F(-5), 10.0));
    TEST_CHECK(almostEqual(f.F(15), 20.0));
}

void test_load_constant() {
    TabulatedFunction f;
    f.LoadConstant(100, -5, 5);
    TEST_CHECK(f.GetNdots() == 1);
    TEST_CHECK(almostEqual(f.F(0), 100.0));
    TEST_CHECK(almostEqual(f.F(-10), 100.0));
}

void test_clear() {
    TabulatedFunction f;
    f.AddPoint(1, 1, 0);
    f.Clear();
    TEST_CHECK(f.GetNdots() == 0);
    TEST_CHECK(std::isnan(f.F(1)));
}

void test_derivative_and_integral() {
    TabulatedFunction f;
    f.AddPoint(-2, 4, 0);
    f.AddPoint(-1, 1, 0);
    f.AddPoint(0, 0, 0);
    f.AddPoint(1, 1, 0);
    f.AddPoint(2, 4, 0);
    f.SetOrder(1);

    const double integral = f.Integrate();
    TEST_CHECK(almostEqual(integral, 16.0 / 3.0));

    f.Derivative();
    TEST_CHECK(almostEqual(f.F(0), 0.0));
    TEST_CHECK(almostEqual(f.F(1), 2.0));

    f.Integral();
    const double y_neg1 = f.F(-1);
    const double y0 = f.F(0);
    const double y1 = f.F(1);
    TEST_CHECK(std::abs((y1 - y0) - (y0 - y_neg1) - 2.0) <= 0.5);
}

void test_from_dump_unsorted_and_deduplication() {
    Dump d{
        .order = 1,
        .trapolation = Trapolation::Linear,
        .x = {3.0, 1.0, 2.0, 2.0},
        .y = {30.0, 10.0, 20.0, 40.0},
        .epoch = {0, 0, 1, 2},
    };

    TabulatedFunction f;
    f.FromDump(d);

    TEST_CHECK(f.GetNdots() == 3);
    TEST_CHECK(almostEqual(f.X[0], 1.0) && almostEqual(f.Y[0], 10.0));
    TEST_CHECK(almostEqual(f.X[1], 2.0) && almostEqual(f.Y[1], 30.0) && f.epoch[1] == 2);
    TEST_CHECK(almostEqual(f.X[2], 3.0) && almostEqual(f.Y[2], 30.0));
}

void test_more_points() {
    TabulatedFunction f;
    f.SetOrder(1);
    f.AddPoint(0, 0, 0);
    f.AddPoint(2, 4, 0);
    f.MorePoints();

    TEST_CHECK(f.GetNdots() == 3);
    TEST_CHECK(almostEqual(f.X[1], 1.0));
    TEST_CHECK(almostEqual(f.Y[1], 2.0));
}

void test_epoch() {
    TabulatedFunction f;
    f.AddPoint(0, 0, 0);
    f.AddPoint(1, 1, 1);
    f.AddPoint(2, 2, 2);
    f.AddPoint(3, 3, 3);

    f.Epoch(2);
    TEST_CHECK(f.GetNdots() == 2);
    TEST_CHECK(f.X[0] == 2 && f.X[1] == 3);

    f.Epoch(10);
    TEST_CHECK(f.GetNdots() == 0);
}

void test_smooth() {
    TabulatedFunction f;
    f.AddPoint(0, 0, 0);
    f.AddPoint(1, 10, 0);
    f.AddPoint(2, 4, 0);

    f.Smooth();
    TEST_CHECK(almostEqual(f.Y[1], 14.0 / 3.0));
}

void test_multiply_and_merge() {
    TabulatedFunction f1;
    f1.AddPoint(0, 2, 0);
    f1.AddPoint(5, 2, 0);

    TabulatedFunction f2;
    f2.AddPoint(0, 0, 0);
    f2.AddPoint(5, 5, 0);

    f1.Multiply(f2);
    TEST_CHECK(almostEqual(f1.F(2.5), 5.0));
    TEST_CHECK(almostEqual(f1.F(5), 10.0));

    TabulatedFunction m1;
    m1.AddPoint(0, 0, 0);
    m1.AddPoint(2, 20, 0);
    TabulatedFunction m2;
    m2.AddPoint(1, 10, 0);
    m2.AddPoint(3, 30, 0);
    m2.AddPoint(2, 22, 1);

    m1.Merge(m2);
    TEST_CHECK(m1.GetNdots() == 4);
    TEST_CHECK(almostEqual(m1.X[2], 2.0) && almostEqual(m1.Y[2], 22.0) && m1.epoch[2] == 1);
}

int main() {
    const std::vector<std::pair<std::string, void (*)()>> tests = {
        {"test_new", test_new},
        {"test_add_point_and_f", test_add_point_and_f},
        {"test_getters", test_getters},
        {"test_orders", test_orders},
        {"test_trapolation_linear", test_trapolation_linear},
        {"test_trapolate_direct", test_trapolate_direct},
        {"test_interpolate_edge_cases_and_panic", test_interpolate_edge_cases_and_panic},
        {"test_calculus_border_cases", test_calculus_border_cases},
        {"test_normalise_border_cases", test_normalise_border_cases},
        {"test_multiply_border_cases", test_multiply_border_cases},
        {"test_can_insert_point_and_expand_border_cases", test_can_insert_point_and_expand_border_cases},
        {"test_draw_ps_border_cases", test_draw_ps_border_cases},
        {"test_dump_and_string", test_dump_and_string},
        {"test_trapolation_opposite", test_trapolation_opposite},
        {"test_trapolation_nearest", test_trapolation_nearest},
        {"test_trapolation_cosine", test_trapolation_cosine},
        {"test_trapolation_shift", test_trapolation_shift},
        {"test_trapolation_min_max", test_trapolation_min_max},
        {"test_trapolation_cosine_extrapolation", test_trapolation_cosine_extrapolation},
        {"test_load_constant", test_load_constant},
        {"test_clear", test_clear},
        {"test_derivative_and_integral", test_derivative_and_integral},
        {"test_from_dump_unsorted_and_deduplication", test_from_dump_unsorted_and_deduplication},
        {"test_more_points", test_more_points},
        {"test_epoch", test_epoch},
        {"test_smooth", test_smooth},
        {"test_multiply_and_merge", test_multiply_and_merge},
    };

    int failed = 0;
    for (const auto& [name, func] : tests) {
        try {
            func();
            std::cout << "[       OK ] " << name << "\n";
        } catch (const std::exception& e) {
            std::cerr << "[  FAILED  ] " << name << ": " << e.what() << "\n";
            failed++;
        }
    }

    if (failed > 0) {
        std::cerr << failed << " test(s) failed.\n";
        return 1;
    }
    std::cout << "All " << tests.size() << " tests passed.\n";
    return 0;
}