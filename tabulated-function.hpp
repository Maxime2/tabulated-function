#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <format>
#include <fstream>
#include <iterator>
#include <limits>
#include <numbers>
#include <numeric>
#include <ranges>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include <nlohmann/json.hpp>

namespace tabulatedfunction {

enum class Trapolation : int {
    Linear   = 0,
    Shift    = 2,
    MinMax   = 3,
    Opposite = 4,
    Nearest  = 5,
    Cosine   = 6,
};

struct Dump {
    int order{1};
    Trapolation trapolation{Trapolation::Linear};
    std::vector<double> x;
    std::vector<double> y;
    std::vector<uint32_t> epoch;
};

inline void to_json(nlohmann::json& j, const Dump& d) {
    j = nlohmann::json{
        {"order", d.order},
        {"trapolation", static_cast<int>(d.trapolation)},
        {"x", d.x},
        {"y", d.y}
    };
    if (!d.epoch.empty()) {
        j["epoch"] = d.epoch;
    }
}

inline void from_json(const nlohmann::json& j, Dump& d) {
    d.order = j.value("order", 1);
    d.trapolation = static_cast<Trapolation>(j.value("trapolation", 0));
    if (j.contains("x") && !j["x"].is_null()) {
        j.at("x").get_to(d.x);
    } else {
        d.x.clear();
    }
    if (j.contains("y") && !j["y"].is_null()) {
        j.at("y").get_to(d.y);
    } else {
        d.y.clear();
    }
    if (j.contains("epoch") && !j["epoch"].is_null()) {
        j.at("epoch").get_to(d.epoch);
    } else {
        d.epoch.clear();
    }
}

class TabulatedFunction {
public:
    mutable double ixmin{0.0};
    mutable double ixmax{0.0};
    mutable double iymin{0.0};
    mutable double iymax{0.0};
    mutable double istep{0.0};
    mutable bool changed{false};

    int Order{1};
    Trapolation trapolation{Trapolation::Linear};
    std::vector<double> X;
    std::vector<double> Y;
    std::vector<uint32_t> epoch;
    std::vector<uint32_t> indices;
    uint32_t nextIndex{1};

    [[nodiscard]] static inline size_t branchless_lower_bound(const double* arr, size_t n, double target) noexcept {
        const double* base = arr;
        while (n > 1) {
            const size_t half = n >> 1;
            base = (base[half] < target) ? (base + half) : base;
            n -= half;
        }
        return static_cast<size_t>(base - arr) + (*base < target ? 1 : 0);
    }

    TabulatedFunction() = default;

    void update_spline() const {
        changed = false;
        const size_t n = X.size();
        if (n == 0) {
            ixmin = 0.0;
            ixmax = 0.0;
            iymin = 0.0;
            iymax = 0.0;
            istep = 0.0;
            return;
        }
        const size_t j = n - 1;
        ixmin = X[0];
        ixmax = X[j];
        const auto [min_it, max_it] = std::ranges::minmax_element(Y);
        iymin = *min_it;
        iymax = *max_it;
        if (j > 0) {
            double min_step = X[1] - X[0];
            for (size_t i = 2; i <= j; ++i) {
                const double diff = X[i] - X[i - 1];
                min_step = std::min(min_step, diff);
            }
            istep = min_step;
        } else {
            istep = 0.0;
        }
    }

    [[nodiscard]] double _interpolate(double xi, int left, int right, Trapolation trap) const {
        switch (trap) {
        case Trapolation::Linear: {
            if (left < 0) left = 0;
            if (left == right || left < 0 || left >= static_cast<int>(X.size())) {
                return Y[right];
            }
            if (Order == 0) {
                return Y[left];
            }
            const double dx = X[right] - X[left];
            if (dx == 0.0) {
                return Y[left];
            }
            const double dy = Y[right] - Y[left];
            return Y[left] + dy * (xi - X[left]) / dx;
        }
        case Trapolation::Opposite: {
            if (left == right || left < 0 || left >= static_cast<int>(X.size())) {
                return GetYmin() + GetYmax() - Y[right];
            }
            const double midpoint = (X[left] + X[right]) * 0.5;
            if (xi <= midpoint) {
                return Y[right];
            }
            return Y[left];
        }
        case Trapolation::Shift: {
            if (left < 0 || left >= static_cast<int>(X.size())) {
                return Y[right];
            }
            return (Y[left] + Y[right]) * 0.5;
        }
        case Trapolation::MinMax: {
            if (left < 0) left = 0;
            if (right >= static_cast<int>(X.size())) right = static_cast<int>(X.size()) - 1;
            if (changed) update_spline();

            std::array<double, 2> v{0.0, 0.0};
            v[0] += std::abs(iymin - Y[left]);
            v[0] += std::abs(iymin - Y[right]);
            v[1] += std::abs(iymax - Y[left]);
            v[1] += std::abs(iymax - Y[right]);

            const double avg = (Y[left] + Y[right]) * 0.5;
            if (v[0] > v[1]) {
                return (iymin + avg) * 0.5;
            }
            return (iymax + avg) * 0.5;
        }
        case Trapolation::Nearest: {
            if (left == right || left < 0 || left >= static_cast<int>(X.size())) {
                return Y[right];
            }
            if (std::abs(xi - X[left]) < std::abs(xi - X[right])) {
                return Y[left];
            }
            return Y[right];
        }
        case Trapolation::Cosine: {
            if (left == right || left < 0 || left >= static_cast<int>(X.size())) {
                return Y[right];
            }
            const double dx = X[right] - X[left];
            if (dx == 0.0) {
                return Y[left];
            }
            const double mu = (xi - X[left]) / dx;
            const double mu2 = (1.0 - std::cos(mu * std::numbers::pi)) / 2.0;
            return Y[left] * (1.0 - mu2) + Y[right] * mu2;
        }
        }
        throw std::runtime_error("unhandled trapolation type");
    }

    [[nodiscard]] double F(double xi) const {
        const size_t l = X.size();
        if (l == 0) [[unlikely]] {
            return std::numeric_limits<double>::quiet_NaN();
        }
        const size_t k = branchless_lower_bound(X.data(), l, xi);
        if (k < l && X[k] == xi) [[likely]] {
            return Y[k];
        }

        const int left = (k == 0) ? 0 : static_cast<int>(k) - 1;
        const int right = (k >= l) ? static_cast<int>(l) - 1 : static_cast<int>(k);

        return _interpolate(xi, left, right, trapolation);
    }

    [[nodiscard]] double Trapolate(double xi, Trapolation customTrapolation) const {
        const size_t l = X.size();
        if (l == 0) [[unlikely]] {
            return std::numeric_limits<double>::quiet_NaN();
        }
        const size_t k = branchless_lower_bound(X.data(), l, xi);
        const bool found = (k < l && X[k] == xi);

        int left = 0;
        int right = 0;
        if (found) {
            right = (k == l - 1) ? static_cast<int>(k) : static_cast<int>(k) + 1;
            left = (k == 0) ? static_cast<int>(k) : static_cast<int>(k) - 1;
        } else {
            right = (k >= l) ? static_cast<int>(l) - 1 : static_cast<int>(k);
            left = (k == 0) ? 0 : static_cast<int>(k) - 1;
        }

        return _interpolate(xi, left, right, customTrapolation);
    }

    void SetOrder(int new_value) noexcept {
        Order = new_value;
        changed = true;
    }

    void SetTrapolation(Trapolation new_value) noexcept {
        trapolation = new_value;
        changed = true;
    }

    double AddPoint(double Xn, double Yn, uint32_t ep = 0) {
        changed = true;
        const size_t n = X.size();

        if (n == 0 || Xn > X[n - 1]) {
            X.push_back(Xn);
            Y.push_back(Yn);
            epoch.push_back(ep);
            indices.push_back(nextIndex++);
            return Yn;
        }

        if (Xn == X[n - 1]) {
            epoch[n - 1] = ep;
            const double old = Y[n - 1];
            Y[n - 1] = Yn;
            indices[n - 1] = nextIndex++;
            return old;
        }

        if (Xn < X[0]) {
            X.insert(X.begin(), Xn);
            Y.insert(Y.begin(), Yn);
            epoch.insert(epoch.begin(), ep);
            indices.insert(indices.begin(), nextIndex++);
            return Yn;
        }

        const size_t idx = branchless_lower_bound(X.data(), n, Xn);
        if (idx < n && X[idx] == Xn) {
            epoch[idx] = ep;
            const double old = Y[idx];
            Y[idx] = Yn;
            indices[idx] = nextIndex++;
            return old;
        }

        X.insert(X.begin() + idx, Xn);
        Y.insert(Y.begin() + idx, Yn);
        epoch.insert(epoch.begin() + idx, ep);
        indices.insert(indices.begin() + idx, nextIndex++);
        return Yn;
    }

    void LoadConstant(double new_Y, double new_xmin, double new_xmax) noexcept {
        ixmin = new_xmin;
        ixmax = new_xmax;
        iymin = new_Y;
        iymax = iymin;
        X = {ixmin};
        Y = {iymin};
        epoch = {0};
        indices = {nextIndex++};
        istep = ixmax - ixmin;
        changed = false;
    }

    void Normalise() noexcept {
        if (changed) update_spline();
        const double ym = std::max(std::abs(iymin), std::abs(iymax));
        if (ym > 0.0) {
            const double inv_ym = 1.0 / ym;
            for (auto& y : Y) {
                y *= inv_ym;
            }
            changed = true;
        }
    }

    std::pair<uint32_t, uint32_t> NormaliseIndices() noexcept {
        if (indices.empty()) {
            nextIndex = 1;
            return {0, 0};
        }
        const auto [min_it, max_it] = std::ranges::minmax_element(indices);
        const uint32_t minIndex = *min_it;
        const uint32_t maxIndex = *max_it;

        for (auto& idx : indices) {
            idx = idx - minIndex + 1;
        }
        nextIndex = maxIndex - minIndex + 2;
        return {1, maxIndex - minIndex + 1};
    }

    void Smooth() noexcept {
        const size_t n = Y.size();
        if (n < 3) return;

        constexpr double inv3 = 1.0 / 3.0;
        double prev = Y[0];
        double curr = Y[1];
        for (size_t i = 1; i < n - 1; ++i) {
            const double next = Y[i + 1];
            Y[i] = (prev + curr + next) * inv3;
            prev = curr;
            curr = next;
        }
        changed = true;
    }

    void Multiply(const TabulatedFunction& by) {
        if (changed) update_spline();
        if (by.changed) by.update_spline();

        if (X.empty()) return;
        if (by.X.empty()) {
            Clear();
            return;
        }

        std::vector<double> sortedX;
        sortedX.reserve(X.size() + by.X.size());
        size_t i = 0, j = 0;
        while (i < X.size() && j < by.X.size()) {
            if (X[i] < by.X[j]) {
                sortedX.push_back(X[i++]);
            } else if (X[i] > by.X[j]) {
                sortedX.push_back(by.X[j++]);
            } else {
                sortedX.push_back(X[i++]);
                j++;
            }
        }
        sortedX.insert(sortedX.end(), X.begin() + i, X.end());
        sortedX.insert(sortedX.end(), by.X.begin() + j, by.X.end());

        const size_t k = sortedX.size();
        std::vector<double> newY(k);
        std::vector<uint32_t> newEpoch(k, 0);
        std::vector<uint32_t> newIndices(k);

        for (size_t idx = 0; idx < k; ++idx) {
            const double x = sortedX[idx];
            newY[idx] = F(x) * by.F(x);
            newIndices[idx] = nextIndex++;
        }

        X = std::move(sortedX);
        Y = std::move(newY);
        epoch = std::move(newEpoch);
        indices = std::move(newIndices);
        changed = true;
    }

    void MultiplyByScalar(double by) noexcept {
        for (auto& y : Y) {
            y *= by;
        }
        changed = true;
    }

    void Assign(const TabulatedFunction& s) {
        ixmin = s.ixmin;
        ixmax = s.ixmax;
        iymin = s.iymin;
        iymax = s.iymax;
        istep = s.istep;
        Order = s.Order;
        trapolation = s.trapolation;
        nextIndex = s.nextIndex;

        X = s.X;
        Y = s.Y;
        epoch = s.epoch;
        indices = s.indices;
        changed = true;
    }

    void Merge(const TabulatedFunction& m) {
        if (m.X.empty()) return;
        if (X.empty()) {
            Assign(m);
            return;
        }

        const size_t capGuess = X.size() + m.X.size();
        std::vector<double> newX;
        std::vector<double> newY;
        std::vector<uint32_t> newEpoch;
        std::vector<uint32_t> newIndices;
        newX.reserve(capGuess);
        newY.reserve(capGuess);
        newEpoch.reserve(capGuess);
        newIndices.reserve(capGuess);

        size_t i = 0, j = 0;
        while (i < X.size() && j < m.X.size()) {
            if (X[i] < m.X[j]) {
                newX.push_back(X[i]);
                newY.push_back(Y[i]);
                newEpoch.push_back(epoch[i]);
                newIndices.push_back(indices[i]);
                i++;
            } else if (X[i] > m.X[j]) {
                newX.push_back(m.X[j]);
                newY.push_back(m.Y[j]);
                newEpoch.push_back(m.epoch[j]);
                newIndices.push_back(m.indices[j]);
                j++;
            } else {
                newX.push_back(m.X[j]);
                newY.push_back(m.Y[j]);
                newEpoch.push_back(m.epoch[j]);
                newIndices.push_back(m.indices[j]);
                i++;
                j++;
            }
        }
        newX.insert(newX.end(), X.begin() + i, X.end());
        newY.insert(newY.end(), Y.begin() + i, Y.end());
        newEpoch.insert(newEpoch.end(), epoch.begin() + i, epoch.end());
        newIndices.insert(newIndices.end(), indices.begin() + i, indices.end());

        newX.insert(newX.end(), m.X.begin() + j, m.X.end());
        newY.insert(newY.end(), m.Y.begin() + j, m.Y.end());
        newEpoch.insert(newEpoch.end(), m.epoch.begin() + j, m.epoch.end());
        newIndices.insert(newIndices.end(), m.indices.begin() + j, m.indices.end());

        X = std::move(newX);
        Y = std::move(newY);
        epoch = std::move(newEpoch);
        indices = std::move(newIndices);
        nextIndex = std::max(nextIndex, m.nextIndex);
        changed = true;
    }

    [[nodiscard]] double Integrate() const noexcept {
        const size_t n = X.size();
        if (n < 2) return 0.0;
        double sum = 0.0;
        size_t i = 0;
        constexpr double inv6 = 1.0 / 6.0;
        for (; i + 2 < n; i += 2) {
            const double h1 = X[i + 1] - X[i];
            const double h2 = X[i + 2] - X[i + 1];
            if (h1 <= 0.0 || h2 <= 0.0) {
                sum += h1 * (Y[i] + Y[i + 1]) * 0.5 + h2 * (Y[i + 1] + Y[i + 2]) * 0.5;
                continue;
            }
            const double inv_h1 = 1.0 / h1;
            const double inv_h2 = 1.0 / h2;
            const double h1_plus_h2 = h1 + h2;
            const double term1 = (2.0 - h2 * inv_h1) * Y[i];
            const double term2 = (h1_plus_h2 * h1_plus_h2 * (inv_h1 * inv_h2)) * Y[i + 1];
            const double term3 = (2.0 - h1 * inv_h2) * Y[i + 2];
            sum += h1_plus_h2 * inv6 * (term1 + term2 + term3);
        }
        if (i + 1 < n) {
            const double h = X[i + 1] - X[i];
            sum += h * (Y[i] + Y[i + 1]) * 0.5;
        }
        return sum;
    }

    void Clear() noexcept {
        X.clear();
        Y.clear();
        epoch.clear();
        indices.clear();
        ixmin = 0.0;
        ixmax = 0.0;
        iymin = 0.0;
        iymax = 0.0;
        istep = 0.0;
        changed = false;
        nextIndex = 1;
    }

    void MorePoints() {
        if (changed) update_spline();
        const size_t numPoints = X.size();
        if (numPoints <= 1) return;

        const size_t newSize = numPoints + (numPoints - 1);
        std::vector<double> newX;
        std::vector<double> newY;
        std::vector<uint32_t> newEpoch;
        std::vector<uint32_t> newIndices;
        newX.reserve(newSize);
        newY.reserve(newSize);
        newEpoch.reserve(newSize);
        newIndices.reserve(newSize);

        newX.push_back(X[0]);
        newY.push_back(Y[0]);
        newEpoch.push_back(epoch[0]);
        newIndices.push_back(indices[0]);

        for (size_t i = 0; i < numPoints - 1; ++i) {
            const double x1 = X[i], x2 = X[i + 1];
            const double midX = (x1 + x2) * 0.5;

            newX.push_back(midX);
            newY.push_back(_interpolate(midX, static_cast<int>(i), static_cast<int>(i + 1), trapolation));
            newEpoch.push_back(epoch[i + 1]);
            newIndices.push_back(nextIndex++);

            newX.push_back(x2);
            newY.push_back(Y[i + 1]);
            newEpoch.push_back(epoch[i + 1]);
            newIndices.push_back(indices[i + 1]);
        }

        X = std::move(newX);
        Y = std::move(newY);
        epoch = std::move(newEpoch);
        indices = std::move(newIndices);
        changed = true;
    }

    void Derivative() {
        const size_t n = X.size();
        if (n == 0) return;
        if (n == 1) {
            Y[0] = 0.0;
            changed = true;
            return;
        }
        std::vector<double> newY(n);
        if (n == 2) {
            const double slope = (Y[1] - Y[0]) / (X[1] - X[0]);
            newY[0] = slope;
            newY[1] = slope;
        } else {
            const double h0 = X[1] - X[0];
            const double h1 = X[2] - X[1];
            newY[0] = -Y[0] * (2 * h0 + h1) / (h0 * (h0 + h1)) + Y[1] * (h0 + h1) / (h0 * h1) - Y[2] * h0 / (h1 * (h0 + h1));

            double hPrev = h0;
            for (size_t i = 1; i < n - 1; ++i) {
                const double hNext = X[i + 1] - X[i];
                newY[i] = -Y[i - 1] * hNext / (hPrev * (hPrev + hNext)) + Y[i] * (hNext - hPrev) / (hPrev * hNext) + Y[i + 1] * hPrev / (hNext * (hPrev + hNext));
                hPrev = hNext;
            }

            hPrev = X[n - 2] - X[n - 3];
            const double hLast = X[n - 1] - X[n - 2];
            newY[n - 1] = Y[n - 3] * hLast / (hPrev * (hPrev + hLast)) - Y[n - 2] * (hPrev + hLast) / (hPrev * hLast) + Y[n - 1] * (hPrev + 2 * hLast) / (hLast * (hPrev + hLast));
        }
        Y = std::move(newY);
        changed = true;
    }

    void Integral() noexcept {
        const size_t n = X.size();
        if (n < 2) return;
        double prevY = Y[0];
        Y[0] = 0.0;
        for (size_t i = 1; i < n; ++i) {
            const double dx = X[i] - X[i - 1];
            const double currY = Y[i];
            Y[i] = Y[i - 1] + (prevY + currY) * 0.5 * dx;
            prevY = currY;
        }
        changed = true;
    }

    [[nodiscard]] bool canInsertPoint(double x) const noexcept {
        if (changed) update_spline();
        if (X.empty()) return true;
        const size_t k = branchless_lower_bound(X.data(), X.size(), x);
        if (k < X.size() && X[k] == x) {
            return false;
        }
        if (k > 0 && x - X[k - 1] < istep) {
            return false;
        }
        if (k < X.size() && X[k] - x < istep) {
            return false;
        }
        return true;
    }

    void Expand(int n) {
        if (changed) update_spline();
        if (X.size() < 2) return;

        const double v1X = ixmin - istep, v1Y = iymax;
        const uint32_t v1Epoch = epoch[0];
        const double v2X = ixmax + istep, v2Y = iymax;
        const uint32_t v2Epoch = epoch.back();

        double midY = (iymin + iymax) * 0.5;
        std::vector<int> indices_list;
        std::vector<int> andices_list;

        for (int step = 0; step < n; ++step) {
            if (changed) {
                update_spline();
                midY = (iymin + iymax) / 2.0;
            }
            if (X.size() < 2) break;

            auto getX = [&](int i) -> double {
                if (i == 0) return v1X;
                if (i <= static_cast<int>(X.size())) return X[i - 1];
                return v2X;
            };
            auto getY = [&](int i) -> double {
                if (i == 0) return v1Y;
                if (i <= static_cast<int>(Y.size())) return Y[i - 1];
                return v2Y;
            };
            auto getEpoch = [&](int i) -> uint32_t {
                if (i == 0) return v1Epoch;
                if (i <= static_cast<int>(epoch.size())) return epoch[i - 1];
                return v2Epoch;
            };

            indices_list.clear();
            andices_list.clear();
            const int total = static_cast<int>(X.size()) + 2;
            for (int i = 0; i < total; ++i) {
                if (getY(i) > midY) {
                    indices_list.push_back(i);
                } else {
                    andices_list.push_back(i);
                }
            }

            if (indices_list.size() > 1) {
                double maxDist = -1.0;
                int bestIdx = -1;
                for (size_t j = 0; j + 1 < indices_list.size(); ++j) {
                    const double dist = getX(indices_list[j + 1]) - getX(indices_list[j]);
                    if (dist > maxDist) {
                        maxDist = dist;
                        bestIdx = static_cast<int>(j);
                    }
                }
                if (bestIdx == -1) break;
                const int idx1 = indices_list[bestIdx], idx2 = indices_list[bestIdx + 1];
                const double midX = (getX(idx1) + getX(idx2)) * 0.5;
                if (canInsertPoint(midX)) {
                    AddPoint(midX, (getY(idx1) + getY(idx2)) * 0.5, getEpoch(idx1));
                }
            }

            if (andices_list.size() > 1) {
                double maxDist = -1.0;
                int bestIdx = -1;
                for (size_t j = 0; j + 1 < andices_list.size(); ++j) {
                    const double dist = getX(andices_list[j + 1]) - getX(andices_list[j]);
                    if (dist > maxDist) {
                        maxDist = dist;
                        bestIdx = static_cast<int>(j);
                    }
                }
                if (bestIdx == -1) break;
                const int idx1 = andices_list[bestIdx], idx2 = andices_list[bestIdx + 1];
                const double midX = (getX(idx1) + getX(idx2)) * 0.5;
                if (canInsertPoint(midX)) {
                    AddPoint(midX, (getY(idx1) + getY(idx2)) * 0.5, getEpoch(idx1));
                }
            }
        }
        changed = true;
    }

    [[nodiscard]] double GetStep() const noexcept {
        if (changed) update_spline();
        return istep;
    }
    [[nodiscard]] double GetXmin() const noexcept {
        if (changed) update_spline();
        return ixmin;
    }
    [[nodiscard]] double GetXmax() const noexcept {
        if (changed) update_spline();
        return ixmax;
    }
    [[nodiscard]] double GetYmin() const noexcept {
        if (changed) update_spline();
        return iymin;
    }
    [[nodiscard]] double GetYmax() const noexcept {
        if (changed) update_spline();
        return iymax;
    }
    [[nodiscard]] size_t GetNdots() const noexcept {
        return X.size();
    }

    std::string to_string() const {
        if (changed) update_spline();
        return std::format(
            "\nTabulated function:\n"
            "\tiOrder: {}; changed: {}\n"
            "\tixmin: {}; ixmax: {}\n"
            "\tiymin: {}; iymax: {}\n"
            "\tistep: {}\n"
            "\tPoints count: {}\n",
            Order, changed, ixmin, ixmax, iymin, iymax, istep, X.size());
    }

    void Epoch(uint32_t targetEpoch) noexcept {
        size_t w = 0;
        const size_t n = X.size();
        for (size_t r = 0; r < n; ++r) {
            if (epoch[r] >= targetEpoch) {
                if (w != r) {
                    X[w] = X[r];
                    Y[w] = Y[r];
                    epoch[w] = epoch[r];
                    indices[w] = indices[r];
                }
                w++;
            }
        }
        if (w != X.size()) {
            X.resize(w);
            Y.resize(w);
            epoch.resize(w);
            indices.resize(w);
            changed = true;
        }
    }

    bool DrawPS(const std::string& path) const {
        std::ofstream ps(path);
        if (!ps.is_open()) {
            return false;
        }

        if (changed) update_spline();

        auto [minIndex, maxIndex] = const_cast<TabulatedFunction*>(this)->NormaliseIndices();

        if (X.empty()) {
            ps << "%!PS\nshowpage\nquit\n";
            return true;
        }

        ps << R"(%!PS
	% This is the color that the grid is drawn in.
/grid_major_color {1 .6 .6} def
/grid_color {.7 1 1} def
/line_color {.5 .5 .5} def
/dot_color {.1 .1 .1} def
/radius 1 def
/set_gray_by_index {
    MaxIdx MinIdx sub dup 0 ne {
        exch MinIdx sub exch div
        1.0 exch sub 0.85 mul
    } {
        pop pop 0.0
    } ifelse
    setgray
} bind def
%% The line width used for the grid.
/grid_major_lw 1.5 def
/grid_lw .5 def

%% Every major-th line is drawn in a different color and thickness.
/major 10 def

%% Usage: dx dy w h gridwh
%% Draw a grid over a supplied width and height
/gridwh {
  4 dict begin
    /h exch def
    /w exch def
    /dy exch def
    /dx exch def
    gsave
        %% Set line width and color
        grid_lw setlinewidth
        grid_color setrgbcolor
        %% draw
        newpath
        %% vertical lines
        dx dx w {
            0 moveto
            0 h rlineto
        } for
        %% horizontal lines
        dy dy h {
            0 exch moveto
            w 0 rlineto
        } for
        stroke
        newpath
        grid_major_lw setlinewidth
        grid_major_color setrgbcolor
        %% every 10th line
        0 dx major mul w {
            0 moveto
            0 h rlineto
        } for
        0 dy major mul h {
            0 exch moveto
            w 0 rlineto
        } for
        stroke
    grestore
  end
} bind def

%% Distance between dimension point and start of witness line
/dimoffs 6 def
%% Distance that the witness line descends past the dimension line.
/dimext 20 def
%% Font for dimensions
/dimfont /Helvetica def
%% Font size
/dimscale 12 def
%% dimension text offset
/dimtextoffs 12 def
%% Dimension color
/dimcol {1 .1 .1 setrgbcolor} def
%% This defines the length of the arrow-head.
/dimhead 30 def

%% Usage: x1 y1 x2 y2 arrow_head x3 y3
%% Sees a line from x1,y1 to x2,y2 and draws an arrow head on the latter.
%% Returns x3 y3, leaving it to the user to draw the line (x1,y1)--(x3,y3).
/_arrow_head {
    9 dict begin
        /y2 exch def /x2 exch def /y1 exch def /x1 exch def
        /dx x2 x1 sub def /dy y2 y1 sub def /ang dy dx atan def
        /len dx dup mul dy dup mul add sqrt def
        /fact dimhead len 0.8 div div def
        gsave
            x2 y2 translate ang rotate
            newpath 0 0 moveto dimhead neg 4 {dup} repeat -.25 mul lineto
            .8 mul 0 lineto .25 mul lineto closepath fill
        grestore
        x2 dx fact mul sub y2 dy fact mul sub %% inside of the arrowhead
    end
} bind def

%% Usage: (text) _align_middle
/_align_middle {
	dimfont findfont dimscale scalefont setfont
	dimcol
    dup %% (text) (text)
    stringwidth pop %% (text) w
    -2 div 0 rmoveto
	show
} bind def

%% Draw a horizontal dimension
%% Usage x1 y1 x2 y2 offs (label) horizontal_dim
/horizontal_dim {
	gsave
	dimcol
    9 dict begin
        /label exch def
        /offs exch def
        /y2 exch def
        /x2 exch def
        /y1 exch def
        /x1 exch def
        /q y1 offs add def
        offs 0 ge {
            /v y1 dimoffs add def
            /w q dimext add def
        } {
            /v y1 dimoffs sub def
            /w q dimext sub def
        } ifelse
        %% Left witness line
        x1 v moveto x1 w lineto stroke
        %% Right witness line
        x2 v moveto x2 w lineto stroke
        %% arrow heads
        x2 q x1 q _arrow_head
        x1 q x2 q _arrow_head
        %% Dimension line
        moveto lineto stroke
        x1 x2 add 2 div q dimtextoffs add moveto label _align_middle
    end
	grestore
} bind def

%% Draw a vertical dimension
%% Usage x1 y1 x2 y2 offs (label) vertical_dim
/vertical_dim {
	gsave
	dimcol
    9 dict begin
        /label exch def
        /offs exch def
        /y2 exch def
        /x2 exch def
        /y1 exch def
        /x1 exch def
        /q x1 offs add def
        offs 0 ge {
            /v x1 dimoffs add def
            /w q dimext add def
        } {
            /v x1 dimoffs sub def
            /w q dimext sub def
        } ifelse
        %% Bottom witness line
        v y1 moveto w y1 lineto stroke
        %% Top witness line
        v y2 moveto w y2 lineto stroke
        %% arrow heads
        q y2 q y1 _arrow_head
        q y1 q y2 _arrow_head
        %% Dimension line
        moveto lineto stroke
        %% Rotated label
        q dimtextoffs sub y1 y2 add 2 div moveto
        gsave 90 rotate label _align_middle grestore
    end
	grestore
} bind def
)";

        ps << "\n/XValues [\n";
        const double xRange = (ixmax != ixmin) ? (ixmax - ixmin) : 1.0;
        std::string buffer;
        buffer.reserve(X.size() * 32);
        for (size_t i = 0; i < X.size(); ++i) {
            std::format_to(std::back_inserter(buffer), " {}\t% {}\n", (X[i] - ixmin) / xRange, i);
        }
        ps << buffer;
        ps << "] def\n";

        ps << "/YValues [\n";
        buffer.clear();
        for (size_t i = 0; i < Y.size(); ++i) {
            std::format_to(std::back_inserter(buffer), " {}\t% {}", Y[i], i);
            if (i > 0 && i < X.size() - 1) {
                const double yPrev = Y[i - 1];
                const double yNext = Y[i + 1];
                const double xPrev = X[i - 1];
                const double xNext = X[i + 1];
                const double xCurr = X[i];
                const double dy = yNext - yPrev;
                const double dx1 = xCurr - xPrev;
                const double dx2 = xNext - xPrev;
                const double val = (dx2 != 0.0) ? (yPrev + dy * dx1 / dx2) : yPrev;
                std::format_to(std::back_inserter(buffer), "\t% interp: {}", val);
            }
            buffer.push_back('\n');
        }
        ps << buffer;
        ps << "] def\n";

        ps << "/ColorValues [\n";
        buffer.clear();
        for (size_t i = 0; i < indices.size(); ++i) {
            std::format_to(std::back_inserter(buffer), " {}\t% {}\n", indices[i], i);
        }
        ps << buffer;
        ps << "] def\n";

        ps << std::format("/MinIdx {} def\n", minIndex);
        ps << std::format("/MaxIdx {} def\n", maxIndex);
        ps << "/Xmin 0 def\n";
        ps << "/Xmax 1 def\n";
        ps << std::format("/Ymin {} def\n", iymin);
        ps << std::format("/Ymax {} def\n", iymax);

        ps << R"(
/Xsize Xmax Xmin sub def
/Ysize Ymax Ymin sub def

/w currentpagedevice /PageSize get 0 get def
/h currentpagedevice /PageSize get 1 get def

w 10 div h 10 div w h gridwh

/Translate { %% x y Translate
	Ymin sub Ysize div h mul
	exch
	Xmin sub Xsize div w mul
	exch 
} bind def

%% lines
1 1 XValues length 1 sub {  %% i
newpath
ColorValues 1 index get set_gray_by_index
XValues
1 index 1 sub
get
YValues
2 index 1 sub
get
Translate
moveto
XValues
1 index
get
YValues
2 index
get
Translate
lineto
stroke
pop
} for

%% dots
newpath
ColorValues 0 get set_gray_by_index
XValues 0 get YValues 0 get
Translate
radius 0 360 arc
stroke
1 1 XValues length 1 sub {  %% i
ColorValues 1 index get set_gray_by_index
XValues
1 index
get
YValues
2 index
get
Translate
radius 0 360 arc
stroke
pop
} for
)";

        ps << std::format("0 5 w 5 10 ({} - {}) horizontal_dim\n", ixmin, ixmax);
        ps << std::format("5 0 5 h 20 ({} - {}) vertical_dim\n", iymin, iymax);
        ps << "\nshowpage\nquit\n";

        return true;
    }

    void FromDump(const tabulatedfunction::Dump& d) {
        Order = d.order;
        trapolation = d.trapolation;

        const size_t n = std::min(d.x.size(), d.y.size());
        if (n == 0) {
            Clear();
            Order = d.order;
            trapolation = d.trapolation;
            return;
        }

        X.assign(d.x.begin(), d.x.begin() + n);
        Y.assign(d.y.begin(), d.y.begin() + n);
        epoch.assign(n, 0);
        const size_t epCount = std::min(n, d.epoch.size());
        std::copy_n(d.epoch.begin(), epCount, epoch.begin());

        if (!std::ranges::is_sorted(X)) {
            struct Point {
                double x, y;
                uint32_t epoch;
            };
            std::vector<Point> pts(n);
            for (size_t i = 0; i < n; ++i) {
                pts[i] = {X[i], Y[i], epoch[i]};
            }
            std::ranges::sort(pts, [](const Point& a, const Point& b) {
                return a.x < b.x;
            });
            for (size_t i = 0; i < n; ++i) {
                X[i] = pts[i].x;
                Y[i] = pts[i].y;
                epoch[i] = pts[i].epoch;
            }
        }

        if (X.size() > 1) {
            size_t k = 0;
            size_t count = 1;
            double sumY = Y[0];
            for (size_t i = 1; i < X.size(); ++i) {
                if (X[i] == X[k]) {
                    sumY += Y[i];
                    count++;
                    if (epoch[i] > epoch[k]) {
                        epoch[k] = epoch[i];
                    }
                } else {
                    Y[k] = sumY / static_cast<double>(count);
                    k++;
                    X[k] = X[i];
                    Y[k] = Y[i];
                    epoch[k] = epoch[i];
                    sumY = Y[k];
                    count = 1;
                }
            }
            Y[k] = sumY / static_cast<double>(count);
            X.resize(k + 1);
            Y.resize(k + 1);
            epoch.resize(k + 1);
        }

        indices.resize(X.size());
        nextIndex = 1;
        for (size_t i = 0; i < indices.size(); ++i) {
            indices[i] = nextIndex++;
        }

        update_spline();
    }

    [[nodiscard]] tabulatedfunction::Dump ToDump() const {
        return tabulatedfunction::Dump{
            .order = Order,
            .trapolation = trapolation,
            .x = X,
            .y = Y,
            .epoch = epoch,
        };
    }

    [[nodiscard]] tabulatedfunction::Dump Dump() const {
        return ToDump();
    }

    [[nodiscard]] std::string ToJSON() const {
        nlohmann::json j;
        to_json(j, ToDump());
        return j.dump();
    }

    bool FromJSON(std::string_view json_str) {
        try {
            const auto j = nlohmann::json::parse(json_str);
            tabulatedfunction::Dump d;
            from_json(j, d);
            FromDump(d);
            return true;
        } catch (...) {
            return false;
        }
    }
};

inline void to_json(nlohmann::json& j, const TabulatedFunction& f) {
    to_json(j, f.ToDump());
}

inline void from_json(const nlohmann::json& j, TabulatedFunction& f) {
    Dump d;
    from_json(j, d);
    f.FromDump(d);
}

} // namespace tabulatedfunction