# Tabulated Function

A high-performance tabulated function interpolation and numerical calculus library implemented in both **Go** and modern **C++23** (header-only).

Originally based on [bravesoftdz/Table-func-lib](https://github.com/bravesoftdz/Table-func-lib).

> [!WARNING]
> **Warning:** The C++ version is not yet tested in real applications and must not be used in production environments.

---

## Features

- **Fast Point Management**: Binary-search point insertion, deduplication, and fast-paths for sequential appending/prepending.
- **Interpolation & Extrapolation Modes**:
  - `Linear`
  - `Shift`
  - `MinMax`
  - `Opposite`
  - `Nearest`
  - `Cosine`
- **Calculus Operations**:
  - Numerical differentiation (3-point finite differences)
  - Definite integration (Simpson's 3-point rule with trapezoidal boundary correction)
  - Indefinite integration (cumulative trapezoidal)
- **Transformation & Manipulation**:
  - In-place 3-point moving average smoothing (`Smooth`)
  - Midpoint expansion and densification (`MorePoints`, `Expand`)
  - Function merging and pointwise multiplication (`Merge`, `Multiply`)
  - Scalar multiplication and normalization
  - Epoch-based point filtering
- **Visualization**: PostScript (`.ps`) vector chart generation (`DrawPS`).
- **Serialization**: State serialization via `Dump` structs and JSON.

---

## Getting Started (C++23)

The C++ implementation is header-only and requires a C++23-compliant compiler (e.g., GCC 13+, Clang 17+, MSVC 19.36+).

### CMake Integration

Include the directory in your `CMakeLists.txt`:

```cmake
add_subdirectory(tabulated-function)
target_link_libraries(your_target PRIVATE tabulatedfunction::tabulated-function)
```

### Quick C++ Example

```cpp
#include "tabulated-function.hpp"
#include <iostream>

int main() {
    tabulatedfunction::TabulatedFunction f;
    f.SetOrder(1); // Linear interpolation

    f.AddPoint(0.0, 0.0);
    f.AddPoint(1.0, 10.0);
    f.AddPoint(2.0, 20.0);

    std::cout << "f(0.5) = " << f.F(0.5) << "\n";       // 5.0
    std::cout << "Integral = " << f.Integrate() << "\n"; // 20.0

    return 0;
}
```

### Building and Running Tests

Use `bootstrap.sh` to initialize vcpkg, build, and run the test suite:

```bash
chmod +x bootstrap.sh
./bootstrap.sh
```

Or manually via CMake:

```bash
cmake -B build -S . -DBUILD_TESTING=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

---

## Getting Started (Go)

### Quick Go Example

```go
package main

import (
    "fmt"
    tabulatedfunction "github.com/Maxime2/tabulated-function"
)

func main() {
    f := tabulatedfunction.New()
    f.AddPoint(0.0, 0.0, 0)
    f.AddPoint(1.0, 10.0, 0)
    f.AddPoint(2.0, 20.0, 0)

    fmt.Println("f(0.5) =", f.F(0.5)) // 5.0
}
```

Run tests:
```bash
go test ./... -v
```

---

## License

This project is licensed under the MIT License.