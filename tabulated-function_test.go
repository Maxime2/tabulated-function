package tabulatedfunction

import (
	"encoding/json"
	"math"
	"os"
	"testing"
)

const float64EqualityThreshold = 1e-9

func almostEqual(a, b float64) bool {
	if math.IsNaN(a) && math.IsNaN(b) {
		return true
	}
	return math.Abs(a-b) <= float64EqualityThreshold
}

func TestNew(t *testing.T) {
	f := New()
	if f == nil {
		t.Fatal("New() returned nil")
	}
	if f.Order != 1 {
		t.Errorf("Expected default order to be 1, got %d", f.Order)
	}
	if f.Trapolation != TrapolationLinear {
		t.Errorf("Expected default trapolation to be TrapolationLinear, got %v", f.Trapolation)
	}
	if f.GetNdots() != 0 {
		t.Errorf("Expected new function to have 0 points, got %d", f.GetNdots())
	}
	if f.changed != false {
		t.Error("Expected new function to have changed=false")
	}
}

func TestAddPointAndF(t *testing.T) {
	f := New()
	f.SetOrder(1) // Use linear for simplicity in this test

	// Test F on empty function
	if !math.IsNaN(f.F(0)) {
		t.Errorf("F(0) on empty function should be NaN, got %v", f.F(0))
	}

	// Add points
	f.AddPoint(1, 10, 0)
	f.AddPoint(3, 30, 0)
	f.AddPoint(2, 20, 0) // Add out of order

	if f.GetNdots() != 3 {
		t.Fatalf("Expected 3 points, got %d", f.GetNdots())
	}

	// Check if sorted
	if !(f.X[0] == 1 && f.X[1] == 2 && f.X[2] == 3) {
		t.Errorf("Points are not sorted correctly: %v", f.X)
	}

	testCases := []struct {
		name     string
		x        float64
		expected float64
	}{
		{"Exact Point 1", 1.0, 10.0},
		{"Exact Point 2", 2.0, 20.0},
		{"Exact Point 3", 3.0, 30.0},
		{"Interpolation", 1.5, 15.0},
		{"Interpolation 2", 2.5, 25.0},
		{"Extrapolation Left", 0.0, 10.0},
		{"Extrapolation Right", 4.0, 30.0},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			y := f.F(tc.x)
			if !almostEqual(y, tc.expected) {
				t.Errorf("F(%v) = %v; want %v", tc.x, y, tc.expected)
			}
		})
	}

	// Test adding an existing point (should overwrite and return the old value)
	// The original point is (2, 20). We add (2, 22).
	oldY := f.AddPoint(2, 22, 1)
	if !almostEqual(oldY, 20.0) {
		t.Errorf("AddPoint(2, 22) returned %v; want 20.0 (old value)", oldY)
	}
	expectedY := 22.0
	if !almostEqual(f.F(2), expectedY) {
		t.Errorf("F(2) after adding existing point = %v; want %v", f.F(2), expectedY)
	}
}

func TestGetters(t *testing.T) {
	f := New()
	f.AddPoint(0, 5, 0)
	f.AddPoint(10, -5, 0)
	f.AddPoint(5, 15, 0)

	// Force update
	_ = f.F(1)

	if !almostEqual(f.GetXmin(), 0) {
		t.Errorf("GetXmin() = %v; want 0", f.GetXmin())
	}
	if !almostEqual(f.GetXmax(), 10) {
		t.Errorf("GetXmax() = %v; want 10", f.GetXmax())
	}
	if !almostEqual(f.GetYmin(), -5) {
		t.Errorf("GetYmin() = %v; want -5", f.GetYmin())
	}
	if !almostEqual(f.GetYmax(), 15) {
		t.Errorf("GetYmax() = %v; want 15", f.GetYmax())
	}
	if f.GetNdots() != 3 {
		t.Errorf("GetNdots() = %v; want 3", f.GetNdots())
	}
	if !almostEqual(f.GetStep(), 5) { // min step is between (0,x) and (5,x) or (5,x) and (10,x)
		t.Errorf("GetStep() = %v; want 5", f.GetStep())
	}
}

func TestOrders(t *testing.T) {
	f := New()
	f.AddPoint(0, 0, 0)
	f.AddPoint(1, 1, 0)
	f.AddPoint(2, 0, 0)

	t.Run("Order 0 (Step function)", func(t *testing.T) {
		f.SetOrder(0)
		if !almostEqual(f.F(0.5), 0) {
			t.Errorf("F(0.5) = %v; want 0", f.F(0.5))
		}
		if !almostEqual(f.F(1.5), 1) {
			t.Errorf("F(1.5) = %v; want 1", f.F(1.5))
		}
		if !almostEqual(f.F(-1), 0) { // Extrapolation
			t.Errorf("F(-1) = %v; want 0", f.F(-1))
		}
		if !almostEqual(f.F(3), 0) { // Extrapolation
			t.Errorf("F(3) = %v; want 0", f.F(3))
		}
	})

	t.Run("Order 1 (Linear)", func(t *testing.T) {
		f.SetOrder(1)
		if !almostEqual(f.F(0.5), 0.5) {
			t.Errorf("F(0.5) = %v; want 0.5", f.F(0.5))
		}
		if !almostEqual(f.F(1.5), 0.5) {
			t.Errorf("F(1.5) = %v; want 0.5", f.F(1.5))
		}
	})
}

func TestTrapolationLinear(t *testing.T) {
	f := New()
	f.SetTrapolation(TrapolationLinear)
	f.AddPoint(0, 0, 0)
	f.AddPoint(2, 2, 0)

	// Test interpolation
	y := f.F(1)
	if !almostEqual(y, 1) {
		t.Errorf("Linear interpolation F(1) = %v; want 1", y)
	}
}

func TestTrapolateDirect(t *testing.T) {
	f := New()
	// Test empty Trapolate
	if !math.IsNaN(f.Trapolate(10, TrapolationLinear)) {
		t.Errorf("Trapolate on empty function should be NaN, got %v", f.Trapolate(10, TrapolationLinear))
	}

	f.AddPoint(0, 10, 0)
	f.AddPoint(5, 20, 0)
	f.AddPoint(10, 30, 0)

	// Exact match at left boundary (k=0)
	if y := f.Trapolate(0, TrapolationShift); !almostEqual(y, 15.0) { // average with right neighbor: (10+20)/2
		t.Errorf("Trapolate exact match k=0 with Shift = %v; want 15.0", y)
	}
	// Exact match at right boundary (k=l-1)
	if y := f.Trapolate(10, TrapolationShift); !almostEqual(y, 25.0) { // average with left neighbor: (20+30)/2
		t.Errorf("Trapolate exact match k=l-1 with Shift = %v; want 25.0", y)
	}
	// Exact match at interior knot (k=1)
	if y := f.Trapolate(5, TrapolationShift); !almostEqual(y, 20.0) { // (f.Y[0]+f.Y[2])/2 = (10+30)/2
		t.Errorf("Trapolate exact interior match with Shift = %v; want 20.0", y)
	}
}

func TestInterpolateEdgeCasesAndPanic(t *testing.T) {
	f := New()
	f.X = []float64{1.0, 1.0} // Artificial dx = 0 to test div-by-zero guard
	f.Y = []float64{10.0, 20.0}

	// Linear with dx == 0 should return f.Y[left]
	if y := f._interpolate(1.0, 0, 1, TrapolationLinear); !almostEqual(y, 10.0) {
		t.Errorf("Linear interpolation with dx=0 gave %v; want 10.0", y)
	}
	// Cosine with dx == 0 should return f.Y[left]
	if y := f._interpolate(1.0, 0, 1, TrapolationCosine); !almostEqual(y, 10.0) {
		t.Errorf("Cosine interpolation with dx=0 gave %v; want 10.0", y)
	}

	// TrapolationMinMax branches
	f2 := New()
	f2.AddPoint(0, 0, 0)
	f2.AddPoint(10, 10, 0)
	_ = f2.Trapolate(5, TrapolationMinMax)
	// Out-of-bounds indices handling in MinMax
	yMinMax := f2._interpolate(5, -1, 5, TrapolationMinMax)
	if math.IsNaN(yMinMax) {
		t.Errorf("MinMax with out of bounds indices produced NaN")
	}

	// Unhandled trapolation type should panic
	defer func() {
		if r := recover(); r == nil {
			t.Errorf("Expected panic on unhandled Trapolation type, but did not panic")
		}
	}()
	_ = f2._interpolate(5, 0, 1, Trapolation(999))
}

func TestCalculusBorderCases(t *testing.T) {
	// 1. Derivative with 0, 1, and 2 points
	f0 := New()
	f0.Derivative()
	if f0.GetNdots() != 0 {
		t.Errorf("Derivative on empty should have 0 points, got %d", f0.GetNdots())
	}

	f1 := New()
	f1.AddPoint(1, 50, 0)
	f1.Derivative()
	if f1.Y[0] != 0 {
		t.Errorf("Derivative on single point should set Y[0] = 0, got %v", f1.Y[0])
	}

	f2 := New()
	f2.AddPoint(0, 10, 0)
	f2.AddPoint(2, 20, 0)
	f2.Derivative()
	if !almostEqual(f2.Y[0], 5.0) || !almostEqual(f2.Y[1], 5.0) {
		t.Errorf("Derivative on 2 points should yield slope 5.0, got (%v, %v)", f2.Y[0], f2.Y[1])
	}

	// 2. Integral with 0, 1, and 2 points
	f0.Integral()
	if f0.GetNdots() != 0 {
		t.Errorf("Integral on empty should have 0 points, got %d", f0.GetNdots())
	}
	f1.Integral()
	if f1.Y[0] != 0 {
		t.Errorf("Integral on single point should be no-op, got %v", f1.Y[0])
	}
	f2Int := New()
	f2Int.AddPoint(0, 4, 0)
	f2Int.AddPoint(3, 4, 0)
	f2Int.Integral() // integral of constant 4 from 0 to 3 should be [0, 12]
	if !almostEqual(f2Int.Y[0], 0) || !almostEqual(f2Int.Y[1], 12.0) {
		t.Errorf("Integral on 2 points = (%v, %v); want (0, 12.0)", f2Int.Y[0], f2Int.Y[1])
	}

	// 3. Integrate with n=0, 1, 2, 3, 4 and non-uniform intervals
	if f0.Integrate() != 0 {
		t.Errorf("Integrate on empty = %v; want 0", f0.Integrate())
	}
	if f1.Integrate() != 0 {
		t.Errorf("Integrate on 1 point = %v; want 0", f1.Integrate())
	}
	// n = 2 (trapezoid test)
	fTrap := New()
	fTrap.AddPoint(0, 2, 0)
	fTrap.AddPoint(4, 6, 0) // Area = 4 * (2+6)/2 = 16
	if !almostEqual(fTrap.Integrate(), 16.0) {
		t.Errorf("Integrate on 2 points = %v; want 16.0", fTrap.Integrate())
	}
	// n = 4 (Simpson's 3-point rule + 1 trapezoid step, non-uniform intervals)
	fNonUniform := New()
	fNonUniform.AddPoint(0, 0, 0)
	fNonUniform.AddPoint(1, 1, 0)  // h1 = 1
	fNonUniform.AddPoint(3, 9, 0)  // h2 = 2
	fNonUniform.AddPoint(4, 16, 0) // h3 = 1 (trapezoid segment: (9+16)/2 * 1 = 12.5)
	val := fNonUniform.Integrate()
	if math.IsNaN(val) || val <= 0 {
		t.Errorf("Integrate on non-uniform grid produced invalid result: %v", val)
	}
}

func TestNormaliseBorderCases(t *testing.T) {
	// All zeros: ym == 0, should not divide by zero or change anything
	fZeros := New()
	fZeros.AddPoint(0, 0, 0)
	fZeros.AddPoint(1, 0, 0)
	fZeros.Normalise()
	if fZeros.Y[0] != 0 || fZeros.Y[1] != 0 {
		t.Errorf("Normalise on all zeros altered values: %v", fZeros.Y)
	}

	// Negative maximum magnitude
	fNeg := New()
	fNeg.AddPoint(0, -2, 0)
	fNeg.AddPoint(1, -10, 0)
	fNeg.Normalise()
	if !almostEqual(fNeg.Y[0], -0.2) || !almostEqual(fNeg.Y[1], -1.0) {
		t.Errorf("Normalise on negative values = %v; want [-0.2, -1.0]", fNeg.Y)
	}
}

func TestMultiplyBorderCases(t *testing.T) {
	empty := New()
	populated := New()
	populated.AddPoint(0, 5, 0)
	populated.AddPoint(1, 10, 0)

	// Empty receiver
	empty.Multiply(populated)
	if empty.GetNdots() != 0 {
		t.Errorf("Multiplying empty receiver resulted in %d points, want 0", empty.GetNdots())
	}

	// Multiplying populated by empty should clear receiver
	populated.Multiply(New())
	if populated.GetNdots() != 0 {
		t.Errorf("Multiplying by empty did not clear receiver, got %d points", populated.GetNdots())
	}
}

func TestCanInsertPointAndExpandBorderCases(t *testing.T) {
	f := New()
	if !f.canInsertPoint(5.0) {
		t.Errorf("canInsertPoint on empty function should return true")
	}

	f.AddPoint(0, 0, 0)
	f.AddPoint(10, 10, 0) // istep = 10
	_ = f.GetStep()

	// Inserting exactly on existing point
	if f.canInsertPoint(0.0) || f.canInsertPoint(10.0) {
		t.Errorf("canInsertPoint on existing knot should return false")
	}
	// Inserting within istep distance
	if f.canInsertPoint(2.0) { // 2.0 - 0.0 < 10.0
		t.Errorf("canInsertPoint within istep distance should return false")
	}

	// Expand with n <= 0
	dotsBefore := f.GetNdots()
	f.Expand(0)
	if f.GetNdots() != dotsBefore {
		t.Errorf("Expand(0) altered point count: got %d, want %d", f.GetNdots(), dotsBefore)
	}
}

func TestDrawPSBorderCases(t *testing.T) {
	// 1. Single point PS rendering (xRange == 0 and yRange == 0)
	fSingle := New()
	fSingle.AddPoint(5, 5, 0)
	tmp1, err := os.CreateTemp("", "test_single_ps_*.ps")
	if err != nil {
		t.Fatalf("Failed to create temp file: %v", err)
	}
	defer os.Remove(tmp1.Name())
	tmp1.Close()

	if err := fSingle.DrawPS(tmp1.Name()); err != nil {
		t.Errorf("DrawPS on single point failed: %v", err)
	}

	// 2. Invalid file path
	err = fSingle.DrawPS("/invalid_nonexistent_directory/file.ps")
	if err == nil {
		t.Errorf("DrawPS on invalid path should return error, got nil")
	}
}

func TestFromDumpAndJSONBorderCases(t *testing.T) {
	// Invalid JSON string
	f := New()
	if err := f.UnmarshalJSON([]byte("invalid json")); err == nil {
		t.Errorf("UnmarshalJSON on invalid json should return error, got nil")
	}

	// String output representation
	str := f.String()
	if len(str) == 0 {
		t.Errorf("String() returned empty string")
	}
}

func TestTrapolationOpposite(t *testing.T) {
	f := New()
	f.SetTrapolation(TrapolationOpposite)
	f.AddPoint(0, 10, 0)
	f.AddPoint(10, 20, 0)

	// Test interpolation: returns the Y value of the opposite point.
	// Midpoint is (0+10)/2 = 5.
	// For xi <= 5, it should return the right point's Y value (20).
	// For xi > 5, it should return the left point's Y value (10).
	t.Run("Interpolation at midpoint", func(t *testing.T) {
		y := f.F(5)
		expected := 20.0 // xi <= midpoint, so return right Y
		if !almostEqual(y, expected) {
			t.Errorf("Opposite interpolation F(5) = %v; want %v", y, expected)
		}
	})
	t.Run("Interpolation closer to right point", func(t *testing.T) {
		y := f.F(6)
		expected := 10.0 // xi > midpoint, so return left Y
		if !almostEqual(y, expected) {
			t.Errorf("Opposite interpolation F(6) = %v; want %v", y, expected)
		}
	})
	t.Run("Interpolation closer to left point", func(t *testing.T) {
		y := f.F(4)
		expected := 20.0 // xi <= midpoint, so return right Y
		if !almostEqual(y, expected) {
			t.Errorf("Opposite interpolation F(4) = %v; want %v", y, expected)
		}
	})

	// Test extrapolation
	// For TrapolationOpposite, extrapolation now also uses the opposite logic.
	// For x=15, it extrapolates from the rightmost point (10, 20).
	// Opposite value should be: Ymin + Ymax - boundary_Y = 10 + 20 - 20 = 10.
	t.Run("Extrapolation", func(t *testing.T) {
		_ = f.GetYmax() // Force update of ymin/ymax
		y := f.F(15)
		expected := 10.0
		if !almostEqual(y, expected) {
			t.Errorf("Opposite extrapolation F(15) = %v; want %v", y, expected)
		}
	})
}

func TestTrapolationNearest(t *testing.T) {
	f := New()
	f.SetTrapolation(TrapolationNearest)
	f.AddPoint(0, 10, 0)
	f.AddPoint(10, 20, 0)

	testCases := []struct {
		x        float64
		expected float64
	}{
		{2, 10},
		{8, 20},
		{5, 20}, // Exactly in middle, implementation returns right
		{-5, 10},
		{15, 20},
	}

	for _, tc := range testCases {
		if y := f.F(tc.x); !almostEqual(y, tc.expected) {
			t.Errorf("Nearest F(%v) = %v; want %v", tc.x, y, tc.expected)
		}
	}
}

func TestTrapolationCosine(t *testing.T) {
	f := New()
	f.SetTrapolation(TrapolationCosine)
	f.AddPoint(0, 0, 0)
	f.AddPoint(10, 100, 0)

	// At midpoint 5, mu=0.5. mu2 = (1 - cos(pi/2))/2 = 0.5. Result should be 50.
	// Midpoint test verifies the cosine transition curve is centered correctly.
	if y := f.F(5); !almostEqual(y, 50) {
		t.Errorf("Cosine F(5) = %v; want 50", y)
	}
}

func TestTrapolationShift(t *testing.T) {
	f := New()
	f.SetTrapolation(TrapolationShift)
	f.AddPoint(0, 10, 0)
	f.AddPoint(2, 20, 0)

	// Interpolation between 0 and 2 should return the midpoint average: (10 + 20) / 2 = 15
	if y := f.F(1); !almostEqual(y, 15.0) {
		t.Errorf("Shift F(1) = %v; want 15.0", y)
	}

	// Extrapolation outside boundary
	if y := f.F(-1); !almostEqual(y, 10.0) {
		t.Errorf("Shift extrapolation left F(-1) = %v; want 10.0", y)
	}
	if y := f.F(3); !almostEqual(y, 20.0) {
		t.Errorf("Shift extrapolation right F(3) = %v; want 20.0", y)
	}
}

func TestTrapolationMinMax(t *testing.T) {
	f := New()
	f.SetTrapolation(TrapolationMinMax)
	f.AddPoint(0, 0, 0)
	f.AddPoint(10, 10, 0)

	// Force update_spline to populate bounds
	_ = f.F(5)

	// Should produce valid numeric value within [iymin, iymax]
	y := f.F(5)
	if math.IsNaN(y) || y < 0 || y > 10 {
		t.Errorf("MinMax F(5) = %v; want value in [0, 10]", y)
	}
}

func TestTrapolationCosineExtrapolation(t *testing.T) {
	f := New()
	f.SetTrapolation(TrapolationCosine)
	f.AddPoint(0, 10, 0)
	f.AddPoint(10, 20, 0)

	// Extrapolation beyond bounds should clamp to boundary points
	if y := f.F(-5); !almostEqual(y, 10.0) {
		t.Errorf("Cosine left extrapolation F(-5) = %v; want 10.0", y)
	}
	if y := f.F(15); !almostEqual(y, 20.0) {
		t.Errorf("Cosine right extrapolation F(15) = %v; want 20.0", y)
	}
}

func TestLoadConstant(t *testing.T) {
	f := New()
	f.LoadConstant(100, -5, 5)
	if f.GetNdots() != 1 {
		t.Fatalf("LoadConstant should create 1 point, got %d", f.GetNdots())
	}
	if !almostEqual(f.F(0), 100) {
		t.Errorf("F(0) = %v; want 100", f.F(0))
	}
	if !almostEqual(f.F(-10), 100) { // Extrapolation
		t.Errorf("F(-10) = %v; want 100", f.F(-10))
	}
}

func TestClear(t *testing.T) {
	f := New()
	f.AddPoint(1, 1, 0)
	f.Clear()
	if f.GetNdots() != 0 {
		t.Errorf("GetNdots() after Clear() = %d; want 0", f.GetNdots())
	}
	if !math.IsNaN(f.F(1)) {
		t.Errorf("F(1) after Clear() should be NaN, got %v", f.F(1))
	}
}

func TestJSON(t *testing.T) {
	f1 := New()
	f1.AddPoint(0, 0, 1)
	f1.AddPoint(1, 1, 2)
	f1.SetOrder(1)
	f1.SetTrapolation(TrapolationLinear)

	jsonData, err := json.Marshal(f1)
	if err != nil {
		t.Fatalf("json.Marshal failed: %v", err)
	}

	var f2 TabulatedFunction
	err = json.Unmarshal(jsonData, &f2)
	if err != nil {
		t.Fatalf("json.Unmarshal failed: %v", err)
	}

	if f1.Trapolation != f2.Trapolation {
		t.Errorf("Trapolation mismatch: original=%v, unmarshaled=%v", f1.Trapolation, f2.Trapolation)
	}
	if f1.Order != f2.Order {
		t.Errorf("Order mismatch: original=%d, unmarshaled=%d", f1.Order, f2.Order)
	}
	if f1.GetNdots() != f2.GetNdots() {
		t.Fatalf("Points count mismatch: original=%d, unmarshaled=%d", f1.GetNdots(), f2.GetNdots())
	}
	for i := range f1.X {
		if !almostEqual(f1.X[i], f2.X[i]) || !almostEqual(f1.Y[i], f2.Y[i]) || f1.epoch[i] != f2.epoch[i] {
			t.Errorf("Point mismatch at index %d: original=(%v,%v,%v) unmarshaled=(%v,%v,%v)", i, f1.X[i], f1.Y[i], f1.epoch[i], f2.X[i], f2.Y[i], f2.epoch[i])
		}
	}
}

func TestDerivativeAndIntegral(t *testing.T) {
	// y = x^2
	f := New()
	f.AddPoint(-2, 4, 0)
	f.AddPoint(-1, 1, 0)
	f.AddPoint(0, 0, 0)
	f.AddPoint(1, 1, 0)
	f.AddPoint(2, 4, 0)
	f.SetOrder(1) // Order doesn't matter for numerical calculus methods

	// Test definite integral: Integral of x^2 from -2 to 2 is 16/3
	integral := f.Integrate()
	// Simpson's rule is exact for quadratics.
	if !almostEqual(integral, 16.0/3.0) {
		t.Errorf("Integral of x^2 from -2 to 2 is %v, analytical is %v", integral, 16.0/3.0)
	}

	// Test Derivative: should be approx y' = 2x
	f.Derivative()
	// 3-point finite difference is exact for quadratics.
	if !almostEqual(f.F(0), 0) {
		t.Errorf("Derivative at F(0) is %v, want 0", f.F(0))
	}
	if !almostEqual(f.F(1), 2.0) {
		t.Errorf("Derivative at F(1) is %v, want 2.0", f.F(1))
	}

	// Test Integral (indefinite): integral of 2x should be y = x^2 + C
	// The trapezoidal rule used is an approximation.
	f.Integral()
	y_neg1 := f.F(-1)
	y0 := f.F(0)
	y1 := f.F(1)
	// Check second difference to verify quadratic shape, independent of integration constant
	if math.Abs((y1-y0)-(y0-y_neg1)-2.0) > 0.5 {
		t.Errorf("Shape of integral is incorrect. Second difference at 0 is %v, want ~2.0", (y1-y0)-(y0-y_neg1))
	}
}

func TestFromDumpUnsortedAndDeduplication(t *testing.T) {
	// Test FromDump with unsorted X points and duplicates
	d := &Dump{
		Order:       1,
		Trapolation: TrapolationLinear,
		X:           []float64{3.0, 1.0, 2.0, 2.0},
		Y:           []float64{30.0, 10.0, 20.0, 40.0},
		Epoch:       []uint32{0, 0, 1, 2},
	}

	f := New()
	f.FromDump(d)

	// After sorting: (1, 10), (2, 20), (2, 40), (3, 30)
	// After deduplicating (2, 20) and (2, 40): Y=(20+40)/2=30, Epoch=max(1,2)=2
	if f.GetNdots() != 3 {
		t.Fatalf("Expected 3 points after deduplication, got %d", f.GetNdots())
	}
	if !almostEqual(f.X[0], 1.0) || !almostEqual(f.Y[0], 10.0) {
		t.Errorf("Point 0 mismatch: got (%v, %v), want (1.0, 10.0)", f.X[0], f.Y[0])
	}
	if !almostEqual(f.X[1], 2.0) || !almostEqual(f.Y[1], 30.0) || f.epoch[1] != 2 {
		t.Errorf("Point 1 mismatch: got (%v, %v, ep %v), want (2.0, 30.0, ep 2)", f.X[1], f.Y[1], f.epoch[1])
	}
	if !almostEqual(f.X[2], 3.0) || !almostEqual(f.Y[2], 30.0) {
		t.Errorf("Point 2 mismatch: got (%v, %v), want (3.0, 30.0)", f.X[2], f.Y[2])
	}
}

func TestJSONEmpty(t *testing.T) {
	f := New()
	data, err := json.Marshal(f)
	if err != nil {
		t.Fatalf("Marshal empty failed: %v", err)
	}
	var f2 TabulatedFunction
	if err := json.Unmarshal(data, &f2); err != nil {
		t.Fatalf("Unmarshal empty failed: %v", err)
	}
	if f2.GetNdots() != 0 {
		t.Errorf("Expected 0 points, got %d", f2.GetNdots())
	}
}

func TestMorePoints(t *testing.T) {
	f := New()
	f.SetOrder(1)
	f.AddPoint(0, 0, 0)
	f.AddPoint(2, 4, 0)

	f.MorePoints()

	if f.GetNdots() != 3 {
		t.Fatalf("MorePoints should have added 1 point, got %d total", f.GetNdots())
	}
	if !almostEqual(f.X[1], 1.0) {
		t.Errorf("New point has X=%v, want 1.0", f.X[1])
	}
	if !almostEqual(f.Y[1], 2.0) {
		t.Errorf("New point has Y=%v, want 2.0", f.Y[1])
	}
}

func TestEpoch(t *testing.T) {
	f := New()
	f.AddPoint(0, 0, 0)
	f.AddPoint(1, 1, 1)
	f.AddPoint(2, 2, 2)
	f.AddPoint(3, 3, 3)

	f.Epoch(2)

	if f.GetNdots() != 2 {
		t.Fatalf("After Epoch(2), expected 2 points, got %d", f.GetNdots())
	}
	if f.X[0] != 2 || f.X[1] != 3 {
		t.Errorf("Remaining points are incorrect: %v", f.X)
	}

	// Purge all points
	f.Epoch(10)
	if f.GetNdots() != 0 {
		t.Errorf("Expected 0 points after Epoch(10), got %d", f.GetNdots())
	}
}

func TestSmooth(t *testing.T) {
	f := New()
	f.AddPoint(0, 0, 0)
	f.AddPoint(1, 10, 0) // A noisy point
	f.AddPoint(2, 4, 0)

	y_before := f.Y[1]
	f.Smooth()
	y_after := f.Y[1]

	if almostEqual(y_before, y_after) {
		t.Errorf("Smooth() did not change the Y value of the middle point. Before: %v, After: %v", y_before, y_after)
	}
	// Based on manual calculation of the 3-point average: (0.0 + 10.0 + 4.0) / 3.0
	expected_y := 14.0 / 3.0
	if !almostEqual(y_after, expected_y) {
		t.Errorf("Smooth() produced %v, want %v", y_after, expected_y)
	}
}

func TestMultiply(t *testing.T) {
	// f1(x) = 2
	f1 := New()
	f1.AddPoint(0, 2, 0)
	f1.AddPoint(5, 2, 0)
	f1.SetOrder(1)

	// f2(x) = x
	f2 := New()
	f2.AddPoint(0, 0, 0)
	f2.AddPoint(5, 5, 0)
	f2.SetOrder(1)

	f1.Multiply(f2) // f1 becomes f1*f2 = 2x

	if !almostEqual(f1.F(2.5), 5.0) {
		t.Errorf("F(2.5) after multiply is %v, want 5.0", f1.F(2.5))
	}
	if !almostEqual(f1.F(5), 10) {
		t.Errorf("F(5) after multiply is %v, want 10.0", f1.F(5))
	}

	// Multiply functions with non-identical X points
	g1 := New()
	g1.SetOrder(1)
	g1.AddPoint(0, 2, 0)
	g1.AddPoint(2, 2, 0)

	g2 := New()
	g2.SetOrder(1)
	g2.AddPoint(1, 3, 0)
	g2.AddPoint(3, 3, 0)

	g1.Multiply(g2)
	// Combined domain should span union [0, 1, 2, 3]
	if g1.GetNdots() != 4 {
		t.Fatalf("Expected 4 points after disjoint Multiply, got %d", g1.GetNdots())
	}
	if !almostEqual(g1.F(1), 6.0) { // 2 * 3
		t.Errorf("F(1) after disjoint Multiply = %v; want 6.0", g1.F(1))
	}
}

func TestMultiplyByScalar(t *testing.T) {
	f := New()
	f.AddPoint(0, 10, 0)
	f.AddPoint(1, 20, 0)
	f.SetOrder(1) // Linear interpolation

	scalar := 2.0
	f.MultiplyByScalar(scalar)

	// Force update_spline to ensure internal min/max are updated
	_ = f.F(0.5)

	if !almostEqual(f.F(0), 20.0) {
		t.Errorf("F(0) after MultiplyByScalar = %v; want 20.0", f.F(0))
	}
	if !almostEqual(f.F(1), 40.0) {
		t.Errorf("F(1) after MultiplyByScalar = %v; want 40.0", f.F(1))
	}
	if !almostEqual(f.F(0.5), 30.0) { // (10*2 + 20*2)/2 = 30
		t.Errorf("F(0.5) after MultiplyByScalar = %v; want 30.0", f.F(0.5))
	}
	if !almostEqual(f.GetYmin(), 20.0) {
		t.Errorf("GetYmin() after MultiplyByScalar = %v; want 20.0", f.GetYmin())
	}
	if !almostEqual(f.GetYmax(), 40.0) {
		t.Errorf("GetYmax() after MultiplyByScalar = %v; want 40.0", f.GetYmax())
	}
}

func TestAssign(t *testing.T) {
	source := New()
	source.AddPoint(1, 10, 1)
	source.AddPoint(2, 20, 2)
	source.SetOrder(1)
	source.SetTrapolation(TrapolationLinear)

	dest := New()
	dest.Assign(source)

	_ = dest.GetXmin()

	if dest.Order != source.Order {
		t.Errorf("Order mismatch: dest=%d, source=%d", dest.Order, source.Order)
	}
	if dest.Trapolation != source.Trapolation {
		t.Errorf("Trapolation mismatch: dest=%v, source=%v", dest.Trapolation, source.Trapolation)
	}
	if dest.GetNdots() != source.GetNdots() {
		t.Fatalf("Points count mismatch: dest=%d, source=%d", dest.GetNdots(), source.GetNdots())
	}
	for i := range source.X {
		if !almostEqual(dest.X[i], source.X[i]) ||
			!almostEqual(dest.Y[i], source.Y[i]) ||
			dest.epoch[i] != source.epoch[i] {
			t.Errorf("Point mismatch at index %d", i)
		}
	}
	// Check internal state variables (will trigger update_spline in dest again, but values should be consistent)
	if !almostEqual(dest.GetXmin(), source.GetXmin()) {
		t.Errorf("Xmin mismatch: dest=%v, source=%v", dest.GetXmin(), source.GetXmin())
	}
	if !almostEqual(dest.GetXmax(), source.GetXmax()) {
		t.Errorf("Xmax mismatch: dest=%v, source=%v", dest.GetXmax(), source.GetXmax())
	}
	if !almostEqual(dest.GetYmin(), source.GetYmin()) {
		t.Errorf("Ymin mismatch: dest=%v, source=%v", dest.GetYmin(), source.GetYmin())
	}
	if !almostEqual(dest.GetYmax(), source.GetYmax()) {
		t.Errorf("Ymax mismatch: dest=%v, source=%v", dest.GetYmax(), source.GetYmax())
	}
	if !almostEqual(dest.GetStep(), source.GetStep()) {
		t.Errorf("Step mismatch: dest=%v, source=%v", dest.GetStep(), source.GetStep())
	}
}

func TestMerge(t *testing.T) {
	f1 := New()
	f1.AddPoint(0, 0, 0)
	f1.AddPoint(2, 20, 0)

	f2 := New()
	f2.AddPoint(1, 10, 0)
	f2.AddPoint(3, 30, 0)
	f2.AddPoint(2, 22, 1) // Point with same X as in f1, will be overwritten

	f1.Merge(f2)

	// Expected points: (0,0), (1,10), (2, 22), (3,30)
	if f1.GetNdots() != 4 {
		t.Fatalf("After merge, expected 4 points, got %d", f1.GetNdots())
	}

	type expectedPoint struct {
		x, y  float64
		epoch uint32
	}
	expectedPoints := []expectedPoint{
		{x: 0, y: 0, epoch: 0},
		{x: 1, y: 10, epoch: 0},
		{x: 2, y: 22, epoch: 1},
		{x: 3, y: 30, epoch: 0},
	}

	for i, ep := range expectedPoints {
		if i >= len(f1.X) {
			t.Fatalf("Missing point at index %d", i)
		}
		if !almostEqual(f1.X[i], ep.x) || !almostEqual(f1.Y[i], ep.y) || f1.epoch[i] != ep.epoch {
			t.Errorf("Point %d mismatch: got (%v, %v, %v), want (%v, %v, %v)", i, f1.X[i], f1.Y[i], f1.epoch[i], ep.x, ep.y, ep.epoch)
		}
	}

	// Merging into empty or merging empty
	empty := New()
	nonEmpty := New()
	nonEmpty.AddPoint(1, 5, 0)
	empty.Merge(nonEmpty)
	if empty.GetNdots() != 1 || !almostEqual(empty.F(1), 5.0) {
		t.Errorf("Merge into empty failed, dots=%d", empty.GetNdots())
	}
	nonEmpty.Merge(New())
	if nonEmpty.GetNdots() != 1 {
		t.Errorf("Merge empty into non-empty altered point count: %d", nonEmpty.GetNdots())
	}
}

func TestDrawPS(t *testing.T) {
	f := New()
	f.AddPoint(0, 0, 0)
	f.AddPoint(1, 1, 0)
	f.AddPoint(2, 0, 0)

	// Create a temporary file for the PS output
	tempFile, err := os.CreateTemp("", "test_draw_ps_*.ps")
	if err != nil {
		t.Fatalf("Failed to create temp file: %v", err)
	}
	defer os.Remove(tempFile.Name()) // Clean up the file after the test
	tempFile.Close()

	err = f.DrawPS(tempFile.Name())
	if err != nil {
		t.Errorf("DrawPS failed: %v", err)
	}

	// Optionally, read the file content to ensure it's not empty
	content, err := os.ReadFile(tempFile.Name())
	if err != nil || len(content) == 0 {
		t.Errorf("Generated PS file is empty or could not be read: %v", err)
	}
}

func TestExpand(t *testing.T) {
	f := New()
	// ymin = 0, ymax = 10
	// Points closer to 10 than 0 are those with Y > 5.
	f.AddPoint(0, 0, 0)   // Not closer
	f.AddPoint(1, 8, 0)   // Closer (Y=8 > 5)
	f.AddPoint(2, 2, 0)   // Not closer
	f.AddPoint(4, 9, 0)   // Closer (Y=9 > 5)
	f.AddPoint(10, 10, 0) // Closer (Y=10 > 5)

	f.Expand(1)

	if f.GetNdots() != 6 {
		t.Fatalf("Expected 6 points, got %d", f.GetNdots())
	}

	y := f.F(7)
	if !almostEqual(y, 9.5) {
		t.Errorf("Expected F(7) to be 9.5, got %v", y)
	}
}

func TestNormalizeIndices(t *testing.T) {
	f := New()
	f.AddPoint(1, 10, 0)
	f.AddPoint(2, 20, 0)
	f.AddPoint(3, 30, 0)

	f.indices[0] = 10
	f.indices[1] = 12
	f.indices[2] = 15
	f.nextIndex = 16

	f.NormaliseIndices()

	if f.indices[0] != 1 {
		t.Errorf("Expected point 0 index to be 1, got %d", f.indices[0])
	}
	if f.indices[1] != 3 {
		t.Errorf("Expected point 1 index to be 3, got %d", f.indices[1])
	}
	if f.indices[2] != 6 {
		t.Errorf("Expected point 2 index to be 6, got %d", f.indices[2])
	}
	if f.nextIndex != 7 {
		t.Errorf("Expected f.nextIndex generator to be updated to 7, got %d", f.nextIndex)
	}
}

func TestEmptyAndSinglePointSafeguards(t *testing.T) {
	// Empty function operations should not panic
	f := New()
	f.Normalise()
	minIdx, maxIdx := f.NormaliseIndices()
	if minIdx != 0 || maxIdx != 0 {
		t.Errorf("NormaliseIndices on empty = (%d, %d); want (0, 0)", minIdx, maxIdx)
	}
	f.Smooth()
	f.MorePoints()
	f.Expand(2)
	if !math.IsNaN(f.Integrate()) && f.Integrate() != 0 {
		t.Errorf("Integrate on empty = %v; want 0", f.Integrate())
	}

	// Single point function
	f.AddPoint(5, 42, 0)
	f.Smooth()     // Should be a no-op (< 3 points)
	f.MorePoints() // Should be a no-op (<= 1 point)
	if f.GetNdots() != 1 {
		t.Fatalf("Expected 1 point, got %d", f.GetNdots())
	}
	// Querying single point
	if !almostEqual(f.F(5), 42.0) {
		t.Errorf("Single point F(5) = %v; want 42.0", f.F(5))
	}
	if !almostEqual(f.F(0), 42.0) { // Boundary extrapolation returns the sole point
		t.Errorf("Single point extrapolation F(0) = %v; want 42.0", f.F(0))
	}
}
