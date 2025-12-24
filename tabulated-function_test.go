package tabulatedfunction

import (
	"encoding/json"
	"math"
	"math/big"
	"testing"
)

const float64EqualityThreshold = 1e-9

func almostEqual(a, b float64) bool {
	if math.IsNaN(a) && math.IsNaN(b) {
		return true
	}
	return math.Abs(a-b) <= float64EqualityThreshold
}

func bf(f float64) *big.Float {
	return big.NewFloat(f)
}

func fl(b *big.Float) float64 {
	if b == nil {
		return math.NaN()
	}
	f, _ := b.Float64()
	return f
}

func TestNew(t *testing.T) {
	f := New()
	if f == nil {
		t.Fatal("New() returned nil")
	}
	if f.iOrder != 3 {
		t.Errorf("Expected default order to be 3, got %d", f.iOrder)
	}
	if f.trapolation != TrapolationSpline {
		t.Errorf("Expected default trapolation to be TrapolationSpline, got %v", f.trapolation)
	}
	if len(f.P) != 0 {
		t.Errorf("Expected new function to have 0 points, got %d", len(f.P))
	}
	if f.changed != false {
		t.Error("Expected new function to have changed=false")
	}
	if f.ixmin == nil {
		t.Error("ixmin is nil")
	}
}

func TestAddPointAndF(t *testing.T) {
	f := New()
	f.SetOrder(1) // Use linear for simplicity in this test

	// Test F on empty function
	if f.F(bf(0)) != nil {
		t.Errorf("F(0) on empty function should be nil, got %v", f.F(bf(0)))
	}

	// Add points
	f.AddPoint(bf(1), bf(10), 0)
	f.AddPoint(bf(3), bf(30), 0)
	f.AddPoint(bf(2), bf(20), 0) // Add out of order

	if len(f.P) != 3 {
		t.Fatalf("Expected 3 points, got %d", len(f.P))
	}

	// Check if sorted
	if !(f.P[0].X.Cmp(bf(1)) == 0 && f.P[1].X.Cmp(bf(2)) == 0 && f.P[2].X.Cmp(bf(3)) == 0) {
		t.Errorf("Points are not sorted correctly: %v", f.P)
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
			y := f.F(bf(tc.x))
			if !almostEqual(fl(y), tc.expected) {
				t.Errorf("F(%v) = %v; want %v", tc.x, y, tc.expected)
			}
		})
	}

	// Test adding an existing point (should average)
	// The original point is (2, 20). We add (2, 22).
	f.AddPoint(bf(2), bf(22), 1)
	expectedY := (20.0 + 22.0) / 2.0
	if !almostEqual(fl(f.F(bf(2))), expectedY) {
		t.Errorf("F(2) after adding existing point = %v; want %v", f.F(bf(2)), expectedY)
	}
}

func TestGetters(t *testing.T) {
	f := New()
	f.AddPoint(bf(0), bf(5), 0)
	f.AddPoint(bf(10), bf(-5), 0)
	f.AddPoint(bf(5), bf(15), 0)

	// Force update
	_ = f.F(bf(1))

	if !almostEqual(fl(f.GetXmin()), 0) {
		t.Errorf("GetXmin() = %v; want 0", f.GetXmin())
	}
	if !almostEqual(fl(f.GetXmax()), 10) {
		t.Errorf("GetXmax() = %v; want 10", f.GetXmax())
	}
	if !almostEqual(fl(f.GetYmin()), -5) {
		t.Errorf("GetYmin() = %v; want -5", f.GetYmin())
	}
	if !almostEqual(fl(f.GetYmax()), 15) {
		t.Errorf("GetYmax() = %v; want 15", f.GetYmax())
	}
	if f.GetNdots() != 3 {
		t.Errorf("GetNdots() = %v; want 3", f.GetNdots())
	}
	if !almostEqual(fl(f.GetStep()), 5) { // min step is between (0,x) and (5,x) or (5,x) and (10,x)
		t.Errorf("GetStep() = %v; want 5", f.GetStep())
	}
}

func TestOrders(t *testing.T) {
	f := New()
	f.AddPoint(bf(0), bf(0), 0)
	f.AddPoint(bf(1), bf(1), 0)
	f.AddPoint(bf(2), bf(0), 0)

	t.Run("Order 0 (Step function)", func(t *testing.T) {
		f.SetOrder(0)
		if !almostEqual(fl(f.F(bf(0.5))), 0) {
			t.Errorf("F(0.5) = %v; want 0", f.F(bf(0.5)))
		}
		if !almostEqual(fl(f.F(bf(1.5))), 1) {
			t.Errorf("F(1.5) = %v; want 1", f.F(bf(1.5)))
		}
		if !almostEqual(fl(f.F(bf(-1))), 0) { // Extrapolation
			t.Errorf("F(-1) = %v; want 0", f.F(bf(-1)))
		}
		// This test exposes a bug in right-side extrapolation for order 0.
		// It returns f.P[l-2].Y instead of f.P[l-1].Y.
		if !almostEqual(fl(f.F(bf(3))), 0) { // Extrapolation
			t.Errorf("F(3) = %v; want 0", f.F(bf(3)))
		}
	})

	t.Run("Order 1 (Linear)", func(t *testing.T) {
		f.SetOrder(1)
		if !almostEqual(fl(f.F(bf(0.5))), 0.5) {
			t.Errorf("F(0.5) = %v; want 0.5", f.F(bf(0.5)))
		}
		if !almostEqual(fl(f.F(bf(1.5))), 0.5) {
			t.Errorf("F(1.5) = %v; want 0.5", f.F(bf(1.5)))
		}
	})

	t.Run("Order 3 (Cubic Spline)", func(t *testing.T) {
		f.SetOrder(3)
		// With natural spline conditions, for these points, the spline should be symmetric around x=1.
		y1 := fl(f.F(bf(0.5)))
		y2 := fl(f.F(bf(1.5)))
		if !almostEqual(y1, y2) {
			t.Errorf("Expected symmetry for cubic spline, F(0.5)=%v, F(1.5)=%v", y1, y2)
		}
		// For these points, the spline in [0,1] is S(x) = -0.5x^3 + 1.5x.
		// So S(0.5) = -0.5*(0.125) + 1.5*0.5 = -0.0625 + 0.75 = 0.6875
		if !almostEqual(fl(f.F(bf(0.5))), 0.6875) {
			t.Errorf("F(0.5) = %v; want 0.6875", f.F(bf(0.5)))
		}
	})
}

// This test exposes a bug in the TrapolationLinear implementation, which fails to interpolate.
func TestTrapolationLinear(t *testing.T) {
	f := New()
	f.SetTrapolation(TrapolationLinear)
	f.AddPoint(bf(0), bf(0), 0)
	f.AddPoint(bf(2), bf(2), 0)

	// Test interpolation
	y := f.F(bf(1))
	if !almostEqual(fl(y), 1) {
		t.Errorf("Linear interpolation F(1) = %v; want 1", y)
	}
}

func TestTrapolationOpposite(t *testing.T) {
	f := New()
	f.SetTrapolation(TrapolationOpposite)
	f.AddPoint(bf(0), bf(10), 0)
	f.AddPoint(bf(10), bf(20), 0)

	// Test interpolation: returns the Y value of the opposite point.
	// Midpoint is (0+10)/2 = 5.
	// For xi <= 5, it should return the right point's Y value (20).
	// For xi > 5, it should return the left point's Y value (10).
	t.Run("Interpolation at midpoint", func(t *testing.T) {
		y := fl(f.F(bf(5)))
		expected := 20.0 // xi <= midpoint, so return right Y
		if !almostEqual(y, expected) {
			t.Errorf("Opposite interpolation F(5) = %v; want %v", y, expected)
		}
	})
	t.Run("Interpolation closer to right point", func(t *testing.T) {
		y := fl(f.F(bf(6)))
		expected := 10.0 // xi > midpoint, so return left Y
		if !almostEqual(y, expected) {
			t.Errorf("Opposite interpolation F(6) = %v; want %v", y, expected)
		}
	})
	t.Run("Interpolation closer to left point", func(t *testing.T) {
		y := fl(f.F(bf(4)))
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
		y := fl(f.F(bf(15)))
		expected := 10.0
		if !almostEqual(y, expected) {
			t.Errorf("Opposite extrapolation F(15) = %v; want %v", y, expected)
		}
	})
}

func TestLoadConstant(t *testing.T) {
	f := New()
	f.LoadConstant(bf(100), bf(-5), bf(5))
	if f.GetNdots() != 1 {
		t.Fatalf("LoadConstant should create 1 point, got %d", f.GetNdots())
	}
	if !almostEqual(fl(f.F(bf(0))), 100) {
		t.Errorf("F(0) = %v; want 100", f.F(bf(0)))
	}
	if !almostEqual(fl(f.F(bf(-10))), 100) { // Extrapolation
		t.Errorf("F(-10) = %v; want 100", f.F(bf(-10)))
	}
}

func TestClear(t *testing.T) {
	f := New()
	f.AddPoint(bf(1), bf(1), 0)
	f.Clear()
	if f.GetNdots() != 0 {
		t.Errorf("GetNdots() after Clear() = %d; want 0", f.GetNdots())
	}
	if f.F(bf(1)) != nil {
		t.Errorf("F(1) after Clear() should be nil, got %v", f.F(bf(1)))
	}
}

func TestJSON(t *testing.T) {
	f1 := New()
	f1.AddPoint(bf(0), bf(0), 1)
	f1.AddPoint(bf(1), bf(1), 2)
	f1.SetOrder(2)

	jsonData, err := json.Marshal(f1)
	if err != nil {
		t.Fatalf("json.Marshal failed: %v", err)
	}

	var f2 TabulatedFunction
	err = json.Unmarshal(jsonData, &f2)
	if err != nil {
		t.Fatalf("json.Unmarshal failed: %v", err)
	}

	if f1.iOrder != f2.iOrder {
		t.Errorf("Order mismatch: original=%d, unmarshaled=%d", f1.iOrder, f2.iOrder)
	}
	if len(f1.P) != len(f2.P) {
		t.Fatalf("Points count mismatch: original=%d, unmarshaled=%d", len(f1.P), len(f2.P))
	}
}

func TestDerivativeAndIntegral(t *testing.T) {
	// y = x^2
	f := New()
	f.AddPoint(bf(-2), bf(4), 0)
	f.AddPoint(bf(-1), bf(1), 0)
	f.AddPoint(bf(0), bf(0), 0)
	f.AddPoint(bf(1), bf(1), 0)
	f.AddPoint(bf(2), bf(4), 0)
	f.SetOrder(3)

	// Test definite integral: Integral of x^2 from -2 to 2 is 16/3
	integral := f.Integrate()
	// The natural cubic spline is an approximation. Its boundary conditions (second derivative is zero at endpoints)
	// do not match the true function y=x^2 (where the second derivative is 2 everywhere).
	// This causes a small, expected error. We relax the tolerance to account for this.
	if math.Abs(fl(integral)-16.0/3.0) > 1e-1 {
		t.Errorf("Spline integral of x^2 from -2 to 2 is %v, analytical is %v", integral, 16.0/3.0)
	}

	// Test Derivative: should be approx y' = 2x
	f.Derivative()
	if !almostEqual(fl(f.F(bf(0))), 0) {
		t.Errorf("Derivative at F(0) is %v, want ~0", f.F(bf(0)))
	}
	if math.Abs(fl(f.F(bf(1)))-2.0) > 0.5 { // Allow tolerance for spline approximation
		t.Errorf("Derivative at F(1) is %v, want ~2.0", f.F(bf(1)))
	}

	// Test Integral (indefinite): should be approx y = x^2 + C
	f.Integral()
	y_neg1 := fl(f.F(bf(-1)))
	y0 := fl(f.F(bf(0)))
	y1 := fl(f.F(bf(1)))
	// Check second difference to verify quadratic shape, independent of integration constant
	if math.Abs((y1-y0)-(y0-y_neg1)-2.0) > 0.5 {
		t.Errorf("Shape of integral is incorrect. Second difference at 0 is %v, want ~2.0", (y1-y0)-(y0-y_neg1))
	}
}

func TestMorePoints(t *testing.T) {
	f := New()
	f.SetOrder(1)
	f.AddPoint(bf(0), bf(0), 0)
	f.AddPoint(bf(2), bf(4), 0)

	f.MorePoints()

	if f.GetNdots() != 3 {
		t.Fatalf("MorePoints should have added 1 point, got %d total", f.GetNdots())
	}
	if !almostEqual(fl(f.P[1].X), 1.0) {
		t.Errorf("New point has X=%v, want 1.0", f.P[1].X)
	}
	if !almostEqual(fl(f.P[1].Y), 2.0) {
		t.Errorf("New point has Y=%v, want 2.0", f.P[1].Y)
	}
}

func TestEpoch(t *testing.T) {
	f := New()
	f.AddPoint(bf(0), bf(0), 0)
	f.AddPoint(bf(1), bf(1), 1)
	f.AddPoint(bf(2), bf(2), 2)
	f.AddPoint(bf(3), bf(3), 3)

	f.Epoch(2)

	if f.GetNdots() != 2 {
		t.Fatalf("After Epoch(2), expected 2 points, got %d", f.GetNdots())
	}
	if f.P[0].X.Cmp(bf(2)) != 0 || f.P[1].X.Cmp(bf(3)) != 0 {
		t.Errorf("Remaining points are incorrect: %v", f.P)
	}
}

func TestSmooth(t *testing.T) {
	f := New()
	f.AddPoint(bf(0), bf(0), 0)
	f.AddPoint(bf(1), bf(10), 0) // A noisy point
	f.AddPoint(bf(2), bf(4), 0)

	y_before := fl(f.P[1].Y)
	f.Smooth()
	y_after := fl(f.P[1].Y)

	if almostEqual(y_before, y_after) {
		t.Error("Smooth() did not change the Y value of the middle point")
	}
	// Based on manual calculation of the implemented formula:
	expected_y := (2.0 + 2.0) / 2.0
	if !almostEqual(y_after, expected_y) {
		t.Errorf("Smooth() produced %v, want %v", y_after, expected_y)
	}
}

func TestMultiply(t *testing.T) {
	// f1(x) = 2
	f1 := New()
	f1.AddPoint(bf(0), bf(2), 0)
	f1.AddPoint(bf(5), bf(2), 0)
	f1.SetOrder(1)

	// f2(x) = x
	f2 := New()
	f2.AddPoint(bf(0), bf(0), 0)
	f2.AddPoint(bf(5), bf(5), 0)
	f2.SetOrder(1)

	f1.Multiply(f2) // f1 becomes f1*f2 = 2x

	if !almostEqual(fl(f1.F(bf(2.5))), 5.0) {
		t.Errorf("F(2.5) after multiply is %v, want 5.0", f1.F(bf(2.5)))
	}
	if !almostEqual(fl(f1.F(bf(5))), 10) {
		t.Errorf("F(5) after multiply is %v, want 10.0", f1.F(bf(5)))
	}
}
