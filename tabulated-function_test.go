package tabulatedfunction

import (
	"testing"
)

func TestAddPoint(t *testing.T) {
	direct := New()
	direct.SetOrder(1)

	direct.AddPoint(-0.3388588053705943, 0.1)
	t.Logf("f.X: %v\n", direct.X)
	t.Logf("f.Y: %v\n", direct.Y)
	t.Logf("f.Cnt: %v\n\n", direct.Cnt)
	direct.AddPoint(1.319795729531707, 0.1)
	t.Logf("f.X: %v\n", direct.X)
	t.Logf("f.Y: %v\n", direct.Y)
	t.Logf("f.Cnt: %v\n\n", direct.Cnt)
	direct.AddPoint(0.3850489958660157, 0.5)
	t.Logf("f.X: %v\n", direct.X)
	t.Logf("f.Y: %v\n", direct.Y)
	t.Logf("f.Cnt: %v\n\n", direct.Cnt)
	direct.AddPoint(0.3850489958660157, 0.1)
	t.Logf("f.X: %v\n", direct.X)
	t.Logf("f.Y: %v\n", direct.Y)
	t.Logf("f.Cnt: %v\n\n", direct.Cnt)
	y := direct.F(-1.0666666666666695)
	t.Logf("  y = %v; expected: 0.1\n", y)
	y = direct.F(0.3850489958660157)
	t.Logf("  y = %v; expected: 0.3\n", y)
	t.Logf("  number of dots: %v; expected: 3\n", direct.GetNdots())
	t.Logf("%v\n", direct)
}
