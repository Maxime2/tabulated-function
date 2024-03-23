package tabulatedfunction

import (
	"testing"
)

func TestAddPoint(t *testing.T) {
	direct := New()
	direct.SetOrder(1)

	direct.AddPoint(-0.3388588053705943, 0.1)
	t.Logf("f.X: %v\n", direct.X)
	t.Logf("f.Y: %v\n\n", direct.Y)
	direct.AddPoint(1.319795729531707, 0.1)
	t.Logf("f.X: %v\n", direct.X)
	t.Logf("f.Y: %v\n\n", direct.Y)
	direct.AddPoint(0.3850489958660157, 0.5)
	t.Logf("f.X: %v\n", direct.X)
	t.Logf("f.Y: %v\n\n", direct.Y)
	y := direct.F(-1.0666666666666695)
	t.Logf("  y = %v", y)
	y = direct.F(0.3850489958660157)
	t.Logf("  y = %v\n", y)
	t.Logf("%v\n", direct)
}
