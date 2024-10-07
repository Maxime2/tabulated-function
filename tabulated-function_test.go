package tabulatedfunction

import (
	"testing"
)

func TestAddPoint(t *testing.T) {
	direct := New()
	direct.SetOrder(1)

	direct.AddPoint(-10.3850489958660157, 0.1, 0)
	direct.AddPoint(0.3850489958660157, 0.5, 1)
	t.Logf("f.P: %v\n", direct.P)
	direct.AddPoint(-0.3388588053705943, 0.1, 2)
	t.Logf("f.P: %v\n", direct.P)
	direct.AddPoint(1.319795729531707, 0.1, 2)
	t.Logf("f.P: %v\n", direct.P)
	direct.DrawPS("test.ps")
	direct.AddPoint(0.3850489958660157, 0.1, 1)
	t.Logf("f.P: %v\n", direct.P)
	y := direct.F(-1.0666666666666695)
	t.Logf("  y = %v; expected: 0.1\n", y)
	y = direct.F(0.3850489958660157)
	t.Logf("  y = %v; expected: 0.3\n", y)
	t.Logf("  number of dots: %v; expected: 4\n", direct.GetNdots())
	t.Logf("%v\n", direct)
	//	direct.DrawPS("test.ps")
	direct.Epoch(2)
	t.Logf("  number of dots: %v; expected: 2\n", direct.GetNdots())
	t.Logf("%v\n", direct)
}
