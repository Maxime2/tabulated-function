package tabulatedfunction

import (
	"testing"
)

func Test_AddPoint(t *testing.T) {
	direct := New()
	direct.SetOrder(1)

	direct.AddPoint(2.0494472732002826, 0.1)
	direct.AddPoint(1.156694013301916, 0.5)
	direct.AddPoint(0.46530775203579466, 0.1)
	direct.AddPoint(-1.1237643368175254, 0.1)
	direct.AddPoint(2.5864746065598427, 0.5)
	t.Logf("%v\n", direct)
	t.Logf("x=2.0494472732002826; y=%v; expected 0.1\n", direct.F(2.0494472732002826))
	t.Logf("x=1.156694013301916; y=%v; expected 0.5\n", direct.F(1.156694013301916))
	t.Logf("x=-1.1237643368175254; y=%v; expected 0.1\n", direct.F(-1.1237643368175254))
	t.Logf("x=2.5864746065598427; y=%v; expected 0.5\n", direct.F(2.5864746065598427))
	t.Logf("x=-2; y=%v; expected 0.1\n", direct.F(-2))
	t.Logf("x=3; y=%v; expected 0.5\n", direct.F(3))
	t.Logf("x=2.4; y=%v; expected 0.361106059244108...\n", direct.F(2.4))

}
