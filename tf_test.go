package tabulatedfunction

import (
	"math/big"
	"testing"
)

func Test_AddPoint(t *testing.T) {
	direct := New()
	direct.SetOrder(1)

	direct.AddPoint(big.NewFloat(2.0494472732002826), big.NewFloat(0.1), 1)
	direct.AddPoint(big.NewFloat(1.156694013301916), big.NewFloat(0.5), 1)
	direct.AddPoint(big.NewFloat(0.46530775203579466), big.NewFloat(0.1), 1)
	direct.AddPoint(big.NewFloat(-1.1237643368175254), big.NewFloat(0.1), 1)
	direct.AddPoint(big.NewFloat(2.5864746065598427), big.NewFloat(0.5), 1)
	t.Logf("%v\n", direct)
	t.Logf("x=2.0494472732002826; y=%v; expected 0.1\n", direct.F(big.NewFloat(2.0494472732002826)))
	t.Logf("x=1.156694013301916; y=%v; expected 0.5\n", direct.F(big.NewFloat(1.156694013301916)))
	t.Logf("x=-1.1237643368175254; y=%v; expected 0.1\n", direct.F(big.NewFloat(-1.1237643368175254)))
	t.Logf("x=2.5864746065598427; y=%v; expected 0.5\n", direct.F(big.NewFloat(2.5864746065598427)))
	t.Logf("x=-2; y=%v; expected 0.1\n", direct.F(big.NewFloat(-2)))
	t.Logf("x=3; y=%v; expected 0.5\n", direct.F(big.NewFloat(3)))
	t.Logf("x=2.4; y=%v; expected 0.361106059244108...\n", direct.F(big.NewFloat(2.4)))

}
