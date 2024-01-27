package tabulatedfunction

import (
	"testing"
)

func AddPoint_test(t *testing.T) {
	direct := New()
	direct.SetOrder(1)

	direct.AddPoint(2.0494472732002826, 0.1)
	direct.AddPoint(1.156694013301916, 0.5)
	direct.AddPoint(0.46530775203579466, 0.1)

}
