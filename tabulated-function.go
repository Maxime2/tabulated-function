package tabulatedfunction

import (
	"cmp"
	"fmt"
	"math"
	"slices"
)

type TFPoint struct {
	b, c, d float64
	X, Y    float64
	Cnt     uint64
	epoch   uint32
}

type TabulatedFunction struct {
	ixmin, ixmax, iymin, iymax float64
	istep                      float64
	iOrder                     int
	changed                    bool
	//
	P []TFPoint
}

// Create
func New() *TabulatedFunction {
	return &TabulatedFunction{
		iOrder:  3,
		changed: false,
	}
}

// splinevalue
func (f *TabulatedFunction) F(xi float64) float64 {
	var k, l int
	var r float64

	l = len(f.P)
	if 1 > l {
		return math.NaN()
	}
	k, found := slices.BinarySearchFunc(f.P, TFPoint{X: xi}, func(a, b TFPoint) int {
		return cmp.Compare(a.X, b.X)
	})
	if found {
		return f.P[k].Y
	}
	if k == l {
		return f.P[l-1].Y
	}
	if k == 0 {
		return f.P[0].Y
	}
	switch f.iOrder {
	case 0:
		return f.P[k-1].Y
	case 1:
		return f.P[k-1].Y + (f.P[k].Y-f.P[k-1].Y)*(xi-f.P[k-1].X)/(f.P[k].X-f.P[k-1].X)
	}
	if f.changed {
		f.update_spline()
	}
	r = xi - f.P[k].X
	return f.P[k].Y + r*(f.P[k].b+r*(f.P[k].c+r*f.P[k].d))
}

func (f *TabulatedFunction) update_spline() {
	var i, j int
	var h, alpha, l, mu, z []float64
	var det, x1, x2, y1, y2 float64

	f.changed = false
	i = len(f.P)
	j = i - 1
	f.ixmin = f.P[0].X
	f.ixmax = f.ixmin
	f.iymin = f.P[0].Y
	f.iymax = f.iymin
	for i = 1; i <= j; i++ {
		if f.P[i].X < f.ixmin {
			f.ixmin = f.P[i].X
		} else if f.P[i].X > f.ixmax {
			f.ixmax = f.P[i].X
		}
		if f.P[i].Y < f.iymin {
			f.iymin = f.P[i].Y
		}
		if f.P[i].Y > f.iymax {
			f.iymax = f.P[i].Y
		}
	}
	if j > 0 {
		f.istep = f.P[1].X - f.P[0].X
		for i = 2; i <= j; i++ {
			if (f.P[i].X - f.P[i-1].X) < f.istep {
				f.istep = f.P[i].X - f.P[i-1].X
			}
		}
	} else {
		f.istep = 0
	}
	for i = 0; i <= j; i++ {
		f.P[i].b = 0
		f.P[i].c = 0
		f.P[i].d = 0
	}
	if f.iOrder == 0 {
		return
	}
	h = make([]float64, j+1)
	for i = 0; i <= j-1; i++ {
		h[i] = f.P[i+1].X - f.P[i].X
	}
	if f.iOrder == 1 {
		for i = 0; i <= j-1; i++ {
			f.P[i].b = (f.P[i+1].Y - f.P[i].Y) / h[i]
		}
		return
	}
	if f.iOrder == 2 {
		for i = 0; i <= j-2; i++ {
			x1 = f.P[i+1].X - f.P[i].X
			x2 = f.P[i+2].X - f.P[i].X
			y1 = f.P[i+1].Y - f.P[i].Y
			y2 = f.P[i+2].Y - f.P[i].Y
			det = x1 * x2 * (x2 - x1)
			det = 1 / det
			f.P[i].b = (y1*x2*x2 - y2*x1*x1) * det
			f.P[i].c = (y2*x1 - y1*x2) * det
		}
		f.P[j-1].b = (f.P[j].Y - f.P[j-1].Y) / h[i]
		return
	}
	alpha = make([]float64, j+1)
	l = make([]float64, j+1)
	mu = make([]float64, j+1)
	z = make([]float64, j+1)
	for i = 1; i <= j-1; i++ {
		alpha[i] = 3/h[i]*(f.P[i+1].Y-f.P[i].Y) - 3/h[i-1]*(f.P[i].Y-f.P[i-1].Y)
	}
	l[0] = 1
	mu[0] = 0
	z[0] = 0
	for i = 1; i <= j-1; i++ {
		l[i] = 2*(f.P[i+1].X-f.P[i-1].X) - h[i-1]*mu[i-1]
		mu[i] = h[i] / l[i]
		z[i] = (alpha[i] - h[i-1]*z[i-1]) / l[i]
	}
	l[j] = 1
	z[j] = 0
	f.P[j].c = 0
	for i = j - 1; i >= 0; i-- {
		f.P[i].c = z[i] - mu[i]*f.P[i+1].c
		f.P[i].b = (f.P[i+1].Y-f.P[i].Y)/h[i] - h[i]*(f.P[i+1].c+2*f.P[i].c)/3
		f.P[i].d = (f.P[i+1].c - f.P[i].c) / 3 / h[i]
	}
}

func (f *TabulatedFunction) SetOrder(new_value int) {
	f.iOrder = new_value
	f.changed = true
}

func (f *TabulatedFunction) AddPoint(Xn, Yn float64, epoch uint32, args ...uint64) float64 {
	var i, l int
	var cnt uint64 = 1
	if len(args) > 0 {
		cnt = args[0]
	}
	f.changed = true
	l = len(f.P)
	if l == 0 {
		f.P = append(f.P, TFPoint{X: Xn, Y: Yn, Cnt: cnt, epoch: epoch})
		return Yn
	}
	i, found := slices.BinarySearchFunc(f.P, TFPoint{X: Xn}, func(a, b TFPoint) int {
		return cmp.Compare(a.X, b.X)
	})
	if found {
		//f.P[i].X = Xn
		if f.P[i].epoch < epoch {
			f.P[i].Y = Yn
			f.P[i].Cnt = 1
			f.P[i].epoch = epoch
			return Yn
		}
		f.P[i].Y = (float64(f.P[i].Cnt)*f.P[i].Y + Yn)
		f.P[i].Cnt += cnt
		f.P[i].Y /= float64(f.P[i].Cnt)
		return f.P[i].Y
	}
	if i == l {
		f.P = append(f.P, TFPoint{X: Xn, Y: Yn, Cnt: cnt, epoch: epoch})
		return Yn
	}
	f.P = slices.Insert(f.P, i, TFPoint{X: Xn, Y: Yn, Cnt: cnt, epoch: epoch})
	return Yn
}

func (f *TabulatedFunction) LoadConstant(new_Y, new_xmin, new_xmax float64) {
	f.ixmin = new_xmin
	f.ixmax = new_xmax
	f.iymin = new_Y
	f.iymax = f.iymin
	f.P = append([]TFPoint{}, TFPoint{X: f.ixmin, Y: f.iymin, Cnt: 1, epoch: 0})
	f.istep = f.ixmax - f.ixmin
	f.changed = false
}

func (f *TabulatedFunction) Normalise() {
	var i int
	var ym float64 = math.Max(math.Abs(f.iymax), math.Abs(f.iymin))

	for i = range f.P {
		f.P[i].Y /= ym
	}
	f.changed = true
}

func (f *TabulatedFunction) Multiply(by *TabulatedFunction) {
	var i int
	var Yt []float64

	if f.changed {
		f.update_spline()
	}
	Yt = make([]float64, len(by.P))
	for i = range by.P {
		Yt[i] = by.P[i].Y * f.F(by.P[i].X)
	}

	for i = range f.P {
		f.P[i].Y *= by.F(f.P[i].X)
	}

	for i = range by.P {
		f.AddPoint(by.P[i].X, Yt[i], by.P[i].epoch, by.P[i].Cnt)
	}
	f.changed = true
}

func (f *TabulatedFunction) MultiplyByFloat64(by float64) {
	if f.changed {
		f.update_spline()
	}
	for i := range f.P {
		f.P[i].Y *= by
		f.P[i].b *= by
		f.P[i].c *= by
		f.P[i].d *= by
	}
	f.iymin *= by
	f.iymax *= by
}

func (f *TabulatedFunction) Assign(s *TabulatedFunction) {
	f.ixmin = s.ixmin
	f.ixmax = s.ixmax
	f.iymin = s.iymin
	f.iymax = s.iymax
	f.istep = s.istep
	f.iOrder = s.iOrder

	f.P = slices.Clone(s.P)
	f.changed = true
}

func (f *TabulatedFunction) Integrate() float64 {
	var i, l int
	var tmp, dif float64

	if f.changed {
		f.update_spline()
	}
	tmp = 0
	l = len(f.P) - 1
	for i = 0; i < l; i++ {
		dif = f.P[i+1].X - f.P[i].X
		tmp += dif * (f.P[i].Y + dif*(f.P[i].b/2+dif*(f.P[i].c/3+dif*f.P[i].d/4)))
	}
	return tmp
}

func (f *TabulatedFunction) Clear() {
	f.P = make([]TFPoint, 0)
	f.ixmin = 0
	f.ixmax = 0
	f.iymin = 0
	f.iymax = 0
	f.istep = 0
	f.changed = false
}

func (f *TabulatedFunction) MorePoints() {
	var i, j, k int
	var Yt []float64

	if f.changed {
		f.update_spline()
	}
	j = len(f.P) - 1
	if j <= 0 {
		return
	}
	i = 2 * j
	Yt = make([]float64, j)
	for k = 0; k < j; k++ {
		Yt[k] = f.F((f.P[k].X + f.P[k+1].X) / 2)
	}

	f.P = append(f.P, make([]TFPoint, j+1)...)

	for k = j; k >= 1; k-- {
		f.P[i] = f.P[k]
		i--
		f.P[i].X = (f.P[k].X + f.P[k-1].X) / 2
		f.P[i].Y = Yt[k-1]
		f.P[i].Cnt = 1
		f.P[i].epoch = f.P[k].epoch
		i--
	}
	f.changed = true
}

func (f *TabulatedFunction) Derivative() {
	var i int
	if f.changed {
		f.update_spline()
	}
	if f.iOrder > 0 {
		f.iOrder--
	}
	for i = range f.P {
		f.P[i].Y = f.P[i].b
		f.P[i].b = 2 * f.P[i].c
		f.P[i].c = 3 * f.P[i].d
		f.P[i].d = 0
	}
}

func (f *TabulatedFunction) Integral() {
	var i, j int
	var acc, prev_acc, r float64

	if f.changed {
		f.update_spline()
	}
	j = len(f.P) - 1
	acc = 0
	if f.iOrder < 3 {
		f.P[0].d = f.P[0].c / 3
		f.P[0].c = f.P[0].b / 2
		f.P[0].b = f.P[0].Y
		f.P[0].Y = 0
		for i = 1; i <= j; i++ {
			r = f.P[i].X - f.P[i-1].X
			f.P[i].d = f.P[i].c / 3
			f.P[i].c = f.P[i].b / 2
			f.P[i].b = f.P[i].Y
			f.P[i].Y = f.P[i-1].Y + r*(f.P[i-1].b+r*(f.P[i-1].c+r*f.P[i-1].d))
		}
		f.iOrder++
	} else {
		prev_acc = f.P[0].d / 4
		f.P[0].d = f.P[0].c / 3
		f.P[0].c = f.P[0].b / 2
		f.P[0].b = f.P[0].Y
		f.P[0].Y = 0
		for i = 1; i <= j; i++ {
			r = f.P[i].X - f.P[i-1].X
			acc = f.P[i].d / 4
			f.P[i].d = f.P[i].c / 3
			f.P[i].c = f.P[i].b / 2
			f.P[i].b = f.P[i].Y
			f.P[i].Y = f.P[i-1].Y + r*(f.P[i-1].b+r*(f.P[i].c+r*(f.P[i-1].d+r*prev_acc)))
			prev_acc = acc
		}
		f.changed = true
	}
}

func (f *TabulatedFunction) GetStep() float64 {
	if f.changed {
		f.update_spline()
	}
	return f.istep
}

func (f *TabulatedFunction) GetXmin() float64 {
	if f.changed {
		f.update_spline()
	}
	return f.ixmin
}

func (f *TabulatedFunction) GetXmax() float64 {
	if f.changed {
		f.update_spline()
	}
	return f.ixmax
}

func (f *TabulatedFunction) GetYmin() float64 {
	if f.changed {
		f.update_spline()
	}
	return f.iymin
}

func (f *TabulatedFunction) GetYmax() float64 {
	if f.changed {
		f.update_spline()
	}
	return f.iymax
}

func (f *TabulatedFunction) GetNdots() int {
	return len(f.P)
}

func (f *TabulatedFunction) String() string {
	s := "\nTabulated function:\n"
	s = fmt.Sprintf("%s\tiOrder: %v; changed: %v\n", s, f.iOrder, f.changed)
	s = fmt.Sprintf("%s\tixmin: %v; ixmax: %v\n", s, f.ixmin, f.ixmax)
	s = fmt.Sprintf("%s\tiymin: %v; iymax: %v\n", s, f.iymin, f.iymax)
	s = fmt.Sprintf("%s\tistep: %v\n", s, f.istep)
	s = fmt.Sprintf("%s\tPoints: %v\n", s, f.P)
	return s
}

func (f *TabulatedFunction) Epoch(epoch uint32) {
	slices.SortFunc(f.P, func(a, b TFPoint) int {
		if a.epoch >= epoch && b.epoch < epoch {
			return -1
		}
		return cmp.Compare(a.X, b.X)
	})
	i := slices.IndexFunc(f.P, func(p TFPoint) bool {
		return p.epoch < epoch
	})
	if i >= 0 {
		f.P = slices.Delete(f.P, i, len(f.P))
		f.changed = true
	}
}
