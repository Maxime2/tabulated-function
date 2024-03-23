package tabulatedfunction

import (
	"fmt"
	"math"
)

type TabulatedFunction struct {
	ixmin, ixmax, iymin, iymax float64
	b, c, d                    []float64
	istep                      float64
	iOrder                     int
	changed                    bool
	//
	X, Y []float64
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
	var i, j, k int
	var r float64

	i = 0
	j = len(f.X) - 1
	if j < i {
		return math.NaN()
	}
	if f.changed {
		f.update_spline()
	}
	if xi > f.ixmax {
		return f.Y[j]
	}
	if xi < f.ixmin {
		return f.Y[0]
	}
	k = 0
	for j >= i {
		k = (i + j) >> 1
		if xi < f.X[k] {
			j = k - 1
		} else if xi > f.X[k] {
			i = k + 1
		} else {
			goto found
		}
	}
	if xi < f.X[k] {
		k--
	}
	if k < 0 {
		return f.Y[0]
	}
found:
	r = xi - f.X[k]
	return f.Y[k] + r*(f.b[k]+r*(f.c[k]+r*f.d[k]))
}

func (f *TabulatedFunction) update_spline() {
	var i, j int
	var h, alpha, l, mu, z []float64
	var det, x1, x2, y1, y2 float64

	f.changed = false
	i = len(f.X)
	j = i - 1
	f.b = make([]float64, i)
	f.c = make([]float64, i)
	f.d = make([]float64, i)
	f.ixmin = f.X[0]
	f.ixmax = f.ixmin
	f.iymin = f.Y[0]
	f.iymax = f.iymin
	for i = 1; i <= j; i++ {
		if f.X[i] < f.ixmin {
			f.ixmin = f.X[i]
		}
		if f.X[i] > f.ixmax {
			f.ixmax = f.X[i]
		}
		if f.Y[i] < f.iymin {
			f.iymin = f.Y[i]
		}
		if f.Y[i] > f.iymax {
			f.iymax = f.Y[i]
		}
	}
	if j > 0 {
		f.istep = f.X[1] - f.X[0]
		for i = 2; i <= j; i++ {
			if (f.X[i] - f.X[i-1]) < f.istep {
				f.istep = f.X[i] - f.X[i-1]
			}
		}
	} else {
		f.istep = 0
	}
	for i = 0; i <= j; i++ {
		f.b[i] = 0
		f.c[i] = 0
		f.d[i] = 0
	}
	if f.iOrder == 0 {
		return
	}
	h = make([]float64, j+1)
	for i = 0; i <= j-1; i++ {
		h[i] = f.X[i+1] - f.X[i]
	}
	if f.iOrder == 1 {
		for i = 0; i <= j-1; i++ {
			f.b[i] = (f.Y[i+1] - f.Y[i]) / h[i]
		}
		return
	}
	if f.iOrder == 2 {
		for i = 0; i <= j-2; i++ {
			x1 = f.X[i+1] - f.X[i]
			x2 = f.X[i+2] - f.X[i]
			y1 = f.Y[i+1] - f.Y[i]
			y2 = f.Y[i+2] - f.Y[i]
			det = x1 * x2 * (x2 - x1)
			det = 1 / det
			f.b[i] = (y1*x2*x2 - y2*x1*x1) * det
			f.c[i] = (y2*x1 - y1*x2) * det
		}
		f.b[j-1] = (f.Y[j] - f.Y[j-1]) / h[i]
		return
	}
	alpha = make([]float64, j+1)
	l = make([]float64, j+1)
	mu = make([]float64, j+1)
	z = make([]float64, j+1)
	for i = 1; i <= j-1; i++ {
		alpha[i] = 3/h[i]*(f.Y[i+1]-f.Y[i]) - 3/h[i-1]*(f.Y[i]-f.Y[i-1])
	}
	l[0] = 1
	mu[0] = 0
	z[0] = 0
	for i = 1; i <= j-1; i++ {
		l[i] = 2*(f.X[i+1]-f.X[i-1]) - h[i-1]*mu[i-1]
		mu[i] = h[i] / l[i]
		z[i] = (alpha[i] - h[i-1]*z[i-1]) / l[i]
	}
	l[j] = 1
	z[j] = 0
	f.c[j] = 0
	for i = j - 1; i >= 0; i-- {
		f.c[i] = z[i] - mu[i]*f.c[i+1]
		f.b[i] = (f.Y[i+1]-f.Y[i])/h[i] - h[i]*(f.c[i+1]+2*f.c[i])/3
		f.d[i] = (f.c[i+1] - f.c[i]) / 3 / h[i]
	}
}

func (f *TabulatedFunction) SetOrder(new_value int) {
	f.iOrder = new_value
	f.changed = true
}

func (f *TabulatedFunction) AddPoint(Xn, Yn float64) {
	var i, k, l int
	f.changed = true
	l = len(f.X)
	if l == 0 {
		f.X = append(f.X, Xn)
		f.Y = append(f.Y, Yn)
		return
	}
	for i = 0; i < l && f.X[i] < Xn; i++ {
	}
	if i == l {
		f.X = append(f.X, Xn)
		f.Y = append(f.Y, Yn)
		return
	}
	if f.X[i] == Xn {
		f.X[i] = Xn
		f.Y[i] = (f.Y[i] + Yn) / 2
		return
	}
	k = l - 1
	f.X = append(f.X, f.X[k])
	copy(f.X[i+1:], f.X[i:k])
	f.X[i] = Xn
	f.Y = append(f.Y, f.Y[k])
	copy(f.Y[i+1:], f.Y[i:k])
	f.Y[i] = Yn
}

func (f *TabulatedFunction) LoadConstant(new_Y, new_xmin, new_xmax float64) {
	f.X = make([]float64, 1)
	f.Y = make([]float64, 1)
	f.b = make([]float64, 1)
	f.c = make([]float64, 1)
	f.d = make([]float64, 1)
	f.ixmin = new_xmin
	f.ixmax = new_xmax
	f.iymin = new_Y
	f.iymax = f.iymin
	f.X[0] = f.ixmin
	f.Y[0] = f.iymin
	f.b[0] = 0
	f.c[0] = 0
	f.d[0] = 0
	f.istep = f.ixmax - f.ixmin
	f.changed = false
}

func (f *TabulatedFunction) Normalise() {
	var i int
	var ym float64
	ym = math.Max(math.Abs(f.iymax), math.Abs(f.iymin))

	for i = range f.Y {
		f.Y[i] /= ym
	}
	f.changed = true
}

func (f *TabulatedFunction) Multiply(by *TabulatedFunction) {
	var i int
	var Yt []float64

	if f.changed {
		f.update_spline()
	}
	Yt = make([]float64, len(by.X))
	for i = range by.X {
		Yt[i] = by.Y[i] * f.F(by.X[i])
	}

	for i = range f.X {
		f.Y[i] *= by.F(f.X[i])
	}

	for i = range by.X {
		f.AddPoint(by.X[i], Yt[i])
	}
	f.changed = true
}

func (f *TabulatedFunction) MultiplyByFloat64(by float64) {
	if f.changed {
		f.update_spline()
	}
	for i := range f.Y {
		f.Y[i] *= by
		f.b[i] *= by
		f.c[i] *= by
		f.d[i] *= by
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

	f.X = append([]float64{}, s.X...)
	f.Y = append([]float64{}, s.Y...)
	f.b = append([]float64{}, s.b...)
	f.c = append([]float64{}, s.c...)
	f.d = append([]float64{}, s.d...)
	f.changed = true
}

func (f *TabulatedFunction) Integrate() float64 {
	var i, l int
	var tmp, dif float64

	if f.changed {
		f.update_spline()
	}
	tmp = 0
	l = len(f.X) - 1
	for i = 0; i < l; i++ {
		dif = f.X[i+1] - f.X[i]
		tmp += dif * (f.Y[i] + dif*(f.b[i]/2+dif*(f.c[i]/3+dif*f.d[i]/4)))
	}
	return tmp
}

func (f *TabulatedFunction) Clear() {
	f.X = make([]float64, 0)
	f.Y = make([]float64, 0)
	f.b = make([]float64, 0)
	f.c = make([]float64, 0)
	f.d = make([]float64, 0)
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
	j = len(f.X) - 1
	if j <= 0 {
		return
	}
	i = 2 * j
	Yt = make([]float64, j)
	for k = 0; k < j; k++ {
		Yt[k] = f.F((f.X[k] + f.X[k+1]) / 2)
	}

	f.X = append(f.X, make([]float64, j+1)...)
	f.Y = append(f.Y, make([]float64, j+1)...)

	for k = j; k >= 1; k-- {
		f.X[i] = f.X[k]
		f.Y[i] = f.Y[k]
		i--
		f.X[i] = (f.X[k] + f.X[k-1]) / 2
		f.Y[i] = Yt[k-1]
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
	for i = range f.X {
		f.Y[i] = f.b[i]
		f.b[i] = 2 * f.c[i]
		f.c[i] = 3 * f.d[i]
		f.d[i] = 0
	}
}

func (f *TabulatedFunction) Integral() {
	var i, j int
	var acc, prev_acc, r float64

	if f.changed {
		f.update_spline()
	}
	j = len(f.X) - 1
	acc = 0
	if f.iOrder < 3 {
		f.d[0] = f.c[0] / 3
		f.c[0] = f.b[0] / 2
		f.b[0] = f.Y[0]
		f.Y[0] = 0
		for i = 1; i <= j; i++ {
			r = f.X[i] - f.X[i-1]
			f.d[i] = f.c[i] / 3
			f.c[i] = f.b[i] / 2
			f.b[i] = f.Y[i]
			f.Y[i] = f.Y[i-1] + r*(f.b[i-1]+r*(f.c[i-1]+r*f.d[i-1]))
		}
		f.iOrder++
	} else {
		prev_acc = f.d[0] / 4
		f.d[0] = f.c[0] / 3
		f.c[0] = f.b[0] / 2
		f.b[0] = f.Y[0]
		f.Y[0] = 0
		for i = 1; i <= j; i++ {
			r = f.X[i] - f.X[i-1]
			acc = f.d[i] / 4
			f.d[i] = f.c[i] / 3
			f.c[i] = f.b[i] / 2
			f.b[i] = f.Y[i]
			f.Y[i] = f.Y[i-1] + r*(f.b[i-1]+r*(f.c[i-1]+r*(f.d[i-1]+r*prev_acc)))
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

func (f *TabulatedFunction) String() string {
	s := "\nTabulated function:\n"
	s = fmt.Sprintf("%s\tiOrder: %v; changed: %v\n", s, f.iOrder, f.changed)
	s = fmt.Sprintf("%s\tixmin: %v; ixmax: %v\n", s, f.ixmin, f.ixmax)
	s = fmt.Sprintf("%s\tiymin: %v; iymax: %v\n", s, f.iymin, f.iymax)
	s = fmt.Sprintf("%s\tistep: %v\n", s, f.istep)
	s = fmt.Sprintf("%s\tb: %v\n", s, f.b)
	s = fmt.Sprintf("%s\tc: %v\n", s, f.c)
	s = fmt.Sprintf("%s\td: %v\n", s, f.d)
	s = fmt.Sprintf("%s\tX: %v\n", s, f.X);
	s = fmt.Sprintf("%s\tY: %v\n", s, f.Y);
	return s
}
