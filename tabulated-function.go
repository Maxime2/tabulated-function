package tabulatedfunction

import (
	"fmt"
	"math/big"
	"os"
	"slices"
)

type Trapolation int

const (
	TrapolationLinear   Trapolation = 0
	TrapolationSpline   Trapolation = 1
	TrapolationShift    Trapolation = 2
	TrapolationMinMax   Trapolation = 3
	TrapolationOpposite Trapolation = 4
)

type TFPoint struct {
	b, c, d *big.Float
	X, Y    *big.Float
	epoch   uint32
}

type TabulatedFunction struct {
	ixmin, ixmax, iymin, iymax *big.Float
	istep                      *big.Float
	iOrder                     int
	changed                    bool
	trapolation                Trapolation
	//
	P []TFPoint
}

// Create
func New() *TabulatedFunction {
	return &TabulatedFunction{
		iOrder:      3,
		trapolation: TrapolationSpline,
		changed:     false,
		ixmin:       new(big.Float),
		ixmax:       new(big.Float),
		iymin:       new(big.Float),
		iymax:       new(big.Float),
		istep:       new(big.Float),
	}
}

// splinevalue
func (f *TabulatedFunction) F(xi *big.Float) *big.Float {
	l := len(f.P)
	if l == 0 {
		return nil
	}
	k, found := slices.BinarySearchFunc(f.P, TFPoint{X: xi}, func(a, b TFPoint) int {
		return a.X.Cmp(b.X)
	})
	if found {
		return new(big.Float).Set(f.P[k].Y)
	}

	var left, right int
	if k >= l {
		right = l - 1
	} else {
		right = k
	}
	if k == 0 {
		left = 0
	} else {
		left = k - 1
	}

	return f._interpolate(xi, left, right, f.trapolation)
}

func (f *TabulatedFunction) Trapolate(xi *big.Float, trapolation Trapolation) *big.Float {
	l := len(f.P)
	if l == 0 {
		return big.NewFloat(1)
	}
	k, found := slices.BinarySearchFunc(f.P, TFPoint{X: xi}, func(a, b TFPoint) int {
		return a.X.Cmp(b.X)
	})

	var left, right int
	if found {
		if k == l-1 {
			right = k
		} else {
			right = k + 1
		}
		if k == 0 {
			left = k
		} else {
			left = k - 1
		}
	} else {
		if k >= l {
			right = l - 1
		} else {
			right = k
		}
		left = k - 1
	}

	return f._interpolate(xi, left, right, trapolation)
}

// _interpolate calculates the value at xi based on the surrounding points.
// It assumes left and right indices are correctly set.
func (f *TabulatedFunction) _interpolate(xi *big.Float, left, right int, trapolation Trapolation) *big.Float {
	switch trapolation {
	case TrapolationLinear:
		// If left==right, it's extrapolation. Return the boundary value.
		if left == right || left < 0 {
			return new(big.Float).Set(f.P[right].Y)
		}
		// Avoid division by zero if points are not distinct on X.
		dx := new(big.Float).Sub(f.P[right].X, f.P[left].X)
		if dx.Sign() == 0 {
			return new(big.Float).Set(f.P[left].Y)
		}
		// f.P[left].Y + (f.P[right].Y-f.P[left].Y)*(xi-f.P[left].X)/dx
		dy := new(big.Float).Sub(f.P[right].Y, f.P[left].Y)
		ratio := new(big.Float).Sub(xi, f.P[left].X)
		ratio.Quo(ratio, dx)
		res := new(big.Float).Mul(dy, ratio)
		res.Add(f.P[left].Y, res)
		return res

		// Similar to Linear then switch to opposite value
	case TrapolationOpposite:
		// If left==right, it's extrapolation. Return the boundary value.
		if left == right || left < 0 {
			// f.GetYmin() + f.GetYmax() - f.P[right].Y
			res := new(big.Float).Add(f.GetYmin(), f.GetYmax())
			res.Sub(res, f.P[right].Y)
			return res
		}

		// Determine if xi is closer to the left or right point.
		// midpoint := (f.P[left].X + f.P[right].X) / 2
		midpoint := new(big.Float).Add(f.P[left].X, f.P[right].X)
		midpoint.Quo(midpoint, big.NewFloat(2))

		if xi.Cmp(midpoint) <= 0 {
			// If closer to the left point, return the right point's Y value.
			return new(big.Float).Set(f.P[right].Y)
		}
		// Otherwise (closer to the right point), return the left point's Y value.
		return new(big.Float).Set(f.P[left].Y)

	case TrapolationSpline:
		// For order 1 (linear), extrapolation should clamp to the boundary value.
		// This makes its extrapolation behavior consistent with TrapolationLinear
		// and fixes an inconsistency where left-side extrapolated linearly while
		// the right-side clamped.
		if f.iOrder == 1 && (left == right || left < 0) {
			return new(big.Float).Set(f.P[right].Y)
		}

		if f.changed {
			f.update_spline()
		}

		// For interpolation, the segment is defined by P[left].
		// For extrapolation, left==right, and we use the coefficients of that boundary point.
		// r := xi - f.P[left].X
		// return f.P[left].Y + r*(f.P[left].b+r*(f.P[left].c+r*f.P[left].d))
		r := new(big.Float).Sub(xi, f.P[left].X)
		term := new(big.Float).Mul(r, f.P[left].d)
		term.Add(f.P[left].c, term)
		term.Mul(r, term)
		term.Add(f.P[left].b, term)
		term.Mul(r, term)
		term.Add(f.P[left].Y, term)
		return term

	case TrapolationShift:
		if left < 0 {
			return new(big.Float).Set(f.P[right].Y)
		}
		res := new(big.Float).Add(f.P[left].Y, f.P[right].Y)
		return res.Quo(res, big.NewFloat(2))

	case TrapolationMinMax:
		var v [2]*big.Float

		if left < 0 {
			left = 0
		}

		if f.changed {
			f.update_spline()
		}

		v[0] = new(big.Float)
		v[1] = new(big.Float)

		diff := new(big.Float).Sub(f.iymin, f.P[left].Y)
		v[0].Add(v[0], diff.Abs(diff))
		diff.Sub(f.iymin, f.P[right].Y)
		v[0].Add(v[0], diff.Abs(diff))

		diff.Sub(f.iymax, f.P[left].Y)
		v[1].Add(v[1], diff.Abs(diff))
		diff.Sub(f.iymax, f.P[right].Y)
		v[1].Add(v[1], diff.Abs(diff))

		avg := new(big.Float).Add(f.P[left].Y, f.P[right].Y)
		avg.Quo(avg, big.NewFloat(2))

		if v[0].Cmp(v[1]) > 0 {
			return avg.Add(f.iymin, avg).Quo(avg, big.NewFloat(2))
		}
		return avg.Add(f.iymax, avg).Quo(avg, big.NewFloat(2))

	}
	// This is unreachable if all Trapolation values are handled.
	// A panic is better than returning a magic number.
	panic("unhandled trapolation type")
}

func (f *TabulatedFunction) update_spline() {
	var i, j int
	var h, alpha, l, mu, z []*big.Float
	var det, x1, x2, y1, y2 *big.Float

	f.changed = false
	i = len(f.P)
	if i == 0 {
		return
	}
	j = i - 1
	f.ixmin.Set(f.P[0].X)
	f.ixmax.Set(f.ixmin)
	f.iymin.Set(f.P[0].Y)
	f.iymax.Set(f.iymin)
	for i = 1; i <= j; i++ {
		if f.P[i].X.Cmp(f.ixmin) < 0 {
			f.ixmin.Set(f.P[i].X)
		} else if f.P[i].X.Cmp(f.ixmax) > 0 {
			f.ixmax.Set(f.P[i].X)
		}
		if f.P[i].Y.Cmp(f.iymin) < 0 {
			f.iymin.Set(f.P[i].Y)
		}
		if f.P[i].Y.Cmp(f.iymax) > 0 {
			f.iymax.Set(f.P[i].Y)
		}
	}
	if j > 0 {
		f.istep.Sub(f.P[1].X, f.P[0].X)
		for i = 2; i <= j; i++ {
			diff := new(big.Float).Sub(f.P[i].X, f.P[i-1].X)
			if diff.Cmp(f.istep) < 0 {
				f.istep.Set(diff)
			}
		}
	} else {
		f.istep.SetInt64(0)
	}
	for i = 0; i <= j; i++ {
		f.P[i].b = new(big.Float)
		f.P[i].c = new(big.Float)
		f.P[i].d = new(big.Float)
	}
	if f.iOrder == 0 {
		return
	}
	h = make([]*big.Float, j+1)
	for i = 0; i <= j-1; i++ {
		h[i] = new(big.Float).Sub(f.P[i+1].X, f.P[i].X)
	}
	if f.iOrder == 1 {
		for i = 0; i <= j-1; i++ {
			f.P[i].b.Sub(f.P[i+1].Y, f.P[i].Y)
			f.P[i].b.Quo(f.P[i].b, h[i])
		}
		return
	}
	if f.iOrder == 2 {
		for i = 0; i <= j-2; i++ {
			x1 = new(big.Float).Sub(f.P[i+1].X, f.P[i].X)
			x2 = new(big.Float).Sub(f.P[i+2].X, f.P[i].X)
			y1 = new(big.Float).Sub(f.P[i+1].Y, f.P[i].Y)
			y2 = new(big.Float).Sub(f.P[i+2].Y, f.P[i].Y)

			// det = x1 * x2 * (x2 - x1)
			det = new(big.Float).Mul(x1, x2)
			det.Mul(det, new(big.Float).Sub(x2, x1))
			// det = 1 / det
			det.Quo(big.NewFloat(1), det)

			// b = (y1*x2*x2 - y2*x1*x1) * det
			t1 := new(big.Float).Mul(y1, x2)
			t1.Mul(t1, x2)
			t2 := new(big.Float).Mul(y2, x1)
			t2.Mul(t2, x1)
			f.P[i].b.Sub(t1, t2).Mul(f.P[i].b, det)

			// c = (y2*x1 - y1*x2) * det
			t1.Mul(y2, x1)
			t2.Mul(y1, x2)
			f.P[i].c.Sub(t1, t2).Mul(f.P[i].c, det)
		}
		f.P[j-1].b.Sub(f.P[j].Y, f.P[j-1].Y).Quo(f.P[j-1].b, h[j-1])
		return
	}
	alpha = make([]*big.Float, j+1)
	l = make([]*big.Float, j+1)
	mu = make([]*big.Float, j+1)
	z = make([]*big.Float, j+1)
	for i = 1; i <= j-1; i++ {
		// alpha[i] = 3/h[i]*(f.P[i+1].Y-f.P[i].Y) - 3/h[i-1]*(f.P[i].Y-f.P[i-1].Y)
		t1 := new(big.Float).Sub(f.P[i+1].Y, f.P[i].Y)
		t1.Mul(t1, big.NewFloat(3)).Quo(t1, h[i])
		t2 := new(big.Float).Sub(f.P[i].Y, f.P[i-1].Y)
		t2.Mul(t2, big.NewFloat(3)).Quo(t2, h[i-1])
		alpha[i] = new(big.Float).Sub(t1, t2)
	}
	l[0] = big.NewFloat(1)
	mu[0] = new(big.Float)
	z[0] = new(big.Float)
	for i = 1; i <= j-1; i++ {
		// l[i] = 2*(f.P[i+1].X-f.P[i-1].X) - h[i-1]*mu[i-1]
		t1 := new(big.Float).Sub(f.P[i+1].X, f.P[i-1].X)
		t1.Mul(t1, big.NewFloat(2))
		t2 := new(big.Float).Mul(h[i-1], mu[i-1])
		l[i] = new(big.Float).Sub(t1, t2)
		// mu[i] = h[i] / l[i]
		mu[i] = new(big.Float).Quo(h[i], l[i])
		// z[i] = (alpha[i] - h[i-1]*z[i-1]) / l[i]
		t3 := new(big.Float).Mul(h[i-1], z[i-1])
		z[i] = new(big.Float).Sub(alpha[i], t3)
		z[i].Quo(z[i], l[i])
	}
	l[j] = big.NewFloat(1)
	z[j] = new(big.Float)
	f.P[j].c = new(big.Float)
	for i = j - 1; i >= 0; i-- {
		// f.P[i].c = z[i] - mu[i]*f.P[i+1].c
		t1 := new(big.Float).Mul(mu[i], f.P[i+1].c)
		f.P[i].c.Sub(z[i], t1)

		// f.P[i].b = (f.P[i+1].Y-f.P[i].Y)/h[i] - h[i]*(f.P[i+1].c+2*f.P[i].c)/3
		t2 := new(big.Float).Sub(f.P[i+1].Y, f.P[i].Y)
		t2.Quo(t2, h[i])
		t3 := new(big.Float).Mul(f.P[i].c, big.NewFloat(2))
		t3.Add(f.P[i+1].c, t3)
		t3.Mul(t3, h[i]).Quo(t3, big.NewFloat(3))
		f.P[i].b.Sub(t2, t3)

		// f.P[i].d = (f.P[i+1].c - f.P[i].c) / 3 / h[i]
		t4 := new(big.Float).Sub(f.P[i+1].c, f.P[i].c)
		t4.Quo(t4, big.NewFloat(3))
		f.P[i].d.Quo(t4, h[i])
	}
}

func (f *TabulatedFunction) SetOrder(new_value int) {
	f.iOrder = new_value
	f.changed = true
}

func (f *TabulatedFunction) SetTrapolation(new_value Trapolation) {
	f.trapolation = new_value
	f.changed = true
}

func (f *TabulatedFunction) AddPoint(Xn, Yn *big.Float, epoch uint32) *big.Float {
	var i int
	f.changed = true

	// Xn = math.Round(Xn*10000) / 10000
	scale := big.NewFloat(10000)
	val := new(big.Float).Mul(Xn, scale)
	half := big.NewFloat(0.5)
	if val.Sign() < 0 {
		half.Neg(half)
	}
	val.Add(val, half)
	intVal := new(big.Int)
	val.Int(intVal)
	XnRounded := new(big.Float).SetInt(intVal)
	XnRounded.Quo(XnRounded, scale)

	i, found := slices.BinarySearchFunc(f.P, TFPoint{X: Xn}, func(a, b TFPoint) int {
		return a.X.Cmp(b.X)
	})
	if found {
		//f.P[i].X = Xn
		f.P[i].epoch = epoch
		f.P[i].Y.Add(f.P[i].Y, Yn).Quo(f.P[i].Y, big.NewFloat(2))
		return new(big.Float).Set(f.P[i].Y)
	}
	f.P = slices.Insert(f.P, i, TFPoint{
		X:     XnRounded,
		Y:     new(big.Float).Set(Yn),
		epoch: epoch,
		b:     new(big.Float),
		c:     new(big.Float),
		d:     new(big.Float),
	})
	return Yn
}

func (f *TabulatedFunction) LoadConstant(new_Y, new_xmin, new_xmax *big.Float) {
	f.ixmin.Set(new_xmin)
	f.ixmax.Set(new_xmax)
	f.iymin.Set(new_Y)
	f.iymax.Set(f.iymin)
	f.P = append([]TFPoint{}, TFPoint{
		X: new(big.Float).Set(f.ixmin), Y: new(big.Float).Set(f.iymin), epoch: 0,
		b: new(big.Float), c: new(big.Float), d: new(big.Float),
	})
	f.istep.Sub(f.ixmax, f.ixmin)
	f.changed = false
}

func (f *TabulatedFunction) Normalise() {
	var i int
	ym := new(big.Float).Abs(f.iymax)
	yminAbs := new(big.Float).Abs(f.iymin)
	if yminAbs.Cmp(ym) > 0 {
		ym.Set(yminAbs)
	}

	for i = range f.P {
		f.P[i].Y.Quo(f.P[i].Y, ym)
	}
	f.changed = true
}

func (f *TabulatedFunction) Smooth() {
	for i := range f.P {
		if i == 0 || i == len(f.P)-1 {
			continue
		}
		// h2 := (f.P[i].X - f.P[i-1].X) + (f.P[i+1].X - f.P[i].X)
		h2 := new(big.Float).Sub(f.P[i].X, f.P[i-1].X)
		h2.Add(h2, new(big.Float).Sub(f.P[i+1].X, f.P[i].X))

		// d0 := (-3*f.P[i-1].Y + 4*f.P[i].Y - f.P[i+1].Y) / h2
		d0 := new(big.Float).Mul(f.P[i-1].Y, big.NewFloat(-3))
		d0.Add(d0, new(big.Float).Mul(f.P[i].Y, big.NewFloat(4)))
		d0.Sub(d0, f.P[i+1].Y)
		d0.Quo(d0, h2)

		// d1 := (f.P[i+1].Y - f.P[i-1].Y) / h2
		d1 := new(big.Float).Sub(f.P[i+1].Y, f.P[i-1].Y)
		d1.Quo(d1, h2)

		// d2 := (f.P[i-1].Y - 4*f.P[i].Y + 3*f.P[i+1].Y) / h2
		d2 := new(big.Float).Sub(f.P[i-1].Y, new(big.Float).Mul(f.P[i].Y, big.NewFloat(4)))
		d2.Add(d2, new(big.Float).Mul(f.P[i+1].Y, big.NewFloat(3)))
		d2.Quo(d2, h2)

		// d := (d0 + d1 + d2) / 3
		d := new(big.Float).Add(d0, d1).Add(d0, d2).Quo(d0, big.NewFloat(3))

		// f1_1 := (d*h2 + 3*f.P[i-1].Y + f.P[i+1].Y) / 4
		f1_1 := new(big.Float).Mul(d, h2)
		f1_1.Add(f1_1, new(big.Float).Mul(f.P[i-1].Y, big.NewFloat(3)))
		f1_1.Add(f1_1, f.P[i+1].Y).Quo(f1_1, big.NewFloat(4))

		// f1_2 := (d*h2 - f.P[i-1].Y - 3*f.P[i+1].Y) / -4
		f1_2 := new(big.Float).Mul(d, h2)
		f1_2.Sub(f1_2, f.P[i-1].Y)
		f1_2.Sub(f1_2, new(big.Float).Mul(f.P[i+1].Y, big.NewFloat(3)))
		f1_2.Quo(f1_2, big.NewFloat(-4))

		// f.P[i].Y = (f1_1 + f1_2) / 2
		f.P[i].Y.Add(f1_1, f1_2).Quo(f.P[i].Y, big.NewFloat(2))
	}
	f.changed = true
}

func (f *TabulatedFunction) Multiply(by *TabulatedFunction) {
	var i int
	var Yt []*big.Float

	if f.changed {
		f.update_spline()
	}
	Yt = make([]*big.Float, len(by.P))
	for i = range by.P {
		Yt[i] = new(big.Float).Mul(by.P[i].Y, f.F(by.P[i].X))
	}

	for i = range f.P {
		f.P[i].Y.Mul(f.P[i].Y, by.F(f.P[i].X))
	}

	for i = range by.P {
		f.AddPoint(by.P[i].X, Yt[i], by.P[i].epoch)
	}
	f.changed = true
}

func (f *TabulatedFunction) MultiplyByScalar(by *big.Float) {
	if f.changed {
		f.update_spline()
	}
	for i := range f.P {
		f.P[i].Y.Mul(f.P[i].Y, by)
		f.P[i].b.Mul(f.P[i].b, by)
		f.P[i].c.Mul(f.P[i].c, by)
		f.P[i].d.Mul(f.P[i].d, by)
	}
	f.iymin.Mul(f.iymin, by)
	f.iymax.Mul(f.iymax, by)
}

func (f *TabulatedFunction) Assign(s *TabulatedFunction) {
	f.ixmin.Set(s.ixmin)
	f.ixmax.Set(s.ixmax)
	f.iymin.Set(s.iymin)
	f.iymax.Set(s.iymax)
	f.istep.Set(s.istep)
	f.iOrder = s.iOrder

	f.P = make([]TFPoint, len(s.P))
	for i, p := range s.P {
		f.P[i] = TFPoint{
			X:     new(big.Float).Set(p.X),
			Y:     new(big.Float).Set(p.Y),
			b:     new(big.Float).Set(p.b),
			c:     new(big.Float).Set(p.c),
			d:     new(big.Float).Set(p.d),
			epoch: p.epoch,
		}
	}
	f.changed = true
}

func (f *TabulatedFunction) Merge(m *TabulatedFunction) {
	for i := range m.P {
		f.AddPoint(m.P[i].X, m.P[i].Y, m.P[i].epoch)
	}
	f.changed = true
}

func (f *TabulatedFunction) Integrate() *big.Float {
	var i, l int
	var tmp, dif *big.Float

	if f.changed {
		f.update_spline()
	}
	tmp = new(big.Float)
	l = len(f.P) - 1
	for i = 0; i < l; i++ {
		// dif = f.P[i+1].X - f.P[i].X
		dif = new(big.Float).Sub(f.P[i+1].X, f.P[i].X)
		// term = f.P[i].Y + dif*(f.P[i].b/2+dif*(f.P[i].c/3+dif*f.P[i].d/4))
		term := new(big.Float).Mul(f.P[i].d, dif)
		term.Quo(term, big.NewFloat(4))
		term.Add(term, new(big.Float).Quo(f.P[i].c, big.NewFloat(3)))
		term.Mul(term, dif)
		term.Add(term, new(big.Float).Quo(f.P[i].b, big.NewFloat(2)))
		term.Mul(term, dif)
		term.Add(term, f.P[i].Y)
		term.Mul(term, dif)
		tmp.Add(tmp, term)
	}
	return tmp
}

func (f *TabulatedFunction) Clear() {
	f.P = make([]TFPoint, 0)
	f.ixmin.SetInt64(0)
	f.ixmax.SetInt64(0)
	f.iymin.SetInt64(0)
	f.iymax.SetInt64(0)
	f.istep.SetInt64(0)
	f.changed = false
}

func (f *TabulatedFunction) MorePoints() {
	if f.changed {
		f.update_spline()
	}

	numPoints := len(f.P)
	if numPoints <= 1 {
		return
	}

	// The new slice will have the original points plus one interpolated point between each pair.
	newSize := numPoints + (numPoints - 1)
	newP := make([]TFPoint, 0, newSize)

	// Add the first point
	newP = append(newP, f.P[0])

	for i := 0; i < numPoints-1; i++ {
		p1 := f.P[i]
		p2 := f.P[i+1]

		// Add interpolated midpoint and then the next original point
		midX := new(big.Float).Add(p1.X, p2.X)
		midX.Quo(midX, big.NewFloat(2))
		newP = append(newP, TFPoint{
			X:     midX,
			Y:     f.F(midX),
			epoch: p2.epoch,
			b:     new(big.Float),
			c:     new(big.Float),
			d:     new(big.Float),
		})
		newP = append(newP, p2)
	}

	f.P = newP
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
		f.P[i].Y.Set(f.P[i].b)
		f.P[i].b.Mul(f.P[i].c, big.NewFloat(2))
		f.P[i].c.Mul(f.P[i].d, big.NewFloat(3))
		f.P[i].d.SetInt64(0)
	}
}

func (f *TabulatedFunction) Integral() {
	var i, j int
	var acc, prev_acc, r *big.Float

	if f.changed {
		f.update_spline()
	}
	j = len(f.P) - 1
	acc = new(big.Float)
	prev_acc = new(big.Float)
	if f.iOrder < 3 {
		f.P[0].d.Quo(f.P[0].c, big.NewFloat(3))
		f.P[0].c.Quo(f.P[0].b, big.NewFloat(2))
		f.P[0].b.Set(f.P[0].Y)
		f.P[0].Y.SetInt64(0)
		for i = 1; i <= j; i++ {
			r = new(big.Float).Sub(f.P[i].X, f.P[i-1].X)
			f.P[i].d.Quo(f.P[i].c, big.NewFloat(3))
			f.P[i].c.Quo(f.P[i].b, big.NewFloat(2))
			f.P[i].b.Set(f.P[i].Y)
			// f.P[i].Y = f.P[i-1].Y + r*(f.P[i-1].b+r*(f.P[i-1].c+r*f.P[i-1].d))
			term := new(big.Float).Mul(r, f.P[i-1].d)
			term.Add(f.P[i-1].c, term)
			term.Mul(r, term)
			term.Add(f.P[i-1].b, term)
			term.Mul(r, term)
			f.P[i].Y.Add(f.P[i-1].Y, term)
		}
		f.iOrder++
	} else {
		prev_acc.Quo(f.P[0].d, big.NewFloat(4))
		f.P[0].d.Quo(f.P[0].c, big.NewFloat(3))
		f.P[0].c.Quo(f.P[0].b, big.NewFloat(2))
		f.P[0].b.Set(f.P[0].Y)
		f.P[0].Y.SetInt64(0)
		for i = 1; i <= j; i++ {
			r = new(big.Float).Sub(f.P[i].X, f.P[i-1].X)
			acc.Quo(f.P[i].d, big.NewFloat(4))
			f.P[i].d.Quo(f.P[i].c, big.NewFloat(3))
			f.P[i].c.Quo(f.P[i].b, big.NewFloat(2))
			f.P[i].b.Set(f.P[i].Y)
			// f.P[i].Y = f.P[i-1].Y + r*(f.P[i-1].b+r*(f.P[i].c+r*(f.P[i-1].d+r*prev_acc)))
			term := new(big.Float).Mul(r, prev_acc)
			term.Add(f.P[i-1].d, term)
			term.Mul(r, term)
			term.Add(f.P[i].c, term)
			term.Mul(r, term)
			term.Add(f.P[i-1].b, term)
			term.Mul(r, term)
			f.P[i].Y.Add(f.P[i-1].Y, term)
			prev_acc.Set(acc)
		}
		f.changed = true
	}
}

func (f *TabulatedFunction) GetStep() *big.Float {
	if f.changed {
		f.update_spline()
	}
	return f.istep
}

func (f *TabulatedFunction) GetXmin() *big.Float {
	if f.changed {
		f.update_spline()
	}
	return f.ixmin
}

func (f *TabulatedFunction) GetXmax() *big.Float {
	if f.changed {
		f.update_spline()
	}
	return f.ixmax
}

func (f *TabulatedFunction) GetYmin() *big.Float {
	if f.changed {
		f.update_spline()
	}
	return f.iymin
}

func (f *TabulatedFunction) GetYmax() *big.Float {
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
		if a.epoch < epoch && b.epoch >= epoch {
			return 1
		}
		return a.X.Cmp(b.X)
	})
	i := slices.IndexFunc(f.P, func(p TFPoint) bool {
		return p.epoch < epoch
	})
	if i >= 0 {
		f.P = slices.Delete(f.P, i, len(f.P))
		f.changed = true
	}
}

// https://github.com/rsmith-nl/ps-lib/blob/main/grid.inc
// https://stackoverflow.com/a/20866012

func (f *TabulatedFunction) DrawPS(path string) error {
	ps, err := os.Create(path)
	if err != nil {
		return err
	}
	defer ps.Close()

	if f.changed {
		f.update_spline()
	}
	fmt.Fprintf(ps, `%%!PS
	%% This is the color that the grid is drawn in.
/grid_major_color {1 .6 .6} def
/grid_color {.7 1 1} def
/line_color {.5 .5 .5} def
/dot_color {.1 .1 .1} def
/radius 1 def
%% The line width used for the grid.
/grid_major_lw 1.5 def
/grid_lw .5 def

%% Every major-th line is drawn in a different color and thickness.
/major 10 def

%% Usage: dx dy w h gridwh
%% Draw a grid over a supplied width and height
/gridwh {
  4 dict begin
    /h exch def
    /w exch def
    /dy exch def
    /dx exch def
    gsave
        %% Set line width and color
        grid_lw setlinewidth
        grid_color setrgbcolor
        %% draw
        newpath
        %% vertical lines
        dx dx w {
            0 moveto
            0 h rlineto
        } for
        %% horizontal lines
        dy dy h {
            0 exch moveto
            w 0 rlineto
        } for
        stroke
        newpath
        grid_major_lw setlinewidth
        grid_major_color setrgbcolor
        %% every 10th line
        0 dx major mul w {
            0 moveto
            0 h rlineto
        } for
        0 dy major mul h {
            0 exch moveto
            w 0 rlineto
        } for
        stroke
    grestore
  end
} bind def

%% Distance between dimension point and start of witness line
/dimoffs 6 def
%% Distance that the witness line descends past the dimension line.
/dimext 20 def
%% Font for dimensions
/dimfont /Alegraya-Regular def
%% Font size
/dimscale 12 def
%% dimension text offset
/dimtextoffs 12 def
%% Dimension color
/dimcol {1 .1 .1 setrgbcolor} def
%% This defines the length of the arrow-head.
/dimhead 30 def

%% Usage: x1 y1 x2 y2 arrow_head x3 y3
%% Sees a line from x1,y1 to x2,y2 and draws an arrow head on the latter.
%% Returns x3 y3, leaving it to the user to draw the line (x1,y1)--(x3,y3).
/_arrow_head {
    9 dict begin
        /y2 exch def /x2 exch def /y1 exch def /x1 exch def
        /dx x2 x1 sub def /dy y2 y1 sub def /ang dy dx atan def
        /len dx dup mul dy dup mul add sqrt def
        /fact dimhead len 0.8 div div def
        gsave
            x2 y2 translate ang rotate
            newpath 0 0 moveto dimhead neg 4 {dup} repeat -.25 mul lineto
            .8 mul 0 lineto .25 mul lineto closepath fill
        grestore
        x2 dx fact mul sub y2 dy fact mul sub %% inside of the arrowhead
    end
} bind def

%% Usage: (text) _align_middle
/_align_middle {
    dup %% (text) (text)
    stringwidth pop %% (text) w
    -2 div 0 rmoveto
	dimfont findfont dimscale scalefont setfont
	dimcol show
} bind def

%% Draw a horizontal dimension
%% Usage x1 y1 x2 y2 offs (label) horizontal_dim
/horizontal_dim {
	gsave
	dimcol
    9 dict begin
        /label exch def
        /offs exch def
        /y2 exch def
        /x2 exch def
        /y1 exch def
        /x1 exch def
        /q y1 offs add def
        offs 0 ge {
            /v y1 dimext add def
            /w y1 offs add dimext add def
        } {
            /v y1 dimext sub def
            /w y1 offs add dimext sub def
        } ifelse
        %% Left witness line
        x1 v moveto x1 w lineto stroke
        %% Right witness line
        x2 v moveto x2 w lineto stroke
        %% arrow heads
        x2 q x1 q _arrow_head
        x1 q x2 q _arrow_head
        %% Dimension line
        moveto lineto stroke
        x1 x2 add 2 div q dimtextoffs add moveto label _align_middle
    end
	grestore
} bind def

%% Draw a vertical dimension
%% Usage x1 y1 x2 y2 offs (label) vertical_dim
/vertical_dim {
	gsave
	dimcol
    9 dict begin
        /label exch def
        /offs exch def
        /y2 exch def
        /x2 exch def
        /y1 exch def
        /x1 exch def
        /q x1 offs add def
        offs 0 ge {
            /v x1 dimext add def
            /w x1 offs add dimext add def
        } {
            /v x1 dimext sub def
            /w x1 offs add dimext sub def
        } ifelse
        %% Bottom witness line
        v y1 moveto w y1 lineto stroke
        %% Top witness line
        v y2 moveto w y2 lineto stroke
        %% arrow heads
        q y2 q y1 _arrow_head
        q y1 q y2 _arrow_head
        %% Dimension line
        moveto lineto stroke
        %% Rotated label
        q dimtextoffs sub y1 y2 add 2 div moveto
        gsave 90 rotate label _align_middle grestore
    end
	grestore
} bind def


`)

	fmt.Fprintf(ps, "/XValues [\n")
	for i, p := range f.P {
		fmt.Fprintf(ps, " %v\t%% %v\n", p.X, i)
	}
	fmt.Fprintf(ps, "] def\n")

	fmt.Fprintf(ps, "/YValues [\n")
	for i, p := range f.P {
		fmt.Fprintf(ps, " %v\t%% %v", p.Y, i)
		if i > 0 && i < len(f.P)-1 {
			// f.P[i-1].Y+(f.P[i+1].Y-f.P[i-1].Y)*(f.P[i].X-f.P[i-1].X)/(f.P[i+1].X-f.P[i-1].X)
			yPrev := f.P[i-1].Y
			yNext := f.P[i+1].Y
			xPrev := f.P[i-1].X
			xNext := f.P[i+1].X
			xCurr := f.P[i].X

			dy := new(big.Float).Sub(yNext, yPrev)
			dx1 := new(big.Float).Sub(xCurr, xPrev)
			dx2 := new(big.Float).Sub(xNext, xPrev)
			val := new(big.Float).Mul(dy, dx1)
			val.Quo(val, dx2)
			val.Add(yPrev, val)
			fmt.Fprintf(ps, " %v", val)
		}
		fmt.Fprintf(ps, "\n")
	}
	fmt.Fprintf(ps, "] def\n")

	fmt.Fprintf(ps, "/Xmin %v dup %v exch sub 0.01 mul abs sub def\n", f.ixmin, f.ixmax)
	fmt.Fprintf(ps, "/Xmax %v dup %v sub 0.01 mul abs add def\n", f.ixmax, f.ixmin)
	fmt.Fprintf(ps, "/Ymin %v dup %v exch sub 0.01 mul abs sub def\n", f.iymin, f.iymax)
	fmt.Fprintf(ps, "/Ymax %v dup %v sub 0.01 mul abs add def\n", f.iymax, f.iymin)

	fmt.Fprintf(ps, `
/Xsize Xmax Xmin sub def
/Ysize Ymax Ymin sub def
`)

	fmt.Fprintf(ps, `
/w currentpagedevice /PageSize get 0 get def
/h currentpagedevice /PageSize get 1 get def

w 10 div h 10 div w h gridwh

/Translate { %% x y Translate
	Ymin sub h mul Ysize div
	exch
	Xmin sub w mul Xsize div
	exch 
} bind def
`)

	fmt.Fprintf(ps, `
%% lines

newpath
line_color setrgbcolor
XValues 0 get YValues 0 get %% X[0] Y[0]
Translate
moveto                      %% move to first point
1 1 XValues length 1 sub {  %% i    push integer i = 1 .. length(XValues)-1 on each iteration
XValues                 %% i XVal    push X array
1 index                 %% i XVal i  copy i from stack
get                     %% i x       get ith X value from array
YValues                 %% i x YVal
2 index                 %% i x YVal i  i is 1 position deeper now, so 2 index instead of 1
get                     %% i x y
Translate
lineto                  %% i    line to next point
pop                     %%      discard index variable
} for
stroke
`)

	fmt.Fprintf(ps, `
%% dots

newpath
dot_color setrgbcolor
XValues 0 get YValues 0 get %% X[0] Y[0]
Translate
radius 0 360 arc           %% draw the first point
stroke
1 1 XValues length 1 sub {  %% i    push integer i = 1 .. length(XValues)-1 on each iteration
XValues                 %% i XVal    push X array
1 index                 %% i XVal i  copy i from stack
get                     %% i x       get ith X value from array
YValues                 %% i x YVal
2 index                 %% i x YVal i  i is 1 position deeper now, so 2 index instead of 1
get                     %% i x y
Translate
radius 0 360 arc        %% i    draw the next point
stroke
pop                     %%      discard index variable
} for
`)

	fmt.Fprintf(ps, `0 5 w 5 10 (%v - %v) horizontal_dim
	`, f.ixmin, f.ixmax)
	fmt.Fprintf(ps, `5 0 5 h 20 (%v - %v) vertical_dim
	`, f.iymin, f.iymax)

	fmt.Fprintf(ps, `

showpage
quit
`)

	return nil
}
