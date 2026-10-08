package tabulatedfunction

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"slices"
)

type Trapolation int

const (
	TrapolationLinear   Trapolation = 0
	TrapolationShift    Trapolation = 2
	TrapolationMinMax   Trapolation = 3
	TrapolationOpposite Trapolation = 4
	TrapolationNearest  Trapolation = 5
	TrapolationCosine   Trapolation = 6
)

type TabulatedFunction struct {
	ixmin, ixmax, iymin, iymax float64
	istep                      float64
	changed                    bool
	Order                      int
	Trapolation                Trapolation
	X                          []float64
	Y                          []float64
	epoch                      []uint32
	indices                    []uint32
	nextIndex                  uint32
}

// Create
func New() *TabulatedFunction {
	return &TabulatedFunction{
		Order:       1,
		Trapolation: TrapolationLinear,
		changed:     false,
		nextIndex:   1,
	}
}

// splinevalue
func (f *TabulatedFunction) F(xi float64) float64 {
	l := len(f.X)
	if l == 0 {
		return math.NaN()
	}
	k, found := slices.BinarySearch(f.X, xi)
	if found {
		return f.Y[k]
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

	return f._interpolate(xi, left, right, f.Trapolation)
}

func (f *TabulatedFunction) Trapolate(xi float64, trapolation Trapolation) float64 {
	l := len(f.X)
	if l == 0 {
		return math.NaN()
	}
	k, found := slices.BinarySearch(f.X, xi)

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
		if k == 0 {
			left = 0
		} else {
			left = k - 1
		}
	}

	return f._interpolate(xi, left, right, trapolation)
}

// _interpolate calculates the value at xi based on the surrounding points.
// It assumes left and right indices are correctly set.
func (f *TabulatedFunction) _interpolate(xi float64, left, right int, trapolation Trapolation) float64 {
	switch trapolation {
	case TrapolationLinear:
		// If left==right, it's extrapolation. Return the boundary value.
		if left == right || left < 0 || left >= len(f.X) {
			return f.Y[right]
		}
		if f.Order == 0 {
			return f.Y[left]
		}
		// Avoid division by zero if points are not distinct on X.
		dx := f.X[right] - f.X[left]
		if dx == 0 {
			return f.Y[left]
		}
		dy := f.Y[right] - f.Y[left]
		return f.Y[left] + dy*(xi-f.X[left])/dx

		// Similar to Linear then switch to opposite value
	case TrapolationOpposite:
		// If left==right, it's extrapolation. Return the boundary value.
		if left == right || left < 0 || left >= len(f.X) {
			return f.GetYmin() + f.GetYmax() - f.Y[right]
		}

		// Determine if xi is closer to the left or right point.
		midpoint := (f.X[left] + f.X[right]) / 2.0

		if xi <= midpoint {
			return f.Y[right]
		}
		return f.Y[left]

	case TrapolationShift:
		if left < 0 || left >= len(f.X) {
			return f.Y[right]
		}
		return (f.Y[left] + f.Y[right]) / 2.0

	case TrapolationMinMax:
		var v [2]float64

		if left < 0 {
			left = 0
		}
		if right >= len(f.X) {
			right = len(f.X) - 1
		}

		if f.changed {
			f.update_spline()
		}

		v[0] = 0
		v[1] = 0

		v[0] += math.Abs(f.iymin - f.Y[left])
		v[0] += math.Abs(f.iymin - f.Y[right])

		v[1] += math.Abs(f.iymax - f.Y[left])
		v[1] += math.Abs(f.iymax - f.Y[right])

		avg := (f.Y[left] + f.Y[right]) / 2.0

		if v[0] > v[1] {
			return (f.iymin + avg) / 2.0
		}
		return (f.iymax + avg) / 2.0

	case TrapolationNearest:
		if left == right || left < 0 || left >= len(f.X) {
			return f.Y[right]
		}
		if math.Abs(xi-f.X[left]) < math.Abs(xi-f.X[right]) {
			return f.Y[left]
		}
		return f.Y[right]

	case TrapolationCosine:
		if left == right || left < 0 || left >= len(f.X) {
			return f.Y[right]
		}
		dx := f.X[right] - f.X[left]
		if dx == 0 {
			return f.Y[left]
		}
		mu := (xi - f.X[left]) / dx
		mu2 := (1 - math.Cos(mu*math.Pi)) / 2.0
		return f.Y[left]*(1.0-mu2) + f.Y[right]*mu2
	}
	// This is unreachable if all Trapolation values are handled.
	// A panic is better than returning a magic number.
	panic("unhandled trapolation type")
}

func (f *TabulatedFunction) update_spline() {
	var i, j int

	f.changed = false
	i = len(f.X)
	if i == 0 {
		f.ixmin = 0
		f.ixmax = 0
		f.iymin = 0
		f.iymax = 0
		f.istep = 0
		return
	}
	j = i - 1
	f.ixmin = f.X[0]
	f.ixmax = f.X[j]
	f.iymin = f.Y[0]
	f.iymax = f.iymin
	for i = 1; i <= j; i++ {
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
			diff := f.X[i] - f.X[i-1]
			if diff < f.istep {
				f.istep = diff
			}
		}
	} else {
		f.istep = 0
	}
}

func (f *TabulatedFunction) SetOrder(new_value int) {
	// With splines removed, Order > 1 has no effect.
	// 0 = step, 1 = linear.
	f.Order = new_value
	f.changed = true
}

func (f *TabulatedFunction) SetTrapolation(new_value Trapolation) {
	f.Trapolation = new_value
	f.changed = true
}

func (f *TabulatedFunction) AddPoint(Xn, Yn float64, epoch uint32) float64 {
	f.changed = true

	n := len(f.X)
	// Fast-path: strictly increasing insertions
	if n == 0 || Xn > f.X[n-1] {
		f.X = append(f.X, Xn)
		f.Y = append(f.Y, Yn)
		f.epoch = append(f.epoch, epoch)
		f.indices = append(f.indices, f.nextIndex)
		f.nextIndex++
		return Yn
	}

	// Fast-path: updating the latest point
	if Xn == f.X[n-1] {
		f.epoch[n-1] = epoch
		old := f.Y[n-1]
		f.Y[n-1] = Yn
		f.indices[n-1] = f.nextIndex
		f.nextIndex++
		return old
	}

	// Fast-path: strictly decreasing insertions (prepend)
	if Xn < f.X[0] {
		f.X = slices.Insert(f.X, 0, Xn)
		f.Y = slices.Insert(f.Y, 0, Yn)
		f.epoch = slices.Insert(f.epoch, 0, epoch)
		f.indices = slices.Insert(f.indices, 0, f.nextIndex)
		f.nextIndex++
		return Yn
	}

	i, found := slices.BinarySearch(f.X, Xn)
	if found {
		f.epoch[i] = epoch
		old := f.Y[i]
		f.Y[i] = Yn
		f.indices[i] = f.nextIndex
		f.nextIndex++
		return old
	}

	f.X = slices.Insert(f.X, i, Xn)
	f.Y = slices.Insert(f.Y, i, Yn)
	f.epoch = slices.Insert(f.epoch, i, epoch)
	f.indices = slices.Insert(f.indices, i, f.nextIndex)
	f.nextIndex++

	return Yn
}

func (f *TabulatedFunction) LoadConstant(new_Y, new_xmin, new_xmax float64) {
	f.ixmin = new_xmin
	f.ixmax = new_xmax
	f.iymin = new_Y
	f.iymax = f.iymin
	f.X = []float64{f.ixmin}
	f.Y = []float64{f.iymin}
	f.epoch = []uint32{0}
	f.indices = []uint32{f.nextIndex}
	f.nextIndex++
	f.istep = f.ixmax - f.ixmin
	f.changed = false
}

func (f *TabulatedFunction) Normalise() {
	if f.changed {
		f.update_spline()
	}

	var i int
	ym := math.Max(math.Abs(f.iymin), math.Abs(f.iymax))

	if ym > 0 {
		for i = range f.Y {
			f.Y[i] /= ym
		}
		f.changed = true
	}
}

func (f *TabulatedFunction) NormaliseIndices() (uint32, uint32) {
	if len(f.indices) == 0 {
		f.nextIndex = 1
		return 0, 0
	}

	minIndex := f.indices[0]
	maxIndex := f.indices[0]
	for _, idx := range f.indices {
		if idx < minIndex {
			minIndex = idx
		}
		if idx > maxIndex {
			maxIndex = idx
		}
	}

	for i := range f.indices {
		f.indices[i] = f.indices[i] - minIndex + 1
	}
	f.nextIndex = maxIndex - minIndex + 2
	return 1, maxIndex - minIndex + 1
}

func (f *TabulatedFunction) Smooth() {
	n := len(f.Y)
	if n < 3 {
		return
	}

	prev := f.Y[0]
	curr := f.Y[1]
	for i := 1; i < n-1; i++ {
		next := f.Y[i+1]
		f.Y[i] = (prev + curr + next) / 3.0
		prev = curr
		curr = next
	}
	f.changed = true
}

func (f *TabulatedFunction) Multiply(by *TabulatedFunction) {
	if f.changed {
		f.update_spline()
	}
	if by.changed {
		by.update_spline()
	}
	if len(f.X) == 0 {
		return
	}
	if len(by.X) == 0 {
		f.Clear()
		return
	}

	// Two-pointer merge over already sorted X slices
	capGuess := len(f.X) + len(by.X)
	sortedX := make([]float64, 0, capGuess)
	i, j := 0, 0
	for i < len(f.X) && j < len(by.X) {
		if f.X[i] < by.X[j] {
			sortedX = append(sortedX, f.X[i])
			i++
		} else if f.X[i] > by.X[j] {
			sortedX = append(sortedX, by.X[j])
			j++
		} else {
			sortedX = append(sortedX, f.X[i])
			i++
			j++
		}
	}
	sortedX = append(sortedX, f.X[i:]...)
	sortedX = append(sortedX, by.X[j:]...)

	k := len(sortedX)
	newY := make([]float64, k)
	newEpoch := make([]uint32, k)
	newIndices := make([]uint32, k)

	for i, x := range sortedX {
		newY[i] = f.F(x) * by.F(x)
		newEpoch[i] = 0
		newIndices[i] = f.nextIndex
		f.nextIndex++
	}

	f.X = sortedX
	f.Y = newY
	f.epoch = newEpoch
	f.indices = newIndices

	f.changed = true
}

func (f *TabulatedFunction) MultiplyByScalar(by float64) {
	for i := range f.Y {
		f.Y[i] *= by
	}
	f.changed = true
}

func (f *TabulatedFunction) Assign(s *TabulatedFunction) {
	f.ixmin = s.ixmin
	f.ixmax = s.ixmax
	f.iymin = s.iymin
	f.iymax = s.iymax
	f.istep = s.istep
	f.Order = s.Order
	f.Trapolation = s.Trapolation
	f.nextIndex = s.nextIndex

	f.X = slices.Clone(s.X)
	f.Y = slices.Clone(s.Y)
	f.epoch = slices.Clone(s.epoch)
	f.indices = slices.Clone(s.indices)
	f.changed = true
}

func (f *TabulatedFunction) Merge(m *TabulatedFunction) {
	if len(m.X) == 0 {
		return
	}
	if len(f.X) == 0 {
		f.Assign(m)
		return
	}

	capGuess := len(f.X) + len(m.X)
	newX := make([]float64, 0, capGuess)
	newY := make([]float64, 0, capGuess)
	newEpoch := make([]uint32, 0, capGuess)
	newIndices := make([]uint32, 0, capGuess)

	i, j := 0, 0
	for i < len(f.X) && j < len(m.X) {
		if f.X[i] < m.X[j] {
			newX = append(newX, f.X[i])
			newY = append(newY, f.Y[i])
			newEpoch = append(newEpoch, f.epoch[i])
			newIndices = append(newIndices, f.indices[i])
			i++
		} else if f.X[i] > m.X[j] {
			newX = append(newX, m.X[j])
			newY = append(newY, m.Y[j])
			newEpoch = append(newEpoch, m.epoch[j])
			newIndices = append(newIndices, m.indices[j])
			j++
		} else {
			newX = append(newX, m.X[j])
			newY = append(newY, m.Y[j])
			newEpoch = append(newEpoch, m.epoch[j])
			newIndices = append(newIndices, m.indices[j])
			i++
			j++
		}
	}
	for ; i < len(f.X); i++ {
		newX = append(newX, f.X[i])
		newY = append(newY, f.Y[i])
		newEpoch = append(newEpoch, f.epoch[i])
		newIndices = append(newIndices, f.indices[i])
	}
	for ; j < len(m.X); j++ {
		newX = append(newX, m.X[j])
		newY = append(newY, m.Y[j])
		newEpoch = append(newEpoch, m.epoch[j])
		newIndices = append(newIndices, m.indices[j])
	}

	f.X = newX
	f.Y = newY
	f.epoch = newEpoch
	f.indices = newIndices
	f.nextIndex = max(f.nextIndex, m.nextIndex)
	f.changed = true
}

func (f *TabulatedFunction) Integrate() float64 {
	n := len(f.X)
	if n < 2 {
		return 0
	}
	var sum float64
	i := 0
	for ; i+2 < n; i += 2 {
		h1 := f.X[i+1] - f.X[i]
		h2 := f.X[i+2] - f.X[i+1]
		if h1 <= 0 || h2 <= 0 {
			sum += h1*(f.Y[i]+f.Y[i+1])/2.0 + h2*(f.Y[i+1]+f.Y[i+2])/2.0
			continue
		}
		term1 := (2.0 - h2/h1) * f.Y[i]
		term2 := ((h1 + h2) * (h1 + h2) / (h1 * h2)) * f.Y[i+1]
		term3 := (2.0 - h1/h2) * f.Y[i+2]
		sum += (h1 + h2) / 6.0 * (term1 + term2 + term3)
	}
	if i+1 < n {
		h := f.X[i+1] - f.X[i]
		sum += h * (f.Y[i] + f.Y[i+1]) / 2.0
	}
	return sum
}

func (f *TabulatedFunction) Clear() {
	f.X = f.X[:0]
	f.Y = f.Y[:0]
	f.epoch = f.epoch[:0]
	f.indices = f.indices[:0]
	f.ixmin = 0
	f.ixmax = 0
	f.iymin = 0
	f.iymax = 0
	f.istep = 0
	f.changed = false
	f.nextIndex = 1
}

func (f *TabulatedFunction) MorePoints() {
	if f.changed {
		f.update_spline()
	}

	numPoints := len(f.X)
	if numPoints <= 1 {
		return
	}

	newSize := numPoints + (numPoints - 1)
	newX := make([]float64, 0, newSize)
	newY := make([]float64, 0, newSize)
	newEpoch := make([]uint32, 0, newSize)
	newIndices := make([]uint32, 0, newSize)

	newX = append(newX, f.X[0])
	newY = append(newY, f.Y[0])
	newEpoch = append(newEpoch, f.epoch[0])
	newIndices = append(newIndices, f.indices[0])

	for i := 0; i < numPoints-1; i++ {
		x1, x2 := f.X[i], f.X[i+1]
		midX := (x1 + x2) / 2.0

		newX = append(newX, midX)
		newY = append(newY, f._interpolate(midX, i, i+1, f.Trapolation))
		newEpoch = append(newEpoch, f.epoch[i+1])
		newIndices = append(newIndices, f.nextIndex)
		f.nextIndex++

		newX = append(newX, x2)
		newY = append(newY, f.Y[i+1])
		newEpoch = append(newEpoch, f.epoch[i+1])
		newIndices = append(newIndices, f.indices[i+1])
	}

	f.X = newX
	f.Y = newY
	f.epoch = newEpoch
	f.indices = newIndices
	f.changed = true
}

func (f *TabulatedFunction) Derivative() {
	n := len(f.X)
	if n == 0 {
		return
	}
	if n == 1 {
		f.Y[0] = 0
		f.changed = true
		return
	}
	newY := make([]float64, n)
	if n == 2 {
		slope := (f.Y[1] - f.Y[0]) / (f.X[1] - f.X[0])
		newY[0] = slope
		newY[1] = slope
	} else {
		// Forward 3-point difference at the left boundary
		h0 := f.X[1] - f.X[0]
		h1 := f.X[2] - f.X[1]
		newY[0] = -f.Y[0]*(2*h0+h1)/(h0*(h0+h1)) + f.Y[1]*(h0+h1)/(h0*h1) - f.Y[2]*h0/(h1*(h0+h1))

		// Central 3-point difference for interior points
		for i := 1; i < n-1; i++ {
			hPrev := f.X[i] - f.X[i-1]
			hNext := f.X[i+1] - f.X[i]
			newY[i] = -f.Y[i-1]*hNext/(hPrev*(hPrev+hNext)) + f.Y[i]*(hNext-hPrev)/(hPrev*hNext) + f.Y[i+1]*hPrev/(hNext*(hPrev+hNext))
		}

		// Backward 3-point difference at the right boundary
		hPrev := f.X[n-2] - f.X[n-3]
		hLast := f.X[n-1] - f.X[n-2]
		newY[n-1] = f.Y[n-3]*hLast/(hPrev*(hPrev+hLast)) - f.Y[n-2]*(hPrev+hLast)/(hPrev*hLast) + f.Y[n-1]*(hPrev+2*hLast)/(hLast*(hPrev+hLast))
	}
	f.Y = newY
	f.changed = true
}

func (f *TabulatedFunction) Integral() {
	n := len(f.X)
	if n < 2 {
		return
	}
	prevY := f.Y[0]
	f.Y[0] = 0
	for i := 1; i < n; i++ {
		dx := f.X[i] - f.X[i-1]
		currY := f.Y[i]
		f.Y[i] = f.Y[i-1] + (prevY+currY)/2.0*dx
		prevY = currY
	}
	f.changed = true
}

func (f *TabulatedFunction) Expand(n int) {
	if f.changed {
		f.update_spline()
	}
	if len(f.X) < 2 {
		return
	}

	v1X, v1Y, v1Epoch := f.ixmin-f.istep, f.iymax, f.epoch[0]
	v2X, v2Y, v2Epoch := f.ixmax+f.istep, f.iymax, f.epoch[len(f.X)-1]

	midY := (f.iymin + f.iymax) / 2.0
	var indices []int
	var andices []int

	for step := 0; step < n; step++ {
		if f.changed {
			f.update_spline()
			midY = (f.iymin + f.iymax) / 2.0
		}
		if len(f.X) < 2 {
			break
		}

		getX := func(i int) float64 {
			if i == 0 {
				return v1X
			} else if i <= len(f.X) {
				return f.X[i-1]
			}
			return v2X
		}
		getY := func(i int) float64 {
			if i == 0 {
				return v1Y
			} else if i <= len(f.Y) {
				return f.Y[i-1]
			}
			return v2Y
		}
		getEpoch := func(i int) uint32 {
			if i == 0 {
				return v1Epoch
			} else if i <= len(f.epoch) {
				return f.epoch[i-1]
			}
			return v2Epoch
		}

		indices = indices[:0]
		andices = andices[:0]
		total := len(f.X) + 2
		for i := 0; i < total; i++ {
			if getY(i) > midY {
				indices = append(indices, i)
			} else {
				andices = append(andices, i)
			}
		}

		if len(indices) > 1 {
			maxDist := -1.0
			bestIdx := -1
			for j := 0; j < len(indices)-1; j++ {
				dist := getX(indices[j+1]) - getX(indices[j])
				if dist > maxDist {
					maxDist = dist
					bestIdx = j
				}
			}
			if bestIdx == -1 {
				break
			}
			idx1, idx2 := indices[bestIdx], indices[bestIdx+1]
			midX := (getX(idx1) + getX(idx2)) / 2.0
			if f.canInsertPoint(midX) {
				f.AddPoint(midX, (getY(idx1)+getY(idx2))/2.0, getEpoch(idx1))
			}
		}

		if len(andices) > 1 {
			maxDist := -1.0
			bestIdx := -1
			for j := 0; j < len(andices)-1; j++ {
				dist := getX(andices[j+1]) - getX(andices[j])
				if dist > maxDist {
					maxDist = dist
					bestIdx = j
				}
			}
			if bestIdx == -1 {
				break
			}
			idx1, idx2 := andices[bestIdx], andices[bestIdx+1]
			midX := (getX(idx1) + getX(idx2)) / 2.0
			if f.canInsertPoint(midX) {
				f.AddPoint(midX, (getY(idx1)+getY(idx2))/2.0, getEpoch(idx1))
			}
		}
	}

	f.changed = true
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
	return len(f.X)
}

func (f *TabulatedFunction) String() string {
	if f.changed {
		f.update_spline()
	}
	return fmt.Sprintf("\nTabulated function:\n"+
		"\tiOrder: %v; changed: %v\n"+
		"\tixmin: %v; ixmax: %v\n"+
		"\tiymin: %v; iymax: %v\n"+
		"\tistep: %v\n"+
		"\tPoints count: %v\n",
		f.Order, f.changed,
		f.ixmin, f.ixmax,
		f.iymin, f.iymax,
		f.istep, len(f.X))
}

func (f *TabulatedFunction) Epoch(epoch uint32) {
	w := 0
	for r := 0; r < len(f.X); r++ {
		if f.epoch[r] >= epoch {
			if w != r {
				f.X[w] = f.X[r]
				f.Y[w] = f.Y[r]
				f.epoch[w] = f.epoch[r]
				f.indices[w] = f.indices[r]
			}
			w++
		}
	}
	if w != len(f.X) {
		f.X = f.X[:w]
		f.Y = f.Y[:w]
		f.epoch = f.epoch[:w]
		f.indices = f.indices[:w]
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
	bw := bufio.NewWriter(ps)
	defer bw.Flush()
	defer ps.Close()

	if f.changed {
		f.update_spline()
	}

	minIndex, maxIndex := f.NormaliseIndices()

	// If there are no points, draw a blank page and exit to avoid errors.
	if len(f.X) == 0 {
		fmt.Fprintf(bw, `%%!PS
showpage
quit
`)
		return nil
	}

	fmt.Fprintf(bw, `%%!PS
	%% This is the color that the grid is drawn in.
/grid_major_color {1 .6 .6} def
/grid_color {.7 1 1} def
/line_color {.5 .5 .5} def
/dot_color {.1 .1 .1} def
/radius 1 def
/set_gray_by_index {
    MaxIdx MinIdx sub dup 0 ne {
        exch MinIdx sub exch div
        1.0 exch sub 0.85 mul
    } {
        pop pop 0.0
    } ifelse
    setgray
} bind def
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
/dimfont /Helvetica def
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
	dimfont findfont dimscale scalefont setfont
	dimcol
    dup %% (text) (text)
    stringwidth pop %% (text) w
    -2 div 0 rmoveto
	show
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
            /v y1 dimoffs add def
            /w q dimext add def
        } {
            /v y1 dimoffs sub def
            /w q dimext sub def
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
            /v x1 dimoffs add def
            /w q dimext add def
        } {
            /v x1 dimoffs sub def
            /w q dimext sub def
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

	fmt.Fprintf(bw, "/XValues [\n")
	xRange := f.ixmax - f.ixmin
	for i, x := range f.X {
		xNorm := 0.0
		if xRange != 0 {
			xNorm = (x - f.ixmin) / xRange
		}
		fmt.Fprintf(bw, " %v\t%% %v\n", xNorm, i)
	}
	fmt.Fprintf(bw, "] def\n")

	fmt.Fprintf(bw, "/YValues [\n")
	for i, y := range f.Y {
		fmt.Fprintf(bw, " %v\t%% %v", y, i)
		if i > 0 && i < len(f.X)-1 {
			yPrev := f.Y[i-1]
			yNext := f.Y[i+1]
			xPrev := f.X[i-1]
			xNext := f.X[i+1]
			xCurr := f.X[i]

			dy := yNext - yPrev
			dx1 := xCurr - xPrev
			dx2 := xNext - xPrev
			if dx2 != 0 {
				val := yPrev + dy*dx1/dx2
				fmt.Fprintf(bw, "\t%% interp: %v", val)
			}
		}
		fmt.Fprintf(bw, "\n")
	}
	fmt.Fprintf(bw, "] def\n")

	fmt.Fprintf(bw, "/ColorValues [\n")
	for i, idx := range f.indices {
		fmt.Fprintf(bw, " %v\t%% %v\n", idx, i)
	}
	fmt.Fprintf(bw, "] def\n")

	fmt.Fprintf(bw, "/MinIdx %v def\n", minIndex)
	fmt.Fprintf(bw, "/MaxIdx %v def\n", maxIndex)

	fmt.Fprintf(bw, "/Xmin 0 def\n")
	fmt.Fprintf(bw, "/Xmax 1 def\n")
	fmt.Fprintf(bw, "/Ymin %v def\n", f.iymin)
	fmt.Fprintf(bw, "/Ymax %v def\n", f.iymax)

	fmt.Fprintf(bw, `
/Xsize Xmax Xmin sub def
/Ysize Ymax Ymin sub dup 0 eq { pop 1.0 } if def
`)

	fmt.Fprintf(bw, `
/w currentpagedevice /PageSize get 0 get def
/h currentpagedevice /PageSize get 1 get def

w 10 div h 10 div w h gridwh

/Translate { %% x y Translate
	Ymin sub Ysize div h mul
	exch
	Xmin sub Xsize div w mul
	exch 
} bind def
`)

	fmt.Fprintf(bw, `
%% lines

1 1 XValues length 1 sub {  %% i    push integer i = 1 .. length(XValues)-1 on each iteration
newpath
ColorValues 1 index get set_gray_by_index
XValues                 %% i XVal
1 index 1 sub           %% i XVal i-1
get                     %% i x_{i-1}
YValues                 %% i x_{i-1} YVal
2 index 1 sub           %% i x_{i-1} YVal i-1
get                     %% i x_{i-1} y_{i-1}
Translate
moveto
XValues                 %% i XVal
1 index                 %% i XVal i
get                     %% i x_i
YValues                 %% i x_i YVal
2 index                 %% i x_i YVal i
get                     %% i x_i y_i
Translate
lineto
stroke
pop                     %% discard index variable
} for
`)

	fmt.Fprintf(bw, `
%% dots

newpath
ColorValues 0 get set_gray_by_index
XValues 0 get YValues 0 get %% X[0] Y[0]
Translate
radius 0 360 arc           %% draw the first point
stroke
1 1 XValues length 1 sub {  %% i    push integer i = 1 .. length(XValues)-1 on each iteration
ColorValues 1 index get set_gray_by_index
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

	fmt.Fprintf(bw, `0 5 w 5 10 (%v - %v) horizontal_dim
	`, f.ixmin, f.ixmax)
	fmt.Fprintf(bw, `5 0 5 h 20 (%v - %v) vertical_dim
	`, f.iymin, f.iymax)

	fmt.Fprintf(bw, `

showpage
quit
`)

	return nil
}

func (f *TabulatedFunction) canInsertPoint(x float64) bool {
	if f.changed {
		f.update_spline()
	}
	if len(f.X) == 0 {
		return true
	}
	k, found := slices.BinarySearch(f.X, x)
	if found {
		return false
	}
	if k > 0 && x-f.X[k-1] < f.istep {
		return false
	}
	if k < len(f.X) && f.X[k]-x < f.istep {
		return false
	}
	return true
}
