package tabulatedfunction

import (
	"cmp"
	"encoding/json"
	"slices"
)

// Dump is a serializable representation of a TabulatedFunction.
type Dump struct {
	Order       int         `json:"order"`
	Trapolation Trapolation `json:"trapolation"`
	X           []float64   `json:"x"`
	Y           []float64   `json:"y"`
	Epoch       []uint32    `json:"epoch,omitempty"`
}

// FromDump restores a tabulated function from a dump.
// It ensures the points are sorted by X before updating internal bounds.
func (f *TabulatedFunction) FromDump(d *Dump) {
	f.Order = d.Order
	f.Trapolation = d.Trapolation

	n := min(len(d.X), len(d.Y))
	if n == 0 {
		f.Clear()
		f.Order = d.Order
		f.Trapolation = d.Trapolation
		return
	}

	f.X = make([]float64, n)
	f.Y = make([]float64, n)
	f.epoch = make([]uint32, n)
	copy(f.X, d.X[:n])
	copy(f.Y, d.Y[:n])
	copy(f.epoch, d.Epoch)

	// Ensure points are sorted by X
	if !slices.IsSorted(f.X) {
		type point struct {
			x, y  float64
			epoch uint32
		}
		pts := make([]point, n)
		for i := 0; i < n; i++ {
			pts[i] = point{x: f.X[i], y: f.Y[i], epoch: f.epoch[i]}
		}
		slices.SortFunc(pts, func(a, b point) int {
			return cmp.Compare(a.x, b.x)
		})
		for i := 0; i < n; i++ {
			f.X[i] = pts[i].x
			f.Y[i] = pts[i].y
			f.epoch[i] = pts[i].epoch
		}
	}

	// Deduplicate points with the same X coordinate
	if len(f.X) > 1 {
		k := 0
		count := 1
		sumY := f.Y[0]
		for i := 1; i < len(f.X); i++ {
			if f.X[i] == f.X[k] {
				sumY += f.Y[i]
				count++
				if f.epoch[i] > f.epoch[k] {
					f.epoch[k] = f.epoch[i]
				}
			} else {
				f.Y[k] = sumY / float64(count)
				k++
				f.X[k] = f.X[i]
				f.Y[k] = f.Y[i]
				f.epoch[k] = f.epoch[i]
				sumY = f.Y[k]
				count = 1
			}
		}
		f.Y[k] = sumY / float64(count)
		f.X = f.X[:k+1]
		f.Y = f.Y[:k+1]
		f.epoch = f.epoch[:k+1]
	}

	f.indices = make([]uint32, len(f.X))
	f.nextIndex = 1
	for i := range f.indices {
		f.indices[i] = f.nextIndex
		f.nextIndex++
	}

	f.update_spline()
}

// Dump generates a serializable dump for a tabulated function.
func (f *TabulatedFunction) Dump() *Dump {
	return &Dump{
		Order:       f.Order,
		Trapolation: f.Trapolation,
		X:           slices.Clone(f.X),
		Y:           slices.Clone(f.Y),
		Epoch:       slices.Clone(f.epoch),
	}
}

// MarshalJSON implements the json.Marshaler interface for TabulatedFunction.
func (f *TabulatedFunction) MarshalJSON() ([]byte, error) {
	return json.Marshal(f.Dump())
}

// UnmarshalJSON implements the json.Unmarshaler interface for TabulatedFunction.
func (f *TabulatedFunction) UnmarshalJSON(bytes []byte) error {
	var dump Dump
	if err := json.Unmarshal(bytes, &dump); err != nil {
		return err
	}

	// The json.Unmarshal call on the parent struct has already allocated
	// a zero-value TabulatedFunction for us. We just need to populate it.
	f.FromDump(&dump)

	return nil
}
