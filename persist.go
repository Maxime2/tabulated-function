package tabulatedfunction

import (
	"encoding/json"
	"math/big"
	"slices"
)

// Dump is a serializable representation of a TabulatedFunction.
type Dump struct {
	Order  int       `json:"order"`
	Points []TFPoint `json:"points"`
}

// FromDump restores a tabulated function from a dump.
// It ensures the points are sorted by X before updating the spline.
func (f *TabulatedFunction) FromDump(d *Dump) {
	f.iOrder = d.Order

	f.P = make([]TFPoint, len(d.Points))
	for i, p := range d.Points {
		f.P[i] = TFPoint{
			epoch: p.epoch,
		}
		if p.X != nil {
			f.P[i].X = new(big.Float).Set(p.X)
		}
		if p.Y != nil {
			f.P[i].Y = new(big.Float).Set(p.Y)
		}
	}

	// Ensure points are sorted, as they may come from an untrusted source.
	slices.SortFunc(f.P, func(a, b TFPoint) int {
		if a.X == nil && b.X == nil {
			return 0
		}
		if a.X == nil {
			return -1
		}
		if b.X == nil {
			return 1
		}
		return a.X.Cmp(b.X)
	})

	if f.ixmin == nil {
		f.ixmin = new(big.Float)
	}
	if f.ixmax == nil {
		f.ixmax = new(big.Float)
	}
	if f.iymin == nil {
		f.iymin = new(big.Float)
	}
	if f.iymax == nil {
		f.iymax = new(big.Float)
	}
	if f.istep == nil {
		f.istep = new(big.Float)
	}

	f.update_spline()
}

// Dump generates a serializable dump for a tabulated function.
func (f *TabulatedFunction) Dump() *Dump {
	points := make([]TFPoint, len(f.P))
	for i, p := range f.P {
		points[i] = TFPoint{
			epoch: p.epoch,
		}
		if p.X != nil {
			points[i].X = new(big.Float).Set(p.X)
		}
		if p.Y != nil {
			points[i].Y = new(big.Float).Set(p.Y)
		}
		if p.b != nil {
			points[i].b = new(big.Float).Set(p.b)
		}
		if p.c != nil {
			points[i].c = new(big.Float).Set(p.c)
		}
		if p.d != nil {
			points[i].d = new(big.Float).Set(p.d)
		}
	}
	return &Dump{
		Order:  f.iOrder,
		Points: points,
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

	// Set defaults for a newly unmarshaled function.
	// `update_spline` is called within FromDump.
	f.trapolation = TrapolationSpline

	return nil
}
