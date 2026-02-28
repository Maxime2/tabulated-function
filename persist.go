package tabulatedfunction

import (
	"encoding/json"
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
	copy(f.P, d.Points)

	// Ensure points are sorted, as they may come from an untrusted source.
	slices.SortFunc(f.P, func(a, b TFPoint) int {
		if a.X < b.X {
			return -1
		}
		if a.X > b.X {
			return 1
		}
		return 0
	})

	f.ixmin = 0
	f.ixmax = 0
	f.iymin = 0
	f.iymax = 0
	f.istep = 0

	f.update_spline()
}

// Dump generates a serializable dump for a tabulated function.
func (f *TabulatedFunction) Dump() *Dump {
	points := make([]TFPoint, len(f.P))
	copy(points, f.P)
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
