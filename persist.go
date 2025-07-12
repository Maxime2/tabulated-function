package tabulatedfunction

import (
	"cmp"
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
	f.P = slices.Clone(d.Points)

	// Ensure points are sorted, as they may come from an untrusted source.
	slices.SortFunc(f.P, func(a, b TFPoint) int {
		return cmp.Compare(a.X, b.X)
	})

	f.update_spline()
}

// Dump generates a serializable dump for a tabulated function.
func (f *TabulatedFunction) Dump() *Dump {
	return &Dump{
		Order:  f.iOrder,
		Points: slices.Clone(f.P),
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
