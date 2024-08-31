package tabulatedfunction

import (
	"encoding/json"
	"slices"
)

// Dump is a tabulated function dump
type Dump struct {
	iOrder int
	P      []TFPoint
}

// FromDump restores a tabulated function from a dump
func (f *TabulatedFunction) FromDump(d *Dump) error {
	f.P = slices.Clone(d.P)
	f.iOrder = d.iOrder
	f.update_spline()

	return nil
}

// Dump() generates dump for a tabulated function
func (f *TabulatedFunction) Dump() *Dump {
	return &Dump{
		iOrder: f.iOrder,
		P:      slices.Clone(f.P),
	}
}

// Marshals to JSON from network
func (f *TabulatedFunction) Marshal() ([]byte, error) {
	return json.Marshal(f.Dump())
}

// Restores tabulated function from a JSON blob
func Unmarshal(bytes []byte) (*TabulatedFunction, error) {
	var dump Dump
	if err := json.Unmarshal(bytes, &dump); err != nil {
		return nil, err
	}
	f := New()
	if err := f.FromDump(&dump); err != nil {
		return nil, err
	}
	return f, nil
}
