package tabulatedfunction

import (
	"encoding/json"
	"errors"
)

// Dump is a tabulated function dump
type Dump struct {
	iOrder int
	X, Y   []float64
}

// FromDump restores a tabulated function from a dump
func (f *TabulatedFunction) FromDump(d *Dump) error {
	if len(d.X) != len(d.Y) {
		return errors.New("incorrect Dump: len(X) != len(Y)")
	}
	l := len(d.X)
	f.X = make([]float64, l)
	f.Y = make([]float64, l)
	for i := range d.X {
		f.X[i] = d.X[i]
		f.Y[i] = d.Y[i]
	}
	f.iOrder = d.iOrder
	f.changed = true

	return nil
}

// Dump() generates dump for a tabulated function
func (f *TabulatedFunction) Dump() *Dump {
	l := len(f.X)
	X := make([]float64, l)
	Y := make([]float64, l)
	return &Dump{
		iOrder: f.iOrder,
		X:      X,
		Y:      Y,
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
