package tabulatedfunction

import (
	"encoding/json"
	"math"
	"slices"
)

// Dump is a serializable representation of a TabulatedFunction.
type Dump struct {
	Order       int         `json:"order"`
	Trapolation Trapolation `json:"trapolation"`
	Points      []TFPoint   `json:"points"`
	Precision   int         `json:"precision"`
}

// FromDump restores a tabulated function from a dump.
// It ensures the points are sorted by X before updating the spline.
func (f *TabulatedFunction) FromDump(d *Dump) {
	f.Order = d.Order
	f.Trapolation = d.Trapolation
	f.Precision = d.Precision

	f.P = make([]TFPoint, len(d.Points))
	copy(f.P, d.Points)
	// Ensure points follow the same rounding logic as AddPoint to prevent
	// numerical instability in spline calculation.
	f.P = make([]TFPoint, 0, len(d.Points))
	for _, p := range d.Points {
		p.X = math.Round(p.X*float64(f.Precision)) / float64(f.Precision)
		f.P = append(f.P, p)
	}

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

	// Deduplicate points with the same X coordinate to avoid division by zero in update_spline.
	if len(f.P) > 1 {
		k := 0
		for i := 1; i < len(f.P); i++ {
			if f.P[i].X == f.P[k].X {
				// Average Y values and take the highest epoch, matching AddPoint behavior.
				f.P[k].Y = (f.P[k].Y + f.P[i].Y) / 2.0
				if f.P[i].Epoch > f.P[k].Epoch {
					f.P[k].Epoch = f.P[i].Epoch
				}
			} else {
				k++
				f.P[k] = f.P[i]
			}
		}
		f.P = f.P[:k+1]
	}

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
		Order:       f.Order,
		Trapolation: f.Trapolation,
		Precision:   f.Precision,
		Points:      points,
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
