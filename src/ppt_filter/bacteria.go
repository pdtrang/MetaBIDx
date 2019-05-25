package ppt_filter

import (
	"fmt"
	"time"
)

type Bacteria struct {
	Signatures *Int64Set
	UpperThreshold float32
	LowerThreshold float32
	Reported bool
	QueryTime time.Duration
}

func NewBacteria(ut float32, lt float32) *Bacteria {
	return &Bacteria {
		Signatures: NewInt64Set(),
		UpperThreshold: ut,
		LowerThreshold: lt,
		Reported: false,
		QueryTime: time.Duration(0),
	}
}

func (b *Bacteria) AddSignature(j int64) {
	if !b.Signatures.Has(j) {
		b.Signatures.Add(j)
	}
}

func (b *Bacteria) PrintBacteria() {
	fmt.Println("Signatures: ", b.Signatures)
	fmt.Println("UpperThreshold: ", b.UpperThreshold)
	fmt.Println("LowerThreshold: ", b.LowerThreshold)
	fmt.Println("Reported: ", b.Reported)
}

func (b *Bacteria) ReachUpperThreshold() bool {
	if float32(b.Signatures.Size()) >= b.UpperThreshold {
		return true
	} else {
		return false
	}
}

func (b *Bacteria) ReachLowerThreshold() bool {
	if float32(b.Signatures.Size()) >= b.LowerThreshold {
		return true
	} else {
		return false
	}	
}