package ppt_filter

import (
	"fmt"
	"time"
	"sync"
)

type Bacteria struct {
	ID uint16
	Signatures *Int64Set
	UpperThreshold float64
	LowerThreshold float64
	Reported bool
	QueryTime time.Duration
	Mutex sync.Mutex
}

func NewBacteria(id uint16, ut float64, lt float64) *Bacteria {
	return &Bacteria {
		ID: id, 
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

func (b *Bacteria) AppendSignature(j int64) {
	b.Signatures.Add(j)
}

func (b *Bacteria) PrintBacteria() {
	fmt.Println("ID", b.ID)
	fmt.Println("Signatures: ", b.Signatures)
	fmt.Println("UpperThreshold: ", b.UpperThreshold)
	fmt.Println("LowerThreshold: ", b.LowerThreshold)
	fmt.Println("Reported: ", b.Reported)
}

func (b *Bacteria) ReachUpperThreshold() bool {
	if float64(b.Signatures.Size()) >= b.UpperThreshold {
		return true
	} else {
		return false
	}
}

func (b *Bacteria) ReachLowerThreshold() bool {
	if float64(b.Signatures.Size()) >= b.LowerThreshold {
		return true
	} else {
		return false
	}	
}