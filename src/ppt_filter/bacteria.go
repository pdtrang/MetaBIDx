package ppt_filter

import (
	"fmt"
)

type Bacteria struct {
	Signatures *StringSet
	Threshold float32
	Reported bool
}

func NewBacteria(t float32) *Bacteria {
	b := &Bacteria {
		Signatures: NewStringSet(),
		Threshold: t,
		Reported: false,
	}

	return b
}

func (b *Bacteria) PrintBacteria() {
	fmt.Println("Signatures: ", b.Signatures)
	fmt.Println("Threshold: ", b.Threshold)
	fmt.Println("Reported: ", b.Reported)
	
}

func (b *Bacteria) AddSignature(sig string) {
	temp := b.Signatures
	temp.Add(sig)
	b.Signatures = temp
}

func (b *Bacteria) ReachThreshold() bool {
	temp := b.Signatures
	if float32(temp.Size()) >= b.Threshold {
		return true
	} else {
		return false
	}
}
