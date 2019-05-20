package ppt_filter

import (
	// "fmt"
)


func FindMajority(gidx map[uint16][][]byte) uint16{
	idx := uint16(0)
	value := 0

	for k, v := range gidx {
		if len(v) > value {
			value = len(v)
			idx = k
		}
	}

	return idx
}

func OneOrNothing(gidx map[uint16][][]byte) (uint16, bool){

	if len(gidx) == 1 {
		for k, _ := range gidx {
			return k, true
		}

	}
	
	return uint16(0), false
}