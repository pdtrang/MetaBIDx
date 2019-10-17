//-----------------------------------------------------------------------------
//
// Create a Bloom filter for each different genome by their signatures.
//
//-----------------------------------------------------------------------------
package ppt_filter

import (
	"fmt"
)

const Dirty = uint16(65534)

//-----------------------------------------------------------------------------
// If all slots have either 0 (clean) or gid, then kmer is unique.
// If so, set these slots to gid.  If not, set them to Dirty.
//-----------------------------------------------------------------------------
func (f *Filter) HashSignature(kmer []byte, is_first_kmer bool, gid uint16, ph int) {
	unique_to_genome := true
	idx := make([]int64, 0)
	for i := 0; i < len(f.HashFunction); i++ {
		j := f.HashFunction[i].SlidingHashKmer(kmer, is_first_kmer)
		idx = append(idx, j)
		if f.table[j] != 0 && f.table[j] != gid {
			unique_to_genome = false
		}
	}

	// fmt.Println(len(idx))
	if unique_to_genome {
		fmt.Print("Phase = ", ph, ", kmer = ", string(kmer), "hash values = ")
		for i := 0; i < len(idx); i++ {
			fmt.Print(idx[i], "\t")
			f.table[idx[i]] = gid
		}
		fmt.Println()
	} else {
		for i := 0; i < len(idx); i++ {
			f.table[idx[i]] = Dirty
		}
	}
}

