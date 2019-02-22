//-----------------------------------------------------------------------------
// Author: Vinhthuy Phan, 2018
//
// Three phases will effectively create a Bloom filter for each different genome.
//
//-----------------------------------------------------------------------------
package ppt_filter

import (
	// "fmt"
)

const Dirty = uint16(65534)

//-----------------------------------------------------------------------------
// If all slots have either 0 (clean) or gid, then kmer is unique.
// If so, set these slots to gid.  If not, set them to Dirty.
//-----------------------------------------------------------------------------
func (f *Filter) Phase1Hash(kmer []byte, is_first_kmer bool, gid uint16) {
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
		for i := 0; i < len(idx); i++ {
			f.table[idx[i]] = gid
		}
	} else {
		for i := 0; i < len(idx); i++ {
			f.table[idx[i]] = Dirty
		}
	}
}

//-----------------------------------------------------------------------------
// If all slots have either 0 (clean) or gid, then kmer is unique.
// If so, save the kmer and set all slots to 0. If not, set them to Dirty.
//-----------------------------------------------------------------------------
func (f *Filter) Phase2Hash(kmer []byte, is_first_kmer bool, gid uint16) [][]byte {
	unique_to_genome := true
	unique_kmers := make([][]byte, 0)
	idx := make([]int64, 0)
	for i := 0; i < len(f.HashFunction); i++ {
		j := f.HashFunction[i].SlidingHashKmer(kmer, is_first_kmer)
		idx = append(idx, j)
		if f.table[j] != 0 && f.table[j] != gid {
			unique_to_genome = false
		}
	}
	if unique_to_genome {
		for i := 0; i < len(idx); i++ {
			f.table[idx[i]] = 0
			unique_kmers = append(unique_kmers, kmer)
		}
	} else {
		for i := 0; i < len(idx); i++ {
			f.table[idx[i]] = Dirty
		}
	}
	return unique_kmers
}

//-----------------------------------------------------------------------------
// Given unique kmers, insert them into the filter.
//-----------------------------------------------------------------------------
func (f *Filter) Phase3Hash(unique_kmers [][]byte, gid uint16) {
	for i := 0; i < len(unique_kmers); i++ {
		for j := 0; j < len(f.HashFunction); j++ {
			k := f.HashFunction[j].HashKmer(unique_kmers[i])
			// fmt.Println()
			// fmt.Println(string(unique_kmers[i]))
			f.table[k] = gid

		}
	}

}

//-----------------------------------------------------------------------------
