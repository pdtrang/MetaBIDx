//-----------------------------------------------------------------------------
//
// Create a Bloom filter for each different genome by their signatures.
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
func (f *Filter) HashSignature(kmer []byte, is_first_kmer bool, gid uint16, ph int, gname string) {
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
	// fmt.Print(seq_name, ", phase = ", ph, ", kmer = ", string(kmer), ", hash values = ")
	if unique_to_genome {
		//fmt.Print(seq_name, ", phase = ", ph, ", kmer = ", string(kmer), ", hash values = ")
		for i := 0; i < len(idx); i++ {
			// fmt.Print(idx[i], "\t")
			f.table[idx[i]] = gid
		}
		// fmt.Print(", Unique")
		// fmt.Println()
	} else {
		for i := 0; i < len(idx); i++ {
			// fmt.Print(idx[i], "\t")
			f.table[idx[i]] = Dirty
		}
		// fmt.Print(", Dirty")
		// fmt.Println()
	}
}

func (f *Filter) HashSignatureWithWindow(kmer []byte, is_first_kmer bool, gid uint16, is_max_num_kmers bool) bool {
	unique_to_genome := true
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
			if is_max_num_kmers {
				f.table[idx[i]] = Dirty
			} else {
				f.table[idx[i]] = gid	
			}
		}
	} else {
		for i := 0; i < len(idx); i++ {
			f.table[idx[i]] = Dirty
		}
	}

	return unique_to_genome
}


