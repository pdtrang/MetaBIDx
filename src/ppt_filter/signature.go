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
func (f *Filter) HashSignature(kmer []byte, is_first_kmer bool, isPrimary bool, gid uint16, ph int, gname string, header string, kmer_pos int) {
	unique_to_genome := true
	idx := make([]int64, 0)
	for i := 0; i < len(f.HashFunction); i++ {
		// j := f.HashFunction[i].SlidingHashKmerModified(kmer, is_first_kmer, isPrimary)
		j := f.HashFunction[i].HashKmer(kmer)
		idx = append(idx, j)
		if f.table[j] != Empty && f.table[j] != gid {
			unique_to_genome = false
		}
	}

	if unique_to_genome {

		
		// fmt.Println("unique", gid, string(kmer), isPrimary, idx)

		for i := 0; i < len(idx); i++ {
			f.table[idx[i]] = gid
		}

		// store all positions of unique kmers in phase 2
		// if the position is already stored, skip it
		if ph == 2 {
			// fmt.Println(string(kmer), isPrimary, is_first_kmer, idx)
			// fmt.Println("Unique kmers", string(kmer), idx, gid)
			_, found := Find(f.Kmer_pos[header], kmer_pos)
		    
			if !found {
				f.Kmer_pos[header] = append(f.Kmer_pos[header], kmer_pos)
				// sort.Ints(f.Kmer_pos[header])	
			}			
		}
		
	} else {
		// fmt.Println("Dirty", gid, string(kmer), idx)
		for i := 0; i < len(idx); i++ {
			f.table[idx[i]] = Dirty
		}
	}
}

// Find takes a slice and looks for an element in it. If found it will
// return it's key, otherwise it will return -1 and a bool of false.
func Find(slice []int, val int) (int, bool) {
    for i, item := range slice {
        if item == val {
            return i, true
        }
    }
    return -1, false
}

