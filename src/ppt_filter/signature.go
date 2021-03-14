//-----------------------------------------------------------------------------
//
// Create a Bloom filter for each different genome by their signatures.
//
//-----------------------------------------------------------------------------
package ppt_filter

import (
	// "fmt" 
	"sync"
)

const Dirty = uint16(65534)

//-----------------------------------------------------------------------------
// If all slots have either 0 (clean) or gid, then kmer is unique.
// If so, set these slots to gid.  If not, set them to Dirty.
//-----------------------------------------------------------------------------
func (f *Filter) HashSignature(kmer []byte, gid uint16, ph int, header string, kmer_pos int, mutex *sync.Mutex) {

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
		for i := 0; i < len(idx); i++ {
			// mut := f.GetLock(idx[i])
			mut := int(idx[i])
			f.lock[mut].Lock()

			f.table[idx[i]] = gid	
			
			f.lock[mut].Unlock()
		}
	} else {

		for i := 0; i < len(idx); i++ {
			// mut := f.GetLock(idx[i])
			mut := int(idx[i])
			f.lock[mut].Lock()

			f.table[idx[i]] = Dirty
			
			f.lock[mut].Unlock()
		}

	}
	

	if unique_to_genome {
		if ph == 2 {
            // store all positions of unique kmers in phase 2
			f.GetPositionofUniqueKmer(kmer_pos, header, mutex)
		}
	}
	
}

func (f *Filter) GetPositionofUniqueKmer(kmer_pos int, header string, mutex *sync.Mutex){
	
	mutex.Lock()
	_, found := Find(f.Kmer_pos[header], kmer_pos)
	
	// if the position is already stored, skip it
	if !found {
		
		f.Kmer_pos[header] = append(f.Kmer_pos[header], kmer_pos)
		
	}	
	mutex.Unlock()

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

