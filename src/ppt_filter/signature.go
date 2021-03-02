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
// func (f *Filter) HashSignature(kmer []byte, is_first_kmer bool, isPrimary bool, gid uint16, ph int, gname string, header string, kmer_pos int) {
func (f *Filter) HashSignature(kmer []byte, gid uint16, ph int, header string, kmer_pos int, mutex *sync.Mutex) {
// func (f *Filter) HashSignature(kmer []byte, gid uint16, ph int, header string, kmer_pos int) {

	unique_to_genome := true
	idx := make([]int64, 0)
	// fmt.Println(string(kmer), gid)

	// mut:= f.GetLock()
	
	for i := 0; i < len(f.HashFunction); i++ {
		// j := f.HashFunction[i].SlidingHashKmerModified(kmer, is_first_kmer, isPrimary)
		j := f.HashFunction[i].HashKmer(kmer)
		idx = append(idx, j)
		if f.table[j] != Empty && f.table[j] != gid {
			unique_to_genome = false
		}
	}

	var wg sync.WaitGroup

	if unique_to_genome {

		for i := 0; i < len(idx); i++ {
			wg.Add(1)
			go func(entry int64, i int, kmer_pos int, header string) {
				defer wg.Done()

				mut := f.GetLock(entry)

				f.lock[mut].Lock()
				f.table[entry] = gid	
				f.lock[mut].Unlock()

				if i == len(f.HashFunction) { 
					// store all positions of unique kmers in phase 2
					if ph == 2 {
						f.GetPositionofUniqueKmer(kmer_pos, header)			
					}
				}
				
			}(idx[i], i, kmer_pos, header)
			
		}

		wg.Wait()

		
		
		
	} else {
		// fmt.Println("Dirty", gid, string(kmer), idx)
		for i := 0; i < len(idx); i++ {
			wg.Add(1)
			go func(entry int64) {
				defer wg.Done()

				mut := f.GetLock(entry)

				f.lock[mut].Lock()
				f.table[entry] = Dirty	
				f.lock[mut].Unlock()

			}(idx[i])
		}

		wg.Wait()
	}
	
}

func (f *Filter) GetPositionofUniqueKmer(kmer_pos int, header string){
	_, found := Find(f.Kmer_pos[header], kmer_pos)
    
    // if the position is already stored, skip it
	if !found {
		f.Kmer_pos[header] = append(f.Kmer_pos[header], kmer_pos)
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

