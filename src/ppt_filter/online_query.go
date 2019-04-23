package ppt_filter

import (
	// "fmt"
	"log"
)

//-----------------------------------------------------------------------------
// For each genome, store all the kmers which are matched.
//-----------------------------------------------------------------------------
func (f *Filter) OnlineQuery(read []byte, bacteria_map map[int64]*Bacteria) {

	kmer_scanner := NewKmerScanner(read, f.K)

	for kmer_scanner.Scan() {

		idx := int64(0)
        j := int64(0)
		for i := 0; i < len(f.HashFunction); i++ {
			j = f.HashFunction[i].SlidingHashKmer(kmer_scanner.Kmer, kmer_scanner.IsFirstKmer)
			
			if int64(f.table[j]) != int64(65534) && int64(f.table[j]) != int64(0) {
				idx = int64(f.table[j])
				
				if !(bacteria_map[idx].Signatures.Has(j)) {
					bacteria_map[idx].Signatures.Add(j)
				}

				// fmt.Println("Updated signatures:", idx, bacteria_map[idx].Signatures)
				if bacteria_map[idx].ReachThreshold() && (bacteria_map[idx].Reported == false) {
					// fmt.Println("Found bacteria ", idx)
					log.Printf("Found bacteria %d", idx)
					bacteria_map[idx] = &Bacteria{bacteria_map[idx].Signatures, bacteria_map[idx].Threshold, true}
				}
			}

		}
		
    }	

}
