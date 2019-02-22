package ppt_filter

import (
	"fmt"
)

//-----------------------------------------------------------------------------
func (f *Filter) Query2(read []byte) uint16{

	kmer_scanner := NewKmerScanner(read, f.K)
	count := make(map[uint16]int)
	max_count, max_loc := 0, uint16(0)
	second_max_loc := uint16(0)

	for kmer_scanner.Scan() {
  
        j := int64(0)
		for i := 0; i < len(f.HashFunction); i++ {
			j = f.HashFunction[i].SlidingHashKmer(kmer_scanner.Kmer, kmer_scanner.IsFirstKmer)
			
			count[f.table[j]]++
			
			if f.table[j] != uint16(0) && f.table[j] != uint16(65534) {
				if count[f.table[j]] > max_count{
					second_max_loc = max_loc

					max_count = count[f.table[j]]
					max_loc = f.table[j]

				}	
			}
		}
    }	

    fmt.Print(count,",")    
	if max_loc != 0 {
		return max_loc
	} else {
		return second_max_loc
	}

}