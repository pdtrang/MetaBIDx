package ppt_filter

import (
	"fmt"
)

//-----------------------------------------------------------------------------
// For each genome, store all the kmers which are matched.
//-----------------------------------------------------------------------------
func (f *Filter) OnlineQuery(read []byte, bacteria_map map[int64]*Bacteria) {

	kmer_scanner := NewKmerScanner(read, f.K)

	for kmer_scanner.Scan() {

		gidx := make([]int64, 0)
        j := int64(0)
		for i := 0; i < len(f.HashFunction); i++ {
			j = f.HashFunction[i].SlidingHashKmer(kmer_scanner.Kmer, kmer_scanner.IsFirstKmer)
			
			if int64(f.table[j]) != int64(65534) && int64(f.table[j]) != int64(0) {
				gidx = append(gidx, int64(f.table[j]))
			}

		}

		if isUniqueGenome(gidx) && (len(gidx) > 0) {
			idx := gidx[0]

			if !(bacteria_map[idx].Signatures.Has(string(kmer_scanner.Kmer))) {
				tempSet := bacteria_map[idx].Signatures
				tempSet.Add(string(kmer_scanner.Kmer))
				bacteria_map[idx] = &Bacteria{tempSet, bacteria_map[idx].Threshold, bacteria_map[idx].Reported}
			}

			if bacteria_map[idx].ReachThreshold() && (bacteria_map[idx].Reported == false) {
				fmt.Println("Found bacteria ", idx)
				bacteria_map[idx] = &Bacteria{bacteria_map[idx].Signatures, bacteria_map[idx].Threshold, true}
			}
		}
    }	

}

func isUniqueGenome(idx []int64) bool {
	for i := 1; i < len(idx); i++ {
        if idx[i] != idx[0] {
            return false
        }
    }
    return true
}
