package ppt_filter

import (
	"fmt"
	"log"
)

//-----------------------------------------------------------------------------
// For each genome, store all the kmers which are matched.
//-----------------------------------------------------------------------------
func (f *Filter) OnlineQuerySingle(read []byte, bacteria_map map[int64]*Bacteria) {

	kmer_scanner := NewKmerScanner(read, f.K)

	for kmer_scanner.Scan() {

		idx := int64(0)
        j := int64(0)
		for i := 0; i < len(f.HashFunction); i++ {
			j = f.HashFunction[i].SlidingHashKmer(kmer_scanner.Kmer, kmer_scanner.IsFirstKmer)
			
			if int64(f.table[j]) != int64(65534) && int64(f.table[j]) != int64(0) {
				idx = int64(f.table[j])
				
				bacteria_map[idx].AddSignature(j)
				
				fmt.Println("Updated signatures:", idx, bacteria_map[idx].Signatures)
				if bacteria_map[idx].ReachThreshold() {
					log.Printf("Found bacteria %d", idx)
				}
			}

		}
		
    }	

}


func (f *Filter) OnlineQueryPair(read_1 []byte, read_2 []byte, bacteria_map map[int64]*Bacteria) {

// 	kmer_scanner := NewKmerScanner(read_1, f.K)
// 	kmer_scanner2 := NewKmerScanner(read_2, f.K)

// 	gidx_1 := make([]int64, 0)
// 	gidx_2 := make([]int64, 0)
// 	for kmer_scanner.Scan() && kmer_scanner2.Scan() {

// 		idx := int64(0)
//         j1 := int64(0)
//         j2 := int64(0)
        
// 		for i := 0; i < len(f.HashFunction); i++ {
// 			j1 = f.HashFunction[i].SlidingHashKmer(kmer_scanner.Kmer, kmer_scanner.IsFirstKmer)
// 			j2 = f.HashFunction[i].SlidingHashKmer(kmer_scanner2.Kmer, kmer_scanner2.IsFirstKmer)

// 			if int64(f.table[j1]) != int64(65534) && int64(f.table[j1]) != int64(0) {
// 				gidx_1 = append(gidx_1, j1)
// 			}

// 			if int64(f.table[j2]) != int64(65534) && int64(f.table[j2]) != int64(0) {
// 				gidx_2 = append(gidx_2, j2)
// 			}


}
