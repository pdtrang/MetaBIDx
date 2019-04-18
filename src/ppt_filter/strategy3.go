package ppt_filter

import (
	"fmt"
)

//-----------------------------------------------------------------------------
// For each genome, store all the kmers which are matched.
//-----------------------------------------------------------------------------
func (f *Filter) Query3(read []byte, bacteria_map map[int64]StringSet) {

	kmer_scanner := NewKmerScanner(read, f.K)

	for kmer_scanner.Scan() {

        j := int64(0)
		for i := 0; i < len(f.HashFunction); i++ {
			j = f.HashFunction[i].SlidingHashKmer(kmer_scanner.Kmer, kmer_scanner.IsFirstKmer)
			
			if int64(f.table[j]) != int64(65534) && int64(f.table[j]) != int64(0) {	

				// fmt.Println(f.table[j], string(kmer_scanner.Kmer))						
				if _, ok := bacteria_map[int64(f.table[j])]; ok {
					var tempSet = NewStringSet()
					*tempSet = bacteria_map[int64(f.table[j])]
					tempSet.Add(string(kmer_scanner.Kmer))
					bacteria_map[int64(f.table[j])] = *tempSet
				} else {
					s := NewStringSet()

					s.Add(string(kmer_scanner.Kmer))
					bacteria_map[int64(f.table[j])] = *s
				}
			}
		}
    }	

}

//-----------------------------------------------------------------------------
// Count the total stored kmers
// If the count is larger than the threshold, report the genome id.
//-----------------------------------------------------------------------------
func (f *Filter) PostQuery3(bacteria_map map[int64]StringSet, threshold float64){

	count := make(map[uint16]int)
	for i := int64(0); i < f.M; i++ {
		count[f.table[i]]++
	}
	var t = NewStringSet()
	fmt.Println("gid\ttotal_matched\ttotal_unique")
	for k, v := range count {
		*t = bacteria_map[int64(k)]
		length_t := len(t.Table)

		if float64(length_t) >= ((float64(v) * threshold) / float64(len(f.HashFunction)) ) {
			fmt.Printf("%d\t%d\t\t%d\n", k, length_t, v/len(f.HashFunction))
		}
		// if number of hash functions is larger than 1, kmer is added to set only 1 time,
		// so we need to divide to the len(f.HashFunction)
		// if storing index, there's no need to do the division.
	}


}
