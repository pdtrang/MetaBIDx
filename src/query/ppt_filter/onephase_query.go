package ppt_filter

import (
	"time"
	// "fmt"
	//"os"
	// "sync"
)

const Dirty = uint16(65534)

// func (f *Filter) OnePhaseQuery(read_1 []byte, read_2 []byte, qual1 []byte, qual2 []byte, start_time time.Time, strategy string, kmer_qual_threshold int) string {
// 	//StartProfile()
// 	//defer Timer()()
// 	return f.OnePhaseMajorityQuery(read_1, read_2, qual1, qual2, header, start_time, kmer_qual_threshold)
// }

//-----------------------------------------------------------------------------
// Majority
//-----------------------------------------------------------------------------
func FindMajority_GID(gidx map[uint16][][]byte) uint16 {
	// fmt.Println("Find majority GID", gidx)
	maxCount := 0
	index := uint16(0)
	total_count := 0 
	for key, element := range gidx {
		if len(element) > maxCount {
			maxCount = len(element)
			index = key
		}	
		total_count += len(element)
	}
	
	if maxCount > total_count/2 {
		// fmt.Println(maxCount, total_count)
		return index
	}
	
	return uint16(0)
}

func (f *Filter) OnePhaseMajorityQuery(read_1 []byte, read_2 []byte, qual1 []byte, qual2 []byte, start_time time.Time, strategy string, kmer_qual_threshold int) string {
	// defer Timer()()
	//fmt.Println("Read ", header)
	// gidx := make(map[uint16][][]byte) // map to keep all the hit kmers for each genome
	gidx := make(map[uint16]int) // map to keep all the hit kmers for each genome

	f.OnePhaseMajorityQueryRead(read_1, qual1, gidx, kmer_qual_threshold)

	if len(read_2) > 0 {
		f.OnePhaseMajorityQueryRead(read_2, qual2, gidx, kmer_qual_threshold)
	}
		
	// idx := FindMajority_GID(gidx)	
	idx := uint16(0)
	maxCount := 0
	total_count := 0
	for key, val := range(gidx) {
		if val > maxCount {
			maxCount = val
			idx = key
		}
		total_count += val
	}

	if maxCount < total_count/2 {
		idx = uint16(0)
	} 

	if idx != uint16(0) {
		return f.Gid[idx]
	} else {
		return "unclassified"
	}
}

func (f *Filter) OnePhaseMajorityQueryRead(read []byte, qual []byte, gidx map[uint16]int, kmer_qual_threshold int) {
	if len(qual) != 0 {
		// fmt.Println("OnePhaseMajQueryRead - func inputs", " read ", string(read), " qual ", string(qual))

		// fmt.Println("OnePhaseMajQueryRead - before loop ", string(kmer_scanner.Seq), string(kmer_scanner.Qual))
		kmer_gid := uint16(0)
		is_valid_kmer := false
		for i := 0; i <= (len(read) - f.K); i++ {
			// check kmer quality 
			// if !isGoodKmer(read, qual, i, f.K, kmer_qual_threshold){
			// 	continue
			// }

			// fmt.Println("OnePhaseMajQueryRead ", string(read), "   kmer: ", string(kmer_scanner.Kmer), "  kmer_qual: ",string(kmer_scanner.Kmer_qual))
			// fmt.Println("OnePhaseMajQueryRead ", i, string(read[i:i+f.K]))
			// continue query if it is a good kmer
			kmer_gid, is_valid_kmer = f.OnePhaseQueryHashKmer(read, qual, i, kmer_qual_threshold)	

			if is_valid_kmer {
				// gidx[kmer_gid] = append(gidx[kmer_gid], read[i:i+f.K])
				gidx[kmer_gid] += 1
			}
		}
	}

}

func CheckMajorityHashValues(gid_map map[uint16]int, num_hash int) (uint16, bool) {
	gid := uint16(0)
	maxCount := 0
	for key, value := range gid_map {
		if key != Dirty {
			if value == num_hash {
				return key, true
			}

			if value > maxCount {
				maxCount = value
				gid = key
			}
		}
	}

	if maxCount > num_hash/2 {
		return gid, true
	} else {
		return gid, false
	}
		

}

func isGoodKmer(read []byte, read_qual []byte, start int, k int, kmer_qual_threshold int) bool {
	total := 0
	for i := start; i < (start + k); i++ {
		if read[i] != 'A' && read[i] != 'C' && read[i] != 'G' && read[i] != 'T' {
			return false
		}
		r := read_qual[i] - 33
		total += int(r)
	}
	mean_qual := total / k
	if mean_qual < kmer_qual_threshold {
		return false
	}
	return true
}

func (f *Filter) OnePhaseQueryHashKmer(read []byte, qual []byte, start int, kmer_qual_threshold int) (uint16, bool) {
	// fmt.Println("\nOnePhaseQueryHashKmer: ", string(read))
	gid_map := make(map[uint16]int)
	for i := 0; i < len(f.HashFunction); i++ {
		// fmt.Println("HashKmer - kmer: ", string(read[start:start + f.K]))
		j := f.HashFunction[i].HashKmer(read, qual, start, f.K, kmer_qual_threshold)
		if j == int64(-1) {
			return uint16(0), false
		} 

		// is Empty
		if f.table[j] == Empty {
			return uint16(0), false
		}

		if _, ok := gid_map[f.table[j]]; ok {
	    	gid_map[f.table[j]] += 1
		} else {
			gid_map[f.table[j]] = 1
		}
	}

	gid := uint16(0)
	is_valid_kmer := false
	gid, is_valid_kmer = CheckMajorityHashValues(gid_map, len(f.HashFunction))

	return gid, is_valid_kmer
	
}
