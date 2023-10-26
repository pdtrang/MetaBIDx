package ppt_filter

import (
	"time"
	//"fmt"
	//"os"
	// "sync"
)

func (f *Filter) OnePhaseQuery(read_1 []byte, read_2 []byte, qual1 []byte, qual2 []byte, header string, start_time time.Time, strategy string, kmer_qual_threshold int) string {
	//StartProfile()
	//defer Timer()()
	return f.OnePhaseMajorityQuery(read_1, read_2, qual1, qual2, header, start_time, kmer_qual_threshold)
}

//////////////////////////////////////////////////////////////
// Majority
//////////////////////////////////////////////////////////////
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

func (f *Filter) OnePhaseMajorityQuery(read_1 []byte, read_2 []byte, qual1 []byte, qual2 []byte, header string, start_time time.Time, kmer_qual_threshold int) string {
	//defer Timer()()
	//fmt.Println("Read ", header)
	gidx := make(map[uint16][][]byte) // map to keep all the hit kmers for each genome

	f.OnePhaseMajorityQueryRead(read_1, qual1, gidx, kmer_qual_threshold)

	// if string(read_2) != "" {
	if len(read_2) > 0 {
		f.OnePhaseMajorityQueryRead(read_2, qual2, gidx, kmer_qual_threshold)
	}
		
	idx := FindMajority_GID(gidx)	

	if idx != uint16(0) {
		//fmt.Println("Read ", header, "|", f.Gid[idx])
		//fmt.Println("Read ", string(read_1), "|", f.Gid[idx])
		//fmt.Println("Read ", header, "|", f.Gid[idx])
		//fmt.Println("Read qual ", string(qual1), string(qual2))
                //text := header + " | " + f.Gid[idx] + "\n"
		//WriteResult(out_filename, text)
		// query_results.Add(header, f.Gid[idx])
		return f.Gid[idx]
	} else {
		//fmt.Println("Read ", header, "|", idx," |unclassified")
		//fmt.Println("Read ", string(read_1), "|", idx," |unclassified")
		//fmt.Println("Read qual ", string(qual1), string(qual2))
		// return 0
                //text := header + "| unclassified\n"
		//WriteResult(out_filename, text)
		// query_results.Add(header, "unclassified")
		return "unclassified"
	}
}

func (f *Filter) OnePhaseMajorityQuery_QuickStop(read_1 []byte, read_2 []byte, qual1 []byte, qual2 []byte, header string, start_time time.Time, kmer_qual_threshold int, query_results SafeMap) {
	gidx1 := make(map[uint16][][]byte)
	idx1 := f.OnePhaseMajorityQueryRead_QuickStop(read_1, qual1, gidx1, kmer_qual_threshold)
	gidx2 := make(map[uint16][][]byte)
	idx2 := f.OnePhaseMajorityQueryRead_QuickStop(read_2, qual2, gidx2, kmer_qual_threshold)
	if idx1 != uint16(0) && idx1 == idx2 {
		// return f.Gid[idx1]
	} else {
		// unclassified
	}
}

func (f *Filter) OnePhaseMajorityQueryRead_QuickStop(read []byte, qual []byte, gidx map[uint16][][]byte, kmer_qual_threshold int) uint16 {
	if len(qual) != 0 {
		kmer_scanner := NewKmerScannerQual(read, f.K, qual)
		
		kmer_gid := uint16(0)
		is_valid_kmer := false
		for kmer_scanner.ScanOneStrand() {
			kmer_gid, is_valid_kmer = f.OnePhaseQueryHashKmer(kmer_scanner.Kmer, kmer_scanner.Kmer_qual, kmer_qual_threshold)
			
			if is_valid_kmer {
				gidx[kmer_gid] = append(gidx[kmer_gid], kmer_scanner.Kmer)
			}
			if len(gidx[kmer_gid]) > ((len(read)-f.K+1)/2) {
				return kmer_gid
				break 
			}
		}
	}
	return uint16(0)
}

func (f *Filter) OnePhaseMajorityQueryRead(read []byte, qual []byte, gidx map[uint16][][]byte, kmer_qual_threshold int) {
	if len(qual) != 0 {

		kmer_scanner := NewKmerScannerQual(read, f.K, qual)

		kmer_gid := uint16(0)
		is_valid_kmer := false
		for kmer_scanner.ScanOneStrand() {
			//fmt.Println(string(kmer_scanner.Kmer), string(kmer_scanner.Kmer_qual))
			kmer_gid, is_valid_kmer = f.OnePhaseQueryHashKmer(kmer_scanner.Kmer, kmer_scanner.Kmer_qual, kmer_qual_threshold)	
			

			if is_valid_kmer {
				gidx[kmer_gid] = append(gidx[kmer_gid], kmer_scanner.Kmer)
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
				// fmt.Println("Valid kmer", gid_map, key, value)
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

func isGoodKmer(kmer_qual []byte, kmer_qual_threshold int) bool {
	//runes := []rune(kmer_qual)
	//total := 0
	//for i := 0; i < len(runes); i++ {
	//	r := runes[i] - 33
	//	total += int(r)
	//}
	//mean_qual := total/len(kmer_qual)
	//if mean_qual < kmer_qual_threshold {
	//	return false
	//}
	//return true

	total := 0
	for i := 0; i < len(kmer_qual); i++ {
		r := kmer_qual[i] - 33
		total += int(r)
	}
	mean_qual := total / len(kmer_qual)
	if mean_qual < kmer_qual_threshold {
		return false
	}
	return true
}

func (f *Filter) OnePhaseQueryHashKmer(kmer []byte, kmer_qual []byte, kmer_qual_threshold int) (uint16, bool) {
	// check kmer quality 
	if !isGoodKmer(kmer_qual, kmer_qual_threshold){
		return uint16(0), false
	}

	// continue query if it is a good kmer
	//fmt.Println(string(kmer), string(kmer_qual))
	gid_map := make(map[uint16]int)
	for i := 0; i < len(f.HashFunction); i++ {
		j := f.HashFunction[i].HashKmer(kmer)
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
