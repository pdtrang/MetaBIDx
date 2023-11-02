package ppt_filter

import (
	"time"
	"fmt"
	//"os"
	// "sync"
)

const Dirty = uint16(65534)

func (f *Filter) OnePhaseQuery(read_1 []byte, read_2 []byte, qual1 []byte, qual2 []byte, header []byte, start_time time.Time, strategy string, kmer_qual_threshold int) string {
	//StartProfile()
	//defer Timer()()
	return f.OnePhaseMajorityQuery(read_1, read_2, qual1, qual2, header, start_time, kmer_qual_threshold)
}

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

func (f *Filter) OnePhaseMajorityQuery(read_1 []byte, read_2 []byte, qual1 []byte, qual2 []byte, header []byte, start_time time.Time, kmer_qual_threshold int) string {
	//defer Timer()()
	//fmt.Println("Read ", header)
	gidx := make(map[uint16][][]byte) // map to keep all the hit kmers for each genome

	f.OnePhaseMajorityQueryRead(read_1, qual1, gidx, kmer_qual_threshold)

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

func (f *Filter) OnePhaseMajorityQueryRead(read []byte, qual []byte, gidx map[uint16][][]byte, kmer_qual_threshold int) {
	if len(qual) != 0 {
		// fmt.Println("OnePhaseMajQueryRead - func inputs", " read ", string(read), " qual ", string(qual))

		kmer_scanner := NewKmerScannerQual(read, f.K, qual)
		// fmt.Println("OnePhaseMajQueryRead - before loop ", string(kmer_scanner.Seq), string(kmer_scanner.Qual))
		kmer_gid := uint16(0)
		is_valid_kmer := false
		for kmer_scanner.ScanOneStrand() {
			// check kmer quality 
			if !isGoodKmer(kmer_scanner.Kmer_qual, kmer_qual_threshold){
				continue
			}

			// fmt.Println("OnePhaseMajQueryRead ", string(read), "   kmer: ", string(kmer_scanner.Kmer), "  kmer_qual: ",string(kmer_scanner.Kmer_qual))
			// continue query if it is a good kmer
			kmer_gid, is_valid_kmer = f.OnePhaseQueryHashKmer(kmer_scanner.Kmer)	

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

func (f *Filter) OnePhaseQueryHashKmer(kmer []byte) (uint16, bool) {
	gid_map := make(map[uint16]int)
	for i := 0; i < len(f.HashFunction); i++ {
		// fmt.Println("HashKmer - kmer: ", string(kmer))
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
