package ppt_filter

import (
	"time"
	// "fmt"
	//"os"
	// "sync"
)

const Dirty = uint16(65534)

//-----------------------------------------------------------------------------
// Majority
//-----------------------------------------------------------------------------
func FindMajority_GID(gidx map[uint16]int) uint16 {
	// fmt.Println("Find majority GID", gidx)
	maxCount := 0
	index := uint16(0)
	total_count := 0 
	for key, val := range gidx {
		if val > maxCount {
			maxCount = val
			index = key
		}	
		total_count += val
	}
	
	if maxCount > total_count/2 {
		// fmt.Println(maxCount, total_count)
		return index
	}
	
	return uint16(0)
}

func (f *FilterInt64) OnePhaseMajorityQuery(read_1 string, read_2 string, qual1 string, qual2 string, start_time time.Time, kmer_qual_threshold int) string {
	// defer Timer()()
	gidx := make(map[uint16]int) // map to keep all the hit kmers for each genome

	f.OnePhaseMajorityQueryRead(read_1, qual1, gidx, kmer_qual_threshold)

	if len(read_2) > 0 {
		f.OnePhaseMajorityQueryRead(read_2, qual2, gidx, kmer_qual_threshold)
	}
		
	idx := FindMajority_GID(gidx)	

	if idx != uint16(0) {
		return f.Gid[idx]
	} else {
		return "unclassified"
	}
}

func (f *FilterInt64) OnePhaseMajorityQueryRead(read string, qual string, gidx map[uint16]int, kmer_qual_threshold int) {
	if len(qual) != 0 {

		kmer_gid := uint16(0)
		is_valid_kmer := false
		for i := 0; i <= (len(read) - f.K); i++ {

			kmer_gid, is_valid_kmer = f.OnePhaseQueryHashKmer(read, qual, i, kmer_qual_threshold)	

			if is_valid_kmer {
				gidx[kmer_gid] += 1
				if gidx[kmer_gid] > (len(read)/2) {
					break
				}
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

func (f *FilterInt64) OnePhaseQueryHashKmer(read string, qual string, start int, kmer_qual_threshold int) (uint16, bool) {	
	gid_map := make(map[uint16]int)
	for i := 0; i < len(f.HashFunction); i++ {
		j := f.HashFunction[i].HashKmerInt64(read, qual, f.K, start, kmer_qual_threshold)
		
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
