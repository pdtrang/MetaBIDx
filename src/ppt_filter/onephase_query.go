package ppt_filter

import (
	"time"
	"fmt"
	//"os"
	// "sync"
)

func (f *Filter) OnePhaseQuery(read_1 []byte, read_2 []byte, bacteria_map map[uint16]*Bacteria, start_time time.Time, strategy string) int {
	return f.OnePhaseMajorityQuery(read_1, read_2, bacteria_map, start_time)

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

func (f *Filter) OnePhaseMajorityQuery(read_1 []byte, read_2 []byte, bacteria_map map[uint16]*Bacteria, start_time time.Time) int {
	gidx := make(map[uint16][][]byte) // map to keep all the hit kmers for each genome

	f.OnePhaseMajorityQueryRead(read_1, gidx)

	if string(read_2) != "" {
		f.OnePhaseMajorityQueryRead(read_2, gidx)
	}
		
	idx := FindMajority_GID(gidx)	

	if idx != uint16(0) {
		fmt.Println("Read ", string(read_1), "|", f.Gid[idx])
		// if bacteria idx has been reported,
		// there is no need to count the signatures
		// if bacteria idx is not reported, 
		// save and count the signatures
		if (bacteria_map[idx].Reported == false) {  
			signatures := make([]int64, 0)
			for j := 0; j < len(gidx[idx]); j++ {
				for i := 0; i < len(f.HashFunction); i++ {
					signatures = append(signatures, f.HashFunction[i].HashKmer(gidx[idx][j]))
					
				}
			}

			return SaveSignatures2(f, signatures, idx, bacteria_map, start_time)
		} else {
			return 0
		}
	} else {
		fmt.Println("Read ", string(read_1), " |unclassified")
		return 0
	}
}


func (f *Filter) OnePhaseMajorityQueryRead(read []byte, gidx map[uint16][][]byte) {
	kmer_scanner := NewKmerScanner(read, f.K)

	kmer_gid := uint16(0)
	is_valid_kmer := false
	for kmer_scanner.ScanOneStrand() {
		
		kmer_gid, is_valid_kmer = f.OnePhaseQueryHashKmer(kmer_scanner.Kmer)	
		

		if is_valid_kmer {
			gidx[kmer_gid] = append(gidx[kmer_gid], kmer_scanner.Kmer)
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
		// fmt.Println("Valid kmer", gid_map, gid, maxCount)
		return gid, true
	} else {
		// fmt.Println("\t\tNot valid kmer", gid_map, gid, maxCount)
		return gid, false
	}
		

}

func (f *Filter) OnePhaseQueryHashKmer(kmer []byte) (uint16, bool) {
	gid_map := make(map[uint16]int)
	for i := 0; i < len(f.HashFunction); i++ {
		j := f.HashFunction[i].HashKmer(kmer)

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
