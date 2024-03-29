package ppt_filter

import (
	"time"
	"fmt"
	//"os"
	// "sync"
)

func (f *FilterInt64) TwoPhaseQuery(read_1 []byte, read_2 []byte, start_time time.Time, strategy string, level string) {
	fmt.Println(strategy, level)
	// if strategy == "majority" {
	// 	f.TwoPhaseMajorityQuery(read_1, read_2, start_time)
	// } else if strategy == "one_hit" {
	// 	f.TwoPhaseOneHitQuery(read_1, read_2, start_time)
	// } else {
	// 	f.TwoPhaseOneOrNothingQuery(read_1, read_2, start_time, level)
	// }

}

// //////////////////////////////////////////////////////////////
// // One Hit
// //////////////////////////////////////////////////////////////
// func (f *FilterInt64) TwoPhaseOneHitQuery(read_1 []byte, read_2 []byte, start_time time.Time) {

// 	idx, is_valid_gid, _ := f.TwoPhaseOneHitQueryRead(read_1)
// 	// idx, is_valid_gid, kmer := f.TwoPhaseOneHitQueryRead(read_1)

// 	if string(read_2) != "" {
// 		if !is_valid_gid && idx == uint16(0) {
// 			idx, is_valid_gid, _ = f.TwoPhaseOneHitQueryRead(read_2)
// 		} 
// 	}
	

// 	if is_valid_gid && idx != uint16(0) {
// 		// if (bacteria_map[idx].Reported == false) {
// 		// 	signatures := make([]int64, 0)

// 		// 	for i := 0; i < len(f.HashFunction); i++ {
// 		// 		signatures = append(signatures, f.HashFunction[i].HashKmer(kmer))

// 		// 	}

// 		// 	return SaveSignatures2(f, signatures, idx, bacteria_map, start_time)	
// 		// } else {
// 		// 	 return 0
// 		// }
// 	} else {
// 		// return 0
// 	}
	
// }

// func (f *FilterInt64) TwoPhaseOneHitQueryRead(read []byte) (uint16, bool, []byte) {
// 	kmer_scanner := NewKmerScanner(read, f.K)

// 	kmer_gid := uint16(0)
// 	is_unique_kmer := false
// 	for kmer_scanner.ScanOneStrand() {
		
// 		kmer_gid, is_unique_kmer = f.TwoPhasesQueryHashKmer(kmer_scanner.Kmer, kmer_scanner.IsFirstKmer)	

// 		if is_unique_kmer {
// 			return kmer_gid, true, kmer_scanner.Kmer
// 		}
// 	}

// 	return uint16(0), false, kmer_scanner.Kmer
// } 

// //////////////////////////////////////////////////////////////
// // Majority
// //////////////////////////////////////////////////////////////
// func FindMajority(gidx map[uint16][][]byte) uint16 {
// 	maxCount := 0
// 	index := uint16(0)
// 	total_count := 0 
// 	for key, element := range gidx {
// 		if len(element) > maxCount {
// 			maxCount = len(element)
// 			index = key
// 		}	
// 		total_count += len(element)
// 	}
	
// 	if maxCount > total_count/2 {
// 		// fmt.Println(maxCount, total_count)
// 		return index
// 	}
	
// 	return uint16(0)
// }

// func (f *FilterInt64) TwoPhaseMajorityQuery(read_1 []byte, read_2 []byte, start_time time.Time)  {
// 	gidx := make(map[uint16][][]byte) // map to keep all the hit kmers for each genome

// 	f.TwoPhasesMajorityQueryRead(read_1, gidx)

// 	if string(read_2) != "" {
// 		f.TwoPhasesMajorityQueryRead(read_2, gidx)
// 	}
		
// 	idx := FindMajority(gidx)	

// 	if idx != uint16(0) {
// 		fmt.Println(string(read_1), f.Gid[idx])
// 		// if bacteria idx has been reported,
// 		// there is no need to count the signatures
// 		// if bacteria idx is not reported, 
// 		// save and count the signatures
// 		// if (bacteria_map[idx].Reported == false) {  
// 		// 	signatures := make([]int64, 0)
// 		// 	for j := 0; j < len(gidx[idx]); j++ {
// 		// 		for i := 0; i < len(f.HashFunction); i++ {
// 		// 			signatures = append(signatures, f.HashFunction[i].HashKmer(gidx[idx][j]))
					
// 		// 		}
// 		// 	}

// 		// 	return SaveSignatures2(f, signatures, idx, bacteria_map, start_time)
// 		// } else {
// 		// 	return 0
// 		// }
// 	} else {
// 		fmt.Println(string(read_1), "unclassified")
// 		// return 0
// 	}
// }


// func (f *FilterInt64) TwoPhasesMajorityQueryRead(read []byte, gidx map[uint16][][]byte) {
// 	kmer_scanner := NewKmerScanner(read, f.K)

// 	kmer_gid := uint16(0)
// 	is_unique_kmer := false
// 	for kmer_scanner.ScanOneStrand() {
		
// 		kmer_gid, is_unique_kmer = f.TwoPhasesQueryHashKmer(kmer_scanner.Kmer, kmer_scanner.IsFirstKmer)	
		

// 		if is_unique_kmer {
// 			gidx[kmer_gid] = append(gidx[kmer_gid], kmer_scanner.Kmer)
// 		}
// 	}

// }

// //////////////////////////////////////////////////////////////
// // One Or Nothing
// //////////////////////////////////////////////////////////////
// func (f *FilterInt64) TwoPhaseOneOrNothingQuery(read_1 []byte, read_2 []byte, start_time time.Time, level string) {
// 	kmers := make([][]byte, 0)
// 	idx := uint16(0)

// 	// idx, is_valid_gid, kmer := f.TwoPhasesOONQueryRead(read_1, &kmers, idx)
// 	idx, is_valid_gid, _ := f.TwoPhasesOONQueryRead(read_1, &kmers, idx)

// 	if string(read_2) != "" {
// 		if is_valid_gid {
// 			// idx, is_valid_gid, kmer = f.TwoPhasesOONQueryRead(read_2, &kmers, idx)
// 			idx, is_valid_gid, _ = f.TwoPhasesOONQueryRead(read_2, &kmers, idx)	
// 		} else {
// 			// return 0
// 		}
// 	}	
// 	// if there is only one gid and that gid is not 0
// 	if is_valid_gid && idx != uint16(0) {
// 		// fmt.Println(string(kmer))
// 		// PrintOnlineResult(f, idx, read_1, read_2, kmer, header_1, header_2, genome_info, level)
		

// 		// if bacteria idx has been reported,
// 		// there is no need to count the signatures
// 		// if bacteria idx is not reported, 
// 		// continue to save and count the signatures
// 		// if (bacteria_map[idx].Reported == false) {
// 		// 	signatures := make([]int64, 0)

// 		// 	for j := 0; j < len(kmers); j++ {
// 		// 		for i := 0; i < len(f.HashFunction); i++ {
// 		// 			signatures = append(signatures, f.HashFunction[i].HashKmer(kmers[j]))

// 		// 		}
// 		// 	}

// 		// 	return SaveSignatures2(f, signatures, idx, bacteria_map, start_time)	
// 		// } else {
// 		// 	return 0
// 		// }
// 	}

// 	// return 0

// }


// func (f *FilterInt64) TwoPhasesOONQueryRead(read []byte, kmers *[][]byte, idx uint16) (uint16, bool, []byte) {
// 	kmer_scanner := NewKmerScanner(read, f.K)
// 	for kmer_scanner.ScanOneStrand() {
// 		is_unique_kmer := true
// 		kmer_gid := uint16(0)

		
// 		kmer_gid, is_unique_kmer = f.TwoPhasesQueryHashKmer(kmer_scanner.Kmer, kmer_scanner.IsFirstKmer)	
		
		  
// 		if is_unique_kmer {
// 			// if it is the first gid queried, or
// 			// if the current id is same as the queried gid
// 			if (idx != uint16(0) && kmer_gid == idx) || (idx == uint16(0) && kmer_gid != uint16(0)) {
// 				*kmers = append(*kmers, kmer_scanner.Kmer)	
// 				return kmer_gid, true, kmer_scanner.Kmer	

// 			} else {
// 				return uint16(0), false, kmer_scanner.Kmer
// 			}
			
// 		} 
		
// 	}	

// 	// no kmers get hit in the query
// 	return uint16(0), true, kmer_scanner.Kmer
// }

// func (f *FilterInt64) TwoPhasesQueryHashKmer(kmer []byte, is_first_kmer bool) (uint16, bool) {

// 	idx := uint16(0)
// 	for i := 0; i < len(f.HashFunction); i++ {
// 		j := f.HashFunction[i].HashKmer(kmer)

// 		// is either Dirty or Empty
// 		if f.table[j] == Dirty || f.table[j] == Empty {
// 			return uint16(0), false
// 		}

// 		// get different gid with the current gid
// 		if idx != Empty && f.table[j] != idx {
// 			return uint16(0), false
// 		}

// 		idx = f.table[j]

// 	}

// 	return idx, true
	
// }
