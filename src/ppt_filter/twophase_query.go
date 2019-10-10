package ppt_filter

import (
	"time"
	"fmt"
	"os"
)

func (f *Filter) TwoPhaseQuery(read_1 []byte, read_2 []byte, bacteria_map map[uint16]*Bacteria, start_time time.Time, strategy string, analysis bool, analysis_fi *os.File, header_1 string, header_2 string, genome_info map[string]string) int {

	if strategy == "majority" {
		return f.TwoPhaseMajorityQuery(read_1, read_2, bacteria_map, start_time, analysis, analysis_fi)
	} else if strategy == "one_hit" {
		return f.TwoPhaseOneHitQuery(read_1, read_2, bacteria_map, start_time, analysis, analysis_fi)
	} else {
		return f.TwoPhaseOneOrNothingQuery(read_1, read_2, bacteria_map, start_time, analysis, analysis_fi, header_1, header_2, genome_info)
	}

}

//////////////////////////////////////////////////////////////
// One Hit
//////////////////////////////////////////////////////////////
func (f *Filter) TwoPhaseOneHitQuery(read_1 []byte, read_2 []byte, bacteria_map map[uint16]*Bacteria, start_time time.Time, analysis bool, analysis_fi *os.File) int {

	idx, is_valid_gid, kmer := f.TwoPhaseOneHitQueryRead(read_1)

	if !is_valid_gid && idx == uint16(0) {
		idx, is_valid_gid, kmer = f.TwoPhaseOneHitQueryRead(read_2)
	} 

	if is_valid_gid {
		if analysis == true {
			_, err := analysis_fi.WriteString(string(read_1)+ "," +string(read_2)+","+string(idx)+"\n")
			if err != nil {
				fmt.Println(err)
				analysis_fi.Close()
			}	
		}	
	} else {
		if analysis == true {
			_, err := analysis_fi.WriteString(string(read_1)+ "," +string(read_2)+",NA\n")
			if err != nil {
				fmt.Println(err)
				analysis_fi.Close()
			}	
		}
	}
	

	if is_valid_gid && idx != uint16(0) {
		if (bacteria_map[idx].Reported == false) {
			signatures := make([]int64, 0)

			for i := 0; i < len(f.HashFunction); i++ {
				signatures = append(signatures, f.HashFunction[i].HashKmer(kmer))

			}

			return SaveSignatures(f, signatures, idx, bacteria_map, start_time)	
		} else {
			 return 0
		}
	} else {
		return 0
	}
	
}

func (f *Filter) TwoPhaseOneHitQueryRead(read []byte) (uint16, bool, []byte) {
	kmer_scanner := NewKmerScanner(read, f.K)

	for kmer_scanner.ScanOneStrand() {
		kmer_gid, is_unique_kmer := f.TwoPhasesQueryHashKmer(kmer_scanner.Kmer, kmer_scanner.IsFirstKmer)  
		
		if is_unique_kmer {
			return kmer_gid, true, kmer_scanner.Kmer
		}
	}

	return uint16(0), false, kmer_scanner.Kmer
} 

//////////////////////////////////////////////////////////////
// Majority
//////////////////////////////////////////////////////////////
func (f *Filter) TwoPhaseMajorityQuery(read_1 []byte, read_2 []byte, bacteria_map map[uint16]*Bacteria, start_time time.Time, analysis bool, analysis_fi *os.File) int {
	gidx := make(map[uint16][][]byte) // map to keep all the hit kmers for each genome

	f.TwoPhasesMajorityQueryRead(read_1, gidx)
	f.TwoPhasesMajorityQueryRead(read_2, gidx)

	idx := FindMajority(gidx)	

	if idx != uint16(0) {

		if analysis == true {
			_, err := analysis_fi.WriteString(string(read_1)+ "," +string(read_2)+","+string(idx)+"\n")
			if err != nil {
				fmt.Println(err)
				analysis_fi.Close()
			}
		}
		

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

			return SaveSignatures(f, signatures, idx, bacteria_map, start_time)
		} else {
			return 0
		}
	} else {
		if analysis == true {
			_, err := analysis_fi.WriteString(string(read_1)+ "," +string(read_2)+",NA\n")
			if err != nil {
				fmt.Println(err)
				analysis_fi.Close()
			}	
		}
		
		return 0
	}
}


func (f *Filter) TwoPhasesMajorityQueryRead(read []byte, gidx map[uint16][][]byte) {
	kmer_scanner := NewKmerScanner(read, f.K)

	for kmer_scanner.ScanOneStrand() {
		kmer_gid, is_unique_kmer := f.TwoPhasesQueryHashKmer(kmer_scanner.Kmer, kmer_scanner.IsFirstKmer)  
		
		if is_unique_kmer {
			gidx[kmer_gid] = append(gidx[kmer_gid], kmer_scanner.Kmer)
		}
	}

}

//////////////////////////////////////////////////////////////
// One Or Nothing
//////////////////////////////////////////////////////////////
func (f *Filter) TwoPhaseOneOrNothingQuery(read_1 []byte, read_2 []byte, bacteria_map map[uint16]*Bacteria, start_time time.Time, analysis bool, analysis_fi *os.File, header_1 string, header_2 string, genome_info map[string]string) int {
	kmers := make([][]byte, 0)
	idx := uint16(0)

	idx, is_valid_gid, kmer := f.TwoPhasesOONQueryRead(read_1, &kmers, idx)
	if is_valid_gid {
		idx, is_valid_gid, kmer = f.TwoPhasesOONQueryRead(read_2, &kmers, idx)	
	} else {
		return 0
	}
	
	// if there is only one gid and that gid is not 0
	if is_valid_gid && idx != uint16(0) {
		PrintOnlineResult(f, idx, read_1, read_2, kmer, bacteria_map, header_1, header_2, genome_info)
		if analysis == true {
			_, err := analysis_fi.WriteString(string(read_1)+ "," +string(read_2)+","+string(idx)+"\n")
			if err != nil {
				fmt.Println(err)
				analysis_fi.Close()
			}	
		}

		// if bacteria idx has been reported,
		// there is no need to count the signatures
		// if bacteria idx is not reported, 
		// continue to save and count the signatures
		if (bacteria_map[idx].Reported == false) {
			signatures := make([]int64, 0)

			for j := 0; j < len(kmers); j++ {
				for i := 0; i < len(f.HashFunction); i++ {
					signatures = append(signatures, f.HashFunction[i].HashKmer(kmers[j]))

				}
			}

			return SaveSignatures(f, signatures, idx, bacteria_map, start_time)	
		} else {
			return 0
		}
	} else {
		if analysis == true {
			_, err := analysis_fi.WriteString(string(read_1)+ "," +string(read_2)+",NA\n")
			if err != nil {
				fmt.Println(err)
				analysis_fi.Close()
			}	
		}
		
	}

	return 0

}


func (f *Filter) TwoPhasesOONQueryRead(read []byte, kmers *[][]byte, idx uint16) (uint16, bool, []byte) {
	kmer_scanner := NewKmerScanner(read, f.K)

	for kmer_scanner.ScanOneStrand() {
		kmer_gid, is_unique_kmer := f.TwoPhasesQueryHashKmer(kmer_scanner.Kmer, kmer_scanner.IsFirstKmer)  
		
		if is_unique_kmer {

			// if it is the first gid queried, or
			// if the current id is same as the queried gid
			if (idx != uint16(0) && kmer_gid == idx) || (idx == uint16(0) && kmer_gid != uint16(0)) {
				*kmers = append(*kmers, kmer_scanner.Kmer)	
				return kmer_gid, true, kmer_scanner.Kmer	

			} else {
				return uint16(0), false, kmer_scanner.Kmer
			}
			
		} 
		
	}	

	// no kmers get hit in the query
	return uint16(0), true, kmer_scanner.Kmer
}


func (f *Filter) TwoPhasesQueryHashKmer(kmer []byte, is_first_kmer bool) (uint16, bool) {

	j := f.HashFunction[0].SlidingHashKmer(kmer, is_first_kmer)

	if f.table[j] == Dirty || f.table[j] == Empty {
		return uint16(0), false
	}
 
	
	return f.table[j], true
	
}