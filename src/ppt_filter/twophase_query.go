package ppt_filter

import (
	"time"
)

func (f *Filter) TwoPhaseQuery(read_1 []byte, read_2 []byte, bacteria_map map[uint16]*Bacteria, start_time time.Time, strategy string) int {

	if strategy == "majority" {
		return f.TwoPhaseMajorityQuery(read_1, read_2, bacteria_map, start_time)
	} else {
		return f.TwoPhaseOneOrNothingQuery(read_1, read_2, bacteria_map, start_time)
	}

}

func (f *Filter) TwoPhaseMajorityQuery(read_1 []byte, read_2 []byte, bacteria_map map[uint16]*Bacteria, start_time time.Time) int {
	gidx := make(map[uint16][][]byte)

	f.TwoPhasesMajorityQueryRead(read_1, gidx)
	f.TwoPhasesMajorityQueryRead(read_2, gidx)

	idx := FindMajority(gidx)	

	if idx != uint16(0) {
		for j := 0; j < len(gidx[idx]); j++ {
			for i := 0; i < len(f.HashFunction); i++ {
				SaveSignatures(f, f.HashFunction[i].HashKmer(gidx[idx][j]), idx, bacteria_map, start_time)
			}
		}

		return 1	
	} else {
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


func (f *Filter) TwoPhaseOneOrNothingQuery(read_1 []byte, read_2 []byte, bacteria_map map[uint16]*Bacteria, start_time time.Time) int {
	kmers := make([][]byte, 0)
	idx := uint16(0)

	idx, is_valid := f.TwoPhasesOONQueryRead(read_1, kmers, idx)
	if is_valid {
		idx, is_valid = f.TwoPhasesOONQueryRead(read_2, kmers, idx)	
	} else {
		return 0
	}
		
	if is_valid {
		if idx != uint16(0) {
			for j := 0; j < len(kmers); j++ {
				for i := 0; i < len(f.HashFunction); i++ {
					SaveSignatures(f, f.HashFunction[i].HashKmer(kmers[j]), idx, bacteria_map, start_time)
				}
			}

			return 1	
		} else {
			return 0
		}
	}


}


func (f *Filter) TwoPhaseOONQueryRead(read []byte, kmers [][]byte, idx uint16) (uint16, bool) {
	kmer_scanner := NewKmerScanner(read, f.K)

	for kmer_scanner.ScanOneStrand() {
		kmer_gid, is_unique_kmer := f.TwoPhasesQueryHashKmer(kmer_scanner.Kmer, kmer_scanner.IsFirstKmer)  
		
		if is_unique_kmer {
			if idx != uint16(0) && kmer_gid == idx {
				kmers = append(kmers, kmer_scanner.Kmer)	
				return kmer_gid, true	
			} else if idx == uint16(0) {
				kmers = append(kmers, kmer_scanner.Kmer)
				return kmer_gid, true
			} else {
				return uint16(0), false
			}
			
		} 
		
	}	

	return uint16(0), true
}


func (f *Filter) TwoPhasesQueryHashKmer(kmer []byte, is_first_kmer bool) (uint16, bool) {

	j := f.HashFunction[0].SlidingHashKmer(kmer, is_first_kmer)

	if f.table[j] == Dirty || f.table[j] == Empty {
		return uint16(0), false
	}
 
	
	return f.table[j], true
	
}