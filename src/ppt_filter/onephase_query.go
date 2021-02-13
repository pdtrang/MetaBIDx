package ppt_filter

import (
	"time"
	//"os"
)


func (f *Filter) OnePhaseQuery(read_1 []byte, read_2 []byte, bacteria_map map[uint16]*Bacteria, start_time time.Time, strategy string) int {
	if strategy == "majority" {
		return f.OnePhaseMajorityQuery(read_1, read_2, bacteria_map, start_time)
	} else {
		return f.OnePhaseOneOrNothingQuery(read_1, read_2, bacteria_map, start_time)
	}	
}

func (f *Filter) OnePhaseMajorityQuery(read_1 []byte, read_2 []byte, bacteria_map map[uint16]*Bacteria, start_time time.Time) int {
	gidx := make(map[uint16][][]byte)

	f.OnePhaseMajorityQueryRead(read_1, gidx)
	f.OnePhaseMajorityQueryRead(read_2, gidx)

	idx := FindMajority(gidx)	

	if idx != uint16(0) {
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
		return 0
	}
}

func (f *Filter) OnePhaseMajorityQueryRead(read []byte, gidx map[uint16][][]byte) {
	kmer_scanner := NewKmerScanner(read, f.K)

	for kmer_scanner.ScanOneStrand() {
		kmer_gid, is_unique_kmer := f.OnePhaseQueryHashKmer(kmer_scanner.Kmer, kmer_scanner.IsFirstKmer)  
		
		if is_unique_kmer {
			gidx[kmer_gid] = append(gidx[kmer_gid], kmer_scanner.Kmer)
		}
	}
}

func (f *Filter) OnePhaseOneOrNothingQuery(read_1 []byte, read_2 []byte, bacteria_map map[uint16]*Bacteria, start_time time.Time) int {
	kmers := make([][]byte, 0)
	idx := uint16(0)

	idx, is_valid_gid := f.OnePhaseOONQueryRead(read_1, &kmers, idx)
	if is_valid_gid {
		idx, is_valid_gid = f.OnePhaseOONQueryRead(read_2, &kmers, idx)	
	} else {
		return 0
	}
	
	if is_valid_gid && idx != uint16(0) {

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
	}

	return 0


}

func (f *Filter) OnePhaseOONQueryRead(read []byte, kmers *[][]byte, idx uint16) (uint16, bool) {
	kmer_scanner := NewKmerScanner(read, f.K)

	for kmer_scanner.ScanOneStrand() {
		kmer_gid, is_unique_kmer := f.OnePhaseQueryHashKmer(kmer_scanner.Kmer, kmer_scanner.IsFirstKmer)  
		
		if is_unique_kmer {

			if (idx != uint16(0) && kmer_gid == idx) || (idx == uint16(0) && kmer_gid != uint16(0)) {
				*kmers = append(*kmers, kmer_scanner.Kmer)	
				return kmer_gid, true	
			} else {
				return uint16(0), false
			}
			
		} 
		
	}	

	return uint16(0), true
}


func (f *Filter) OnePhaseQueryHashKmer(kmer []byte, is_first_kmer bool) (uint16, bool) {
	idx := uint16(0)
	for i := 0; i < len(f.HashFunction); i++ {
		j := f.HashFunction[i].SlidingHashKmer(kmer, is_first_kmer)

		if f.table[j] == Dirty || f.table[j] == Empty {
			return uint16(0), false
		}
		idx = f.table[j]

	}

	return idx, true
}

