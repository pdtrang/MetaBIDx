package ppt_filter

import (
	"time"
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

	f.OnePhaseQueryRead(read_1, gidx)
	f.OnePhaseQueryRead(read_2, gidx)

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

func (f *Filter) OnePhaseOneOrNothingQuery(read_1 []byte, read_2 []byte, bacteria_map map[uint16]*Bacteria, start_time time.Time) int {
	gidx := make(map[uint16][][]byte)

	f.OnePhaseQueryRead(read_1, gidx)
	f.OnePhaseQueryRead(read_2, gidx)

	idx, is_gid := OneOrNothing(gidx)	

	if is_gid == true {
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

func (f *Filter) OnePhaseQueryRead(read []byte, gidx map[uint16][][]byte) {
	kmer_scanner := NewKmerScanner(read, f.K)

	for kmer_scanner.ScanOneStrand() {
		kmer_gid, is_unique_kmer := f.OnePhaseQueryHashKmer(kmer_scanner.Kmer, kmer_scanner.IsFirstKmer)  
		
		if is_unique_kmer {
			gidx[kmer_gid] = append(gidx[kmer_gid], kmer_scanner.Kmer)
		}
	}
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

