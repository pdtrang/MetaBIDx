package ppt_filter

import (
	// "fmt"
)

//-----------------------------------------------------------------------------
// This scanner is designed to work with the hash functions.
//-----------------------------------------------------------------------------

type KmerScanner struct {
	Seq         []byte
	Qual		[]byte
	Kmer        []byte
	Kmer_qual	[]byte
	K           int
	I           int
	IsFirstKmer bool
	Restarted   bool // when encountered non A,C,G,T character, must compute kmer
	IsPrimary   bool
}

//-----------------------------------------------------------------------------
func NewKmerScanner(seq []byte, k int) *KmerScanner {
	return &KmerScanner{
		Seq:         seq,
		K:           k,
		I:           0,
		IsFirstKmer: true,
		Restarted:   false,
		IsPrimary:   true,
	}
}

//-----------------------------------------------------------------------------
func NewKmerScannerQual(seq []byte, k int, qual []byte) *KmerScanner {
	return &KmerScanner{
		Seq:         seq,
		Qual: 		 qual,
		K:           k,
		I:           0,
		IsPrimary:   true,
	}
}

//-----------------------------------------------------------------------------
func (s *KmerScanner) ScanOneStrand() bool {

	if s.I >= len(s.Seq)-s.K+1 || s.K > len(s.Seq) {
		return false
	} 

	// status := true
	// for i := s.I; i < s.K+s.I; i++ {
	// 	if s.Seq[i] != 'A' && s.Seq[i] != 'C' && s.Seq[i] != 'G' && s.Seq[i] != 'T' {
	// 		status = false
	// 		break
	// 	}
	// }

	// if status {
	// 	s.Kmer = s.Seq[s.I : s.K+s.I]
	// 	s.Kmer_qual = s.Qual[s.I : s.K+s.I]	
	// 	// fmt.Println("ScanOneStrand - ", string(s.Kmer), string(s.Kmer_qual))
	// } else {
	// 	s.Kmer = []byte("")
	// 	s.Kmer_qual = []byte("")
	// }

	s.Kmer = s.Seq[s.I : s.K+s.I]
	s.Kmer_qual = s.Qual[s.I : s.K+s.I]	

	s.I++
	return true

}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Legacy code below
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

func base_to_num(base byte) int {
	if base == 'A' {
		return 0
	} else if base == 'C' {
		return 1
	} else if base == 'G' {
		return 2
	} else if base == 'T' {
		return 3
	} else {
		panic("Unknown character: " + string(base))
	}
}

//-----------------------------------------------------------------------------
func base_to_num_rc(base byte) int {
	if base == 'A' {
		return 3
	} else if base == 'C' {
		return 2
	} else if base == 'G' {
		return 1
	} else if base == 'T' {
		return 0
	} else {
		panic("Unknown character: " + string(base))
	}
}

//-----------------------------------------------------------------------------
// return the decimal value of a k-mer (A=0,C=1,G=2,T=3)
//-----------------------------------------------------------------------------
func kmer_to_dec(kmer []byte) int {
	value := 0
	exp := 1
	for i := len(kmer) - 1; i >= 0; i-- {
		if kmer[i] == 'A' {
			value += 0 * exp
		} else if kmer[i] == 'C' {
			value += 1 * exp
		} else if kmer[i] == 'G' {
			value += 2 * exp
		} else if kmer[i] == 'T' {
			value += 3 * exp
		} else {
			panic("Unknown character: " + string(kmer[i]))
		}
		exp *= 4
	}
	return value
}

//-----------------------------------------------------------------------------
// return the decimal value of the reverse complement of a k-mer (A=0,C=1,G=2,T=3)
//-----------------------------------------------------------------------------
func kmer_to_dec_rc(kmer []byte) int {
	value := 0
	exp := 1
	for i := 0; i < len(kmer); i++ {
		if kmer[i] == 'A' {
			value += 3 * exp
		} else if kmer[i] == 'C' {
			value += 2 * exp
		} else if kmer[i] == 'G' {
			value += 1 * exp
		} else if kmer[i] == 'T' {
			value += 0 * exp
		} else {
			panic("Unknown character: " + string(kmer[i]))
		}
		exp *= 4
	}
	return value
}

//-----------------------------------------------------------------------------
func DecToKmer(x int, K int) string {
	y := make([]byte, K)
	for i := K - 1; i >= 0; i-- {
		base := x % 4
		switch base {
		case 0:
			y[i] = 'A'
		case 1:
			y[i] = 'C'
		case 2:
			y[i] = 'G'
		case 3:
			y[i] = 'T'
		}
		x = (x - base) >> 2
	}
	return string(y)
}

//-----------------------------------------------------------------------------
func (s *KmerScanner) ReverseComplement(dna []byte) []byte {
	r := make([]byte, len(dna))
	var c byte
	for i := 0; i < len(dna); i++ {
		c = dna[len(dna)-i-1]
		if c == 'A' {
			r[i] = 'T'
		} else if c == 'C' {
			r[i] = 'G'
		} else if c == 'G' {
			r[i] = 'C'
		} else if c == 'T' {
			r[i] = 'A'
		} else {
			panic("Unknown character: " + string(c))
		}
	}
	return r
}

//-----------------------------------------------------------------------------
func ReverseComplement(s string) string {
	r := make([]byte, len(s))
	var c byte
	for i := 0; i < len(s); i++ {
		c = s[len(s)-i-1]
		if c == 'A' {
			r[i] = 'T'
		} else if c == 'C' {
			r[i] = 'G'
		} else if c == 'G' {
			r[i] = 'C'
		} else if c == 'T' {
			r[i] = 'A'
		} else {
			panic("Unknown character: " + string(c))
		}
	}
	return string(r)
}
