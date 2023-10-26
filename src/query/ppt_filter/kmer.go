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
	Kmer_loc      int // current location of Kmer
	Kmer        []byte
	Kmer_rc 	[]byte
	Kmer_qual	[]byte
	K           int
	I           int
	SWindow		int
	WindowPos	int
	IsFirstKmer bool
	Restarted   bool // when encountered non A,C,G,T character, must compute kmer
	IsPrimary   bool
}

//-----------------------------------------------------------------------------
func NewKmerScanner(seq []byte, k int) *KmerScanner {
	return &KmerScanner{
		Seq:         seq,
		Kmer_loc:      0,
		K:           k,
		I:           0,
		SWindow: 	 1,
		WindowPos:	 0,
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
		Kmer_loc:      0,
		K:           k,
		I:           0,
		SWindow: 	 1,
		WindowPos:	 0,
		IsFirstKmer: true,
		Restarted:   false,
		IsPrimary:   true,
	}
}

//-----------------------------------------------------------------------------
func NewKmerScannerSkip(seq []byte, k int, swindow int) *KmerScanner {
	return &KmerScanner{
		Seq:         seq,
		K:           k,
		I:           0,
		SWindow: 	 swindow,
		WindowPos: 	 0,
		IsFirstKmer: true,
		Restarted:   false,
		IsPrimary:   true,
	}
}

//-----------------------------------------------------------------------------
func NewKmerScannerAtIndex(seq []byte, k int, swindow int, index int) *KmerScanner {
	return &KmerScanner{
		Seq:         seq,
		K:           k,
		I:           index,
		SWindow: 	 swindow,
		WindowPos: 	 0,
		IsFirstKmer: true,
		Restarted:   false,
		IsPrimary:   true,
	}
}

func (s *KmerScanner) Scan() bool {
	return s.ScanBothStrands()
}


//-----------------------------------------------------------------------------
// (1) Scan k-mers from the primary strand from left to right, then
// (2) Scan k-mers from the complementary strand from right to left.
// Skip k-mers that contain characters other than A, C, G, T.
//-----------------------------------------------------------------------------
func (s *KmerScanner) ScanBothStrands() bool {
	if s.IsPrimary {
		if s.I >= len(s.Seq)-s.K+1 || s.K > len(s.Seq) {
			s.I = len(s.Seq) - s.K
			s.IsFirstKmer = true
			s.Restarted = false
			s.IsPrimary = false
			// do not return false because we need to go to complementary strand.
		} else {
			if s.I == 0 || s.Restarted {
				s.IsFirstKmer = true
				s.Restarted = false
			} else {
				s.IsFirstKmer = false
			}
			for i := s.I; i < s.K+s.I; i++ {
				if s.Seq[i] != 'A' && s.Seq[i] != 'C' && s.Seq[i] != 'G' && s.Seq[i] != 'T' {
					s.I = i + 1
					s.Restarted = true
					return s.ScanBothStrands()
				}
			}
			s.Kmer = s.Seq[s.I : s.K+s.I]
			// fmt.Println("Primary", string(s.Kmer), s.I)
			s.Kmer_loc = s.I
			s.I++
			return true
		}
	}
	if s.I < 0 || s.K > len(s.Seq) {
		s.I = 0
		s.IsFirstKmer = false
		s.Restarted = false
		s.IsPrimary = true
		return false
	} else {
		if s.I == len(s.Seq)-s.K || s.Restarted {
			s.IsFirstKmer = true
			s.Restarted = false
		} else {
			s.IsFirstKmer = false
		}
		for i := s.K + s.I - 1; i >= s.I; i-- {
			if s.Seq[i] != 'A' && s.Seq[i] != 'C' && s.Seq[i] != 'G' && s.Seq[i] != 'T' {
				s.I = i - s.K
				s.Restarted = true
				return s.ScanBothStrands()
			}
		}
		s.Kmer = s.ReverseComplement(s.Seq[s.I : s.K+s.I])
		// fmt.Println("Reverse", string(s.Kmer), s.I)
		s.Kmer_loc = s.I
		s.I--
		return true
	}
}

func (s *KmerScanner) ScanOne() bool {
	
	if s.I >= len(s.Seq)-s.K+1 || s.K > len(s.Seq) {
		s.I = len(s.Seq) - s.K
		s.IsFirstKmer = true
		s.Restarted = false
		//s.IsPrimary = false
		return false
		// do not return false because we need to go to complementary strand.
	} else {
		if s.I == 0 || s.Restarted {
			s.IsFirstKmer = true
			s.Restarted = false
		} else {
			s.IsFirstKmer = false
		}
		for i := s.I; i < s.K+s.I; i++ {
			if s.Seq[i] != 'A' && s.Seq[i] != 'C' && s.Seq[i] != 'G' && s.Seq[i] != 'T' {
				s.I = i + 1
				s.Restarted = true
				return s.ScanOne()
			}
		}
		s.Kmer = s.Seq[s.I : s.K+s.I]
		s.Kmer_qual = s.Qual[s.I : s.K+s.I]
		s.Kmer_loc = s.I
		
		s.I++
		return true
	}
	
}

func (s *KmerScanner) ScanBothStrandsModified() bool {
	if s.IsPrimary {
		if s.I >= len(s.Seq)-s.K+1 || s.K > len(s.Seq) {
			s.I = len(s.Seq) - s.K
			s.IsFirstKmer = true
			s.Restarted = false
			// s.IsPrimary = false
			// do not return false because we need to go to complementary strand.
			return false
		} else {
			if s.I == 0 || s.Restarted {
				s.IsFirstKmer = true
				s.Restarted = false
			} else {
				s.IsFirstKmer = false
			}
			for i := s.I; i < s.K+s.I; i++ {
				if s.Seq[i] != 'A' && s.Seq[i] != 'C' && s.Seq[i] != 'G' && s.Seq[i] != 'T' {
					s.I = i + 1
					s.Restarted = true
					return s.ScanBothStrandsModified()
				}
			}
			s.Kmer = s.Seq[s.I : s.K+s.I]
			// fmt.Println("Primary", string(s.Kmer), s.I)
			s.Kmer_loc = s.I
			s.I++
			s.IsPrimary = false
			return true
		}
	}

	s.Kmer_rc = s.ReverseComplement(s.Kmer)
	s.IsPrimary = true
	return true
}

func (s *KmerScanner) ScanOneStrand() bool {

	if s.I >= len(s.Seq)-s.K+1 || s.K > len(s.Seq) {
		s.I = len(s.Seq) - s.K
		s.IsFirstKmer = true
		s.Restarted = false
		return false
	} else {
		if s.I == 0 || s.Restarted {
			s.IsFirstKmer = true
			s.Restarted = false
		} else {
			s.IsFirstKmer = false
		}
		for i := s.I; i < s.K+s.I; i++ {
			if s.Seq[i] != 'A' && s.Seq[i] != 'C' && s.Seq[i] != 'G' && s.Seq[i] != 'T' {
				s.I = i + 1
				s.Restarted = true
				return s.ScanOneStrand()
			}
		}
		s.Kmer = s.Seq[s.I : s.K+s.I]
		s.Kmer_qual = s.Qual[s.I : s.K+s.I]
		s.Kmer_loc = s.I

		s.I++
		return true
	}
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
