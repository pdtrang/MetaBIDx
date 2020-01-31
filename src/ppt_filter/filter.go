//-----------------------------------------------------------------------------
// Author: Vinhthuy Phan, 2018
//-----------------------------------------------------------------------------
package ppt_filter

import (
	"bufio"
	"encoding/binary"
	"encoding/gob"
	"fmt"
	"log"
	"os"
	"path"
	"unsafe"
)

const Unused = uint16(65535)

type Filter struct {
	M            int64
	K            int
	HashFunction []*LinearHash
	table        []uint16
	Gid			 map[uint16]string
	Kmer_pos     map[string][]int
	N_phases	 int
}

//-----------------------------------------------------------------------------
// m: size of hash table.
// k: length of kmers
//-----------------------------------------------------------------------------

func NewFilter(m int64, k int, num_hashes int, n_phases int) *Filter {
	//num_hashes := 2
	f := &Filter{
		M:     m,
		K:     k,
		table: make([]uint16, m),
		Gid: make(map[uint16]string),
		Kmer_pos: make(map[string][]int),
		N_phases: n_phases,
	}
	f.HashFunction = make([]*LinearHash, num_hashes)
	for i := 0; i < num_hashes; i++ {
		f.HashFunction[i] = NewLinearHash(m)
		f.HashFunction[i].SetK(k)
	}
	return f
}

//-----------------------------------------------------------------------------
func (f *Filter) IsEmpty() bool {
	return len(f.table) == 0
}

//-----------------------------------------------------------------------------
func (f *Filter) Summarize() {
	fmt.Println("Number of hash functions: ", len(f.HashFunction))
	fmt.Println("Kmer length:              ", f.K)
	fmt.Println("Table size:               ", f.M)
	count := make(map[uint16]int)
	for i := int64(0); i < f.M; i++ {
		count[f.table[i]]++
	}
	for k, v := range count {
		fmt.Printf("%d\t%d\n", k, v)
	}

	for header, pos := range f.Kmer_pos {
		fmt.Println(header, pos)
	}
}

//-----------------------------------------------------------------------------
func (f *Filter) RemoveUnusedKmers(gid uint16, seq []byte, pos_array []int) {

	// fmt.Println("pos array:", pos_array)
	// keep the first kmer, start from i = 1
	for p := 1; p < len(pos_array); p++ {
		kmer := seq[pos_array[p] : f.K+pos_array[p]]

		unique_to_genome := true
		idx := make([]int64, 0)
		for i := 0; i < len(f.HashFunction); i++ {
			j := f.HashFunction[i].HashKmer(kmer)
			idx = append(idx, j)
			if f.table[j] != 65535 && f.table[j] != gid {
				unique_to_genome = false
				break
			}
		}

		unique_to_genome_rc := true
		kmer_rc := []byte(ReverseComplement(string(kmer)))
		idx_rc := make([]int64, 0)
		for i := 0; i < len(f.HashFunction); i++ {
			j := f.HashFunction[i].HashKmer(kmer_rc)
			idx_rc = append(idx_rc, j)
			if f.table[j] != 65535 && f.table[j] != gid {
				unique_to_genome_rc = false
				break
			}
		}

		if unique_to_genome {
			for i := 0; i < len(idx); i++ {
				f.table[idx[i]] = Unused	
			} 
		} 

		if unique_to_genome_rc {
			for i := 0; i < len(idx_rc); i++ {
				f.table[idx_rc[i]] = Unused	
			} 
		}

	}

	
}
//-----------------------------------------------------------------------------
func (f *Filter) SetGid(gid uint16, seq []byte, pos_array []int) (int, int) {
	count := 0
	count_rc := 0
	for p := 0; p < len(pos_array); p++ {
		kmer := seq[pos_array[p] : f.K+pos_array[p]]

		unique_to_genome := true
		idx := make([]int64, 0)
		for i := 0; i < len(f.HashFunction); i++ {
			j := f.HashFunction[i].HashKmer(kmer)
			idx = append(idx, j)
			if f.table[j] != 65535 && f.table[j] != gid {
				unique_to_genome = false
				break
			}
		}

		if unique_to_genome {
			count += 1
			for i := 0; i < len(idx); i++ {
				f.table[idx[i]] = gid	
			}
		}

		unique_to_genome_rc := true
		kmer_rc := []byte(ReverseComplement(string(kmer)))
		idx_rc := make([]int64, 0)
		for i := 0; i < len(f.HashFunction); i++ {
			j := f.HashFunction[i].HashKmer(kmer_rc)
			idx_rc = append(idx_rc, j)
			if f.table[j] != 65535 && f.table[j] != gid {
				unique_to_genome_rc = false
				break
			}
		}

		if unique_to_genome_rc {
			count_rc += 1
			for i := 0; i < len(idx_rc); i++ {
				f.table[idx_rc[i]] = gid	
			}
		}

	}
	return count, count_rc
}


//-----------------------------------------------------------------------------
func (f *Filter) GetNumberOfUniqueKmers() {
	fmt.Println("Number of hash functions: ", len(f.HashFunction))
	fmt.Println("Kmer length:              ", f.K)
	fmt.Println("Table size:               ", f.M)
	count := make(map[uint16]int)
	for i := int64(0); i < f.M; i++ {
		count[f.table[i]]++
	}
	for k, v := range count {
		fmt.Printf("%s,%d,%.1f\n", f.Gid[k], v, float64(v) / float64(len(f.HashFunction)))
	}
}

func (f *Filter) CountSignature()  map[uint16]int {
	count := make(map[uint16]int)
	for i := int64(0); i < f.M; i++ {
		count[f.table[i]]++
	}

	return count
}

//-----------------------------------------------------------------------------
func (f *Filter) Show() {
	fmt.Println("Number of hash functions: ", len(f.HashFunction))
	fmt.Println("Kmer length:              ", f.K)
	for i := 0; i < len(f.HashFunction); i++ {
		fmt.Println("Hash function", i)
		f.HashFunction[i].Show()
	}
	fmt.Println("Filter table")
	for i := int64(0); i < f.M; i++ {
		fmt.Printf("%d:%d", i, f.table[i])
	}
	fmt.Println("")

}

//-----------------------------------------------------------------------------
func _load_table_alone(filename string, length int64) []uint16 {
	f, err := os.Open(filename)
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()
	numBytes := uint(unsafe.Sizeof(uint16(0)))
	v := make([]uint16, length)

	scanner := bufio.NewScanner(f)
	scanner.Split(bufio.ScanBytes)
	for i, b := 0, uint(0); scanner.Scan(); b++ {
		if b == numBytes {
			b, i = 0, i+1
		}
		v[i] += uint16(scanner.Bytes()[0]) << (b * 8)
	}
	return v
}

//-----------------------------------------------------------------------------
func _save_table_alone(s []uint16, filename string) {
	f, err := os.Create(filename)
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	err = binary.Write(w, binary.LittleEndian, s)
	if err != nil {
		log.Fatal(err)
	}
	w.Flush()
}

//-----------------------------------------------------------------------------
func (f *Filter) SaveFilterGob(fn string) {
	file, err := os.Create(fn)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()
	encoder := gob.NewEncoder(file)
	encoder.Encode(f)
}

//-----------------------------------------------------------------------------
func LoadFilterGob(fn string) *Filter {
	file, err := os.Open(fn)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()
	decoder := gob.NewDecoder(file)
	filter := &Filter{}
	err = decoder.Decode(filter)
	return filter
}

//-----------------------------------------------------------------------------
func (f *Filter) Save(fn string) {
	f.SaveFilterGob(fn)
	_save_table_alone(f.table, path.Join(fn+".table"))
}

//-----------------------------------------------------------------------------
func Load(fn string) *Filter {
	filter := LoadFilterGob(fn)
	filter.table = _load_table_alone(fn+".table", filter.M)
	return filter
}

//-----------------------------------------------------------------------------
func LoadFilter(fn string) * Filter {
	filter := LoadFilterGob(fn)
	// filter.N_phases = n_phases
	filter.table = make([]uint16, filter.M)
	filter.Gid = make(map[uint16]string)
	// fmt.Println(filter)
	return filter
}

//-----------------------------------------------------------------------------
func (f *Filter) PrintHashFunction() {
	for i := range f.HashFunction {
		fmt.Println("Function ", i, ": ", f.HashFunction[i])
	}
}

