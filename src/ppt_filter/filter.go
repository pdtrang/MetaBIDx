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

type Filter struct {
	M            int64
	K            int
	HashFunction []*LinearHash
	table        []uint16
	Gid			 map[uint16]string
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
		fmt.Printf("%s,%d,%f\n", f.Gid[k], v, float64(v) / float64(len(f.HashFunction)))
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
