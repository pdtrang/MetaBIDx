package ppt_filter

import (
	"bufio"
	"encoding/binary"
	"encoding/gob"
	// "encoding/json" // Encoding and Decoding Package
	"fmt"
	"log"
	"os"
	"path"
	"unsafe"
	// "io/ioutil"
	"sync"
	// "build_index/utils"
	// "time"
)

const Empty = uint16(0)

const Unused = uint16(65535)

type FilterInt64 struct {
	M            int64
	K            int
	HashFunction []*LinearHashInt64
	table        []uint16
	Gid			 map[uint16]string // map gids and strains/species names
	Gid_header   map[uint16][]string // map gids and sequence headers (for query)
	SeqLength    map[string]int // map of sequence length
	// Kmer_pos     map[string][]int // map of position of unique kmers in each sequence
	N_phases	 int
	Total_signatures map[uint16]int //map of total signatures of each bacteria
	NumOfLocks	int
	lock map[int]*sync.Mutex
}

//-----------------------------------------------------------------------------
// m: size of hash table.
// k: length of kmers
//-----------------------------------------------------------------------------

func NewFilterInt64(m int64, k int, num_hashes int, n_phases int, nlocks int) *FilterInt64 {
	//num_hashes := 2
	f := &FilterInt64{
		M:     m,
		K:     k,
		table: make([]uint16, m),
		Gid: make(map[uint16]string),
		Gid_header: make(map[uint16][]string),
		SeqLength: make(map[string]int),
		// Kmer_pos: make(map[string][]int),
		N_phases: n_phases,
		Total_signatures: make(map[uint16]int),
		NumOfLocks: nlocks,
	}
	f.HashFunction = make([]*LinearHashInt64, num_hashes)
	// fmt.Println("Generate random hash functions")
	for i := 0; i < num_hashes; i++ {
		f.HashFunction[i] = NewLinearHashInt64(m)
		f.HashFunction[i].SetKInt64(k)
	}
	f.lock = make(map[int]*sync.Mutex, f.NumOfLocks)
	for i := 0; i <= f.NumOfLocks; i++ {
		f.lock[i] = new(sync.Mutex)
	}

	return f
}

//-----------------------------------------------------------------------------
func (f *FilterInt64) ResetTotalSignatures(){
	f.Total_signatures = make(map[uint16]int)
}

//-----------------------------------------------------------------------------
func (f *FilterInt64) InitLocks() {
	f.lock = make(map[int]*sync.Mutex, f.NumOfLocks)
	for i := 0; i <= f.NumOfLocks; i++ {
		f.lock[i] = new(sync.Mutex)
	}
}

//-----------------------------------------------------------------------------
func (f *FilterInt64) Summarize() {
	f.CountSignature()
	fmt.Println("Number of hash functions: ", len(f.HashFunction))
	fmt.Println("Kmer length:              ", f.K)
	fmt.Println("Table size:               ", f.M)
	fmt.Println("Gid")
	for key, value := range f.Gid {
	    fmt.Println(key,":",value)
	}
	fmt.Println("Total_signatures")
	for key, value := range f.Total_signatures {
	    fmt.Println(key,":",value)
	}
}


func (f *FilterInt64) GetLock(entry int64) int {
	div := int(f.M) / f.NumOfLocks

	if div == 0 {
		return 0
	}

	mut := int(entry) / div
	if mut > f.NumOfLocks {
		mut = f.NumOfLocks
	}

	return mut
}

//-----------------------------------------------------------------------------
func (f *FilterInt64) CountSignature() {
	f.ResetTotalSignatures()
	for i := int64(0); i < f.M; i++ {
		f.Total_signatures[f.table[i]] += 1
	}

}

//-----------------------------------------------------------------------------
func (f *FilterInt64) Show() {
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
// func _load_hashfunction(fn string) []*LinearHash {
//     // read file
//     f, err := ioutil.ReadFile(fn)
//     if err != nil {
//       fmt.Print(err)
//     }

//     var data []*LinearHash
//     err = json.Unmarshal(f, &data)
//     if err != nil {
//         fmt.Println("error:", err)
//     }
    

//     return data
// }

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
// func _save_kmerpos_to_json(data map[string][]int, fn string) {

//     // Marshal the map into a JSON string.
//     // saveData, err := json.Marshal(data)   
//     // if err != nil {
//     //     fmt.Println(err.Error())
//     //     return
//     // }
     
//     // jsonStr := string(saveData)
//     // fmt.Println("The JSON data is:")
//     // fmt.Println(jsonStr)

//     file, _ := json.MarshalIndent(data, "", " ")
 
//     _ = ioutil.WriteFile(fn, file, 0644)

// }

//-----------------------------------------------------------------------------
// func _save_kmerpos_to_binary(data map[string][]int, fn string) {

//     file, err := os.Create(fn)
//     if err != nil {
//         log.Fatal(err)
//     }
//     defer file.Close()
//     encoder := gob.NewEncoder(file)
//     encoder.Encode(data)

// }

//-----------------------------------------------------------------------------
// func _save_hashfunction_to_json(data []*LinearHash, fn string) {

//     // Marshal the map into a JSON string.
//     // saveData, err := json.Marshal(data)   
//     // if err != nil {
//     //     fmt.Println(err.Error())
//     //     return
//     // }
     
//     // jsonStr := string(saveData)
//     // fmt.Println("The JSON data is:")
//     // fmt.Println(jsonStr)

//     file, _ := json.MarshalIndent(data, "", " ")
 
//     _ = ioutil.WriteFile(fn, file, 0644)

// }

//-----------------------------------------------------------------------------
func (f *FilterInt64) SaveFilterGob(fn string) {
	for i := range f.HashFunction {
		f.HashFunction[i] = ResetLinearHashInt64(f.HashFunction[i], f.K)
	}

	file, err := os.Create(fn)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()
	encoder := gob.NewEncoder(file)
	encoder.Encode(f)
}

//-----------------------------------------------------------------------------
func LoadFilterGob(fn string) *FilterInt64 {
	file, err := os.Open(fn)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()
	decoder := gob.NewDecoder(file)
	filter := &FilterInt64{}
	err = decoder.Decode(filter)
	return filter
}

//-----------------------------------------------------------------------------
func (f *FilterInt64) Save(fn string) {
	f.SaveFilterGob(fn)
	_save_table_alone(f.table, path.Join(fn+".table"))
}

//-----------------------------------------------------------------------------
// load the table
// func Load(fn string) *Filter {
// 	// defer utils.TimeConsume(time.Now(), "Load filter: ")
// 	filter := LoadFilterGob(fn)
// 	filter.table = _load_table_alone(fn+".table", filter.M)
// 	return filter
// }

//-----------------------------------------------------------------------------
// no need to load the table, just initialize it
func LoadFilter(fn string) * FilterInt64 {
	filter := LoadFilterGob(fn)
	// filter.table = make([]uint16, filter.M)
	filter.table = _load_table_alone(fn+".table", filter.M)
	filter.InitLocks()
	filter.Gid_header = make(map[uint16][]string)
	return filter
}

//-----------------------------------------------------------------------------
func (f *FilterInt64) PrintHashFunction() {
	for i := range f.HashFunction {
		fmt.Println("Function ", i, ": ", f.HashFunction[i])
	}
}

