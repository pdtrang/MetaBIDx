package ppt_filter

import (
	"bufio"
	"encoding/binary"
	"encoding/gob"
	"encoding/json" // Encoding and Decoding Package
	"fmt"
	"log"
	"os"
	"path"
	"unsafe"
	"io/ioutil"
	"sync"
	"../utils"
	"time"
)

const Unused = uint16(65535)

type Filter struct {
	M            int64
	K            int
	HashFunction []*LinearHash
	table        []uint16
	Gid			 map[uint16]string // map gids and strains/species names
	Gid_header   map[uint16][]string // map gids and sequence headers (for query)
	SeqLength    map[string]int // map of sequence length
	Kmer_pos     map[uint16][]int // map of position of unique kmers in each sequence
	N_phases	 int
	Total_signatures map[uint16]int //map of total signatures of each bacteria
	NumOfLocks	int
	lock map[int]*sync.Mutex
}

//-----------------------------------------------------------------------------
// m: size of hash table.
// k: length of kmers
//-----------------------------------------------------------------------------

func NewFilter(m int64, k int, num_hashes int, n_phases int, nlocks int) *Filter {
	//num_hashes := 2
	f := &Filter{
		M:     m,
		K:     k,
		table: make([]uint16, m),
		Gid: make(map[uint16]string),
		Gid_header: make(map[uint16][]string),
		SeqLength: make(map[string]int),
		Kmer_pos: make(map[uint16][]int),
		N_phases: n_phases,
		Total_signatures: make(map[uint16]int),
		NumOfLocks: nlocks,
	}
	f.HashFunction = make([]*LinearHash, num_hashes)
	for i := 0; i < num_hashes; i++ {
		f.HashFunction[i] = NewLinearHash(m)
		f.HashFunction[i].SetK(k)
	}
	f.lock = make(map[int]*sync.Mutex, f.NumOfLocks)
	for i := 0; i <= f.NumOfLocks; i++ {
		f.lock[i] = new(sync.Mutex)
	}

	return f
}

//-----------------------------------------------------------------------------
func (f *Filter) IsEmpty() bool {
	return len(f.table) == 0
}

//-----------------------------------------------------------------------------
func (f *Filter) InitLocks() {
	f.lock = make(map[int]*sync.Mutex, f.NumOfLocks)
	for i := 0; i <= f.NumOfLocks; i++ {
		f.lock[i] = new(sync.Mutex)
	}
}

//-----------------------------------------------------------------------------
func (f *Filter) Summarize() {
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

func (f *Filter) SetTable(table []uint16) {
	f.table = table
}

func (f *Filter) GetLock(entry int64) int {
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
func (f *Filter) RemoveUnusedKmers(gid uint16, seq []byte, pos_array []int) {

	// fmt.Println("Remove Unused Kmers")
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
			}
		}

		if unique_to_genome {
			fmt.Println("remove unique on main", string(kmer), idx)
			for i := 0; i < len(idx); i++ {
				f.table[idx[i]] = Unused	
			} 
		} 

		if unique_to_genome_rc {
			fmt.Println("remove unique on rc", string(kmer_rc), idx_rc)
			for i := 0; i < len(idx_rc); i++ {
				f.table[idx_rc[i]] = Unused	
			} 
		}

	}

	
}

//-----------------------------------------------------------------------------
func (f *Filter) RemoveAllKmers() {
	for i := int64(0); i < f.M; i++ {
		f.table[i] = Empty
	}
}

//-----------------------------------------------------------------------------
func (f *Filter) SetGidWithMutex(gid uint16, seq []byte, pos_array []int, temp_table []uint16, mutex *sync.Mutex) (int, int) {
	// fmt.Println("Set GID")
	// fmt.Println("SetGid", len(pos_array))
	count := 0
	count_rc := 0
	mutex.Lock()
	for p := 0; p < len(pos_array); p++ {
		// fmt.Println(gid, pos_array[p])
		kmer := seq[pos_array[p] : f.K+pos_array[p]]

		unique_to_genome := true
		idx := make([]int64, 0)
		for i := 0; i < len(f.HashFunction); i++ {
			j := f.HashFunction[i].HashKmer(kmer)
			idx = append(idx, j)
			if f.table[j] != Empty && f.table[j] != gid {
				unique_to_genome = false
				break
			}
		}

		if unique_to_genome {
			count += 1
			for i := 0; i < len(idx); i++ {
				temp_table[idx[i]] = gid	
			}
			// fmt.Println("main", string(kmer), idx)
		}

		unique_to_genome_rc := true
		kmer_rc := []byte(ReverseComplement(string(kmer)))
		idx_rc := make([]int64, 0)
		for i := 0; i < len(f.HashFunction); i++ {
			j := f.HashFunction[i].HashKmer(kmer_rc)
			idx_rc = append(idx_rc, j)
			if f.table[j] != Empty && f.table[j] != gid {
				unique_to_genome_rc = false
				break
			}
		}

		if unique_to_genome_rc {
			count_rc += 1
			for i := 0; i < len(idx_rc); i++ {
				temp_table[idx_rc[i]] = gid	
			}
			// fmt.Println("rc", string(kmer_rc), idx_rc)
		}

	}
	mutex.Unlock()
	return count, count_rc
}

//-----------------------------------------------------------------------------
func (f *Filter) SetGid(gid uint16, seq []byte, pos_array []int, temp_table []uint16) (int, int) {
	// fmt.Println("Set GID")
	// fmt.Println("SetGid", len(pos_array))
	count := 0
	count_rc := 0
	for p := 0; p < len(pos_array); p++ {
		// fmt.Println(gid, pos_array[p])
		kmer := seq[pos_array[p] : f.K+pos_array[p]]

		unique_to_genome := true
		idx := make([]int64, 0)
		for i := 0; i < len(f.HashFunction); i++ {
			j := f.HashFunction[i].HashKmer(kmer)
			idx = append(idx, j)
			if f.table[j] != Empty && f.table[j] != gid {
				unique_to_genome = false
				break
			}
		}

		if unique_to_genome {
			count += 1
			for i := 0; i < len(idx); i++ {
				temp_table[idx[i]] = gid	
			}
			// fmt.Println("main", string(kmer), idx)
		}

		unique_to_genome_rc := true
		kmer_rc := []byte(ReverseComplement(string(kmer)))
		idx_rc := make([]int64, 0)
		for i := 0; i < len(f.HashFunction); i++ {
			j := f.HashFunction[i].HashKmer(kmer_rc)
			idx_rc = append(idx_rc, j)
			if f.table[j] != Empty && f.table[j] != gid {
				unique_to_genome_rc = false
				break
			}
		}

		if unique_to_genome_rc {
			count_rc += 1
			for i := 0; i < len(idx_rc); i++ {
				temp_table[idx_rc[i]] = gid	
			}
			// fmt.Println("rc", string(kmer_rc), idx_rc)
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

func (f *Filter) CountSignature() {
	for i := int64(0); i < f.M; i++ {
		f.Total_signatures[f.table[i]] += 1
	}

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
func _load_kmerpos(fn string) map[string][]int {
    // read file
    f, err := ioutil.ReadFile(fn)
    if err != nil {
      fmt.Print(err)
    }

    var data map[string][]int
    err = json.Unmarshal(f, &data)
    if err != nil {
        fmt.Println("error:", err)
    }
    

    return data
}

//-----------------------------------------------------------------------------
func _load_binary_kmerpos(fn string) map[uint16][]int {
    
    file, err := os.Open(fn)
    if err != nil {
        log.Fatal(err)
    }
    defer file.Close()
    decoder := gob.NewDecoder(file)
    var data map[uint16][]int
    err = decoder.Decode(&data)

    // fmt.Println("kmer pos: ", data)

    return data
}

//-----------------------------------------------------------------------------
func _load_hashfunction(fn string) []*LinearHash {
    // read file
    f, err := ioutil.ReadFile(fn)
    if err != nil {
      fmt.Print(err)
    }

    var data []*LinearHash
    err = json.Unmarshal(f, &data)
    if err != nil {
        fmt.Println("error:", err)
    }
    

    return data
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
func _save_kmerpos_to_json(data map[string][]int, fn string) {

    // Marshal the map into a JSON string.
    // saveData, err := json.Marshal(data)   
    // if err != nil {
    //     fmt.Println(err.Error())
    //     return
    // }
     
    // jsonStr := string(saveData)
    // fmt.Println("The JSON data is:")
    // fmt.Println(jsonStr)

    file, _ := json.MarshalIndent(data, "", " ")
 
    _ = ioutil.WriteFile(fn, file, 0644)

}

//-----------------------------------------------------------------------------
func _save_kmerpos_to_binary(data map[uint16][]int, fn string) {

    file, err := os.Create(fn)
    if err != nil {
        log.Fatal(err)
    }
    defer file.Close()
    encoder := gob.NewEncoder(file)
    encoder.Encode(data)

}

//-----------------------------------------------------------------------------
func _save_hashfunction_to_json(data []*LinearHash, fn string) {

    // Marshal the map into a JSON string.
    // saveData, err := json.Marshal(data)   
    // if err != nil {
    //     fmt.Println(err.Error())
    //     return
    // }
     
    // jsonStr := string(saveData)
    // fmt.Println("The JSON data is:")
    // fmt.Println(jsonStr)

    file, _ := json.MarshalIndent(data, "", " ")
 
    _ = ioutil.WriteFile(fn, file, 0644)

}

//-----------------------------------------------------------------------------
func (f *Filter) SaveFilterGob(fn string) {
	for i := range f.HashFunction {
		f.HashFunction[i] = ResetLinearHash(f.HashFunction[i], f.K)
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
	// _save_kmerpos_to_json(f.Kmer_pos, path.Join(fn+".json"))
	_save_kmerpos_to_binary(f.Kmer_pos, path.Join(fn+"_kmerpos.bin"))
	// _save_hashfunction_to_json(f.HashFunction, path.Join(fn+"_hf.json"))
}

//-----------------------------------------------------------------------------
func (f *Filter) SaveReducedFilter(fn string) {
	f.SaveFilterGob(fn)
	_save_table_alone(f.table, path.Join(fn+".table"))
	// _save_kmers_bases_to_json(f.Kmers_bases, path.Join(fn+"_bases.json"))
}

func LoadReducedFilter(fn string) * Filter {
	filter := LoadFilterGob(fn)
	filter.table = _load_table_alone(fn+".table", filter.M)
	// filter.Kmers_bases = _load_kmers_bases(fn+"_bases.json")
	return filter
}


//-----------------------------------------------------------------------------
// load the table
func Load(fn string) *Filter {
	defer utils.TimeConsume(time.Now(), "Load filter: ")
	filter := LoadFilterGob(fn)
	filter.table = _load_table_alone(fn+".table", filter.M)
	return filter
}

//-----------------------------------------------------------------------------
// no need to load the table, just initialize it
func LoadFilter(fn string) * Filter {
	filter := LoadFilterGob(fn)
	// filter.table = make([]uint16, filter.M)
	filter.table = _load_table_alone(fn+".table", filter.M)
	filter.InitLocks()
	filter.Gid_header = make(map[uint16][]string)
	// filter.Kmer_pos = _load_kmerpos(fn+".json")
	filter.Kmer_pos = _load_binary_kmerpos(fn+"_kmerpos.bin")
	// filter.HashFunction = _load_hashfunction(fn+"_hf.json")
	return filter
}

//-----------------------------------------------------------------------------
func (f *Filter) AddLengthInfo(data map[string]int) {
	f.SeqLength = data
}

//-----------------------------------------------------------------------------
func (f *Filter) PrintHashFunction() {
	for i := range f.HashFunction {
		fmt.Println("Function ", i, ": ", f.HashFunction[i])
	}
}

