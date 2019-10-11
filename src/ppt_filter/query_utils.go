package ppt_filter

import (
	"fmt"
	"time"
	"log"
	"os"
	"../utils"
	"strings"
	"encoding/csv"
	"io"
)

func SaveSignatures(f *Filter, signatures []int64, idx uint16, bacteria_map map[uint16]*Bacteria, start_time time.Time) int {
	bac_found := 0

	for i := 0; i < len(signatures); i++ {
		bacteria_map[idx].AddSignature(signatures[i])

		if bacteria_map[idx].ReachUpperThreshold() && bacteria_map[idx].Reported == false {
			elapsed := time.Since(start_time)
			log.Printf("Found [%s], elapsed: %s ", f.Gid[idx], elapsed)
			log.Printf("Threshold: ", bacteria_map[idx].UpperThreshold, bacteria_map[idx].LowerThreshold)
			bacteria_map[idx].Reported = true
			bacteria_map[idx].QueryTime = elapsed
			bac_found = 1
		}
	}

	return bac_found
}

func IsExactSubstring(fasta_file string, substring string) bool {

	fa, err := os.Open(fasta_file)
    if err != nil {
        panic(err)
    }
    fa_scanner := NewFastaScanner(fa)
    rc_substring := RevComp(substring)
    for fa_scanner.Scan() {
    	if strings.Contains(string(fa_scanner.Seq), substring) || strings.Contains(string(fa_scanner.Seq), rc_substring) {
    		return true
    	} 
    }

    return false
}

// analysis utils
func PrintOnlineResult(f *Filter, idx uint16, read_1 []byte, read_2 []byte, kmer []byte, bacteria_map map[uint16]*Bacteria, header_1 string, header_2 string, genome_info map[string]string) {
	// fmt.Println("-------------------------")	

	var header_parts []string
	if strings.Contains(string(read_1), string(kmer)) || strings.Contains(RevComp(string(read_1)), string(kmer)) {
		fmt.Println(header_1)
		fmt.Println("Read 1: ", string(read_1))
		header_parts = strings.Split(header_1, ":")	
	} else if strings.Contains(string(read_2), string(kmer)) || strings.Contains(RevComp(string(read_2)), string(kmer)) {
		fmt.Println(header_2)
		fmt.Println("Read 2: ", string(read_2))
		header_parts = strings.Split(header_2, ":")
	} else {
		fmt.Println("Kmer is not in reads.")
	}

	gid := strings.Replace(header_parts[0], ".fna","",-1)
	gid = strings.Replace(gid, "@", "", -1)
	if val, ok := genome_info[gid]; ok {
		if val != f.Gid[idx] {
			fmt.Println("False Positive")			
		} else {
			fmt.Println("True Positive")
		}
	} else {
		fmt.Println("Can not find", idx, "in genome_info.")
	}

	fmt.Println("Kmer: ", string(kmer))
	fmt.Println("True: ", genome_info[gid])
	true_fasta := "/backup2/dpham2/mende_metagenomics_data/new_groupRef_2/"+genome_info[gid]+".fa"
	fmt.Println("Kmer in True genome:", IsExactSubstring(true_fasta, string(kmer)))
	fmt.Println("Predicted", f.Gid[idx])
	predicted_fasta := "/backup2/dpham2/mende_metagenomics_data/new_groupRef_2/"+f.Gid[idx]+".fa"
	fmt.Println("Kmer in Predicted genome:", IsExactSubstring(predicted_fasta, string(kmer)))
	// fmt.Println("Predicted strain: ", f.Gid[idx])
	// fmt.Println("Number of signature found: ", bacteria_map[idx].Signatures.Size()+1) 
	// fmt.Println("Threshold: ", bacteria_map[idx].UpperThreshold, bacteria_map[idx].LowerThreshold)
	// fmt.Println("-------------------------")	
}

// analysis utils
func RevComp(s string) (string){
	bases := map[string]string{"A": "T", "T":"A", "C": "G", "G":"C"}
	rc_s := ""
	for i := range(s) {
		rc_s = bases[string(s[i])] + rc_s
	}
	// fmt.Println(rc_s)

	return rc_s
}

func LoadGenomeInfo() map[string]string{
	file := "/backup2/dpham2/mende_metagenomics_data/scripts/new_groupRef2_genome_reference.csv"

	csvfile, err := os.Open(file)
	if err != nil {
		panic(err)
	}

	r := csv.NewReader(csvfile)

	genome_info := make(map[string]string)
	// Iterate through the record
	for {
		record, err := r.Read()
		if err == io.EOF{
			break
		}
		if err != nil {
			log.Fatal(err)
		}

		id := record[0]
		id = strings.Replace(id, ".1", "", -1)
		id = strings.Replace(id, ".2", "", -1)
		name := record[1]

		genome_info[id] = name
	}

	return genome_info

}

func ComputeAverageQueryTime(bacteria_map map[uint16]*Bacteria, num_bacteria int) time.Duration {
	
    t := time.Duration(0)
	
	sum := float64(0)
	for _, b := range bacteria_map {
		if b.Reported == true {
			sum += float64(b.QueryTime)
		}
	}

	t = time.Duration(sum/float64(num_bacteria))*time.Nanosecond

	return t
}

func ComputeAverageQueryTimeAll(bacteria_map map[uint16]*Bacteria, start_time time.Time) time.Duration {
	
    t := time.Duration(0)
	
	count := 0
	sum := float64(0)
	for _, b := range bacteria_map {
		if b.ReachLowerThreshold() == true {
			if b.ReachUpperThreshold() == false {
				b.QueryTime = time.Since(start_time)
			}
			sum += float64(b.QueryTime)
			count += 1
		}
		
	}

	t = time.Duration(sum/float64(count))*time.Nanosecond

	return t
}


func SaveQueryResult(f *Filter, bacteria_map map[uint16]*Bacteria, num_bacteria int, fn string, start_time time.Time) {
	fi, err := os.Create(fn)
    if err != nil {
        fmt.Println(err)
        return
    }

    lt_fn := strings.Replace(fn, ".txt", "_lowthreshold.txt", -1)
    lt_fi, err := os.Create(lt_fn)
    if err != nil {
        fmt.Println(err)
        return
    }

	if num_bacteria > 0 {

    	// compute avg query time
    	t := ComputeAverageQueryTime(bacteria_map, num_bacteria)
    	t_all := ComputeAverageQueryTimeAll(bacteria_map, start_time)
    	fmt.Printf("Average query time = %s\n", t)
    	fmt.Printf("Average query time of all bacteria = %s\n", t_all)
    	// s := "# Reported bacteria \n"
    	// _, err = fi.WriteString(s)
	    // if err != nil {
	    //     fmt.Println(err)
	    //     fi.Close()
	    //     return
	    // }

	    for k, b := range bacteria_map {
	    	if b.Reported == true {
	    		s := f.Gid[k] + "\t" + b.QueryTime.String()+ "\n"
	    		_, err := fi.WriteString(s)
			    if err != nil {
			        fmt.Println(err)
			        fi.Close()
			        return
			    }	
	    	}
	        
	    }

	    // Save unreported bacteria
	    if num_bacteria < len(bacteria_map) {
	    	SaveLowThresholdBacteria(f, bacteria_map, start_time, lt_fi)
	    }

	    s := "# Average query time of high-threshold bacteria = " + t.String() + "\n"
    	_, err = fi.WriteString(s)
	    if err != nil {
	        fmt.Println(err)
	        fi.Close()
	        return
	    }

	    s = "# Average query time of all bacteria = " + t_all.String() + "\n"
    	_, err = fi.WriteString(s)
	    if err != nil {
	        fmt.Println(err)
	        fi.Close()
	        return
	    }

    } else {
    	fmt.Println("No bacteria found.")
    	// Save unreported bacteria
	    SaveLowThresholdBacteria(f, bacteria_map, start_time, fi)
    }

    utils.SaveMemUsage(fi)
    
}

func SaveLowThresholdBacteria(f *Filter, bacteria_map map[uint16]*Bacteria, start_time time.Time, fi *os.File) {

	// s := "# Low-threshold bacteria: " + "\n"
	// _, err := fi.WriteString(s)
 //    if err != nil {
 //        fmt.Println(err)
 //        fi.Close()
 //        return
 //    }

	for k, b := range bacteria_map {
    	if b.ReachLowerThreshold() == true && b.ReachUpperThreshold() == false {
    		b.QueryTime = time.Since(start_time)
    		s := f.Gid[k] + "\t" + b.QueryTime.String() + "\n"
    		_, err := fi.WriteString(s)
		    if err != nil {
		        fmt.Println(err)
		        fi.Close()
		        return
		    }	
    	}
    
	}

}

func PrintLowThresholdBacteria(f *Filter, bacteria_map map[uint16]*Bacteria) {
	fmt.Println("Low-threshold bacteria:")
	for i, b := range bacteria_map {
		if b.Reported == false && b.Signatures.Size() > 0 {
			fmt.Println(f.Gid[i], b.Signatures.Size())
		}
	}

}