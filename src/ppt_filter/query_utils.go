package ppt_filter

import (
	"fmt"
	"time"
	"log"
	"os"
	"../utils"
	"strings"
)

func SaveSignatures(f *Filter, signatures []int64, idx uint16, bacteria_map map[uint16]*Bacteria, start_time time.Time) int {
	bac_found := 0

	for i := 0; i < len(signatures); i++ {
		bacteria_map[idx].AddSignature(signatures[i])

		//if bacteria_map[idx].ReachUpperThreshold() && bacteria_map[idx].Reported == false {
		if bacteria_map[idx].Reported == false {
			elapsed := time.Since(start_time)
			log.Printf("Found [%s], elapsed: %s ", f.Gid[idx], elapsed)
			bacteria_map[idx].Reported = true
			bacteria_map[idx].QueryTime = elapsed
			bac_found = 1
		}
	}

	return bac_found
}

func PrintOnlineResult(f * Filter, idx uint16, read_1 []byte, read_2 []byte, kmer []byte) {
	fmt.Println("-------------------------")
	fmt.Println("Signature found!")
	fmt.Println("Read 1: ", string(read_1))
	fmt.Println("Read 2: ", string(read_2))
	fmt.Println("Signature from filter: ", string(kmer))
	fmt.Println(f.Gid[idx])
	fmt.Println("-------------------------")	
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

    // save all reported bacteria
	// for k, b := range bacteria_map {
 //    	if b.Reported == true {
 //    		s := f.Gid[k] +  "\n"
 //    		_, err := fi.WriteString(s)
	// 	    if err != nil {
	// 	        fmt.Println(err)
	// 	        fi.Close()
	// 	        return
	// 	    }	
 //    	}
        
 //    }    

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