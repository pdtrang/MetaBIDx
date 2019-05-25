package ppt_filter

import (
	"fmt"
	"time"
	"log"
	"os"
	"../utils"
)

func SaveSignatures(f *Filter, signatures []int64, idx uint16, bacteria_map map[uint16]*Bacteria, start_time time.Time) int {
	bac_found := 0

	for i := 0; i < len(signatures); i++ {
		bacteria_map[idx].AddSignature(signatures[i])

		if bacteria_map[idx].ReachUpperThreshold() && bacteria_map[idx].Reported == false {
			elapsed := time.Since(start_time)
			log.Printf("Found [%s], elapsed: %s ", f.Gid[idx], elapsed)
			bacteria_map[idx].Reported = true
			bacteria_map[idx].QueryTime = elapsed
			bac_found = 1
		}
	}

	return bac_found
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

func ComputeAverageQueryTimeAll(bacteria_map map[uint16]*Bacteria) time.Duration {
	
    t := time.Duration(0)
	
	count := 0
	sum := float64(0)
	for _, b := range bacteria_map {
		if b.ReachLowerThreshold() == true {
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

	if num_bacteria > 0 {

    	// compute avg query time
    	t := ComputeAverageQueryTime(bacteria_map, num_bacteria)
    	t_all := ComputeAverageQueryTimeAll(bacteria_map)
    	fmt.Printf("Average query time = %s\n", t)
    	s := "# Reported bacteria \n"
    	_, err = fi.WriteString(s)
	    if err != nil {
	        fmt.Println(err)
	        fi.Close()
	        return
	    }

	    for k, b := range bacteria_map {
	    	if b.Reported == true {
	    		s = f.Gid[k] + "\t" + b.QueryTime.String()+ "\n"
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
	    	SaveLowThresholdBacteria(f, bacteria_map, start_time, fi)
	    }

	    s = "# Average query time of reported bacteria = " + t.String() + " | All bacteria = " + t_all.String() + "\n"
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

	s := "# Low-threshold bacteria: " + "\n"
	_, err := fi.WriteString(s)
    if err != nil {
        fmt.Println(err)
        fi.Close()
        return
    }

	for k, b := range bacteria_map {
    	if b.ReachLowerThreshold() == true && b.ReachUpperThreshold() == false {
    		s = ">" + f.Gid[k] + "\t" + b.QueryTime.String() + "\n"
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