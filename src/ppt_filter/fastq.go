package ppt_filter

import (
    "bufio"
    "io"
    "log"
)

type FastqScanner struct {
    Header     string
    NextHeader string
    Seq        []byte
    Qual       string
    Finished   bool
    Scanner    *bufio.Scanner
}

//-----------------------------------------------------------------------------
func NewFastqScanner(r io.Reader) *FastqScanner {
    scanner := &FastqScanner{Scanner: bufio.NewScanner(r), Finished: false}
    return scanner
}

//-----------------------------------------------------------------------------
func (s *FastqScanner) Scan() bool {
    if s.Finished {
        return false
    }
    var line []byte
    var flag bool
    // 1. Read Fasta header
    if s.NextHeader == "" {
        for flag = s.Scanner.Scan(); flag; flag = s.Scanner.Scan() {
            line = s.Scanner.Bytes()
            if len(line)==0 { continue }
            if line[0] == '@' {
                s.Header = string(line)
                break
            }
        }
        if err := s.Scanner.Err(); err != nil {
            log.Fatal(err)
        }
        // This happens when the file has no FASTA sequence
        if flag == false {
            return false
        }
    } else {
        s.Header = s.NextHeader
    }

    // 2. Read Fastq sequence
    var seq []byte
    for flag = s.Scanner.Scan(); flag; flag = s.Scanner.Scan() {
        line = s.Scanner.Bytes()
        if len(line)==0 { continue }
        if line[0] == '@' {
            s.NextHeader = string(line)
            break
        } 
        if line[0] == '+' {         
            break
        }       
        seq = line
        
    }
    // s.Seq = string(seq)
    s.Seq = seq
    if err := s.Scanner.Err(); err != nil {
        log.Fatal(err)
    }
    s.Finished = !flag
    return true
}
