package ppt_filter

import (
    "bufio"
    "io"
    "log"
)

type FastqScanner struct {
    Header     []byte
    NextHeader []byte
    Seq        []byte
    Qual       []byte
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
//	defer Timer()()
    if s.Finished {
        return false
    }
    var line []byte
    var flag bool
    // 1. Read Fasta header
    if len(s.NextHeader) == 0 {
        for flag = s.Scanner.Scan(); flag; flag = s.Scanner.Scan() {
            line = s.Scanner.Bytes()
            if len(line)==0 { continue }
            if line[0] == '@' {
                // s.Header = make([]byte, len(line))
                // copy(s.Header, line)
                s.Header = line
                break
            }
        }
        if err := s.Scanner.Err(); err != nil {
            log.Fatal(err)
        }
        // This happens when the file has no sequence
        if flag == false {
            return false
        }
    } else {
        // s.Header = s.NextHeader
        s.Header = make([]byte, len(s.NextHeader))
        copy(s.Header, s.NextHeader)
    }

    // 2. Read Fastq sequence
    var seq []byte
    for flag = s.Scanner.Scan(); flag; flag = s.Scanner.Scan() {
        line = s.Scanner.Bytes()
        if len(line)==0 { continue }
        if line[0] == '@' {
            s.NextHeader = line
            // s.NextHeader = make([]byte, len(line))
            // copy(s.NextHeader, line)
            break
        }
        if line[0] == '+' {
            break
        }
        // seq = line
        copy(seq, line)

    }

    // 3. Read Quality
    var qual []byte
    for flag = s.Scanner.Scan(); flag; flag = s.Scanner.Scan() {
       line = s.Scanner.Bytes()
        if len(line)==0 { continue }
        if line[0] == '@' {
            s.NextHeader = line
            // s.NextHeader = make([]byte, len(line))
            // copy(s.NextHeader, line)
            break
        }
        if line[0] == '+' {
            break
        }

        if line[0] == 'A' || line[0] == 'T' || line[0] == 'G' || line[0] == 'C' || line[0] == 'N' {
            break
        }
        // qual = line
        copy(qual, line)

    }
    s.Seq = make([]byte, len(seq))
    s.Qual = make([]byte, len(qual))
    copy(s.Seq, seq)
    copy(s.Qual, qual)
    // s.Seq = seq
    // s.Qual = qual
    if err := s.Scanner.Err(); err != nil {
        log.Fatal(err)
    }
    s.Finished = !flag
    return true
}
