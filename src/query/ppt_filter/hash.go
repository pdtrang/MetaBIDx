//-----------------------------------------------------------------------------
// Author: Vinhthuy Phan, 2018
//-----------------------------------------------------------------------------
package ppt_filter

import (
    "fmt"
    "math/rand"
    "time"
)

//-----------------------------------------------------------------------------
// Integer hashing:
//  h(x) = ((ax+b) mod p) mod m
// String hashing:
//  h(x) = ((x_0*e^{k-1} + x_1*e^{k-2} + ... + x_{k-2}*e + x_{k-1}*1) mod p) mod m
//-----------------------------------------------------------------------------
type LinearHash struct {
    A               int64 // randomly selected
    B               int64 // randomly selected
    P               int64 // large prime
    M               int64    // universe/hash table size
    K               int
    Base            int64
    Term0           int64
    Term0_rc        int64
    PrevValue       int64
    PrevValue_rc    int64
    Exponents       []int64
}

//-----------------------------------------------------------------------------
func NewLinearHash(m int64) *LinearHash {
    rand.Seed(time.Now().UTC().UnixNano())

    // temporarily: 2^61 - 1
    var p int64 = 2305843009213693951
    // generate random numbers in [2, p-1]
    a := rand.Int63()%(p-2) + 2
    b := rand.Int63()%(p-2) + 2

    return &LinearHash{
        // A: big.NewInt(23),
        // B: big.NewInt(17),
        // P: big.NewInt(97),
        A:              int64(a),
        B:              int64(b),
        P:              int64(p),
        M:              m,
        Term0:          int64(0),
        Term0_rc:       int64(0),
        PrevValue:      int64(0),
        PrevValue_rc:   int64(0),
    }
}

//-----------------------------------------------------------------------------
func ResetLinearHash(linear_hash *LinearHash, k int) *LinearHash {

    return &LinearHash{
        // A: big.NewInt(23),
        // B: big.NewInt(17),
        // P: big.NewInt(97),
        A:              linear_hash.A,
        B:              linear_hash.B,
        P:              linear_hash.P,
        M:              linear_hash.M,
        K:              k,
        Base:           linear_hash.Base,
        Term0:          int64(0),
        Term0_rc:       int64(0),
        PrevValue:      int64(0),
        PrevValue_rc:   int64(0),
        Exponents:      linear_hash.Exponents,
    }
}

//-----------------------------------------------------------------------------
func (h *LinearHash) SetK(k int) {
    h.K = k
    h.Exponents = make([]int64, k)
    h.Base = rand.Int63n(65536-4) + 4
    b := int64(1)
    for i := k - 1; i >= 0; i-- {
        h.Exponents[i] = 1
        h.Exponents[i] = (h.Exponents[i] * b) % h.P
        b = (h.Base * b) % h.P
    }
}


//-----------------------------------------------------------------------------
func (h *LinearHash) ComputeKmer(kmer []byte) int64 {
    if len(kmer) != h.K {
        panic("Unmatched k-mer length")
    }
    var base int64
    value := int64(0)
    for i := 0; i < len(kmer); i++ {
        if kmer[i] == 'A' {
            base = 0
        } else if kmer[i] == 'C' {
            base = 1
        } else if kmer[i] == 'G' {
            base = 2
        } else if kmer[i] == 'T' {
            base = 3
        } else {
            // fmt.Println(string(kmer))
            panic("ComputeKmer" + string(kmer) + "Unknown character: " + string(kmer[i]))
        }
        cur_term := (base * h.Exponents[i]) % h.P
        value = (value + cur_term) % h.P
        if i == 0 {
            h.Term0 = cur_term
        }
    }
    h.PrevValue = value
    return value % h.M
}


//-----------------------------------------------------------------------------
func (h *LinearHash) HashKmer(kmer []byte) int64 {
    // fmt.Println("HashKmer func: ", string(kmer))
    if len(kmer) != h.K {
        panic("Unmatched k-mer length")
    }    

    return h.HashInt64(h.ComputeKmer(kmer))
}

//-----------------------------------------------------------------------------
func (h *LinearHash) HashInt64(x int64) int64 {
    value := (h.A * x + h.B) % h.P
    return value % h.M
}


//-----------------------------------------------------------------------------
func (h *LinearHash) Show() {
    fmt.Println("A: ", h.A)
    fmt.Println("B: ", h.B)
    fmt.Println("P: ", h.P)
    fmt.Println("M: ", h.M)
    fmt.Println("K: ", h.K)
    fmt.Println("Exponents:")
    for i := 0; i < h.K; i++ {
        fmt.Println("\t", i, h.Exponents[i])
    }
}

//-----------------------------------------------------------------------------
