//-----------------------------------------------------------------------------
// Author: Vinhthuy Phan, 2018
//-----------------------------------------------------------------------------
package ppt_filter

import (
    "fmt"
    "math/big"
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
    A               *big.Int // randomly selected
    B               *big.Int // randomly selected
    P               *big.Int // large prime
    M               int64    // universe/hash table size
    K               int
    Base            *big.Int
    Term0           *big.Int
    Term0_rc        *big.Int
    PrevValue       *big.Int
    PrevValue_rc    *big.Int
    Exponents       []*big.Int
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
        A:              big.NewInt(a),
        B:              big.NewInt(b),
        P:              big.NewInt(p),
        M:              m,
        Term0:          big.NewInt(0),
        Term0_rc:       big.NewInt(0),
        PrevValue:      big.NewInt(0),
        PrevValue_rc:   big.NewInt(0),
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
        Term0:          big.NewInt(0),
        Term0_rc:       big.NewInt(0),
        PrevValue:      big.NewInt(0),
        PrevValue_rc:   big.NewInt(0),
        Exponents:      linear_hash.Exponents,
    }
}

//-----------------------------------------------------------------------------
func NewLinearHashSetK(m int64, k int) *LinearHash {
    rand.Seed(time.Now().UTC().UnixNano())

    // temporarily: 2^61 - 1
    var pint int64 = 2305843009213693951
    p := big.NewInt(pint)
    // generate random numbers in [2, p-1]
    a := big.NewInt((rand.Int63()%(pint-2) + 2))
    b := big.NewInt((rand.Int63()%(pint-2) + 2))
    // SetK
    Expo := make([]*big.Int, k)
    tbase := big.NewInt(rand.Int63n(65536-4) + 4)
    tmp := big.NewInt(1)
    for i := k - 1; i >= 0; i-- {
        Expo[i] = big.NewInt(1)
        Expo[i].Mul(Expo[i], tmp)
        Expo[i].Mod(Expo[i], p)
        tmp = tmp.Mul(tbase, tmp)
    }

    return &LinearHash{
        // A: big.NewInt(23),
        // B: big.NewInt(17),
        // P: big.NewInt(97),
        A:              a,
        B:              b,
        P:              p,
        M:              m,
        Term0:          big.NewInt(0),
        Term0_rc:       big.NewInt(0),
        PrevValue:      big.NewInt(0),
        PrevValue_rc:   big.NewInt(0),
        Base:           tbase,
        K:              k,
        Exponents:      Expo,
    }
}

//-----------------------------------------------------------------------------
func (h *LinearHash) SetK(k int) {
    h.K = k
    h.Exponents = make([]*big.Int, k)
    h.Base = big.NewInt(rand.Int63n(65536-4) + 4)
    // h.Base = big.NewInt(4)
    b := big.NewInt(1)
    for i := k - 1; i >= 0; i-- {
        h.Exponents[i] = big.NewInt(1)
        h.Exponents[i].Mul(h.Exponents[i], b)
        h.Exponents[i].Mod(h.Exponents[i], h.P)
        b = b.Mul(h.Base, b)
    }
}

//-----------------------------------------------------------------------------
func (h *LinearHash) ComputeKmer(read []byte, start int, k int) int64 {
    // fmt.Println("---Compute Kmer: ", string(kmer))
    // if len(kmer) != h.K {
    //     panic("Unmatched k-mer length")
    // }
    fmt.Println("ComputeKmer - read ", string(read))
    fmt.Println("ComputeKmer - kmer ", string(read[start:start+k]))
    var base *big.Int
    value := big.NewInt(0)
    for i := start; i <= (start + k - 1); i++ {
        // fmt.Println(i, string(read[i]))
        if read[i] == 'A' {
            base = big.NewInt(0)
        } else if read[i] == 'C' {
            base = big.NewInt(1)
        } else if read[i] == 'G' {
            base = big.NewInt(2)
        } else if read[i] == 'T' {
            base = big.NewInt(3)
        } else {
            // fmt.Println(string(kmer))
            panic("ComputeKmer: " + string(read) + " ------ Unknown character: " + string(read[i]))
        }
        cur_term := big.NewInt(0)
        cur_term.Mul(base, h.Exponents[i-start])
        cur_term.Mod(cur_term, h.P)
        value.Add(value, cur_term)
        value.Mod(value, h.P)
        if (i-start) == 0 {
            h.Term0 = cur_term
        }
    }
    h.PrevValue = value
    return value.Int64() % h.M
}

//-----------------------------------------------------------------------------
func (h *LinearHash) HashKmer(read []byte, start int, k int) int64 {
    // fmt.Println("HashKmer func: ", string(kmer))
    // if len(kmer) != h.K {
    //     panic("Unmatched k-mer length")
    // }    

    return h.HashInt64(h.ComputeKmer(read, start, k))
}

//-----------------------------------------------------------------------------
func (h *LinearHash) HashInt64(x int64) int64 {
    value := big.NewInt(0)
    value.Mul(h.A, big.NewInt(x))
    value.Add(value, h.B)
    value.Mod(value, h.P)
    return value.Int64() % h.M
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
