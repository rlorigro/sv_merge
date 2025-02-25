#pragma once

#include "MurmurHash3.hpp"
#include "MurmurHash2.hpp"

#include <type_traits>
#include <iostream>
#include <vector>
#include <bitset>
#include <string>
#include <deque>
#include <cmath>
#include <array>
#include <bit>

using std::runtime_error;
using std::is_integral;
using std::to_string;
using std::vector;
using std::bitset;
using std::string;
using std::deque;
using std::array;
using std::pow;
using std::cerr;
using std::rotr;
using std::rotl;


namespace sv_merge {


template<class T> class BinarySequence {
public:
    /// Attributes ///
    vector <T> sequence;
    int32_t length;

    static const array<char,4> index_to_base;
    static const array<uint16_t,128> base_to_index;

    /// Methods ///
    BinarySequence();
    BinarySequence(const BinarySequence& s);
    template <class T2> explicit BinarySequence(const T2& s);
    [[nodiscard]] BinarySequence<T> substr(size_t start, size_t n) const;
    void encode(const string& s);
    void push_back(char c);
    void shift(char c);
    void reverse_complement();
    void clear();

    // For short sequences it may be important to reserve so that default exponential allocator does not overallocate
    void reserve(size_t bp_length);
    [[nodiscard]] bool empty() const;
    void to_string(string& s) const;
    [[nodiscard]] size_t get_byte_length() const;
    [[nodiscard]] size_t size() const;

    void print_as_bits() const;

    [[nodiscard]] T get_random_base() const;
};


template <class T> void BinarySequence<T>::clear(){
    sequence.clear();
    length = 0;
}


template <class T> void BinarySequence<T>::reserve(size_t bp_length){
    // 4 bases per 8bits
    auto n_bases_per_word = sizeof(T)*4;

    size_t words_needed = int(ceil(double(bp_length) / double(n_bases_per_word)));

    sequence.reserve(words_needed);
}


template <class T> bool BinarySequence<T>::empty() const{
    return length == 0;
}


template <class T> void BinarySequence<T>::print_as_bits() const {
    for (auto& item: sequence){
        cerr << bitset<sizeof(T)*8>(item) << ' ';
    }
    cerr << '\n';
}


template <class T> bool operator==(const BinarySequence<T>& a, const BinarySequence<T>& b)
{
    return a.sequence == b.sequence;
}


template <class T> const array<char,4> BinarySequence<T>::index_to_base = {'A','C','G','T'};
template <class T> const array<uint16_t,128> BinarySequence<T>::base_to_index = {
        4,4,4,4,4,4,4,4,4,4,      // 0
        4,4,4,4,4,4,4,4,4,4,      // 10
        4,4,4,4,4,4,4,4,4,4,      // 20
        4,4,4,4,4,4,4,4,4,4,      // 30
        4,4,4,4,4,4,4,4,4,4,      // 40
        4,4,4,4,4,4,4,4,4,4,      // 50
        4,4,4,4,4,0,4,1,4,4,      // 60  A = 65, C = 67
        4,2,4,4,4,4,4,4,4,4,      // 70  G = 71
        4,4,4,4,3,4,4,4,4,4,      // 80  T = 84
        4,4,4,4,4,4,4,4,4,4,      // 90
        4,4,4,4,4,4,4,4,4,4,      // 100
        4,4,4,4,4,4,4,4,4,4,      // 110
        4,4,4,4,4,4,4,4           // 120
};


template<class T> BinarySequence<T>::BinarySequence():
        length(0)
{}


template<class T> T BinarySequence<T>::get_random_base() const{
    return rand() % 4;
}


template <class T> BinarySequence<T>::BinarySequence(const BinarySequence<T>& s):
        sequence(s.sequence),
        length(s.length)
{}


template<class T> template <class T2> BinarySequence<T>::BinarySequence(const T2& s):
        length(0)
{
    static_assert(is_integral<T>::value, "ERROR: provided type for BinarySequence is not integer");

    for (auto& c: s){
        push_back(c);
    }
}


template <class T> void BinarySequence<T>::push_back(char c){
    T bits = base_to_index.at(c);

    if (bits == 4){
        // throw runtime_error("ERROR: non ACGT character encountered in sequence: " + string(1,c) + " (ord=" + std::to_string(int(c)) + ")");
        bits = get_random_base();
    }

    uint8_t shift_size = (2*length) % (sizeof(T)*8);

    // If we have reached the beginning of a new word, append the word vector with 0
    if (shift_size == 0){
        sequence.emplace_back(0);
    }
    // Otherwise shift the last word over to prepare to insert another base
    else {
        bits <<= shift_size;
    }

    sequence.back() |= bits;
    length++;
}


template <class T> void BinarySequence<T>::encode(const string& s){
    auto n_bases_per_word = sizeof(T)*4;
    size_t words_needed = int(ceil(double(s.length()) / double(n_bases_per_word)));

    sequence.reserve(words_needed);

    for (auto& c: s) {
        T bits = base_to_index.at(c);

        if (bits == 4){
            // throw runtime_error("ERROR: non ACGT character encountered in sequence: " + string(1,c) + " (ord=" + std::to_string(int(c)) + ")");
            bits = get_random_base();
        }

        uint8_t shift_size = (2*length) % (sizeof(T)*8);

        // If we have reached the beginning of a new word, append the word vector with 0
        if (shift_size == 0){
            sequence.emplace_back(0);
        }
        // Otherwise shift the last word over to prepare to insert another base
        else {
            bits <<= shift_size;
        }

        sequence.back() |= bits;
        length++;
    }
}


/// A function to push a new base but not alter the length of the sequence (as a fixed length queue would)
/// \tparam T
/// \param c
template <class T> void BinarySequence<T>::shift(char c){
    T bits = base_to_index.at(c);

    if (bits == 4){
        // throw runtime_error("ERROR: non ACGT character encountered in sequence: " + string(1,c) + " (ord=" + std::to_string(int(c)) + ")");
        bits = get_random_base();
    }

    for (size_t i=0; i<sequence.size(); i++){
        if (i > 0) {
            T mask = 3;
//            cerr << bitset<sizeof(T)*8>(mask) << '\n';
//            cerr << '\n';

            T leftover = sequence[i] & mask;
//            cerr << bitset<sizeof(T)*8>(sequence[i]) << '\n';
//            cerr << bitset<sizeof(T)*8>(leftover) << '\n';

            leftover <<= sizeof(T)*8 - 2;
//            cerr << bitset<sizeof(T)*8>(leftover) << '\n';
//            cerr << '\n';

//            cerr << bitset<sizeof(T)*8>(sequence[i-1]) << '\n';
//            cerr << bitset<sizeof(T)*8>(T(pow(2,sizeof(T)*8 - 2) - 1)) << '\n';
            sequence[i-1] &= T(pow(2,sizeof(T)*8 - 2) - 1);
//            cerr << bitset<sizeof(T)*8>(sequence[i-1]) << '\n';

            sequence[i-1] |= leftover;
//            cerr << bitset<sizeof(T)*8>(sequence[i-1]) << '\n';
//
//            cerr << '\n';
//            cerr << '\n';

        }
        sequence[i] >>= 2;
    }

    uint8_t shift_size = (2*(length-1)) % (sizeof(T)*8);

    bits <<= shift_size;

    sequence.back() |= bits;
}


template<class T> void BinarySequence<T>::to_string(string& s) const{
    if (sequence.empty()){
        return;
    }

    s.clear();
    s.reserve(length);

    T mask = 3;

    size_t n_bits = sizeof(T)*8;
    size_t n_remainder_bases = length - (n_bits*(sequence.size() - 1))/2;

    for (size_t i=0; i<sequence.size(); i++){
        T word = sequence[i];
        size_t l = n_bits/2;

        if (i == sequence.size() - 1){
            l = n_remainder_bases;
        }

        // Iterate/consume the word and produce bases (chars)
        for (size_t j=0; j<l; j++){
//            cerr << bitset<sizeof(T)*8>(word) << '\n';

            auto index = word & mask;
            s += index_to_base[index];
            word >>= 2;
        }
    }
}


template<class T> void BinarySequence<T>::reverse_complement(){
    if (sequence.empty()){
        return;
    }

    // We iterate the existing bit vector backwards, and a new one forwards.
    // Length doesn't change. Remainder bases/bits don't change.

    // Initialize a temp sequence of the same length as this sequence, all zeros
    vector<T> temp_sequence(sequence.size(), 0);

    // 0b00000011 or however many bits are used in T int
    T mask = 3;

    int64_t n_bits = sizeof(T)*8;
    int64_t n_remainder_bases = length - (n_bits*(sequence.size() - 1))/2;

    // As we iterate backwards, the integer bounds occur in different places for the forward and reverse vectors
    // so they are tracked separately
    int64_t i_temp = 0;
    int64_t j_temp = 0;

    for (int64_t i=sequence.size() - 1; i>=0; i--){
        T word = sequence[i];
        int64_t l = n_bits/2;

        if (i == sequence.size() - 1){
            l = n_remainder_bases;
        }

        // Iterate the input word and update the forward oriented vector with the complement
        for (int64_t j=0; j<(n_bits/2); j++){
            i_temp = j_temp / (n_bits/2);

            word = rotl(word, 2);

            if ((n_bits/2) - j <= l) {
                temp_sequence[i_temp] |= 3 - word & mask;
                temp_sequence[i_temp] = rotr(temp_sequence[i_temp], 2);

                j_temp++;
            }
        }

        // Only the last word needs to be shifted because it has dangling 00s
        if (i_temp == temp_sequence.size() - 1) {
            temp_sequence[i_temp] >>= (n_bits/2 - n_remainder_bases)*2;
        }
    }

    sequence = std::move(temp_sequence);
}


template<class T> BinarySequence<T> BinarySequence<T>::substr(size_t start, size_t n) const{
    BinarySequence<T> s;
    size_t n_bits = sizeof(T)*8;

    // Reserve exactly how much space is needed for the sequence
    auto n_bases_per_word = sizeof(T)*4;
    size_t words_needed = int(ceil(double(n) / double(n_bases_per_word)));

    s.sequence.reserve(words_needed);

    if (start + n > length) {
        throw runtime_error("ERROR: cannot extract substr greater than size of BinarySequence: start: " + std::to_string(start) + " n: " + std::to_string(n) + " length: " + std::to_string(length));
    }

    if (empty()) {
        return s;
    }

    T mask = 3;

    size_t i_start = (start*2) / n_bits;
    size_t i_stop = ((start + n)*2) / n_bits;

    size_t stop = start + n;

    // cerr << "start: " << start << '\n';
    // cerr << "n: " << n << '\n';
    // cerr << "n_bits: " << n_bits << '\n';
    // cerr << "i_start: " << i_start << '\n';
    // cerr << "i_stop: " << i_stop << '\n';
    // cerr << "start: " << start << '\n';
    // cerr << "stop: " << stop << '\n';

    for (size_t i=i_start; i<i_stop + 1; i++){
        T word = sequence[i];
        size_t a = 0;
        size_t b = n_bits/2;

        // cerr << "i: " << i << '\n';
        // cerr << "a: " << a << '\n';
        // cerr << "b: " << b << '\n';

        if (i == i_stop){
            b = stop % (n_bits/2);
        }

        if (i == i_start){
            a = start % (n_bits/2);
        }

        // Iterate/consume the word and produce bases (chars)
        for (size_t j=0; j<b; j++){
            // cerr << bitset<sizeof(T)*8>(word) << '\n';

            if (j >= a) {
                auto index = word & mask;
                s.push_back(index_to_base[index]);
            }

            word >>= 2;
        }
    }

    return s;
}


template <class T> size_t BinarySequence<T>::get_byte_length() const{
    size_t bit_length = size_t(length)*2;
    size_t byte_length = bit_length/8 + (bit_length % 8 != 0);

    return byte_length;
}


template <class T> size_t BinarySequence<T>::size() const{
    return size_t(length);
}




}


namespace std {
template<>
class hash<sv_merge::BinarySequence<uint64_t> > {
public:
    size_t operator()(const sv_merge::BinarySequence<uint64_t>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.get_byte_length()), 14741);
    }
};


template<>
class hash<sv_merge::BinarySequence<int64_t> > {
public:
    size_t operator()(const sv_merge::BinarySequence<int64_t>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.get_byte_length()), 14741);
    }
};


template<>
class hash<sv_merge::BinarySequence<uint32_t> > {
public:
    size_t operator()(const sv_merge::BinarySequence<uint32_t>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.get_byte_length()), 14741);
    }
};


template<>
class hash<sv_merge::BinarySequence<int32_t> > {
public:
    size_t operator()(const sv_merge::BinarySequence<int32_t>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.get_byte_length()), 14741);
    }
};


template<>
class hash<sv_merge::BinarySequence<uint16_t> > {
public:
    size_t operator()(const sv_merge::BinarySequence<uint16_t>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.get_byte_length()), 14741);
    }
};


template<>
class hash<sv_merge::BinarySequence<int16_t> > {
public:
    size_t operator()(const sv_merge::BinarySequence<int16_t>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.get_byte_length()), 14741);
    }
};


template<>
class hash<sv_merge::BinarySequence<uint8_t> > {
public:
    size_t operator()(const sv_merge::BinarySequence<uint8_t>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.get_byte_length()), 14741);
    }
};


template<>
class hash<sv_merge::BinarySequence<int8_t> > {
public:
    size_t operator()(const sv_merge::BinarySequence<int8_t>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.get_byte_length()), 14741);
    }
};


template<>
class hash<sv_merge::BinarySequence<__uint128_t> > {
public:
    size_t operator()(const sv_merge::BinarySequence<__uint128_t>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.get_byte_length()), 14741);
    }
};


template<>
class hash<sv_merge::BinarySequence<__int128_t> > {
public:
    size_t operator()(const sv_merge::BinarySequence<__int128_t>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.get_byte_length()), 14741);
    }
};


}
