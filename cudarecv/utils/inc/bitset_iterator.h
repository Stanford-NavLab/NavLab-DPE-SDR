#ifndef __INC_UTILS_BITSET_ITERATOR_H_
#define __INC_UTILS_BITSET_ITERATOR_H_

#include <bitset>
#include <vector>

template <bool flag, class IfTrue, class IfFalse>
struct IF;

template <class IfTrue, class IfFalse>
struct IF<true, IfTrue, IfFalse> {
    typedef IfTrue val;
};

template <class IfTrue, class IfFalse>
struct IF<false, IfTrue, IfFalse> {
    typedef IfFalse val;
};

template <std::size_t N, bool is_const>
class bitset_iterator {
    private:
        typedef std::bitset<N> Bitset;
        typedef typename IF<is_const, const Bitset, Bitset>::val
            qBitset;

        qBitset* B;
        std::size_t n;

    public:
        typedef std::random_access_iterator_tag iterator_category;
        typedef bool value_type;
        typedef std::ptrdiff_t difference_type;
        typedef typename IF<is_const, const bool, bool>::val* pointer;
        typedef typename IF<is_const,
                bool,
                typename Bitset::reference>::val reference;

        bitset_iterator() : B(), n() { }
        bitset_iterator(qBitset& b, std::size_t sz)
            : B(&b), n(sz) { }
        bitset_iterator(const bitset_iterator<N, false>& x)
            : B(x.B), n(x.n) { }

        bitset_iterator& operator=(const bitset_iterator& x) {
            B = x.B;
            n = x.n;
            return *this;
        }

    public:
        reference operator*() const { return (*B)[n]; }
        reference operator[](std::ptrdiff_t x) const {
            return (*B)[n + x];
        }

        bitset_iterator& operator++() { ++n; return *this; }
        bitset_iterator operator++(int) {
            ++n;
            return bitset_iterator(*B, n-1);
        }
        bitset_iterator& operator--() { --n; return *this; }
        bitset_iterator operator--(int) {
            --n;
            return bitset_iterator(*B, n+1);
        }

        bitset_iterator operator+(std::ptrdiff_t x) const {
            return bitset_iterator(*B, n + x);
        }
        bitset_iterator& operator+=(std::ptrdiff_t x) {
            n += x;
            return *this;
        }
        bitset_iterator operator-(std::ptrdiff_t x) const {
            return bitset_iterator(*B, n - x);
        }
        bitset_iterator& operator-=(std::ptrdiff_t x) {
            n -= x;
            return *this;
        }

    public:
        friend bool operator==(bitset_iterator x,
                bitset_iterator y) {
            return x.B == y.B && x.n == y.n;
        }
        friend bool operator!=(bitset_iterator x,
                bitset_iterator y) {
            return !(x == y);
        }

        friend bool operator<(bitset_iterator x,
                bitset_iterator y) {
            return x.n < y.n;
        }
        friend bool operator>(bitset_iterator x,
                bitset_iterator y) {
            return y < x;
        }
        friend bool operator<=(bitset_iterator x,
                bitset_iterator y) {
            return !(y < x);
        }
        friend bool operator>=(bitset_iterator x,
                bitset_iterator y) {
            return !(x < y);
        }

        friend std::ptrdiff_t operator-(bitset_iterator x,
                bitset_iterator y) {
            return x.n - y.n;
        }
        friend bitset_iterator operator+(std::ptrdiff_t n1,
                bitset_iterator x) {
            return bitset_iterator(*x.B, x.n + n1);
        }
};


template <std::size_t N>
bitset_iterator<N, true>
begin(const std::bitset<N>& b) {
    return bitset_iterator<N, true>(b, 0);
}

template <std::size_t N>
bitset_iterator<N, true>
end(const std::bitset<N>& b) {
    return bitset_iterator<N, true>(b, N);
}

template <std::size_t N>
bitset_iterator<N, false>
begin(std::bitset<N>& b) {
    return bitset_iterator<N, false>(b, 0);
}

template <std::size_t N>
bitset_iterator<N, false>
end(std::bitset<N>& b) {
    return bitset_iterator<N, false>(b, N);
}

#endif
