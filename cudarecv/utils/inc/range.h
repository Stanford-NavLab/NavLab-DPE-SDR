#ifndef __INC_UTILS_RANGE_H_
#define __INC_UTILS_RANGE_H_

#include <vector>

namespace dsp { namespace utils {

template <long int T_begin, long int step, long int T_end>
class range_class {
public:
    class iterator {
        friend class range_class;

    public:
        long int operator*() const { return i_; }
        const iterator &operator++() {
            i_ += step;
            return *this;
        }
        iterator operator++(int) {
            iterator copy(*this);
            i_ += step;
            return copy;
        }

        bool operator ==(const iterator &other) const
        {
            if (other.isEnd_)
            {
                if (other.i_ > i_) return false; else return true;
            }
            if (isEnd_)
            {
                if (i_ < other.i_) return true; else return false;
            }
            return i_ == other.i_;
        }
        bool operator !=(const iterator &other) const {
            if (other.isEnd_)
            {
                if (other.i_ > i_) return true; else return false;
            }
            if (isEnd_)
            {
                if (i_ < other.i_) return false; else return true;
            }
            return i_ != other.i_;
        }

    protected:
        iterator(long int start, bool isEnd = false) : i_(start), isEnd_(isEnd) {}

    private:
        long i_;
        bool isEnd_;
    };

    iterator begin() const { return iterator(T_begin); }
    iterator end() const { return iterator(T_end, true); }
};

template <long int T_begin, long int T_step, long int T_end>
const range_class<T_begin, T_step, T_end> range = range_class<T_begin, T_step, T_end>();

template<typename T>
class ranged_class {
public:
    class iterator {
        friend class range;
    public:
        T operator *() const { return i_; }
        const iterator &operator ++() { i_ += step_; return *this; }
        iterator operator ++(int) { iterator copy(*this); i_ += step_; return copy; }

        bool operator ==(const iterator &other) const
        {
            if (other.isEnd_)
            {
                if (other.i_ > i_) return false; else return true;
            }
            if (isEnd_)
            {
                if (i_ < other.i_) return true; else return false;
            }
            return i_ == other.i_;
        }
        bool operator !=(const iterator &other) const {
            if (other.isEnd_)
            {
                if (other.i_ > i_) return true; else return false;
            }
            if (isEnd_)
            {
                if (i_ < other.i_) return false; else return true;
            }
            return i_ != other.i_;
        }

        iterator(T start, T step, bool isEnd = false) : i_ (start), step_ (step), isEnd_(isEnd) { }

        typedef T value_type;
        typedef T& reference;
        typedef T difference_type;
        typedef T* pointer;
        typedef std::input_iterator_tag iterator_category;

    private:
        T i_;
        T step_;
        bool isEnd_;
    };

    iterator begin() const { return begin_; }
    iterator end() const { return end_; }

public:
    ranged_class(T begin, T step, T end) : begin_(begin, step), end_(end, step, true) {}

    std::vector<T> toVector()
    {
        return std::vector<T>(begin(), end());
    }

private:
    iterator begin_;
    iterator end_;
};

template<typename T>
ranged_class<T> ranged(T begin, T step, T end)
{
    return ranged_class<T>(begin, step, end);
}

} } // namespace utils; namespace dsp

#endif
