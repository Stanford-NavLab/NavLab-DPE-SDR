#ifndef __INC_UTILS_ITERATOR_REPEATER_H_
#define __INC_UTILS_ITERATOR_REPEATER_H_

namespace dsp { namespace utils {

namespace details {

template <typename T>
class IteratorRepeater {
public:
    typename std::iterator_traits<T>::value_type operator*() const { return *internalItr_; }
    const IteratorRepeater& operator++() {
        counter_++;
        if (counter_ % repeatCount_ == 0)
            internalItr_++;
        counter_ %= repeatCount_;
        return *this;
    }

    IteratorRepeater operator+(const int steps) {
        counter_ += steps;
        int numIncr = counter_ / repeatCount_;
        counter_ %= repeatCount_;
        internalItr_ += numIncr;
        return *this;
    }

    bool operator==(const IteratorRepeater &other) const {
        return internalItr_ == other.internalItr_ && counter_ == other.counter_;
    }
    bool operator!=(const IteratorRepeater &other) const {
        return internalItr_ != other.internalItr_ || counter_ != other.counter_;
    }

    bool operator==(const T &other) const { return internalItr_ == other && counter_ == 0; }
    bool operator!=(const T &other) const { return internalItr_ != other || counter_ != 0; }

    IteratorRepeater(T start, int repeatCount)
        : internalItr_(start), counter_(0), repeatCount_(repeatCount)
    {
    }

    IteratorRepeater(const IteratorRepeater& other)
    {
        internalItr_ = other.internalItr_;
        counter_ = other.counter_;
        repeatCount_ = other.repeatCount_;
    }

    IteratorRepeater& operator=(const IteratorRepeater& other)
    {
        internalItr_ = other.internalItr_;
        counter_ = other.counter_;
        repeatCount_ = other.repeatCount_;
        return *this;
    }

    /* Type definitions */
    typedef typename std::iterator_traits<T>::value_type value_type;

private:
    T internalItr_;
    int counter_;
    int repeatCount_;
};

}

template <typename T>
details::IteratorRepeater<T> iteratorRepeater(T start, int repeatCount) {
    return details::IteratorRepeater<T>(start, repeatCount);
}

} } // namespace utils; namespace dsp

#endif
