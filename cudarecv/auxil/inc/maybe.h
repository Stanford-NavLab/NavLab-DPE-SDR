#ifndef _INC_MAYBE_H_
#define _INC_MAYBE_H_

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <functional>
#include <vector>
#include <tuple>
#include "function_traits.h"

#ifdef _OPENMP
    #include <omp.h>
#endif

namespace auxil {

namespace maybe_details {

#define _unused(x) ((void)(x))
class NothingType {
    uint8_t not_used;
public:
    NothingType() { _unused(not_used); }
};

}   // namespace maybe_details

extern maybe_details::NothingType Nothing;

template <typename T>
union MaybeUnion
{
    uint8_t _nothing;
    T v;
    MaybeUnion() : _nothing(0) {}
    ~MaybeUnion() {}
};

template <typename T>
struct identity
{
  typedef T type;
};

template <typename T>
class Maybe {
private:
    MaybeUnion<T> val;
    bool exists;

public:
    Maybe() : exists(false) { }
    Maybe(const T& bla) : exists(true) { val.v = bla; }
    Maybe(const Maybe<T>& other) : exists(other.exists)
    {
        if (exists)
        {
            val.v = other.val.v;
        }
    }
    ~Maybe() { if (exists) (&val.v)->~T(); }
    Maybe(const maybe_details::NothingType nothing) : exists(false) { }

    Maybe& operator = (const Maybe& other)
    {
        if(exists && (!other.exists))
            (&val.v)->~T();
        exists = other.exists;
        if (exists)
        {
            val.v = other.val.v;
        }
        return *this;
    }

    const T& operator + () const {
        if (!exists)
        {
            // Trigger an exception
            throw std::logic_error("Unwrapping empty optional value");
            // We'll never reach here
            for (;;) { }
            return val.v;
        }
        return val.v;
    }

    const Maybe<T> operator >> (std::function<Maybe<T>()> other) const {
        if (exists)
        {
            return *this;
        }
        else
        {
            return other();
        }
    }

    template<class Cont>
    const typename function_traits<Cont>::result_type operator >>= (Cont other) const {
        if (exists)
        {
            return other(+*this);
        }
        else
        {
            return Nothing;
        }
    }

    bool operator ! () const { return !exists; }
    operator bool() const { return exists; }
};

template<class T> struct remove_maybe {};
template<class R> struct remove_maybe<Maybe<R>> { typedef R type; };

template<class T, class Func>
std::vector<typename remove_maybe<typename function_traits<Func>::result_type>::type>
parallelFlatMap(const std::vector<T> input, Func func) {
    typedef typename remove_maybe<typename function_traits<Func>::result_type>::type ResultType;
    std::vector<ResultType> output;

    #if defined(HAS_OPENMP)
        #pragma omp parallel for
    #endif
    for (int i = 0; i <= input.size(); i++) {
        auto maybeVal = func(input[i]);
        if (maybeVal) {
            #if defined(HAS_OPENMP)
                #pragma omp critical
            #endif
            {
                output.push_back(+maybeVal);
            }
        }
    }

    return output;
}

}   // namespace auxil

#endif
