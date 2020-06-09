#ifndef __INC_UTILS_FUNCTIONAL_H_
#define __INC_UTILS_FUNCTIONAL_H_

#include <vector>
#include <numeric>
#include "function_traits.h"
#include "arithmetic.h"

namespace dsp { namespace utils {

struct identity {
    template<typename U>
    constexpr auto operator()(U&& v) const noexcept
        -> decltype(std::forward<U>(v))
    {
        return std::forward<U>(v);
    }
};

template<class T, class Func>
std::vector<typename function_traits<Func>::result_type>
map(const std::vector<T> input, Func func) {
    typedef typename function_traits<Func>::result_type ResultType;
    std::vector<ResultType> output(input.size());
    std::transform(input.begin(), input.end(), output.begin(), func);
    return output;
}

template<class T>
T sum(std::vector<T> inputs) {
    return std::accumulate(inputs.begin() + 1, inputs.end(), inputs[0], [](T a, T b) { return a + b; });
}

} }     // namespace utils; namespace dsp

#endif
