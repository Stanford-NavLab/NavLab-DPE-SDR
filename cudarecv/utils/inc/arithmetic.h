#ifndef __INC_UTILS_ARITHMETIC_H_
#define __INC_UTILS_ARITHMETIC_H_

#include <cmath>
#include <vector>
#include <complex>
#include <type_traits>
#include <cassert>
#include <algorithm>

namespace dsp { namespace utils {

template<class T, class U, class Enable = void> struct MoreGeneralType { typedef int value; };
template<class T> struct MoreGeneralType <T, double, typename std::enable_if<std::is_integral<T>::value, T>::type> { typedef double value; };
template<class T> struct MoreGeneralType <double, T, typename std::enable_if<std::is_integral<T>::value, T>::type> { typedef double value; };
template<> struct MoreGeneralType <float, float> { typedef float value; };
template<> struct MoreGeneralType <double, double> { typedef double value; };
template<> struct MoreGeneralType <int, int> { typedef int value; };
template<> struct MoreGeneralType <float, double> { typedef double value; };
template<> struct MoreGeneralType <double, float> { typedef double value; };
template<class K> struct MoreGeneralType<std::complex<K>, std::complex<K>> { typedef std::complex<K> value; };
template<class R, class B> struct MoreGeneralType <std::complex<B>, R> { typedef std::complex<typename MoreGeneralType<B, R>::value> value; };
template<class T, class B> struct MoreGeneralType <T, std::complex<B>> { typedef std::complex<typename MoreGeneralType<B, T>::value> value; };

template<class T> struct is_complex : std::integral_constant<bool, false> {};
template<class R> struct is_complex <std::complex<R>> : std::integral_constant<bool, true> {};

template<class T>
struct is_number : std::integral_constant<bool,
      std::is_integral<T>::value ||
      std::is_floating_point<T>::value ||
      is_complex<T>::value> {};

template<bool A, bool B> struct both : std::integral_constant<bool, false> {};
template<> struct both<true, true> : std::integral_constant<bool, true> {};

template<bool A> struct negate {};
template<> struct negate<true> : std::integral_constant<bool, false> {};
template<> struct negate<false> : std::integral_constant<bool, true> {};

// scalar times vector

template<class T, class Q>
typename std::enable_if<both<is_number<Q>::value, std::is_same<typename MoreGeneralType<T, Q>::value, T>::value>::value, std::vector<T>>::type
operator* (const Q c, std::vector<T> A)
{
    std::transform(A.begin(), A.end(), A.begin(), std::bind1st(std::multiplies<T>(), c)) ;
    return A;
}

template<class T, class Q>
typename std::enable_if<both<is_number<Q>::value, std::is_same<typename MoreGeneralType<T, Q>::value, Q>::value>::value, std::vector<Q>>::type
operator* (const Q c, std::vector<T> A)
{
    // type-extend the vector first. e.g. Q = complex and T = double
    std::vector<Q> converted(A.begin(), A.end());
    std::transform(converted.begin(), converted.end(), converted.begin(), std::bind1st(std::multiplies<Q>(), c)) ;
    return converted;
}

template<class T, class Q>
typename std::enable_if<is_number<Q>::value, std::vector<typename MoreGeneralType<T, Q>::value>>::type
operator* (std::vector<T> A, const Q c)
{
    return c * A;
}

// vector-vector element-wise product

template<class T, class U>
typename std::enable_if<std::is_same<typename MoreGeneralType<T, U>::value, T>::value && negate<std::is_same<T, U>::value>::value, std::vector<T>>::type
operator* (std::vector<T> A, std::vector<U> B)
{
    assert(A.size() == B.size());
    std::transform(A.begin(), A.end(), B.begin(), A.begin(), std::multiplies<T>()) ;
    return A;
}

template<class T, class U>
typename std::enable_if<std::is_same<typename MoreGeneralType<T, U>::value, U>::value && negate<std::is_same<T, U>::value>::value, std::vector<U>>::type
operator* (std::vector<T> A, std::vector<U> B)
{
    assert(A.size() == B.size());
    std::transform(A.begin(), A.end(), B.begin(), B.begin(), std::multiplies<U>()) ;
    return B;
}

template<class T, class U>
typename std::enable_if<std::is_same<T, U>::value, std::vector<U>>::type
operator* (std::vector<T> A, std::vector<U> B)
{
    assert(A.size() == B.size());
    std::transform(A.begin(), A.end(), B.begin(), A.begin(), std::multiplies<T>()) ;
    return A;
}

// vector-vector element-wise sum

template<class T, class U>
typename std::enable_if<std::is_same<typename MoreGeneralType<T, U>::value, T>::value && negate<std::is_same<T, U>::value>::value, std::vector<T>>::type
operator+ (std::vector<T> A, std::vector<U> B)
{
    assert(A.size() == B.size());
    std::transform(A.begin(), A.end(), B.begin(), A.begin(), std::plus<T>()) ;
    return A;
}

template<class T, class U>
typename std::enable_if<std::is_same<typename MoreGeneralType<T, U>::value, U>::value && negate<std::is_same<T, U>::value>::value, std::vector<U>>::type
operator+ (std::vector<T> A, std::vector<U> B)
{
    assert(A.size() == B.size());
    std::transform(A.begin(), A.end(), B.begin(), B.begin(), std::plus<U>()) ;
    return B;
}

template<class T, class U>
typename std::enable_if<std::is_same<T, U>::value, std::vector<U>>::type
operator+ (std::vector<T> A, std::vector<U> B)
{
    assert(A.size() == B.size());
    std::transform(A.begin(), A.end(), B.begin(), A.begin(), std::plus<T>()) ;
    return A;
}

// vector conversions

/** \brief Casts an input vector to a desired type.
 *  \param input Incoming vector.
 *  \return Output vector whose elements are of the stipulated type.
 */
template<typename T, typename R>
std::vector<T> vector_cast(std::vector<R> input) {
    return std::vector<T>(input.begin(), input.end());
}

/** \brief Converts a real vector of size N to a complex vector of size N/2.
 */
template<typename T>
std::vector<std::complex<T>> reshapeToComplexVector(std::vector<T> input) {
    std::vector<std::complex<T>> comp(input.size() / 2);
    for (int i = 0; i < input.size() / 2; i++) {
        comp[i] = std::complex<T>(input[i * 2], input[i * 2 + 1]);
    }
    return comp;
}

/** \brief Converts a complex vector of size N to a real vector of size 2N.
 */
template <typename T>
std::vector<T> reshapeToRealVector(std::vector<std::complex<T>> input) {
    std::vector<T> flattened(input.size() * 2);
    for (int i = 0; i < input.size(); i++) {
        flattened[i * 2] = input[i].real();
        flattened[i * 2 + 1] = input[i].imag();
    }
    return flattened;
}

} }     // namespace utils; namespace dsp

namespace std {

template <bool B, typename T = void>
struct disable_if {
    typedef T type;
};

template <typename T>
struct disable_if<true, T> {
};

// extending functions to operate on vectors

template<typename T>
vector<T> exp(vector<T> A)
{
    transform(A.begin(), A.end(), A.begin(), [](T x) { return exp(x); });
    return A;
}

template<typename T>
vector<T> conj(vector<T> A)
{
    transform(A.begin(), A.end(), A.begin(), [](T x) { return conj(x); });
    return A;
}

template<class T>
typename std::enable_if<std::is_same<decltype(abs(T())), T>::value, vector<T>>::type
abs(vector<T> A)
{
    transform(A.begin(), A.end(), A.begin(), [](T x) { return abs(x); });
    return A;
}

template<class T>
typename std::disable_if<std::is_same<decltype(abs(T())), T>::value, vector<decltype(abs(T()))>>::type
abs(vector<T> A)
{
    typedef vector<decltype(abs(T()))> ReturnType;
    ReturnType B(A.size());
    transform(A.begin(), A.end(), B.begin(), [](T x) { return abs(x); });
    return B;
}

}

#endif
