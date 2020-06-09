
#ifndef INC__ERRORHANDLER_H_
#define INC__ERRORHANDLER_H_

#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cuComplex.h>
#include <cufft.h>
#include <iostream>


/**
 * cufft error handling
 */

#define cuCheckHelper(stmt, failStmt)                               \
    do {                                                            \
        cudaError_t err = stmt;                                     \
        if (err != cudaSuccess){                                    \
            std::cerr << "Error: " << #stmt << std::endl            \
                      << "in file: " << __FILE__ << " and line: "   \
                      << __LINE__ << std::endl                      \
                      << cudaGetErrorString(err) << std::endl ;     \
            failStmt;                                               \
        }                                                           \
    } while (0)

#define cuCheck(stmt)       cuCheckHelper(stmt, return -1)
#define cuCheckV(stmt)      cuCheckHelper(stmt, )

#define cuCheckModHelper(stmt, failStmt)                            \
    do {                                                            \
        cudaError_t err = stmt;                                     \
        if (err != cudaSuccess){                                    \
            std::cerr << "[" << ModuleName << "] Error: " << #stmt  \
                      << "\nin file: " << __FILE__ << " and line: " \
                      << __LINE__ << "\n"                           \
                      << cudaGetErrorString(err) << std::endl;      \
            failStmt;                                               \
        }                                                           \
    } while (0)

#define cuCheckMSt(stmt)    cuCheckModHelper(stmt, Stop(); return -1)
#define cuCheckM(stmt)      cuCheckModHelper(stmt, return -1)
#define cuCheckMSp(stmt)    cuCheckModHelper(stmt, ret = -1)

static const char *_cudaGetErrorEnum(cufftResult error) {
    switch (error) {
        case CUFFT_SUCCESS:
            return "CUFFT_SUCCESS";
        case CUFFT_INVALID_PLAN:
            return "CUFFT_INVALID_PLAN";
        case CUFFT_ALLOC_FAILED:
            return "CUFFT_ALLOC_FAILED";
        case CUFFT_INVALID_TYPE:
            return "CUFFT_INVALID_TYPE";
        case CUFFT_INVALID_VALUE:
            return "CUFFT_INVALID_VALUE";
        case CUFFT_INTERNAL_ERROR:
            return "CUFFT_INTERNAL_ERROR";
        case CUFFT_EXEC_FAILED:
            return "CUFFT_EXEC_FAILED";
        case CUFFT_SETUP_FAILED:
            return "CUFFT_SETUP_FAILED";
        case CUFFT_INVALID_SIZE:
            return "CUFFT_INVALID_SIZE";
        case CUFFT_UNALIGNED_DATA:
            return "CUFFT_UNALIGNED_DATA";
        case CUFFT_INCOMPLETE_PARAMETER_LIST:
            return "CUFFT_INCOMPLETE_PARAMETER_LIST";
        case CUFFT_INVALID_DEVICE:
            return "CUFFT_INVALID_DEVICE";
        case CUFFT_PARSE_ERROR:
            return "CUFFT_PARSE_ERROR";
        case CUFFT_NO_WORKSPACE:
            return "CUFFT_NO_WORKSPACE";
        case CUFFT_NOT_IMPLEMENTED:
            return "CUFFT_NOT_IMPLEMENTED";
        case CUFFT_LICENSE_ERROR:
            return "CUFFT_LICENSE_ERROR";
        case CUFFT_NOT_SUPPORTED:
            return "CUFFT_NOT_SUPPORTED";
        default:
            return "<unknown>";
    }
}

#define cufftCheckModHelper(stmt, failStmt)                         \
    do {                                                            \
        cufftResult err = stmt;                                     \
        if (err != CUFFT_SUCCESS){                                  \
            std::cerr << "[" << ModuleName << "] Error: " << #stmt  \
                      << "\nin file: " << __FILE__ << " and line: " \
                      << __LINE__ << "\n"                           \
                      << _cudaGetErrorEnum(err) << std::endl;       \
            failStmt;                                               \
        }                                                           \
    } while (0)

#define cufftCheckMSt(stmt)    cufftCheckModHelper(stmt, Stop(); return -1)
#define cufftCheckM(stmt)      cufftCheckModHelper(stmt, return -1)
#define cufftCheckMSp(stmt)    cufftCheckModHelper(stmt, ret = -1)

#define ptrCheckHelper(ptr, failStmt)                                   \
    do {                                                                \
        if (ptr == NULL){                                               \
            std::cerr << "Error: Malloc Failed: " << #ptr << std::endl  \
                      << "in file: " << __FILE__ << " and line: "       \
                      << __LINE__ << std::endl;                         \
            failStmt;                                                   \
        }                                                               \
    } while (0)

#define ptrCheck(ptr)       ptrCheckHelper(ptr, return -1)

#define errCheckHelper(stmt, failStmt)                               \
    do {                                                            \
        int err = stmt;                                     \
        if (err < 0){                                    \
            std::cerr << "Error: " << #stmt << std::endl            \
                      << "in file: " << __FILE__ << " and line: "   \
                      << __LINE__ << std::endl;                     \
            failStmt;                                               \
        }                                                           \
    } while (0)

#define errCheck(stmt)      errCheckHelper(stmt, return -1)

#define errCheckModHelper(stmt, failStmt)                           \
    do {                                                            \
        int err = stmt;                                             \
        if (err < 0){                                               \
            std::cerr << "[" << ModuleName << "] Error: "           \
                      << #stmt << std::endl                         \
                      << "in file: " << __FILE__ << " and line: "   \
                      << __LINE__ << std::endl;                     \
            failStmt;                                               \
        }                                                           \
    } while (0)

#define errCheckMod(stmt)      errCheckModHelper(stmt, return -1)
#define errCheckModSt(stmt)    errCheckModHelper(stmt, Stop(); return -1)




/**
 * cuBLAS error handling
 */


static const char *_cudaBLASGetErrorEnum(cublasStatus_t error) {
    switch (error) {
        case CUBLAS_STATUS_SUCCESS:
            return "CUBLAS_STATUS_SUCCESS";
        case CUBLAS_STATUS_NOT_INITIALIZED:
            return "CUBLAS_STATUS_NOT_INITIALIZED";
        case CUBLAS_STATUS_ALLOC_FAILED:
            return "CUBLAS_STATUS_ALLOC_FAILED";
        case CUBLAS_STATUS_INVALID_VALUE:
            return "CUBLAS_STATUS_INVALID_VALUE";
        case CUBLAS_STATUS_ARCH_MISMATCH:
            return "CUBLAS_STATUS_ARCH_MISMATCH";
        case CUBLAS_STATUS_MAPPING_ERROR:
            return "CUBLAS_STATUS_MAPPING_ERROR";
        case CUBLAS_STATUS_EXECUTION_FAILED:
            return "CUBLAS_STATUS_EXECUTION_FAILED";
        case CUBLAS_STATUS_INTERNAL_ERROR:
            return "CUBLAS_STATUS_INTERNAL_ERROR";
        case CUBLAS_STATUS_NOT_SUPPORTED:
            return "CUBLAS_STATUS_NOT_SUPPORTED";
        case CUBLAS_STATUS_LICENSE_ERROR:
            return "CUBLAS_STATUS_LICENSE_ERROR";
        default:
            return "<unknown>";
    }
}

#define cuBLASCheckModHelper(stmt, failStmt)                            \
    do {                                                            \
        cublasStatus_t err = stmt;                                     \
        if (err != CUBLAS_STATUS_SUCCESS){                                    \
            std::cerr << "[" << ModuleName << "] Error: " << #stmt  \
                      << "\nin file: " << __FILE__ << " and line: " \
                      << __LINE__ << "\n"                           \
                      << _cudaBLASGetErrorEnum(err) << std::endl;      \
            failStmt;                                               \
        }                                                           \
    } while (0)

#define cuBLASCheckMSt(stmt)    cuBLASCheckModHelper(stmt, Stop(); return -1)
#define cuBLASCheckM(stmt)      cuBLASCheckModHelper(stmt, return -1)
#define cuBLASCheckMSp(stmt)    cuBLASCheckModHelper(stmt, ret = -1)



// Macro definition
#define cuAssert( X ) if ( !(X) ) { printf( "Thread %d:%d failed assert at %s:%d!", blockIdx.x, threadIdx.x, __FILE__, __LINE__ ); return; }



/**
 *  Misc helpers
 */
#define funcLogHelper(msg, failStmt)									\
do {																	\
		std::cout << "In " << __func__ << ", " << msg << std::endl;		\
		failStmt;														\
	} while (0)

#define funcLogHelperM(msg)			funcLogHelper(msg, return -1)
#define funcLogHelperSp(msg)		funcLogHelper(msg, ret = -1)



#endif
