
#ifndef INC__DSP_H_
#define INC__DSP_H_

#include <iostream>
//#include "flow.h"
//#include "scalarflow.h"
//#include "module.h"

/** \brief Namespace containing data structures for DSP Flow.
*/
namespace dsp {

    /** \todo Add Typedefs
     * Plan to support float and fixed point types.
     * Plan for different Q format fixed point types.
     */

    const unsigned short VECTORLENGTH_ANY = 0;
    const uint8_t ANY = 0;



    /** \brief The type of quantity this value represents.
     *  Note that this does not represent the meaning of the value. For
     *  example, an SNR value will be of type RATIO or RATIO_DB.
     */
    enum ValueType_t : uint8_t {
        VALUETYPE_ANY = 0,
        VALUE,          /**< Can be any type of value, i.e. ADC samples. */
        VALUE_CMPX,     /**< Can be any type of value, i.e. ADC samples. */
        RATIO,          /**< Name is self-explanatory. */
        RATIO_DB,       /**< = 10 log (Ratio) or 20 log (Ratio) */
        FREQUENCY_HZ,   /**< duh... */
        FREQUENCY_RAD,  /**< = FrequencyHz * 2 * pi */
        PHASE,          /**< in Radians */
        MAGNITUDE,      /**< Kinda like Value */
        RS_CORR_OUT,    /**< RSCorrOut Struct */
        SS_CORR_OUT,    /**< Signal Synchronous SSCorrOut struct. */
        CHANNEL,		/**< Parameters for tracking a specific PRN in receiver channel. */
        STATE,			/**< An array representing a state column vector. */
        COVARIANCE,		/**< A row-major matrix representing a covariance. */
        FUNCTION_PTR,	/**< Function pointer. */
        EPHEMS,			/**< Ephemeris parameters. */
        GRID			/**< A special type for the manifold grid. */
    };

    /** \brief Memory Location*/
    enum MemLoc_t : uint8_t {
        HOST = 0,
        CUDA_DEVICE = 1
    };

    /** \brief Datatype */
    enum DataType_t : uint8_t {
        DATATYPE_ANY = 0,
        UNDEFINED_t,
        FLOAT_t,
        DOUBLE_t,
        FIXED_Q15_t,    /**< 15 bits decimal with 1 sign bit stored in an int16 datatype. */
        FIXED_Q31_t,    /**< 31 bits decimal with 1 sign bit stored in an int32 datatype. */
        FIXED_I15Q16_t, /**< 15 bits integer, 1 sign bit, 16 bits decimal, stored in int32. */
        CHAR_t,
        STRING_t,
        INT_t,
        BOOL_t,
        CUFFTCOMP_t		/**< Two elements, x and y, both are double. */
    };

    typedef struct Param {
        void *Ptr;
        DataType_t Datatype;
        unsigned int Capacity;  /**< In bytes. */
        unsigned int Size;      /**< In bytes. */
    } Param;

    /** \brief Port data structure for I/O of modules.
     *
     * All member fields other than *Data are fixed at the initialization
     * stage. They are not expected to change dynamically. Modules will treat
     * contents in *Data as how they are described in these fields.
     */
    typedef struct Port {
        char Name[32];

        /** \brief Datatype of the data in this port. */
        DataType_t Datatype;

        /** \brief Actual value = *Data * 10 ^ Exponent
         * probably more important for fixed point types
         * should set this value to 0 when using float?
         */
        signed char Exponent;

        /** \brief The quantity that this value represents.
         *  For example: frequency, phase, ratio, dB etc.
         *  This will allow for sanity checks in flow diagram.
         *  Note that this does not represent the meaning of the value. For
         *  example, an SNR value will be of type RATIO or RATIO_DB.
         */
        ValueType_t ValueType;

        /** \brief Location of memory pointer */
        MemLoc_t MemLoc;

        /** \brief Length on the incoming data array.
         *  1=scalar, otherwise, length of array.
         */
        unsigned short VectorLength;

        /** \brief Pointer to actual data.
         * Can be any datatype but typecast to (void*) to be generic for all
         * modules
         */
        void *Data;

        /** \brief Aux value for any module specific implementation.
         *  This value can only be set once during initialization stage.
         */
        int AuxValue;
    } Port;

    /** \brief Expected Port data structure for I/O of modules.
     *
     * For error checking
     */
    typedef struct ExpectedPort {
        char Name[32];

        /** \brief Datatype of the data in this port. */
        DataType_t Datatype;

        /** \brief The quantity that this value represents.
         *  For example: frequency, phase, ratio, dB etc.
         *  This will allow for sanity checks in flow diagram.
         *  Note that this does not represent the meaning of the value. For
         *  example, an SNR value will be of type RATIO or RATIO_DB.
         */
        ValueType_t ValueType;

        /** \brief Length on the incoming data array.
         *  1=scalar, 0=vector of any length > 1, otherwise, length of array.
         */
        unsigned short VectorLength;

    } ExpectedPort;


    /** \brief Pointer to function pointer for EKF motion model (f) */
	typedef void (*fFunc_t)(double*, double*, double*, double*);

	/** \brief Pointer to function pointer for EKF measurement model (h) */
	typedef void (*hFunc_t)(double*, double*, double*);


	/**
	 * Class Declarations for Flows and Modules
	 */

	// Template Flow class
    class Flow; 
    
    // Flow helper functions that is friended with Flow class
    class FlowMgr;

	// Specific Flows
    class DPEFlow;

    
    
	// Template Module class
    class Module;

	// Specific Modules
	class DPInit;
	
    class SampleBlock;
    
    class BatchCorrScores;

    class BatchCorrManifold;

    class Correlator;

    class LoopFilter;

    class PhaseDiscriminator;

    class DelayDiscriminator;

    class LPF;

    class LockDetector;

    class cuEKF;

    class cuChanMgr;

    class SatPos;

    class ScalarControl;

    class DPControl;

    class DataLogger;

    class Acquisition;

    class AcqControl;

    class SampleServer;
}

#endif
