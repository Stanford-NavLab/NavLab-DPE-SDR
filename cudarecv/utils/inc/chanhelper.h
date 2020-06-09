#ifndef __INC_UTILS_CHANHELPER_H_
#define __INC_UTILS_CHANHELPER_H_

namespace dsp { namespace utils {

	/** \brief Correlation sums. */
	template <typename T>
	struct CorrOut {
		T iE;           /**< Early */
		T qE;           /**< Early */
		T iP;           /**< Prompt */
		T qP;           /**< Prompt */
		T iL;           /**< Late */
		T qL;           /**< Late */
	};

	/** \brief Receiver Synchronous Correlation sums for 1ms */
	template <typename T>
	struct RSCorrOut {
		unsigned char PRN;
		CorrOut<T> Corr;
	};

	/** \brief Signal Synchronous correlation sums for 1 ms */
	template <typename T>
	struct SSCorrOut {
		unsigned char PRN;
		unsigned char NumCorr;  /**< Number of correlations in 1ms. */
		CorrOut<T> Corr[3];     /**< Max 3 possible correlations in 1 ms. */
	};

	/** \brief Channel parameters */
	struct Channel {
		double rc; /**< Code phase (# of chips into sequence). */
		double ri; /**< Carrier phase. */
		double fc; /**< Code frequency (Doppler effects). */
		double fi; /**< Carrier frequency (Doppler effects). */
		unsigned long cp; /**< Number of complete PRN codes elapsed since start of data. */
		unsigned long cp_timestamp; /**< Number of complete PRN codes elapsed between start of data and reference subframe (first acquired). */
	};

}
}

#endif
