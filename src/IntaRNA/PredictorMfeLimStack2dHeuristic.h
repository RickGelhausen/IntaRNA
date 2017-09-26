
#ifndef INTARNA_PREDICTORMFELIMSTACK2DHEURISTIC_H_
#define INTARNA_PREDICTORMFELIMSTACK2DHEURISTIC_H_

#include "IntaRNA/PredictorMfe2dHeuristic.h"
#include "IntaRNA/Interaction.h"

#include <boost/numeric/ublas/matrix.hpp>

namespace IntaRNA {

/**
 * Memory efficient interaction predictor that uses a heuristic to
 * find the mfe or a close-to-mfe interaction.
 *
 * To this end, for each interaction start i1,i2 only the optimal right side
 * interaction with boundaries j1,j2 is considered in the recursion instead of
 * all possible interaction ranges.
 *
 * Helices (continues stackings of base pairs) are restricted to a given length.
 *
 * This yields a quadratic time and space complexity.
 *
 * @author Martin Mann
 *
 */
class PredictorMfeLimStack2dHeuristic: public PredictorMfe2dHeuristic {

protected:

	//! matrix type to hold the mfe energies and boundaries for interaction site starts
	typedef PredictorMfe2dHeuristic::E2dMatrix E2dMatrix;

public:

	/**
	 * Constructs a predictor and stores the energy and output handler
	 *
	 * @param energy the interaction energy handler
	 * @param output the output handler to report mfe interactions to
	 * @param predTracker the prediction tracker to be used or NULL if no
	 *         tracking is to be done; if non-NULL, the tracker gets deleted
	 *         on this->destruction.
	 */
	PredictorMfeLimStack2dHeuristic( const InteractionEnergy & energy
							, OutputHandler & output
							, PredictionTracker * predTracker );

	virtual ~PredictorMfeLimStack2dHeuristic();

protected:

	//! access to the interaction energy handler of the super class
	using PredictorMfe2dHeuristic::energy;

	//! access to the output handler of the super class
	using PredictorMfe2dHeuristic::output;

	//! access to the list of reported interaction ranges of the super class
	using PredictorMfe2dHeuristic::reportedInteractions;

	//! energy of all interaction hybrids starting in i1,i2
	using PredictorMfe2dHeuristic::hybridE;

	//! maximal number of continues stackings within an interaction
	size_t maxHelixLength;

	// TODO overwrite "predict()" to enable correct verbose output

protected:

	/**
	 * Computes all entries of the hybridE matrix
	 */
	virtual
	void
	fillHybridE();

	/**
	 * Fills a given interaction (boundaries given) with the according
	 * hybridizing base pairs.
	 * @param interaction IN/OUT the interaction to fill
	 */
	virtual
	void
	traceBack( Interaction & interaction );


};

} // namespace

#endif /* INTARNA_PREDICTORMFELIMSTACK2DHEURISTIC_H_ */
