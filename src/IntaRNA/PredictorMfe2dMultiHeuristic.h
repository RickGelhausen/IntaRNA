#ifndef INTARNA_PredictorMfe2dMultiHeuristic_H
#define INTARNA_PredictorMfe2dMultiHeuristic_H

#include "IntaRNA/PredictorMfe.h"
#include "IntaRNA/Interaction.h"

#include <boost/numeric/ublas/matrix.hpp>

namespace IntaRNA {
/**
 * Predictor for RNAup-like computation, i.e. full DP-implementation without
 * seed-heuristic using a 2D matrix
 *
 * @author Martin Mann & Rick Gelhausen
 *
 */
class PredictorMfe2dMultiHeuristic: public PredictorMfe {

protected:

	/**
	 * Describes the currently best interaction found for a left interaction
	 * boundary i1,i2
	 */
	class BestInteraction {
	public:

		//! init data
		BestInteraction( const E_type E=E_INF, const size_t j1=RnaSequence::lastPos, const size_t j2=RnaSequence::lastPos )
				: E(E), j1(j1), j2(j2)
		{}

	public:
		//! energy of the interaction
		E_type E;
		//! right end of the interaction in seq1
		size_t j1;
		//! right end of the interaction in seq2
		size_t j2;
	};

	//! matrix type to hold the mfe energies and boundaries for interaction site starts
	typedef boost::numeric::ublas::matrix<BestInteraction> E2dMatrix;

public:

	/**
	 * Constructs a predictor and stores the energy and output handler
	 *
	 * @param energy the interaction energy handler
	 * @param output the output handler to report mfe interactions to
	 * @param predTracker the prediction tracker to be used or NULL if no
	 *         tracking is to be done; if non-NULL, the tracker gets deleted
	 *         on this->destruction.
	 * @param allowES where ES-terms are to be considered
	 */
	PredictorMfe2dMultiHeuristic( const InteractionEnergy & energy
			, OutputHandler & output
			, PredictionTracker * predTracker
			, const AllowES allowES);

	virtual ~PredictorMfe2dMultiHeuristic();

	/**
	 * Computes the mfe for the given sequence ranges (i1-j1) in the first
	 * sequence and (i2-j2) in the second sequence and reports it to the output
	 * handler.
	 *
	 * @param r1 the index range of the first sequence interacting with r2
	 * @param r2 the index range of the second sequence interacting with r1
	 * @param outConstraint constrains the interactions reported to the output handler
	 *
	 */
	virtual
	void
	predict( const IndexRange & r1 = IndexRange(0,RnaSequence::lastPos)
			, const IndexRange & r2 = IndexRange(0,RnaSequence::lastPos)
			, const OutputConstraint & outConstraint = OutputConstraint() );

protected:

	//! access to the interaction energy handler of the super class
	using PredictorMfe::energy;

	//! access to the output handler of the super class
	using PredictorMfe::output;

	//! access to the list of reported interaction ranges of the super class
	using PredictorMfe::reportedInteractions;

	// TODO provide all data structures as arguments to make predict() call threadsafe

	//! energy of all interaction hybrids starting in i1,i2
	E2dMatrix hybridE;

	//! Auxillary Matrix
	//! Composed of the ES2 values for a fixed value in S1 and hybridE of the remaining part.
	E2dMatrix hybridO;

	//! defines where ES-terms are considered
	AllowES allowES;

protected:

	/**
	 * computes all entries of the hybridE matrix
	 */
	virtual
	void
	fillHybridE( );

	/**
	 * Recurse into HybridO to find the index k2 for which k1 returns the minimal energy contribution.
	 * @return index k2
	 */
	virtual
	size_t traceHybridO(const size_t i1, const size_t j1,
						const size_t i2, const size_t j2) const;

	/**
	 * Fills a given interaction (boundaries given) with the according
	 * hybridizing base pairs.
	 * @param interaction IN/OUT the interaction to fill
	 */
	virtual
	void
	traceBack( Interaction & interaction );

	/**
	 * Identifies the next best interaction with an energy equal to or higher
	 * than the given interaction. The new interaction will not overlap any
	 * index range stored in reportedInteractions.
	 *
	 * @param curBest IN/OUT the current best interaction to be replaced with one
	 *        of equal or higher energy not overlapping with any reported
	 *        interaction so far; an interaction with energy E_INF is set, if
	 *        there is no better interaction left
	 */
	virtual
	void
	getNextBest( Interaction & curBest );
};

} // namespace

#endif //INTARNA_PredictorMfe2dMultiHeuristic_H
