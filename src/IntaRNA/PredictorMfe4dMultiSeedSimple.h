
#ifndef INTARNA_PREDICTORMFE4DMULTISEEDSIMPLE_H
#define INTARNA_PREDICTORMFE4DMULTISEEDSIMPLE_H

#include "IntaRNA/PredictorMfe4d.h"
#include "IntaRNA/SeedConstraint.h"
#include "IntaRNA/SeedHandlerIdxOffset.h"

namespace IntaRNA {

/**
 * Predictor for exact multi-site interaction prediction with seed constraint using a
 * 4D matrix.
 * 
 * The recursion is a straight-forward extension with O(n^6) time consumption.
 *
 * This enables non-overlapping suboptimal enumeration.
 *
 * @author Martin Mann
 * @author Rick Gelhausen
 *
 */
class PredictorMfe4dMultiSeedSimple: public PredictorMfe4d {


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
	 * @param seedConstraint the seed constraint to be used for seed identification
	 */
	PredictorMfe4dMultiSeedSimple( const InteractionEnergy & energy
						, OutputHandler & output
						, PredictionTracker * predTracker
						, const AllowES allowES
						, SeedHandler * seedHandlerInstance
						);

	/**
	 * destruction
	 */
	virtual ~PredictorMfe4dMultiSeedSimple();

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
	using PredictorMfe4d::energy;

	//! access to the output handler of the super class
	using PredictorMfe4d::output;

	//! access to the list of reported interaction ranges of the super class
	using PredictorMfe4d::reportedInteractions;

	//! energy of all interaction hybrids computed by the recursion with indices
	//! hybridE(i1,i2)->(w1,w2), with interaction start i1 (seq1) and i2 (seq2) and
	//! interaction end j1=i1+w1 and j2=j2+w2. Interactions do not necessarily
	//! contain a seed interaction.
	//! NOTE: hybridE(i1,i2)==NULL if not complementary(seq1[i1],seq2[i2])
	using PredictorMfe4d::hybridE;


	//! the seed handler (with idx offset)
	SeedHandlerIdxOffset seedHandler;

	//! energy of all interaction hybrids that contain a seed interaction.
	//! they are computed by the recursion with indices
	//! hybridE_seed(i1,i2)->(w1,w2), with interaction start i1 (seq1) and i2 (seq2) and
	//! interaction end j1=i1+w1 and j2=j2+w2.
	//! NOTE: hybridE_seed(i1,i2)==NULL if not complementary(seq1[i1],seq2[i2])
	E4dMatrix hybridE_seed;

	//! defines where ES-terms are considered
    AllowES allowES;


protected:

	/**
		 * disables tracker update for mfe updates
		 *
		 * @param i1 the index of the first sequence interacting with i2
		 * @param j1 the index of the first sequence interacting with j2
		 * @param i2 the index of the second sequence interacting with i1
		 * @param j2 the index of the second sequence interacting with j1
		 * @param energy ignored
		 * @param isHybridE ignored
		 */
		virtual
		void
		updateOptima( const size_t i1, const size_t j1
				, const size_t i2, const size_t j2
				, const E_type energy
				, const bool isHybridE );

		/**
	 * Removes all temporary data structures and resets the predictor
	 */
	void
	clear();

	/**
	 * computes all entries of the hybridE matrix
	 */
	void
	fillHybridE_seed( );

	/**
	 * Fills a given interaction (boundaries given) with the according
	 * hybridizing base pairs.
	 * @param interaction IN/OUT the interaction to fill
	 */
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
	void
	getNextBest( Interaction & curBest );
};


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

inline
void
PredictorMfe4dMultiSeedSimple::
updateOptima( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2
		, const E_type energy
		, const bool isHybridE )
{
	// temporarily disable tracker
	PredictionTracker * curPredTracker = this->predTracker;
	this->predTracker = NULL;
	// update optimum information
	PredictorMfe4d::updateOptima(i1,j1,i2,j2,energy,isHybridE);
	// reenable tracker
	this->predTracker = curPredTracker;
}

//////////////////////////////////////////////////////////////////////////


} // namespace
#endif /* INTARNA_PREDICTORMFE4DMULTISEEDSIMPLE_H */
