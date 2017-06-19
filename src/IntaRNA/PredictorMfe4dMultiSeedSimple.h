
#ifndef INTARNA_PREDICTORMFE4DMULTISEED_H_
#define INTARNA_PREDICTORMFE4DMULTISEED_H_

#include "IntaRNA/PredictorMfe4d.h"
#include "IntaRNA/SeedConstraint.h"
#include "IntaRNA/SeedHandlerIdxOffset.h"

namespace IntaRNA {

/**
 * Predictor for exact interaction prediction with seed constraint using a
 * 4D matrix.
 *
 * This enables non-overlapping suboptimal enumeration.
 *
 * @author Martin Mann
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
						, const SeedConstraint & seedConstraint
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

	//! energy of all multi-mode interaction hybrids that contain a seeded
	//! interaction right of the gap and no seed left of it.
	//! they are computed by the recursion with indices
	//! hybridE_seed(i1,i2)->(w1,w2), with interaction start i1 (seq1) and i2 (seq2) and
	//! interaction end j1=i1+w1 and j2=j2+w2.
	//! NOTE: hybridE_seed(i1,i2)==NULL if not complementary(seq1[i1],seq2[i2])
	E4dMatrix hybridE_multi;


protected:

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

};

} // namespace
#endif /* INTARNA_PREDICTORMFE4DMULTISEED_H_ */
