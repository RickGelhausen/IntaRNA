#ifndef INTARNA_PredictorMfe4dMultiPlus_H_
#define INTARNA_PredictorMfe4dMultiPlus_H_

#include "IntaRNA/PredictorMfe4d.h"

#include <boost/numeric/ublas/matrix.hpp>

namespace IntaRNA {
/**
 * Predictor for RNAup-like computation, i.e. full DP-implementation without
 * seed-heuristic using a 4D matrix
 *
 * @author Martin Mann & Rick Gelhausen
 *
 */
class PredictorMfe4dMultiPlus: public PredictorMfe4d {

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
    PredictorMfe4dMultiPlus( const InteractionEnergy & energy
            , OutputHandler & output
            , PredictionTracker * predTracker
            , const AllowES allowES);

    virtual ~PredictorMfe4dMultiPlus();

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

    // TODO provide all data structures as arguments to make predict() call threadsafe

    //! energy of all interaction hybrids computed by the recursion with indices
    //! hybridE(i1,i2)->(w1,w2), with interaction start i1 (seq1) and i2 (seq2) and
    //! ineraction end j1=i1+w1 and j2=j2+w2
    //! NOTE: hybridE(i1,i2)==NULL if not complementary(seq1[i1],seq2[i2])
    using PredictorMfe4d::hybridE;

    //! Auxillary Matrix
    //! Composed of the ES2 values for a fixed value in S1 and hybridE of remaining part.
    E4dMatrix hybridO;


    //! defines where ES-terms are considered
    AllowES allowES;

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
    fillHybridE( );

    /**
     * Recurse into HybridO to find the index k2 for which k1 returns the minimal energy contribution.
     * @return std::pair<size_t, E_type> index k2 and energy contribution
     */
	size_t traceHybridO(const size_t i1, const size_t j1,
						  const size_t i2, const size_t j2) const;

    /**
     * Fills a given interaction (boundaries given) with the according
     * hybridizing base pairs.
     * @param interaction IN/OUT the interaction to fill
     */
    void
    traceBack( Interaction & interaction );

};

} // namespace
#endif /* INTARNA_PredictorMfe4dMultiPlus_H_ */