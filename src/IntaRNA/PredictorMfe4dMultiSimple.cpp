#include "IntaRNA/PredictorMfe4dMultiSimple.h"

namespace IntaRNA {

//////////////////////////////////////////////////////////

PredictorMfe4dMultiSimple::
PredictorMfe4dMultiSimple( const InteractionEnergy & energy
        , OutputHandler & output
        , PredictionTracker * predTracker
        , const AllowES allowES_)
        : PredictorMfe4d(energy,output, predTracker)
        , allowES( allowES_ )
{
}


////////////////////////////////////////////////////////////////////////////

PredictorMfe4dMultiSimple::
~PredictorMfe4dMultiSimple()
{
    // clean up
    this->clear();
}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfe4dMultiSimple::
predict( const IndexRange & r1
        , const IndexRange & r2
        , const OutputConstraint & outConstraint
)
{
    VLOG(2) <<"predicting mfe interactions in O(n^4) space and O(n^5) time...";
    // measure timing
    TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

#if IN_DEBUG_MODE
    // check indices
    if (!(r1.isAscending() && r2.isAscending()) )
        throw std::runtime_error("PredictorMfe4d::predict("+toString(r1)+","+toString(r2)+") is not sane");
#endif

    // clear data
    clear();

    // setup index offset
    energy.setOffset1(r1.from);
    energy.setOffset2(r2.from);

    // resize matrix
    hybridE.resize( std::min( energy.size1()
            , (r1.to==RnaSequence::lastPos?energy.size1()-1:r1.to)-r1.from+1 )
            , std::min( energy.size2()
                    , (r2.to==RnaSequence::lastPos?energy.size2()-1:r2.to)-r2.from+1 ) );

    size_t w1, w2;


    bool i1blocked, i1or2blocked, skipw1w2;
    // initialize 3rd and 4th dimension of the matrix
    for (size_t i1=0; i1<hybridE.size1(); i1++) {
        // check if i1 is blocked for interaction
        i1blocked = !energy.isAccessible1(i1);
        for (size_t i2=0; i2<hybridE.size2(); i2++) {
            // check whether i1 or i2 is blocked for interaction
            i1or2blocked = i1blocked || !energy.isAccessible2(i2);

            // check if i1 and i2 are not blocked and can form a base pair
            if ( ! i1or2blocked
                 && energy.areComplementary( i1, i2 ))
            {
                // create new 2d matrix for different interaction site widths
                hybridE(i1,i2) = new E2dMatrix(
                        /*w1 = */ std::min(energy.getAccessibility1().getMaxLength(), hybridE.size1()-i1 ),
                        /*w2 = */ std::min(energy.getAccessibility2().getMaxLength(), hybridE.size2()-i2 ));
            } else {
                // reduce memory consumption and avoid computation for this start index combination
                hybridE(i1,i2) = NULL;
            }
        }
    }
    // initialize mfe interaction for updates
    initOptima( outConstraint );

    // fill matrix
    fillHybridE( );

    // report mfe interaction
    reportOptima( outConstraint );
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe4dMultiSimple::
clear()
{
    // delete 3rd and 4th dimension of the matrix
    for (E4dMatrix::iterator1 ijEntry = hybridE.begin1(); ijEntry != hybridE.end1(); ijEntry++) {
        if (*ijEntry != NULL) {
            // delete 2d matrix for current ij
            delete (*ijEntry);
            *ijEntry = NULL;
        }
    }
    // clear matrix, free data
    hybridE.clear();
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe4dMultiSimple::
fillHybridE( ) {
    // global vars to avoid reallocation
    size_t i1, i2, j1, j2, w1, w2, k1, k2;

    //////////  COMPUTE HYBRIDIZATION ENERGIES ONLY INCLUDING ES-TERMS  ////////////

    // current minimal value
    E_type curMinE = E_INF;
    // iterate increasingly over all window sizes w1 (seq1) and w2 (seq2)
    for (w1 = 0; w1 < energy.getAccessibility1().getMaxLength(); w1++) {
        for (w2 = 0; w2 < energy.getAccessibility2().getMaxLength(); w2++) {
            // iterate over all window starts i1 (seq1) and i2 (seq2)
            // TODO PARALLELIZE THIS DOUBLE LOOP ?!
            for (i1 = 0; i1 + w1 < hybridE.size1(); i1++) {
                for (i2 = 0; i2 + w2 < hybridE.size2(); i2++) {

                    // check if left boundary is complementary
                    if (hybridE(i1, i2) == NULL) {
                        // interaction not possible: nothing to do, since no storage reserved
                        continue;
                    }

                    // and widths are possible (ie available within data structure)
                    if (hybridE(i1, i2)->size1() <= w1 || hybridE(i1, i2)->size2() <= w2) {
                        // interaction not possible: nothing to do, since no storage reserved
                        continue;
                    }


                    // get window ends j (seq1) and l (seq2)
                    j1 = i1 + w1;
                    j2 = i2 + w2;

                    // check if right boundary is complementary
                    if (hybridE(j1, j2) == NULL) {
                        // not complementary -> ignore this entry
                        (*hybridE(i1, i2))(w1, w2) = E_INF;
                        continue;
                    }

                    // compute entry, since (i1,i2) complementary
                    {
                        // init
                        curMinE = E_INF;

                        // either interaction initiation
                        if (w1 == 0 && w2 == 0) {
                            // interaction initiation
                            curMinE = energy.getE_init();
                        }
                        else
                            // or loop
                            // or ES-gap
                        if (w1 > 0 && w2 > 0)
                        {
                            // interior loop case

                            // check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
                            for (k1 = std::min(j1, i1 + energy.getMaxInternalLoopSize1() + 1); k1 > i1; k1--) {
                                for (k2 = std::min(j2, i2 + energy.getMaxInternalLoopSize2() + 1); k2 > i2; k2--) {
                                    // check if (k1,k2) are complementary
                                    if (hybridE(k1, k2) != NULL && hybridE(k1, k2)->size1() > (j1 - k1) &&
                                        hybridE(k1, k2)->size2() > (j2 - k2))
                                    {
                                        // update minE
                                        curMinE = std::min(curMinE,
                                                           (energy.getE_interLeft(i1, k1, i2, k2)
                                                            + (*hybridE(k1, k2))(j1 - k1, j2 - k2)));
                                    }
                                }
                            }

                            // Multiloop cases = ES-gap

                            // Both-sided structure
                            if (allowES == ES_both) {
                                for (k1 = j1; k1 > i1 + InteractionEnergy::minDistES; k1--) {
                                for (k2 = j2; k2 > i2 + InteractionEnergy::minDistES; k2--) {
                                    if (hybridE(k1, k2) != NULL
                                        && hybridE(k1, k2)->size1() > (j1 - k1)
                                        && hybridE(k1, k2)->size2() > (j2 - k2))
                                    {
                                        // update minE
                                        curMinE = std::min(curMinE,
                                                           (energy.getE_multi(i1, k1, i2, k2, InteractionEnergy::ES_multi_mode::ES_multi_both)
                                                            + (*hybridE(k1, k2))(j1 - k1, j2 - k2)
                                                           ));
                                    }
                                    }
                                }
                            }

                            // Structure in S1
                            if (allowES == ES_target || allowES == ES_xorQueryTarget) {
                                for (k1 = j1; k1 > i1 + InteractionEnergy::minDistES; k1--) {
                                for (k2 = std::min(j2, i2 + energy.getMaxInternalLoopSize2() + 1);
                                     k2 > i2; k2--) {
                                    if (hybridE(k1, k2) != NULL
                                        && hybridE(k1, k2)->size1() > (j1 - k1)
                                        && hybridE(k1, k2)->size2() > (j2 - k2))
                                    {
                                        // update minE
                                        curMinE = std::min(curMinE,
                                                           (energy.getE_multi(i1, k1, i2, k2, InteractionEnergy::ES_multi_mode::ES_multi_1only)
                                                            + (*hybridE(k1, k2))(j1 - k1, j2 - k2)
                                                           ));
                                    }
                                    }
                                }
                            }

                            // Structure in S2
                            if (allowES == ES_query || allowES == ES_xorQueryTarget) {
                                for (k1 = std::min(j1, i1 + energy.getMaxInternalLoopSize1() + 1); k1 > i1; k1--) {
                                for (k2 = j2; k2 > i2 + InteractionEnergy::minDistES; k2--) {
                                    if (hybridE(k1, k2) != NULL
                                        && hybridE(k1, k2)->size1() > (j1 - k1)
                                        && hybridE(k1, k2)->size2() > (j2 - k2))
                                    {
                                        // update minE
                                        curMinE = std::min(curMinE,
                                                           (energy.getE_multi(i1, k1, i2, k2, InteractionEnergy::ES_multi_mode::ES_multi_2only)
                                                            + (*hybridE(k1, k2))(j1 - k1, j2 - k2)
                                                           ));
                                    }
                                }
                                }
                            }
                        }
                        // store value
                        (*hybridE(i1, i2))(w1, w2) = curMinE;

                        // update mfe if needed
                        updateOptima(i1, j1, i2, j2, (*hybridE(i1, i2))(w1, w2), true);
                        continue;
                    }
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe4dMultiSimple::
traceBack( Interaction & interaction )
{
    // check if something to trace
    if (interaction.basePairs.size() < 2) {
        return;
    }

#if IN_DEBUG_MODE
    // sanity checks
    if ( ! interaction.isValid() ) {
        throw std::runtime_error("PredictorMfe4d::traceBack() : given interaction not valid");
    }
    if ( interaction.basePairs.size() != 2 ) {
        throw std::runtime_error("PredictorMfe4d::traceBack() : given interaction does not contain boundaries only");
    }
#endif

    // check for single interaction
    if (interaction.basePairs.at(0).first == interaction.basePairs.at(1).first) {
        // delete second boundary (identical to first)
        interaction.basePairs.resize(1);
        // update done
        return;
    }

    // ensure sorting
    interaction.sort();
    // get indices in hybridE for boundary base pairs
    size_t	i1 = energy.getIndex1(interaction.basePairs.at(0)),
            j1 = energy.getIndex1(interaction.basePairs.at(1)),
            i2 = energy.getIndex2(interaction.basePairs.at(0)),
            j2 = energy.getIndex2(interaction.basePairs.at(1));

    // the currently traced value for i1-j1, i2-j2
    E_type curE = (*hybridE(i1,i2))(j1-i1,j2-i2);

    // trace back
    // do until only right boundary is left over
    while( i1 != j1 ) {
        // temp variables
        size_t k1,k2;
        bool traceNotFound = true;

        // check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
        for (k1=std::min(j1,i1+energy.getMaxInternalLoopSize1()+1); traceNotFound && k1>i1; k1--) {
		for (k2=std::min(j2,i2+energy.getMaxInternalLoopSize2()+1); traceNotFound && k2>i2; k2--) {
			// check if (k1, k2) are complementary
			if (hybridE(k1,k2) != NULL
				&& hybridE(k1,k2)->size1() > (j1-k1)
				&& hybridE(k1,k2)->size2() > (j2-k2))
			{

				if ( E_equal( curE,
							  (energy.getE_interLeft(i1,k1,i2,k2)
							   + (*hybridE(k1,k2))(j1-k1,j2-k2)
							  ) ) )
				{
					// stop searching
					traceNotFound = false;
					// store splitting base pair
					interaction.basePairs.push_back( energy.getBasePair(k1,k2) );
					// trace right part of split
					i1=k1;
					i2=k2;
					curE = (*hybridE(i1,i2))(j1-i1,j2-i2);
				}
			}
		}
        }

        // Structure in both
        if (traceNotFound && allowES == ES_both) {
            for (k1 = j1; traceNotFound && k1 > i1 + InteractionEnergy::minDistES; k1--) {
            for (k2 = j2; traceNotFound && k2 > i2 + InteractionEnergy::minDistES; k2--) {
				if (hybridE(k1, k2) != NULL
					&& hybridE(k1, k2)->size1() > (j1 - k1)
					&& hybridE(k1, k2)->size2() > (j2 - k2))
				{
					if (E_equal(curE,
						  (energy.getE_multi(i1, k1, i2, k2, InteractionEnergy::ES_multi_mode::ES_multi_both)
						   + (*hybridE(k1, k2))(j1 - k1, j2 - k2)
						  ))) {
					  // stop searching
					  traceNotFound = false;
					  // store splitting base pair
					  interaction.basePairs.push_back(energy.getBasePair(k1, k2));
					  // store gap information
					  if (interaction.gap == NULL) { interaction.gap = new Interaction::Gap(); }
					  interaction.gap->energy += energy.getE_multi(i1, k1, i2, k2, InteractionEnergy::ES_multi_mode::ES_multi_both);
					  Interaction::BasePair bpLeft = energy.getBasePair(i1,i2);
					  interaction.gap->gaps1.push_back( IndexRange(bpLeft.first,interaction.basePairs.rbegin()->first) );
					  interaction.gap->gaps2.push_back( IndexRange(interaction.basePairs.rbegin()->second,bpLeft.second) );
					  // trace right part of split
					  i1 = k1;
					  i2 = k2;
					  curE = (*hybridE(i1, i2))(j1 - i1, j2 - i2);
					}
				}
            }
            }
        }

        // Structure in S1
        if (traceNotFound && (allowES == ES_target || allowES == ES_xorQueryTarget)) {
            for (k1 = j1; traceNotFound && k1 > i1 + InteractionEnergy::minDistES; k1--) {
            for (k2 = std::min(j2, i2 + energy.getMaxInternalLoopSize2() + 1); traceNotFound && k2 > i2; k2--) {
				if (hybridE(k1, k2) != NULL
					&& hybridE(k1, k2)->size1() > (j1 - k1)
					&& hybridE(k1, k2)->size2() > (j2 - k2))
				{
					if (E_equal(curE,
						  (energy.getE_multi(i1, k1, i2, k2, InteractionEnergy::ES_multi_mode::ES_multi_1only)
						   + (*hybridE(k1, k2))(j1 - k1, j2 - k2)
						  ))) {
					  // stop searching
					  traceNotFound = false;
					  // store splitting base pair
					  interaction.basePairs.push_back(energy.getBasePair(k1, k2));
					  // store gap information
					  if (interaction.gap == NULL) { interaction.gap = new Interaction::Gap(); }
					  interaction.gap->energy += energy.getE_multi(i1, k1, i2, k2, InteractionEnergy::ES_multi_mode::ES_multi_1only);
					  Interaction::BasePair bpLeft = energy.getBasePair(i1,i2);
					  interaction.gap->gaps1.push_back( IndexRange(bpLeft.first,interaction.basePairs.rbegin()->first) );

					  // trace right part of split
					  i1 = k1;
					  i2 = k2;
					  curE = (*hybridE(i1, i2))(j1 - i1, j2 - i2);
					}
				}
			}
			}
        }

    // Structure in S2
        if (traceNotFound && (allowES == ES_query || allowES == ES_xorQueryTarget)) {
            for (k1 = std::min(j1, i1 + energy.getMaxInternalLoopSize1() + 1); traceNotFound && k1 > i1; k1--) {
            for (k2 = j2; traceNotFound && k2 > i2 + InteractionEnergy::minDistES; k2--) {
				if (hybridE(k1, k2) != NULL
					&& hybridE(k1, k2)->size1() > (j1 - k1)
					&& hybridE(k1, k2)->size2() > (j2 - k2))
				{
					if (E_equal(curE,
						  (energy.getE_multi(i1, k1, i2, k2, InteractionEnergy::ES_multi_mode::ES_multi_2only)
						   + (*hybridE(k1, k2))(j1 - k1, j2 - k2)
						  ))) {
					  // stop searching
					  traceNotFound = false;
					  // store splitting base pair
					  interaction.basePairs.push_back(energy.getBasePair(k1, k2));
					  // store gap information
					  if (interaction.gap == NULL) { interaction.gap = new Interaction::Gap(); }
					  interaction.gap->energy += energy.getE_multi(i1, k1, i2, k2, InteractionEnergy::ES_multi_mode::ES_multi_2only);
					  Interaction::BasePair bpLeft = energy.getBasePair(i1,i2);
					  interaction.gap->gaps2.push_back( IndexRange(interaction.basePairs.rbegin()->second,bpLeft.second) );
					  // trace right part of split
					  i1 = k1;
					  i2 = k2;
					  curE = (*hybridE(i1, i2))(j1 - i1, j2 - i2);
					}
				}
            }
            }
        }
		// final sanity check
		assert(!traceNotFound);
    }

    // sort final interaction (to make valid) (faster than calling sort())
    if (interaction.basePairs.size() > 2) {
        Interaction::PairingVec & bps = interaction.basePairs;
        // check if last added base pair is a duplicate of the boundary
        if (bps.rbegin()->first == j1) {
            // remove last base pair
            bps.erase( bps.begin()+bps.size()-1 );
        }
        // shift all added base pairs to the front
        for (size_t i=2; i<bps.size(); i++) {
            bps.at(i-1).first = bps.at(i).first;
            bps.at(i-1).second = bps.at(i).second;
        }
        // set last to j1-j2
        (*bps.rbegin()) = energy.getBasePair( j1, j2 );
    }

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe4dMultiSimple::
getNextBest( Interaction & curBest )
{
    // store current best energy.second
    const E_type curBestE = curBest.energy;

    // overwrite energy for update
    curBest.energy = E_INF;
    curBest.basePairs.resize(2);

    // TODO replace index iteration with something based on ranges from reportedInteractions

    // identify cell with next best non-overlapping interaction site
    // iterate (decreasingly) over all left interaction starts
    E2dMatrix * curTable = NULL;
    IndexRange r1,r2;
	size_t d1 = 0, d2 = 0;  // temp vars to deals with possible interaction lengths
    E_type curE = E_INF;
    for (r1.from=hybridE.size1(); r1.from-- > 0;) {

        // ensure interaction site start is not covered
        if (reportedInteractions.first.covers(r1.from)) {
            continue;
        }

        for (r2.from=hybridE.size2(); r2.from-- > 0;) {

            // ensure interaction site start is not covered
            if (reportedInteractions.second.covers(r2.from)) {
                continue;
            }
            // check if left boundary is complementary
            if (hybridE(r1.from,r2.from) == NULL) {
                // interaction not possible: nothing to do, since no storage reserved
                continue;
            }

			// access energy table for left-most interaction base pair
            curTable = hybridE(r1.from,r2.from);

			// iterate over all available interaction site lengths in seq1
			for (d1 = 0; d1<curTable->size1(); d1++) {

				// set according right interaction boundary in seq1
				r1.to = r1.from + d1;

                // check of overlapping
                if (reportedInteractions.first.overlaps(r1)) {
                    // stop since all larger sites will overlap as well
                    break;;
                }

				// iterate over all available interaction site lengths in seq2
				for (d2=0; d2<curTable->size2(); d2++) {

					// set according right interaction boundary in seq2
					r2.to = r2.from + d2;

                    // check of overlapping
                    if (reportedInteractions.second.overlaps(r2)) {
                        // stop since all larger sites will overlap as well
                        break;;
                    }

                    // get overall energy of entry
					curE = energy.getE( r1.from, r1.to, r2.from, r2.to, (*curTable)(d1,d2));

                    // skip sites with energy too low
                    // or higher than current best found so far
                    if (  curE< curBestE || curE >= curBest.energy ) {
                        continue;
                    }

                    //// FOUND THE NEXT BETTER SOLUTION
                    // overwrite current best found so far
                    curBest.energy = curE;
                    curBest.basePairs[0] = energy.getBasePair( r1.from, r2.from );
                    curBest.basePairs[1] = energy.getBasePair( r1.to, r2.to );

                } // j2
            } // j1


        } // i2
    } // i1

}

} // namespace