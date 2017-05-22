
#include "IntaRNA/PredictorMfe4dMultiSeed.h"

#include <stdexcept>
#include <algorithm>

namespace IntaRNA {

////////////////////////////////////////////////////////////////////////////

PredictorMfe4dMultiSeed::
PredictorMfe4dMultiSeed( const InteractionEnergy & energy
					, OutputHandler & output
					, PredictionTracker * predTracker
					, const AllowES allowES_
					, const SeedConstraint & seedConstraint
					)
 : PredictorMfe4d(energy,output,predTracker)
	, seedHandler(energy,seedConstraint)
	, hybridE_seed(0,0)
	, allowES( allowES_ )
	, hybridE_multi(0,0)
{
}


////////////////////////////////////////////////////////////////////////////

PredictorMfe4dMultiSeed::
~PredictorMfe4dMultiSeed()
{
	// clean up
	this->clear();
}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfe4dMultiSeed::
predict( const IndexRange & r1
		, const IndexRange & r2
		, const OutputConstraint & outConstraint
		)
{
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) <<"predicting mfe multi-side interactions with seed in O(n^4) space and O(n^5) time..."; }
	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

#if IN_DEBUG_MODE
	// check indices
	if (!(r1.isAscending() && r2.isAscending()) )
		throw std::runtime_error("PredictorMfe4dMultiSeed::predict("+toString(r1)+","+toString(r2)+") is not sane");
#endif

	// clear data
	clear();

	// setup index offset
	energy.setOffset1(r1.from);
	energy.setOffset2(r2.from);
	seedHandler.setOffset1(r1.from);
	seedHandler.setOffset2(r2.from);

	const size_t hybridEsize1 = std::min( energy.size1()
			, (r1.to==RnaSequence::lastPos?energy.size1()-1:r1.to)-r1.from+1 );
	const size_t hybridEsize2 = std::min( energy.size2()
			, (r2.to==RnaSequence::lastPos?energy.size2()-1:r2.to)-r2.from+1 );


	// compute seed interactions for whole range
	// and check if any seed possible
	if (seedHandler.fillSeed( 0, hybridEsize1-1, 0, hybridEsize2-1 ) == 0) {
		// trigger empty interaction reporting
		initOptima(outConstraint);
		reportOptima(outConstraint);
		// stop computation
		return;
	}


	// resize matrix
	hybridE.resize( hybridEsize1, hybridEsize2 );
	hybridE_seed.resize( hybridEsize1, hybridEsize2 );
	hybridE_multi.resize( hybridEsize1, hybridEsize2 );

	size_t w1, w2;


	bool i1blocked, i1or2blocked, skipw1w2;
	// initialize 3rd and 4th dimension of the matrix
	for (size_t i1=0; i1<hybridEsize1; i1++) {
		// check if i1 is blocked for interaction
		i1blocked = !energy.isAccessible1(i1);
		for (size_t i2=0; i2<hybridEsize2; i2++) {
			// check whether i1 or i2 is blocked for interaction
			i1or2blocked = i1blocked || !energy.isAccessible2(i2);

			// check if i1 and i2 are not blocked and can form a base pair
			if ( ! i1or2blocked
				&& energy.areComplementary( i1, i2 ))
			{
				// create new 2d matrix for different interaction site widths
				hybridE(i1,i2) = new E2dMatrix(
					/*w1 = */ std::min(energy.getAccessibility1().getMaxLength(), hybridEsize1-i1 ),
					/*w2 = */ std::min(energy.getAccessibility2().getMaxLength(), hybridEsize2-i2 ));
				hybridE_seed(i1,i2) = new E2dMatrix( hybridE(i1,i2)->size1(), hybridE(i1,i2)->size2() );
				hybridE_multi(i1,i2) = new E2dMatrix( hybridE(i1,i2)->size1(), hybridE(i1,i2)->size2() );

			} else {
				// reduce memory consumption and avoid computation for this start index combination
				hybridE(i1,i2) = NULL;
				hybridE_seed(i1,i2) = NULL;
				hybridE_multi(i1,i2) = NULL;
			}
		}
	}

	// init mfe without seed condition
	OutputConstraint tmpOutConstraint(1, outConstraint.reportOverlap, outConstraint.maxE, outConstraint.deltaE);
	initOptima( tmpOutConstraint );

	// fill matrix
	fillHybridE( );

	// check if any interaction possible
	// if not no seed-containing interaction is possible neither
	if (!(this->mfeInteractions.begin()->energy < tmpOutConstraint.maxE)) {
		// stop computation since no favorable interaction found
		reportOptima(tmpOutConstraint);
		return;
	}

	// initialize mfe interaction with seed for updates
	initOptima( outConstraint );
	// fill final matrix
	fillHybridE_seed( );

	// report mfe interaction
	reportOptima( outConstraint );

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe4dMultiSeed::
clear()
{
	// delete 3rd and 4th dimension of the matrix
	for (E4dMatrix::iterator1 iRows = hybridE_seed.begin1(); iRows != hybridE_seed.end1(); iRows++) {
		for (E4dMatrix::iterator2 ijEntry = iRows.begin(); ijEntry != iRows.end(); ijEntry++) {
			// delete 2d matrix for current ij
			INTARNA_CLEANUP(*ijEntry);
		}
	}
	// clear matrix, free data
	hybridE_seed.clear();

	// delete 3rd and 4th dimension of the matrix
	for (E4dMatrix::iterator1 iRows = hybridE_multi.begin1(); iRows != hybridE_multi.end1(); iRows++) {
		for (E4dMatrix::iterator2 ijEntry = iRows.begin(); ijEntry != iRows.end(); ijEntry++) {
			// delete 2d matrix for current ij
			INTARNA_CLEANUP(*ijEntry);
		}
	}
	// clear matrix, free data
	hybridE_multi.clear();

	// clean up super class data structures
	PredictorMfe4d::clear();
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe4dMultiSeed::
fillHybridE_seed( )
{

	// global vars to avoid reallocation
	size_t i1,i2,j1,j2,w1,w2,k1,k2;

	//////////  FIRST ROUND : COMPUTE HYBRIDIZATION ENERGIES ONLY  ////////////

	// current minimal value
	E_type curMinE = E_INF, curMinE_multi = E_INF;
	// iterate increasingly over all window sizes w1 (seq1) and w2 (seq2)
	// minimal size == number of seed base pairs
	for (w1=0; w1<energy.getAccessibility1().getMaxLength(); w1++) {
	for (w2=0; w2<energy.getAccessibility2().getMaxLength(); w2++) {
		// iterate over all window starts i1 (seq1) and i2 (seq2)
		// TODO PARALLELIZE THIS DOUBLE LOOP ?!
		for (i1=0; i1+w1<hybridE_seed.size1(); i1++) {
		for (i2=0; i2+w2<hybridE_seed.size2(); i2++) {
			// check if left boundary is complementary
			if (hybridE_seed(i1,i2) == NULL) {
				// interaction not possible: nothing to do, since no storage reserved
				continue;
			}

			// check if widths are possible (ie available within data structure)
			if ( w1 >= hybridE_seed(i1,i2)->size1() || w2 >= hybridE_seed(i1,i2)->size2() ) {
				// interaction not possible: nothing to do, since no storage reserved
				continue;
			}
			assert(hybridE(i1,i2)->size1() > w1);
			assert(hybridE(i1,i2)->size2() > w2);
			// check if no interaction without seed possible -> none with neither
			if ( E_isINF((*hybridE(i1,i2))(w1,w2))
					// check if seed base pairs fit not into interaction window
					|| w1+1 < seedHandler.getConstraint().getBasePairs()
					|| w2+1 < seedHandler.getConstraint().getBasePairs() )
			{
				// ignore this entry
				(*hybridE_seed(i1,i2))(w1,w2) = E_INF;
				(*hybridE_multi(i1,i2))(w1,w2) = E_INF;
				continue;
			}

			// get window ends j1 (seq1) and j2 (seq2)
			j1=i1+w1;
			j2=i2+w2;

			// compute entry hybridE_seed ///////////////////////////
			curMinE = E_INF;
			curMinE_multi = E_INF;

			// base case = incorporate mfe seed starting at (i1,i2)
			//             + interaction on right side up to (j1,j2)
			// -> check if widths are wide enough to contain a seed
			if ( E_isNotINF( seedHandler.getSeedE(i1,i2) ) ) {
				// decode right mfe boundary
				k1 = i1+seedHandler.getSeedLength1(i1,i2)-1;
				k2 = i2+seedHandler.getSeedLength2(i1,i2)-1;
				// compute overall energy of seed+singleSide
				if ( k1 <= j1 && k2 <= j2
						&& hybridE(k1,k2) != NULL
						&& j1-k1 < hybridE(k1,k2)->size1()
						&& j2-k2 < hybridE(k1,k2)->size2()
						&& E_isNotINF((*hybridE(k1,k2))(j1-k1,j2-k2)))
				{
					// check extension of single-side interaction
					curMinE = seedHandler.getSeedE(i1,i2) + (*hybridE(k1,k2))(j1-k1,j2-k2);
				}
				// compute overall energy of seed+multiSide
				// TODO shrink k1<j1 checks to skip more entries
				if ( k1 <= j1 && k2 <= j2
						&& hybridE_multi(k1,k2) != NULL
						&& j1-k1 < hybridE_multi(k1,k2)->size1()
						&& j2-k2 < hybridE_multi(k1,k2)->size2()
						&& E_isNotINF((*hybridE_multi(k1,k2))(j1-k1,j2-k2)))
				{
					// check extension of multi-side interaction
					curMinE = std::min( curMinE, seedHandler.getSeedE(i1,i2) + (*hybridE_multi(k1,k2))(j1-k1,j2-k2) );
				}
			}

			// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
			// where k1..j1 contains a seed
			for (k1=std::min(j1-seedHandler.getConstraint().getBasePairs()+1,i1+energy.getMaxInternalLoopSize1()+1); k1>i1; k1--) {
			for (k2=std::min(j2-seedHandler.getConstraint().getBasePairs()+1,i2+energy.getMaxInternalLoopSize2()+1); k2>i2; k2--) {
				// check if (k1,k2) are valid left boundaries including a seed
				if ( hybridE_seed(k1,k2) != NULL
						&& j1-k1 < hybridE_seed(k1,k2)->size1()
						&& j2-k2 < hybridE_seed(k1,k2)->size2()
						&& E_isNotINF( (*hybridE_seed(k1,k2))(j1-k1,j2-k2) ) )
				{
					curMinE = std::min( curMinE,
							(energy.getE_interLeft(i1,k1,i2,k2)
									+ (*hybridE_seed(k1,k2))(j1-k1,j2-k2) )
						);
				}
				// check if (k1,k2) are valid left boundaries including a multi-side
				if ( hybridE_multi(k1,k2) != NULL
						&& j1-k1 < hybridE_multi(k1,k2)->size1()
						&& j2-k2 < hybridE_multi(k1,k2)->size2()
						&& E_isNotINF( (*hybridE_multi(k1,k2))(j1-k1,j2-k2) ) )
				{
					curMinE_multi = std::min( curMinE_multi,
							(energy.getE_interLeft(i1,k1,i2,k2)
									+ (*hybridE_multi(k1,k2))(j1-k1,j2-k2) )
						);
				}
			}
			}

			///////////////// compute multi-side cases /////////////////

            // Both-sided structure seeded interaction on the right
            if (allowES == ES_both) {
                for (k1 = j1; k1 > i1 + InteractionEnergy::minDistES; k1--) {
                    for (k2 = j2; k2 > i2 + InteractionEnergy::minDistES; k2--) {
                        if (hybridE_seed(k1, k2) != NULL
                        		&& hybridE_seed(k1, k2)->size1() > (j1 - k1)
                        		&& hybridE_seed(k1, k2)->size2() > (j2 - k2))
                        {
							// update minE
                            curMinE_multi = std::min(curMinE_multi,
                                               (energy.getE_multi(i1, k1, i2, k2, InteractionEnergy::ES_multi_mode::ES_multi_both)
                                                + (*hybridE_seed(k1, k2))(j1 - k1, j2 - k2)
                                               ));
                        }
                    }
                }
            }

            // Structure in S1 seeded interaction on the right
            if (allowES == ES_target || allowES == ES_xorQueryTarget) {
                for (k1 = j1; k1 > i1 + InteractionEnergy::minDistES; k1--) {
                    for (k2 = std::min(j2, i2 + energy.getMaxInternalLoopSize2() + 1);
                         k2 > i2; k2--) {
                        if (hybridE_seed(k1, k2) != NULL
                        		&& hybridE_seed(k1, k2)->size1() > (j1 - k1)
                        		&& hybridE_seed(k1, k2)->size2() > (j2 - k2))
                        {
							// update minE
                            curMinE_multi = std::min(curMinE_multi,
                                               (energy.getE_multi(i1, k1, i2, k2, InteractionEnergy::ES_multi_mode::ES_multi_1only)
                                                + (*hybridE_seed(k1, k2))(j1 - k1, j2 - k2)));
                        }
                    }
                }
            }


            // Structure in S2 with seeded interaction on the right
            if (allowES == ES_query || allowES == ES_xorQueryTarget) {
                for (k1 = std::min(j1, i1 + energy.getMaxInternalLoopSize1() + 1); k1 > i1; k1--) {
                    for (k2 = j2; k2 > i2 + InteractionEnergy::minDistES; k2--) {
                        if (hybridE_seed(k1, k2) != NULL
                        		&& hybridE_seed(k1, k2)->size1() > (j1 - k1)
                        		&& hybridE_seed(k1, k2)->size2() > (j2 - k2))
                        {
							// update minE
                            curMinE_multi = std::min(curMinE_multi,
                                               (energy.getE_multi(i1, k1, i2, k2, InteractionEnergy::ES_multi_mode::ES_multi_2only)
                                                + (*hybridE_seed(k1, k2))(j1 - k1, j2 - k2)
                                               ));
                        }
                    }
                }
            }


			// store computed value
			(*hybridE_seed(i1,i2))(w1,w2) = curMinE;
			(*hybridE_multi(i1,i2))(w1,w2) = curMinE_multi;

			// update mfe if needed (call super class)
			if (E_isNotINF(curMinE)) {
				updateOptima( i1,j1,i2,j2, curMinE, true );
			}

		}
		}
	}
	}

}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfe4dMultiSeed::
traceBack( Interaction & interaction )
{
	// check if something to trace
	if (interaction.basePairs.size() < 2) {
		return;
	}

#if IN_DEBUG_MODE
	// sanity checks
	if ( ! interaction.isValid() ) {
		throw std::runtime_error("PredictorMfe4dMultiSeed::traceBack() : given interaction not valid");
	}
	if ( interaction.basePairs.size() != 2 ) {
		throw std::runtime_error("PredictorMfe4dMultiSeed::traceBack() : given interaction does not contain boundaries only");
	}
#endif


	// check for single interaction (start==end)
	if (interaction.basePairs.begin()->first == interaction.basePairs.rbegin()->first) {
		// delete second boundary (identical to first)
		interaction.basePairs.resize(1);
		// update done
		return;
	}

	// ensure sorting
	interaction.sort();
	// get indices in hybridE_seed for boundary base pairs
	size_t	i1 = energy.getIndex1(interaction.basePairs.at(0)),
			j1 = energy.getIndex1(interaction.basePairs.at(1)),
			i2 = energy.getIndex2(interaction.basePairs.at(0)),
			j2 = energy.getIndex2(interaction.basePairs.at(1)),
			k1, k2;

	// the currently traced value for i1-j1, i2-j2
	E_type curE = (*hybridE_seed(i1,i2))(j1-i1,j2-i2);

	// trace back
	// do until only right boundary is left over
	bool doLocalTrace = true;
	bool traceInESeed = true; // if false: traces in hybridE_multi
	while( i1 != j1 ) {

		// check if we still have to find the seed
		if (doLocalTrace) {

			// check seed extended to the right
			if ( traceInESeed && E_isNotINF( seedHandler.getSeedE(i1,i2) ) ) {

				// right boundary of seed
				k1 = i1 + seedHandler.getSeedLength1(i1,i2) -1;
				k2 = i2 + seedHandler.getSeedLength2(i1,i2) -1;

				// check if correct trace with single-side extension
				if ( hybridE_seed(k1,k2) != NULL
						&& E_equal( curE, seedHandler.getSeedE(i1,i2) + (*hybridE(k1,k2))(j1-k1,j2-k2) ) )
				{
					// store seed information
					interaction.setSeedRange(
									energy.getBasePair(i1,i2),
									energy.getBasePair(k1,k2),
									energy.getE(i1,k1,i2,k2,seedHandler.getSeedE(i1,i2))+energy.getE_init());
					// trace back seed base pairs
					seedHandler.traceBackSeed( interaction, i1, i2 );
					// continue after seed
					i1 = k1;
					i2 = k2;
					curE = (*hybridE(k1,k2))(j1-k1,j2-k2);
					doLocalTrace = false;
					continue;
				}

				// check if correct trace with multi-side extension
				if ( hybridE_multi(k1,k2) != NULL
						&& E_equal( curE, seedHandler.getSeedE(i1,i2) + (*hybridE_multi(k1,k2))(j1-k1,j2-k2) ) )
				{
					// store seed information
					interaction.setSeedRange(
									energy.getBasePair(i1,i2),
									energy.getBasePair(k1,k2),
									energy.getE(i1,k1,i2,k2,seedHandler.getSeedE(i1,i2))+energy.getE_init());
					// trace back seed base pairs
					seedHandler.traceBackSeed( interaction, i1, i2 );
					// add right most base pair of seed
					interaction.basePairs.push_back( energy.getBasePair(k1,k2) );
					// move seed to gap information
					if (interaction.gap == NULL) { interaction.gap = new Interaction::Gap(); }
					interaction.gap->seeds.push_back( *(interaction.seed) );
					// continue after seed
					i1 = k1;
					i2 = k2;
					curE = (*hybridE_multi(k1,k2))(j1-k1,j2-k2);
					traceInESeed = false;
					continue;
				}
			}
			// check all interval splits if no trace already found
			if ( doLocalTrace )
			{
				// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
				// where k1..j1 contains a seed
				bool traceNotFound = true;
				for (k1=std::min(j1-seedHandler.getConstraint().getBasePairs()+1,i1+energy.getMaxInternalLoopSize1()+1); traceNotFound && k1>i1; k1--) {
				for (k2=std::min(j2-seedHandler.getConstraint().getBasePairs()+1,i2+energy.getMaxInternalLoopSize2()+1); traceNotFound && k2>i2; k2--) {
					// check if (k1,k2) are valid left boundaries including a seed
					if ( 	traceInESeed
							&& hybridE_seed(k1,k2) != NULL
							&& j1-k1 < hybridE_seed(k1,k2)->size1()
							&& j2-k2 < hybridE_seed(k1,k2)->size2()
							&& E_isNotINF( (*hybridE_seed(k1,k2))(j1-k1,j2-k2) ) )
					{
						// check if correct split
						if (E_equal ( curE,
								(energy.getE_interLeft(i1,k1,i2,k2)
										+ (*hybridE_seed(k1,k2))(j1-k1,j2-k2) )
								) )
						{
							// update trace back boundary
							i1=k1;
							i2=k2;
							curE= (*hybridE_seed(k1,k2))(j1-k1,j2-k2);
							// stop search splits
							traceNotFound = false;
							// store splitting base pair
							interaction.basePairs.push_back( energy.getBasePair(k1,k2) );
							continue;
						}
					} else
					// check if (k1,k2) are valid left boundaries including a seed
					if ( 	!traceInESeed
							&& hybridE_multi(k1,k2) != NULL
							&& j1-k1 < hybridE_multi(k1,k2)->size1()
							&& j2-k2 < hybridE_multi(k1,k2)->size2()
							&& E_isNotINF( (*hybridE_multi(k1,k2))(j1-k1,j2-k2) ) )
					{
						// check if correct split
						if (E_equal ( curE,
								(energy.getE_interLeft(i1,k1,i2,k2)
										+ (*hybridE_multi(k1,k2))(j1-k1,j2-k2) )
								) )
						{
							// update trace back boundary
							i1=k1;
							i2=k2;
							curE= (*hybridE_multi(k1,k2))(j1-k1,j2-k2);
							// stop search splits
							traceNotFound = false;
							// store splitting base pair
							interaction.basePairs.push_back( energy.getBasePair(k1,k2) );
							continue;
						}
					}
				} // k2
				} // k1

				///////////////  multi-side trace  ///////////////////

				// Structure in both
				if (!traceInESeed && traceNotFound && allowES == ES_both) {
					for (k1 = j1; traceNotFound && k1 > i1 + InteractionEnergy::minDistES; k1--) {
						for (k2 = j2; traceNotFound && k2 > i2 + InteractionEnergy::minDistES; k2--) {
							if (hybridE_seed(k1, k2) != NULL
									&& hybridE_seed(k1, k2)->size1() > (j1 - k1)
									&& hybridE_seed(k1, k2)->size2() > (j2 - k2))
							{
								if (E_equal(curE,
											(energy.getE_multi(i1, k1, i2, k2, InteractionEnergy::ES_multi_mode::ES_multi_both)
											 + (*hybridE_seed(k1, k2))(j1 - k1, j2 - k2)
											))) {
									// stop searching
									traceNotFound = false;
									traceInESeed = true;
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
									curE = (*hybridE_seed(i1, i2))(j1 - i1, j2 - i2);
									continue;
								}
							}
						}
					}
				}

				// Structure in S1
				if (!traceInESeed && traceNotFound && (allowES == ES_target || allowES == ES_xorQueryTarget)) {
					for (k1 = j1; traceNotFound && k1 > i1 + InteractionEnergy::minDistES; k1--) {
						for (k2 = std::min(j2, i2 + energy.getMaxInternalLoopSize2() + 1); traceNotFound && k2 > i2; k2--) {
							if (hybridE_seed(k1, k2) != NULL
									&& hybridE_seed(k1, k2)->size1() > (j1 - k1)
									&& hybridE_seed(k1, k2)->size2() > (j2 - k2))
							{
								if (E_equal(curE,
											(energy.getE_multi(i1, k1, i2, k2, InteractionEnergy::ES_multi_mode::ES_multi_1only)
											 + (*hybridE_seed(k1, k2))(j1 - k1, j2 - k2)
											))) {
									// stop searching
									traceNotFound = false;
									traceInESeed = true;
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
									curE = (*hybridE_seed(i1, i2))(j1 - i1, j2 - i2);
									continue;
								}
							}
						}
					}
				}

				// Structure in S2
				if (!traceInESeed && traceNotFound && (allowES == ES_query || allowES == ES_xorQueryTarget)) {
					for (k1 = std::min(j1, i1 + energy.getMaxInternalLoopSize1() + 1); traceNotFound && k1 > i1; k1--) {
						for (k2 = j2; traceNotFound && k2 > i2 + InteractionEnergy::minDistES; k2--) {
							if (hybridE_seed(k1, k2) != NULL
									&& hybridE_seed(k1, k2)->size1() > (j1 - k1)
									&& hybridE_seed(k1, k2)->size2() > (j2 - k2))
							{
								if (E_equal(curE,
											(energy.getE_multi(i1, k1, i2, k2, InteractionEnergy::ES_multi_mode::ES_multi_2only)
											 + (*hybridE_seed(k1, k2))(j1 - k1, j2 - k2)
											))) {
									// stop searching
									traceNotFound = false;
									traceInESeed = true;
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
									curE = (*hybridE_seed(i1, i2))(j1 - i1, j2 - i2);
									continue;
								}
							}
						}
					}
				}


				// final sanity check
				assert(!traceNotFound);

			} // local trace
		}
		// seed was already traced, do "normal" interaction trace
		else {
			// create temporary data structure to be filed
			Interaction rightSide( *interaction.s1, *interaction.s2 );
			rightSide.basePairs.push_back( energy.getBasePair(i1,i2) );
			rightSide.basePairs.push_back( energy.getBasePair(j1,j2) );
			// call traceback of super class
			PredictorMfe4d::traceBack( rightSide );
			// copy base pairs (excluding last)
			for (size_t i=0; i+1<rightSide.basePairs.size(); i++) {
				interaction.basePairs.push_back( rightSide.basePairs.at(i) );
			}
			i1 = j1;
			i2 = j2;
			// stop traceback
			break;
		}
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

} // namespace

////////////////////////////////////////////////////////////////////////////

