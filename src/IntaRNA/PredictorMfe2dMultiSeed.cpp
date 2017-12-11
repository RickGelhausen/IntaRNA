
#include "PredictorMfe2dMultiSeed.h"

namespace IntaRNA {

//////////////////////////////////////////////////////////////////////////

PredictorMfe2dMultiSeed::
PredictorMfe2dMultiSeed(
		const InteractionEnergy & energy
		, OutputHandler & output
		, PredictionTracker * predTracker
		, const AllowES allowES_
		, SeedHandler * seedHandlerInstance
)
		: PredictorMfe2dMulti(energy,output,predTracker, allowES_)
		, hybridE_pq_seed(0,0)
		, seedHandler(seedHandlerInstance)
{
	assert( seedHandler.getConstraint().getBasePairs() > 1 );
}

//////////////////////////////////////////////////////////////////////////

PredictorMfe2dMultiSeed::
~PredictorMfe2dMultiSeed()
{
}

//////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dMultiSeed::
predict( const IndexRange & r1
		, const IndexRange & r2
		, const OutputConstraint & outConstraint )
{
#if INTARNA_MULITHREADING
#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) <<"predicting mfe multi-side interactions with seed in O(n^2) space and O(n^5) time..."; }
	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

	// suboptimal setup check
	if (outConstraint.reportMax>1 && outConstraint.reportOverlap != OutputConstraint::ReportOverlap::OVERLAP_BOTH) {
		throw std::runtime_error("PredictorMfe2dMultiSeed : the enumeration of non-overlapping suboptimal interactions is not supported in this prediction mode");
	}

#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!(r1.isAscending() && r2.isAscending()) )
		throw std::runtime_error("PredictorMfe2dMultiSeed::predict("+toString(r1)+","+toString(r2)+") is not sane");
#endif

	// setup index offset
	energy.setOffset1(r1.from);
	energy.setOffset2(r2.from);
	seedHandler.setOffset1(r1.from);
	seedHandler.setOffset2(r2.from);

	const size_t hybridE_pqsize1 = std::min( energy.size1()
			, (r1.to==RnaSequence::lastPos?energy.size1()-1:r1.to)-r1.from+1 );
	const size_t hybridE_pqsize2 = std::min( energy.size2()
			, (r2.to==RnaSequence::lastPos?energy.size2()-1:r2.to)-r2.from+1 );


	// compute seed interactions for whole range
	// and check if any seed possible
	if (seedHandler.fillSeed( 0, hybridE_pqsize1-1, 0, hybridE_pqsize2-1 ) == 0) {
		// trigger empty interaction reporting
		initOptima(outConstraint);
		reportOptima(outConstraint);
		// stop computation
		return;
	}

	// resize matrix
	hybridE_pq.resize( hybridE_pqsize1, hybridE_pqsize2 );
	hybridE_pq_seed.resize( hybridE_pqsize1, hybridE_pqsize2 );
	hybridO.resize( hybridE_pqsize1, hybridE_pqsize2 );

	// initialize mfe interaction for updates
	initOptima( outConstraint );

	// for all right ends j1
	for (size_t j1 = hybridE_pq.size1(); j1-- > 0; ) {
		// check if j1 is accessible
		if (!energy.isAccessible1(j1))
			continue;
		// iterate over all right ends j2
		for (size_t j2 = hybridE_pq.size2(); j2-- > 0; ) {
			// check if j2 is accessible
			if (!energy.isAccessible2(j2))
				continue;
			// check if base pair (j1,j2) possible
			if (!energy.areComplementary( j1, j2 ))
				continue;

			// compute hybridE_pq_seed and update mfe via PredictorMfe2d::updateOptima()
			fillHybridE_seed( j1, j2, 0, 0, outConstraint );
		}
	}

	// report mfe interaction
	reportOptima( outConstraint );

}

//////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dMultiSeed::
fillHybridE_seed( const size_t j1, const size_t j2, const size_t i1min, const size_t i2min
				, const OutputConstraint & outConstraint  )
{

	// compute hybridE_pq
	fillHybridE( j1, j2, outConstraint, i1min, i2min );

	assert(i1min <= j1);
	assert(i2min <= j2);
	assert(hybridErange.r1.from <= i1min);
	assert(hybridErange.r2.from <= i2min);
	assert(j1==hybridErange.r1.to);
	assert(j2==hybridErange.r2.to);
	assert(j1<hybridE_pq.size1());
	assert(j2<hybridE_pq.size2());

	// check if it is possible to have a seed ending on the right at (j1,j2)
	if (std::min(j1-i1min,j2-i2min)+1 < seedHandler.getConstraint().getBasePairs()) {
		// no seed possible, abort table computation
		return;
	}

	// global vars to avoid reallocation
	size_t i1,i2,k1,k2,w1,w2;

	// get i1/i2 index boundaries for computation
	const IndexRange i1range( std::max(hybridErange.r1.from,i1min), j1+1-seedHandler.getConstraint().getBasePairs() );
	const IndexRange i2range( std::max(hybridErange.r2.from,i2min), j2+1-seedHandler.getConstraint().getBasePairs() );

	//////////  COMPUTE HYBRIDIZATION ENERGIES (WITH SEED)  ////////////

	// current minimal value
	E_type curMinE = E_INF;
	E_type curMinE_seed = E_INF;
	E_type curMinO = E_INF;
	// iterate over all window starts i1 (seq1) and i2 (seq2)
	// TODO PARALLELIZE THIS DOUBLE LOOP ?!

	for (i1=hybridErange.r1.to+1; i1-- > hybridErange.r1.from; ) {

		// screen for left boundaries i2 in seq2
		for (i2=hybridErange.r2.to+1; i2-- > hybridErange.r2.from; ) {


			// check if this cell is to be computed (!=E_INF)
			if( E_isNotINF( hybridE_pq(i1,i2) ) ) {
				///////////////////////////////////////////////////
				// hybridO(i1,i2) computation
				///////////////////////////////////////////////////

				if (allowES != ES_target)
					// compute entry, since (i1,i2) complementary
				{
					// init
					curMinO = E_INF;


					for (k2 = j2; k2 > i2 + InteractionEnergy::minDistES; k2--) {
						if (E_isNotINF(hybridE_pq_seed(i1, k2))) {
							curMinO = std::min(curMinO,
											   energy.getE_multiRight(i1, i2, k2)
											   + hybridE_pq_seed(i1, k2));
						}
					}

					hybridO(i1, i2) = curMinO;
				}

				// compute entry
				curMinE = hybridE_pq(i1, i2);
				curMinE_seed = E_INF;

				// base case = incorporate mfe seed starting at (i1,i2)
				//             + interaction on right side up to (p,q)
				if (E_isNotINF(seedHandler.getSeedE(i1, i2))) {
					// decode right mfe boundary
					k1 = i1 + seedHandler.getSeedLength1(i1, i2) - 1;
					k2 = i2 + seedHandler.getSeedLength2(i1, i2) - 1;
					// compute overall energy of seed+singleSide
					if (k1 <= j1 && k2 <= j2 && E_isNotINF(hybridE_pq(k1, k2))) {
						curMinE_seed = seedHandler.getSeedE(i1, i2) + hybridE_pq(k1, k2);
					}
				}

				// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
				// where k1..j1 contains a seed
				for (k1 = std::min(j1 + 1 - seedHandler.getConstraint().getBasePairs(),
								   i1 + energy.getMaxInternalLoopSize1() + 1); k1 > i1; k1--) {
					for (k2 = std::min(j2 + 1 - seedHandler.getConstraint().getBasePairs(),
									   i2 + energy.getMaxInternalLoopSize2() + 1); k2 > i2; k2--) {

						///////////////////////////////////////////////////////////////////
						// hybridE(i1,i2) computation
						///////////////////////////////////////////////////////////////////

						if (E_isNotINF(hybridE_pq(k1, k2))) {
							curMinE = std::min(curMinE,
											   (energy.getE_interLeft(i1, k1, i2, k2)
												+ hybridE_pq(k1, k2))
							);

						}

						///////////////////////////////////////////////////////////////////
						// hybridE_seed(i1,i2) computation
						///////////////////////////////////////////////////////////////////

						if (E_isNotINF(hybridE_pq_seed(k1, k2))) {
							curMinE_seed = std::min(curMinE_seed,
													(energy.getE_interLeft(i1, k1, i2, k2)
													 + hybridE_pq_seed(k1, k2))
							);
						}

					} // k2
				} // k1

				////////////////////////////////////////////////////////////////
				// hybridE(i1,i2) Multiloop cases = ES-gap
				////////////////////////////////////////////////////////////////

				// Both-sided structure
				if (allowES == ES_both) {
					for (k1 = j1; k1 > i1 + InteractionEnergy::minDistES; k1--) {
						if (E_isNotINF(hybridO(k1, i2))) {
							// update minE
							curMinE = std::min(curMinE,
											   (energy.getE_multiLeft(i1, k1, i2,
																	  InteractionEnergy::ES_multi_mode::ES_multi_both)
												+ hybridO(k1, i2)
											   ));
						}
					}
				}

				// Structure in S1
				if (allowES == ES_target || allowES == ES_xorQueryTarget) {
					for (k1 = j1; k1 > i1 + InteractionEnergy::minDistES; k1--) {
						for (k2 = std::min(j2, i2 + energy.getMaxInternalLoopSize2() + 1); k2 > i2; k2--) {
							if (E_isNotINF(hybridE_pq_seed(k1, k2))) {
								// update minE
								curMinE = std::min(curMinE,
												   (energy.getE_multi(i1, k1, i2, k2,
																	  InteractionEnergy::ES_multi_mode::ES_multi_1only)
													+ hybridE_pq_seed(k1, k2)
												   ));
							}
						}
					}
				}

				// Structure in S2
				if (allowES == ES_query || allowES == ES_xorQueryTarget) {
					for (k1 = std::min(j2, i1 + energy.getMaxInternalLoopSize1() + 1); k1 > i1; k1--) {
						if (E_isNotINF(hybridO(k1, i2))) {
							// update minE
							curMinE = std::min(curMinE,
											   (energy.getE_multiLeft(i1, k1, i2,
																	  InteractionEnergy::ES_multi_mode::ES_multi_2only)
												+ hybridO(k1, i2)
											   ));
						}
					}
				}

				// store value
				hybridE_pq(i1, i2) = curMinE;
				hybridE_pq_seed(i1, i2) = curMinE_seed;

				// update mfe if needed (call super class)
				if (E_isNotINF(curMinE_seed)) {
					PredictorMfe2dMulti::updateOptima(i1, j1, i2, j2, curMinE_seed, true);
				}
			}
		}
	}

}

////////////////////////////////////////////////////////////////////////////

size_t
PredictorMfe2dMultiSeed::
traceHybridO( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2 ) const
{
	E_type curE = hybridO(i1,i2);

	size_t k2;
	for (k2 = j2; k2 > i2 + InteractionEnergy::minDistES; k2--) {
		if (E_isNotINF(hybridE_pq_seed(i1, k2)))
		{
			if ( E_equal(curE, energy.getE_multiRight(i1, i2, k2)
							   + hybridE_pq_seed(i1, k2)))
			{
				return k2;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dMultiSeed::
traceBack( Interaction & interaction, const OutputConstraint & outConstraint  )
{
	// check if something to trace
	if (interaction.basePairs.size() < 2) {
		return;
	}

#if INTARNA_IN_DEBUG_MODE
	// sanity checks
	if ( ! interaction.isValid() ) {
		throw std::runtime_error("PredictorMfe2dMultiSeed::traceBack() : given interaction not valid");
	}
	if ( interaction.basePairs.size() != 2 ) {
		throw std::runtime_error("PredictorMfe2dMultiSeed::traceBack() : given interaction does not contain boundaries only");
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
			j2 = energy.getIndex2(interaction.basePairs.at(1))
	;

#if INTARNA_IN_DEBUG_MODE
	// check if intervals are larger enough to contain a seed
	if (std::min(j1-i1,j2-i2)+1 < seedHandler.getConstraint().getBasePairs()) {
		// no seed possible, abort computation
		throw std::runtime_error("PredictorMfe2dMultiSeed::traceBack() : given boundaries "+toString(interaction)+" can not hold a seed of "+toString(seedHandler.getConstraint().getBasePairs())+" base pairs");
	}
#endif


	// refill submatrices of mfe interaction
	fillHybridE_seed( j1, j2, i1, i2, outConstraint );

	// the currently traced value for i1-j1, i2-j2
	E_type curE = hybridE_pq_seed(i1,i2);

	// trace back
	// temp variables
	size_t k1,k2;
	bool traceInESeed = true; // if false: traces hybridE

	while( (j1-i1) > 1 ) {

		// check if we still have to find the seed
		if (traceInESeed) {

			// check base case == seed only
			if (E_isNotINF(seedHandler.getSeedE(i1, i2))) {

				// right boundary of seed
				k1 = i1 + seedHandler.getSeedLength1(i1, i2) - 1;
				k2 = i2 + seedHandler.getSeedLength2(i1, i2) - 1;

				// check if correct trace single-side extension
				if (E_isNotINF(hybridE_pq(k1, k2))
					&& E_equal(curE, seedHandler.getSeedE(i1, i2) + hybridE_pq(k1, k2))) {
					// store seed information
					interaction.setSeedRange(
							energy.getBasePair(i1, i2),
							energy.getBasePair(k1, k2),
							energy.getE(i1, k1, i2, k2, seedHandler.getSeedE(i1, i2)) + energy.getE_init());
					// trace back seed base pairs
					seedHandler.traceBackSeed(interaction, i1, i2);
					// continue after seed
					i1 = k1;
					i2 = k2;
					curE = hybridE_pq(k1, k2);
					traceInESeed = false;
					continue;
				}
			} // seed extension

			// check interior loop extension of hybridE_seed()
			// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
			bool traceNotFound = true;
			for (k1 = std::min(j1 - seedHandler.getConstraint().getBasePairs() + 1,
							   i1 + energy.getMaxInternalLoopSize1() + 1); traceNotFound && k1 > i1; k1--) {
			for (k2 = std::min(j2 - seedHandler.getConstraint().getBasePairs() + 1,
							   i2 + energy.getMaxInternalLoopSize2() + 1); traceNotFound && k2 > i2; k2--) {
				// check if (k1,k2) are valid left boundaries including a seed
				// check if correct split
				if (E_isNotINF(hybridE_pq_seed(k1, k2))
					&& E_equal (curE,
							 (energy.getE_interLeft(i1, k1, i2, k2)
							  + hybridE_pq_seed(k1, k2))
				)) {
					// stop searching
					traceNotFound = false;
					// store splitting base pair
					interaction.basePairs.push_back(energy.getBasePair(k1, k2));
					// update trace back boundary
					i1 = k1;
					i2 = k2;
					curE = hybridE_pq_seed(k1, k2);
				}
			} // k2
			} // k1
		} // traceInESeed

		// trace in hybridE()
		else {
			// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
			bool traceNotFound = true;
			for (k1 = std::min(j1, i1 + energy.getMaxInternalLoopSize1() + 1); traceNotFound && k1 > i1; k1--) {
			for (k2 = std::min(j2, i2 + energy.getMaxInternalLoopSize2() + 1); traceNotFound && k2 > i2; k2--) {
				// check if (k1,k2) are valid left boundaries including a seed
				// check if correct split
				if (E_isNotINF(hybridE_pq(k1, k2))
					&& E_equal (curE,
							 (energy.getE_interLeft(i1, k1, i2, k2)
							  + hybridE_pq(k1, k2))
				)) {
					// stop searching
					traceNotFound = false;
					// store splitting base pair
					interaction.basePairs.push_back(energy.getBasePair(k1, k2));
					// update trace back boundary
					i1 = k1;
					i2 = k2;
					curE = hybridE_pq(k1, k2);
					continue;
				}
			} // k2
			} // k1

			///////////////  multi-side trace  ///////////////////

			// Both-sided structure
			if (traceNotFound && allowES == ES_both) {
				for (k1 = j1; traceNotFound && k1 > i1 + InteractionEnergy::minDistES; k1--) {
					if (E_equal(curE,
								(energy.getE_multiLeft(i1, k1, i2, InteractionEnergy::ES_multi_mode::ES_multi_both)
								 + hybridO(k1, i2)
								))) {
						// stop searching
						traceNotFound = false;
						traceInESeed = true;
						// Determine k2 based on k1
						k2 = traceHybridO(k1, j1, i2, j2);
						E_type E_multiRight = energy.getE_multiRight(i1, i2, k2);
						// store splitting base pair
						interaction.basePairs.push_back(energy.getBasePair(k1, k2));
						// store gap information
						if (interaction.gap == NULL) { interaction.gap = new Interaction::Gap(); }
						interaction.gap->energy += energy.getE_multiLeft(i1, k1, i2, InteractionEnergy::ES_multi_mode::ES_multi_both) + E_multiRight;
						Interaction::BasePair bpLeft = energy.getBasePair(i1,i2);
						interaction.gap->gaps1.insert( IndexRange(bpLeft.first+1,interaction.basePairs.rbegin()->first-1) );
						interaction.gap->gaps2.insert( IndexRange(interaction.basePairs.rbegin()->second+1,bpLeft.second-1) );
						// move seed to gap information
						interaction.gap->seeds.push_back( *(interaction.seed) );
						// trace right part of split
						i1 = k1;
						i2 = k2;
						curE = hybridE_pq_seed(k1, k2);
						continue;
					}
				}
			}

			// Structure in S1
			if (traceNotFound && (allowES == ES_target || allowES == ES_xorQueryTarget)) {
				for (k1 = j1; traceNotFound && k1 > i1 + InteractionEnergy::minDistES; k1--) {
				for (k2 = std::min(j2, i2 + energy.getMaxInternalLoopSize2() + 1); traceNotFound && k2 > i2; k2--) {
					if (E_equal(curE,
								(energy.getE_multi(i1, k1, i2, k2, InteractionEnergy::ES_multi_mode::ES_multi_1only)
								 + hybridE_pq_seed(k1, k2)
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
						interaction.gap->gaps1.insert( IndexRange(bpLeft.first+1,interaction.basePairs.rbegin()->first-1) );
						// move seed to gap information
						interaction.gap->seeds.push_back( *(interaction.seed) );
						// trace right part of split
						i1 = k1;
						i2 = k2;
						curE = hybridE_pq_seed(k1, k2);
						continue;
					}
				}
				}
			}

			// Structure in S2
			if (traceNotFound && (allowES == ES_query || allowES == ES_xorQueryTarget)) {
				for (k1 = std::min(j1, i1 + energy.getMaxInternalLoopSize1() + 1); traceNotFound && k1 > i1; k1--) {
					if (E_equal(curE,
								(energy.getE_multiLeft(i1, k1, i2, InteractionEnergy::ES_multi_mode::ES_multi_2only)
								 + hybridO(k1, i2)
								))) {
						// stop searching
						traceNotFound = false;
						traceInESeed = true;
						// Determine k2 based on k1
						k2 = traceHybridO(k1, j1, i2, j2);
						E_type E_multiRight = energy.getE_multiRight(i1, i2, k2);
						// store splitting base pair
						interaction.basePairs.push_back(energy.getBasePair(k1, k2));
						// store gap information
						if (interaction.gap == NULL) { interaction.gap = new Interaction::Gap(); }
						interaction.gap->energy += energy.getE_multiLeft(i1, k1, i2, InteractionEnergy::ES_multi_mode::ES_multi_2only) + E_multiRight;
						Interaction::BasePair bpLeft = energy.getBasePair(i1,i2);
						interaction.gap->gaps2.insert( IndexRange(interaction.basePairs.rbegin()->second+1,bpLeft.second-1) );
						// move seed to gap information
						interaction.gap->seeds.push_back( *(interaction.seed) );
						// trace right part of split
						i1 = k1;
						i2 = k2;
						curE = hybridE_pq_seed(k1, k2);
						continue;
					}
				}
			}

			// final sanity check
			assert(!traceNotFound);
		} // local trace
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
PredictorMfe2dMultiSeed::
getNextBest( Interaction & curBest )
{
	throw std::runtime_error("PredictorMfe2dMultiSeed::getNextBest() : This prediction mode does not support non-overlapping suboptimal interaction enumeration.");
}

//////////////////////////////////////////////////////////////////////////




} // namespace
