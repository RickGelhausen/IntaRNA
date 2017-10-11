#include "PredictorMfe2dMultiHeuristicSeed.h"

namespace IntaRNA {

////////////////////////////////////////////////////////////////////////////


// TODO: FIX CONSTRUCTOR
PredictorMfe2dMultiHeuristicSeed::
PredictorMfe2dMultiHeuristicSeed( const InteractionEnergy & energy
		, OutputHandler & output
		, PredictionTracker * predTracker
		, const AllowES allowES_
		, const SeedConstraint & seedConstraint
)
		: PredictorMfe2dMultiHeuristic(energy,output,predTracker, allowES_)
		, seedHandler( energy, seedConstraint )
		, hybridE_seed(0,0)
{
}


////////////////////////////////////////////////////////////////////////////

PredictorMfe2dMultiHeuristicSeed::
~PredictorMfe2dMultiHeuristicSeed()
{
	// clean up
}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dMultiHeuristicSeed::
predict( const IndexRange & r1
		, const IndexRange & r2
		, const OutputConstraint & outConstraint )
{
#if INTARNA_MULITHREADING
#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) <<"predicting mfe multi-side interactions with seed heuristically in O(n^2) space and O(n^3) time..."; }
	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!(r1.isAscending() && r2.isAscending()) )
		throw std::runtime_error("PredictorMfe2dMultiHeuristicSeed::predict("+toString(r1)+","+toString(r2)+") is not sane");
#endif


	// set index offset
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
	hybridO.resize( hybridEsize1, hybridEsize2 );

	// temp vars
	size_t i1,i2,w1,w2;

	// init hybridE matrix
	bool isValidCell = true;
	for (i1=0; i1<hybridE.size1(); i1++) {
	for (i2=0; i2<hybridE.size2(); i2++) {

		// check if positions can form interaction
		if (	energy.isAccessible1(i1)
				&& energy.isAccessible2(i2)
				&& energy.areComplementary(i1,i2) )
		{
			// set to interaction initiation with according boundary
			hybridE(i1,i2) = BestInteraction(energy.getE_init(), i1, i2);
		} else {
			// set to infinity, ie not used
			hybridE(i1,i2) = BestInteraction(E_INF, RnaSequence::lastPos, RnaSequence::lastPos);
		}

		hybridE_multi(i1,i2) = BestInteraction(E_INF, RnaSequence::lastPos, RnaSequence::lastPos);
		// init seed data
		hybridE_seed(i1,i2) = BestInteraction(E_INF, RnaSequence::lastPos, RnaSequence::lastPos);
	} // i2
	} // i1

	for (i1=0; i1<hybridO.size1(); i1++) {
	for (i2=0; i2<hybridO.size2(); i2++) {
		hybridO(i1,i2) =  BestInteraction(E_INF, RnaSequence::lastPos, RnaSequence::lastPos);
	} // i2
	} // i1

	// init mfe without seed condition
	OutputConstraint tmpOutConstraint(1, outConstraint.reportOverlap, outConstraint.maxE, outConstraint.deltaE);
	initOptima( tmpOutConstraint );

	// compute hybridization energies WITHOUT seed condition
	// sets also -energy -hybridE
	// -> no hybrid update since updateOptima overwritten
	PredictorMfe2dMultiHeuristic::fillHybridE();

	// check if any interaction possible
	// if not no seed-containing interaction is possible neither
	if (!(this->mfeInteractions.begin()->energy < tmpOutConstraint.maxE)) {
		// stop computation since no favorable interaction found
		reportOptima(tmpOutConstraint);
		return;
	}


	// init mfe for later updates
	initOptima( outConstraint );

	// compute entries
	// current minimal value
	E_type curE = E_INF, curEtotal = E_INF, curCellEtotal = E_INF;
	BestInteraction * curCell = NULL;
	BestInteraction * curCellMulti = NULL;
	BestInteraction * curCellO = NULL;
	const BestInteraction * rightExt = NULL;

	// iterate (decreasingly) over all left interaction starts
	for (i1=hybridE_seed.size1(); i1-- > 0;) {
	for (i2=hybridE_seed.size2(); i2-- > 0;) {

		// check if left side can pair
		if (E_isINF(hybridE(i1,i2).E)) {
			continue;
		}
		// direct cell access
		curCell  = &(hybridE_seed(i1,i2));
		curCellMulti = &(hybridE_multi(i1,i2));
		curCellO = &(hybridO(i1, i2));

		// reset temporary variables
		curCellEtotal = E_INF;

		// HybridO computation
		if (allowES != ES_target) {
			for (w2 = 1; i2 + InteractionEnergy::minDistES + w2 < hybridE_seed.size2(); w2++) {
				// direct cell access (const)
				rightExt = &(hybridE_seed(i1, i2 + w2));
				// check if right side can pair
				if (E_isINF(rightExt->E)) {
					continue;
				}
				// check if interaction length is within boundary
				if ((rightExt->j1 + 1 - i1) > energy.getAccessibility1().getMaxLength()
					|| (rightExt->j2 + 1 - i2) > energy.getAccessibility2().getMaxLength()) {
					continue;
				}
//				LOG(DEBUG) << "HybridO computation \n"
//						   << "i1: " << i1 << "\n"
//						   << "j1: " << rightExt->j1 << "\n"
//						   << "i2: " << i2 << "\n"
//						   << "j1: " << rightExt->j2 << "\n";
				// compute energy for this loop sizes
				curE = energy.getE_multiRight( i1, i2, i2 + w2) + rightExt->E;
				// check if this combination yields better energy
				curEtotal = energy.getE(i1, rightExt->j1, i2, rightExt->j2, curE);
				if (curEtotal < curCellO->E) {
					// update current best for this left boundary
					// copy right boundary
					*curCellO = *rightExt;
					// set new energy
					curCellO->E = curE;
				}
			} // w2
		}

		///////////////////////////////////////////////////////////////////
		// check all extensions of interactions CONTAINING a seed already
		///////////////////////////////////////////////////////////////////

		// TODO PARALLELIZE THIS DOUBLE LOOP ?!
		// iterate over all loop sizes w1 (seq1) and w2 (seq2)
		for (w1=1; w1-1 <= energy.getMaxInternalLoopSize1(); w1++) {
		for (w2=1; w2-1 <= energy.getMaxInternalLoopSize2(); w2++) {
			//(&(hybridE_seed(i1, i2)) != NULL)
			if (  (i1+w1 < hybridE_seed.size1())
			    && (i2+w2 < hybridE_seed.size2())) {
				// direct cell access to right side end of loop (seed has to be to the right of it)
				rightExt = &(hybridE_seed(i1 + w1, i2 + w2));
				// check if right side of loop can pair
				if (E_isINF(rightExt->E)) {
					continue;
				}
				// check if interaction length is within boundary
				if ((rightExt->j1 + 1 - i1) > energy.getAccessibility1().getMaxLength()
					|| (rightExt->j2 + 1 - i2) > energy.getAccessibility2().getMaxLength()) {
					continue;
				}
//				LOG(DEBUG) << "IL case 1 computation \n"
//						   << "i1: " << i1 << "\n"
//						   << "j1: " << rightExt->j1 << "\n"
//						   << "i2: " << i2 << "\n"
//						   << "j1: " << rightExt->j2 << "\n";
				// compute energy for this loop sizes
				curE = energy.getE_interLeft(i1, i1 + w1, i2, i2 + w2) + rightExt->E;
				// check if this combination yields better energy
				curEtotal = energy.getE(i1, rightExt->j1, i2, rightExt->j2, curE);
				if (curEtotal < curCellEtotal) {
					// update current best for this left boundary
					// copy right boundary
					*curCell = *rightExt;
					// set new energy
					curCell->E = curE;
					// store overall energy
					curCellEtotal = curEtotal;
				}
			}
			// &(hybridE_multi(i1, i2)) != NULL)
			if (  (i1+w1 < hybridE_multi.size1())
				 && (i2+w2 < hybridE_multi.size2())) {

				// direct cell access to right side end of loop (seed has to be to the right of it)
				rightExt = &(hybridE_multi(i1 + w1, i2 + w2));
				// check if right side of loop can pair
				if (E_isINF(rightExt->E)) {
					continue;
				}
				// check if interaction length is within boundary
				if ((rightExt->j1 + 1 - i1) > energy.getAccessibility1().getMaxLength()
					|| (rightExt->j2 + 1 - i2) > energy.getAccessibility2().getMaxLength()) {
					continue;
				}
//				LOG(DEBUG) << "IL case 2 computation \n"
//						   << "i1: " << i1 << "\n"
//						   << "j1: " << rightExt->j1 << "\n"
//						   << "i2: " << i2 << "\n"
//						   << "j1: " << rightExt->j2 << "\n";
				// compute energy for this loop sizes
				curE = energy.getE_interLeft(i1, i1 + w1, i2, i2 + w2) + rightExt->E;
				// check if this combination yields better energy
				curEtotal = energy.getE(i1, rightExt->j1, i2, rightExt->j2, curE);
				if (curEtotal < curCellEtotal) {
					// update current best for this left boundary
					// copy right boundary
					*curCellMulti = *rightExt;
					// set new energy
					curCellMulti->E = curE;
					// store overall energy
					curCellEtotal = curEtotal;
				}
			}
		} // w2
		} // w1

		// Multiloop cases = ES-gap

		// Both-sided structure
		if (allowES == ES_both) {
			for (w1 = 1; i1 + InteractionEnergy::minDistES + w1 < hybridE_seed.size1(); w1++) {
				// direct cell access (const)
				rightExt = &(hybridO(i1 + w1, i2));
				// check if right side can pair
				if (E_isINF(rightExt->E)) {
					continue;
				}
				// check if interaction length is within boundary
				if ((rightExt->j1 + 1 - i1) > energy.getAccessibility1().getMaxLength()
					|| (rightExt->j2 + 1 - i2) > energy.getAccessibility2().getMaxLength()) {
					continue;
				}
//				LOG(DEBUG) << "Both computation \n"
//						   << "i1: " << i1 << "\n"
//						   << "j1: " << rightExt->j1 << "\n"
//						   << "i2: " << i2 << "\n"
//						   << "j1: " << rightExt->j2 << "\n";
				// compute energy for this loop sizes
				curE = energy.getE_multiLeft(i1, i1 + w1, i2,
											 InteractionEnergy::ES_multi_mode::ES_multi_both) + rightExt->E;
				// check if this combination yields better energy
				curEtotal = energy.getE(i1, rightExt->j1, i2, rightExt->j2, curE);
				if (curEtotal < curCellEtotal) {
					// update current best for this left boundary
					// copy right boundary
					*curCellMulti = *rightExt;
					// set new energy
					curCellMulti->E = curE;
					// store total energy to avoid recomputation
					curCellEtotal = curEtotal;
				}

			} // w1
		}

		// Structure in S1
		if (allowES == ES_target || allowES == ES_xorQueryTarget) {
			for (w1 = 1; i1 + InteractionEnergy::minDistES + w1 < hybridE_seed.size1(); w1++) {
			for (w2 = 1; w2 - 1 <= energy.getMaxInternalLoopSize2() && i2 + w2 < hybridE_seed.size2(); w2++) {
				// direct cell access (const)
				rightExt = &(hybridE_seed(i1 + w1, i2 + w2));
				// check if right side can pair
				if (E_isINF(rightExt->E)) {
					continue;
				}
				// check if interaction length is within boundary
				if ((rightExt->j1 + 1 - i1) > energy.getAccessibility1().getMaxLength()
					|| (rightExt->j2 + 1 - i2) > energy.getAccessibility2().getMaxLength()) {
					continue;
				}
//				LOG(DEBUG) << "Target computation \n"
//						   << "i1: " << i1 << "\n"
//						   << "j1: " << rightExt->j1 << "\n"
//						   << "i2: " << i2 << "\n"
//						   << "j1: " << rightExt->j2 << "\n";
				// compute energy for this loop sizes
				curE = energy.getE_multi(i1, i1 + w1, i2, i2 + w2,
										 InteractionEnergy::ES_multi_mode::ES_multi_1only) + rightExt->E;
				// check if this combination yields better energy
				curEtotal = energy.getE(i1, rightExt->j1, i2, rightExt->j2, curE);
				if (curEtotal < curCellEtotal) {
					// update current best for this left boundary
					// copy right boundary
					*curCellMulti = *rightExt;
					// set new energy
					curCellMulti->E = curE;
					// store total energy to avoid recomputation
					curCellEtotal = curEtotal;
				}

			} // w2
			} // w1
		}

		// Structure in S2
		if (allowES == ES_query || allowES == ES_xorQueryTarget) {
			for (w1 = 1; w1 - 1 <= energy.getMaxInternalLoopSize1() && i1 + w1 < hybridO.size1(); w1++) {
				// direct cell access (const)
				rightExt = &(hybridO(i1 + w1, i2));
				// check if right side can pair
				if (E_isINF(rightExt->E)) {
					continue;
				}
				// check if interaction length is within boundary
				if ((rightExt->j1 + 1 - i1) > energy.getAccessibility1().getMaxLength()
					|| (rightExt->j2 + 1 - i2) > energy.getAccessibility2().getMaxLength()) {
					continue;
				}
//				LOG(DEBUG) << "Query computation \n"
//						   << "i1: " << i1 << "\n"
//						   << "j1: " << rightExt->j1 << "\n"
//						   << "i2: " << i2 << "\n"
//						   << "j1: " << rightExt->j2 << "\n";
				// compute energy for this loop sizes
				curE = energy.getE_multiLeft(i1, i1 + w1, i2,
											 InteractionEnergy::ES_multi_mode::ES_multi_2only) + rightExt->E;
				// check if this combination yields better energy
				curEtotal = energy.getE(i1, rightExt->j1, i2, rightExt->j2, curE);
				if (curEtotal < curCellEtotal) {
					// update current best for this left boundary
					// copy right boundary
					*curCellMulti = *rightExt;
					// set new energy
					curCellMulti->E = curE;
					// store total energy to avoid recomputation
					curCellEtotal = curEtotal;
				}

			} // w1
		}

		///////////////////////////////////////////////////////////////////
		// check if seed is starting here
		///////////////////////////////////////////////////////////////////
		// LOG(DEBUG) << "SeedHandler.getSeedE(i1,i2): " << seedHandler.getSeedE(i1,i2);
		// check if seed is possible for this left boundary
		if ( E_isNotINF( seedHandler.getSeedE(i1,i2) ) ) {
			// get right extension
			w1 = seedHandler.getSeedLength1(i1,i2)-1;
			w2 = seedHandler.getSeedLength2(i1,i2)-1;
//			LOG(DEBUG) << "Entered SeedHandling! " << i1 << ", " << i2;
			// compute overall energy of seed+singleside
			if (  i1+w1 < hybridE.size1()
				&& i2+w2 < hybridE.size2()) {
//				LOG(DEBUG) << "Entered SeedHandling case single!";
				rightExt = &(hybridE(i1 + w1, i2 + w2));
				// check if right side can pair
				if (E_isINF(rightExt->E)) {
					continue;
				}
				// check if interaction length is within boundary
				if ((rightExt->j1 + 1 - i1) > energy.getAccessibility1().getMaxLength()
					|| (rightExt->j2 + 1 - i2) > energy.getAccessibility2().getMaxLength()) {
					continue;
				}
//				LOG(DEBUG) << "Seed case 1 computation \n"
//						   << "i1: " << i1 << "\n"
//						   << "j1: " << rightExt->j1 << "\n"
//						   << "i2: " << i2 << "\n"
//						   << "j1: " << rightExt->j2 << "\n";
				// get energy of seed interaction with best right extension
				curE = seedHandler.getSeedE(i1, i2) + rightExt->E;
				// check if this combination yields better energy
				curEtotal = energy.getE(i1, rightExt->j1, i2, rightExt->j2, curE);
				if (curEtotal < curCellEtotal) {
					// update current best for this left boundary
					// copy right boundary
					*curCell = *rightExt;
					// set new energy
					curCell->E = curE;
					// store total energy
					curCellEtotal = curEtotal;
				}
			}

			// compute overall energy of seed+multiside
			if ( i1+w1 < hybridE_multi.size1()
			    && i2+w2 < hybridE_multi.size2()) {
//				LOG(DEBUG) << "Entered SeedHandling case Multi!";
				rightExt = &(hybridE_multi(i1 + w1, i2 + w2));
				// check if right side can pair
				if (E_isINF(rightExt->E)) {
					continue;
				}
				// check if interaction length is within boundary
				if ((rightExt->j1 + 1 - i1) > energy.getAccessibility1().getMaxLength()
					|| (rightExt->j2 + 1 - i2) > energy.getAccessibility2().getMaxLength()) {
					continue;
				}
//				LOG(DEBUG) << "Seed case 2 computation \n"
//						   << "i1: " << i1 << "\n"
//						  << "j1: " << rightExt->j1 << "\n"
//						  << "i2: " << i2 << "\n"
//						  << "j1: " << rightExt->j2 << "\n";
				// get energy of seed interaction with best right extension
				curE = seedHandler.getSeedE(i1, i2) + rightExt->E;
				// check if this combination yields better energy
				curEtotal = energy.getE(i1, rightExt->j1, i2, rightExt->j2, curE);

				if (curEtotal < curCellEtotal) {
					// update current best for this left boundary
					// copy right boundary
					*curCell = *rightExt;
					// set new energy
					curCell->E = curE;
					// store total energy
					curCellEtotal = curEtotal;
				}
			}
		}

//		LOG(DEBUG) << "At updateOptima: \n"
//				   << "i1 : " << i1 << "\n"
//				   << "j1 : " << curCell->j1 << "\n"
//				   << "j1 - Multi : " << curCellMulti->j1 << "\n"
//				   << "i2 : " << i2 << "\n"
//				   << "j2 : " << curCell->j2 << "\n"
//		           << "j2 - Multi : " << curCellMulti->j2 << "\n"
//				   << "CurCellETotal : " << curCellEtotal << "\n";

		// update mfe if needed (call superclass update routine)
		PredictorMfe2dMultiHeuristic::updateOptima( i1,curCell->j1, i2,curCell->j2, curCellEtotal, false );

	} // i2
	} // i1


	// report mfe interaction
	reportOptima( outConstraint );

}

////////////////////////////////////////////////////////////////////////////
size_t
PredictorMfe2dMultiHeuristicSeed::
traceHybridO( const size_t i1, const size_t j1
		, const size_t i2, const size_t j2 ) const
{
	E_type curE = hybridO(i1,i2).E;
	const BestInteraction * curCell = NULL;

	size_t k2;
	for (k2 = j2; k2 > i2 + InteractionEnergy::minDistES; k2--) {
		curCell = &(hybridE_seed(i1,k2));
		if (curCell->j1 == j1 && curCell->j2 == j2 &&
			E_equal(curE, energy.getE_multiRight(i1, i2, k2)
						  + curCell->E))
		{
			return k2;
		}
	}
	throw std::runtime_error("PredictorMfe2dMultiHeuristicSeed::traceHybridO() : can not trace k2 ");
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dMultiHeuristicSeed::
traceBack( Interaction & interaction )
{
	LOG(DEBUG) << "Entered traceback!";
	// check if something to trace
	if (interaction.basePairs.size() < 2) {
		return;
	}

#if INTARNA_IN_DEBUG_MODE
	// sanity checks
	if ( ! interaction.isValid() ) {
		throw std::runtime_error("PredictorMfe2dMultiHeuristicSeed::traceBack() : given interaction not valid");
	}
	if ( interaction.basePairs.size() != 2 ) {
		throw std::runtime_error("PredictorMfe2dMultiHeuristicSeed::traceBack() : given interaction does not contain boundaries only");
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
			i2 = energy.getIndex2(interaction.basePairs.at(0));
	const size_t j1 = energy.getIndex1(interaction.basePairs.at(1));
	const size_t j2 = energy.getIndex2(interaction.basePairs.at(1));


	// the currently traced value for i1-j1, i2-j2
	E_type curE = hybridE_seed(i1,i2).E;
	assert( hybridE_seed(i1,i2).j1 == j1 );
	assert( hybridE_seed(i1,i2).j2 == j2 );
	assert( i1 <= j1 );
	assert( i2 <= j2 );
	assert( j1 < hybridE_seed.size1() );
	assert( j2 < hybridE_seed.size2() );

	// trace back
	// temp variables
	size_t k1,k2;
	bool doLocalTrace = true;
	bool traceInESeed = true; // if false: traces hybridE_multi
	// do until only right boundary is left over
	while( (j1-i1) > 1 ) {
		const BestInteraction * curCell = NULL;

		// check if we still have to find the seed
        if (doLocalTrace) {

			// check seed extended to the right
			if (traceInESeed && E_isNotINF( seedHandler.getSeedE(i1,i2))) {

				k1 = i1+seedHandler.getSeedLength1(i1,i2)-1;
				k2 = i2+seedHandler.getSeedLength2(i1,i2)-1;

				// check if correct trace with single-side extension
				curCell = &(hybridE(k1,k2));
				if ( E_equal( curE, (seedHandler.getSeedE(i1,i2)+curCell->E) )) {
					// store seed information
					interaction.setSeedRange(
							energy.getBasePair(i1, i2),
							energy.getBasePair(k1, k2),
							energy.getE(i1, k1, i2, k2, seedHandler.getSeedE(i1, i2)) + energy.getE_init());
					// traceback seed base pairs (excludes right most = (k1,k2))
					seedHandler.traceBackSeed(interaction, i1, i2);
					// update position to point after seed interaction
					i1 = k1;
					i2 = k2;
					curE = curCell->E;
					doLocalTrace = false;
					continue;
				}

				// check if correct trace with multi-side extension
				curCell = &(hybridE_multi(k1,k2));
				if ( E_equal( curE, (seedHandler.getSeedE(i1,i2)+(&(hybridE_multi(k1,k2)))->E) )) {
					// store seed information
					interaction.setSeedRange(
							energy.getBasePair(i1, i2),
							energy.getBasePair(k1, k2),
							energy.getE(i1, k1, i2, k2, seedHandler.getSeedE(i1, i2)) + energy.getE_init());
					// traceback seed base pairs (excludes right most = (k1,k2))
					seedHandler.traceBackSeed(interaction, i1, i2);
					// update position to point after seed interaction
					i1 = k1;
					i2 = k2;
					curE = curCell->E;
					traceInESeed = false;
					continue;
				}
			}
		}
	    // check all interval splits if no trace already found
	    if ( doLocalTrace ) {

			// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
			bool traceNotFound = true;
			for (k1=std::min(j1-seedHandler.getConstraint().getBasePairs()+1,i1+energy.getMaxInternalLoopSize1()+1); traceNotFound && k1>i1; k1--) {
			for (k2=std::min(j2-seedHandler.getConstraint().getBasePairs()+1,i2+energy.getMaxInternalLoopSize2()+1); traceNotFound && k2>i2; k2--) {
				if ( traceInESeed
					&& j1-k1 < hybridE_seed.size1()
					&& j2-k2 < hybridE_seed.size2()
					&& (E_isNotINF((&(hybridE_seed(k1,k2)))->E))) {
					// temp access to current cell
					curCell = &(hybridE_seed(k1,k2));
					// check if right boundary is equal (part of the heuristic)
					if ( curCell->j1 == j1 && curCell->j2 == j2 &&
						 // and energy is the source of curE
						 E_equal( curE, (energy.getE_interLeft(i1,k1,i2,k2) + curCell->E ) ) )
					{
						// stop searching
						traceNotFound = false;
						// store splitting base pair
						if ( k1 < j1 ) {
							interaction.basePairs.push_back( energy.getBasePair(k1,k2) );
						}
						// trace right part of split
						i1=k1;
						i2=k2;
						curE = curCell->E;
						continue;
					}
				} else
					// check if (k1,k2) are valid left boundaries including a seed
				if ( !traceInESeed
					&& j1-k1 < hybridE_multi.size1()
					&& j2-k2 < hybridE_multi.size2()
					&& (E_isNotINF((&(hybridE_multi(k1,k2)))->E))) {

					curCell = &(hybridE_multi(k1,k2));
					// check if right boundary is equal (part of the heuristic)
					if ( curCell->j1 == j1 && curCell->j2 == j2 &&
						 // and energy is the source of curE
						 E_equal( curE, (energy.getE_interLeft(i1,k1,i2,k2) + curCell->E ) ) )
					{
						// stop searching
						traceNotFound = false;
						// store splitting base pair
						if ( k1 < j1 ) {
							interaction.basePairs.push_back( energy.getBasePair(k1,k2) );
						}
						// trace right part of split
						i1=k1;
						i2=k2;
						curE = curCell->E;
						continue;
					}
				}

			} // k2
			} // k1

			///////////////  multi-side trace  ///////////////////

			// Both-sided structure
			if (!traceInESeed && traceNotFound && allowES == ES_both) {
				for (k1 = j1; traceNotFound && k1 > i1 + InteractionEnergy::minDistES; k1--) {
					// temp access to current cell
					curCell = &(hybridO(k1, i2));
					// check if right boundary is equal (part of the heuristic)
					if (curCell->j1 == j1 && curCell->j2 == j2 &&
						// and energy is the source of curE
						E_equal(curE, (energy.getE_multiLeft(i1, k1, i2,
															 InteractionEnergy::ES_multi_mode::ES_multi_both) +
									   curCell->E))) {
						// stop searching
						traceNotFound = false;
						traceInESeed = true;
						E_type E_multiRIght = energy.getE_multiRight(i1, i2, k2);
						// Determine k2 based on k1
						k2 = traceHybridO(k1, j1, i2, j2);
						// store splitting base pair if not last one of interaction range
						if (k1 < j1) {
							interaction.basePairs.push_back(energy.getBasePair(k1, k2));
						}

						// store gap information
						if (interaction.gap == NULL) { interaction.gap = new Interaction::Gap(); }
						interaction.gap->energy += energy.getE_multiLeft(i1, k1, i2, InteractionEnergy::ES_multi_mode::ES_multi_both) + E_multiRIght;
						Interaction::BasePair bpLeft = energy.getBasePair(i1,i2);
						interaction.gap->gaps1.insert( IndexRange(bpLeft.first, interaction.basePairs.rbegin()->first) );
						interaction.gap->gaps2.insert( IndexRange(interaction.basePairs.rbegin()->second,bpLeft.second) );
						// trace right part of split
						i1 = k1;
						i2 = k2;
						curE = hybridE_seed(k1, k2).E;
						continue;
					}
				}
			}

			// Structure in S1
			if (!traceInESeed && traceNotFound && (allowES == ES_target || allowES == ES_xorQueryTarget)) {
				for (k1 = j1; traceNotFound && k1 > i1 + InteractionEnergy::minDistES; k1--) {
				for (k2 = std::min(j2, i2 + energy.getMaxInternalLoopSize2() + 1); traceNotFound && k2 > i2; k2--) {
					// temp access to current cell
					curCell = &(hybridE_seed(k1, k2));
					// check if right boundary is equal (part of the heuristic)
					if (curCell->j1 == j1 && curCell->j2 == j2 &&
						// and energy is the source of curE
						E_equal(curE, (energy.getE_multi(i1, k1, i2, k2,
														 InteractionEnergy::ES_multi_mode::ES_multi_1only) +
									   curCell->E))) {
						// stop searching
						traceNotFound = false;
						traceInESeed = true;

						// store splitting base pair if not last one of interaction range
						if (k1 < j1) {
							interaction.basePairs.push_back(energy.getBasePair(k1, k2));
						}

						// store gap information
						if (interaction.gap == NULL) { interaction.gap = new Interaction::Gap(); }
						interaction.gap->energy += energy.getE_multi(i1, k1, i2, k2,
																	 InteractionEnergy::ES_multi_mode::ES_multi_1only);
						Interaction::BasePair bpLeft = energy.getBasePair(i1, i2);
						interaction.gap->gaps1.insert(
								IndexRange(bpLeft.first, interaction.basePairs.rbegin()->first));

						// trace right part of split
						i1 = k1;
						i2 = k2;
						curE = curCell->E;
						continue;
					}
					}
				}
			}

			// Structure in S2
			if (!traceInESeed && traceNotFound && (allowES == ES_query || allowES == ES_xorQueryTarget)) {
				for (k1 = std::min(j1, i1 + energy.getMaxInternalLoopSize1() + 1); traceNotFound && k1 > i1; k1--) {
					// temp access to current cell
					curCell = &(hybridO(k1, i2));
					// check if right boundary is equal (part of the heuristic)
					if (curCell->j1 == j1 && curCell->j2 == j2 &&
						// and energy is the source of curE
						E_equal(curE, (energy.getE_multiLeft(i1, k1, i2,
															 InteractionEnergy::ES_multi_mode::ES_multi_2only) +
									   curCell->E))) {
						// stop searching
						traceNotFound = false;
						// Determine k2 based on k1
						k2 = traceHybridO(k1, j1, i2, j2);
						E_type E_multiRight = energy.getE_multiRight(i1, i2, k2);
						// store splitting base pair if not last one of interaction range
						if (k1 < j1) {
							interaction.basePairs.push_back(energy.getBasePair(k1, k2));
						}
						// store gap information
						if (interaction.gap == NULL) { interaction.gap = new Interaction::Gap(); }
						interaction.gap->energy += energy.getE_multiLeft(i1, k1, i2, InteractionEnergy::ES_multi_mode::ES_multi_2only) + E_multiRight;
						Interaction::BasePair bpLeft = energy.getBasePair(i1,i2);
						interaction.gap->gaps2.insert( IndexRange(interaction.basePairs.rbegin()->second,bpLeft.second) );
						// trace right part of split
						i1 = k1;
						i2 = k2;
						curE = hybridE_seed(k1, k2).E;
						continue;
					}
				}


				// final sanity check
				assert(!traceNotFound);
			} // local trace
		}
	       // seed was already traced, do "normal" interaction trace
	    else {
			// create temporary data structure to be filed
			Interaction bpsRight(*(interaction.s1), *(interaction.s2) );
			bpsRight.basePairs.push_back( energy.getBasePair(i1,i2) );
			bpsRight.basePairs.push_back( energy.getBasePair(j1,j2) );
			PredictorMfe2dMultiHeuristic::traceBack( bpsRight );
			// copy remaining base pairs
			Interaction::PairingVec & bps = bpsRight.basePairs;
			// copy all base pairs excluding the right most
			for (size_t i=0; i+1<bps.size(); i++) {
				interaction.basePairs.push_back( bps.at(i) );
			}
			i1 = j1;
			i2 = j2;
			// stop traceback
			break;
		}


	}
#if INTARNA_IN_DEBUG_MODE
	if ( (j2-i2) > 1 ) {
		throw std::runtime_error("PredictorMfe2dMultiHeuristicSeed::traceBack() : trace leaves ji<j2 : "+toString(i2)+"<"+toString(j2));
	}
#endif

	// sort final interaction (to make valid) (faster than calling sort())
	if (interaction.basePairs.size() > 2) {
		Interaction::PairingVec & bps = interaction.basePairs;
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
PredictorMfe2dMultiHeuristicSeed::
getNextBest( Interaction & curBest )
{

	const E_type curBestE = curBest.energy;

	// TODO replace index iteration with something based on ranges from reportedInteractions

	// identify cell with next best non-overlapping interaction site
	// iterate (decreasingly) over all left interaction starts
	size_t i1,i2;
	BestInteraction * curBestCell = NULL;
	E_type curBestCellE = E_INF;
	Interaction::BasePair curBestCellStart;
	BestInteraction * curCell = NULL;
	E_type curCellE = E_INF;
	IndexRange r1,r2;
	for (i1=hybridE_seed.size1(); i1-- > 0;) {
		// ensure interaction site start is not covered
		if (reportedInteractions.first.covers(i1)) {
			continue;
		}
		for (i2=hybridE_seed.size2(); i2-- > 0;) {
			// ensure interaction site start is not covered
			if (reportedInteractions.second.covers(i2)) {
				continue;
			}
			// direct cell access
			curCell = &(hybridE_seed(i1,i2));
			// check if left side can pair
			if (E_isINF(curCell->E))
			{
				continue;
			}
			// get overall energy of the interaction
			curCellE = energy.getE(i1,curCell->j1,i2,curCell->j2,curCell->E);
			// or energy is too low to be considered
			// or energy is higher than current best found so far
			if (curCellE < curBestE || curCellE >= curBestCellE )
			{
				continue;
			}
			// ensure site is not overlapping
			r1.from = i1;
			r1.to = curCell->j1;
			if ( reportedInteractions.first.overlaps( r1 )) {
				continue;
			}
			r2.from = i2;
			r2.to = curCell->j2;
			if ( reportedInteractions.second.overlaps( r2 )) {
				continue;
			}
			//// FOUND THE NEXT BETTER SOLUTION
			// overwrite current best found so far
			curBestCell = curCell;
			curBestCellE = curCellE;
			curBestCellStart.first = i1;
			curBestCellStart.second = i2;

		} // i2
	} // i1

	// overwrite curBest
	curBest.energy = curBestCellE;
	curBest.basePairs.resize(2);
	if (E_isNotINF(curBestCellE)) {
		curBest.basePairs[0] = energy.getBasePair( curBestCellStart.first, curBestCellStart.second );
		curBest.basePairs[1] = energy.getBasePair( curBestCell->j1, curBestCell->j2 );
	}

}

////////////////////////////////////////////////////////////////////////////

} // namespace