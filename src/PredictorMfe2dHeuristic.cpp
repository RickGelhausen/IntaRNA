
#include "PredictorMfe2dHeuristic.h"

#include <stdexcept>

// #include <omp.h>
#include <smmintrin.h>
#include <vector>

////////////////////////////////////////////////////////////////////////////

PredictorMfe2dHeuristic::
PredictorMfe2dHeuristic( const InteractionEnergy & energy, OutputHandler & output )
 : PredictorMfe(energy,output)
	, hybridE( 0,0 )
{
}


////////////////////////////////////////////////////////////////////////////

PredictorMfe2dHeuristic::
~PredictorMfe2dHeuristic()
{
	// clean up
}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dHeuristic::
predict( const IndexRange & r1
		, const IndexRange & r2
		, const OutputConstraint & outConstraint
		)
{

	VLOG(2) <<"predicting mfe interactions heuristically in O(n^2) space and time...";
	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

#if IN_DEBUG_MODE
	// check indices
	if (!(r1.isAscending() && r2.isAscending()) )
		throw std::runtime_error("PredictorMfe2dHeuristic::predict("+toString(r1)+","+toString(r2)+") is not sane");
#endif


	// set index offset
	energy.setOffset1(r1.from);
	energy.setOffset2(r2.from);

	// resize matrix
	hybridE.resize( std::min( energy.size1()
						, (r1.to==RnaSequence::lastPos?energy.size1()-1:r1.to)-r1.from+1 )
				, std::min( energy.size2()
						, (r2.to==RnaSequence::lastPos?energy.size2()-1:r2.to)-r2.from+1 ) );

	// temp vars
	size_t i1,i2,w1,w2;

	// init matrix
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
//		// alternative init obsolete due to default constructor used in matrix cell creation
//		} else {
//			// set to infinity, ie not used
//			hybridE(i1,i2) = BestInteraction(E_INF, RnaSequence::lastPos, RnaSequence::lastPos);
		}

	} // i2
	} // i1

	// init mfe for later updates
	initOptima( outConstraint );

	// compute table and update mfeInteraction
	fillHybridE();

	// trace back and output handler update
	reportOptima( outConstraint );

}


////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dHeuristic::
fillHybridE()
{
	// std::cout << "foo" << std::endl;
	// compute entries
	// current minimal value
	E_type curE = E_INF, curEtotal = E_INF, curCellEtotal = E_INF;
	__m128 curE_SSE = {E_INF, E_INF, E_INF, E_INF};
	std::vector<float> curEtotal_SSE(4,E_INF);
	__m128 curCellEtotal_SSE = {E_INF, E_INF, E_INF, E_INF};
	std::vector<size_t> rightExt_j1_SSE(4,0);
	std::vector<size_t> rightExt_j2_SSE(4,0);
	std::vector<float> result(4,0);
	E_type energy0 = E_INF;
	E_type energy1 = E_INF;
	E_type energy2 = E_INF;
	E_type energy3 = E_INF;
	


	
	size_t i1,i2,w1,w2;
	BestInteraction * curCell = NULL;
	const BestInteraction * rightExt = NULL;
	const BestInteraction * rightExt_SSE1 = NULL;
	const BestInteraction * rightExt_SSE2  = NULL;
	const BestInteraction * rightExt_SSE3  = NULL;
	// std::cout << __LINE__ << std::endl;
	
	// iterate (decreasingly) over all left interaction starts
	std::vector< BestInteraction* > curCellVector(4, NULL);
	std::vector< E_type > curCellEtotalVector(4, 0.0);
	size_t loopCounter = 0;
	for (i1=hybridE.size1(); i1-- > 0;) {
	for (i2=hybridE.size2(); i2-- > 0;) {
		// direct cell access
		curCell = &(hybridE(i1,i2));
		curCellVector[0] = &(hybridE(i1,i2));
		curCellVector[1] = &(hybridE(i1,i2));
		curCellVector[2] = &(hybridE(i1,i2));
		curCellVector[3] = &(hybridE(i1,i2));
		
		// check if left side can pair
		if (E_isINF(curCell->E)) {
			continue;
		}
		
		// current
		curCellEtotal = energy.getE(i1,curCell->j1,i2,curCell->j2,curCell->E);
		curCellEtotalVector[0] = curCellEtotal;
		curCellEtotalVector[1] = curCellEtotal;
		curCellEtotalVector[2] = curCellEtotal;
		curCellEtotalVector[3] = curCellEtotal;
		
		// TODO PARALLELIZE THIS DOUBLE LOOP ?!
		// iterate over all loop sizes w1 (seq1) and w2 (seq2) (minus 1)
		// #pragma omp parallel
		for (w1=1; w1-1 <= energy.getMaxInternalLoopSize1() && i1+w1<hybridE.size1(); w1++) {
		for (w2=1; w2-1 <= energy.getMaxInternalLoopSize2() -4 && i2+w2<hybridE.size2() -4; w2 += 4) {
			loopCounter++;
			// direct cell access (const)
			rightExt = &(hybridE(i1+w1,i2+w2));
			rightExt_SSE1 = &(hybridE(i1+w1,i2+w2+1));
			rightExt_SSE2 = &(hybridE(i1+w1,i2+w2+2));
			rightExt_SSE3 = &(hybridE(i1+w1,i2+w2+3));
			// std::cout << __LINE__ << std::endl;
			if (E_isINF(rightExt->E) && E_isINF(rightExt_SSE1->E) && E_isINF(rightExt_SSE2->E) && E_isINF(rightExt_SSE3->E)) {
				continue;
			}
			// check if right side can pair
			if (E_isINF(rightExt->E)) {
				energy0 = E_INF;
				rightExt_j1_SSE[0] = i1;
				rightExt_j2_SSE[0] = i2;
			} else {
				energy0 = energy.getE_interLeft(i1,i1+w1,i2,i2+w2) + rightExt->E;
				rightExt_j1_SSE[0] = rightExt->j1;
				rightExt_j2_SSE[0] = rightExt->j2;
			}
			if (E_isINF(rightExt_SSE1->E)) {
				energy1 = E_INF;
				rightExt_j1_SSE[1] = i1;
				rightExt_j2_SSE[1] = i2;
			} else {
				energy1 = energy.getE_interLeft(i1,i1+w1,i2,i2+w2+1) + rightExt_SSE1->E;
				rightExt_j1_SSE[1] = rightExt_SSE1->j1;
				rightExt_j2_SSE[1] = rightExt_SSE1->j2;
			}
			if (E_isINF(rightExt_SSE2->E)) {
				energy2 = E_INF;
				rightExt_j1_SSE[2] = i1;
				rightExt_j2_SSE[2] = i2;
			} else {
				energy2 = energy.getE_interLeft(i1,i1+w1,i2,i2+w2+2) + rightExt_SSE2->E;
				rightExt_j1_SSE[2] = rightExt_SSE2->j1;
				rightExt_j2_SSE[2] = rightExt_SSE2->j2;
			}
			if (E_isINF(rightExt_SSE3->E)) {
				energy3 = E_INF;
				rightExt_j1_SSE[3] = i1;
				rightExt_j2_SSE[3] = i2;
			} else {
				energy3 = energy.getE_interLeft(i1,i1+w1,i2,i2+w2+3) + rightExt_SSE3->E;
				rightExt_j1_SSE[3] = rightExt_SSE3->j1;
				rightExt_j2_SSE[3] = rightExt_SSE3->j2;
			}
			// compute energy for this loop sizes
			// check if this combination yields better energy
			curE_SSE = _mm_setr_ps(energy0, energy1, energy2, energy3);
			 energy.getE_SSE(i1, rightExt_j1_SSE, i2, rightExt_j2_SSE, curE_SSE, curEtotal_SSE);
			// curEtotal = energy.getE(i1,rightExt->j1,i2,rightExt->j2,curE);
			// #pragma omp critical

			if (curEtotal_SSE[0] < curCellEtotalVector[0] && !E_isINF(energy0)) {
				*(curCellVector[0]) = *rightExt;
				curCellVector[0]->E = energy0;
				curCellEtotalVector[0] = curEtotal_SSE[0];
			}
			
			if (curEtotal_SSE[1] < curCellEtotalVector[1] && !E_isINF(energy1)) {
				*(curCellVector[1]) = *rightExt_SSE1;
				curCellVector[1]->E = energy1;
				curCellEtotalVector[1] = curEtotal_SSE[1];
			}
			
			if (curEtotal_SSE[2] < curCellEtotalVector[2] && !E_isINF(energy2)) {
				*(curCellVector[2]) = *rightExt_SSE2;
				curCellVector[2]->E = energy2;
				curCellEtotalVector[2] = curEtotal_SSE[2];
			}
			
			if (curEtotal_SSE[3] < curCellEtotalVector[3] && !E_isINF(energy3)) {
				*(curCellVector[3]) = *rightExt_SSE3;
				curCellVector[3]->E = energy3;
				curCellEtotalVector[3] = curEtotal_SSE[3];
			}
		} // w2
		
		// check eventually missing remainders
		if (loopCounter * 4 < w2) {
			for (w2=loopCounter; w2-1 <= energy.getMaxInternalLoopSize2() && i2+w2<hybridE.size2(); w2 += 1) {
				// direct cell access (const)
				rightExt = &(hybridE(i1+w1,i2+w2));
				// check if right side can pair
				if (E_isINF(rightExt->E)) {
					continue;
				}
				// compute energy for this loop sizes
				curE = energy.getE_interLeft(i1,i1+w1,i2,i2+w2) + rightExt->E;
				// check if this combination yields better energy
				curEtotal = energy.getE(i1,rightExt->j1,i2,rightExt->j2,curE);
				if ( curEtotal < curCellEtotal )
				{
					// update current best for this left boundary
					// copy right boundary
					*curCell = *rightExt;
					// set new energy
					curCell->E = curE;
					// store total energy to avoid recomputation
					curCellEtotal = curEtotal;
				}
			}
		}
		
		loopCounter = 0;
		} // w1
		
		// merge values from the two loops

		for (size_t i = 0; i < 4; ++i) {
			if (curCellEtotalVector[i] < curCellEtotal) {
				*curCell = *(curCellVector[i]);
				curCell->E = curCellVector[i]->E;
				curCellEtotal = curCellEtotalVector[i];
			}
		}

		// update mfe if needed
		updateOptima( i1,curCell->j1, i2,curCell->j2, curCellEtotal, false );
		

	} // i2
	} // i1

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dHeuristic::
traceBack( Interaction & interaction )
{
	// check if something to trace
	if (interaction.basePairs.size() < 2) {
		return;
	}

#if IN_DEBUG_MODE
	// sanity checks
	if ( ! interaction.isValid() ) {
		throw std::runtime_error("PredictorMfe2dHeuristic::traceBack() : given interaction not valid");
	}
	if ( interaction.basePairs.size() != 2 ) {
		throw std::runtime_error("PredictorMfe2dHeuristic::traceBack() : given interaction does not contain boundaries only");
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
	E_type curE = hybridE(i1,i2).E;
	assert( hybridE(i1,i2).j1 == j1 );
	assert( hybridE(i1,i2).j2 == j2 );
	assert( i1 <= j1 );
	assert( i2 <= j2 );
	assert( j1 < hybridE.size1() );
	assert( j2 < hybridE.size2() );

	// trace back
	// temp variables
	size_t k1,k2;
	// do until only right boundary is left over
	while( (j1-i1) > 1 ) {
		const BestInteraction * curCell = NULL;
		bool traceNotFound = true;
		// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
		for (k1=std::min(j1,i1+energy.getMaxInternalLoopSize1()+1); traceNotFound && k1>i1; k1--) {
		for (k2=std::min(j2,i2+energy.getMaxInternalLoopSize2()+1); traceNotFound && k2>i2; k2--) {
			// temp access to current cell
			curCell = &(hybridE(k1,k2));
			// check if right boundary is equal (part of the heuristic)
			if ( curCell->j1 == j1 && curCell->j2 == j2 &&
					// and energy is the source of curE
					E_equal( curE, (energy.getE_interLeft(i1,k1,i2,k2) + curCell->E ) ) )
			{
				// stop searching
				traceNotFound = false;
				// store splitting base pair if not last one of interaction range
				if ( k1 < j1 ) {
					interaction.basePairs.push_back( energy.getBasePair(k1,k2) );
				}
				// trace right part of split
				i1=k1;
				i2=k2;
				curE = curCell->E;
			}
		}
		}
		assert( !traceNotFound );
	}
#if IN_DEBUG_MODE
	if ( (j2-i2) > 1 ) {
		throw std::runtime_error("PredictorMfe2dHeuristic::traceBack() : trace leaves ji<j2 : "+toString(i2)+"<"+toString(j2));
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
PredictorMfe2dHeuristic::
getNextBest( Interaction & curBest )
{

	// get original
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
	for (i1=hybridE.size1(); i1-- > 0;) {
		// ensure interaction site start is not covered
		if (reportedInteractions.first.covers(i1)) {
			continue;
		}
		for (i2=hybridE.size2(); i2-- > 0;) {
			// ensure interaction site start is not covered
			if (reportedInteractions.second.covers(i2)) {
				continue;
			}
			// direct cell access
			curCell = &(hybridE(i1,i2));
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

