
#include "IntaRNA/PredictorMfe2dLimStackHeuristic.h"

#include <stdexcept>

namespace IntaRNA {

////////////////////////////////////////////////////////////////////////////

PredictorMfe2dLimStackHeuristic::
PredictorMfe2dLimStackHeuristic( const InteractionEnergy & energy
		, OutputHandler & output
		, PredictionTracker * predTracker
		, const HelixConstraint & helixConstraint
		)
	: PredictorMfe2dHeuristic(energy,output,predTracker)
		, helixHandler( energy, helixConstraint )
{
}


////////////////////////////////////////////////////////////////////////////

PredictorMfe2dLimStackHeuristic::
~PredictorMfe2dLimStackHeuristic()
{
	// clean up
}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dLimStackHeuristic::
predict( const IndexRange & r1
		, const IndexRange & r2
		, const OutputConstraint & outConstraint
)
{
#if INTARNA_MULITHREADING
#pragma omp critical(intarna_omp_logOutput)
#endif
	// TODO: Change
	{ VLOG(2) <<"predicting mfe interactions heuristically in O(n^2) space and time..."; }
	// measure timing
	TIMED_FUNC_IF(timerObj,VLOG_IS_ON(9));

#if INTARNA_IN_DEBUG_MODE
	// check indices
	if (!(r1.isAscending() && r2.isAscending()) )
		throw std::runtime_error("PredictorMfe2dLimStackHeuristic::predict("+toString(r1)+","+toString(r2)+") is not sane");
#endif


	// set index offset
	energy.setOffset1(r1.from);
	energy.setOffset2(r2.from);
	helixHandler.setOffset1(r1.from);
	helixHandler.setOffset2(r2.from);

	const size_t hybridEsize1 = std::min( energy.size1()
			, (r1.to==RnaSequence::lastPos?energy.size1()-1:r1.to)-r1.from+1 );
	const size_t hybridEsize2 = std::min( energy.size2()
			, (r2.to==RnaSequence::lastPos?energy.size2()-1:r2.to)-r2.from+1 );

	if (helixHandler.fillHelix( 0, hybridEsize1-1, 0, hybridEsize2-1 ) == 0) {
		// trigger empty interaction reporting
		initOptima(outConstraint);
		reportOptima(outConstraint);
		// stop computation
		return;
	}

	// resize matrix
	hybridE.resize( hybridEsize1, hybridEsize2 );

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
		} else {
			// set to infinity, ie not used
			hybridE(i1,i2) = BestInteraction(E_INF, RnaSequence::lastPos, RnaSequence::lastPos);
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
PredictorMfe2dLimStackHeuristic::
fillHybridE()
{
	// compute entries
	// current minimal value
	E_type curE = E_INF, curEtotal = E_INF, curCellEtotal = E_INF;
	size_t i1,i2,h1,h2,w1,w2;
	BestInteraction * curCell = NULL;
	const BestInteraction * rightExt = NULL;
	// iterate (decreasingly) over all left interaction starts
	for (i1=hybridE.size1(); i1-- > 0;) {
	for (i2=hybridE.size2(); i2-- > 0;) {
		// direct cell access
		curCell = &(hybridE(i1,i2));
		// check if left side can pair
		if (E_isINF(curCell->E)) {
			continue;
		}

		// E_init initialization
		curCellEtotal = energy.getE(i1,curCell->j1,i2,curCell->j2,curCell->E);

		// check if helix is possible for this left boundary
		if (E_isNotINF( helixHandler.getHelixE(i1,i2) ) ) {
			// get right extension
			h1 = helixHandler.getHelixLength1(i1,i2)-1; assert(i1+h1 < hybridE.size1());
			h2 = helixHandler.getHelixLength2(i1,i2)-1; assert(i2+h2 < hybridE.size2());

			curE = helixHandler.getHelixE(i1,i2) + energy.getE_init();

			// check if this combination yields better energy
			curEtotal = energy.getE(i1, i1+h1, i2, i2+h2, curE);
			if ( curEtotal < curCellEtotal )
			{
				// TODO: Check this "right" boundary
				// update current best for this left boundary
				// set right boundary
				curCell->j1 = i1+h1;
				curCell->j2 = i2+h2;
				// set new energy
				curCell->E = curE;
				// store total energy to avoid recomputation
				curCellEtotal = curEtotal;
			}

			/////////// extend stacking with bulge and further interaction

			// iterate over all loop sizes w1 (seq1) and w2 (seq2) (minus 1)
			for (w1=1; w1-1 <= energy.getMaxInternalLoopSize1() && i1+h1+w1<hybridE.size1(); w1++) {
			for (w2=1; w2-1 <= energy.getMaxInternalLoopSize2() && i2+h2+w2<hybridE.size2(); w2++) {

				// skip stacking to ensure enumerating only bulges and interior loops
				if (w1 == 1 && w2 == 1) {
					continue;
				}

				// direct cell access (const)
				rightExt = &(hybridE(i1+h1+w1,i2+h2+w2));
				// check if right side can pair
				if (E_isINF(rightExt->E)) {
					continue;
				}
				// check if interaction length is within boundary
				if ( (rightExt->j1 +1 -i1) > energy.getAccessibility1().getMaxLength()
					 || (rightExt->j2 +1 -i2) > energy.getAccessibility2().getMaxLength() )
				{
					continue;
				}

				// compute energy for this loop sizes
				curE = helixHandler.getHelixE(i1,i2) + energy.getE_interLeft(i1+h1,i1+h1+w1,i2+h2,i2+h2+w2) + rightExt->E;

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

			} // w2
			} // w1
		} // helix

		// update mfe if needed
		updateOptima( i1,curCell->j1, i2,curCell->j2, curCellEtotal, false );

	} // i2
	} // i1

}

////////////////////////////////////////////////////////////////////////////

void
PredictorMfe2dLimStackHeuristic::
traceBack( Interaction & interaction )
{
	// check if something to trace
	if (interaction.basePairs.size() < 2) {
		return;
	}

#if INTARNA_IN_DEBUG_MODE
	// sanity checks
	if ( ! interaction.isValid() ) {
		throw std::runtime_error("PredictorMfeLimStack2dHeuristic::traceBack() : given interaction not valid");
	}
	if ( interaction.basePairs.size() != 2 ) {
		throw std::runtime_error("PredictorMfeLimStack2dHeuristic::traceBack() : given interaction does not contain boundaries only");
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
	size_t h1,h2,k1,k2;
	// do until only right boundary is left over
	while( (j1-i1) > 1 ) {
		const BestInteraction * curCell = NULL;
		bool traceNotFound = true;

		assert(E_isNotINF(helixHandler.getHelixE(i1,i2)));
		h1 = helixHandler.getHelixLength1(i1,i2)-1; assert(h1 < hybridE.size1());
		h2 = helixHandler.getHelixLength2(i1,i2)-1; assert(h2 < hybridE.size2());

		// init case
		if ( E_equal(curE, helixHandler.getHelixE(i1,i2) + energy.getE_init() ) ) {
			// stop searching
			traceNotFound = false;
			// traceback helix base pairs ( excluding right most = (k1,k2))
			helixHandler.traceBackHelix(interaction, i1, i2);
			// trace right part of split
			i1=i1+h1;
			i2=i2+h2;
			curE = 0.0;
		}

		// check all combinations of decompositions into (i1,i2)..(k1,k2)-(j1,j2)
		for (size_t w1=1; traceNotFound && w1-1 <= energy.getMaxInternalLoopSize1() && i1+h1+w1<hybridE.size1(); w1++) {
		for (size_t w2=1; traceNotFound && w2-1 <= energy.getMaxInternalLoopSize2() && i2+h2+w2<hybridE.size2(); w2++) {

			// skip stacking
			if (w1 == 1 && w2 == 1) {
				continue;
			}

			k1 = i1+h1+w1;
			k2 = i2+h2+w2;

			// temp access to current cell
			curCell = &(hybridE(k1,k2));
			// check if right boundary is equal (part of the heuristic)
			if ( curCell->j1 == j1 && curCell->j2 == j2 &&
				 // and energy is the source of curE
				 E_equal( curE, (helixHandler.getHelixE(i1,i2) + energy.getE_interLeft(i1+h1,k1,i2+h2,k2) + curCell->E ) ) )
			{
					// stop searching
					traceNotFound = false;
					// store helix base pairs
					helixHandler.traceBackHelix( interaction, i1, i2 );
					// store splitting base pair if not last one of interaction range
					if ( k1 < j1 ) {
						interaction.basePairs.push_back( energy.getBasePair(k1,k2) );
					}
					// trace right part of split
					i1=k1;
					i2=k2;
					curE = curCell->E;
				}
			} // w1
			} // w2

		assert( !traceNotFound );
	}
#if INTARNA_IN_DEBUG_MODE
	if ( (j2-i2) > 1 ) {
		throw std::runtime_error("PredictorMfeLimStack2dHeuristic::traceBack() : trace leaves ji<j2 : "+toString(i2)+"<"+toString(j2));
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

} // namespace
