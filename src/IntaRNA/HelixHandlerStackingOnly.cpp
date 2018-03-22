#include "IntaRNA/HelixHandlerStackingOnly.h"

namespace IntaRNA {

/////////////////////////////////////////////////////////////////////////////

size_t
HelixHandlerStackingOnly::
fillHelix(const size_t i1min, const size_t i1max, const size_t i2min, const size_t i2max)
{

#if INTARNA_IN_DEBUG_MODE
	if ( i1min > i1max ) throw std::runtime_error("HelixHandlerStackingOnly::fillHelix: i1min("+toString(i1min)+") > i1max("+toString(i1max)+")");
	if ( i2min > i2max ) throw std::runtime_error("HelixHandlerStackingOnly::fillHelix: i2min("+toString(i2min)+") > i2max("+toString(i2max)+")");
	if ( i1max > energy.size1() ) throw std::runtime_error("HelixHandlerStackingOnly::fillHelix: i1max("+toString(i1max)+") > energy.size1("+toString(energy.size1())+")");
	if ( i2max > energy.size2() ) throw std::runtime_error("HelixHandlerStackingOnly::fillHelix: i2max("+toString(i2max)+") > energy.size2("+toString(energy.size2())+")");
	if ( helixConstraint.getMinBasePairs() > helixConstraint.getMaxBasePairs() )
		throw std::runtime_error("HelixHandlerStackingOnly::fillHelix: bpMin("+toString(helixConstraint.getMinBasePairs()) +") > bpMax("+toString(helixConstraint.getMaxBasePairs())+")");
#endif

	helix.resize( i1max-i1min+1, i2max-i2min+1 );

	// store index offset due to restricted matrix size generation
	offset1 = i1min;
	offset2 = i2min;

	// temporary variables
	size_t i1, i2, curBP, j1, j2;

	size_t  helixCountNotInf = 0, helixCount = 0;

	// fill for all start indices
	// in increasing index order
	for (i1=i1min; i1 < i1max+1; i1++ ) {
	for (i2=i1min; i2 < i2max+1; i2++ ) {

		// count possible helices
		helixCount++;

		// init according to no helix interaction
		helix(i1-offset1,i2-offset2) = HelixMatrix::value_type( E_INF, 0 );

		// skip non-complementary left helix boundaries
		if (!energy.areComplementary(i1,i2)) {
			continue; // go to next helixE index
		}

		E_type leftHelixE = 0.0;
		// screen over all possible base pair combinations (starting at 2)
		for (curBP=2; curBP < helixConstraint.getMaxBasePairs()+1 && (i1+curBP-1-offset1<helix.size1())
					  											  && (i2+curBP-1-offset2<helix.size2())
					  											  && energy.areComplementary(i1+curBP-1,i2+curBP-1); curBP++) {

			j1 = i1+curBP-1;
			j2 = i2+curBP-1;

			// if next base pair is not complementary skip to next (i1,i2)
			if (!energy.areComplementary(j1,j2)) {
				break;
			}

			if (curBP == 2) {
				// energy for initial stacking
				leftHelixE = energy.getE_interLeft(i1,j1, i2,j2);
			} else {
				leftHelixE += energy.getE_interLeft(j1-1,j1, j2-1,j2);
			}

			// check whether this helix has best energy and is within lower boundary
			if (leftHelixE < helix(i1,i2).first && curBP >= helixConstraint.getMinBasePairs() && leftHelixE != 0.0 ) {
				helix(i1-offset1, i2-offset2) = HelixMatrix::value_type( leftHelixE, encodeHelixLength(curBP, curBP) );
			}
		}

		// Ensures that the helixCount is only inceased for the mfe helix.
		if (E_isNotINF(helix(i1-offset1, i2-offset2).first)) {
			// count possible helices
			helixCountNotInf++;
		}

	} // i2
	} // i1

#if INTARNA_MULITHREADING
#pragma omp critical(intarna_omp_logOutput)
#endif
	{ VLOG(2) << "valid helices = " << helixCountNotInf << " (" << (helixCountNotInf/helixCount) << "% of start index combinations)"; }

	return helixCountNotInf;
}

//////////////////////////////////////////////////////////////////////////

void
HelixHandlerStackingOnly::
traceBackHelix( Interaction & interaction
		, const size_t i1_
		, const size_t i2_)
{

	// get boundaries
	size_t 	  i1 = i1_-offset1
	, i2 = i2_-offset2;

	// Get base pair length of best interaction
	size_t numberOfBP = decodeHelixLength1(helix(i1,i2).second);

	// Check whether minimum requirement met to avoid unexpected behavior
	if (numberOfBP < helixConstraint.getMinBasePairs())
		return;

	// trace helices
	// trace each helix base pair (excluding right most)
	for (size_t bp = 0; bp < numberOfBP-1; bp++) {
		if (i1 != i1_-offset1) {
			interaction.basePairs.push_back(energy.getBasePair(i1+offset1, i2+offset2));
		}
		i1++;
		i2++;
	}

}

} // namespace