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
	helixE_rec.resize( HelixIndex({{
				   (HelixRecMatrix::index)(helix.size1())
				   , (HelixRecMatrix::index)(helix.size2())
				   , (HelixRecMatrix::index)(helixConstraint.getMaxBasePairs()+1)}}) );

	std::fill_n(helixE_rec.data(), helixE_rec.num_elements(), E_INF);

	// store index offset due to restricted matrix size generation
	offset1 = i1min;
	offset2 = i2min;

	// temporary variables
	size_t i1, i2, curBP, u1, u2, j1, j2, u1p, u2p, k1, k2, u1best, u2best, bestBP;
	E_type curE, bestE;

	size_t  helixCountNotInf = 0, helixCount = 0;

	// fill for all start indeices
	// in decreasing index order
	for (i1=i1max+1; i1-- > i1min;) {
	for (i2=i2max+1; i2-- > i2min;) {

		// count possible helices
		helixCount++;

		// init according to no helix interaction
		helix(i1-offset1,i2-offset2) = HelixMatrix::value_type( E_INF, 0 );

		// skip non-complementary left helix boundaries
		if (!energy.areComplementary(i1,i2)) {
			continue; // go to next helixE index
		}

		// Calculate energy for all different numbers of base pairs (bpMin to bpMax)
		for (curBP=2; curBP < helixConstraint.getMaxBasePairs()+1 && (i1+curBP-1-offset1)<helix.size1()
					  											  && (i2+curBP-1-offset2)<helix.size2()
																  && energy.areComplementary(i1+curBP-1,i2+curBP-1); curBP++) {

			// get right helix boundaries
			j1 = i1+curBP-1;
			j2 = i2+curBP-1;

			// init current helix energy
			curE = E_INF;

			// base case: only left and right base pair present
			if ( curBP==2 ) {
				// energy for stacking
				curE = energy.getE_interLeft(i1,j1,i2,j2);

			} else {

				k1 = i1+1;
				k2 = i2+1;

				// check if split pair is complementary
				// and recursed entry is < E_INF
				if (! ( energy.areComplementary(k1,k2) && E_isNotINF( getHelixE( k1-offset1, k2-offset2, curBP-1) ) ) ) {
					continue; // not complementary -> skip
				}

				// update mfe for split at k1,k2
				curE = std::min( curE,
								 energy.getE_interLeft(i1,k1,i2,k2)
								 + getHelixE( k1-offset1, k2-offset2, curBP-1)
				);
			} // more than two base pairs

			// store helix energy
			setHelixE( i1-offset1, i2-offset2, curBP, curE );

		} // curBP

		// find best combination in helix for i1,i2,bp
		bestE = E_INF;
		bestBP = 0;

		// Check all base pair combinations to find the best (> minBP)
		for ( size_t bp = helixConstraint.getMinBasePairs(); bp < helixConstraint.getMaxBasePairs()+1
															 && i1+bp-1-offset1 < helix.size1()
															 && i2+bp-1-offset2 < helix.size2(); bp++) {
			// get right helix boundaries
			j1 = i1+bp-1;
			j2 = i2+bp-1;

			// energy for current number of base pairs
			curE = energy.getE(i1, j1, i2, j2, getHelixE(i1-offset1,i2-offset2, bp)) + energy.getE_init();

			// check if better than what is known so far
			if ( curE < bestE ) {
				bestE = curE;
				bestBP = bp;
			}
		} // bp

		// reduce bestE to hybridization energy only (init+loops)
		if (E_isNotINF( bestE )) {
			// get helix hybridization loop energies only
			bestE = getHelixE( i1-offset1, i2-offset2, bestBP);
			// count true helix
			helixCountNotInf++;
		}

		// store best (mfe) helix for all u1/u2
		helix( i1-offset1, i2-offset2 ) = HelixMatrix::value_type( bestE
				, E_isINF(bestE)?0:bestBP);

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
	size_t 	  i1 = i1_
	, i2 = i2_;

	// Get base pair length of best interaction
	size_t numberOfBP = getHelixBP(i1_,i2_);

	std::cout << offset1 << " " << offset2 << std::endl;
	// Check whether this trace is within the allowed boundaries
	if (numberOfBP >= helixConstraint.getMinBasePairs() && numberOfBP <= helixConstraint.getMaxBasePairs()) {
		// trace helices
		// trace each helix base pair (excluding right most)
		for (size_t bp = 0; bp < numberOfBP; bp++) {
			if (i1 != i1_) {
				//TODO: ? interaction.basePairs.push_back(energy.getBasePair(i1+offset1, i2+offset2)
				interaction.basePairs.push_back(energy.getBasePair(i1, i2));
			}
			i1++;
			i2++;
		}
	}
}

} // namespace