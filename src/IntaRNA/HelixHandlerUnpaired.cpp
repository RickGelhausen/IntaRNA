#include "IntaRNA/HelixHandlerUnpaired.h"

namespace IntaRNA {

/////////////////////////////////////////////////////////////////////////////

size_t
HelixHandlerUnpaired::
fillHelix(const size_t i1min, const size_t i1max, const size_t i2min, const size_t i2max)
{

#if INTARNA_IN_DEBUG_MODE
	if ( i1min > i1max ) throw std::runtime_error("HelixHandlerUnpaired::fillHelix: i1min("+toString(i1min)+") > i1max("+toString(i1max)+")");
	if ( i2min > i2max ) throw std::runtime_error("HelixHandlerUnpaired::fillHelix: i2min("+toString(i2min)+") > i2max("+toString(i2max)+")");
	if ( i1max > energy.size1() ) throw std::runtime_error("HelixHandlerUnpaired::fillHelix: i1max("+toString(i1max)+") > energy.size1("+toString(energy.size1())+")");
	if ( i2max > energy.size2() ) throw std::runtime_error("HelixHandlerUnpaired::fillHelix: i2max("+toString(i2max)+") > energy.size2("+toString(energy.size2())+")");
	if ( helixConstraint.getMinBasePairs() > helixConstraint.getMaxBasePairs() )
		throw std::runtime_error("HelixHandlerUnpaired::fillHelix: bpMin("+toString(helixConstraint.getMinBasePairs()) +") > bpMax("+toString(helixConstraint.getMaxBasePairs())+")");
#endif


	helix.resize( i1max-i1min+1, i2max-i2min+1 );
	helixE_rec.resize( HelixIndex({{
				   (HelixRecMatrix::index)(helix.size1())
				   , (HelixRecMatrix::index)(helix.size2())
				   , (HelixRecMatrix::index)(getConstraint().getMaxBasePairs()+1)
				   , (HelixRecMatrix::index)(getConstraint().getMaxUnpaired()+1)
				   , (HelixRecMatrix::index)(getConstraint().getMaxUnpaired()+1)}}) );

	// store index offset due to restricted matrix size generation
	offset1 = i1min;
	offset2 = i2min;

	// temporary variables
	size_t i1, i2, curBP, bestBP, u1, u2, j1, j2, u1p, u2p, k1, k2, u1best, u2best;
	E_type curE, bestE;

	size_t  helixCountNotInf = 0, helixCount = 0;

	// fill for all start indeices
	// in decreasing index order
	for (i1=i1max+1; i1-- > i1min;) {
	for (i2=i2max+1; i2-- > i2min;) {

		// count possible helices
		helixCount++;

		// init according to no helix interaction
		helix(i1-offset1,i2-offset2) = HelixMatrix::value_type( E_INF, 0, 0 );

		// skip non-complementary left helix boundaries
		if (!energy.areComplementary(i1,i2)) {
			continue; // go to next helixE index
		}

		// Calculate energy for all different numbers of base pairs (bpMin to bpMax)
		for (curBP=2; curBP < getConstraint().getMaxBasePairs()+1
					  && (i1+curBP-1-offset1) < helix.size1()
					  && (i2+curBP-1-offset2) < helix.size2(); curBP++) {

			// for feasible unpaired in seq1 in increasing order
			for (u1 = 0; u1 < helixE_rec.shape()[3] && (i1 + curBP - 1 + u1 - offset1) < helix.size1(); u1++) {

			// for feasible unpaired in seq2 in increasing order
			for (u2 = 0; u2 < helixE_rec.shape()[4] - u1 && (i2 + curBP - 1 + u2 - offset2) < helix.size2(); u2++) {

				// get right helix boundaries
				j1 = i1 + curBP - 1 + u1;
				j2 = i2 + curBP - 1 + u2;

				// init current helix energy
				curE = E_INF;

				// check if right boundary is complementary
				if (energy.areComplementary(j1, j2)) {

					// base case: only left and right base pair present
					if (curBP == 2) {
						// energy for stacking/bulge/interior depending on u1/u2
						curE = energy.getE_interLeft(i1, j1, i2, j2);
//						LOG(DEBUG) << "i1, i2, curBP, u1, u2: " << i1 << " " << i2 << " " << curBP << " " << u1 << " " << u2;
//						LOG(DEBUG) << "curE: " << curE;
						} else {
						// split helix recursively into all possible leading interior loops
						// i1 .. i1+u1p+1 .. j1
						// i2 .. i2+u2p+1 .. j2
						for (u1p = 1 + std::min(u1, energy.getMaxInternalLoopSize1()); u1p-- > 0;) {
						for (u2p = 1 + std::min(u2, energy.getMaxInternalLoopSize2()); u2p-- > 0;) {

							k1 = i1 + u1p + 1;
							k2 = i2 + u2p + 1;
							// check if split pair is complementary
							// and recursed entry is < E_INF
							if (!(energy.areComplementary(k1, k2) && E_isNotINF(
									getHelixE(k1 - offset1, k2 - offset2, curBP - 1, u1 - u1p, u2 - u2p)))) {
								continue; // not complementary -> skip
							}

//							LOG(DEBUG) << "i1, i2, curBP, u1, u2: " << i1 << " " << i2 << " " << curBP << " " << u1 << " " << u2;
//							LOG(DEBUG) << "curE: " << energy.getE_interLeft(i1, k1, i2, k2)
													  + getHelixE(k1 - offset1, k2 - offset2, curBP - 1, u1 - u1p, u2 - u2p);
							// update mfe for split at k1,k2
							curE = std::min(curE,
											energy.getE_interLeft(i1, k1, i2, k2)
											+ getHelixE(k1 - offset1, k2 - offset2, curBP - 1, u1 - u1p, u2 - u2p)
							);

						} // u2p
						} // u1p
					} // more than two base pairs
				} // (j1, j2) complementary

				// store helix energy
				setHelixE(i1 - offset1, i2 - offset2, curBP, u1, u2, curE);
			} // u2
			} // u1


		} // curBP

		// TODO: Possible runtime improvement by checking this while calculating
		// find best unpaired combination in helx for i1,i2,bp
		u1best = 0;
		u2best = 0;
		bestBP = 2;
		bestE = E_INF;

//		LOG(DEBUG) << "FindBEST: (i1,i2) " << i1 << " " << i2;
		// Calculate energy for all different numbers of base pairs
		// Ensuring minimum number of base pairs here
		for (curBP = helixConstraint.getMinBasePairs(); curBP < helixConstraint.getMaxBasePairs() + 1
						&& (i1 + curBP - 1 - offset1) < helix.size1()
						&& (i2 + curBP - 1 - offset2) < helix.size2(); curBP++) {

			// for feasible unpaired in seq1 in increasing order
			for (u1 = 0; u1 < helixE_rec.shape()[3] && (i1 + curBP - 1 + u1 - offset1) < helix.size1(); u1++) {
			// for feasible unpaired in seq2 in increaing order
			for (u2 = 0; u2 < helixE_rec.shape()[4]-u1 && (i2 + curBP - 1 + u2 - offset2) < helix.size2(); u2++) {

				// get right helix boundaries
				j1 = i1 + curBP - 1 + u1;
				j2 = i2 + curBP - 1 + u2;

				// get overall interaction energy
				curE = energy.getE(i1, j1, i2, j2, getHelixE(i1 - offset1, i2 - offset2, curBP, u1, u2)) +
					   energy.getE_init();

//				LOG(DEBUG) << "i1, i2, curBP, u1, u2, curE " << i1 << " " << i2 << " " << curBP << " " << u1 << " " << u2 << " " << curE;
				// check if better than what is known so far
				if (curE < bestE) {
					bestE = curE;
					u1best = u1;
					u2best = u2;
					bestBP = curBP;
				}
			} // u2
			} // u1
		} // curBP
		// reduce bestE to hybridization energy only (init+loops)
		if (E_isNotINF(bestE)) {
			// get helix hybridization loop energies only
			bestE = getHelixE(i1 - offset1, i2 - offset2, bestBP, u1best, u2best);
			// count true helix
			helixCountNotInf++;
		}

		// store best (mfe) helix for all u1/u2
		helix(i1 - offset1, i2 - offset2) = HelixMatrix::value_type(bestE,
																	E_isINF(bestE) ? 0 : encodeHelixLength(
																			bestBP + u1best,
																			bestBP + u2best),
																	E_isINF(bestE) ? 0 : bestBP);


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
HelixHandlerUnpaired::
traceBackHelix( Interaction & interaction
			, const size_t i1_
			, const size_t i2_
			, const size_t bestBP
			, const size_t u1_
			, const size_t u2_)
{
	// get boundaries
	size_t 	  i1 = i1_
			, i2 = i2_
			, u1max = u1_
			, u2max = u2_
			, uMax = helixConstraint.getMaxUnpaired()
			, u1, u2
			, k1, k2
			;

	// get energy of provided seed
	E_type curE = getHelixE(i1_,i2_,bestBP,u1_,u2_);

	// trace helices
	// trace each helix base pair (excluding right most)
	for ( size_t bp=1+bestBP; bp-- > 2; ) {
		// base case: only left and right base pair present
		if (bp==2) {
			// add left base pair if not left helix boundary
			if (i1 != i1_) {
				interaction.basePairs.push_back( energy.getBasePair( i1+offset1, i2+offset2 ) );
			}

		} else {
			// split helix recursively into all possible leading interior loops
			// i1 .. i1+u1p+1 .. j1
			// i2 .. i2+u2p+1 .. j2
			bool traceNotFound = true;
			for (u1=1+u1max; traceNotFound && u1-- > 0;) {
			for (u2=1+u2max; traceNotFound && u2-- > 0;) {
				// TODO umax might be useless here
				// check if overall number of unpaired is not exceeded
				if (u1+u2 > uMax) {
					continue;
				}
				k1 = i1+u1+1;
				k2 = i2+u2+1;

				// check if valid trace
				if ( E_isNotINF( getHelixE( k1, k2, bp-1, u1max-u1, u2max-u2) ) ) {

					// check if correct trace
					if ( E_equal( curE, energy.getE_interLeft(i1+offset1, k1+offset1, i2+offset2, k2+offset2)
										+ getHelixE( k1, k2, bp-1, u1max-u1, u2max-u2 )) )
					{
						// store left base pair if not left helix boundary
						if (i1 != i1_) {
							interaction.basePairs.push_back( energy.getBasePair(i1+offset1, i2+offset2) );
						}
						// store next energy value in trace
						curE = getHelixE( k1, k2, bp-1, u1max-u1, u2max-u2 );
						// reset for next trace step
						i1 = k1;
						i2 = k2;
						// update boundaries for unpaired positions to reduce trace effort
						u1max -= u1;
						u2max -= u2;
						uMax -= (u1+u2);
						// mark trace step done
						traceNotFound = false;
					}
				}

			} // u2
			} // u1
			assert( !traceNotFound ); // sanity check
		} // more than two base pairs

	} // bp
}


} // namespace


