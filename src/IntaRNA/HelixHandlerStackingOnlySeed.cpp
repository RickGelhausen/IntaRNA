#include "IntaRNA/HelixHandlerStackingOnly.h"

namespace IntaRNA {

/////////////////////////////////////////////////////////////////////////////

size_t
HelixHandlerStackingOnly::
fillHelixSeed(const size_t i1min, const size_t i1max, const size_t i2min, const size_t i2max)
{
	helixSeed.resize( i1max-i1min+1, i2max-i2min+1 );

	// store index offset due to restricted matrix size generation
	offset1 = i1min;
	offset2 = i2min;

	// temporary variables
	size_t i1, i2,j1,j2, seedStart1, seedStart2, seedEnd1, seedEnd2, bestL1, bestL2, possibleBasePairs;
	size_t  helixCountNotInf = 0, helixCount = 0;

	E_type bestE, curE, bestRawE, rawE;

	// fill for all start indices
	// in decreasing index order
	for (i1=i1max+1; i1-- > i1min;) {
	for (i2=i2max+1; i2-- > i2min;) {

		// count possible helices
		helixCount++;

		// init according to no helix interaction
		helixSeed(i1 - offset1, i2 - offset2) = HelixMatrix::value_type(E_INF, 0);

		// skip non-complementary left helix boundaries
		if (!energy.areComplementary(i1, i2)) {
			continue; // go to next helixSeedE index
		}

		// TODO: THIS might make the boundary conditions useless
		// Check if a seed can fit given the left boundaries
		// Note: If seedHandler allows unpaired positions this check is not enough, check happens in loop
		if (std::min(helixSeed.size1()-i1+offset1, helixSeed.size2()-i2+offset2) < seedHandler->getConstraint().getBasePairs()) {
			continue;
		} else {
			// Seed fits, check how many bases are possible around
			possibleBasePairs = std::min(std::min(helixSeed.size1()-i1+offset1, helixSeed.size2()-i2+offset2), helixConstraint.getMaxBasePairs())-seedHandler->getConstraint().getBasePairs();
		}

		// Initialuze variables
		curE = E_INF;
		rawE = E_INF;
		bestE = E_INF;
		bestRawE = E_INF;
		bestL1 = 0;
		bestL2 = 0;

		// screen over all possible leading and trailing base pair combinations
		for (size_t leadingBP=0; leadingBP <= possibleBasePairs && (i1+leadingBP-offset1) < helixSeed.size1()
								 								&& (i2+leadingBP-offset2) < helixSeed.size2(); leadingBP++) {
			// If leading base pairs exist and helixE = E_INF -> skip to the next leadingBP
			if (leadingBP != 0 && E_isINF(getHelixE(i1-offset1,i2-offset2,leadingBP+1))) {
				continue;
			}

			// the start positions for the seed
			seedStart1 = i1+leadingBP;
			seedStart2 = i2+leadingBP;

			// If no seed is possible here, skip to next leading base pair number
			if (E_isINF(seedHandler->getSeedE(seedStart1, seedStart2))) {
				continue;
			}

			// the end positions of the seed
			seedEnd1 = seedStart1+seedHandler->getSeedLength1(seedStart1, seedStart2)-1;
			seedEnd2 = seedStart2+seedHandler->getSeedLength2(seedStart1, seedStart2)-1;

			// Run over all trailing base pairs
			for (size_t trailingBP=0; trailingBP <= possibleBasePairs - leadingBP
									  && (seedEnd1+trailingBP-offset1) < helixSeed.size1()
									  && (seedEnd2+trailingBP-offset2) < helixSeed.size2(); trailingBP++) {

				// If trailing base pairs exist and helixE = E_INF -> skip to the next leadingBP
				if (trailingBP != 0 && E_isINF(getHelixE(seedEnd1-offset1,seedEnd2-offset2,trailingBP+1))) {
					break;
				}

				j1 = seedEnd1+trailingBP;
				j2 = seedEnd2+trailingBP;

				// ensure that ED-values are within the boundaries (default 999)
				if (energy.getED1(i1, j1) > helixConstraint.getMaxED()
					|| energy.getED2(i2, j2) > helixConstraint.getMaxED()) {
					continue;
				}

				// energy value without dangling ends and E_init
				rawE = getHelixE(i1-offset1,i2-offset2,leadingBP+1) + seedHandler->getSeedE(seedStart1, seedStart2) + getHelixE(seedEnd1-offset1,seedEnd2-offset2,trailingBP+1);

				// energy value
				curE = energy.getE(i1,j1,i2,j2, rawE) + energy.getE_init();

				if (helixConstraint.noED())
					curE -= (energy.getED1(i1,j1) + energy.getED2(i2,j2));

				if (curE < bestE && !E_equal(curE, bestE)) {
					bestE = curE;
					bestRawE = rawE;
					bestL1 = leadingBP + seedHandler->getSeedLength1(seedStart1,seedStart2) + trailingBP;
					bestL2 = leadingBP + seedHandler->getSeedLength2(seedStart1,seedStart2) + trailingBP;
				}
			} // trailingBP
		} // leadingBP

		// reduce bestE to hybridization energy only (init+loops)
		if (E_isNotINF( bestE )) {
			// overwrite all helices with too high energy -> infeasible start interactions
			if (bestE > helixConstraint.getMaxE()) {
				bestE = E_INF;
			} else {
				// get helix hybridization loop energies only
				bestE = bestRawE;
				// count true helix
				helixCountNotInf++;
			}
		}
		// store best (mfe) helix for left boundary i1, i2
		helixSeed(i1-offset1, i2-offset2) = HelixMatrix::value_type(bestE, E_isINF(bestE) ? 0: encodeHelixSeedLength(bestL1,bestL2));

	} // i2
	} // i1

	return helixCountNotInf;
}

//////////////////////////////////////////////////////////////////////////

void
HelixHandlerStackingOnly::
traceBackHelixSeed( Interaction & interaction
		, const size_t i1_
		, const size_t i2_)
{
	size_t i1 = i1_
		 , i2 = i2_
		 , seedStart1, seedEnd1
	     , seedStart2, seedEnd2
		 , j1,j2;

	bool traceNotFound = true;

	E_type curE = getHelixSeedE(i1_,i2_);

	// No traceback possible for current boundary
	if (E_isINF(curE)) {
		return;
	}

	// Calculate how many base pairs are possible allongside the seed.
	// Note: If seedHandler allows unpaired positions this check is not enough, check happens in loop
	size_t possibleBasePairs = std::min(std::min(helixSeed.size1()-i1 +offset1, helixSeed.size2()-i2+offset2), helixConstraint.getMaxBasePairs())-seedHandler->getConstraint().getBasePairs();

	// screen over all possible leading and trailing base pair combinations
	for (size_t leadingBP=0; traceNotFound
							 && leadingBP <= possibleBasePairs
							 && i1 + leadingBP-offset1 < helixSeed.size1()
							 && i2 + leadingBP-offset2 < helixSeed.size2(); leadingBP++) {

		// If leading base pairs exist and helixE = E_INF -> skip to the next leadingBP
		if (leadingBP != 0 && E_isINF(getHelixE(i1-offset1,i2-offset2,leadingBP+1))) {
			continue;
		}

		// start positions of the seed
		seedStart1 = i1 + leadingBP;
		seedStart2 = i2 + leadingBP;

		// Check whether seed is possible for this starting position
		if (E_isINF(seedHandler->getSeedE(seedStart1, seedStart2))) {
			continue;
		}

		// end positions of the seed
		seedEnd1 = seedStart1+seedHandler->getSeedLength1(seedStart1,seedStart2)-1;
		seedEnd2 = seedStart2+seedHandler->getSeedLength2(seedStart1,seedStart2)-1;

		// Trailing base pairs
		for (size_t trailingBP = 0; traceNotFound
									&& trailingBP <= possibleBasePairs - leadingBP
									&& seedEnd1+trailingBP-offset1 < helixSeed.size1()
									&& seedEnd2+trailingBP-offset2 < helixSeed.size2(); trailingBP++) {

			// If trailing base pairs exist and helixE = E_INF -> skip to the next leadingBP
			if (trailingBP != 0 && E_isINF(getHelixE(seedEnd1-offset1,seedEnd2-offset2,trailingBP+1))) {
				break;
			}
			j1 = seedEnd1 + trailingBP;
			j2 = seedEnd2 + trailingBP;

			if (E_equal(curE,getHelixE(i1-offset1, i2-offset2, leadingBP+1)
							  + seedHandler->getSeedE(seedStart1, seedStart2)
							  + getHelixE(seedEnd1-offset1, seedEnd2-offset2, trailingBP+1)))
			{
				// Trace the first part if existing
				if (leadingBP != 0) {
					traceBackHelix(interaction, i1-offset1, i2-offset2, leadingBP+1);
					interaction.basePairs.push_back(energy.getBasePair(seedStart1, seedStart2));
				}
				// Trace the seed
				seedHandler->traceBackSeed(interaction,seedStart1,seedStart2);

				// Trace the last part if existing
				if (trailingBP != 0) {
					interaction.basePairs.push_back(energy.getBasePair(seedEnd1, seedEnd2));
					traceBackHelix(interaction, seedEnd1-offset1, seedEnd2-offset2, trailingBP+1);
				}
				traceNotFound = false;
			}


		} // trailing
	} // leading
	assert(!traceNotFound);

} // traceback

} // namespace