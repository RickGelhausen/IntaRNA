#include "IntaRNA/HelixHandlerUnpaired.h"

namespace IntaRNA {

/////////////////////////////////////////////////////////////////////////////

size_t
HelixHandlerUnpaired::
fillHelixSeed(const size_t i1min, const size_t i1max, const size_t i2min, const size_t i2max)
{
//	LOG(DEBUG) << " ";
//	LOG(DEBUG) << "Start FILLHELIXSEED!";

	helixSeed.resize( i1max-i1min+1, i2max-i2min+1 );

	// store index offset due to restricted matrix size generation
	offset1 = i1min;
	offset2 = i2min;

	// temporary variables
	size_t i1, i2, seedStart1, seedStart2, seedEnd1, seedEnd2, j1, j2, bestL1, bestL2, possibleBasePairs, curL1, curL2;
	size_t  helixCountNotInf = 0, helixCount = 0;

	E_type curE_withED, curE, bestE_withED, bestE;

	// fill for all start indices
	// in decreasing index order
	for (i1=i1max+1; i1-- > i1min;) {
	for (i2=i2max+1; i2-- > i2min;) {
//		LOG(DEBUG) << "i1, i2: " << i1 << " " << i2;

		// count possible helices
		helixCount++;

		// init according to no helix interaction
		helixSeed(i1 - offset1, i2 - offset2) = HelixSeedMatrix::value_type(E_INF, 0);

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
		curE_withED = E_INF;
		bestE = E_INF;
		bestE_withED = E_INF;
		bestL1 = 0;
		bestL2 = 0;

		// screen over all possible leading and trailing base pair combinations
		for (size_t leadingBP=0; leadingBP <= possibleBasePairs
								 && (i1+leadingBP-offset1) < helixSeed.size1()
								 && (i2+leadingBP-offset2) < helixSeed.size2(); leadingBP++) {

			// If leading base pairs exist and helixE = E_INF -> skip to the next leadingBP
			if (leadingBP != 0 && E_isINF(getHelixE(i1-offset1,i2-offset2,leadingBP+1))) {
				continue;
			}

			// the start positions for the seed
			if (leadingBP != 0) {
				seedStart1 = i1 + getHelixLength1(i1 - offset1, i2 - offset2, leadingBP + 1) - 1;
				seedStart2 = i2 + getHelixLength2(i1 - offset1, i2 - offset2, leadingBP + 1) - 1;
			} else {
				seedStart1 = i1;
				seedStart2 = i2;
			}

			// check whether the right boundaries are broken
			if (seedStart1-offset1 >= helixSeed.size1() || seedStart2-offset2 >= helixSeed.size2()) {
				continue;
			}

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
				if (trailingBP != 0 && E_isINF(getHelixE(seedEnd1 - offset1, seedEnd2 - offset2, trailingBP + 1))) {
					continue;
				}

				// right boundary
				if (trailingBP != 0) {
					j1 = seedEnd1+getHelixLength1(seedEnd1-offset1,seedEnd2-offset2,trailingBP+1)-1;
					j2 = seedEnd2+getHelixLength2(seedEnd1-offset1,seedEnd2-offset2,trailingBP+1)-1;
				} else {
					j1 = seedEnd1;
					j2 = seedEnd2;
				}
//				LOG_IF(i1==0 && i2==0, DEBUG) << "trailingBP: " << trailingBP+1;
//				LOG_IF(i1==0 && i2==0, DEBUG) << "j1, j2: " << j1 << " " << j2;
				// check whether the right boundaries are broken
				if (j1-offset1 >= helixSeed.size1()
					|| j2-offset2 >= helixSeed.size2()) {
					continue;
				}

				// ensure that ED-values are within the boundaries (default 999)
				if (energy.getED1(i1, j1) > helixConstraint.getMaxED()
					|| energy.getED2(i2, j2) > helixConstraint.getMaxED()) {
					continue;
				}

				// current lengths
				curL1 = seedHandler->getSeedLength1(seedStart1, seedStart2);
				curL2 = seedHandler->getSeedLength2(seedStart1, seedStart2);

				// Add leadingBP length contribution
				if (leadingBP != 0) {
					curL1 += getHelixLength1(i1-offset1, i2-offset2, leadingBP+1)-1;
					curL2 += getHelixLength2(i1-offset1, i2-offset2, leadingBP+1)-1;
				}
				// Add trailingBP length contribution
				if (trailingBP != 0) {
					curL1 += getHelixLength1(seedEnd1-offset1, seedEnd2-offset2, trailingBP+1)-1;
					curL2 += getHelixLength2(seedEnd1-offset1, seedEnd2-offset2, trailingBP+1)-1;
				}

				// check the lengths
				if (bestL1 != curL1 || bestL2 != curL2) {
					// energy without contributions
					curE = getHelixE(i1 - offset1, i2 - offset2, leadingBP + 1)
						   + seedHandler->getSeedE(seedStart1, seedStart2)
						   + getHelixE(seedEnd1 - offset1, seedEnd2 - offset2, trailingBP + 1);
					// energy value
					curE_withED = energy.getE(i1,j1,i2,j2, curE) + energy.getE_init();

					// If no ED-values are wanted, remove them
					if (helixConstraint.noED())
						curE_withED -= (energy.getED1(i1,j1) + energy.getED2(i2, j2));

		//				LOG(DEBUG) << "curE: " << curE;
					if (curE_withED < bestE_withED) {
		//					LOG(DEBUG) << "leading/trailing: " << leadingBP << " " << trailingBP << " with bestE: " << curE_withED << "(" << curE << ") beating last bestE: " << bestE_withED;
						bestE_withED = curE_withED;
						bestE = curE;
						bestL1 = curL1;
						bestL2 = curL2;
		//					LOG(DEBUG) << "lengths: " << bestL1 << " " << bestL2;
					}
				}
			} // trailingBP
		} // leadingBP

		// Ensures that the helixCount is only increased for the mfe helix.
		if (E_isNotINF(bestE_withED)) {
			// overwrite all helices with too high energy -> infeasible start interactions
			if (bestE_withED > helixConstraint.getMaxE()) {
				bestE = E_INF;
			} else {
//				// energy without dangling ends and other contributions
//				bestE -= energy.getE_init();
//				if (!helixConstraint.noED()) {
//					bestE-= (energy.getED1(i1, i1+bestL1-1) + energy.getED2(i2, i2+bestL2-1));
//				}
				// count true helix
				helixCountNotInf++;
			}
		}
//		LOG(DEBUG) << "(i1,i2): BestL, bestE: (" << i1 <<", " << i2 << ") " << bestL1 << " " << bestL2 << " E: " << bestE;
		helixSeed(i1-offset1, i2-offset2) = HelixSeedMatrix::value_type(bestE, E_isINF(bestE) ? 0: encodeHelixSeedLength(bestL1,bestL2));


	} // i2
	} // i1

	return helixCountNotInf;
}

//////////////////////////////////////////////////////////////////////////

void
HelixHandlerUnpaired::
traceBackHelixSeed( Interaction & interaction
		, const size_t i1_
		, const size_t i2_)
{
//	LOG(DEBUG) << "Traceback Start!";
	size_t i1 = i1_
	, i2 = i2_
	, seedStart1, seedEnd1
	, seedStart2, seedEnd2
	, j1, j2;

	bool traceNotFound = true;

	E_type curE = getHelixSeedE(i1_,i2_);

	// No traceback possible for current boundary
	if (E_isINF(curE)) {
		return;
	}

	// TODO: Check if this work when seed allows unpaired bases
	// Calculate how many base pairs are possible allongside the seed.
	// Note: If seedHandler allows unpaired positions this check is not enough, check happens in loop
	size_t possibleBasePairs = std::min(std::min(helixSeed.size1()-i1 +offset1, helixSeed.size2()-i2+offset2), helixConstraint.getMaxBasePairs())-seedHandler->getConstraint().getBasePairs();
//	LOG(DEBUG) << "possibleBasePairs: " << possibleBasePairs;


//	LOG(DEBUG) << "curE: "<< curE;
	// screen over all possible leading and trailing base pair combinations
	for (size_t leadingBP=0; traceNotFound
							 && leadingBP <= possibleBasePairs
							 && i1 + leadingBP-offset1 < helixSeed.size1()
							 && i2 + leadingBP-offset2 < helixSeed.size2(); leadingBP++) {

		// If leading base pairs exist and helixE = E_INF -> skip to the next leadingBP
		if (leadingBP != 0 && E_isINF(getHelixE(i1-offset1,i2-offset2,leadingBP+1))) {
			continue;
		}

		// the start positions for the seed
		if (leadingBP != 0) {
			seedStart1 = i1 + getHelixLength1(i1 - offset1, i2 - offset2, leadingBP + 1) - 1;
			seedStart2 = i2 + getHelixLength2(i1 - offset1, i2 - offset2, leadingBP + 1) - 1;
		} else {
			seedStart1 = i1;
			seedStart2 = i2;
		}
		// check whether the right boundaries are broken
		if (seedStart1-offset1 >= helixSeed.size1() || seedStart2-offset2 >= helixSeed.size2()) {
			continue;
		}

		// Check whether seed is possible for this starting position
		if (E_isINF(seedHandler->getSeedE(seedStart1, seedStart2))) {
			continue;
		}

		// end positions of the seed
		seedEnd1 = seedStart1 + seedHandler->getSeedLength1(seedStart1, seedStart2) - 1;
		seedEnd2 = seedStart2 + seedHandler->getSeedLength2(seedStart1, seedStart2) - 1;

		// Trailing base pairs
		for (size_t trailingBP = 0; traceNotFound
									&& trailingBP <= possibleBasePairs - leadingBP
									&& (seedEnd1 + trailingBP - offset1) < helixSeed.size1()
									&& (seedEnd2 + trailingBP - offset2) < helixSeed.size2(); trailingBP++) {


			// If trailing base pairs exist and helixE = E_INF -> skip to the next leadingBP
			if (trailingBP != 0 && E_isINF(getHelixE(seedEnd1 - offset1, seedEnd2 - offset2, trailingBP + 1))) {
				continue;
			}

			// right boundary
			if (trailingBP != 0) {
				j1 = seedEnd1+getHelixLength1(seedEnd1-offset1,seedEnd2-offset2,trailingBP+1)-1;
				j2 = seedEnd2+getHelixLength2(seedEnd1-offset1,seedEnd2-offset2,trailingBP+1)-1;
			} else {
				j1 = seedEnd1;
				j2 = seedEnd2;
			}

			// check whether the right boundaries are broken
			if (j1-offset1 >= helixSeed.size1()
				|| j2-offset2 >= helixSeed.size2()) {
				continue;
			}

//			LOG_IF(leadingBP==0 && trailingBP==2,DEBUG) << getHelixE(i1 - offset1, i2 - offset2, leadingBP + 1)
//						  + seedHandler->getSeedE(seedStart1, seedStart2)
//						  + getHelixE(seedEnd1 - offset1, seedEnd2 - offset2, trailingBP + 1);

			if (E_equal(curE, getHelixE(i1 - offset1, i2 - offset2, leadingBP + 1)
							  + seedHandler->getSeedE(seedStart1, seedStart2)
							  + getHelixE(seedEnd1 - offset1, seedEnd2 - offset2, trailingBP + 1))) {
//				LOG(DEBUG) << "Trace Found!";
//				LOG(DEBUG) << "i1, i2: " << i1 << " " << i2;
//				LOG(DEBUG) << "leadinBP+1: " << leadingBP+1;
//				LOG(DEBUG) << "seedStarts: " << seedStart1 << " " << seedStart2;
//				LOG(DEBUG) << "seedEnds: " << seedEnd1 << " " << seedEnd2;
//				LOG(DEBUG) << "trailingBP+1: " << trailingBP+1;
//
//				LOG(DEBUG) << "Energy: " << getHelixE(i1 - offset1, i2 - offset2, leadingBP + 1)
//						   << " + " << seedHandler->getSeedE(seedStart1, seedStart2)
//						   << " + " << getHelixE(seedEnd1 - offset1, seedEnd2 - offset2, trailingBP + 1);
				// Trace the first part if existing
				if (leadingBP != 0) {
					traceBackHelix(interaction, i1 - offset1, i2 - offset2, leadingBP + 1);
					interaction.basePairs.push_back(energy.getBasePair(seedStart1, seedStart2));
				}
//				interaction.clear();
				// Trace the seed
				seedHandler->traceBackSeed(interaction, seedStart1, seedStart2);
//				LOG(DEBUG) << interaction;

				// Trace the last part if existing
				if (trailingBP != 0) {
					interaction.basePairs.push_back(energy.getBasePair(seedEnd1, seedEnd2));
//					LOG(DEBUG) << interaction;
					traceBackHelix(interaction, seedEnd1 - offset1, seedEnd2 - offset2, trailingBP + 1);
//					LOG(DEBUG) << interaction;
				}
				traceNotFound = false;
			}
		} // trailing
	} // leading
	assert(!traceNotFound);

} // traceback

} // namespace