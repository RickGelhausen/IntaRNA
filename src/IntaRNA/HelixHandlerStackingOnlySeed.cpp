#include "IntaRNA/HelixHandlerStackingOnly.h"

namespace IntaRNA {

/////////////////////////////////////////////////////////////////////////////

size_t
HelixHandlerStackingOnly::
fillHelixSeed(const size_t i1min, const size_t i1max, const size_t i2min, const size_t i2max)
{
//	LOG(DEBUG) << "";
//	LOG(DEBUG) << "FILLHELIXSEED!";
	helixSeed.resize( i1max-i1min+1, i2max-i2min+1 );

//	LOG(DEBUG) << "i1min, i1max, i2min, i2max: " << i1min << " " << i1max << " " << i2min << " " << i2max;
	// store index offset due to restricted matrix size generation
	offset1 = i1min;
	offset2 = i2min;

	// temporary variables
	size_t i1, i2, seedStart1, seedStart2, seedEnd1, seedEnd2, bestL1, bestL2, possibleBasePairs;
	size_t  helixCountNotInf = 0, helixCount = 0;

	E_type curE, tmpE;

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

//		LOG(DEBUG) << "i1, i2, possibleBasePairs: " << i1 << " " << i2 << " " << possibleBasePairs;
		// Initialuze variables
		curE = E_INF;
		bestL1 = 0;
		bestL2 = 0;

		// screen over all possible leading and trailing base pair combinations
		for (size_t leadingBP=0; leadingBP <= possibleBasePairs && (i1+leadingBP-offset1) < helixSeed.size1()
								 								&& (i2+leadingBP-offset2) < helixSeed.size2(); leadingBP++) {
//			LOG(DEBUG) << "leadingBP: " << leadingBP;
			// TODO: THESE CONDITIONS SHOULD CONTAIN CONTINUES FOR UNPAIRED VARIANT
			// If leading base pairs exist and helixE = E_INF -> skip to the next leadingBP
//			LOG(DEBUG) << "HelixE LEAD i1, i2, leadingBP+1: " << i1-offset1 << " " << i2-offset2 << " " << leadingBP+1;
//			LOG(DEBUG) << "HelixE LEAD: " << getHelixE(i1-offset1,i2-offset2,leadingBP+1);
			if (leadingBP != 0)
				if (E_isINF(getHelixE(i1-offset1,i2-offset2,leadingBP+1))) {
					continue;
				}

			// the start positions for the seed
			seedStart1 = i1+leadingBP;
			seedStart2 = i2+leadingBP;

//			LOG(DEBUG) << "SeedStarts: " << seedStart1 << " " << seedStart2;
			// If no seed is possible here, skip to next leading base pair number
			if (E_isINF(seedHandler->getSeedE(seedStart1, seedStart2))) {
				continue;
			}
//			LOG(DEBUG) << "SeedE: " << seedHandler->getSeedE(seedStart1, seedStart2);
			// the end positions of the seed
			seedEnd1 = seedStart1+seedHandler->getSeedLength1(seedStart1, seedStart2)-1;
			seedEnd2 = seedStart2+seedHandler->getSeedLength2(seedStart1, seedStart2)-1;

//			LOG(DEBUG) << "SeedEnds: " << seedEnd1 << " " << seedEnd2;
			// Run over all trailing base pairs
			for (size_t trailingBP=0; trailingBP <= possibleBasePairs - leadingBP
									  && (seedEnd1+trailingBP-offset1) < helixSeed.size1()
									  && (seedEnd2+trailingBP-offset2) < helixSeed.size2(); trailingBP++) {

//				LOG(DEBUG) << "trailingBP " << trailingBP;
//				LOG(DEBUG) << "HelixE Trail seedEnd1, seedEnd2, trailingBP+1: " << seedEnd1-offset1 << " " << seedEnd2-offset2 << " " << trailingBP+1;
//				LOG(DEBUG) << "HelixE Trail: " << getHelixE(seedEnd1-offset1, seedEnd2-offset2, trailingBP+1);
				// If trailing base pairs exist and helixE = E_INF -> skip to the next leadingBP
				if (trailingBP != 0)
					if (E_isINF(getHelixE(seedEnd1-offset1,seedEnd2-offset2,trailingBP+1))) {
						break;
					}

				tmpE = getHelixE(i1-offset1,i2-offset2,leadingBP+1) + seedHandler->getSeedE(seedStart1, seedStart2) + getHelixE(seedEnd1-offset1,seedEnd2-offset2,trailingBP+1);
				if (tmpE < curE) {
					curE = tmpE;
					bestL1 = leadingBP + seedHandler->getSeedLength1(seedStart1,seedStart2) + trailingBP;
					bestL2 = leadingBP + seedHandler->getSeedLength2(seedStart1,seedStart2) + trailingBP;
				}
//				LOG(DEBUG) << "-";
			} // trailingBP
		} // leadingBP

//		LOG(DEBUG) << "i1, i2, curE, bestL1, bestL2: " << i1 << " " << i2 << " " << curE << " " << bestL1 << " " << bestL2;
//		LOG(DEBUG) << "OFFSET: i1, i2 "<< i1-offset1 << " " << i2-offset2;
		helixSeed(i1-offset1, i2-offset2) = HelixMatrix::value_type(curE, E_isINF(curE) ? 0: encodeHelixSeedLength(bestL1,bestL2));
//		LOG(DEBUG) << "-----------------------------------------------------------------------------------------------------------";
		// Ensures that the helixCount is only increased for the mfe helix.
		if (E_isNotINF(curE)) {
			helixCountNotInf++;
		}
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
	     , seedStart2, seedEnd2;

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

	// screen over all possible leading and trailing base pair combinations
	for (size_t leadingBP=0; traceNotFound
							 && leadingBP <= possibleBasePairs
							 && i1 + leadingBP-offset1 < helixSeed.size1()
							 && i2 + leadingBP-offset2 < helixSeed.size2(); leadingBP++) {

		// If leading base pairs exist and helixE = E_INF -> skip to the next leadingBP
		if (leadingBP != 0)
			if (E_isINF(getHelixE(i1-offset1,i2-offset2,leadingBP+1))) {
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
		for (size_t trailingBP = 0; trailingBP <= possibleBasePairs - leadingBP
									&& seedEnd1+trailingBP-offset1 < helixSeed.size1()
									&& seedEnd2+trailingBP-offset2 < helixSeed.size2(); trailingBP++) {

			// If trailing base pairs exist and helixE = E_INF -> skip to the next leadingBP
			if (trailingBP != 0)
				if (E_isINF(getHelixE(seedEnd1-offset1,seedEnd2-offset2,trailingBP+1))) {
					break;
				}

			if (E_equal(curE, getHelixE(i1-offset1, i2-offset2, leadingBP+1)
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