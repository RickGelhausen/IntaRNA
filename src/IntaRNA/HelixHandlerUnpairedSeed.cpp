#include "IntaRNA/HelixHandlerUnpaired.h"

namespace IntaRNA {

/////////////////////////////////////////////////////////////////////////////

size_t
HelixHandlerUnpaired::
fillHelixSeed(const size_t i1min, const size_t i1max, const size_t i2min, const size_t i2max)
{
//	LOG(DEBUG) << "";
//	LOG(DEBUG) << "FILLHELIX SEED";
	helixSeed.resize( i1max-i1min+1, i2max-i2min+1 );

	// store index offset due to restricted matrix size generation
	offset1 = i1min;
	offset2 = i2min;

	// temporary variables
	size_t i1, i2, seedStart1, seedStart2, seedEnd1, seedEnd2, u1L, u1T, u2L, u2T, bestL1, bestL2, possibleBasePairs;
	size_t  helixCountNotInf = 0, helixCount = 0;

	E_type curE, tmpE;

	// fill for all start indices
	// in decreasing index order
	for (i1=i1max+1; i1-- > i1min;) {
	for (i2=i2max+1; i2-- > i2min;) {

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
		bestL1 = 0;
		bestL2 = 0;

		// screen over all possible leading and trailing base pair combinations
		for (size_t leadingBP=0; leadingBP <= possibleBasePairs
								 && (i1+leadingBP-offset1) < helixSeed.size1()
								 && (i2+leadingBP-offset2) < helixSeed.size2(); leadingBP++) {

			// Check all leading unpaired bases combinations
			for (u1L = 0; u1L < getConstraint().getMaxUnpaired()+1 && (i1+leadingBP+u1L-offset1) < helixSeed.size1(); u1L++) {
			for (u2L = 0; u2L < getConstraint().getMaxUnpaired()+1-u1L && (i2+leadingBP+u2L-offset2) < helixSeed.size2(); u2L++) {

				// If leading base pairs exist and helixE = E_INF -> skip to the next leadingBP
				if (leadingBP != 0) {
					if (E_isINF(getHelixE(i1-offset1,i2-offset2,leadingBP+1,u1L, u2L))) {
						continue;
					}
				} else {
					// ensure no dangling unpaired bases are considered
					if (u1L != 0 || u2L != 0)
						break;
				}

				// the start positions for the seed
				seedStart1 = i1+leadingBP+u1L;
				seedStart2 = i2+leadingBP+u2L;

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

//					LOG(DEBUG) << "TrailingBP: " << trailingBP;
//					LOG(DEBUG) << "helixSizes: " << helixSeed.size1() << " " << helixSeed.size2();
//					LOG(DEBUG) << "seedEnd+trailingBP: " << seedEnd1+trailingBP << " " << seedEnd2+trailingBP;
					// check all trailing unpaired bases combinations
					for (u1T = 0; u1T < getConstraint().getMaxUnpaired()+1 - (u1L+u2L) && (seedEnd1+trailingBP+u1T-offset1) < helixSeed.size1(); u1T++) {
					for (u2T = 0; u2T < getConstraint().getMaxUnpaired()+1 - (u1L+u2L+u1T) && (seedEnd2+trailingBP+u2T-offset2) < helixSeed.size2(); u2T++) {

//						LOG(DEBUG) << "u1T, u2T: " << u1T << " "<< u2T;
						// If trailing base pairs exist and helixE = E_INF -> skip to the next leadingBP
						if (trailingBP != 0) {
							if (E_isINF(getHelixE(seedEnd1 - offset1, seedEnd2 - offset2, trailingBP + 1, u1T, u2T))) {
								continue;
							}
						} else {
							// ensure no dangling unpaired bases are considered
							if (u1T != 0 || u2T != 0)
								break;
						}
//						LOG(DEBUG) << "i1, i2: " << i1 << " " << i2;
//						LOG(DEBUG) << "leadingBP+1: " << leadingBP+1;
//						LOG(DEBUG) << "u1L, u2L: " << u1L << " " << u2L;
//						LOG(DEBUG) << "seedStart1, seedStart2: " << seedStart1 << " " << seedStart2;
//						LOG(DEBUG) << "seedEnd1, seedEnd2: " << seedEnd1 << " " << seedEnd2;
//						LOG(DEBUG) << "trailingBP+1: " << trailingBP+1;
//						LOG(DEBUG) << "u1T, u2T: " << u1T << " " << u2T;
//						LOG(DEBUG) << "Energie: " << getHelixE(i1 - offset1, i2 - offset2, leadingBP + 1, u1L, u2L)
//								   << " + "<< seedHandler->getSeedE(seedStart1, seedStart2)
//								   << " + "<< getHelixE(seedEnd1 - offset1, seedEnd2 - offset2, trailingBP + 1, u1T, u2T);

						tmpE = getHelixE(i1 - offset1, i2 - offset2, leadingBP + 1, u1L, u2L)
							   + seedHandler->getSeedE(seedStart1, seedStart2)
							   + getHelixE(seedEnd1 - offset1, seedEnd2 - offset2, trailingBP + 1, u1T, u2T);

						if (tmpE < curE) {
							curE = tmpE;
							bestL1 = leadingBP + u1L + seedHandler->getSeedLength1(seedStart1, seedStart2) + trailingBP + u1T;
							bestL2 = leadingBP + u2L + seedHandler->getSeedLength2(seedStart1, seedStart2) + trailingBP + u2T;
 						}
					} // u2T
					} // u1T
				} // trailingBP

			} // u2L
			} // u1L
		} // leadingBP

		helixSeed(i1-offset1, i2-offset2) = HelixSeedMatrix::value_type(curE, E_isINF(curE) ? 0: encodeHelixSeedLength(bestL1,bestL2));

		// Ensures that the helixCount is only increased for the mfe helix.
		if (E_isNotINF(curE)) {
//			LOG(DEBUG) << "i1, i2: " << i1-offset1 << " " << i2-offset2 << " " << curE;
			helixCountNotInf++;
		}
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
	size_t i1 = i1_
	, i2 = i2_
	, seedStart1, seedEnd1
	, seedStart2, seedEnd2
	, u1L, u2L, u1T, u2T;

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
		for (u1L = 0; u1L < getConstraint().getMaxUnpaired()+1 && traceNotFound && (i1+leadingBP+u1L-offset1) < helixSeed.size1(); u1L++) {
		for (u2L = 0; u2L < getConstraint().getMaxUnpaired()+1-u1L && traceNotFound && (i2+leadingBP+u2L-offset2) < helixSeed.size2(); u2L++) {
			// If leading base pairs exist and helixE = E_INF -> skip to the next leadingBP
			if (leadingBP != 0) {
				if (E_isINF(getHelixE(i1-offset1,i2-offset2,leadingBP+1,u1L, u2L))) {
					continue;
				}
			} else {
				// ensure no dangling unpaired bases are considered
				if (u1L != 0 || u2L != 0)
					break;
			}

			// start positions of the seed
			seedStart1 = i1 + leadingBP + u1L;
			seedStart2 = i2 + leadingBP + u2L;

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
										&& seedEnd1 + trailingBP - offset1 < helixSeed.size1()
										&& seedEnd2 + trailingBP - offset2 < helixSeed.size2(); trailingBP++) {

				for (u1T = 0; u1T < getConstraint().getMaxUnpaired()+1 - (u1L+u2L) && traceNotFound && (seedEnd1+trailingBP+u1T-offset1) < helixSeed.size1(); u1T++) {
				for (u2T = 0; u2T < getConstraint().getMaxUnpaired()+1 - (u1L + u2L + u1T) && traceNotFound && (seedEnd2 + trailingBP + u2T - offset2) < helixSeed.size2(); u2T++) {

					// If trailing base pairs exist and helixE = E_INF -> skip to the next leadingBP
					if (trailingBP != 0) {
						if (E_isINF(getHelixE(seedEnd1 - offset1, seedEnd2 - offset2, trailingBP + 1, u1T, u2T))) {
							continue;
						}
					} else {
						// ensure no dangling unpaired bases are considered
						if (u1T != 0 || u2T != 0)
							break;
					}

					if (E_equal(curE, getHelixE(i1 - offset1, i2 - offset2, leadingBP + 1, u1L, u2L)
									  + seedHandler->getSeedE(seedStart1, seedStart2)
									  + getHelixE(seedEnd1 - offset1, seedEnd2 - offset2, trailingBP + 1, u1T, u2T))) {
//						LOG(DEBUG) << "Found Trace!";
//						LOG(DEBUG) << "i1, i2: " << i1 << " " << i2;
//						LOG(DEBUG) << "LeadingBP: " << leadingBP;
//						LOG(DEBUG) << "u1L, u2L: " << u1L << " " << u2L;
//						LOG(DEBUG) << "seedStart1, seedStart2: " << seedStart1 << " " << seedStart2;
//						LOG(DEBUG) << "TrailingBP: " << trailingBP;
//						LOG(DEBUG) << "u1T, u2T: " << u1T << " " << u2T;
//						LOG(DEBUG) << "seedEnd1, seedEnd2: " << seedEnd1 << " " << seedEnd2;
						// Trace the first part if existing
						if (leadingBP != 0) {
							traceBackHelix(interaction, i1 - offset1, i2 - offset2, leadingBP + 1,u1L,u2L);
							interaction.basePairs.push_back(energy.getBasePair(seedStart1, seedStart2));
						}
						// Trace the seed
						seedHandler->traceBackSeed(interaction, seedStart1, seedStart2);

						// Trace the last part if existing
						if (trailingBP != 0) {
							interaction.basePairs.push_back(energy.getBasePair(seedEnd1, seedEnd2));
							traceBackHelix(interaction, seedEnd1 - offset1, seedEnd2 - offset2, trailingBP + 1,u1T,u2T);
						}
						traceNotFound = false;
					}
				} // u2T
				} // u1T
			} // trailing

		} // u2L
		} // u1L
	} // leading
	assert(!traceNotFound);

} // traceback

} // namespace