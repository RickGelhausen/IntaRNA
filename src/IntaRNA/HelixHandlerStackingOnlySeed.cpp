#include "IntaRNA/HelixHandlerStackingOnly.h"

namespace IntaRNA {

/////////////////////////////////////////////////////////////////////////////

size_t
HelixHandlerStackingOnly::
fillHelixSeed(const size_t i1min, const size_t i1max, const size_t i2min, const size_t i2max)
{

	// Fill Seed before computation
	if (seedHandler->fillSeed(i1min,i1max, i2min,i2max) == 0) {
		// No seed possible -> no helix with seed possible
		return 0;
	}

	helixSeed.resize( i1max-i1min+1, i2max-i2min+1 );

	// store index offset due to restricted matrix size generation
	offset1 = i1min;
	offset2 = i2min;

	// temporary variables
	size_t i1, i2, j1, j2, seedStart1, seedStart2, seedEnd1, seedEnd2, bestTrailingBP, possibleBasePairs;
	size_t  helixCountNotInf = 0, helixCount = 0;

	E_type leadingE, trailingE, bestTrailingE, totalEnergy;

	// fill for all start indices
	// in increasing index order
	for (i1=i1min; i1 < i1max+1; i1++ ) {
	for (i2=i1min; i2 < i2max+1; i2++ ) {
		// count possible helices
		helixCount++;

		// init according to no helix interaction
		helixSeed(i1 - offset1, i2 - offset2) = HelixMatrix::value_type(E_INF, 0);

		// skip non-complementary left helix boundaries
		if (!energy.areComplementary(i1, i2)) {
			continue; // go to next helixSeedE index
		}

		// Check if a seed can fit given the left boundaries
		// Note: Do it here to initialize all entries correctly
		if (std::min(helixSeed.size1()-i1, helixSeed.size2()-i2) < seedHandler->getConstraint().getBasePairs()) {
			continue;
		} else {
			// Seed fits, check how many bases are possible around
			possibleBasePairs = std::min(std::min(helixSeed.size1()-i1, helixSeed.size2()-i2), helixConstraint.getMaxBasePairs())-seedHandler->getConstraint().getBasePairs();
		}

		leadingE = 0.0;
		// TODO: Further tests for boundaries should be redundant with the strong condition beforehand
		// screen over all possible leading and trailing base pair combinations
		for (size_t leadingBP=0; leadingBP <= possibleBasePairs && i1 + leadingBP < helixSeed.size1()
																&& i2 + leadingBP < helixSeed.size2(); leadingBP++) {

			seedStart1 = i1 + leadingBP;
			seedStart2 = i2 + leadingBP;

			// Check whether seed is possible for this starting position
			// TODO: If no seed is possible from here, there should never be a possible seed anymore (so break should be alright)
			// TODO: Might need offset for this
			if (E_isINF(seedHandler->getSeedE(seedStart1, seedStart2))) {
				break;
			}

			// Update energy for the leading base pairs
			if (leadingBP > 0) {
				// TODO: THis might be handled by the seedCondition already
				if (!energy.areComplementary(seedStart1, seedStart2)) {
					break;
				}
				leadingE += energy.getE_interLeft(seedStart1 - 1, seedStart1, seedStart2 - 1, seedStart2);
			}

			seedEnd1 = seedStart1+seedHandler->getSeedLength1(seedStart1,seedStart2)-1;
			seedEnd2 = seedStart2+seedHandler->getSeedLength2(seedStart1,seedStart2)-1;

			trailingE = 0.0;
			bestTrailingE = 0.0;
			bestTrailingBP = 0;
			for (size_t trailingBP = 0; trailingBP <= possibleBasePairs - leadingBP
										&& seedEnd1+trailingBP < helixSeed.size1()
										&& seedEnd2+trailingBP < helixSeed.size2(); trailingBP++) {
				j1 = seedEnd1 + trailingBP;
				j2 = seedEnd2 + trailingBP;

				if (trailingBP > 0) {
					if (!energy.areComplementary(j1, j2)) {
						break;
					}
					trailingE += energy.getE_interLeft(j1 - 1, j1, j2 - 1, j2);
				}

				// only keep the best trailing energy, in order to calculate the energy in first loop
				// TODO: bestTrailingE might be redundant
				if (trailingE < bestTrailingE) {
					bestTrailingE = trailingE;
					bestTrailingBP = trailingBP;
				}
			}
			// Check whether this energy is the overall best so far
			// Done here to avoid problems when there are no trailingBP
			totalEnergy = leadingE + seedHandler->getSeedE(seedStart1,seedStart2) + bestTrailingE;
			if ( totalEnergy < helixSeed(i1,i2).first) {
				// Found new working helix with seed
				helixCountNotInf++;
				// Lengths
				size_t helixLength1 = leadingBP+seedHandler->getSeedLength1(seedStart1,seedStart2)+bestTrailingBP;
				size_t helixLength2 = leadingBP+seedHandler->getSeedLength2(seedStart1,seedStart2)+bestTrailingBP;
				// Creating new entry for helixSeed matrix
				helixSeed(i1,i2) = HelixMatrix::value_type(totalEnergy,
														   encodeHelixSeedLength(helixLength1,helixLength2));
			}

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
	     , seedStart2, seedEnd2
	     , bestTrailingBP
	     , j1, j2;

	bool traceNotFound = true;

	E_type curE = getHelixSeedE(i1_,i2_);

	E_type leadingE, trailingE, bestTrailingE, totalEnergy;

	// Calculate how many base pairs are possible allongside the seed.
	size_t possibleBasePairs = std::min(std::min(helixSeed.size1()-i1, helixSeed.size2()-i2), helixConstraint.getMaxBasePairs())-seedHandler->getConstraint().getBasePairs();


	leadingE = 0.0;
	// screen over all possible leading and trailing base pair combinations
	for (size_t leadingBP=0; traceNotFound
							 && leadingBP <= possibleBasePairs
							 && i1+ leadingBP < helixSeed.size1()
							 && i2+ leadingBP < helixSeed.size2(); leadingBP++) {

		seedStart1 = i1 + leadingBP;
		seedStart2 = i2 + leadingBP;

		// Check whether seed is possible for this starting position
		// TODO: If no seed is possible from here, there should never be a possible seed anymore (so break should be alright)
		// TODO: Might need offset for this
		if (E_isINF(seedHandler->getSeedE(seedStart1, seedStart2))) {
			break;
		}

		// Update energy for the leading base pairs
		if (leadingBP > 0) {
			if (!energy.areComplementary(seedStart1, seedStart2)) {
				break;
			}
			leadingE += energy.getE_interLeft(seedStart1 - 1, seedStart1, seedStart2 - 1, seedStart2);
		}

		seedEnd1 = seedStart1+seedHandler->getSeedLength1(seedStart1,seedStart2)-1;
		seedEnd2 = seedStart2+seedHandler->getSeedLength2(seedStart1,seedStart2)-1;

		trailingE = 0.0;
		bestTrailingE = 0.0;
		bestTrailingBP = 0;
		for (size_t trailingBP = 0; trailingBP <= possibleBasePairs - leadingBP
									&& seedEnd1+trailingBP < helixSeed.size1()
									&& seedEnd2+trailingBP < helixSeed.size2(); trailingBP++) {


			j1 = seedEnd1 + trailingBP;
			j2 = seedEnd2 + trailingBP;

			if (trailingBP > 0) {
				if (!energy.areComplementary(j1, j2)) {
					break;
				}
				trailingE += energy.getE_interLeft(j1 - 1, j1, j2 - 1, j2);
			}

			// only keep the best trailing energy, in order to calculate the energy in first loop
			// TODO: bestTrailingE might be redundant
			if (trailingE < bestTrailingE) {
				bestTrailingE = trailingE;
				bestTrailingBP = trailingBP;
			}
		} // trailing
		// Check whether this energy is the overall best so far
		// Done here to avoid problems when there are no trailingBP
		totalEnergy = leadingE + seedHandler->getSeedE(seedStart1,seedStart2) + bestTrailingE;
		if (E_equal(totalEnergy, curE)) {
			// Add leading bases
			for (size_t l = 0; l < leadingBP; l++) {
				interaction.basePairs.push_back( energy.getBasePair(i1+l, i2+l) );
			}

			// Add seed base pairs
			seedHandler->traceBackSeed(interaction, seedStart1, seedStart2);

			// Add trailing base pairs
			for (size_t l = 0; l < bestTrailingBP; l++) {
				interaction.basePairs.push_back( energy.getBasePair(seedEnd1+l, seedEnd2+l));
			}
			// Finish traceback
			traceNotFound = false;
		}
	} // leading
	assert(!traceNotFound);

} // traceback

} // namespace