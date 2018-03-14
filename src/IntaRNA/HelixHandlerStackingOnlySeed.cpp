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
	size_t i1, i2, curBP, j1, j2, seedStart1, seedStart2, seedEnd1, seedEnd2, bestTrailingBP;
	size_t  helixCountNotInf = 0, helixCount = 0;

	E_type leadingE, trailingE, bestTrailingE, totalEnergy;

	// Calculate how many base pairs are possible allongside the seed.
	size_t possibleBasePairs = getConstraint().getMaxBasePairs()-seedHandler.getConstraint().getBasePairs();

	// fill for all start indeices
	// in decreasing index order
	for (i1=0; i1 < i1max; i1++ ) {
	for (i2=0; i2 < i2max; i2++ ) {

		// count possible helices
		helixCount++;

		// init according to no helix interaction
		helixSeed(i1 - offset1, i2 - offset2) = HelixMatrix::value_type(E_INF, 0);

		// skip non-complementary left helix boundaries
		if (!energy.areComplementary(i1, i2)) {
			continue; // go to next helixSeedE index
		}

		leadingE = 0.0;
		 // screen over all possible leading and trailing base pair combinations
		for (size_t leadingBP=0; leadingBP <= possibleBasePairs; leadingBP++) {

			seedStart1 = i1 + leadingBP;
			seedStart2 = i2 + leadingBP;

			// Check whether seed is possible for this starting position
			// TODO: If no seed is possible from here, there should never be a possible seed anymore (so break should be alright)
			if (E_isINF(seedHandler.getSeedE(seedStart1, seedStart2))) {
				break;
			}

			// Update energy for the leading base pairs
			if (leadingBP > 0)
				leadingE += energy.getE_interLeft(seedStart1-1,seedStart1, seedStart2-1,seedStart2);


			seedEnd1 = seedStart1+seedHandler.getSeedLength1(seedStart1,seedStart2);
			seedEnd2 = seedStart2+seedHandler.getSeedLength2(seedStart1,seedStart2);

			trailingE = 0.0;
			bestTrailingE = 0.0;
			bestTrailingBP = 0;
			for (size_t trailingBP = 0; trailingBP <= possibleBasePairs - leadingBP; trailingBP++) {

				j1 = seedEnd1 + trailingBP;
				j2 = seedEnd2 + trailingBP;

				if (trailingBP > 0);
					trailingE += energy.getE_interLeft(j1-1,j1,j2-1,j2);

				// only keep the best trailing energy, in order to calculate the energy in first loop
				// TODO: bestTrailingE might be redundant
				if (trailingE < bestTrailingE) {
					bestTrailingE = trailingE;
					bestTrailingBP = trailingBP;
				}
			}
			// Check whether this energy is the overall best so far
			// Done here to avoid problems when there are no trailingBP
			totalEnergy = leadingE + seedHandler.getSeedE(seedStart1,seedStart2) + bestTrailingE;
			if ( totalEnergy < helixSeed(i1,i2).first) {
				size_t helixLength1 = leadingBP+seedHandler.getSeedLength1(seedStart1,seedStart2)+bestTrailingBP;
				size_t helixLength2 = leadingBP+seedHandler.getSeedLength2(seedStart1,seedStart2)+bestTrailingBP;
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

}

} // namespace