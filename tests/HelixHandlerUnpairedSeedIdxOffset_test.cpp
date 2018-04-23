
#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/RnaSequence.h"
#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/HelixConstraint.h"
#include "IntaRNA/HelixHandlerUnpaired.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/SeedHandlerMfe.h"

using namespace IntaRNA;

TEST_CASE( "HelixHandlerUnpairedSeed", "[HelixHandlerUnpaired]" ) {

	SECTION("HelixSeed: Case 1 - all complementary", "[HelixHandlerUnpaired]") {

		RnaSequence r1("r1", "GGGGG");
		RnaSequence r2("r2", "CCCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 2, 0, false);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
		SeedConstraint sC(3, 0, 0, 0, 0, AccessibilityDisabled::ED_UPPER_BOUND, IndexRangeList(""), IndexRangeList(""),
						  "");

		SeedHandlerMfe sH(energy, sC);
		HelixHandlerUnpaired hhU(energy, hC);

		sH.fillSeed(0, energy.size1() - 1, 0, energy.size2() - 1);
		hhU.setSeedHandler(sH);

		hhU.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		REQUIRE(hhU.fillHelixSeed(0, energy.size1() - 1, 0, energy.size2() - 1) == 9);

		// (0,0)
		REQUIRE(hhU.getHelixSeedE(0, 0) == -3);
		REQUIRE(hhU.getHelixSeedLength1(0, 0) == 4);
		REQUIRE(hhU.getHelixSeedLength2(0, 0) == 4);

		// (0,1)
		REQUIRE(hhU.getHelixSeedE(0, 1) == -3);
		REQUIRE(hhU.getHelixSeedLength1(0, 1) == 4);
		REQUIRE(hhU.getHelixSeedLength2(0, 1) == 4);

		// (0,2)
		REQUIRE(hhU.getHelixSeedE(0, 2) == -2);
		REQUIRE(hhU.getHelixSeedLength1(0, 2) == 3);
		REQUIRE(hhU.getHelixSeedLength2(0, 2) == 3);

		// (0,3)
		REQUIRE(hhU.getHelixSeedE(0, 3) == E_INF);
		REQUIRE(hhU.getHelixSeedLength1(0, 3) == 0);
		REQUIRE(hhU.getHelixSeedLength2(0, 3) == 0);

		// (1,1)
		REQUIRE(hhU.getHelixSeedE(1, 1) == -3);
		REQUIRE(hhU.getHelixSeedLength1(1, 1) == 4);
		REQUIRE(hhU.getHelixSeedLength2(1, 1) == 4);

		// (2,2)
		REQUIRE(hhU.getHelixSeedE(2, 2) == -2);
		REQUIRE(hhU.getHelixSeedLength1(2, 2) == 3);
		REQUIRE(hhU.getHelixSeedLength2(2, 2) == 3);

		// (4,4)
		REQUIRE(hhU.getHelixSeedE(4, 4) == E_INF);
		REQUIRE(hhU.getHelixSeedLength1(4, 4) == 0);
		REQUIRE(hhU.getHelixSeedLength2(4, 4) == 0);


		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Case (0,0)
		Interaction interaction(r1, r2);
		hhU.traceBackHelixSeed(interaction, 0, 0);

		REQUIRE(interaction.basePairs.size() == 2);

		// First / last base pair of helixSeed
		REQUIRE(interaction.basePairs.begin()->first == 1);
		REQUIRE(interaction.basePairs.begin()->second == 3);

		REQUIRE(interaction.basePairs.rbegin()->first == 2);
		REQUIRE(interaction.basePairs.rbegin()->second == 2);

		// Case (4,4)
		interaction.clear();
		hhU.traceBackHelixSeed(interaction, 4, 4);

		REQUIRE(interaction.basePairs.size() == 0);
	}
}