
#include "catch.hpp"

#undef NDEBUG
#define protected public

#include "IntaRNA/RnaSequence.h"
#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/HelixConstraint.h"
#include "IntaRNA/HelixHandlerStackingOnly.h"
#include "IntaRNA/HelixHandlerIdxOffset.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/SeedHandlerMfe.h"


using namespace IntaRNA;

TEST_CASE( "HelixHandlerIdxOffset for StackingOnlySeed", "[HelixHandlerIdxOffset]") {

	SECTION("HelixSeed with Offset: Case 1 - offset 1", "[HelixHandlerIdxOffset]") {

		RnaSequence r1("r1", "GGGGGG");
		RnaSequence r2("r2", "GCCCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2,4,0);

		// seedBP / seedMaxUP / seedTMaxUP / seedQMaxUP / seedMaxE / seedMaxED / seedTRange / seedQRange / seedTQ
		SeedConstraint sC(3,0,0,0,0
				, AccessibilityDisabled::ED_UPPER_BOUND
				, IndexRangeList("")
				, IndexRangeList("")
				, "");

		HelixHandler * hhS = new HelixHandlerStackingOnly(energy, hC);
        SeedHandler * sH = new SeedHandlerMfe(energy, sC);
		SeedHandlerIdxOffset sHIO(sH);

		sH->fillSeed(0,energy.size1()-1, 0,energy.size2()-1);
		hhS->setSeedHandler(sH);

		HelixHandlerIdxOffset hhIO(hhS);

		// Set the offsets
		sHIO.setOffset1(1);
		sHIO.setOffset2(1);
		hhIO.setOffset1(1);
		hhIO.setOffset2(1);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////   FILLHELIXSEED  //////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
//		REQUIRE(sHIO.fillSeed( 0,energy.size1()-sHIO.getOffset1()-1, 0,energy.size2()-sHIO.getOffset2()-1 ) > 0);
		REQUIRE(hhIO.fillHelixSeed( 0,energy.size1()-hhIO.getOffset1()-1, 0,energy.size2()-hhIO.getOffset2()-1 ) == 9);

//		// (0,0)
//		REQUIRE(hhIO.getHelixSeedE(0,0) == -3);
//		REQUIRE(hhIO.getHelixSeedLength1(0,0) == 4);
//		REQUIRE(hhIO.getHelixSeedLength2(0,0) == 4);
//
//		// (0,1)
//		REQUIRE(hhIO.getHelixSeedE(0,1) == -3);
//		REQUIRE(hhIO.getHelixSeedLength1(0,1) == 4);
//		REQUIRE(hhIO.getHelixSeedLength2(0,1) == 4);
//
//		// (0,2)
//		REQUIRE(hhIO.getHelixSeedE(0,2) == -2);
//		REQUIRE(hhIO.getHelixSeedLength1(0,2) == 3);
//		REQUIRE(hhIO.getHelixSeedLength2(0,2) == 3);
//
//		// (0,3)
//		REQUIRE(hhIO.getHelixSeedE(0,3) == E_INF);
//		REQUIRE(hhIO.getHelixSeedLength1(0,3) == 0);
//		REQUIRE(hhIO.getHelixSeedLength2(0,3) == 0);
//
//		// (1,1)
//		REQUIRE(hhIO.getHelixSeedE(1,1) == -3);
//		REQUIRE(hhIO.getHelixSeedLength1(1,1) == 4);
//		REQUIRE(hhIO.getHelixSeedLength2(1,1) == 4);
//
//		// (2,2)
//		REQUIRE(hhIO.getHelixSeedE(2,2) == -2);
//		REQUIRE(hhIO.getHelixSeedLength1(2,2) == 3);
//		REQUIRE(hhIO.getHelixSeedLength2(2,2) == 3);
//
//		// (4,4)
//		REQUIRE(hhIO.getHelixSeedE(4,4) == E_INF);
//		REQUIRE(hhIO.getHelixSeedLength1(4,4) == 0);
//		REQUIRE(hhIO.getHelixSeedLength2(4,4) == 0);


		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////   TRACEBACK   ///////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////

//		// Case (0,0)
//		Interaction interaction(r1,r2);
//		hhIO.traceBackHelixSeed(interaction,0,0);
//
//		REQUIRE(interaction.basePairs.size() == 2);
//
//		// First / last base pair of helixSeed
//		REQUIRE(interaction.basePairs.begin()->first == 1);
//		REQUIRE(interaction.basePairs.begin()->second == 3);
//
//		REQUIRE(interaction.basePairs.rbegin()->first == 2);
//		REQUIRE(interaction.basePairs.rbegin()->second == 2);
//
//		// Case (4,4)
//		interaction.clear();
//		hhIO.traceBackHelixSeed(interaction,4,4);
//
//		REQUIRE(interaction.basePairs.size() == 0);
	}

}