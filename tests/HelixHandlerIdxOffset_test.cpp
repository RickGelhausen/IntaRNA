
#include "catch.hpp"

#undef NDEBUG
#define protected public

#include "IntaRNA/RnaSequence.h"
#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/HelixConstraint.h"
#include "IntaRNA/HelixHandlerStackingOnly.h"
#include "IntaRNA/HelixHandlerIdxOffset.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/Interaction.h"
#include "IntaRNA/ReverseAccessibility.h"


using namespace IntaRNA;

TEST_CASE( "HelixHandlerIdxOffset", "[HelixHandlerIdxOffset]") {

	SECTION("getter", "[HelixHandlerIdxOffset]") {

		RnaSequence r1("r1", "GGGGA");
		RnaSequence r2("r2", "ACCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 10, 0);
		HelixHandlerIdxOffset hhIO(new HelixHandlerStackingOnly(energy, hC));

		// Initial offset of 0
		REQUIRE(hhIO.getOffset1() == 0);
		REQUIRE(hhIO.getOffset2() == 0);

		// Set offset
		hhIO.setOffset1(1);
		hhIO.setOffset2(4);
		REQUIRE(hhIO.getOffset1() == 1);
		REQUIRE(hhIO.getOffset2() == 4);

		// get Constraints
		REQUIRE(hhIO.getConstraint().getMinBasePairs() == 2);
		REQUIRE(hhIO.getConstraint().getMaxBasePairs() == 10);
		REQUIRE(hhIO.getConstraint().getMaxUnpaired() == 0);

	}

	SECTION("fillHelix StackOnly offset 1", "[HelixHandlerIdxOffset]") {

		RnaSequence r1("r1", "GGGAAGG");
		RnaSequence r2("r2", "CCAACCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2,4,0);
		HelixHandlerIdxOffset hhIO(new HelixHandlerStackingOnly(energy, hC));

		// Set the offsets
		hhIO.setOffset1(1);
		hhIO.setOffset2(1);

		// possible helices
		REQUIRE(hhIO.fillHelix(0, energy.size1()-hhIO.getOffset1()-1, 0, energy.size2()-hhIO.getOffset2()-1)==4);

		// All Working
		// (0,0)
		REQUIRE(hhIO.getHelixE(0, 0) == -1);
		REQUIRE(hhIO.getHelixLength1(0, 0) == 2);
		REQUIRE(hhIO.getHelixLength2(0, 0) == hhIO.getHelixLength1(0, 0));

		// (0,4)
		REQUIRE(hhIO.getHelixE(0, 4) == -1);
		REQUIRE(hhIO.getHelixLength1(0, 4) == 2);
		REQUIRE(hhIO.getHelixLength2(0, 4) == hhIO.getHelixLength1(0, 4));

		// (0,4)
		REQUIRE(hhIO.getHelixE(4, 0) == -1);
		REQUIRE(hhIO.getHelixLength1(4, 0) == 2);
		REQUIRE(hhIO.getHelixLength2(4, 0) == hhIO.getHelixLength1(4, 0));

		// (0,4)
		REQUIRE(hhIO.getHelixE(4, 4) == -1);
		REQUIRE(hhIO.getHelixLength1(4, 4) == 2);
		REQUIRE(hhIO.getHelixLength2(4, 4) == hhIO.getHelixLength1(4, 4));

		// Not Working
		// (3,4)
		REQUIRE(hhIO.getHelixE(3, 4) == E_INF);
		REQUIRE(hhIO.getHelixLength1(3, 4) == 0);
		REQUIRE(hhIO.getHelixLength2(3, 4) == hhIO.getHelixLength1(3, 4));

		// (1,1)
		REQUIRE(hhIO.getHelixE(1, 1) == E_INF);
		REQUIRE(hhIO.getHelixLength1(1, 1) == 0);
		REQUIRE(hhIO.getHelixLength2(1, 1) == hhIO.getHelixLength1(1, 1));
	}

	SECTION("TracebackHelix StackOnly offset1", "[HelixHandlerIdxOffset]") {

		RnaSequence r1("r1", "GGGAAGG");
		RnaSequence r2("r2", "CCAACCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 0);
		HelixHandlerIdxOffset hhIO(new HelixHandlerStackingOnly(energy, hC));

		// Set the offsets
		hhIO.setOffset1(1);
		hhIO.setOffset2(1);

		// possible helices
		REQUIRE(hhIO.fillHelix(0, energy.size1() - hhIO.getOffset1() - 1, 0, energy.size2() - hhIO.getOffset2() - 1) == 4);

		// Case (0,0)
		Interaction interaction1(r1,r2);
		hhIO.traceBackHelix(interaction1,0,0);

		REQUIRE(interaction1.basePairs.size() == 1);
		REQUIRE(interaction1.basePairs.begin()->first == 3);
		REQUIRE(interaction1.basePairs.begin()->second == 3);

		REQUIRE(interaction1.basePairs.rbegin()->first == 3);
		REQUIRE(interaction1.basePairs.rbegin()->second == 3);
	}
}