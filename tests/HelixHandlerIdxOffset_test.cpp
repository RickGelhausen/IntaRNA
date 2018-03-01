
#include "catch.hpp"

#undef NDEBUG

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

	SECTION("fillHelix", "[HelixHandlerIdxOffset]") {

		RnaSequence r1("r1", "GGGGG");
		RnaSequence r2("r2", "CCCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2,4,0);
		HelixHandlerIdxOffset hhIO(new HelixHandlerStackingOnly(energy, hC));

		// Set the offsets
		hhIO.setOffset1(1);
		hhIO.setOffset2(2);

		// possible helices
		//REQUIRE(hhIO.fillHelix(0, energy.size1()-1, 0, energy.size2()-1)==6);
	}
}