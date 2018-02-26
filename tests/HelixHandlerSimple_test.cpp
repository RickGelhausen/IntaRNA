
#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/RnaSequence.h"
#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/HelixConstraint.h"
#include "IntaRNA/HelixHandlerSimple.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/ReverseAccessibility.h"

using namespace IntaRNA;

TEST_CASE( "HelixHandlerSimple", "[HelixHandlerSimple]") {


	SECTION( "getter", "[HelixHandlerSimple]") {

		RnaSequence r1("r1","GGGGG");
		RnaSequence r2("r2","CCCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy( acc1, racc );

		HelixConstraint hC(2,10,0);
		HelixHandlerSimple hhS(energy, hC);

		REQUIRE( &hhS.getInteractionEnergy() == &energy);
		REQUIRE( &hhS.getConstraint() == &hC);

	}

	SECTION( "fillHelix", "[HelixHandlerSimple]") {

		RnaSequence r1("r1","GGGGG");
		RnaSequence r2("r2","CCCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy( acc1, racc );

		HelixConstraint hC(2,5,0);
		HelixHandlerSimple hhS(energy, hC);

		REQUIRE(hhS.fillHelix(0, energy.size1()-1, 0, energy.size2()-1) == 2);
	}
}