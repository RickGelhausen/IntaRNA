
#include <IntaRNA/AccessibilityDisabled.h>
#include <IntaRNA/InteractionEnergyBasePair.h>
#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/HelixConstraint.h"
#include "IntaRNA/HelixHandlerSimple.h"

using namespace IntaRNA;

TEST_CASE( "HelixHandlerSimple", "[HelixHandlerSimple]") {

	RnaSequence r1("r1","AAAAAAAA");
	RnaSequence r2("r2","UUUUUUUU");
	AccessibilityDisabled acc1(r1, 0, NULL);
	AccessibilityDisabled acc2(r2, 0, NULL);
	ReverseAccessibility racc(acc2);
	InteractionEnergyBasePair energy( acc1, racc );

	SECTION( "getter", "[HelixHandlerSimple]") {

		HelixConstraint hC(2,10,0);
		HelixHandlerSimple hhS(energy, hC);



	}
}