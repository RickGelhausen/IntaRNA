// TODO: REMOVE this for real test
#define protected public

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

		HelixConstraint hC(2,4,0);
		HelixHandlerSimple hhS(energy, hC);


		// When counting all non-inf values for helices -> 29
		// When only counting the optimal helices that are non-inf -> 16
		// (0,0) -> 4 ; (1,0) -> 4 ; (2,0) -> 3 ; (3,0) -> 2
		// (0,1) -> 4 ; (1,1) -> 4 ; (2,1) -> 3 ; (3,1) -> 2
		// (0,2) -> 3 ; (1,2) -> 3 ; (2,2) -> 3 ; (3,2) -> 2
		// (0,3) -> 2 ; (1,3) -> 2 ; (2,3) -> 2 ; (3,3) -> 2
		REQUIRE(hhS.fillHelix(0, energy.size1()-1, 0, energy.size2()-1) == 16);


		std::cout << "Hello" << std::endl;
		std::cout << energy.size1() << std::endl;
		std::cout << energy.size2() << std::endl;
		for (size_t i1=0; i1 < energy.size1(); i1++){
		for (size_t i2=0; i2 < energy.size2(); i2++){
			std::cout << "(i1,i2): " << "(" << i1 << "," << i2 << ")" << std::endl;
			for (size_t bp=2; bp < 5; bp++) {
				std::cout << "-------------------------" << std::endl;
				std::cout << "     bpE  " << bp << " : " << hhS.getHelixE(i1,i2,bp) << std::endl;
			}
			std::cout << "" << std::endl;
			std::cout << "     Optima bpL  " << " : " << hhS.getHelixLength1(i1,i2) << std::endl;
			std::cout << "" << std::endl;
		}
		}
		REQUIRE(hhS.getHelixLength1(2,2) == 3);
	}
}
