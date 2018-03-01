
#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/RnaSequence.h"
#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/HelixConstraint.h"
#include "IntaRNA/HelixHandlerStackingOnly.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/Interaction.h"
#include "IntaRNA/ReverseAccessibility.h"


using namespace IntaRNA;

TEST_CASE( "HelixHandlerStackingOnly", "[HelixHandlerStackingOnly]") {


	SECTION("getter", "[HelixHandlerStackingOnly]") {

		RnaSequence r1("r1", "GGGGG");
		RnaSequence r2("r2", "CCCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 10, 0);
		HelixHandlerStackingOnly hhS(energy, hC);

		REQUIRE(&hhS.getInteractionEnergy() == &energy);
		REQUIRE(&hhS.getConstraint() == &hC);

	}

	SECTION("fillHelix 1: Everything is complementary", "[HelixHandlerStackingOnly]") {

		// Case 1 Perfect sequence
		RnaSequence r1("r1", "GGGGG");
		RnaSequence r2("r2", "CCCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 0);
		HelixHandlerStackingOnly hhS(energy, hC);


		// When counting all non-inf values for helices -> 29
		// When only counting the optimal helices that are non-inf -> 16
		// (0,0) -> 4 ; (1,0) -> 4 ; (2,0) -> 3 ; (3,0) -> 2
		// (0,1) -> 4 ; (1,1) -> 4 ; (2,1) -> 3 ; (3,1) -> 2
		// (0,2) -> 3 ; (1,2) -> 3 ; (2,2) -> 3 ; (3,2) -> 2
		// (0,3) -> 2 ; (1,3) -> 2 ; (2,3) -> 2 ; (3,3) -> 2
		REQUIRE(hhS.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1) == 16);


		// All optimal combinations
		// (0,0)
		REQUIRE(hhS.getHelixE(0, 0) == -3);
		REQUIRE(hhS.getHelixLength1(0, 0) == 4);
		REQUIRE(hhS.getHelixLength2(0, 0) == hhS.getHelixLength1(0, 0));

		// (0,1)
		REQUIRE(hhS.getHelixE(0, 1) == -3);
		REQUIRE(hhS.getHelixLength1(0, 1) == 4);
		REQUIRE(hhS.getHelixLength2(0, 1) == hhS.getHelixLength1(0, 1));

		// (0,2)
		REQUIRE(hhS.getHelixE(0, 2) == -2);
		REQUIRE(hhS.getHelixLength1(0, 2) == 3);
		REQUIRE(hhS.getHelixLength2(0, 2) == hhS.getHelixLength1(0, 2));

		// (0,3)
		REQUIRE(hhS.getHelixE(0, 3) == -1);
		REQUIRE(hhS.getHelixLength1(0, 3) == 2);
		REQUIRE(hhS.getHelixLength2(0, 3) == hhS.getHelixLength1(0, 3));

		// (1,0)
		REQUIRE(hhS.getHelixE(1, 0) == -3);
		REQUIRE(hhS.getHelixLength1(1, 0) == 4);
		REQUIRE(hhS.getHelixLength2(1, 0) == hhS.getHelixLength1(1, 0));

		// (1,1)
		REQUIRE(hhS.getHelixE(1, 1) == -3);
		REQUIRE(hhS.getHelixLength1(1, 1) == 4);
		REQUIRE(hhS.getHelixLength2(1, 1) == hhS.getHelixLength1(1, 1));

		// (1,2)
		REQUIRE(hhS.getHelixE(1, 2) == -2);
		REQUIRE(hhS.getHelixLength1(1, 2) == 3);
		REQUIRE(hhS.getHelixLength2(1, 2) == hhS.getHelixLength1(1, 2));

		// (1,3)
		REQUIRE(hhS.getHelixE(1, 3) == -1);
		REQUIRE(hhS.getHelixLength1(1, 3) == 2);
		REQUIRE(hhS.getHelixLength2(1, 3) == hhS.getHelixLength1(1, 3));

		// (2,0)
		REQUIRE(hhS.getHelixE(2, 0) == -2);
		REQUIRE(hhS.getHelixLength1(2, 0) == 3);
		REQUIRE(hhS.getHelixLength2(2, 0) == hhS.getHelixLength1(2, 0));

		// (2,1)
		REQUIRE(hhS.getHelixE(2, 1) == -2);
		REQUIRE(hhS.getHelixLength1(2, 1) == 3);
		REQUIRE(hhS.getHelixLength2(2, 1) == hhS.getHelixLength1(2, 1));

		// (2,2)
		REQUIRE(hhS.getHelixE(2, 2) == -2);
		REQUIRE(hhS.getHelixLength1(2, 2) == 3);
		REQUIRE(hhS.getHelixLength2(2, 2) == hhS.getHelixLength1(2, 2));

		// (2,3)
		REQUIRE(hhS.getHelixE(2, 3) == -1);
		REQUIRE(hhS.getHelixLength1(2, 3) == 2);
		REQUIRE(hhS.getHelixLength2(2, 3) == hhS.getHelixLength1(2, 3));

		// (3,0)
		REQUIRE(hhS.getHelixE(3, 0) == -1);
		REQUIRE(hhS.getHelixLength1(3, 0) == 2);
		REQUIRE(hhS.getHelixLength2(3, 0) == hhS.getHelixLength1(3, 0));

		// (3,1)
		REQUIRE(hhS.getHelixE(3, 1) == -1);
		REQUIRE(hhS.getHelixLength1(3, 1) == 2);
		REQUIRE(hhS.getHelixLength2(3, 1) == hhS.getHelixLength1(3, 1));

		// (3,2)
		REQUIRE(hhS.getHelixE(3, 2) == -1);
		REQUIRE(hhS.getHelixLength1(3, 2) == 2);
		REQUIRE(hhS.getHelixLength2(3, 2) == hhS.getHelixLength1(3, 2));

		// (3,3)
		REQUIRE(hhS.getHelixE(3, 3) == -1);
		REQUIRE(hhS.getHelixLength1(3, 3) == 2);
		REQUIRE(hhS.getHelixLength2(3, 3) == hhS.getHelixLength1(3, 3));


		// Several cases that do not allow a helix

		// (4,4)
		REQUIRE(hhS.getHelixE(4, 4) == E_INF);
		REQUIRE(hhS.getHelixLength1(4, 4) == 0);
		REQUIRE(hhS.getHelixLength2(4, 4) == hhS.getHelixLength1(4, 4));

		// (4,3)
		REQUIRE(hhS.getHelixE(0, 4) == E_INF);
		REQUIRE(hhS.getHelixLength1(0, 4) == 0);
		REQUIRE(hhS.getHelixLength2(0, 4) == hhS.getHelixLength1(0, 4));

	}



	SECTION("traceBackHelix 1: Everything is complementary", "[HelixHandlerStackingOnly]") {

		RnaSequence r1("r1", "GGGGG");
		RnaSequence r2("r2", "CCCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 0);
		HelixHandlerStackingOnly hhS(energy, hC);
		// fill the helix
		hhS.fillHelix(0, energy.size1()-1, 0, energy.size2()-1);

		// Case (0,0)
		//////////////////////
		Interaction interaction1(r1,r2);

		hhS.traceBackHelix(interaction1, 0, 0);


		// First / last base pair of helix
		REQUIRE(interaction1.basePairs.begin()->first == 0);
		REQUIRE(interaction1.basePairs.begin()->second == 4);

		REQUIRE(interaction1.basePairs.rbegin()->first == 2);
		REQUIRE(interaction1.basePairs.rbegin()->second == 2);



		// Case (2,1)
		//////////////////////

		Interaction interaction2(r1,r2);

		hhS.traceBackHelix(interaction2, 2, 1);


		// First / last base pair of helix
		REQUIRE(interaction1.basePairs.begin()->first == 0);
		REQUIRE(interaction1.basePairs.begin()->second == 4);

		REQUIRE(interaction1.basePairs.rbegin()->first == 2);
		REQUIRE(interaction1.basePairs.rbegin()->second == 2);


	}


	SECTION("fillHelix 2: Sequence 1 contains an A", "[HelixHandlerStackingOnly]") {
		// Case 2 - sequence containing an "A" to disrupt perfect stacking
		RnaSequence r1("r1", "GGGAGG");
		RnaSequence r2("r2", "CCCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 0);
		HelixHandlerStackingOnly hhS(energy, hC);
		REQUIRE(hhS.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1) == 12);


		// All optimal combinations
		// (0,0)
		REQUIRE(hhS.getHelixE(0, 0) == -2);
		REQUIRE(hhS.getHelixLength1(0, 0) == 3);
		REQUIRE(hhS.getHelixLength2(0, 0) == hhS.getHelixLength1(0, 0));

		// (0,1)
		REQUIRE(hhS.getHelixE(0, 1) == -2);
		REQUIRE(hhS.getHelixLength1(0, 1) == 3);
		REQUIRE(hhS.getHelixLength2(0, 1) == hhS.getHelixLength1(0, 1));

		// (0,2)
		REQUIRE(hhS.getHelixE(0, 2) == -2);
		REQUIRE(hhS.getHelixLength1(0, 2) == 3);
		REQUIRE(hhS.getHelixLength2(0, 2) == hhS.getHelixLength1(0, 2));

		// (0,3)
		REQUIRE(hhS.getHelixE(0, 3) == -1);
		REQUIRE(hhS.getHelixLength1(0, 3) == 2);
		REQUIRE(hhS.getHelixLength2(0, 3) == hhS.getHelixLength1(0, 3));

		// (1,0)
		REQUIRE(hhS.getHelixE(1, 0) == -1);
		REQUIRE(hhS.getHelixLength1(1, 0) == 2);
		REQUIRE(hhS.getHelixLength2(1, 0) == hhS.getHelixLength1(1, 0));

		// (1,1)
		REQUIRE(hhS.getHelixE(1, 1) == -1);
		REQUIRE(hhS.getHelixLength1(1, 1) == 2);
		REQUIRE(hhS.getHelixLength2(1, 1) == hhS.getHelixLength1(1, 1));

		// (1,2)
		REQUIRE(hhS.getHelixE(1, 2) == -1);
		REQUIRE(hhS.getHelixLength1(1, 2) == 2);
		REQUIRE(hhS.getHelixLength2(1, 2) == hhS.getHelixLength1(1, 2));

		// (1,3)
		REQUIRE(hhS.getHelixE(1, 3) == -1);
		REQUIRE(hhS.getHelixLength1(1, 3) == 2);
		REQUIRE(hhS.getHelixLength2(1, 3) == hhS.getHelixLength1(1, 3));

		// (4,0)
		REQUIRE(hhS.getHelixE(4, 0) == -1);
		REQUIRE(hhS.getHelixLength1(4, 0) == 2);
		REQUIRE(hhS.getHelixLength2(4, 0) == hhS.getHelixLength1(4, 0));

		// (4,1)
		REQUIRE(hhS.getHelixE(4, 1) == -1);
		REQUIRE(hhS.getHelixLength1(4, 1) == 2);
		REQUIRE(hhS.getHelixLength2(4, 1) == hhS.getHelixLength1(4, 1));

		// (4,2)
		REQUIRE(hhS.getHelixE(4, 2) == -1);
		REQUIRE(hhS.getHelixLength1(4, 2) == 2);
		REQUIRE(hhS.getHelixLength2(4, 2) == hhS.getHelixLength1(4, 2));

		// (4,3)
		REQUIRE(hhS.getHelixE(4, 3) == -1);
		REQUIRE(hhS.getHelixLength1(4, 3) == 2);
		REQUIRE(hhS.getHelixLength2(4, 3) == hhS.getHelixLength1(4, 3));


		// Not viable cases
		// (3,0) Not complementary
		REQUIRE(hhS.getHelixE(3,0) == E_INF);
		REQUIRE(hhS.getHelixLength1(3,0) == 0);
		REQUIRE(hhS.getHelixLength2(3,0) == hhS.getHelixLength1(3,0));

		// (5,4) no Helix possible (both sides too short)
		REQUIRE(hhS.getHelixE(5,4) == E_INF);
		REQUIRE(hhS.getHelixLength1(5,4) == 0);
		REQUIRE(hhS.getHelixLength2(5,4) == hhS.getHelixLength1(5,4));

		// (5,3) no helix possible
		REQUIRE(hhS.getHelixE(5,3) == E_INF);
		REQUIRE(hhS.getHelixLength1(5,3) == 0);
		REQUIRE(hhS.getHelixLength2(5,3) == hhS.getHelixLength1(5,3));

	}

	SECTION("traceBackHelix 2: Sequence 1 contains an A", "[HelixHandlerStackingOnly]") {

		RnaSequence r1("r1", "GGGAGG");
		RnaSequence r2("r2", "CCCCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 0);
		HelixHandlerStackingOnly hhS(energy, hC);
		// fill the helix
		hhS.fillHelix(0, energy.size1()-1, 0, energy.size2()-1);

		// Case (0,0)
		//////////////////////
		Interaction interaction1(r1,r2);
		hhS.traceBackHelix(interaction1, 0, 0);


		// First / last base pair of helix
		REQUIRE(interaction1.basePairs.begin()->first == 0);
		REQUIRE(interaction1.basePairs.begin()->second == 4);

		REQUIRE(interaction1.basePairs.rbegin()->first == 1);
		REQUIRE(interaction1.basePairs.rbegin()->second == 3);



		// Case (2,1)
		//////////////////////

		Interaction interaction2(r1,r2);

		hhS.traceBackHelix(interaction2, 2, 1);

		REQUIRE(interaction2.basePairs.size() == 0);

	}

	SECTION("fillHelix 3: A wall of A's disrupts the possible helices", "[HelixHandlerStackingOnly]") {
		// Case 2 - sequence containing an "A"-wall to disrupt perfect stacking
		RnaSequence r1("r1", "GGGAAGG");
		RnaSequence r2("r2", "CCAACCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 0);
		HelixHandlerStackingOnly hhS(energy, hC);
		REQUIRE(hhS.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1) == 9);

		REQUIRE_FALSE(energy.areComplementary(5,4));
		// All optimal combinations
		// (0,0)
		REQUIRE(hhS.getHelixE(0, 0) == -2);
		REQUIRE(hhS.getHelixLength1(0, 0) == 3);
		REQUIRE(hhS.getHelixLength2(0, 0) == hhS.getHelixLength1(0, 0));

		// (0,1)
		REQUIRE(hhS.getHelixE(0, 1) == -1);
		REQUIRE(hhS.getHelixLength1(0, 1) == 2);
		REQUIRE(hhS.getHelixLength2(0, 1) == hhS.getHelixLength1(0, 1));

		// (0,5)
		REQUIRE(hhS.getHelixE(0, 5) == -1);
		REQUIRE(hhS.getHelixLength1(0, 5) == 2);
		REQUIRE(hhS.getHelixLength2(0, 5) == hhS.getHelixLength1(0, 5));

		// (1,0)
		REQUIRE(hhS.getHelixE(1, 0) == -1);
		REQUIRE(hhS.getHelixLength1(1, 0) == 2);
		REQUIRE(hhS.getHelixLength2(1, 0) == hhS.getHelixLength1(1, 0));

		// (1,1)
		REQUIRE(hhS.getHelixE(1, 1) == -1);
		REQUIRE(hhS.getHelixLength1(1, 1) == 2);
		REQUIRE(hhS.getHelixLength2(1, 1) == hhS.getHelixLength1(1, 1));

		// (1,5)
		REQUIRE(hhS.getHelixE(1, 5) == -1);
		REQUIRE(hhS.getHelixLength1(1, 5) == 2);
		REQUIRE(hhS.getHelixLength2(1, 5) == hhS.getHelixLength1(1, 5));

		// (5,0)
		REQUIRE(hhS.getHelixE(5, 0) == -1);
		REQUIRE(hhS.getHelixLength1(5, 0) == 2);
		REQUIRE(hhS.getHelixLength2(5, 0) == hhS.getHelixLength1(5, 0));

		// (5,1)
		REQUIRE(hhS.getHelixE(5, 1) == -1);
		REQUIRE(hhS.getHelixLength1(5, 1) == 2);
		REQUIRE(hhS.getHelixLength2(5, 1) == hhS.getHelixLength1(5, 1));

		// (5,5)
		REQUIRE(hhS.getHelixE(5, 5) == -1);
		REQUIRE(hhS.getHelixLength1(5, 5) == 2);
		REQUIRE(hhS.getHelixLength2(5, 5) == hhS.getHelixLength1(5, 5));



		// Not viable cases
		// (3,0) Not complementary
		REQUIRE(hhS.getHelixE(3,0) == E_INF);
		REQUIRE(hhS.getHelixLength1(3,0) == 0);
		REQUIRE(hhS.getHelixLength2(3,0) == hhS.getHelixLength1(3,0));

		// (5,4) no Helix possible (both sides too short)
		REQUIRE(hhS.getHelixE(2,2) == E_INF);
		REQUIRE(hhS.getHelixLength1(2,2) == 0);
		REQUIRE(hhS.getHelixLength2(2,2) == hhS.getHelixLength1(2,2));

		// (5,3) no helix possible
		REQUIRE(hhS.getHelixE(5,3) == E_INF);
		REQUIRE(hhS.getHelixLength1(5,3) == 0);
		REQUIRE(hhS.getHelixLength2(5,3) == hhS.getHelixLength1(5,3));

	}

	SECTION("traceBackHelix 3: A wall of A's disrupts the possible helices", "[HelixHandlerStackingOnly]") {

		RnaSequence r1("r1", "GGGAAGG");
		RnaSequence r2("r2", "CCAACCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 0);
		HelixHandlerStackingOnly hhS(energy, hC);
		// fill the helix
		hhS.fillHelix(0, energy.size1()-1, 0, energy.size2()-1);

		// Case (0,0)
		//////////////////////
		Interaction interaction1(r1,r2);
		hhS.traceBackHelix(interaction1, 0, 0);

		REQUIRE(interaction1.basePairs.size() == 2);
		// First / last base pair of helix
		REQUIRE(interaction1.basePairs.begin()->first == 0);
		REQUIRE(interaction1.basePairs.begin()->second == 6);

		REQUIRE(interaction1.basePairs.rbegin()->first == 1);
		REQUIRE(interaction1.basePairs.rbegin()->second == 5);

		// Case (2,1)
		//////////////////////

		Interaction interaction2(r1,r2);
		hhS.traceBackHelix(interaction2, 2, 1);
		REQUIRE(interaction2.basePairs.size() == 0);

		// Case (5,5)
		//////////////////////

		Interaction interaction3(r1,r2);
		hhS.traceBackHelix(interaction2, 5, 5);

		REQUIRE(interaction2.basePairs.size() == 1);
		// First / last base pair of helix
		REQUIRE(interaction2.basePairs.begin()->first == 5);
		REQUIRE(interaction2.basePairs.begin()->second == 1);

		REQUIRE(interaction2.basePairs.rbegin()->first == 5);
		REQUIRE(interaction2.basePairs.rbegin()->second == 1);


	}

	SECTION("fillHelix 4: No interaction possible", "[HelixHandlerStackingOnly]") {
		// Case 4 -NO HELIX POSSIBLE
		RnaSequence r1("r1", "AAAAAAA");
		RnaSequence r2("r2", "AAAAAAA");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 0);
		HelixHandlerStackingOnly hhS(energy, hC);

		// No helix is possibles
		REQUIRE(hhS.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1) == 0);

		// NO POSSIBLE HELICES
		// (2,2)
		REQUIRE(hhS.getHelixE(2, 2) == E_INF);
		REQUIRE(hhS.getHelixLength1(2, 2) == 0);
		REQUIRE(hhS.getHelixLength2(2, 2) == hhS.getHelixLength1(2, 2));

		// (0,3)
		REQUIRE(hhS.getHelixE(0, 3) == E_INF);
		REQUIRE(hhS.getHelixLength1(0, 3) == 0);
		REQUIRE(hhS.getHelixLength2(0, 3) == hhS.getHelixLength1(0, 3));

	}

	SECTION("traceBackHelix 4: No interaction possible") {
		// Case 4 -NO HELIX POSSIBLE
		RnaSequence r1("r1", "AAAAAAA");
		RnaSequence r2("r2", "AAAAAAA");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hC(2, 4, 0);
		HelixHandlerStackingOnly hhS(energy, hC);

		// No helix is possibles
		REQUIRE(hhS.fillHelix(0, energy.size1() - 1, 0, energy.size2() - 1) == 0);

		// Case (0,0)
		//////////////////////
		Interaction interaction1(r1,r2);

		hhS.traceBackHelix(interaction1, 0, 0);

		REQUIRE(interaction1.basePairs.size() == 0);

		// Case (0,0)
		//////////////////////
		Interaction interaction2(r1,r2);

		hhS.traceBackHelix(interaction2, 1, 3);

		REQUIRE(interaction2.basePairs.size() == 0);

	}

}