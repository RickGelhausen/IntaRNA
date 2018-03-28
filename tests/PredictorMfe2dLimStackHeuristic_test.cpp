
#include "catch.hpp"

#undef NDEBUG
#define protected public

#include "IntaRNA/RnaSequence.h"
#include "IntaRNA/AccessibilityDisabled.h"
#include "IntaRNA/HelixConstraint.h"
#include "IntaRNA/HelixHandlerStackingOnly.h"
#include "IntaRNA/InteractionEnergyBasePair.h"
#include "IntaRNA/Interaction.h"
#include "IntaRNA/ReverseAccessibility.h"
#include "IntaRNA/PredictorMfe2dLimStackHeuristic.h"
#include "IntaRNA/OutputHandlerInteractionList.h"

using namespace IntaRNA;

TEST_CASE( "PredictorMfe2dLimStackHeuristc", "[PredictorMfe2dLimStackHeuristic]") {

	SECTION("getter", "[PredictorMfe2dLimStackHeuristic]") {

		RnaSequence r1("r1", "GGGAAGG");
		RnaSequence r2("r2", "CCAACCC");
		AccessibilityDisabled acc1(r1, 0, NULL);
		AccessibilityDisabled acc2(r2, 0, NULL);
		ReverseAccessibility racc(acc2);
		InteractionEnergyBasePair energy(acc1, racc);

		HelixConstraint hc(2, 4, 0);

		OutputHandlerInteractionList out(1);

		PredictorMfe2dLimStackHeuristic pLSH(energy, out, NULL, hc);

		IndexRange idx1(0,r1.lastPos);
		IndexRange idx2(0,r2.lastPos);
		OutputConstraint outC(1,OutputConstraint::OVERLAP_SEQ2,0,100);

		pLSH.predict(idx1,idx2,outC);

		REQUIRE_FALSE(out.empty());
		REQUIRE(out.reported() == 1);

		const Interaction * interaction((*out.begin()));
		REQUIRE(interaction->basePairs.begin()->first == 0);
		REQUIRE(interaction->basePairs.begin()->second == 6);

		REQUIRE(interaction->basePairs.rbegin()->first == 6);
		REQUIRE(interaction->basePairs.rbegin()->second == 0);
	}

}