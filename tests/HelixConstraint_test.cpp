
#include "catch.hpp"

#undef NDEBUG

#include "IntaRNA/HelixConstraint.h"

using namespace IntaRNA;

TEST_CASE( "HelixConstraint", "[HelixConstraint]" ) {

	SECTION( "getter", "[HelixConstraint]" ) {

		size_t minBP= 2, maxBP = 10, maxUP=2, maxED=0;
		bool withED=false;

		HelixConstraint hC( minBP, maxBP, maxUP, maxED, withED);

		// check data access
		REQUIRE( hC.getMinBasePairs() == 2);
		REQUIRE( hC.getMaxBasePairs() == 10);
		REQUIRE( hC.getMaxUnpaired() == 2);
		REQUIRE( hC.getMaxLength1() == 12);
		REQUIRE( hC.getMaxLength2() == 12);
		REQUIRE( hC.getMaxED() == 0);
		REQUIRE( hC.withED() == false);
	}


}