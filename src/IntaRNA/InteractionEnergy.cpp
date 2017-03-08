
#include "IntaRNA/InteractionEnergy.h"

#include <cmath>

#include <iostream>

namespace IntaRNA {

////////////////////////////////////////////////////////////////////////////

const size_t InteractionEnergy::minDistES = 3;

////////////////////////////////////////////////////////////////////////////

InteractionEnergy::
EnergyContributions
InteractionEnergy::
getE_contributions( const Interaction & interaction ) const
{

	// temporary access to range indices
	const size_t i1 = interaction.basePairs.begin()->first;
	const size_t i2 = getAccessibility2().getReversedIndex(interaction.basePairs.begin()->second);
	const size_t j1 = interaction.basePairs.rbegin()->first;
	const size_t j2 = getAccessibility2().getReversedIndex(interaction.basePairs.rbegin()->second);

	// fill contribution data structure
	EnergyContributions contr;
	contr.init = getE_init();
	contr.ED1 = getED1( i1, j1 );
	contr.ED2 = getED2( i2, j2 );
	contr.dangleLeft = (getE_danglingLeft( i1, i2 )*getPr_danglingLeft(i1,j1,i2,j2));
	contr.dangleRight = (getE_danglingRight( j1, j2 )*getPr_danglingRight(i1,j1,i2,j2));
	contr.endLeft = getE_endLeft( i1, i2 );
	contr.endRight = getE_endRight( j1, j2 );
	// compute loop energy
	contr.loops = interaction.energy
					- contr.init
					- contr.ED1
					- contr.ED2
					- contr.dangleLeft
					- contr.dangleRight
					- contr.endLeft
					- contr.endRight
					;

	// final contribution distribution
	return contr;
}

////////////////////////////////////////////////////////////////////////////


} // namespace
