#include "IntaRNA/HelixHandler.h"
#include "IntaRNA/HelixHandlerStackingOnly.h"

namespace IntaRNA {

HelixHandler* HelixHandler::getHelixHandler(const InteractionEnergy &energy,
											const HelixConstraint &helixConstraint,
											SeedHandler * const seedHandler) {

	if (helixConstraint.getMaxUnpaired() == 0) {
		return new HelixHandlerStackingOnly(energy, helixConstraint, seedHandler);
	} else {
		INTARNA_NOT_IMPLEMENTED("HelixHandlerUnpaired!!");
	}
}

}