
#ifndef INTARNA_HELIXCONSTRAINT_H
#define INTARNA_HELIXCONSTRAINT_H

#include "IntaRNA/general.h"
#include "IntaRNA/IndexRangeList.h"

#include <cstddef>
#include <iostream>

namespace IntaRNA {

/**
 * Encodes helix constraints to be used for interaction prediction.
 *
 * @author Rick Gelhausen
 *
 */
class HelixConstraint {

public:

	/**
	 *  Constructor
	 *
	 * @param bpMin minimal number of base pairs a helix is allowed to have (>= 2)
	 * @param bpMax maximal number of base pairs a helix is allowed to have (>= bpMin)
	 * @param maxUnpaired maximal number of unpaired bases
	 */
	HelixConstraint( const size_t bpMin
				, const size_t bpMax
				, const size_t maxUnpaired
				);

	virtual ~HelixConstraint();

	/**
	 *  Provides the minimum number of base pairs allowed within an helix (>=2)
	 *
	 * @return the minimum number of base pairs an helix is allowed to have (>=2)
	 */
	size_t
	getMinBasePairs() const;

	/**
	 *  Provides the maximum number of base pairs allowed within an helix (>=bpMin)
	 *
	 * @return the maximum number of base pairs an helix is allowed to have (>=bpMin)
	 */
	size_t
	getMaxBasePairs() const;

	/**
	 * Provides the maximal number of unpaired bases allowed for the helix
	 *
	 * @return the maximal number of unpaired bases allowed for the helix
	 */
	size_t
	getMaxUnpaired() const;

	// TODO: Not sure whether the true max length, but surely larger
	/**
     * Provides the maximal length of the helix in seq1
     * @return the maximal length of the helix in seq1
     */
	size_t
	getMaxLength1() const;

	/**
	 * Provides the maximal length of the helix in seq2
	 * @return the maximal length of the helix in seq2
	 */
	size_t
	getMaxLength2() const;

	/**
	 * Prints the helix constraint details to stream
	 * @param out the ostream to write to
	 * @param c the object to add
	 * @return the altered stream out
	 */
	friend std::ostream& operator<<(std::ostream& out, const HelixConstraint& c);

protected:

	//! the minimal number of base pairs allowed in the helix (>=2)
	size_t bpMin;

	//! the maximal number of base pairs allowed in the helix (>=bpMin)
	size_t bpMax;

	//! the maximally allowed number of unpaired bases in the helix
	size_t maxUnpaired;
};


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

inline
HelixConstraint::HelixConstraint(
		const size_t bpMin_
		, const size_t bpMax_
		, const size_t maxUnpaired_)
	:
		bpMin(bpMin_)
	  , bpMax(bpMax_)
      , maxUnpaired(maxUnpaired_) // exclude too large boundaries
{
	if (bpMin < 2) throw std::runtime_error("HelixConstraint() : base pair number ("+toString(bpMin)+") < 2");
}

/////////////////////////////////////////////////////////////////////////////

inline
HelixConstraint::~HelixConstraint() {
}

/////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixConstraint::
getMinBasePairs() const {
	return bpMin;
}

/////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixConstraint::
getMaxBasePairs() const {
	return bpMax;
}

/////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixConstraint::
getMaxUnpaired() const {
	return maxUnpaired;
}


/////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixConstraint::
getMaxLength1() const {
	return getMaxBasePairs() + getMaxUnpaired();
}

/////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixConstraint::
getMaxLength2() const {
	return getMaxBasePairs() + getMaxUnpaired();
}

/////////////////////////////////////////////////////////////////////////////

inline
std::ostream&
operator<<(std::ostream& out, const HelixConstraint& c)
{
	out <<"HelixConstraint( bpMin="<<c.getMinBasePairs()
			<<", bpMax="<<c.getMaxBasePairs()
		    <<", upMax="<<c.getMaxUnpaired()
		    <<")";
	return out;

}

/////////////////////////////////////////////////////////////////////////////

} // namespace

#endif //INTARNA_HELIXCONSTRAINT_H
