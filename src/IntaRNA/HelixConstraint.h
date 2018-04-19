
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
	 * @param minBP minimal number of base pairs a helix is allowed to have (>= 2)
	 * @param maxBP maximal number of base pairs a helix is allowed to have (>= bpMin)
	 * @param maxUnpaired maximal number of unpaired bases
	 * @param maxED maximal allowed energy for ED-values
	 * @parem withED whether or not ED values are added in the computation of helices
	 */
	HelixConstraint( const size_t minBP
				, const size_t maxBP
				, const size_t maxUnpaired
				, const E_type maxED
			    , const bool withED
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
	 *  Provides the maximum number of base pairs allowed within an helix (>=minBP)
	 *
	 * @return the maximum number of base pairs an helix is allowed to have (>=maxBP)
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
	 * Provides the maximally allowed ED value for each seed subsequence to be considered
	 * @return the maximally allowed ED value (per sequence) for a helix
	 */
	E_type
	getMaxED() const;

	/**
	 * Whether or not ED- values are used in the helix computation
	 * @return true / false
	 */
	bool
	withED() const;

	/**
	 * Provides the minimum interior loop length depending on the maximum allowed number of unpaired bases
	 * @return minimum interior loop length
	 */
	size_t
	getMinInternalLoopSize() const;

	/**
	 * Prints the helix constraint details to stream
	 * @param out the ostream to write to
	 * @param c the object to add
	 * @return the altered stream out
	 */
	friend std::ostream& operator<<(std::ostream& out, const HelixConstraint& c);

protected:

	//! the minimal number of base pairs allowed in the helix (>=2)
	size_t minBP;

	//! the maximal number of base pairs allowed in the helix (>=maxBP)
	size_t maxBP;

	//! the maximally allowed number of unpaired bases in the helix
	size_t maxUnpaired;

	//! the maximal ED value (per sequence) allowed for a helix
	E_type maxED;

	//! decision variable, whether ED-values are used in the computation of helix energies or not
	bool wED;
};


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

inline
HelixConstraint::HelixConstraint(
		const size_t minBP_
		, const size_t maxBP_
		, const size_t maxUnpaired_
		, const E_type maxED_
		, const bool withED_)
	:
		minBP(minBP_)
	  , maxBP(maxBP_)
      , maxUnpaired(maxUnpaired_)
	  , maxED(maxED_)
	  , wED(withED_)
{
	if (minBP < 2) throw std::runtime_error("HelixConstraint() : minimal base pair number ("+toString(minBP)+") < 2");
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
	return minBP;
}

/////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixConstraint::
getMaxBasePairs() const {
	return maxBP;
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
E_type
HelixConstraint::
getMaxED() const {
	return maxED;
}

/////////////////////////////////////////////////////////////////////////////

inline
bool
HelixConstraint::
withED() const {
	return wED;
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
size_t
HelixConstraint::
getMinInternalLoopSize() const {
	return getMaxUnpaired();
}

/////////////////////////////////////////////////////////////////////////////

inline
std::ostream&
operator<<(std::ostream& out, const HelixConstraint& c)
{
	out <<"HelixConstraint( minBP="<<c.getMinBasePairs()
			<<", maxBP="<<c.getMaxBasePairs()
		    <<", maxUP="<<c.getMaxUnpaired()
			<<", maxED="<<c.getMaxED()
			<<", withED="<<c.withED()
		    <<")";
	return out;

}

/////////////////////////////////////////////////////////////////////////////

} // namespace

#endif //INTARNA_HELIXCONSTRAINT_H
