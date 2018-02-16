
#ifndef INTARNA_HELIXCONSTRAINT_H
#define INTARNA_HELIXCONSTRAINT_H

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
	 * @param maxUnpaired1 the maximal number of unpaired bases within seq1
	 * 		  allowed within an helix
	 * @param maxUnpaired2 the maximal number of unpaired bases within seq2
	 * 		  allowed within an helix
	 */
	HelixConstraint( const size_t bpMin
				, const size_t bpMax
				, const size_t maxUnpaired1
			 	, const size_t maxUnpaired2
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
	 * Provides the maximal number of unpaired bases within the first sequence
	 * of an helix
	 *
	 * @return the maximal number of unpaired bases within the first sequence
	 *         which an helix is allowed to have
	 */
	size_t
	getMaxUnpaired1() const;

	/**
	 * Provides the maximal number of unpaired bases within the second sequence
	 * of an helix
	 *
	 * @return the maximal number of unpaired bases within the second sequence
	 *         an helix is allowed to have
	 */
	size_t
	getMaxUnpaired2() const;

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

	//! the maximally allowed number of unpaired bases in helix seq1
	size_t maxUnpaired1;

	//! the maximally allowed number of unpaired bases in helix seq2
	size_t maxUnpaired2;
};


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

inline
HelixConstraint::HelixConstraint(
		const size_t bpMin_
		, const size_t bpMax_
		, const size_t maxUnpaired1_
		, const size_t maxUnpaired2_)
	:
		bpMin(bpMin_)
	  , bpMax(bpMax_)
      , maxUnpaired1(maxUnpaired1_) // exclude too large boundaries
      , maxUnpaired2(maxUnpaired2_)
{
	// TODO: Check if 0 helix allowed
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
getMaxUnpaired1() const {
	return maxUnpaired1;
}

/////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixConstraint::
getMaxUnpaired2() const {
	return maxUnpaired2;
}

/////////////////////////////////////////////////////////////////////////////
// TODO: Maybe change this. Might not be the right max length
inline
size_t
HelixConstraint::
getMaxLength1() const {
	return getBasePairs() + getMaxUnpaired1();
}

/////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixConstraint::
getMaxLength2() const {
	return getBasePairs() + getMaxUnpaired2();
}

/////////////////////////////////////////////////////////////////////////////

inline
std::ostream&
operator<<(std::ostream& out, const HelixConstraint& c)
{
	out <<"HelixConstraint( bpMin="<<c.getMinBasePairs()
			<<", bpMax="<<c.getMaxBasePairs()
		    <<", up1="<<c.getMaxUnpaired1()
		    <<", up2="<<c.getMaxUnpaired2()
		    <<")";
	return out;

}

/////////////////////////////////////////////////////////////////////////////

} // namespace

#endif //INTARNA_HELIXCONSTRAINT_H
