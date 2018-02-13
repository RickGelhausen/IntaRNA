
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
	 * @param maxBp maximal number of base pairs a helix is allowed to have (>= 2)
	 * @param maxUnpairedOverall the maximal (summed) number of unpaired bases
	 * 		  within both seq1 and seq2 allowed within an helix
	 * @param maxUnpaired1 the maximal number of unpaired bases within seq1
	 * 		  allowed within an helix
	 * @param maxUnpaired2 the maximal number of unpaired bases within seq2
	 * 		  allowed within an helix
	 */
	HelixConstraint( const size_t maxBp
				, const size_t maxUnpairedOverall
				, const size_t maxUnpaired1
				, const size_t maxUnpaired2
				);

	virtual ~HelixConstraint();

	/**
	 *  Provides the maximum number of base pairs allowed within an helix
	 *
	 * @return the maximum number of base pairs an helix is allowed to have
	 */
	size_t
	getMaxBasePairs() const;

	/**
	 * Provides the overall maximal number of unpaired bases within an helix
	 *
	 * @return the overall maximal number of unpaired bases which
	 *         an helix is allowed to have
	 */
	size_t
	getMaxUnpairedOverall() const;

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
	 * Prints the helix constraint details to stream
	 * @param out the ostream to write to
	 * @param c the object to add
	 * @return the altered stream out
	 */
	friend std::ostream& operator<<(std::ostream& out, const HelixConstraint& c);

protected:

	//! the maximal number of base pairs allowed in the helix
	size_t maxBp;

	//! the overall summed maximally allowed number of unpaired bases in an helix
	size_t maxUnpairedOverall;

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
		const size_t maxBp_
		, const size_t maxUnpairedOverall_
		, const size_t maxUnpaired1_
		, const size_t maxUnpaired2_)
	:
		maxBp(maxBp_)
      , maxUnpairedOverall(maxUnpairedOverall_)
      , maxUnpaired1(std::min(maxUnpaired1_, maxUnpairedOverall_)) // exclude too large boundaries
      , maxUnpaired2(std::min(maxUnpaired2_, maxUnpairedOverall_))
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
getMaxBasePairs() const {
	return maxBp;
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
inline
size_t
HelixConstraint::
getMaxUnpairedOverall() const {
	return maxUnpairedOverall;
}

/////////////////////////////////////////////////////////////////////////////

inline
std::ostream&
operator<<(std::ostream& out, const HelixConstraint& c)
{
	out <<"HelixConstraint( maxBp="<<c.getMaxBasePairs()
			<<", up="<<c.getMaxUnpairedOverall()
		    <<", up1="<<c.getMaxUnpaired1()
		    <<", up2="<<c.getMaxUnpaired2()
		    <<")";
	return out;

}

/////////////////////////////////////////////////////////////////////////////

} // namespace

#endif //INTARNA_HELIXCONSTRAINT_H
