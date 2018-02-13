#ifndef INTARNA_HELIXHANDLER_H
#define INTARNA_HELIXHANDLER_H

namespace IntaRNA {

class HelixHandler {

public:
	/**
	 * Constructor
	 * @param energy the energy function to be used
	 */
	HelixHandler( const InteractionEnergy & energy );

	/**
	 * destructor
	 */
	virtual ~HelixHandler();

	/**
	 * Access to the underlying interaction energy function
	 * @return the underlying energy function
	 */
	virtual
	const IntaractionEnergy&
	getInteractionEnergy() const;

	/**
	 * Compute the helix matrix for the given interval boundaries
	 * @param i1 the first index of seq1 that might interact
	 * @param j1 the last index of seq1 that might interact
	 * @param i2 the first index of seq2 that might interact
	 * @param j2 the last index of seq2 that might interact
	 * @return
	 */
	virtual
	//size_t
	fillHelix( const size_t i1, const size_t j1, const size_t i2, const size_t j2 );

	/**
	 * Access to the mfe of any helix with left-most base pair (i1, i2)
	 * @param i1 the left most interaction base of seq1
	 * @param i2 the left most interaction base of seq2
	 * @return the mfe of any helix starting at (i1,i2) or E_INF if none possible
	 */
	virtual
	E_type getHelixE( const size_t i1, const size_t i2 ) const;

	/**
	 * Access to the length in seq1 of the mfe helix with left-most base pair (i1,i2)
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interacting base of seq2
	 * @return the length in seq1 of the mfe helix starting at (i1,i2) or 0 if none possible
	 */
	virtual
	size_t
	getHelixLength1( const size_t i1, const size_t i2 ) const;

	/**
	 * Access to the length in seq2 of the mfe helix with left-most base pair (i1,i2)
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interacting base of seq2
	 * @return the length in seq2 of the mfe helix starting at (i1,i2) or 0 if none possible
	 */
	virtual
	size_t
	getHelixLength2( const size_t i1, const size_t i2 ) const;

protected:

	//! the used energy function
	const InteractionEnergy& energy;


};

}

#endif //INTARNA_HELIXHANDLER_H
