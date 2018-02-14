#ifndef INTARNA_HELIXHANDLER_H
#define INTARNA_HELIXHANDLER_H

#include "IntaRNA/InteractionEnergy.h"
#include "IntaRNA/HelixConstraint.h"

#include <boost/multi_array.hpp>

#include <boost/numeric/ublas/matrix.hpp>

namespace IntaRNA {

class HelixHandler {

public:

	//! 6D matrix type to hold the mfe energies for helix interactions
	//! of the ranges i1..(i1+bp+u1-1) and i2..(i2+bp+u2-1), with
	//! i1,i2 = the start index of the helix in seq1/2
	//! bp = the number of base pairs within the helix
	//! bpInbetween = the number of base pairs enclosed by left and right base pair, ie == (bp-2)
	//! u1/u2 = the number of unpaired positions within the helix
	//! using the index [i1][i2][bpInbetween][u1][u2][bp] or a SeedIndex object
	typedef boost::multi_array<E_type,6> HelixRecMatrix;

	//! defines the helix data {{ i1, i2, bpInbetween, u1, u2, bp }} to acces elemetns of
	//! the HelixRecMatrix
	typedef boost::array<HelixRecMatrix::index, 6> HelixIndex;

	//! matrix to store the helix information for each helix left side (i1, i2)
	//! it holds both the enrgy (first) as well as the length of the helix using
	//! the length combination of encodeHelixLength()
	typedef boost::numeric::ublas::matrix< std::pair<E_type, size_t> > HelixMatrix;

public:

	/**
	 * Constructor
	 * @param energy the energy function to be used
	 */
	HelixHandler(
			const InteractionEnergy & energy
			, const HelixConstraint & helixConstraint
			);

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
	 * Access to the underlying helix constraint
	 * @return the used helix constraint
	 */
	virtual
	const HelixConstraint&
	getConstraint() const;

	/**
	 * Compute the helix matrix for the given interval boundaries
	 * @param i1 the first index of seq1 that might interact
	 * @param j1 the last index of seq1 that might interact
	 * @param i2 the first index of seq2 that might interact
	 * @param j2 the last index of seq2 that might interact
	 * @return
	 */
	virtual
	size_t
	fillHelix( const size_t i1, const size_t j1, const size_t i2, const size_t j2 );

	/**
	 * Identifies the base pairs of the mfe helix interaction starting at i1,i2
	 * and writes them to the provided container
	 *
	 * Note: the right most base pair is excluded!
	 *
	 * @param interaction the container to add the base pairs too
	 * @param i1 the start of the helix in seq1
	 * @param i2 the start of the helix in seq2
	 */
	virtual
	void
	traceBackHelix( Interaction & interaction, const size_t i1, const size_t i2);

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

	//! the helix constraint to be applied
	const HelixConstraint & helixConstraint;

	//! the recursion data for the computation of a helix interaction
	//! i1..(i1+bpInbetween+u1-1) and i2..(i2+bpInbetween+u2-1)
	//! using the indexing [i1][i2][bpInbetween][u1][u2][bp]
	HelixRecMatrix helixE_rec;

	//! the helix mfe information for helix starting at (i1, i2)
	HelixMatrix helix;

	//! offset for seq1 indices for the current matrices
	size_t offset1;

	//! offset for seq2 indices for the current matrices
	size_t offset2;


	/**
	 *
	 * @param i1
	 * @param i2
	 * @param bpInbetween
	 * @param u1
	 * @param u2
	 * @param bp
	 * @return
	 */
	E_type
	getHelixE( const size_t i1, const size_t i2, const size_t bpInbetween, const size_t u1, const size_t u2, const size_t bp );

};

}

#endif //INTARNA_HELIXHANDLER_H
