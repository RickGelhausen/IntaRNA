
#ifndef INTARNA_HELIXHANDLERSTACKINGONLY_H_
#define INTARNA_HELIXHANDLERSTACKINGONLY_H_

#include "IntaRNA/InteractionEnergy.h"
#include "IntaRNA/HelixConstraint.h"

#include <boost/multi_array.hpp>

#include <boost/numeric/ublas/matrix.hpp>

namespace IntaRNA {

/**
 * Handler to provide helix interaction information
 * Simple Version without allowing unpaired regions
 */

class HelixHandlerStackingOnly {

public:
	//! 4D matrix type to hold the mfe energies for helix interactions
	//! of the ranges i1..(i1+bp-1) and i2..(i2+bp-1), with
	//! i1,i2 = the start index of the helix in seq1/2
	//! bp = the current number of base pairs in the helix (< maxBP)
	//! using the index [i1][i2][bpMax][bp] or a HelixIndex object
	typedef boost::multi_array<E_type,3> HelixRecMatrix;

	//! defines the helix data {{ i1, i2, bp }} to acces elements of
	//! the HelixRecMatrix
	typedef boost::array<HelixRecMatrix::index, 3> HelixIndex;

	//! matrix to store the helix information for each helix left side (i1, i2)
	//! it holds both the energy (first) as well as the length of the helix using
	//! the length combination of encodeHelixLength()
	typedef boost::numeric::ublas::matrix< std::pair<E_type, size_t> > HelixMatrix;

public:

	/**
	 * Constructor
	 * @param energy the energy function to be used
	 */
	HelixHandlerStackingOnly(
			const InteractionEnergy & energy
			, const HelixConstraint & helixConstraint
	);

	/**
	 * destructor
	 */
	virtual ~HelixHandlerStackingOnly();

	/**
	 * Access to the underlying interaction energy function
	 * @return the underlying energy function
	 */
	virtual
	const InteractionEnergy&
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
	 * @return number of valid helices
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
	//! bp: the number of allowed bases
	//! i1..(i1+bp-1) and i2..(i2+bp-1)
	//! using the indexing [i1][i2][bp]
	HelixRecMatrix helixE_rec;

	//! the helix mfe information for helix starting at (i1, i2)
	HelixMatrix helix;

	//! offset for seq1 indices for the current matrices
	size_t offset1;

	//! offset for seq2 indices for the current matrices
	size_t offset2;


	/**
	 * Provides the helix energy during recursion
	 *
	 * @param i1 the helix left end in seq 1 (index including offset)
	 * @param i2 the helix left end in seq 2 (index including offset)
	 * @param bp the current number of base pairs in the helix
	 *
	 * @return the energy of the according helix
	 */
	E_type
	getHelixE( const size_t i1, const size_t i2, const size_t bp );

	/**
	 * Get the number of bp in the mfe helix
	 * @param i1 the helix left end in seq 1
	 * @param i2 the helix left end in seq 2
	 *
	 * @return the number of base pairs in the mfe helix for (i1,i2)
  	 */
	size_t
	getHelixBP( const size_t i1, const size_t i2 ) const;

	/**
	 * Fills the helix energy during recursion
	 *
	 * @param i1 the helix left end in seq 1 (index including offset)
	 * @param i2 the helix left end in seq 2 (index including offset)
	 * @param bp the current number of base pairs in the helix
	 * @param E the energy value to be set
	 */
	void
	setHelixE( const size_t i1, const size_t i2, const size_t bp, const E_type E );

};

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

inline
HelixHandlerStackingOnly::HelixHandlerStackingOnly(
		const InteractionEnergy & energy
		, const HelixConstraint & helixConstraint
)
		:
		energy(energy)
		, helixConstraint(helixConstraint)
		, helixE_rec( HelixIndex({{ 0, 0, 0 }}))
		, helix()
		, offset1(0)
		, offset2(0)
{
}

////////////////////////////////////////////////////////////////////////////

inline
HelixHandlerStackingOnly::~HelixHandlerStackingOnly()
{
}

////////////////////////////////////////////////////////////////////////////

inline
const HelixConstraint&
HelixHandlerStackingOnly::getConstraint() const
{
	return helixConstraint;
}

////////////////////////////////////////////////////////////////////////////

inline
const InteractionEnergy&
HelixHandlerStackingOnly::getInteractionEnergy() const
{
	return energy;
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
HelixHandlerStackingOnly::
getHelixE(const size_t i1, const size_t i2) const
{
	return helix(i1-offset1, i2-offset2).first;
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerStackingOnly::
getHelixBP(const size_t i1, const size_t i2) const
{
	return helix(i1-offset1, i2-offset2).second;
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerStackingOnly::
getHelixLength1(const size_t i1, const size_t i2) const
{
	return helix(i1-offset1, i2-offset2).second;
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerStackingOnly::
getHelixLength2(const size_t i1, const size_t i2) const
{
	return helix(i1-offset1, i2-offset2).second;
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
HelixHandlerStackingOnly::
getHelixE(const size_t i1, const size_t i2, const size_t bp )
{
// return helixE_rec[i1][i2][bpInbetween][u1][u2][bp]
	return helixE_rec( HelixIndex({{
										   (HelixRecMatrix::index) i1
										   , (HelixRecMatrix::index) i2
										   , (HelixRecMatrix::index) bp}}) );
}

////////////////////////////////////////////////////////////////////////////

inline
void
HelixHandlerStackingOnly::
setHelixE(const size_t i1, const size_t i2, const size_t bp, const E_type E)
{
//	helixE_rec[i1][i2][bpMax][bp] = E
	helixE_rec( HelixIndex({{
									(HelixRecMatrix::index) i1
									, (HelixRecMatrix::index) i2
									, (HelixRecMatrix::index) bp}}) ) = E;
}


} // namespace
#endif /* HELIXHANDLERSTACKINGONLY_H_ */
