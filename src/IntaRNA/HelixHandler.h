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
	//! of the ranges i1..(i1+bpMax+u1-1) and i2..(i2+bpMax+u2-1), with
	//! i1,i2 = the start index of the helix in seq1/2
	//! bpMax = the maximal number of base pairs within the helix (>=bpMin)
	//! bpInbetween = the number of base pairs enclosed by left and right base pair, ie bpInBetween <= bpMax-2
	//! u1/u2 = the number of unpaired positions within the helix
	//! using the index [i1][i2][bpInbetween][u1][u2][bpMax] or a HelixIndex object
	typedef boost::multi_array<E_type,6> HelixRecMatrix;

	//! defines the helix data {{ i1, i2, bpInbetween, u1, u2, bpMax }} to acces elements of
	//! the HelixRecMatrix
	typedef boost::array<HelixRecMatrix::index, 6> HelixIndex;

	//! matrix to store the helix information for each helix left side (i1, i2)
	//! it holds both the energy (first) as well as the length of the helix using
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
	//! bp: the number of allowed bases
	//! i1..(i1+bpInbetween+u1-1) and i2..(i2+bpInbetween+u2-1)
	//! using the indexing [i1][i2][bp][bpInbetween][u1][u2]
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
	 * @param bp the number of currently allowed base pairs
	 * @param bpInbetween the number of helix base pairs enclosed by the left-
	 * 		  and right-most base pair, ie. bpInBetween <= bpMax-2
	 * @param u1 the number of unpaired bases within seq 1
	 * @param u2 the number of unpaired bases within seq 2
	 *
	 * @return the energy of the according helix
	 */
	E_type
	getHelixE( const size_t i1, const size_t i2, const size_t bp, const size_t bpInbetween
				, const size_t u1, const size_t u2 );


	/**
	 * Fills the helix energy during recursion
	 *
	 * @param i1 the helix left end in seq 1 (index including offset)
	 * @param i2 the helix left end in seq 2 (index including offset)
	 * @param bp the number of currently allowed base pairs bpMin <= bp <= bpMax
	 * @param bpInbetween the number of helix base pairs enclosed by the left-
	 * 		  and right-most base pair, ie bpInBetween <= bp-2
	 * @param u1 the number of unpaired bases within seq 1
	 * @param u2 the number of unpaired bases within seq 2
	 * @param E the energy value to be set
	 */
	void
	setHelixE( const size_t i1, const size_t i2, const size_t bp, const size_t bpInbetween
				, const size_t u1, const size_t u2, const E_type E );

	/**
	 * Encodes the helix lengths into one number
	 * @param l1 the length of the helix in seq1
	 * @param l2 the length of the helix in seq2
	 * @return the combined encoding = (l1 + l2*(max_l1+1))
	 */
	size_t
	encodeHelixLength( const size_t l1, const size_t l2 ) const;

	/**
	 * Decodes the length of the helix within sequence 1 from an encoding
	 * generated with encodeHelixLength
	 * @param code the lengths encoding
	 * @return the length of the helix in seq1
	 */
	size_t
	decodeHelixLength1( const size_t code ) const;


	/**
	 * Decodes the length of the helix within sequence 2 from an encoding
	 * generated with encodeHelixLength
	 * @param code the lengths encoding
	 * @return the length of the helix in seq2
	 */
	size_t
	decodeHelixLength2( const size_t code ) const;

	// TODO: Might have to change this
	/**
	 * Fills a given interaction with the according
	 * hybridizing base pairs of the provided helix interaction TODO:??? (excluding the right mosot helix base pair)
	 *
	 * @param interaction IN/OUT the interaction to fill
	 * @param i1 the helix left end in seq 1 (index including offset)
	 * @param i2 the helix left end in seq 2 (index including offset)
	 * @param bpMin the minimal number of base pairs allowed within the helix
	 * @param bpMax the maximal number of base pairs allowed within the helix
	 * @param u1 the number of unpaired bases within seq 1
	 * @param u2 the number of unpaired bases within seq 2
	 */
	void
	traceBackHelix( Interaction & interaction
			, const size_t i1, const size_t i2, const size_t bpMin, const size_t bpMax
			, const size_t u1, const size_t u2 );
};


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

inline
HelixHandler::HelixHandler(
		const InteractionEnergy & energy
		, const HelixConstraint & helixConstraint
)
		:
		energy(energy)
		, helixConstraint(helixConstraint)
		, helixE_rec( HelixIndex({{ 0, 0, 0, 0, 0, 0 }}))
		, helix()
		, offset1(0)
		, offset2(0)
{
}

////////////////////////////////////////////////////////////////////////////

inline
HelixHandler::~HelixHandler()
{
}

////////////////////////////////////////////////////////////////////////////

inline
const HelixConstraint&
HelixHandler::getConstraint() const
{
	return helixConstraint;
}

////////////////////////////////////////////////////////////////////////////

inline
const InteractionEnergy&
HelixHandler::getInteractionEnergy() const
{
	return energy;
}

////////////////////////////////////////////////////////////////////////////

// TODO: REWORK THIS
inline
void
HelixHandler::traceBackHelix(Interaction &interaction
		, const size_t i1
		, const size_t i2
)
{
#if INTARNA_IN_DEBUG_MODE
	if ( i1 < offset1 ) throw std::runtime_error("HelixHandler::traceBackHelix(i1="+toString(i1)+") is out of range (>"+toString(offset1)+")");
	if ( i1-offset1 >= helix.size1() ) throw std::runtime_error("HelixHandler::traceBackHelix(i1="+toString(i1)+") is out of range (<"+toString(helix.size1()+offset1)+")");
	if ( i2 < offset2 ) throw std::runtime_error("HelixHandler::traceBackHelix(i2="+toString(i2)+") is out of range (>"+toString(offset2)+")");
	if ( i2-offset2 >= helix.size2() ) throw std::runtime_error("HelixHandler::traceBackHelix(i2="+toString(i2)+") is out of range (<"+toString(helix.size2()+offset2)+")");
	if ( E_isINF( getSeedE(i1,i2) ) ) throw std::runtime_error("HelixHandler::traceBackHelix(i1="+toString(i1)+",i2="+toString(i2)+") no helix known (E_INF)");
	if ( i1+getSeedLength1(i1,i2)-1-offset1 >= helix.size1() ) throw std::runtime_error("HelixHandler::traceBackHelix(i1="+toString(i1)+") helix length ("+toString(getHelixLength1(i1,i2))+") exceeds of range (<"+toString(helix.size1()+offset1)+")");
	if ( i2+getSeedLength2(i1,i2)-1-offset2 >= helix.size2() ) throw std::runtime_error("HelixHandler::traceBackHelix(i2="+toString(i2)+") helix length ("+toString(getHelixLength2(i1,i2))+") exceeds of range (<"+toString(helix.size2()+offset2)+")");
#endif
	// get number of base pairs allowed within the helix
	const size_t maxHelixBps = getConstraint().getMaxBasePairs();
	const size_t minHelixBps = getConstraint().getMinBasePairs();

	// trace back the according helix
	traceBackHelix( interaction, i1-offset1, i2-offset2
			, minHelixBps
			, maxHelixBps
			, getHelixLength1(i1,i2)-maxHelixBps-2
			, getHelixLength2(i1,i2)-maxHelixBps-2 );
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
HelixHandler::
getHelixE(const size_t i1, const size_t i2) const
{
	return helix(i1-offset1, i2-offset2).first;
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandler::
getHelixLength1(const size_t i1, const size_t i2) const
{
	return decodeHelixLength1(seed(i1-offset1, i2-offset2).second);
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandler::
getHelixLength2(const size_t i1, const size_t i2) const
{
	return decodeHelixLength2(seed(i1-offset1, i2-offset2).second);
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
HelixHandler::
getHelixE(const size_t i1, const size_t i2, const size_t bp, const size_t bpInbetween
		, const size_t u1, const size_t u2)
{
// return helixE_rec[i1][i2][bpInbetween][u1][u2][bp]
	return helixE_rec( HelixIndex({{
										   (HelixRecMatrix::Index) i1
										   , (HelixRecMatrix::Index) i2
										   , (HelixRecMatrix::Index) bp
										   , (HelixRecMatrix::Index) bpInbetween
										   , (HelixRecMatrix::Index) u1
										   , (HelixRecMatrix::Index) u2}}) );
}

////////////////////////////////////////////////////////////////////////////

inline
void
HelixHandler::
setHelixE(const size_t i1, const size_t i2, const size_t bp, const size_t bpInbetween
		, const size_t u1, const size_t u2, const E_type E)
{
//	helixE_rec[i1][i2][bpInbetween][u1][u2][bp] = E
	helixE_rec( HelixIndex({{
									(HelixRecMatrix::index) i1
									, (HelixRecMatrix::index) i2
									, (HelixRecMatrix::index) bp
									, (HelixRecMatrix::index) bpInbetween
									, (HelixRecMatrix::index) u1
									, (HelixRecMatrix::index) u2}}) );
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandler::
encodeHelixLength(const size_t l1, const size_t l2) const
{
	return l1 + l2 * (helixConstraint.getMaxLength1()+1);
}

inline
size_t
HelixHandler::
decodeHelixLength1(const size_t code) const
{
	return code % (helixConstraint.getMaxLength1()+1);
}

inline
size_t
HelixHandler::
decodeHelixLength2(const size_t code) const
{
	return code / (helixConstraint.getMaxLength1()+1);
}


} // namespace

#endif //INTARNA_HELIXHANDLER_H