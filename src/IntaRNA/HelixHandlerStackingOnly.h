
#ifndef INTARNA_HELIXHANDLERSTACKINGONLY_H_
#define INTARNA_HELIXHANDLERSTACKINGONLY_H_

#include "IntaRNA/InteractionEnergy.h"
#include "IntaRNA/HelixConstraint.h"
#include "IntaRNA/HelixHandler.h"

#include <boost/numeric/ublas/matrix.hpp>

namespace IntaRNA {

/**
 * Handler to provide helix interaction information
 * Simple Version without allowing unpaired regions
 */

class HelixHandlerStackingOnly : public HelixHandler {

public:

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
			, SeedHandler * const seedHandler = NULL
	);

	/**
	 * destructor
	 */
	virtual ~HelixHandlerStackingOnly();


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
	 * Compute the helix matrix, containing a seed, for the given interval boundaries
	 * @param i1 the first index of seq1 that might interact
	 * @param j1 the last index of seq1 that might interact
	 * @param i2 the first index of seq2 that might interact
	 * @param j2 the last index of seq2 that might interact
	 * @return
	 */
	virtual
	size_t
	fillHelixSeed( const size_t i1, const size_t j1, const size_t i2, const size_t j2 );

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
	 * Identifies the base pairs of the mfe helix interaction, containing a seed, starting at i1,i2
	 * and writes them to the provided container
	 *
	 * @param interaction the container to add the base pairs too
	 * @param i1 the start of the helix in seq1
	 * @param i2 the start of the helix in seq2
	 */
	virtual
	void
	traceBackHelixSeed( Interaction & interaction, const size_t i1, size_t i2);

	/**
	 * Access to the mfe of any helix with left-most base pair (i1, i2)
	 * @param i1 the left most interaction base of seq1
	 * @param i2 the left most interaction base of seq2
	 * @return the mfe of any helix starting at (i1,i2) or E_INF if none possible
	 */
	virtual
	E_type getHelixE( const size_t i1, const size_t i2 ) const;

	/**
	 * Access to the mfe of any helix, containing a seed, with left-most base pair (i1, i2)
	 * @param i1 the left most interaction base of seq1
	 * @param i2 the left most interaction base of seq2
	 * @return the mfe of any helix starting at (i1,i2) or E_INF if none possible
	 */
	virtual
	E_type getHelixSeedE( const size_t i1, const size_t i2 ) const;

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


	/**
	 * Access to the length in seq1 of the mfe helix, containing a seed, with left-most base pair (i1,i2)
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interacting base of seq2
	 * @return the length in seq1 of the mfe helix starting at (i1,i2) or 0 if none possible
	 */
	virtual
	size_t
	getHelixSeedLength1( const size_t i1, const size_t i2 ) const;

	/**
	 * Access to the length in seq2 of the mfe helix, containing a seed, with left-most base pair (i1,i2)
	 * @param i1 the left most interacting base of seq1
	 * @param i2 the left most interacting base of seq2
	 * @return the length in seq2 of the mfe helix starting at (i1,i2) or 0 if none possible
	 */
	virtual
	size_t
	getHelixSeedLength2( const size_t i1, const size_t i2 ) const;

	/**
	 * Set the seedHandler in order to compute helixSeed
	 * @param seedHandler seedHandler to be used in the helix computation
	 */
	void setSeedHandler(SeedHandler * const seedHandler);

protected:
	/**
	 * Encodes the seed lengths into one number
	 * @param l1 the length of the seed in seq1
	 * @param l2 the length of the seed in seq2
	 * @return the combined encoding = (l1 + l2*(max_l1+1))
	 */
	size_t
	encodeHelixLength( const size_t l1, const size_t l2 ) const;

	/**
	 * Decodes the length of the seed within sequence 1 from an encoding
	 * generated with encodeSeedLength()
	 * @param code the lengths encoding
	 * @return the length of the seed in seq1
	 */
	size_t
	decodeHelixLength1( const size_t code ) const;

	/**
	 * Decodes the length of the seed within sequence 2 from an encoding
	 * generated with encodeSeedLength()
	 * @param code the lengths encoding
	 * @return the length of the seed in seq2
	 */
	size_t
	decodeHelixLength2( const size_t code ) const;

	/**
	 * Encodes the seed lengths into one number
	 * @param l1 the length of the seed in seq1
	 * @param l2 the length of the seed in seq2
	 * @return the combined encoding = (l1 + l2*(max_l1+1))
	 */
	size_t
	encodeHelixSeedLength( const size_t l1, const size_t l2 ) const;

	/**
	 * Decodes the length of the seed within sequence 1 from an encoding
	 * generated with encodeSeedLength()
	 * @param code the lengths encoding
	 * @return the length of the seed in seq1
	 */
	size_t
	decodeHelixSeedLength1( const size_t code ) const;

	/**
	 * Decodes the length of the seed within sequence 2 from an encoding
	 * generated with encodeSeedLength()
	 * @param code the lengths encoding
	 * @return the length of the seed in seq2
	 */
	size_t
	decodeHelixSeedLength2( const size_t code ) const;

protected:

	//! the helix mfe information for helix starting at (i1, i2)
	HelixMatrix helix;

	//! the helix mfe information for helix with seed starting at (i1, i2)
	HelixMatrix helixSeed;

	//! offset for seq1 indices for the current matrices
	size_t offset1;

	//! offset for seq2 indices for the current matrices
	size_t offset2;


	SeedHandler * seedHandler;
};

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

inline
HelixHandlerStackingOnly::HelixHandlerStackingOnly(
		const InteractionEnergy & energy
		, const HelixConstraint & helixConstraint
		, SeedHandler * const seedHandler
)
		:
		HelixHandler(energy, helixConstraint, seedHandler)
		, seedHandler(seedHandler)
		, helix()
		, helixSeed()
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
E_type
HelixHandlerStackingOnly::
getHelixE(const size_t i1, const size_t i2) const
{
	return helix(i1-offset1, i2-offset2).first;
}

////////////////////////////////////////////////////////////////////////////

inline
E_type
HelixHandlerStackingOnly::
getHelixSeedE(const size_t i1, const size_t i2) const
{
	return helixSeed(i1-offset1, i2-offset2).first;
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerStackingOnly::
getHelixLength1(const size_t i1, const size_t i2) const
{
	return decodeHelixLength1(helix(i1-offset1, i2-offset2).second);
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerStackingOnly::
getHelixLength2(const size_t i1, const size_t i2) const
{
	return decodeHelixLength2(helix(i1-offset1, i2-offset2).second);
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerStackingOnly::
getHelixSeedLength1(const size_t i1, const size_t i2) const
{
	return decodeHelixSeedLength1(helixSeed(i1-offset1, i2-offset2).second);
}

////////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerStackingOnly::
getHelixSeedLength2(const size_t i1, const size_t i2) const
{
	return decodeHelixSeedLength2(helixSeed(i1-offset1, i2-offset2).second);
}

///////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerStackingOnly::
encodeHelixLength( const size_t l1, const size_t l2 ) const
{
	return l1 + l2*(helixConstraint.getMaxLength1()+1);
}

//////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerStackingOnly::
decodeHelixLength1( const size_t code ) const
{
	return code % (helixConstraint.getMaxLength1()+1);
}

//////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerStackingOnly::
decodeHelixLength2( const size_t code ) const
{
	return code / (helixConstraint.getMaxLength1()+1);
}

///////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerStackingOnly::
encodeHelixSeedLength( const size_t l1, const size_t l2 ) const
{
	return l1 + l2*(helixConstraint.getMaxLength1() + seedHandler->getConstraint().getMaxLength1()+1);
}

//////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerStackingOnly::
decodeHelixSeedLength1( const size_t code ) const
{
	return code % (helixConstraint.getMaxLength1() + seedHandler->getConstraint().getMaxLength1()+1);
}

//////////////////////////////////////////////////////////////////////////

inline
size_t
HelixHandlerStackingOnly::
decodeHelixSeedLength2( const size_t code ) const
{
	return code / (helixConstraint.getMaxLength1() + seedHandler->getConstraint().getMaxLength1()+1);
}

//////////////////////////////////////////////////////////////////////////

inline
void
HelixHandlerStackingOnly::setSeedHandler(SeedHandler *const seedHandler) {
	this->seedHandler = seedHandler;
}

//////////////////////////////////////////////////////////////////////////

} // namespace
#endif /* HELIXHANDLERSTACKINGONLY_H_ */
