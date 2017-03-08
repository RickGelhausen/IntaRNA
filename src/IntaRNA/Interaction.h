
#ifndef INTARNA_INTERACTION_H_
#define INTARNA_INTERACTION_H_

#include <vector>
#include <utility>

#include "IntaRNA/general.h"
#include "IntaRNA/RnaSequence.h"
#include "IntaRNA/InteractionRange.h"
#include "IntaRNA/IndexRangeList.h"

namespace IntaRNA {

// dummy declaration to avoid inclusion loop
class InteractionRange;


/**
 * Data structure to describe and store an Interaction
 */
class Interaction {

public:

	//! type of a base pair index encoding
	typedef std::pair<size_t,size_t> BasePair;

	//! type of a vector encoding base pair indices that are interacting
	typedef std::vector<BasePair> PairingVec;

	/**
	 * Provides the seed information of an Interaction.
	 */
	class Seed {
	public:
		//! left-most base pair of seed
		BasePair bp_i;
		//! right-most base pair of seed
		BasePair bp_j;
		//! overall energy of seed
		E_type energy;

        //! default constructor
        Seed()
                : bp_i(), bp_j(), energy(0.0)
        {}
	};

	/**
	 * Provides information on multi-side gaps
	 */
	class Gap {
	public:
		//! index ranges of the intramolecularly structured gap in seq1
		IndexRangeList gaps1;
		//! index ranges of the intramolecularly structured gap in seq2
		IndexRangeList gaps2;
		//! summed energy terms that score all gaps
		E_type energy;
		//! list of all but the right most (in seq1)
		std::list< Seed > seeds;

		//! default constructor
		Gap()
				: gaps1(), gaps2(), energy(0.0), seeds()
		{}
	};
public:

	//! the first interaction partner
	const RnaSequence * s1;

	//! the second interaction partner
	const RnaSequence * s2;

	//! interacting indices
	PairingVec basePairs;

	//! energy of the interaction (can be NaN)
	E_type energy;

	//! optional: seed information
	Seed * seed;

	//! optional: multi-side gap information
	Gap * gap;

	/**
	 * construction
	 * @param s1 the sequence of the first interaction partner
	 * @param s2 the sequence of the second interaction partner
	 */
	Interaction( const RnaSequence & s1, const RnaSequence & s2 );

	/**
	 * copy construction from interaction
	 *
	 * @param toCopy the interaction range to make this a copy of
	 */
	Interaction( const Interaction & toCopy );

	/**
	 * construction from interaction range, ie. forming base pairs at the
	 * beginning and end of the interaction range.
	 *
	 * NOTE: the ends have to be complementary nucleotides
	 *
	 * @param range the interaction range to be used for initialization
	 */
	Interaction( const InteractionRange & range );

	/**
	 * destruction
	 */
	virtual ~Interaction();

	/**
	 * checks whether or not the current interaction is non-empty, nested and
	 * sorted.
	 * @return true if the interaction encoding is non-empty, nested and sorted;
	 *         false otherwise
	 */
	bool
	isValid() const;

	/**
	 * Checks whether or not the interaction contains no base pairs.
	 * @return true if the interaction encodes no base pair; false otherwise
	 */
	bool
	isEmpty() const;

	/**
	 * Sorts the stored interacting base pairs.
	 */
	void
	sort();

	/**
	 * Clears all interaction information.
	 */
	void
	clear();

	/**
	 * Sets the seedRange member according to the given data
	 * @param bp_left left most base pair in seed
	 * @param bp_right right most base pair in seed
	 * @param energy full energy of the seed interaction
	 */
	void
	setSeedRange( const BasePair bp_left, const BasePair bp_right, const E_type energy );


	/**
	 * Creates an interaction with one base pair for each interaction range
	 * boundary.
	 *
	 * NOTE: the ends have to be complementary nucleotides
	 *
	 * NOTE: the interaction energy is reset too.
	 *
	 * @param toCopy the interaction to copy
	 * @return the altered object (*this)
	 */
	Interaction &
	operator= ( const Interaction & toCopy );

	/**
	 * Creates an interaction with one base pair for each interaction range
	 * boundary.
	 *
	 * NOTE: the ends have to be complementary nucleotides
	 *
	 * NOTE: the interaction energy is reset too.
	 *
	 * @param range the interaction range to copy
	 * @return the altered object (*this)
	 */
	Interaction &
	operator= ( const InteractionRange & range );


	/**
	 * Compares if an interaction has larger energy that a given value
	 *
	 * @param energy the energy to compare to
	 * @param hasLargerE the interaction that is supposed to have larger energy
	 * @return (energy < (hasLargerE.energy-precisionDelta))
	 */
	static
	bool
	compareEnergy( const E_type& energy, const Interaction & hasLargerE );

	/**
	 * Prints the interacting base pairs to stream
	 * @param out the ostream to write to
	 * @param i the Interaction object to add
	 * @return the altered stream out
	 */
	friend std::ostream& operator<<(std::ostream& out, const Interaction& i);


	/**
	 * Produces a dot-bar encoding of the interaction in the form
	 *
	 *    startId1 dot-bar-1 & startId2 dot-bar-2
	 *
	 * where dot-bar-1/2 encodes all positions of the first/second sequence
	 * enclosed by the first and last base pair of the interaction.
	 * Here, a '|' bar encodes an intermolecular base pair and a '.' dot
	 * encodes unpaired positions.
	 *
	 * Note, within this encoding, the '&' separated parts can be swapped if
	 * needed yielding still a valid encoding (for the swapped sequences).
	 *
	 * @param i interaction to encode
	 * @return the dot-bracket encoding
	 */
	static
	std::string
	dotBar( const Interaction & i );

	/**
	 * Produces a dot-bracket encoding of the interaction in VRNA style in the
	 * form
	 *
	 *   dot-bracket-1 & dot-bracket-2
	 *
	 * where dot-bracket-1/2 encodes all positions of the first/second sequence
	 * enclosed by the first and last base pair of the interaction.
	 * In dot-bracket-1, a closing '(' parenthesis encodes an intermolecular
	 * base pair and a '.' dot encodes unpaired positions. In dot-bracket-2,
	 * closing ')' parentheses are used to mark base pairs.
	 *
	 * @param i interaction to encode
	 * @param symOpen the symbol to be used for base pairs in dot-bracket-1
	 * @param symClose the symbol to be used for base pairs in dot-bracket-2
	 *
	 * @return the VRNA-styled dot-bracket encoding
	 */
	static
	std::string
	dotBracket( const Interaction & i, const char symOpen = '(', const char symClose = ')');


protected:

	template < typename bpIterator >
	static
	std::string
	dotSomething( bpIterator bpBegin, const bpIterator bpEnd, const bool handleSeq1, const char bpSymbol );


};

/**
 * Prints the interacting base pair to stream
 * @param out the ostream to write to
 * @param bp the Interaction base pair object to add
 * @return the altered stream out
 */
std::ostream& operator<<(std::ostream& out, const Interaction::BasePair& bp);



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


inline
Interaction::Interaction( const RnaSequence & s1, const RnaSequence & s2 )
:
	s1(&s1)
	, s2(&s2)
	, basePairs()
	, energy( std::numeric_limits<E_type>::signaling_NaN() )
	, seed( NULL )
    , gap( NULL )
{
}

////////////////////////////////////////////////////////////////////////////

inline
Interaction::Interaction( const Interaction & toCopy )
:
    s1(toCopy.s1)
    , s2(toCopy.s2)
    , basePairs(toCopy.basePairs)
    , energy( toCopy.energy )
    , seed( toCopy.seed == NULL ? NULL : new Seed() )
    , gap( toCopy.gap == NULL ? NULL : new Gap() )
{
    // copy seed if needed
    if (seed != NULL) {
        (*seed) = *(toCopy.seed);
    }
    // copy gap if needed
    if (gap != NULL) {
        (*gap) = *(toCopy.gap);
    }
}

////////////////////////////////////////////////////////////////////////////

inline
Interaction::Interaction( const InteractionRange & range )
:
	s1(NULL)
	, s2(NULL)
	, basePairs()
	, energy( std::numeric_limits<E_type>::signaling_NaN() )
	, seed( NULL )
{
	// init data
	this->operator =( range );
}

////////////////////////////////////////////////////////////////////////////

inline
Interaction::~Interaction()
{
    INTARNA_CLEANUP(seed);
    INTARNA_CLEANUP(gap);
}

////////////////////////////////////////////////////////////////////////////

inline
bool
Interaction::
isEmpty() const
{
	return basePairs.size()==0;
}

////////////////////////////////////////////////////////////////////////////

inline
void
Interaction::
sort()
{
	// sort based on STL functionalities for vector and pair
	std::sort(basePairs.begin(), basePairs.end());
}

////////////////////////////////////////////////////////////////////////////

inline
void
Interaction::
clear()
{
	// clear interaction base pairing information
	basePairs.clear();
	// clear energy
	energy = std::numeric_limits<E_type>::signaling_NaN();
	// undo seed information
    INTARNA_CLEANUP(seed);
    INTARNA_CLEANUP(gap);

}

////////////////////////////////////////////////////////////////////////////

inline
bool
Interaction::
compareEnergy( const E_type & energy, const Interaction & hasLargerE )
{
	return energy < hasLargerE.energy && !(E_equal(energy,hasLargerE.energy));
}

////////////////////////////////////////////////////////////////////////////

inline
std::ostream&
operator<<(std::ostream& out, const Interaction::BasePair& bp)
{
	out <<"("<<bp.first<<"-"<<bp.second<<")";
	return out;
}

////////////////////////////////////////////////////////////////////////////

inline
std::ostream&
operator<<(std::ostream& out, const Interaction& i)
{
	for (int p=0; p<i.basePairs.size(); p++) {
		out <<(p==0?"":",") <<i.basePairs.at(p);
	}
	return out;
}

////////////////////////////////////////////////////////////////////////////

template < typename bpIterator >
inline
std::string
Interaction::
dotSomething( bpIterator bp, const bpIterator bpEnd, const bool handleSeq1, const char bpSymbol )
{
	// stream to compile output
	std::stringstream dotBracket;
	// compile dotBracket1
	for (bpIterator bpLast = bp; bp!=bpEnd; bp++) {
		// fill unpaired up to current bp in seq1
		for ( size_t u = handleSeq1 ? (bp->first - bpLast->first) : (bp->second - bpLast->second); u>1; u-- ) {
			dotBracket <<'.';
		}
		// add base pair
		dotBracket <<bpSymbol;
		// update last bp for next round
		bpLast = bp;
	}
	// return compiled string
	return dotBracket.str();
}

////////////////////////////////////////////////////////////////////////////

} // namespace


#endif /* INTERACTION_H_ */
