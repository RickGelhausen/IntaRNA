
#include "IntaRNA/OutputHandlerText.h"

#include <sstream>
#include <iomanip>

#if INTARNA_MULITHREADING
	#include <omp.h>
#endif

namespace IntaRNA {

////////////////////////////////////////////////////////////////////////////

OutputHandlerText::
OutputHandlerText(
		std::ostream & out,
		const InteractionEnergy & energy,
		const size_t flankingLength_,
		const bool detailedOutput
		)
 :
	out(out)
	, energy(energy)
	, flankingLength(flankingLength_)
	, detailedOutput(detailedOutput)
{
	if (flankingLength < 10) {
		throw std::runtime_error("OutputHandlerText::OutputHandlerText : flankingLength ("+toString(flankingLength_)+") has to be at least 9");
	}
}

////////////////////////////////////////////////////////////////////////////

OutputHandlerText::
~OutputHandlerText()
{
	out.flush();
}

////////////////////////////////////////////////////////////////////////////

void
OutputHandlerText::
add( const Interaction & i )
{
#if INTARNA_IN_DEBUG_MODE
	// debug checks
	if ( i.basePairs.size() > 0 && ! i.isValid() ) {
		throw std::runtime_error("OutputHandlerText::add() : given interaction is not valid : "+toString(i));
	}
#endif

	// special handling if no base pairs present
	if (i.basePairs.size() == 0) {
		// ensure outputs do not intervene
#if INTARNA_MULITHREADING
		#pragma omp critical(intarna_omp_outputStreamUpdate)
#endif
		{
			out <<"\n"
				<<"no favorable interaction for "<<i.s1->getId() <<" and "<<i.s2->getId()<<std::endl;
		} // omp critical(intarna_omp_outputStreamUpdate)
		return;
	}

	// get interaction start/end per sequence
	const size_t i1 = i.basePairs.begin()->first;
	const size_t j1 = i.basePairs.rbegin()->first;
	const size_t i2 = i.basePairs.begin()->second;
	const size_t j2 = i.basePairs.rbegin()->second;

	// decompose printing into several rows
	std::ostringstream pos1tag;
	std::ostringstream s1Unbound;
	std::ostringstream s1Bound;
	std::ostringstream pairing;
	std::ostringstream s2Bound;
	std::ostringstream s2Unbound;
	std::ostringstream pos2tag;

	// left flanking
	pos1tag	<<std::setw(flankingLength+3+1) <<'|';
	pos2tag	<<std::setw(flankingLength+3+1) <<'|';
	// unbound region s1
	s1Unbound.width(flankingLength+3);
	s1Unbound <<std::right;
	if (i1 < flankingLength) {
		// full sequence prefix
		s1Unbound <<("5'-" + i.s1->asString().substr( 0, i1));
	} else {
		// prefix shortening needed
		s1Unbound <<("5'-"
					+ i.s1->asString().substr( 0,3)
					+ "..."
					+ i.s1->asString().substr( (size_t)std::max(0,(int)i1+6-(int)flankingLength), flankingLength-6)
					);
	}
	// bound region
	s1Bound.width(flankingLength+3); s1Bound <<' ';
	pairing.width(flankingLength+3); pairing <<' ';
	s2Bound.width(flankingLength+3); s2Bound <<' ';
	// unbound region s2
	s2Unbound.width(flankingLength+3);
	s2Unbound <<std::right; // <<("3'-" + i.s2->asString().substr( (size_t)std::max(0,(int)i2-(int)flankingLength), std::min(i2,flankingLength) ));
//	if (i2 < flankingLength) {
//		// full sequence prefix
//		s2Unbound <<("3'-" + i.s2->asString().substr( 0, i2));
//	} else {
//		// prefix shortening needed
//		s2Unbound <<("3'-"
//					+ i.s2->asString().substr( 0,3)
//					+ "..."
//					+ i.s2->asString().substr( (size_t)std::max(0,(int)i2+6-(int)flankingLength), flankingLength-6)
//					);
//	}
	if (i2+flankingLength > i.s2->size()) {
		// add remaining sequence
		s2Unbound <<("3'-" + reverse(i.s2->asString().substr(i2+1)));
	} else {
		// shortening needed
		s2Unbound <<("3'-"
					+ reverse(i.s2->asString().substr(i.s2->size()-3))
					+ "..."
					+ reverse(i.s2->asString().substr(i2+1, flankingLength-6))
					);
	}

	// interaction start left
	Interaction::PairingVec::const_iterator leftBP = i.basePairs.begin();
	char nt1 = i.s1->asString().at(leftBP->first);
	char nt2 = i.s2->asString().at(leftBP->second);
	// print start
	s1Unbound.width(1);	 s1Unbound <<std::left <<' ';
	s1Bound.width(1);	 s1Bound   <<std::left <<i.s1->asString().at(leftBP->first);
	pairing.width(1);	 pairing   <<std::left;
	if ((nt1=='G' && nt2=='U') || (nt1=='U' && nt2=='G')) {
		pairing   <<':';
	} else {
		pairing   <<'|';
	}
	s2Bound.width(1);	 s2Bound   <<std::left <<i.s2->asString().at(leftBP->second);
	s2Unbound.width(1);	 s2Unbound <<std::left <<' ';

	// iterate loops in interaction region
	Interaction::PairingVec::const_iterator curBP = i.basePairs.begin();
	size_t loop1=0, loop2=0, loop=0, interactionLength = 1;
	for (++curBP; curBP != i.basePairs.end(); ++curBP, ++leftBP) {
		// handle internal loop region
		// get specific loop lengths
		loop1 = curBP->first - leftBP->first -1;
		loop2 = leftBP->second - curBP->second -1;
		loop = std::max(loop1,loop2);
		// print unbound loop regions
		if (loop>0) {
			// gap handling in s1
			if (i.gap != NULL && i.gap->gaps1.covers(leftBP->first +1)) {
				// TODO replace with dot-bracket of subseq's mfe
				// mark gap
				const int leftBrackets = (int)loop/2;
				pos1tag.width(leftBrackets); pos1tag <<'<';
				pos1tag.width(loop-leftBrackets); pos1tag <<'>';
			} else {
				pos1tag.width(loop); pos1tag <<' ';
			}
			// unbound region s1
			s1Unbound.width(loop);
			if (loop1 > 0) {
				s1Unbound <<i.s1->asString().substr( leftBP->first +1, loop1 );
			} else {
				s1Unbound <<' ';
			}
			// bound region
			s1Bound.width(loop); s1Bound <<' ';
			pairing.width(loop); pairing <<' ';
			s2Bound.width(loop); s2Bound <<' ';
			// unbound region s2
			s2Unbound.width(loop);
			if (loop2 > 0) {
				s2Unbound <<reverse(i.s2->asString().substr( curBP->second +1, loop2 ));
			} else {
				s2Unbound <<' ';
			}
			// gap handling in s2
			if (i.gap != NULL && i.gap->gaps2.covers(leftBP->second +1)) {
				// TODO replace with dot-bracket of subseq's mfe
				// mark gap
				const int leftBrackets = (int)loop/2;
				pos2tag.width(leftBrackets); pos2tag <<'<';
				pos2tag.width(loop-leftBrackets); pos2tag <<'>';
			} else {
				pos2tag.width(loop); pos2tag <<' ';
			}
		}
		interactionLength += loop;

		// print current base pair (right end of internal loop)
		nt1 = i.s1->asString().at(curBP->first);
		nt2 = i.s2->asString().at(curBP->second);
		pos1tag << (curBP->first == i.basePairs.rbegin()->first ? "|" : " "); // interaction end tag
		s1Unbound.width(1);	 s1Unbound <<' ';
		s1Bound.width(1);	 s1Bound   <<nt1;
		pairing.width(1);
		if ((nt1=='G' && nt2=='U') || (nt1=='U' && nt2=='G')) {
			pairing   <<':';
		} else {
			pairing   <<'|';
		}
		s2Bound.width(1);	 s2Bound   <<nt2;
		s2Unbound.width(1);	 s2Unbound <<' ';
		pos2tag << (curBP->second == i.basePairs.rbegin()->second ? "|" : " "); // interaction end tag
		interactionLength++;
	}

	// flanking right
	// unbound region s1
	if (j1+flankingLength > i.s1->size()) {
		// add remaining sequence
		s1Unbound <<i.s1->asString().substr(j1+1);
	} else {
		// shortening needed
		s1Unbound <<i.s1->asString().substr(j1+1, flankingLength-6)
					<<"..."
					<<i.s1->asString().substr(i.s1->size()-3);
	}
	s1Unbound <<"-3'";
	// unbound region s2
//	if (j2+flankingLength > i.s2->size()) {
//		// add remaining sequence
//		s2Unbound <<i.s2->asString().substr(j2+1);
//	} else {
//		// shortening needed
//		s2Unbound <<i.s2->asString().substr(j2+1, flankingLength-6)
//					<<"..."
//					<<i.s2->asString().substr(i.s2->size()-3);
//	}
	if (j2 < flankingLength) {
		// full sequence prefix
		s2Unbound <<reverse(i.s2->asString().substr( 0, j2));
	} else {
		// prefix shortening needed
		s2Unbound 	<<reverse(i.s2->asString().substr( (size_t)std::max(0,(int)j2+6-(int)flankingLength), flankingLength-6))
					<<"..."
					<<reverse(i.s2->asString().substr( 0,3))
					;
	}
	s2Unbound <<"-5'";

	// position information (first and last interacting positions in sequence)
	// start indexing with 1 instead of 0
	std::ostringstream pos1;
	pos1	<<std::setw(flankingLength+3+1) <<std::right <<(i1+1);
	size_t numberSize = numStringLength(i1);
	// put right (3') end only if not overlapping with 5' start position (left)
	if (i1+1 < j1 && numberSize+2 < flankingLength+3+interactionLength) {
		pos1 <<std::setw(interactionLength -2) <<' '
				<<std::setw(1) <<std::left <<(j1+1);
	}
	std::ostringstream pos2;
	numberSize = numStringLength(i.s2->size()-i2-1);
	// put left (3') end only if not overlapping with 5' start position (right)
	pos2	<<std::setw(flankingLength+3+1);
	if (i2==j2 || numberSize+2 < (flankingLength+3+interactionLength)) {
		pos2	 <<std::right <<(i2+1);
	} else {
		pos2	 <<' ';
	}
	if (i2 > j2) {
		pos2 <<std::setw(interactionLength - 2) <<' '
				<<std::setw(1) <<std::left <<(j2+1);
	}

	// get individual energy contributions
	InteractionEnergy::EnergyContributions contr = energy.getE_contributions(i);

	// ensure outputs do not intervene
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_outputStreamUpdate)
#endif
	{
		// print full interaction to output stream
		out <<'\n'
			// get ID of s1
			<<i.s1->getId() <<'\n'
			// get position in s1
			<<pos1.str() <<'\n'
			<<pos1tag.str() <<'\n'
			// print collected interaction stuff
			<<s1Unbound.str() <<'\n'
			<<s1Bound.str() <<'\n'
			<<pairing.str() <<'\n'
			<<s2Bound.str() <<'\n'
			<<s2Unbound.str() <<'\n'
			// get position in reversed index order for s2
			<<pos2tag.str() <<'\n'
			<<pos2.str() <<'\n'
			// get ID of s2
			<<i.s2->getId() <<'\n'
			;

		if (detailedOutput) {
			out
				// interaction range
				<<"\n"
				<<"interaction seq1   = "<<(i.basePairs.begin()->first +1)<<" -- "<<(i.basePairs.rbegin()->first +1) <<'\n'
				<<"interaction seq2   = "<<(i.basePairs.rbegin()->second +1)<<" -- "<<(i.basePairs.begin()->second +1) <<'\n'
				;
		} // detailed
			// print energy
		out
			<<"\n"
			<<"interaction energy = "<<i.energy <<" kcal/mol\n"
			;

		if (detailedOutput) {
			E_type Egaps = (i.gap == NULL) ? 0.0 : i.gap->energy;
			out
					<<"  = E(init)        = "<<contr.init<<'\n'
					<<"  + E(loops)       = "<<(contr.loops-Egaps)<<'\n'
					<<"  + E(dangleLeft)  = "<<contr.dangleLeft<<'\n'
					<<"  + E(dangleRight) = "<<contr.dangleRight<<'\n'
					<<"  + E(endLeft)     = "<<contr.endLeft<<'\n'
					<<"  + E(endRight)    = "<<contr.endRight<<'\n'
					<<"  + ED(seq1)       = "<<contr.ED1<<'\n'
					<<"    from Pu(seq1)  = "<<std::exp(-contr.ED1/energy.getRT())<<'\n'
					<<"  + ED(seq2)       = "<<contr.ED2<<'\n'
					<<"    from Pu(seq2)  = "<<std::exp(-contr.ED2/energy.getRT())<<'\n'
					;
			if (i.gap != NULL) {
				out
						<<"  + E(gaps)        = "<<Egaps<<'\n'
						;
			}

			// print seed information if available
			if (i.seed != NULL) {
				if (i.gap != NULL) {
					// print all seeds left of last gap
					for (auto is = i.gap->seeds.begin(); is!=i.gap->seeds.end(); is++) {
						out <<"\n";
						seedToStream( *(is) );
					}
				}
				// print seed right of last gap
				out <<"\n";
				seedToStream( *(i.seed) );
			} // seed

			// print multi-side gap information if available
			if (i.gap != NULL) {
				out
						<<"\n"
						<<"gaps seq1   = "<<i.gap->gaps1.shift(1,i.s1->size()) <<'\n'
						<<"gaps seq2   = "<<i.gap->gaps2.shift(1,i.s1->size()) <<'\n'
						;
			} // gap

		} // detailed
	} // omp critical(intarna_omp_outputStreamUpdate)

}


////////////////////////////////////////////////////////////////////////////

	void
	OutputHandlerText::
	seedToStream( const Interaction::Seed & seed )
	{
		out
				<<"seed seq1   = "<<(seed.bp_i.first +1)<<" -- "<<(seed.bp_j.first +1) <<'\n'
				<<"seed seq2   = "<<(seed.bp_j.second +1)<<" -- "<<(seed.bp_i.second +1) <<'\n'
				<<"seed energy = "<<(seed.energy)<<" kcal/mol\n"
				<<"seed ED1    = "<<energy.getED1( seed.bp_i.first, seed.bp_j.first )<<" kcal/mol\n"
				<<"seed ED2    = "<<energy.getAccessibility2().getAccessibilityOrigin().getED( seed.bp_j.second, seed.bp_i.second )<<" kcal/mol\n"
				<<"seed Pu1    = "<<std::exp(-(energy.getED1( seed.bp_i.first, seed.bp_j.first ))/energy.getRT())<<'\n'
				<<"seed Pu2    = "<<std::exp(-(energy.getAccessibility2().getAccessibilityOrigin().getED( seed.bp_j.second, seed.bp_i.second ))/energy.getRT())<<'\n'
				;

	}

////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////

} // namespace

