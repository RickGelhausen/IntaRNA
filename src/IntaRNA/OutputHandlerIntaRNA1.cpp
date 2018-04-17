
#include "IntaRNA/OutputHandlerIntaRNA1.h"

#include <sstream>
#include <iomanip>

#if INTARNA_MULITHREADING
	#include <omp.h>
#endif

namespace IntaRNA {

////////////////////////////////////////////////////////////////////////////

OutputHandlerIntaRNA1::
OutputHandlerIntaRNA1(
		std::ostream & out
		, const InteractionEnergy & energy
		, const bool detailedOutput
		)
 :
	detailedOutput(detailedOutput)
	, out(out)
	, energy(energy)
	, initialOutputDone(false)
	, printSeparator(false)
{
}

////////////////////////////////////////////////////////////////////////////

OutputHandlerIntaRNA1::
~OutputHandlerIntaRNA1()
{
	out.flush();
}

////////////////////////////////////////////////////////////////////////////

void
OutputHandlerIntaRNA1::
add( const Interaction & i )
{
#if INTARNA_IN_DEBUG_MODE
	// debug checks
	if ( i.basePairs.size() > 0 && ! i.isValid() ) {
		throw std::runtime_error("OutputHandlerIntaRNA1::add() : given interaction is not valid : "+toString(i));
	}
#endif

	// special handling if no base pairs present
	if (i.basePairs.size() == 0) {
		// no output for empty interactions
		return;
	}

	// count report
	reportedInteractions++;

	if (!initialOutputDone) {
		// ensure outputs do not intervene
		std::stringstream outTmp;
		if (printSeparator) {
			outTmp <<"\n=========================\n"
				<<'\n';


		}
		// write sequences in FASTA to outTmp
		outTmp
		<<">" <<energy.getAccessibility1().getSequence().getId() <<"\n"
		<<(detailedOutput ? energy.getAccessibility1().getSequence().asString()+"\n" : "")
		<<">" <<energy.getAccessibility2().getSequence().getId() <<"\n"
		<<(detailedOutput ? energy.getAccessibility2().getAccessibilityOrigin().getSequence().asString()+"\n" : "")
		;
		initialOutputDone = true;

	#if INTARNA_MULITHREADING
		#pragma omp critical(intarna_omp_outputStreamUpdate)
	#endif
		{
			out <<outTmp.str();
		} // omp critical(intarna_omp_outputStreamUpdate)
	}

	// get interaction start/end per sequence
	const size_t i1 = i.basePairs.begin()->first;
	const size_t j1 = i.basePairs.rbegin()->first;
	const size_t i2 = i.basePairs.begin()->second;
	const size_t j2 = i.basePairs.rbegin()->second;

	// determine maximal length of sequence to the left of first base pair
	const size_t flankingLength = 3 + std::max(i1,i.s2->size()-i2-1);

	// decompose printing into several rows
	std::ostringstream s1Unbound;
	std::ostringstream s1Bound;
	std::ostringstream s2Bound;
	std::ostringstream s2Unbound;

	// left flanking
	// unbound region s1
	s1Unbound.width(flankingLength);
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
	s1Bound.width(flankingLength); s1Bound <<' ';
	s2Bound.width(flankingLength); s2Bound <<' ';
	// unbound region s2
	s2Unbound.width(flankingLength);
	s2Unbound <<std::right;
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
			// unbound region s1
			s1Unbound.width(loop);
			if (loop1 > 0) {
				s1Unbound <<i.s1->asString().substr( leftBP->first +1, loop1 );
			} else {
				s1Unbound <<' ';
			}
			// bound region
			s1Bound.width(loop); s1Bound <<' ';
			s2Bound.width(loop); s2Bound <<' ';
			// unbound region s2
			s2Unbound.width(loop);
			if (loop2 > 0) {
				s2Unbound <<reverse(i.s2->asString().substr( curBP->second +1, loop2 ));
			} else {
				s2Unbound <<' ';
			}
		}
		interactionLength += loop;

		// print current base pair (right end of internal loop)
		nt1 = i.s1->asString().at(curBP->first);
		nt2 = i.s2->asString().at(curBP->second);
		s1Unbound.width(1);	 s1Unbound <<' ';
		s1Bound.width(1);	 s1Bound   <<nt1;
		s2Bound.width(1);	 s2Bound   <<nt2;
		s2Unbound.width(1);	 s2Unbound <<' ';
		interactionLength++;
	}

	// flanking right
	// unbound region s1
	// add remaining sequence
	s1Unbound <<i.s1->asString().substr(j1+1)
			<<"-3'";
	// full sequence prefix
	s2Unbound <<reverse(i.s2->asString().substr( 0, j2))
			<<"-5'";

	{
		// ensure outputs do not intervene
		std::stringstream outTmp;
		// print full interaction to output stream
		outTmp <<'\n'
			// print collected interaction stuff
			<<s1Unbound.str() <<'\n'
			<<s1Bound.str() <<'\n'
			<<s2Bound.str() <<'\n'
			<<s2Unbound.str() <<'\n'
			;

		if (detailedOutput) {
			// get individual energy contributions
			InteractionEnergy::EnergyContributions contr = energy.getE_contributions(i);
				// print interaction details
			outTmp	<<'\n'
				<<"positions(target)     : "<<(i.basePairs.begin()->first +1)<<" -- "<<(i.basePairs.rbegin()->first +1) <<'\n'
				<<"positions seed(target): "<<(i.seed!=NULL?toString(i.seed->bp_i.first +1):"?")<<" -- "<<(i.seed!=NULL?toString(i.seed->bp_j.first +1):"?") <<'\n'
				<<"positions with dangle(target): "<<(i.basePairs.begin()->first +1)<<" -- "<<(i.basePairs.rbegin()->first +1) <<'\n'
				<<"positions(ncRNA)      : "<<(i.basePairs.rbegin()->second +1)<<" -- "<<(i.basePairs.begin()->second +1) <<'\n'
				<<"positions seed(ncRNA) : "<<(i.seed!=NULL?toString(i.seed->bp_j.second +1):"?")<<" -- "<<(i.seed!=NULL?toString(i.seed->bp_i.second +1):"?") <<'\n'
				<<"positions with dangle(ncRNA): "<<(i.basePairs.rbegin()->second +1)<<" -- "<<(i.basePairs.begin()->second +1) <<'\n'
				<<"ED target need: "<<contr.ED1 <<" kcal/mol"<<'\n'
				<<"ED ncRNA  need: "<<contr.ED2 <<" kcal/mol"<<'\n'
				<<"hybrid energy : "<<(i.energy-contr.ED1-contr.ED2) <<" kcal/mol"<<'\n'
				<<"\n"
				<<"energy: "<<i.energy <<" kcal/mol\n"
				;
		} else {
			// normal minimal information output
			outTmp	<<'\n'
				<<"energy: "<<i.energy <<" kcal/mol\n"
				;
		}

		// ensure outputs do not intervene
	#if INTARNA_MULITHREADING
		#pragma omp critical(intarna_omp_outputStreamUpdate)
	#endif
		{
			out <<outTmp.str();
		} // omp critical(intarna_omp_outputStreamUpdate)
	}

}


////////////////////////////////////////////////////////////////////////////

void
OutputHandlerIntaRNA1::
add( const InteractionRange & range )
{
	add( Interaction(range) );
}

////////////////////////////////////////////////////////////////////////////

void
OutputHandlerIntaRNA1::
addSeparator (const bool yesNo )
{
	printSeparator = yesNo;
	if (initialOutputDone) {
#if INTARNA_MULITHREADING
		#pragma omp critical(intarna_omp_logOutput)
#endif
		{ LOG(INFO) <<"OutputHandlerIntaRNA1::addSeparator() called but initial output already done..."; }
	}
}

////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////

} // namespace

