
#include "IntaRNA/OutputHandlerCsv.h"

#if INTARNA_MULITHREADING
	#include <omp.h>
#endif

#include <boost/algorithm/string.hpp>

namespace IntaRNA {

////////////////////////////////////////////////////////////////////////

std::map<OutputHandlerCsv::ColType,std::string> OutputHandlerCsv::colType2string;

////////////////////////////////////////////////////////////////////////

OutputHandlerCsv::OutputHandlerCsv(
		  std::ostream & out
		, const InteractionEnergy & energy
		, const ColTypeList colOrder
		, const std::string& colSep
		, const bool printHeader
		)
 :	out(out)
	, energy(energy)
	, colOrder(colOrder)
	, colSep(colSep)
{
	// init mapping of coltypes to string
	initColType2string();

	// print CSV header of column names
	if (printHeader) {
		// ensure outputs do not intervene
#if INTARNA_MULITHREADING
		#pragma omp critical(intarna_omp_outputStreamUpdate)
#endif
		{
			out <<getHeader(colOrder,colSep);
		}
	}
}

////////////////////////////////////////////////////////////////////////

OutputHandlerCsv::~OutputHandlerCsv()
{
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_outputStreamUpdate)
#endif
	{
		// force output
		out.flush();
	}
}

////////////////////////////////////////////////////////////////////////

void
OutputHandlerCsv::
add( const Interaction & i )
{
#if INTARNA_IN_DEBUG_MODE
	// debug checks
	if ( i.basePairs.size() > 0 && ! i.isValid() ) {
		throw std::runtime_error("OutputHandlerCsv::add() : given interaction is not valid : "+toString(i));
	}
#endif

	// special handling if no base pairs present
	if (i.basePairs.size() == 0) {
		return;
	}

	// get interaction start/end per sequence
	const size_t i1 = i.basePairs.begin()->first;
	const size_t j1 = i.basePairs.rbegin()->first;
	const size_t i2 = i.basePairs.begin()->second;
	const size_t j2 = i.basePairs.rbegin()->second;

	// get individual energy contributions
	InteractionEnergy::EnergyContributions contr = energy.getE_contributions(i);

	// ensure outputs do not intervene
	{
		std::stringstream outTmp;

		for (auto col = colOrder.begin(); col != colOrder.end(); col++) {
			// print separator if needed
			if (col != colOrder.begin()) {
				outTmp <<colSep;
			}
			// print this column information
			switch ( *col ) {

			case id1:
				// ensure no colSeps are contained
				outTmp <<boost::replace_all_copy(energy.getAccessibility1().getSequence().getId(), colSep, "_");
				break;

			case id2:
				// ensure no colSeps are contained
				outTmp <<boost::replace_all_copy(energy.getAccessibility2().getSequence().getId(), colSep, "_");
				break;

			case seq1:
				outTmp <<energy.getAccessibility1().getSequence().asString();
				break;

			case seq2:
				outTmp <<energy.getAccessibility2().getAccessibilityOrigin().getSequence().asString();
				break;

			case subseq1:
				outTmp <<energy.getAccessibility1().getSequence().asString().substr(i1, j1-i1+1);
				break;

			case subseq2:
				outTmp <<energy.getAccessibility2().getAccessibilityOrigin().getSequence().asString().substr(j2, i2-j2+1);
				break;

			case subseqDP:
				outTmp <<energy.getAccessibility1().getSequence().asString().substr(i1, j1-i1+1)
					<<'&'
					<<energy.getAccessibility2().getAccessibilityOrigin().getSequence().asString().substr(j2, i2-j2+1);
				break;

			case subseqDB:
				outTmp <<(i1+1)
					<<energy.getAccessibility1().getSequence().asString().substr(i1, j1-i1+1)
					<<'&'
					<<(j2+1)
					<<energy.getAccessibility2().getAccessibilityOrigin().getSequence().asString().substr(j2, i2-j2+1);
				break;

			case start1:
				outTmp <<(i1+1);
				break;

			case end1:
				outTmp <<(j1+1);
				break;

			case start2:
				outTmp <<(j2+1);
				break;

			case end2:
				outTmp <<(i2+1);
				break;

			case hybridDP:
				outTmp <<Interaction::dotBracket( i );
				break;

			case hybridDB:
				outTmp <<Interaction::dotBar( i );
				break;

			case E:
				outTmp <<i.energy;
				break;

			case ED1:
				outTmp <<contr.ED1;
				break;

			case ED2:
				outTmp <<contr.ED2;
				break;

			case Pu1:
				outTmp <<std::exp( - contr.ED1 / energy.getRT() );
				break;

			case Pu2:
				outTmp <<std::exp( - contr.ED2 / energy.getRT() );
				break;

			case E_init:
				outTmp <<contr.init;
				break;

			case E_loops:
				outTmp <<contr.loops;
				break;

			case E_dangleL:
				outTmp <<contr.dangleLeft;
				break;

			case E_dangleR:
				outTmp <<contr.dangleRight;
				break;

			case E_endL:
				outTmp <<contr.endLeft;
				break;

			case E_endR:
				outTmp <<contr.endRight;
				break;

			case E_hybrid:
				outTmp <<(i.energy - contr.ED1 - contr.ED2);
				break;

			case E_norm:
				outTmp <<(i.energy / std::log( energy.size1() * energy.size2() ) );
				break;

			case E_hybridNorm:
				outTmp <<( (i.energy - contr.ED1 - contr.ED2) / std::log( energy.size1() * energy.size2() ) );
				break;

			case seedStart1:
				if (i.seed == NULL) {
					outTmp <<std::numeric_limits<E_type>::signaling_NaN();
				} else {
					outTmp <<(i.seed->bp_i.first+1);
				}
				break;

			case seedEnd1:
				if (i.seed == NULL) {
					outTmp <<std::numeric_limits<E_type>::signaling_NaN();
				} else {
					outTmp <<(i.seed->bp_j.first+1);
				}
				break;

			case seedStart2:
				if (i.seed == NULL) {
					outTmp <<std::numeric_limits<E_type>::signaling_NaN();
				} else {
					outTmp <<(i.seed->bp_j.second+1);
				}
				break;

			case seedEnd2:
				if (i.seed == NULL) {
					outTmp <<std::numeric_limits<E_type>::signaling_NaN();
				} else {
					outTmp <<(i.seed->bp_i.second+1);
				}
				break;

			case seedE:
				if (i.seed == NULL) {
					outTmp <<std::numeric_limits<E_type>::signaling_NaN();
				} else {
					outTmp <<i.seed->energy;
				}
				break;

			case seedED1:
				if (i.seed == NULL) {
					outTmp <<std::numeric_limits<E_type>::signaling_NaN();
				} else {
					outTmp <<energy.getED1( i.seed->bp_i.first, i.seed->bp_j.first );
				}
				break;

			case seedED2:
				if (i.seed == NULL) {
					outTmp <<std::numeric_limits<E_type>::signaling_NaN();
				} else {
					outTmp <<energy.getAccessibility2().getAccessibilityOrigin().getED( i.seed->bp_j.second, i.seed->bp_i.second );
				}
				break;

			case seedPu1:
				if (i.seed == NULL) {
					outTmp <<std::numeric_limits<E_type>::signaling_NaN();
				} else {
					outTmp <<std::exp( - energy.getED1( i.seed->bp_i.first, i.seed->bp_j.first ) / energy.getRT() );
				}
				break;

			case seedPu2:
				if (i.seed == NULL) {
					outTmp <<std::numeric_limits<E_type>::signaling_NaN();
				} else {
					outTmp <<std::exp( - energy.getAccessibility2().getAccessibilityOrigin().getED( i.seed->bp_j.second, i.seed->bp_i.second ) / energy.getRT() );
				}
				break;

			default : throw std::runtime_error("OutputHandlerCsv::add() : unhandled ColType '"+colType2string[*col]+"'");
			}
		}
		outTmp <<'\n';
	#if INTARNA_MULITHREADING
		#pragma omp critical(intarna_omp_outputStreamUpdate)
	#endif
		{
			out << outTmp.str();
		} // omp critical(intarna_omp_outputStreamUpdate)
	}

}

////////////////////////////////////////////////////////////////////////

} // namespace

