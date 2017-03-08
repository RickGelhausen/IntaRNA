
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
	// force output
	out.flush();
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
#if INTARNA_MULITHREADING
	#pragma omp critical(intarna_omp_outputStreamUpdate)
#endif
	{

		for (auto col = colOrder.begin(); col != colOrder.end(); col++) {
			// print separator if needed
			if (col != colOrder.begin()) {
				out <<colSep;
			}
			// print this column information
			switch ( *col ) {

			case id1:
				// ensure no colSeps are contained
				out <<boost::replace_all_copy(energy.getAccessibility1().getSequence().getId(), colSep, "_");
				break;

			case id2:
				// ensure no colSeps are contained
				out <<boost::replace_all_copy(energy.getAccessibility2().getSequence().getId(), colSep, "_");
				break;

			case seq1:
				out <<energy.getAccessibility1().getSequence().asString();
				break;

			case seq2:
				out <<energy.getAccessibility2().getAccessibilityOrigin().getSequence().asString();
				break;

			case subseq1:
				out <<energy.getAccessibility1().getSequence().asString().substr(i1, j1-i1+1);
				break;

			case subseq2:
				out <<energy.getAccessibility2().getAccessibilityOrigin().getSequence().asString().substr(j2, i2-j2+1);
				break;

			case subseqDP:
				out <<energy.getAccessibility1().getSequence().asString().substr(i1, j1-i1+1)
					<<'&'
					<<energy.getAccessibility2().getAccessibilityOrigin().getSequence().asString().substr(j2, i2-j2+1);
				break;

			case subseqDB:
				out <<(i1+1)
					<<energy.getAccessibility1().getSequence().asString().substr(i1, j1-i1+1)
					<<'&'
					<<(j2+1)
					<<energy.getAccessibility2().getAccessibilityOrigin().getSequence().asString().substr(j2, i2-j2+1);
				break;

			case start1:
				out <<(i1+1);
				break;

			case end1:
				out <<(j1+1);
				break;

			case start2:
				out <<(j2+1);
				break;

			case end2:
				out <<(i2+1);
				break;

			case hybridDP:
				out <<Interaction::dotBracket( i );
				break;

			case hybridDB:
				out <<Interaction::dotBar( i );
				break;

			case E:
				out <<i.energy;
				break;

			case ED1:
				out <<contr.ED1;
				break;

			case ED2:
				out <<contr.ED2;
				break;

			case Pu1:
				out <<std::exp( - contr.ED1 / energy.getRT() );
				break;

			case Pu2:
				out <<std::exp( - contr.ED2 / energy.getRT() );
				break;

			case E_init:
				out <<contr.init;
				break;

			case E_loops:
			{
				E_type e_loops = contr.loops;
				if (i.gap != NULL) {
					e_loops -= i.gap->energy;
				}
				out <<e_loops;
			}
				break;

			case E_dangleL:
				out <<contr.dangleLeft;
				break;

			case E_dangleR:
				out <<contr.dangleRight;
				break;

			case E_endL:
				out <<contr.endLeft;
				break;

			case E_endR:
				out <<contr.endRight;
				break;

			case E_gaps:
				if (i.gap == NULL) {
					out <<std::numeric_limits<E_type>::signaling_NaN();
				} else {
					out <<i.gap->energy;
				}
				break;

			case seedStart1:
				if (i.seed == NULL) {
					out <<std::numeric_limits<E_type>::signaling_NaN();
				} else {
					out <<(i.seed->bp_i.first+1);
					// print info for remaining (multi-side) seeds
					if (i.gap != NULL) {
						for (auto s = i.gap->seeds.rbegin(); s!=i.gap->seeds.rend(); s++) {
							out <<','<<(s->bp_i.first+1);
						}
					}
				}
				break;

			case seedEnd1:
				if (i.seed == NULL) {
					out <<std::numeric_limits<E_type>::signaling_NaN();
				} else {
					out <<(i.seed->bp_j.first+1);
					// print info for remaining (multi-side) seeds
					if (i.gap != NULL) {
						for (auto s = i.gap->seeds.rbegin(); s!=i.gap->seeds.rend(); s++) {
							out <<','<<(s->bp_j.first+1);
						}
					}
				}
				break;

			case seedStart2:
				if (i.seed == NULL) {
					out <<std::numeric_limits<E_type>::signaling_NaN();
				} else {
					out <<(i.seed->bp_j.second+1);
					// print info for remaining (multi-side) seeds
					if (i.gap != NULL) {
						for (auto s = i.gap->seeds.rbegin(); s!=i.gap->seeds.rend(); s++) {
							out <<','<<(s->bp_j.second+1);
						}
					}
				}
				break;

			case seedEnd2:
				if (i.seed == NULL) {
					out <<std::numeric_limits<E_type>::signaling_NaN();
				} else {
					out <<(i.seed->bp_i.second+1);
					// print info for remaining (multi-side) seeds
					if (i.gap != NULL) {
						for (auto s = i.gap->seeds.rbegin(); s!=i.gap->seeds.rend(); s++) {
							out <<','<<(s->bp_i.second+1);
						}
					}
				}
				break;

			case seedE:
				if (i.seed == NULL) {
					out <<std::numeric_limits<E_type>::signaling_NaN();
				} else {
					out <<i.seed->energy;
					// print info for remaining (multi-side) seeds
					if (i.gap != NULL) {
						for (auto s = i.gap->seeds.rbegin(); s!=i.gap->seeds.rend(); s++) {
							out <<','<<(s->energy);
						}
					}
				}
				break;

			case seedED1:
				if (i.seed == NULL) {
					out <<std::numeric_limits<E_type>::signaling_NaN();
				} else {
					out <<energy.getED1( i.seed->bp_i.first, i.seed->bp_j.first );
					// print info for remaining (multi-side) seeds
					if (i.gap != NULL) {
						for (auto s = i.gap->seeds.rbegin(); s!=i.gap->seeds.rend(); s++) {
							out <<','<<energy.getED1( s->bp_i.first, s->bp_j.first );
						}
					}
				}
				break;

			case seedED2:
				if (i.seed == NULL) {
					out <<std::numeric_limits<E_type>::signaling_NaN();
				} else {
					out <<energy.getAccessibility2().getAccessibilityOrigin().getED( i.seed->bp_j.second, i.seed->bp_i.second );
					// print info for remaining (multi-side) seeds
					if (i.gap != NULL) {
						for (auto s = i.gap->seeds.rbegin(); s!=i.gap->seeds.rend(); s++) {
							out <<','<<energy.getAccessibility2().getAccessibilityOrigin().getED( s->bp_j.second, s->bp_i.second );
						}
					}
				}
				break;

			case seedPu1:
				if (i.seed == NULL) {
					out <<std::numeric_limits<E_type>::signaling_NaN();
				} else {
					out <<std::exp( - energy.getED1( i.seed->bp_i.first, i.seed->bp_j.first ) / energy.getRT() );
					// print info for remaining (multi-side) seeds
					if (i.gap != NULL) {
						for (auto s = i.gap->seeds.rbegin(); s!=i.gap->seeds.rend(); s++) {
							out <<','<<std::exp( - energy.getED1( s->bp_i.first, s->bp_j.first ) / energy.getRT() );
						}
					}
				}
				break;

			case seedPu2:
				if (i.seed == NULL) {
					out <<std::numeric_limits<E_type>::signaling_NaN();
				} else {
					out <<std::exp( - energy.getAccessibility2().getAccessibilityOrigin().getED( i.seed->bp_j.second, i.seed->bp_i.second ) / energy.getRT() );
					// print info for remaining (multi-side) seeds
					if (i.gap != NULL) {
						for (auto s = i.gap->seeds.rbegin(); s!=i.gap->seeds.rend(); s++) {
							out <<','<<std::exp( - energy.getAccessibility2().getAccessibilityOrigin().getED( s->bp_j.second, s->bp_i.second ) / energy.getRT() );
						}
					}
				}
				break;

			case gaps1:
				if (i.gap != NULL) {
					out <<i.gap->gaps1.shift(1,i.s1->size());
				}
				break;

			case gaps2:
				if (i.gap != NULL) {
					out <<i.gap->gaps2.shift(1,i.s2->size());
				}
				break;

			default : throw std::runtime_error("OutputHandlerCsv::add() : unhandled ColType '"+colType2string[*col]+"'");
			}
		}
		out <<'\n';
	} // omp critical(intarna_omp_outputStreamUpdate)

}

////////////////////////////////////////////////////////////////////////

} // namespace

