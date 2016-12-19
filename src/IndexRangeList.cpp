
#include "IndexRangeList.h"

#include <algorithm>

//////////////////////////////////////////////////////////////////////

const boost::regex IndexRangeList::regex("^(\\d|([123456789]\\d*)-(\\d|[123456789]\\d*),)*(\\d|[123456789]\\d*)-(\\d|[123456789]\\d*)$");

//////////////////////////////////////////////////////////////////////

IndexRangeList::IndexRangeList()
:
list()
{
}

//////////////////////////////////////////////////////////////////////

IndexRangeList::IndexRangeList( const std::string & stringEncoding )
:
list()
{
	fromString(stringEncoding);
}

//////////////////////////////////////////////////////////////////////

IndexRangeList::IndexRangeList( const IndexRangeList & toCopy )
:
list(toCopy.list)
{
}

//////////////////////////////////////////////////////////////////////

IndexRangeList::~IndexRangeList()
{
}

//////////////////////////////////////////////////////////////////////

bool
IndexRangeList::
covers( const size_t index ) const
{
	// quick check
	if (list.empty()) {
		return false;
	}

	// find first range that with begin > index
	const_iterator r = std::upper_bound( list.begin(), list.end(), IndexRange(index,std::numeric_limits<size_t>::max()) );
	if ( r == list.begin() ) {
		return false;
	} else {
		// go to preceding range and check if <= the end of the blocked range
		return index <= (--r)->to;
	}
}

//////////////////////////////////////////////////////////////////////

bool
IndexRangeList::
covers( const size_t from, const size_t to ) const
{
	return covers( IndexRange(from, to) );
}

//////////////////////////////////////////////////////////////////////

bool
IndexRangeList::
covers( const IndexRange & range ) const
{
	// quick check
	if (list.empty()) {
		return false;
	}

	// find first range that with begin > index
	const_iterator r = std::upper_bound( list.begin(), list.end(), range );

	bool isCovered = false;
	if ( r != list.end() ) {
		// check succeeding range and check
		isCovered = r->from <= range.from && range.to <= r->to;
	}
	if ( !isCovered && r != list.begin() ) {
		// go to preceding range and check
		--r;
		isCovered = r->from <= range.from && range.to <= r->to;
	}

	return isCovered;
}

//////////////////////////////////////////////////////////////////////

bool
IndexRangeList::
overlaps( const IndexRange& range ) const
{
#if IN_DEBUG_MODE
	if (!range.isAscending())  {
		throw std::runtime_error("IndexRangeList::overlaps("+toString(range)+") range is not ascending");
	}
#endif

	// quick check
	if (list.empty()) {
		return false;
	}

	// find first range that with begin > range.from
	const_iterator r = std::upper_bound( list.begin(), list.end(), range );

	bool isNotOverlapping = true;

	// check if not overlapping with successor
	if (isNotOverlapping && r != list.end()) {
		// check for overlap
		isNotOverlapping = range.to < r->from;
	}

	// check if not overlapping with predecessor
	if (isNotOverlapping && r != list.begin()) {
		// go to predecessor
		--r;
		// check for overlap
		isNotOverlapping = r->to < range.from;
	}

	return !isNotOverlapping;

}

//////////////////////////////////////////////////////////////////////

void
IndexRangeList::
push_back( const IndexRange& range )
{
#if IN_DEBUG_MODE
	if (!range.isAscending())  {
		throw std::runtime_error("IndexRangeList::push_back("+toString(range)+") range is not ascending");
	}
	if (!list.empty() && list.rbegin()->from >= range.from) {
		throw std::runtime_error("IndexRangeList::push_back("+toString(range)+") violates order given last range = "+toString(*(list.rbegin())));
	}
#endif
	if (!list.empty() && list.rbegin()->to >= range.from) {
		NOTIMPLEMENTED("IndexRangeList::push_back() : overlapping ranges are currently not supported")
	}
	// sorting should be OK (in debug mode.. ;) )
	list.push_back( range );
}

//////////////////////////////////////////////////////////////////////

IndexRangeList::iterator
IndexRangeList::
insert( const IndexRange& range )
{
#if IN_DEBUG_MODE
	if (!range.isAscending())  {
		throw std::runtime_error("IndexRangeList::insert("+toString(range)+") range is not ascending");
	}
#endif
	// add first member to list
	if (list.empty()) {
		list.push_back( range );
		return begin();
	} else
	// insert accordingly
	{
		// find first range that with begin > i
		List::iterator r = std::upper_bound( list.begin(), list.end(), range );
		if (r != list.end()) {
			// check for overlap
			if (range.to >= r->from) {
				NOTIMPLEMENTED("IndexRangeList::insert() : overlapping ranges are currently not supported")
			}
		}
		// check if already existing (predecessor)
		if (r != list.begin()){
			--r;
			// check if already present
			if (*r == range) {
				// return iterator to already present element
				return r;
			}
			// check for overlap
			if (r->to >= range.from) {
				NOTIMPLEMENTED("IndexRangeList::insert() : overlapping ranges are currently not supported")
			}
			++r;
		}
		// insert accordingly preserving sorting
		return list.insert( r, range );
	}
}

//////////////////////////////////////////////////////////////////////

IndexRangeList::iterator IndexRangeList::erase( IndexRangeList::iterator i ) { return list.erase( i ); }

//////////////////////////////////////////////////////////////////////

IndexRangeList::iterator IndexRangeList::begin() { return list.begin(); }

//////////////////////////////////////////////////////////////////////

IndexRangeList::iterator IndexRangeList::end() { return list.end(); }

//////////////////////////////////////////////////////////////////////

IndexRangeList::const_iterator IndexRangeList::begin() const { return list.begin(); }

//////////////////////////////////////////////////////////////////////

IndexRangeList::const_iterator IndexRangeList::end() const { return list.end(); }

//////////////////////////////////////////////////////////////////////

IndexRangeList::reverse_iterator IndexRangeList::rbegin() { return list.rbegin(); }

//////////////////////////////////////////////////////////////////////

IndexRangeList::reverse_iterator IndexRangeList::rend() { return list.rend(); }

//////////////////////////////////////////////////////////////////////

IndexRangeList::const_reverse_iterator IndexRangeList::rbegin() const { return list.rbegin(); }

//////////////////////////////////////////////////////////////////////

IndexRangeList::const_reverse_iterator IndexRangeList::rend() const { return list.rend(); }

//////////////////////////////////////////////////////////////////////

bool IndexRangeList::empty() const { return list.empty(); }

//////////////////////////////////////////////////////////////////////

size_t IndexRangeList::size() const { return list.size(); }

//////////////////////////////////////////////////////////////////////

void IndexRangeList::clear() { return list.clear(); }

//////////////////////////////////////////////////////////////////////

IndexRangeList
IndexRangeList::
shift( const int indexShift, const size_t indexMax ) const
{
	IndexRangeList l;
	IndexRange r2;
	for (IndexRangeList::const_iterator r=begin(); r!=end(); r++) {
		// skip ranges leaving the valid interval
		if ( (indexShift+((int)r->to)) < 0 || (int)indexMax < (indexShift+((int)r->from))) {
			continue;
		}
		// get shifted range boundaries
		r2.from = (size_t)std::max(0,indexShift+((int)r->from));
		r2.to = (size_t)std::min((int)indexMax,indexShift+((int)r->to));
		// store shifted and cut range
		l.insert(r2);
	}

	// final updated range list
	return l;
}

//////////////////////////////////////////////////////////////////////

IndexRangeList &
IndexRangeList::
reverse( const size_t seqLength )
{
	// reverse each entry
	size_t tmpFrom;
	for (IndexRangeList::iterator r=begin(); r!=end(); r++) {
#if IN_DEBUG_MODE
		// check if reversal is possible
		if (r->from >= seqLength || r->to >= seqLength) throw std::runtime_error("IndexRangeList::reverse("+toString(seqLength)+") = range "+toString(*r)+" exceeds seqLength");
#endif
		// reverse boundaries
		tmpFrom = r->from;
		r->from = seqLength -1 - r->to;
		r->to = seqLength -1 - tmpFrom;
	}
	// reverse order of list entries
	list.reverse();
	// return access to altered element
	return *this;
}

//////////////////////////////////////////////////////////////////////

IndexRangeList
IndexRangeList::
reverse( const size_t seqLength ) const
{
	// create copy
	IndexRangeList tmp(*this);
	// reverse and return copy
	return tmp.reverse( seqLength );
}

//////////////////////////////////////////////////////////////////////

std::ostream& operator<<(std::ostream& out, const IndexRangeList& l)
{
	// output according to regex
	for (IndexRangeList::const_iterator i=l.begin(); i!=l.end(); i++)
		out <<(i==l.begin()?"":",") <<*i;
	return out;
}

//////////////////////////////////////////////////////////////////////

void
IndexRangeList::
fromString( const std::string & stringEncoding )
{
	// clear current data
	this->clear();
	// check if something to be parsed
	if (!stringEncoding.empty()) {
		// check if parsable
		if( ! boost::regex_match(stringEncoding, IndexRangeList::regex, boost::match_perl) ) {
			throw std::runtime_error("IndexRangeList::fromString("+stringEncoding+") uses no valid index range string encoding matching '"+regex.str()+"'");
		}
		// find split position
		size_t startPos = 0, splitPos = std::string::npos;
		while (startPos != splitPos) {
			splitPos = stringEncoding.find(',',startPos);
			// insert interval
			this->insert( IndexRange(stringEncoding.substr(startPos,splitPos-(splitPos==std::string::npos?0:startPos))));
			// update start of next interval encoding to parse
			startPos = splitPos + (splitPos != std::string::npos ? 1 : 0);
		}
	}
}

//////////////////////////////////////////////////////////////////////

