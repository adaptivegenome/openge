// ***************************************************************************
// Sort.h (c) 2009 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 10 October 2011 (DB)
// ---------------------------------------------------------------------------
// Provides sorting functionality.
// ***************************************************************************

#ifndef ALGORITHMS_SORT_H
#define ALGORITHMS_SORT_H

#include "api/api_global.h"
#include "api/BamAlignment.h"
#include <cassert>
#include <algorithm>
#include <functional>
#include <string>
#include <vector>

namespace BamTools {
namespace Algorithms {

/*! \struct BamTools::Algorithms::Sort
    \brief Provides classes & methods related to sorting BamAlignments
*/
struct API_EXPORT Sort {

    //! Provides explicit values for specifying desired sort ordering
    enum Order { AscendingOrder = 0
               , DescendingOrder
               };

    /*! \fn template<typename ElemType> static inline bool sort_helper(const Sort::Order& order, const ElemType& lhs, const ElemType& rhs)
        \internal

        Determines necessary STL function object depending on requested Sort::Order
    */
    template<typename ElemType>
    static inline bool sort_helper(const Sort::Order& order, const ElemType& lhs, const ElemType& rhs) {
        switch ( order ) {
            case ( Sort::AscendingOrder  ) : { std::less<ElemType> comp;    return comp(lhs, rhs); }
            case ( Sort::DescendingOrder ) : { std::greater<ElemType> comp; return comp(lhs, rhs); }
            default : BT_ASSERT_UNREACHABLE;
        }
        return false; // <-- unreachable
    }

    //! Base class for our sorting function objects
    typedef std::less<BamAlignment> AlignmentSortBase;
    typedef std::less<BamAlignment *> AlignmentPtrSortBase;

    /*! \struct BamTools::Algorithms::Sort::ByName
        \brief Function object for comparing alignments by name

        Default sort order is Sort::AscendingOrder.

        \code
            std::vector<BamAlignment> a;

            // sort by name, in ascending order (the following two lines are equivalent):
            std::sort( a.begin(), a.end(), Sort::ByName() );
            std::sort( a.begin(), a.end(), Sort::ByName(Sort::AscendingOrder) );

            // OR sort in descending order
            std::sort( a.begin(), a.end(), Sort::ByName(Sort::DescendingOrder) );
        \endcode
    */
    struct ByName : public AlignmentSortBase, public AlignmentPtrSortBase {

        // ctor
        ByName(const Sort::Order& order = Sort::AscendingOrder)
            : m_order(order)
        { }

        // comparison function
        bool operator()(const BamTools::BamAlignment& lhs, const BamTools::BamAlignment& rhs) const {
            return sort_helper(m_order, lhs.getName(), rhs.getName());
        }
        bool operator()(const BamTools::BamAlignment * lhs, const BamTools::BamAlignment * rhs) const {
            return sort_helper(m_order, lhs->getName(), rhs->getName());
        }

        // used by BamMultiReader internals
        static inline bool UsesCharData(void) { return true; }

        // data members
        private:
            const Sort::Order& m_order;
    };

    /*! \struct BamTools::Algorithms::Sort::ByPosition
        \brief Function object for comparing alignments by position

        Default sort order is Sort::AscendingOrder.

        \code
            std::vector<BamAlignment> a;

            // sort by position, in ascending order (the following two lines are equivalent):
            std::sort( a.begin(), a.end(), Sort::ByPosition() );
            std::sort( a.begin(), a.end(), Sort::ByPosition(Sort::AscendingOrder) );

            // OR sort in descending order
            std::sort( a.begin(), a.end(), Sort::ByPosition(Sort::DescendingOrder) );
        \endcode
    */
    struct ByPosition : public AlignmentSortBase, public AlignmentPtrSortBase {

        // ctor
        ByPosition(const Sort::Order& order = Sort::AscendingOrder)
            : m_order(order)
        { }

        // comparison function
        bool operator()(const BamTools::BamAlignment& lhs, const BamTools::BamAlignment& rhs) const {
            // force unmapped aligmnents to end

            if ( lhs.getRefID() == -1 ) return false;
            if ( rhs.getRefID() == -1 ) return true;
            
            if(lhs.getRefID() != rhs.getRefID())
                return sort_helper(m_order, lhs.getRefID(), rhs.getRefID());
            if ( lhs.getPosition() != rhs.getPosition() )
                return sort_helper(m_order, lhs.getPosition(), rhs.getPosition());
            if ( lhs.IsReverseStrand() != rhs.IsReverseStrand() )
                return lhs.IsReverseStrand() ? false: true;
            if ( lhs.getName() != rhs.getName() )
                return sort_helper<std::string>(m_order, lhs.getName(), rhs.getName());
            if ( lhs.getAlignmentFlag() != rhs.getAlignmentFlag() )
                return sort_helper(m_order, lhs.getAlignmentFlag(), rhs.getAlignmentFlag());
            return sort_helper(m_order, &lhs, &rhs);
        }
        bool operator()(const BamTools::BamAlignment * lhs, const BamTools::BamAlignment * rhs) const {
            return operator()(*lhs, *rhs);
        }


        // used by BamMultiReader internals
        static inline bool UsesCharData(void) { return false; }

        // data members
        private:
            Sort::Order m_order;
    };

    /*! \struct BamTools::Algorithms::Sort::ByTag
        \brief Function object for comparing alignments by tag value

        Default sort order is Sort::AscendingOrder.

        \code
            std::vector<BamAlignment> a;

            // sort by edit distance, in ascending order (the following two lines are equivalent):
            std::sort( a.begin(), a.end(), Sort::ByTag<int>("NM") );
            std::sort( a.begin(), a.end(), Sort::ByTag<int>("NM", Sort::AscendingOrder) );

            // OR sort in descending order
            std::sort( a.begin(), a.end(), Sort::ByTag<int>("NM", Sort::DescendingOrder) );
        \endcode
    */
    template<typename T>
    struct ByTag : public AlignmentSortBase {

        // ctor
        ByTag(const std::string& tag,
              const Sort::Order& order = Sort::AscendingOrder)
            : m_tag(tag)
            , m_order(order)
        { }

        // comparison function
        bool operator()(const BamTools::BamAlignment& lhs, const BamTools::BamAlignment& rhs) {

            // force alignments without tag to end
            T lhsTagValue;
            T rhsTagValue;
            if ( !lhs.GetTag(m_tag, lhsTagValue) ) return false;
            if ( !rhs.GetTag(m_tag, rhsTagValue) ) return true;

            // otherwise compare on tag values
            return sort_helper(m_order, lhsTagValue, rhsTagValue);
        }

        // used by BamMultiReader internals
        static inline bool UsesCharData(void) { return true; }

        // data members
        private:
            std::string m_tag;
            Sort::Order m_order;
    };

    /*! \struct BamTools::Algorithms::Sort::Unsorted
        \brief Placeholder function object

        This function object exists purely to allow for dropping a "do not care" ordering
        into methods, containers, etc that are designed to work with the other sorting objects.

        \code
            std::set<BamAlignment, Sort::ByName>;   // STL set, ordered on alignment name
            std::set<BamAlignment, Sort::Unsorted>; // STL set, unsorted (but probably insertion order)
        \endcode
    */
    struct Unsorted : public AlignmentSortBase {

        // comparison function
        inline bool operator()(const BamTools::BamAlignment&, const BamTools::BamAlignment&) {
            return false;   // returning false tends to retain insertion order
        }
    };
};

} // namespace Algorithms
} // namespace BamTools

#endif // ALGORITHMS_SORT_H
