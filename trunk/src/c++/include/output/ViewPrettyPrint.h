#ifndef OUTPUT_VIEW_PRETTYPRINT_H
#define OUTPUT_VIEW_PRETTYPRINT_H

#include "output/Graph.h"
#include "output/Misc.h"
#include "output/RenderingRules.h"
#include "output/Table.h"
#include "util/OutlineCntr.h"

namespace SAGE   {
namespace OUTPUT {

class ViewPrettyPrint
{
  public:

    ///
    /// Renders to a file.
    template<typename T>
    static void renderToFile(const T & element, const std::string & filename)
    {
      std::ofstream o;
      o.open(filename.c_str());

      o << render(element) << std::flush;
    }

    ///
    /// PrettyPrint rendering.
    template<typename T>
    static std::string render(const T & element)
    {
      return internal_render(element, UTIL::OutlineCntr());
    }

  private:

  /// @name Pretty-print functions by type
  //@{

    ///
    /// No pretty-print function for given type.
    ///
    template<typename T>
    static std::string internal_render(const T & element, const UTIL::OutlineCntr & _cntr)
    {
      BOOST_STATIC_ASSERT((sizeof(T) == 0)); // No pretty-print function defined for this type!
      throw std::exception(); // Never happens--only here to make compiler shut up about missing return value.
    }
    
    ///
    /// Section
    ///
    static std::string internal_render(const Section & section, const UTIL::OutlineCntr & _cntr);
    
    ///
    /// Section
    ///
    static std::string internal_render(const List & list_, const UTIL::OutlineCntr & _cntr, bool is_top_level_list = true, HasBulletType::BulletType e = HasBulletType::NUMBERED);
    
    ///
    /// BasicValue<V>
    ///
    template<typename V>
    static std::string internal_render(const BasicValue<V> & val, const UTIL::OutlineCntr & _cntr) 
    { 
      return default_rules.render(val);
    }
    
    ///
    /// NamedValue<V>
    ///
    template<typename V>
    static std::string internal_render(const NamedValue<V> & val, const UTIL::OutlineCntr & _cntr)
    {
      std::ostringstream s; 
      
      s << (val.getTitle() == "" ? "" : val.getTitle() + ": ") << default_rules.render(val) << std::endl; 
      
      return s.str();
    }


    ///
    /// Table
    ///
    static std::string internal_render(const Table & table, const UTIL::OutlineCntr & _cntr);

    class TableRenderer
    {
      public:
      
        explicit TableRenderer(const Table & table);
         
        TableRenderer(const TableRenderer & other);
      
        std::string render(const UTIL::OutlineCntr & _cntr) const;
         
        /// Text sequence for delimiting between table cells (ie: " | ").
        static std::string delimiter;

        /// Text sequence for delimiting between table cells in a horizontal line (ie: "-+-")
        static std::string horizontal_delimiter;

        /// Character used in generating horizontal divider lines (ie: '-').
        static char horizontal_char;

      private:
      
      /// @name Precalculating output stuff
      //@{
  
        ///
        /// Calculates the maximum number of digits to left of the decimal place for each
        /// column (by looking at all the doubles in each column).
        void calculateLeftLengths() const;

        ///
        /// Calculates the rendered text for each cell. This assumes calculateLeftLengths()
        /// has already been invoked.
        void calculateRenderedCells() const;

        ///
        /// Figures out how wide every column has to be. This assumes calculateRenderedCells()
        /// has already been invoked.
        void calculateMaxWidths() const;
    
        ///
        /// Renders the given cell.
        std::string renderCell(const TableRow & row, size_t column_idx, const RenderingRules & rules) const;

      //@}
  
      /// @name Output generation
      //@{

        ///
        /// Returns the header (or blank if there's no title).
        std::string generateHeader() const;

        ///
        /// Returns a blank line (for output formatting).
        std::string generateBlankLine() const;

        ///
        /// Returns a horizontal bar (for output formatting).
        std::string generateHBar() const;

        ///
        /// Returns the column headers.
        ///
        /// NOTE: Assumes calculateMaxWidths() has been invoked.
        std::string generateColumns() const;

        ///
        /// Returns all the rwos (for output formatting).
        ///
        /// NOTE: Assumes calculateMaxWidths() has been invoked.
        std::string generateRows() const;

        ///
        /// Returns the given row (for output formatting).
        ///
        /// NOTE: Assumes calculateMaxWidths() has been invoked.
        std::string generateRow(size_t row_idx, const TableRow & row, const RenderingRules & table_rules) const;

        ///
        /// Returns the table footer (comments included).
        std::string generateFooter() const;
    
      //@}

        TableRenderer & operator= (const TableRenderer & other); // Disallowed
        
        const Table & my_table;

        /// For every i'th column, indicates the maximum 'left-length' of the Double's in that column.
        mutable std::vector<size_t> my_left_lengths;

        mutable std::vector<std::vector<std::string> > my_rendered_cells;
        
        /// For storing column widths
        mutable std::vector<size_t> my_widths;
            
    };

  //@}
    
};

} // End namespace OUTPUT
} // End namespace SAGE


namespace std {

/// \brief For directing an element's contents to an output stream.
/// 
/// This is equivalent to invoking toPrettyPrint(0) on the object.
template<typename SELF_TYPE,
         typename T0,
         typename T1,
         typename T2,
         typename T3, 
         typename T4,
         typename T5,
         typename T6,
         typename T7,
         typename T8,
         typename T9>

inline std::ostream & operator<< (std::ostream & o, const SAGE::OUTPUT::Element<SELF_TYPE, T0, T1, T2, T3, T4, T5, T6, T7, T8, T9> & e)
{
  return (o << SAGE::OUTPUT::ViewPrettyPrint::render((const SELF_TYPE &)e));
}

} // End namespace std

#endif
