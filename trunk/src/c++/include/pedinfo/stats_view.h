#ifndef STATS_VIEW_H
#define STATS_VIEW_H
//============================================================================
// File:     stats_view.h
//                                                                          
// Author:   Dan Baechle
//                                                                          
// History:  Initial version: 10/00
//                                                                          
// Notes:    Supercedes stat_view.h
// 
//           Declares the following classes -
//              base_stats_viewer
//              base_trait_stats_viewer
//              binary_stats_viewer
//              cont_stats_viewer
//              cmpd_stats_viewer
//              general_stats_viewer
//              ped_stats_viewer
//              summary_stats_viewer
//              mp_stats_viewer
//              log_histogram
//                                                                          
// Copyright (c) 2000 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include <iostream>
#include <string>
#include <math.h>
#include "numerics/sinfo.h"
#include "numerics/histogram.h"
#include "pedinfo/stats.h"

namespace SAGE {
namespace RPED {

//----------------------------------------------------------------------------
//  Class:    base_stats_viewer
//                                                                          
//  Purpose:  Base class for all other "stats_viewer" classes in this file.
//                                                                          
//----------------------------------------------------------------------------
//
class base_stats_viewer
{
  protected:
    // Constructor/destructor.
    base_stats_viewer(std::ostream& o);
  
    static const size_t LINE_SIZE        = 78;        // Must be an even number for header
                                                      // centering to work correctly.
    static const size_t FIELD_SIZE       = 10;        // For a 'count' value.
    static const size_t FIELD_SIZE_CORR  = 7;         // For a 'correlation' value.
  
    void         outer_line()                             const;   // ====== ... =====
    void         single_inner_line()                      const;   // |----- ... ----|
    void         double_inner_line()                      const;   // |===== ... ====|
    void         major_header(std::string title)          const;
    void         header_body(std::string title)           const;
    void         header_line(std::string line)            const;
    std::string  pair_label(pg::pair_type type)           const;
    void         sib_header()                             const;
    
    // - Provides uniform formatting for floating point values.
    //    prec is number of digits after the decimal point.  Used where value >= 1.
    //    sig_dig is number of significant digits.  Used where value < 1.
    //
    //    Scientific notation used as a last resort to fit number into space allowed.
    //     Width should be >= 7 + prec to accommodate scientific notation.
    //     Width should be >= 7 + sig_dig - 1 to accomodate scientific notation.
    //
    void         display_double(double value, unsigned int width,
                                unsigned int prec = 2, unsigned int sig_dig = 2)  const; 
    
    // Data members.
    std::ostream&     my_o;
};

//----------------------------------------------------------------------------
//  Class:    base_trait_stats_viewer
//                                                                          
//  Purpose:  Display data contained in a Base_trait_stats object.
//                                                                          
//----------------------------------------------------------------------------
//
class base_trait_stats_viewer : public base_stats_viewer
{
  public:
    // Constructor/destructor.
    base_trait_stats_viewer(std::ostream& o, const Base_trait_stats& base_trait_data);
    
    void          view()                 const;
    
  private:
    static const size_t FIELD_SIZE_BIN = 15;  
    static const size_t FIELD_SIZE_PED    = 19;
    static const size_t FIELD_SIZE_LABEL  = 16;
  
    void          header()               const;
    void          data_line()            const;
    void          sib_line()             const;
    void          ped_size_line()        const;
    void          nuclear_families()     const;
  
    // Data members.
    const Base_trait_stats&      my_base_trait_stats;
};

//----------------------------------------------------------------------------
//  Class:    binary_stats_viewer
//                                                                          
//  Purpose:  Display a Binary_trait_stats object.
//                                                                          
//----------------------------------------------------------------------------
//
class binary_stats_viewer : public base_stats_viewer
{
  public:
    typedef Binary_trait_stats     bt;
  
    // Constructor/destructor.
    binary_stats_viewer(std::ostream& o, const Binary_trait_stats& binary_stats);
    
    void view()                const;
    
  private:
    static const size_t FIELD_SIZE_PED  = 10;
    
    void          ind_data()                                  const;
    std::string   build_title()                               const;
    void          ind_header()                                const;
    std::string   gender_label(bt::gender g)                  const;
    std::string   founder_status_label(bt::founder_status f)  const;
    void          ind_gender_data_line(bt::gender g)          const;
    void          ind_blank_data_line()                       const;
    void          ind_founder_data_line(bt::founder_status f) const;
    void          ind_gender_total_line()                     const;
    void          ind_founder_total_line()                    const;
    void          pair_data()                                 const;
    void          pair_header()                               const;
    void          pair_data_line(pg::pair_type t)             const;
    
    // Data members.
    const Binary_trait_stats&    my_binary_stats;
    base_trait_stats_viewer          my_base_trait_stats_viewer;
};

//----------------------------------------------------------------------------
//  Class:    cont_stats_viewer
//                                                                          
//  Purpose:  Display a Cont_trait_stats object.
//                                                                          
//----------------------------------------------------------------------------
//
class cont_stats_viewer : public base_stats_viewer
{
  public:
    typedef Cont_trait_stats     ct;
  
    // Constructor/destructor.
    cont_stats_viewer(std::ostream& o, const Cont_trait_stats& cont_stats);
    
    void view()                const;
    
  private:
    void          ind_data()                                  const;
    std::string   build_title()                               const;
    void          ind_header()                                const;
    std::string   gender_label(ct::gender g)                  const;
    std::string   founder_status_label(ct::founder_status f)  const;
    void          ind_gender_data_line(ct::gender g)          const;
    void          ind_blank_data_line()                       const;
    void          ind_founder_data_line(ct::founder_status f) const;
    void          pair_data()                                 const;
    void          pair_data_alt()                             const;
    void          pair_header()                               const;
    void          pair_header_alt()                           const;
    void          pair_data_line(pg::pair_type t)             const;
    void          pair_data_half_line(pg::pair_type t)        const;
    void          pair_data_line_alt(pg::pair_type t)         const;
    
    // Data members.
    const Cont_trait_stats&    my_cont_stats;
    base_trait_stats_viewer        my_base_trait_stats_viewer;
};

//----------------------------------------------------------------------------
//  Class:    cmpd_stats_viewer
//                                                                          
//  Purpose:  Display a Cmpd_trait_stats object.
//                                                                          
//----------------------------------------------------------------------------
//
class cmpd_stats_viewer : public base_stats_viewer
{
  public:
    typedef Cmpd_trait_stats     cp;
  
    static const size_t FIELD_SIZE_PED = 12;
  
    // Constructor/destructor.
    cmpd_stats_viewer(std::ostream& o, const Cmpd_trait_stats& cmpd_stats);
    
    void          view()                                   const;
    
  private:  
    std::string   build_title()                            const;
    void          pair_ind_data()                          const;
    void          pair_ind_header()                        const;
    void          pair_half_line(pg::pair_type t)          const;
    size_t        ind_half_line(int i)                     const;
    void          total_ind_half_line(size_t total)        const;
    
    // Data members.
    const Cmpd_trait_stats&    my_cmpd_stats;
    base_trait_stats_viewer        my_base_trait_stats_viewer;
};

//----------------------------------------------------------------------------
//  Class:    general_stats_viewer
//                                                                          
//  Purpose:  Parent class of ped_stats_viewer and summary_stats_viewer.
//                                                                          
//----------------------------------------------------------------------------
//
class general_stats_viewer : public base_stats_viewer
{
  public:
    void  multiple_mates() const;
    void  cons_pairs() const;

  protected:
    static const size_t FIELD_SIZE_PED = 12;
  
    // Constructor/destructor.
    general_stats_viewer(std::ostream& o, const General_Stats& gen_stats);
    
    void view(std::string title)                           const;
    
    void          iv_bit_line()                            const;

    // Essentially the same info as sib_line()
    //void          nuc_line()                               const;
    void          sib_line()                               const;
    void          sub_loop_header()                        const;
    void          sub_loop_line()                          const;
    void          pair_ind_data()                          const;
    void          pair_ind_header()                        const;
    void          pair_half_line(pg::pair_type t)          const;
    size_t        ind_half_line(int i)                     const;
    void          total_ind_half_line(size_t total)        const;
    
    static bool  wrap_not_needed(const ostringstream& mate_line, const string& temp);
    
    // Data members.
    const General_Stats&  my_gen_stats;
};

//----------------------------------------------------------------------------
//  Class:    ped_stats_viewer
//                                                                          
//  Purpose:  Display individual pedigree data.
//                                                                          
//----------------------------------------------------------------------------
//
class ped_stats_viewer : public general_stats_viewer
{
  public:
    // Constructor/destructor.
    ped_stats_viewer(std::ostream& o, const Ped_stats& ped_stats, std::string ped_name);

    void          view()                        const;   
    
  private:
    std::string   build_title()                 const;
    
    // Data members.
    std::string   my_ped_name;
};

//----------------------------------------------------------------------------
//  Class:    summary_stats_viewer
//                                                                          
//  Purpose:  Display collective pedigree data contained in an MP_stats object.
//                                                                          
//----------------------------------------------------------------------------
//
class summary_stats_viewer : public general_stats_viewer
{
  public:
    // Constructor/destructor.
    summary_stats_viewer(std::ostream& o, const MP_stats& multi_stats);

    void          view()                        const;   
};

//----------------------------------------------------------------------------
//  Class:    mp_stats_viewer
//                                                                          
//  Purpose:  Display an total MP_stats object.
//                                                                          
//----------------------------------------------------------------------------
//
class mp_stats_viewer : public base_stats_viewer
{
  public:
    // Constructor/destructor.
    mp_stats_viewer(std::ostream& o, const MP_stats& multi_stats);
    
    void view(bool each = false, bool suppress_general = false) const;
    
  private:
    static const size_t FIELD_SIZE_BIN = 15;
    
    void          ped_size_header()       const;
    void          ped_size_line()         const;
    void          histograms_header()     const;
    void          histograms()            const;
    std::string   build_title()           const;
  
    // Data members.
    const MP_stats&      my_mp_stats;
};

//----------------------------------------------------------------------------
//  Class:    log_histogram
//                                                                          
//  Purpose:  A simple histogram whose bin boundaries vary logarithmically.  The
//            logarithmic base can be specified.
//            Number of bins can be arbitrarily limited, in which case the upper
//            bound of the last bin is as big as is needed to accomodate the data
//            from the underlying "conventional" histogram from which it is constructed.
//                                                                          
//----------------------------------------------------------------------------
//
class log_histogram
{
  friend class  mp_stats_viewer;
  friend class  base_trait_stats_viewer;
  
  public:
    size_t                     size()                 const;
    size_t                     operator[](size_t i)   const;
    std::pair<size_t, size_t>  boundaries(size_t i)   const; 
    bool                       limited()              const;
    size_t                     offset()               const;   // Returns number of empty bins
                                                               // before first non-empty bin.
  
  private:
    // Constructor/destructor.
    log_histogram(const Histogram& histogram, unsigned int base = 2, size_t size_limit = 6);
    size_t  prelim_size(const Histogram& histogram)   const;
    
    // Data members.
    std::vector<size_t>  my_counts;
    unsigned int         my_base;        // Logarithmic base.
    size_t               my_size;        // Number of bins.
    bool                 my_limited;     // Necessary number of bins exceeds that allowed by size limit.
};

} // End namespace RPED
} // End namespace SAGE

#include "pedinfo/stats_view.ipp"

#endif
