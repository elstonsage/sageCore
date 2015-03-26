#ifndef LODPAL_OUTPUT_H
#define LODPAL_OUTPUT_H

//****************************************************************************
//* File:      lodpal_out.h                                                  *
//*                                                                          *
//* Author:    Kevin Jacobs & Yeunjoo Song                                   *
//*                                                                          *
//* History:   Version 0.0 Initial implementation                kbj         *
//*                    1.0 diagnostic outfile added.             yjs Jun. 01 *
//*                                                                          *
//* Notes:     This header file defines ARPTest output classes.              *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "lodpal/arp_base_analysis.h"

namespace SAGE   {
namespace LODPAL {

//////////////////////////////////////////////////////////////////////////

class lodpal_test_viewer
{
  public:

    lodpal_test_viewer(ostream& output)
      : my_output(output), my_max_marker_name(10),
        my_max_ped_name(8), my_max_ind_name(5)
    { 
      set_wide_output(false);
      set_pval_scientific_notation(false);
    }

    lodpal_test_viewer(ostream& output, bool wide_out, bool pval_sci)
      : my_output(output), my_max_marker_name(10),
        my_max_ped_name(8), my_max_ind_name(5)
    { 
      set_wide_output(wide_out);
      set_pval_scientific_notation(pval_sci);
    }

    virtual ~lodpal_test_viewer() {}

    virtual void print_header(const ARP_base_analysis&)  = 0;
    virtual void print_results(const ARP_base_analysis&) = 0;
    virtual void print_footer(const ARP_base_analysis&)  = 0;

    virtual void print_skip_line(const ARP_base_analysis&);

    void       print_title_heading(const ARP_base_analysis& test, bool split=true);
    void       print_model_summary(const ARP_base_analysis&);

    ostream&   output_stream()        { return my_output; }

    void       set_wide_output(bool w)    { my_wide_output   = w; }
    bool       wide_output()    const     { return my_wide_output; }

    void       set_pval_scientific_notation(bool w) { my_pval_scientific_notation = w; }
    bool       get_pval_scientific_notation() const { return my_pval_scientific_notation; }
    
  protected:

    void       set_format_params(const ARP_base_analysis&);

    ostream&                  my_output;

    size_t                    my_max_marker_name;
    size_t                    my_max_ped_name;
    size_t                    my_max_ind_name;

    bool                      my_wide_output;
    bool                      my_pval_scientific_notation;
};

class lodpal_test_textfile : public lodpal_test_viewer
{
  public:

    lodpal_test_textfile(ostream& output)
      : lodpal_test_viewer(output)
    { }

    lodpal_test_textfile(ostream& output, bool wide_out, bool pval_sci)
      : lodpal_test_viewer(output, wide_out, pval_sci)
    { }

    virtual void print_header(const ARP_base_analysis&);
    virtual void print_results(const ARP_base_analysis&);
    virtual void print_footer(const ARP_base_analysis&);

  protected:

    void    print_table_heading(const ARP_base_analysis&);

    void    print_double_line(const ARP_base_analysis&);
    void    print_single_line(const ARP_base_analysis&);

    void    print_table_heading_single_line(const ARP_base_analysis& test);
    void    print_param_names(const ARP_base_analysis&, bool lambda);

    void    print_param_double_line(const ARP_base_analysis&);
    void    print_param_single_line(const ARP_base_analysis&);

    bool    is_valid_for_emp_p_value(const ARP_base_analysis&) const;
};

class lodpal_test_xfile : public lodpal_test_viewer
{
  public:

    lodpal_test_xfile(ostream& output);
    lodpal_test_xfile(ostream& output, bool wide_out, bool pval_sci);

    virtual void print_header(const ARP_base_analysis&);
    virtual void print_results(const ARP_base_analysis&);
    virtual void print_footer(const ARP_base_analysis&);

  protected:

    struct sub_type_info
    {
      string         marker_name;
      double         marker_distance;
      double         lod_score;
      size_t         fsib_pair_count;
      size_t         hsib_pair_count;
      size_t         other_pair_count;
      size_t         pair_count;

      vector<double> param_estimates;
      vector<double> first_derivative;

      size_t         function_evaluations;
      size_t         last_error;
    };
    
    void    print_table_heading(const ARP_base_analysis&, bool full = true);

    void    print_double_line(const ARP_base_analysis&);
    void    print_single_line(const ARP_base_analysis&, bool full = true);

    void    print_mm_result(const ARP_base_analysis&);
    void    print_mf_result(const ARP_base_analysis&);
    void    print_ff_result(const ARP_base_analysis&);
    
    vector<sub_type_info> my_mm_pair_result;
    vector<sub_type_info> my_mf_pair_result;
    vector<sub_type_info> my_ff_pair_result;
};

class lodpal_test_csvfile : public lodpal_test_viewer
{
  public:

    lodpal_test_csvfile(ostream& output)
      : lodpal_test_viewer(output) 
    { printed = false; }

    lodpal_test_csvfile(ostream& output, bool wide_out, bool pval_sci)
      : lodpal_test_viewer(output, wide_out, pval_sci) 
    { printed = false; }

    virtual void print_header(const ARP_base_analysis&);
    virtual void print_results(const ARP_base_analysis&);
    virtual void print_footer(const ARP_base_analysis&);

  private:

    bool printed;
};

class lodpal_test_diagfile : public lodpal_test_viewer
{
  public:

    lodpal_test_diagfile(ostream& output)
      : lodpal_test_viewer(output)
    { }

    lodpal_test_diagfile(ostream& output, bool wide_out, bool pval_sci)
      : lodpal_test_viewer(output, wide_out, pval_sci)
    { }

    virtual void print_header(const ARP_base_analysis&);
    virtual void print_results(const ARP_base_analysis&);
    virtual void print_footer(const ARP_base_analysis&);

  protected:

    void    print_final_result_summary(const ARP_base_analysis&);
    void    print_param_estimates(const ARP_base_analysis&);
    void    print_var_cov_matrix(const ARP_base_analysis&);
    void    print_histogram(const ARP_base_analysis&);
    void    print_individual_lod_table(const ARP_base_analysis&);
    void    print_type_regend(const ARP_base_analysis&);

    void    print_double_line(const ARP_base_analysis&);
    void    print_single_line(const ARP_base_analysis&);
};

//////////////////////////////////////////////////////////////////////////

} // end of namespace LODPAL
} // end of namespace SAGE

#endif
