#ifndef  RELTEST_ANALYSIS_H
#define  RELTEST_ANALYSIS_H

//==========================================================================
//  File:       analysis.h
//
//  Author:     Qing Sun & Yeunjoo Song
//
//  History:    Version 1.0
//                      1.01  Updated for marker bug           gcw 00-11-03
//                      1.1   Added nonparametric estimation   yjs 01-09-13 
//                      1.11  Took out output from analysis    yjs 01-09-27
//                      2.0   Updated to new libraries         yjs  Jul. 03
//
//  Notes:      This class is the main driver performing analysis.
//
//  Copyright (c) 2001 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "ibd/basic_storage_ibd.h"
#include "ibd/exact_ibd_analysis.h"
#include "reltest/parser.h"
#include "reltest/l2_procedure.h"

namespace SAGE
{

namespace RELTEST
{

class reltest_analysis
{
  public:

    reltest_analysis(reltest_data& rtfile, const reltest_parser& rp, cerrorstream& err = sage_cerr);

    ~reltest_analysis();

    bool run_analysis(ostream& sum_out, ostream& nuc_out, ostream& det_out);

    void dump_pairs(ostream& o=cout)                             const;
    void dump_current_pairtype(putative_type p, ostream& o=cout) const;

    const reltest_parser*        get_parser()        const;
    const RPED::RefMultiPedigree*      get_multipedigree() const;

    const vector<putative_pair>& get_ptt_pairs(putative_type p)  const;
    const vector<putative_pair>& get_misc_pairs(putative_type p) const;

    const putative_type&         get_current_pairtype()    const;

    double                       get_chrom_length(size_t)  const;
    double                       get_total_genome_length() const;
    double                       get_total_map_points()    const;

    double                       get_cutpoints(cutpoint_type c)          const;
    double                       get_adjusted_cutpoints(cutpoint_type c) const;

    double                       get_AMIC()           const;
    double                       get_Var_Yj()         const;
    double                       get_Var_Yjp()        const;

    const vector<solution_type>& get_solution_Yj()    const;
    const vector<solution_type>& get_solution_Yjp()   const;

    size_t                       get_picked_Yj()      const;
    size_t                       get_picked_Yjp()     const;

    double                       get_mu_Yj()          const;
    double                       get_mu_Yjp()         const;

    double                       get_Yj_mean()        const;
    double                       get_standard_error() const;

  private:

    // ANALYSIS
    //
    bool  build();         //preparations before starting analysis.
    bool  do_analysis();

    // return pair name giving a pair type
    //
    string get_pair_name(putative_type type) const;

    // PAIR-GENERATION FUNCTIONS
    //
    void  generate_pairs();

    // if any one in a pair has missing data at all markers over the genome,
    //  then this pair is uninformative pair and it is removed.
    void  remove_uninformative_pairs();

    // return the number of markers at which an individual has missing data.
    int   uninformative_markers(ind_id)  const;

    void  pair_analysis(); //do analyses for each specified putative pair type

    IBD*  get_pair_ibds(const region& r, const FPED::Subpedigree& subped, const putative_pair& p);

    bool  update_pair_stats(int pair_index, IBD* ibd, const region& r); 
    void  stats();
    void  reclassify_pairs();

    // find misclassified pairs.
    void  get_misclassified_pairs();

    // L2_error_procedure
    void  nonparametric_estimation();
    void  adjust_cutpoints(bool is_Yj);
    void  encode_params_trial_vector(vector< std::pair<double, double> >& init_theta);
    void  initialize_maxfun(Maxfun& maxfun);
    void  run_optimum_maxfun(Maxfun& maxfun);
    void  test_deviation();

    void  do_new_L2_procedure(vector<solution_type>& solution, bool is_Yj);

    bool  check_mean_validity(size_t max_pos, bool is_Yj);
    pair<double, size_t>  estimate_shift(bool is_Yj, vector<solution_type>&);

    //
    // MEMBER DATA
    //
    reltest_data*                     my_input;

    const reltest_parser*             my_parser;
    const RPED::RefMultiPedigree*     my_multipedigree;

    vector< vector<putative_pair> >   my_ptt_pairs;  //all putative pairs
    vector< vector<putative_pair> >   my_misc_pairs; //all reclassified pairs

    putative_type   my_current_pairtype;           //putative type of the pairs currently
                                                   //being analyzed.

    vector<double>  my_chrom_length;               //the lengths of all regions in genome.

    double          my_total_genome_length;        //total genome length used in analysis
    double          my_total_map_points;           //total map points used to generate IBD

    double          my_cutpoints[Call];            //cutpoints
    double          my_adjusted_cutpoints[Call];   //adjusted_cutpoints

    double          my_AMIC;                       //average marker information content
    double          my_Var_Yj;                     //variance of Yj  over all k regions
    double          my_Var_Yjp;                    //variance of Yjp over all k regions

    exact_ibd_analysis*    my_exact_ibd_analysis;

    vector<solution_type>  my_solution_Yj;
    vector<solution_type>  my_solution_Yjp;

    size_t                 my_picked_Yj;
    size_t                 my_picked_Yjp;

    double                 my_mu_Yj;
    double                 my_mu_Yjp;

    double                 my_Yj_mean;
    double                 my_standard_error;

    bool                   my_unexpected_error;

    cerrorstream           errors;
};

#include "reltest/analysis.ipp"

} // end of namespace RELTEST

} // end of namespace SAGE

#endif
