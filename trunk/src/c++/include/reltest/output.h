#ifndef RELTEST_OUTPUT_H
#define RELTEST_OUTPUT_H

//==========================================================================
//  File:       reltest_output.h
//
//  Author:     Qing Sun & Yeunjoo Song
//
//  History:    Version 1.0
//                      1.1 Updated to print new nonparametric result &
//                          default SAGE title.                - yjs 09/2001
//                      1.2 histogram & summary file merged.   - yjs 09/2001
//                      1.3 new output detailed added.         - yjs 10/2001
//                      2.0 updated to new libraries.            yjs Jul. 03
//
//  Notes:      This header file defines the interfaces to output files
//              related to RELTEST program. All output files should be
//              handled here. There are currently 4 output files from 
//              RELTEST program:  
//                         a.  histogram of Yj : sib_stats.histogram
//                         b.  histogram of Yjp: parent_child_stats.histogram
//                         c.     summary file : summary.out
//                         d.   auxiliary file : nucFam.inf
//
//              New: a, b, c files types above got merged into .sum file.
//                   nucFam.inf got renamed as .fam file.
//                   a new detailed file .det file added.        yjs Oct. 03
//
//  Copyright (c) 1999 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "reltest/analysis.h"

namespace SAGE
{

namespace RELTEST
{

class reltest_output
{
  public:

    reltest_output(const reltest_analysis* p,
                   ostream& out, ostream& fam, ostream& det, cerrorstream& err = sage_cerr);
   ~reltest_output()
    {}

     //If Yj, generate histogram for Yj.
     //else generate histogram for Yj*.
    size_t print_histogram(ostream& out, bool Yj, bool all=false);

    bool print_summary_file(ostream& o);     //Generate the summary file.
    bool print_detailed_file(ostream& o);    //Generate the summary file.
    void print_header(ostream& o, bool d);   //Generate the header of the summary & detailed files.
    void write_sib_info_file(ostream& o);    //Generate auxiliary file, "reltest.fam"

  private:

    //functions
    bool   init();                            //Set up output files.
    void   get_format_params();               //Get "max_ped_name" & "max_ind_name"
    string pair_name(putative_type) const;  

    //output files
    void print_subtitle(ostream& o);
    void print_body    (ostream& o, const vector<putative_pair>&);
    void print_trailer (ostream& o);

    void print_subtitle_det(ostream& o);
    void print_body_det    (ostream& o, const vector<putative_pair>&);

    void graph_hist    (ostream& o, int count);

    //data members
    //
    const reltest_analysis*  my_ra;
    
    putative_type      my_current_pairtype;

    bool               my_file_ready;
    size_t             max_ped_name;
    size_t             max_ind_name;

    cerrorstream       errors;
};

} // end of namespace RELTEST

} // end of namespace SAGE

#endif
