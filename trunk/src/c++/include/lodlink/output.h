#ifndef LODLINK_OUTPUT_H
#define LODLINK_OUTPUT_H
//============================================================================
// File:      output.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   2/27/3 - created.                                   djb
//                                                                          
// Notes:     Common data, classes and functions for output routines.
//                                                                          
// Copyright (c) 2003 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include "maxfun/sub_model.h"

using std::ostream;
using std::string;
using std::vector;

namespace SAGE
{

namespace LODLINK
{

class col_hdr;
class header;

extern const int  PRC1;   // theta, alpha
extern const int  PRC2;   // lod scores, chi square stat
extern const int  PRC3;   // p-values
extern const int  PRC4;   // variance-covariance matrix elements

extern const char U_CHR;

// - Default sizes.
//
extern const int  SPACE_SZ;
extern const int  LOCUS_SZ;
extern const int  RECOM_SZ;
extern const int  LOD_SZ;
extern const int  CHI_SZ;
extern const int  P_VALUE_SZ;
extern const int  P_VALUE_UB_SZ;
extern const int  PROP_SZ;
extern const int  NATURAL_LOG_SZ;
extern const int  POST_SZ;
extern const int  OPTION_SZ;
extern const int  VAR_COV_SZ;
extern const int  GENOTYPE_SZ;

// - Labels (Trailing underscore indicates that name is a abbreviated
//   version of the content).
//
extern const string  LOCUS;
extern const string  BLANK;
extern const string  RECOM;
extern const string  MALE_RC;
extern const string  FEMALE_RC;
extern const string  FRACT_IN;
extern const string  ESTIMATES_IN;
extern const string  SMALL_RANGE;
extern const string  LARGE_RANGE;
extern const string  LOD;
extern const string  SCORE;
extern const string  CHI;
extern const string  SQUARE;
extern const string  STAT;
extern const string  P_VALUE;
extern const string  UPPER;
extern const string  BOUND;
extern const string  EST;
extern const string  LINK_;
extern const string  PROP;
extern const string  POSTERIOR;
extern const string  PROBABILITY;
extern const string  NATURAL_LOG;
extern const string  MAXIMUM;
extern const string  LIKELIHOOD;
extern const string  TYPE_AA;
extern const string  TYPE_AB;
extern const string  TYPE_BB;

extern const string  LOD_SCORES;
extern const string  LOD_SCORES_FAM_;
extern const string  MLE;
extern const string  LOD_SCORE_LINKAGE_;
extern const string  LOD_SCORE_LINKAGE_TEST;
extern const string  USING_RECOM_SMALL_;
extern const string  USING_RECOM_LARGE_;
extern const string  SEX_SPECIFIC_;
extern const string  SEX_AVERAGED_;
extern const string  CLEVES_ELSTON_;
extern const string  SMITHS_;
extern const string  FARAWAYS_;
extern const string  MORTONS_;
extern const string  M_F_;
extern const string  POSTERIORS_FAM_;
extern const string  EVIDENCE_;
extern const string  MATRICES_;
extern const string  PARAM_ORDER_A_;
extern const string  PARAM_ORDER_M_F_;
extern const string  PARAM_ORDER_A_A_;
extern const string  PARAM_ORDER_M_F_A_;
extern const string  GENOTYPES;

// - Column headers (Trailing underscore indicates that name differs
//   content of first line).
//
extern const col_hdr  _LOCUS;
extern const col_hdr  _TWO_LINE_LOCUS;
extern const col_hdr  _THREE_LINE_LOCUS;
extern const col_hdr  _RECOM_SMALL_;
extern const col_hdr  _RECOM_LARGE_;
extern const col_hdr  _M_RECOM_SMALL_;
extern const col_hdr  _M_RECOM_LARGE_;
extern const col_hdr  _F_RECOM_SMALL_;
extern const col_hdr  _F_RECOM_LARGE_;
extern const col_hdr  _LOD;
extern const col_hdr  _CHI;
extern const col_hdr  _P_VALUE;
extern const col_hdr  _P_VALUE_UB_;
extern const col_hdr  _EST_LINKAGE_PROP_;
extern const col_hdr  _POSTERIOR;
extern const col_hdr  _LIKELIHOOD_;
extern const col_hdr  _TYPE_AA;
extern const col_hdr  _TYPE_AB;
extern const col_hdr  _TYPE_BB;

// - Meta headers (Same conventions as for column headers.)
//
extern const col_hdr  _MLE_AVE_;
extern const col_hdr  _MLE_SEX_SPEC_;
extern const col_hdr  _LOD_SCORE_;
extern const col_hdr  _CLEVES_ELSTON_;
extern const col_hdr  _SMITHS_;
extern const col_hdr  _FARAWAYS_;
extern const col_hdr  _MORTONS_SS_;
extern const col_hdr  _MORTONS_NON_SS_;
extern const col_hdr  _LODS_AVE_;
extern const col_hdr  _LODS_AVE_FAM_;
extern const col_hdr  _LODS_SEX_SPEC_;
extern const col_hdr  _LODS_SEX_SPEC_FAM_;
extern const col_hdr  _CLEVES_ELSTON_VC_;
extern const col_hdr  _SMITHS_NON_SS_VC_;
extern const col_hdr  _SMITHS_SS_VC_;
extern const col_hdr  _FARAWAYS_NON_SS_VC_;
extern const col_hdr  _FARAWAYS_SS_VC_;
extern const col_hdr  _NON_SS_LOD_SCORE_VC_;
extern const col_hdr  _SS_LOD_SCORE_VC_;

void  write_double(ostream& out, double d, bool p_value = false);

//----------------------------------------------------------------------------
//  Class:    col_hdr
//                                                                          
//  Purpose:  represent a column header of three lines.
//                                                                          
//----------------------------------------------------------------------------
//
class col_hdr
{
  public:
    struct labels
    {
      labels(const string& one, const string& two, const string& three); 
      
      vector<const string*>  data;
    };

  col_hdr(size_t label_width, size_t spacer_width, const labels& ls);
  
  size_t  lw() const;
  size_t  sw() const;
  size_t  tw() const;
  const string&  operator[](size_t i) const;

  private:
    size_t  my_label_width;
    size_t  my_spacer_width;
    const labels  my_labels;
};


//----------------------------------------------------------------------------
//  Class:    header
//                                                                          
//  Purpose:  represents an entire header consisting of three lines and
//            and optional underline.
//                                                                          
//----------------------------------------------------------------------------
//
class header
{
  public:
    friend ostream&  operator<<(ostream& out, const header& h);
  
    header(size_t offset = 0, size_t underline = '0');
    
    size_t  offset() const;
    size_t  col_w(size_t i) const;
    size_t  spc_w(size_t i) const;
    bool    underline() const;
    const col_hdr&  operator[](size_t index) const;
    
    void  set_offset(size_t offset);
    void  set_underline(char underline);
    void  add_col(const col_hdr& col);

  private:
    size_t  my_offset;
    char    my_underline;             // '0' means no underline.
    vector<const col_hdr*>  my_cols;
};

#include "lodlink/output.ipp"
}
}

#endif

