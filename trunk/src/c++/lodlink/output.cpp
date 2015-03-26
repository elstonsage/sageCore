//============================================================================
// File:      output.cpp                      
//                                                                          
// Author:    Dan Baechle                                     
//                                                                          
// History:   2/27/3 - created.                         djb
//                                                                          
// Notes:     Implementation classes and functions common to output routines.   
//               
//                                                                          
// Copyright (c) 2003 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "lodlink/output.h"
#include "util/StringUtils.h"
#include "output/Output.h"

using namespace std;

namespace SAGE
{

namespace LODLINK
{

const int  PRC1 = 4;   // theta, alpha
const int  PRC2 = 3;   // lod scores, chi square stat
const int  PRC3 = 4;   // p-values, ln likelihoods
const int  PRC4 = 2;   // variance-covariance matrix elements

extern const char U_CHR = '-';

// - Default sizes.
//
const int  SPACE_SZ        = 3;
const int  LOCUS_SZ        = 15;
const int  RECOM_SZ        = 12;
const int  LOD_SZ          = 12;  
const int  CHI_SZ          = 12;  
const int  P_VALUE_SZ      = 12;  
const int  P_VALUE_UB_SZ   = 12;  
const int  PROP_SZ         = 7;
const int  NATURAL_LOG_SZ  = 11;
const int  POST_SZ         = 11;
const int  OPTION_SZ       = 41;
const int  VAR_COV_SZ      = 12;
const int  GENOTYPE_SZ     = 12;

// - Labels
//
const string  LOCUS        = "Locus";
const string  BLANK        = "";
const string  RECOM        = "Recom";
const string  MALE_RC      = "Male Rc";
const string  FEMALE_RC    = "Female Rc";
const string  FRACT_IN     = "Fract in";
const string  ESTIMATES_IN = "Estimates in ";
const string  SMALL_RANGE  = "[0, .5]";
const string  LARGE_RANGE  = "[0, 1]";
const string  LOD          = "Lod";
const string  SCORE        = "Score";
const string  CHI          = "Chi";
const string  SQUARE       = "Square";
const string  STAT         = "Stat";
const string  P_VALUE      = "P-Value";
const string  UPPER        = "Upper";
const string  BOUND        = "Bound";
const string  EST          = "Est";
const string  LINK_        = "Linkage";
const string  PROP         = "Prop";
const string  POSTERIOR    = "Posterior";
const string  PROBABILITY  = "Probability";
const string  NATURAL_LOG  = "Natural Log";
const string  MAXIMUM      = "Maximum";
const string  LIKELIHOOD   = "Likelihood";
extern const col_hdr  _TYPE_AA;extern const col_hdr  _TYPE_AA;extern const col_hdr  _TYPE_AA;extern const col_hdr  _TYPE_AA;extern const col_hdr  _TYPE_AA;const string  TYPE_AA      = "Type AA";
extern const col_hdr  _TYPE_AB;extern const col_hdr  _TYPE_AB;extern const col_hdr  _TYPE_AB;extern const col_hdr  _TYPE_AB;extern const col_hdr  _TYPE_AB;const string  TYPE_AB      = "Type AB";
extern const col_hdr  _TYPE_BB;extern const col_hdr  _TYPE_BB;extern const col_hdr  _TYPE_BB;extern const col_hdr  _TYPE_BB;extern const col_hdr  _TYPE_BB;const string  TYPE_BB      = "Type BB";

const string  LOD_SCORES             = "Lod Scores";
const string  LOD_SCORES_FAM_        = "Lod Scores By Constituent Pedigree";
const string  MLE                    = "MLE";
const string  LOD_SCORE_LINKAGE_     = "Lod Score and Linkage Test";
const string  LOD_SCORE_LINKAGE_TEST = "Lod Score Linkage Test";
const string  SEX_SPECIFIC_          = "Sex-Specific Recombination Fractions";
const string  SEX_AVERAGED_          = "Non-Sex-Specific Recombination Fractions";
const string  USING_RECOM_SMALL_     = "Using Recom. Fract. in [0, .5]";
const string  USING_RECOM_LARGE_     = "Using Recom. Fract. in [0, 1]";
const string  CLEVES_ELSTON_         = "Cleves-Elston Linkage Test";
const string  SMITHS_                = "Smith's Homogeneity Test";
const string  FARAWAYS_              = "Faraway's Linkage Test";
const string  MORTONS_               = "Morton's Homogeneity Test";
const string  M_F_                   = "Male Recom - 1st Row, Female Recom - 2nd Row";
const string  POSTERIORS_FAM_        = "Posterior Probabilities by Constituent Pedigree";
const string  EVIDENCE_              = "** Indicates Evidence for Linkage";
const string  MATRICES_              = "Variance-Covariance Matrices";
const string  PARAM_ORDER_A_         = "Parameter Order (Avg Recomb)";
const string  PARAM_ORDER_M_F_       = "Parameter Order (Male Recomb, Female Recomb)";
const string  PARAM_ORDER_A_A_       = "Parameter Order (Avg Recomb, Alpha)";
const string  PARAM_ORDER_M_F_A_     = "Parameter Order (Male Recomb, Female Recomb, Alpha)";
const string  GENOTYPES              = "Individual Genotype Probabilities";

typedef col_hdr::labels  lbs;

// - Column headers.
//
const col_hdr  _LOCUS(LOCUS_SZ, SPACE_SZ, lbs(LOCUS, BLANK, BLANK));
const col_hdr  _TWO_LINE_LOCUS(LOCUS_SZ, SPACE_SZ, lbs(BLANK, LOCUS, BLANK));
const col_hdr  _THREE_LINE_LOCUS(LOCUS_SZ, SPACE_SZ, lbs(BLANK, BLANK, LOCUS));
const col_hdr  _RECOM_SMALL_(RECOM_SZ, SPACE_SZ, lbs(RECOM, FRACT_IN, SMALL_RANGE));
const col_hdr  _RECOM_LARGE_(RECOM_SZ, SPACE_SZ, lbs(RECOM, FRACT_IN, LARGE_RANGE));
const col_hdr  _M_RECOM_SMALL_(RECOM_SZ, SPACE_SZ, lbs(MALE_RC, FRACT_IN, SMALL_RANGE));
const col_hdr  _M_RECOM_LARGE_(RECOM_SZ, SPACE_SZ, lbs(MALE_RC, FRACT_IN, LARGE_RANGE));
const col_hdr  _F_RECOM_SMALL_(RECOM_SZ, SPACE_SZ, lbs(FEMALE_RC, FRACT_IN, SMALL_RANGE));
const col_hdr  _F_RECOM_LARGE_(RECOM_SZ, SPACE_SZ, lbs(FEMALE_RC, FRACT_IN, LARGE_RANGE));
const col_hdr  _LOD(LOD_SZ, SPACE_SZ, lbs(LOD, SCORE, BLANK));
const col_hdr  _CHI(CHI_SZ, SPACE_SZ, lbs(CHI, SQUARE, STAT));
const col_hdr  _P_VALUE(P_VALUE_SZ, SPACE_SZ, lbs(P_VALUE, BLANK, BLANK));
const col_hdr  _P_VALUE_UB_(P_VALUE_UB_SZ, SPACE_SZ, lbs(P_VALUE, UPPER, BOUND));
const col_hdr  _EST_LINKAGE_PROP_(PROP_SZ, SPACE_SZ, lbs(EST, LINK_, PROP));
const col_hdr  _POSTERIOR(POST_SZ, SPACE_SZ, lbs(BLANK,POSTERIOR, PROBABILITY));
const col_hdr  _LIKELIHOOD_(NATURAL_LOG_SZ, SPACE_SZ, lbs(NATURAL_LOG, MAXIMUM, LIKELIHOOD));
const col_hdr  _TYPE_AA(GENOTYPE_SZ, SPACE_SZ, lbs(BLANK, BLANK, TYPE_AA));
const col_hdr  _TYPE_AB(GENOTYPE_SZ, SPACE_SZ, lbs(BLANK, BLANK, TYPE_AB));
const col_hdr  _TYPE_BB(GENOTYPE_SZ, SPACE_SZ, lbs(BLANK, BLANK, TYPE_BB));


// - Meta headers (Same conventions as for column headers.)
//
const col_hdr  _MLE_AVE_(MLE.size(), (2 * (RECOM_SZ + SPACE_SZ)) - MLE.size(),
                         lbs(MLE, BLANK, BLANK));
const col_hdr  _MLE_SEX_SPEC_(MLE.size(), 4 * (RECOM_SZ + SPACE_SZ) - MLE.size(),
                              lbs(MLE, BLANK, BLANK));
const col_hdr  _LOD_SCORE_(USING_RECOM_SMALL_.size(), 0,
                           lbs(LOD_SCORE_LINKAGE_, USING_RECOM_SMALL_, BLANK));
const col_hdr  _CLEVES_ELSTON_(USING_RECOM_LARGE_.size(), 0,
                               lbs(CLEVES_ELSTON_, BLANK, BLANK));
const col_hdr  _SMITHS_(USING_RECOM_SMALL_.size(), 0,
                        lbs(SMITHS_, BLANK, BLANK));
const col_hdr  _FARAWAYS_(USING_RECOM_SMALL_.size(), 0,
                          lbs(FARAWAYS_, BLANK, BLANK));                        
const col_hdr  _MORTONS_SS_(SEX_SPECIFIC_.size(), 0,
                            lbs(MORTONS_, SEX_SPECIFIC_, BLANK));
const col_hdr  _MORTONS_NON_SS_(SEX_AVERAGED_.size(), 0,
                                lbs(MORTONS_, SEX_AVERAGED_, BLANK));
const col_hdr  _LODS_AVE_(SEX_AVERAGED_.size(), 0,
                          lbs(LOD_SCORES, SEX_AVERAGED_, BLANK));
const col_hdr  _LODS_AVE_FAM_(SEX_AVERAGED_.size(), 0,
                          lbs(LOD_SCORES_FAM_, SEX_AVERAGED_, BLANK));
const col_hdr  _LODS_SEX_SPEC_(M_F_.size(), 0,
                               lbs(LOD_SCORES, SEX_SPECIFIC_, M_F_));
const col_hdr  _LODS_SEX_SPEC_FAM_(M_F_.size(), 0,
                                   lbs(LOD_SCORES_FAM_, SEX_SPECIFIC_, M_F_));
const col_hdr  _CLEVES_ELSTON_VC_(PARAM_ORDER_M_F_.size(), 0,
                                  lbs(CLEVES_ELSTON_, MATRICES_, PARAM_ORDER_M_F_));
const col_hdr  _SMITHS_NON_SS_VC_(PARAM_ORDER_A_A_.size(), 0,
                                  lbs(SMITHS_, MATRICES_, PARAM_ORDER_A_A_));
const col_hdr  _SMITHS_SS_VC_(PARAM_ORDER_M_F_A_.size(), 0,
                              lbs(SMITHS_, MATRICES_, PARAM_ORDER_M_F_A_));
const col_hdr  _FARAWAYS_NON_SS_VC_(PARAM_ORDER_A_A_.size(), 0,
                                  lbs(FARAWAYS_, MATRICES_, PARAM_ORDER_A_A_));
const col_hdr  _FARAWAYS_SS_VC_(PARAM_ORDER_M_F_A_.size(), 0,
                              lbs(FARAWAYS_, MATRICES_, PARAM_ORDER_M_F_A_));
const col_hdr  _NON_SS_LOD_SCORE_VC_(PARAM_ORDER_A_.size(), 0,
                                     lbs(LOD_SCORE_LINKAGE_TEST, MATRICES_, PARAM_ORDER_A_));
const col_hdr  _SS_LOD_SCORE_VC_(PARAM_ORDER_M_F_.size(), 0,
                                 lbs(LOD_SCORE_LINKAGE_TEST, MATRICES_, PARAM_ORDER_M_F_));                                     

void  
write_double(ostream& out, double d, bool p_value)
{
  if(SAGE::isnan(d))
  {
    out << " --- ";
  }
  else if(p_value && d < 1.0e-20)
  {
    out << "< 1.0e-20";
  }
  else
  {
    out << OUTPUT::default_rules.render(OUTPUT::Double(d));
  }
}

ostream&
operator<<(ostream& out, const header& h)
{
  ios::fmtflags old_flags = out.flags();
  
  out << left;
  out << setfill(' ');
  for(size_t line = 0; line < 3; ++line)
  {
    out << setw(h.my_offset) << "";
    
    for(size_t col = 0; col < h.my_cols.size(); ++ col)
    {
      out << setw((h.my_cols[col])->lw()) << (*(h.my_cols[col]))[line]
          << setw((h.my_cols[col])->sw()) << "";
    }
    
    out << endl;
  }
  
  // - Underline.
  //
  if(h.underline())
  {
    for(size_t col = 0; col < h.my_cols.size(); ++ col)
    {
      out << setw(h.my_offset) << "";
      out << setfill(h.my_underline) << setw((h.my_cols[col])->lw()) << "";
      out << setfill(' ') << setw((h.my_cols[col])->sw()) << "";
    }
    
    out << endl;
  }
  
  out.flags(old_flags);
  
  return  out;
}

}
}


