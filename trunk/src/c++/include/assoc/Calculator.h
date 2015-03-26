#ifndef ASSOC_CALCULATOR_H
#define ASSOC_CALCULATOR_H

#include <string>
#include <fstream>
#include "numerics/log_double.h"
#include "maxfun/maxfun.h"
#include "ibd/prior_ibd.h"
#include "assoc/Datatypes.h"
#include "assoc/MatrixDefs.h"
#include "assoc/MemberCovariateCalculator.h"
#include "assoc/Configuration.h"

namespace SAGE  {
namespace ASSOC {

class Calculator
{
  public:
    Calculator(ostream& messages, Configuration& config, MAXFUN::ParameterMgr& mgr,
               const FPED::Multipedigree& fped, const Sampledata& sampledata,
               MemberCovariateCalculator& mcc,
               MFSUBMODELS::Transformation::Facade& facade, const string& model_name,
               bool alt, cerrorstream& errors = sage_cerr);

    MemberCovariateCalculator& getMcc() { return my_mcc; }
    const MemberCovariateCalculator&  getMcc() const { return my_mcc; }
    //const vector<vector<vector<double> > >&  getSharedEffects() {return  my_shared_effects;}

    double calculateLh(MAXFUN::ParameterMgr&);
    
    const map<string, FortranMatrix<double> >&   getSharedEffects() const;
    const FortranMatrix<double>&   getRandomEffect() const;
    const vector<FPED::MemberConstPointer>&  getMemberLookup() const;  // shared effects member index (not mped index) to member pointer
    const map<FPED::MemberConstPointer, size_t>&  getIndexLookup() const;  // member pointer to shared effects member index     

  private:
    typedef vector<vector<size_t> >  IndVarCounts;

    // CALCULATOR CONSTRUCTION
    /// Sets up the PointDensity's.
    void setup_pds();
      void  createIndividualPds(map<size_t, PointDensity>& phenotype_pds);
      void  addPolygenicEffect(map<size_t, PointDensity>& phenotype_pds, size_t& term_idx);
      void  calculateSharedPolygenicEffects();
      void  addFamilyEffect(map<size_t, PointDensity>& phenotype_pds, size_t& term_idx, IndVarCounts& ind_var_counts);
        void  enumerateFamilyMembers(vector<FPED::MemberConstPointer>& valid_members, FPED::FamilyConstIterator fam_itr);
        void  setPointDensity(PointDensity& pd, const vector<FPED::MemberConstPointer>& valid_members,
                              size_t term_idx, size_t group_idx, const string& generic_name);
        void  generatePhiTerms(map<size_t, PointDensity>& phenotype_pds, 
                               const vector<FPED::MemberConstPointer>& valid_members, 
                               IndVarCounts& ind_var_counts, size_t term_idx, size_t group_idx);
      void  addMaritalEffect(map<size_t, PointDensity>& phenotype_pds, size_t& term_idx, IndVarCounts& ind_var_counts);
        void  enumerateMaritalMembers(vector<FPED::MemberConstPointer>& valid_members, FPED::FamilyConstIterator fam_itr);
      void  addSiblingEffect(map<size_t, PointDensity>& phenotype_pds, size_t& term_idx, IndVarCounts& ind_var_counts);
        void  enumerateSibshipMembers(vector<FPED::MemberConstPointer>& valid_members, FPED::FamilyConstIterator fam_itr);
      void  addUserEffects(map<size_t, PointDensity>& phenotype_pds, size_t& term_idx, IndVarCounts& ind_var_counts);
        void  addSingleUserEffect(map<size_t, PointDensity>& phenotype_pds, size_t& term_idx, 
                                  IndVarCounts& ind_var_counts, MAXFUN::ParameterConstIterator p_iter);
          void  enumerate_category_members(vector<FPED::MemberConstPointer>& valid_members, 
                                                  size_t trait_index, size_t category_index);
      void  addResidualVariances(map<size_t, PointDensity>& phenotype_pds, const IndVarCounts& ind_var_counts);
        void  findMaximumEffectCounts(IdxVector& max_effect_counts, const IndVarCounts& ind_var_counts);
      void  addIndividualTerms(const map<size_t, PointDensity>& phenotype_pds);
    
    /// Reshuffles the term indices for greater efficiency.
    void reorganize_pds();

    /// Builds all the necessary structures (PointDensityMatrice's) to evaluate a likelihood.
    void build_pdms();
 
    
    // LIKELIHOOD CALCULATION
    /// Sets up the coefficient of the first term of integration (the constant).
    void set_phenotype_coefficients();

    /// Populates all the PDMs with the current estimates from the mgr.
    log_double populate_matrices(MAXFUN::ParameterMgr& mgr);

    /// Populates the given H matrix, using the reverse_lookup and applicable_pdms list.
    void populate_h_matrix(PointDensityMatrix& H_matrix, const IdxVector& reverse_lookup, const IdxSet& applicable_pdms) const;

    /// Populates a correctly sized R matrix from an H matrix.
    void populate_r_matrix(const PointDensityMatrix& H_matrix, PointDensityMatrix& R_matrix) const;

    // debugging only
    void dump_pdms(const PDMVector& pdms) const;

    void dump_matrix(const string& name, const Matrix& matrix) const;
    
    // - For calculating independent residuals
    //
    void setUpIndices();
    void  dumpMemberLookups() const;
    void initializeSharedEffects();
    void initializeRandomEffect();
    void populateEffectMatrix(const string& effect, vector<FPED::MemberConstPointer> valid_members);
    void dumpSharedEffects() const;
    
    //===================================================================
    /// Data members:
    //===================================================================

    ostream&      my_messages;
    cerrorstream  my_errors;
    
    const Configuration&  my_config;
    const string          my_model_name;
    bool                  alt_covariates;
    
    MAXFUN::ParameterMgr&  my_mgr;
    
    const FPED::Multipedigree&  my_fped;
    const Sampledata&           my_sampledata;
    
    MemberCovariateCalculator&  my_mcc;
    
    size_t           my_term_count;      // The number of terms of integration
    TermNameVector   my_term_names;      // A vector where for each index i:
                                         //   i is the index of a term-of-integration
                                         //   the corresponding string is the user-friendly name for that term    
    PointDensities       my_pds;
    IdxVector            my_phenotype_terms;   // A vector where for each index i:
                                               //   i is the mpindex of an individual
                                               //   the corresponding size_t is the index of the adjusted phenotype 
                                               //     point density term my_pdms.
    TermLookupTable        my_lookups;    // A vector of sets of indices.


    PDMVector              my_pdms;            // A vector of PointDensityMatrix's.
    PDMWithLookupVector    my_h_matrices;      // A vector of H matrices (where each index is the term-of-integration)
    IdxVector              my_r_matrices;      // A vector where for each index i:
                                               //    i is the index of a term-of-integration
                                               //    the corresponding value is the index of the R matrix 
                                               //      for that term-of-integration
    IdxSet     my_constant_matrices;    // A set where each value i is a PDM that corresponds to a matrix with 
                                        // one cell for the constant
                                        
    // - number of variance effects shared by any two members (member i x member j x variance component)
    //   indices are specific to this data structure.
    map<string, FortranMatrix<double> >   my_shared_effects;   // effect name, matrix of shared effects
    FortranMatrix<double>  my_random_effect;
    vector<FPED::MemberConstPointer>  my_member_lookup;  // shared effects member index (not mped index) to member pointer
    map<FPED::MemberConstPointer, size_t>  my_index_lookup;  // member pointer to shared effects member index 
};

#include "assoc/Calculator.ipp"

} 
} 

#endif
