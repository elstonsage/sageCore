//#define DEBUG_CALC 1

//=======================================================================
//
//  File:  Calculator.cpp
//
//  Author:  Stephen Gross
//
//  Copyright 2002 R. C. Elston
//=======================================================================

#include "assoc/Calculator.h"
#include "numerics/normal_pdf.h"

namespace SAGE  {
namespace ASSOC {


//=======================================================================
//  Calculator() 
//=======================================================================
Calculator::Calculator(ostream& messages, Configuration& config,
                       MAXFUN::ParameterMgr& mgr, const FPED::Multipedigree& fped,
                       const Sampledata& sampledata,
                       MemberCovariateCalculator& mcc,
                       MFSUBMODELS::Transformation::Facade& facade, const string& model_name,
                       bool alt, cerrorstream& errors)
      : my_messages(messages), my_errors(errors), my_config(config), my_model_name(model_name),  
        alt_covariates(alt), my_mgr(mgr), my_fped(fped), my_sampledata(sampledata), 
        my_mcc(mcc /*config, mgr, sampledata, facade*/)
{ 
  //size_t  component_count = my_mgr.getParamCount("Variance components");

  // - For calculating independent residuals.
  //
  setUpIndices();
  //dumpMemberLookups();
  initializeSharedEffects();
  initializeRandomEffect();
  setup_pds();
  //dumpSharedEffects();
  //print_matrix(my_random_effect, cout);
  
  reorganize_pds();
  build_pdms();
}

void
Calculator::dumpSharedEffects() const
{
  cout << endl;

  map<string, FortranMatrix<double> >::const_iterator  e_iter = my_shared_effects.begin();
  map<string, FortranMatrix<double> >::const_iterator  e_end_iter = my_shared_effects.end();
  for(; e_iter != e_end_iter; ++e_iter)
  {
    cout << e_iter->first << endl;
    print_matrix(e_iter->second, cout);
  }
}

void
Calculator::initializeRandomEffect()
{
  size_t  member_count = my_member_lookup.size();
  
  my_random_effect.resize_fill(member_count, member_count, 0.0);
  for(size_t i = 0; i < member_count; ++i)
  {
    for(size_t j = 0; j < member_count; ++j)
    {
      if(i == j)
      {
        my_random_effect(i, j) = 1.0;
      }
    }
  }
}

void 
Calculator::initializeSharedEffects()
{
  size_t  member_count = my_member_lookup.size();
  MAXFUN::ParameterConstIterator  vc_iter     = my_mgr.getParamBegin("Variance components");
  MAXFUN::ParameterConstIterator  vc_end_iter = my_mgr.getParamEnd("Variance components");
  for(; vc_iter != vc_end_iter; ++vc_iter)
  {
    string  param_name = vc_iter->getName();
    if(! (param_name == "0.5 * Polygenic" || param_name == "Random"))
    {
      FortranMatrix<double>  effect;
      effect.resize_fill(member_count, member_count, 0.0);
      my_shared_effects.insert(make_pair(param_name, effect));
    }  
  }  
}

void 
Calculator::populateEffectMatrix(const string& effect, vector<FPED::MemberConstPointer> valid_members)
{
  FortranMatrix<double>& effect_matrix = my_shared_effects[effect];
  vector<FPED::MemberConstPointer>::const_iterator  row_iter = valid_members.begin();
  vector<FPED::MemberConstPointer>::const_iterator  row_end_iter = valid_members.end();
  for(; row_iter != row_end_iter; ++row_iter)
  {
    vector<FPED::MemberConstPointer>::const_iterator  col_iter = valid_members.begin();
    vector<FPED::MemberConstPointer>::const_iterator  col_end_iter = valid_members.end();
    for(; col_iter != col_end_iter; ++col_iter)
    {
      effect_matrix(my_index_lookup[*row_iter], my_index_lookup[*col_iter]) = 1.0;
    }
  }
}

// - Populate member lookup containers.
//
void
Calculator::setUpIndices()
{
  size_t  member_count = my_sampledata.getTotalIndividualCount();
  size_t  informative_index = 0;
  for(size_t mped_index = 0; mped_index < member_count; ++mped_index)
  {
    if(my_sampledata.isValid(mped_index))
    {
      FPED::MemberConstPointer  member_ptr = &(my_fped.member_index(mped_index));
      my_member_lookup.push_back(member_ptr);
      my_index_lookup.insert(make_pair(member_ptr, informative_index++));
    }
  }
}

void
Calculator::dumpMemberLookups() const
{
  cout << endl;
  size_t  valid_member_count = my_member_lookup.size();
  for(size_t i = 0; i < valid_member_count; ++i)
  {
    FPED::MemberConstPointer  member_ptr = my_member_lookup[i];
    cout << "shared effects index  "<< i << ", mp index " << member_ptr->mpindex() << " -> pedigree " 
         << member_ptr->pedigree()->name() << "  member " << member_ptr->name() << endl;
  }
  
  cout << endl;
  size_t  total_member_count = my_fped.member_count();
  for(size_t i = 0; i < total_member_count; ++i)
  {
    FPED::MemberConstPointer  member_ptr = &(my_fped.member_index(i));
    map<FPED::MemberConstPointer, size_t>::const_iterator lookup_result = my_index_lookup.find(member_ptr);
    size_t  shared_effects_index = lookup_result == my_index_lookup.end() ? static_cast<size_t>(-1) : lookup_result->second;
    cout << "mped index " << i << ", shared_effects index " << shared_effects_index << endl;
  }
}


//=======================================================================
//  setup_pds()
//=======================================================================
void
Calculator::setup_pds()
{
  my_pds.clear();
  
  // Clear and set up the first named term-of-integration.
  my_term_names.resize(1, "Constant");

  // Term-of-integration counter:
  // NOTE: We start at 1, not 0, since 0 is for the constant portion.
  size_t  term_idx = 1;

  // A vector where for each index i:
  //   i is the mpindex of an individual
  //   the corresponding value is a vector where for each index j:
  //     j is the groupindex of a variance component
  //     the corresponding value is the COUNT of the individual's phenotypic 
  //     additive terms specific to that variance component
  //
  // For instance, let's say individual X has two family effects. Then IndVarCounts[X][FamilyEffectIdx] = 2.
  IndVarCounts  ind_var_counts(my_sampledata.getTotalIndividualCount(), 
                              vector<size_t> (my_mgr.getParamCount("Variance components"), 0));

  // A map to temporarily hold the point density terms for the effect-specific terms.
  map<size_t, PointDensity> phenotype_pds;

  createIndividualPds(phenotype_pds);

  if(my_config.modelHasFixedEffect(my_model_name, alt_covariates, Configuration::POLYGENIC))
  {
    addPolygenicEffect(phenotype_pds, term_idx);
    calculateSharedPolygenicEffects(); 
  }
    
  if(my_config.modelHasFixedEffect(my_model_name, alt_covariates, Configuration::FAMILY))
    addFamilyEffect(phenotype_pds, term_idx, ind_var_counts);
    
  if(my_config.modelHasFixedEffect(my_model_name, alt_covariates, Configuration::SIBLING))
    addSiblingEffect(phenotype_pds, term_idx, ind_var_counts);
    
  if(my_config.modelHasFixedEffect(my_model_name, alt_covariates, Configuration::MARITAL))
    addMaritalEffect(phenotype_pds, term_idx, ind_var_counts);
    
  addUserEffects(phenotype_pds, term_idx, ind_var_counts);
  
  my_term_count = term_idx;
  addResidualVariances(phenotype_pds, ind_var_counts);
  addIndividualTerms(phenotype_pds);
}


void
Calculator::calculateSharedPolygenicEffects()
{
  FortranMatrix<double>& effect_matrix = my_shared_effects["Polygenic"];

  size_t  pedigree_count = my_fped.pedigree_count();
  for(size_t p = 0; p < pedigree_count; ++p)
  {
    const SAGE::FPED::Pedigree&  pedigree = my_fped.pedigree_index(p);
    RefPriorIBD  rpi(pedigree);
    
    rpi.compute(pedigree);
    size_t  member_count = pedigree.member_count();
    for(size_t i = 0; i < member_count; ++i)
    {
      for(size_t j = 0; j < member_count; ++j)
      {
        double  f1 = rpi.f(i, j, 1);
        double  f2 = rpi.f(i, j, 2);
        double  coefficient_of_relationship = f2 + 0.5 * f1;
        
        size_t  shared_effects_index_i = my_index_lookup[&(pedigree.member_index(i))];
        size_t  shared_effects_index_j = my_index_lookup[&(pedigree.member_index(j))];
        
        effect_matrix(shared_effects_index_i, shared_effects_index_j) = coefficient_of_relationship;
      }
    }
  }
}


/*
void
Calculator::dump_shared_effects() const
{
  size_t  member_count = my_shared_effects.size();
  for(size_t i = 0; i < member_count; ++i)
  {
    for(size_t j = 0; j < i; ++j)
    {
      string  member1 = my_fped.member_index(i).pedigree()->name() + ":" + my_fped.member_index(i).name() + " ";
      string  member2 = my_fped.member_index(j).pedigree()->name() + ":" + my_fped.member_index(j).name() + " ";
      
      cout << member1 << ", " << member2 << endl;      
      
      size_t  component_count = my_mgr.getParamCount("Variance components");
      for(size_t vc = 0; vc < component_count; ++vc)
      {
        string  component_name = my_mgr.getParameter("Variance components", vc).getName();
        
        cout << component_name << "  " << my_shared_effects[i][j][vc] << ", ";
      }
      
      cout << "\n" << endl;
    }
  }
}
*/


void
Calculator::addIndividualTerms(const map<size_t, PointDensity>& phenotype_pds)
{
  my_phenotype_terms.resize(my_sampledata.getTotalIndividualCount(), size_t(-1));
  
  map<size_t, PointDensity>::const_iterator  pd_iter = phenotype_pds.begin();
  map<size_t, PointDensity>::const_iterator  pd_end_iter = phenotype_pds.end();
  for(; pd_iter != pd_end_iter; ++pd_iter)
  {
    my_pds.push_back(pd_iter->second);
    my_phenotype_terms[pd_iter->first] = my_pds.size() - 1;
  }
}


void
Calculator::addResidualVariances(map<size_t, PointDensity>& phenotype_pds, const IndVarCounts& ind_var_counts)
{
  IdxVector  max_counts(my_mgr.getParamCount("Variance components"), 0);

  findMaximumEffectCounts(max_counts, ind_var_counts);

  // Loop across individuals and add residual variance.
  map<size_t, PointDensity>::iterator  pd_iter     = phenotype_pds.begin();
  map<size_t, PointDensity>::iterator  pd_end_iter = phenotype_pds.end();
  for(; pd_iter != pd_end_iter; ++pd_iter)
  {
    // Loop across variance components.
    size_t  variance_component_count = max_counts.size();
    for(size_t var = 0; var < variance_component_count; ++var)
    {
      size_t difference = max_counts[var] - ind_var_counts[pd_iter->first][var];
      
      for(size_t x = 0; x < difference; ++x)
        pd_iter->second.var_idxs.push_back(var);
    }
  }
}


// Find the maximum count amoung individuals for each variance component, ie
// what is the largest number of families, sibships, etc that an individual
// belongs to.
//
void
Calculator::findMaximumEffectCounts(IdxVector& max_effect_counts, const IndVarCounts& ind_var_counts)
{
  size_t  individual_count = ind_var_counts.size();
  size_t  component_count  = max_effect_counts.size();
  for(size_t i = 0; i < individual_count; ++i)
    for(size_t var = 0; var < component_count; ++var)
      if(ind_var_counts[i][var] > max_effect_counts[var])
        max_effect_counts[var] = ind_var_counts[i][var];
}


void
Calculator::createIndividualPds(map<size_t, PointDensity>& phenotype_pds)
{
  PointDensity  phenotype_pd;

  // Set up a CoeffPair for a coefficient of the constant.
  // Note: During maximization, this value will be set, so the -999 doesn't matter for the moment.
  phenotype_pd.coeff_pairs.push_back(CoeffPair(-999, 0));

  // Add in the random effect.
  phenotype_pd.var_idxs.push_back(my_mgr.getParameter("Variance components", "Random").getGroupIndex());

  // Add in a point density term for every valid individual.
  size_t  individual_count = my_sampledata.getTotalIndividualCount();
  for(size_t i = 0; i < individual_count; ++i)
  {
    if(my_sampledata.isValid(i))
    {
      phenotype_pd.name = "Adj. pheno. for " + 
                          my_sampledata.getIndividual(i).pedigree()->name() + 
                          ":" + my_sampledata.getIndividual(i).name();
      phenotype_pds[my_sampledata.getIndividual(i).mpindex()] = phenotype_pd;
    }
  }
}

void
Calculator::addUserEffects(map<size_t, PointDensity>& phenotype_pds, 
                                    size_t& term_idx, IndVarCounts& ind_var_counts)
{
  MAXFUN::ParameterConstIterator  p_iter     = my_mgr.getParamBegin("Variance components");
  MAXFUN::ParameterConstIterator  p_end_iter = my_mgr.getParamEnd("Variance components");  
  for(; p_iter != p_end_iter; ++p_iter)
  {
    string  param_name = p_iter->getName();
    
    if(! (fixed_effect(param_name) || param_name == "0.5 * Polygenic"))
    {
      addSingleUserEffect(phenotype_pds, term_idx, ind_var_counts, p_iter);
    }
  }                 
}


void  
Calculator::addSingleUserEffect(map<size_t, PointDensity>& phenotype_pds, size_t& term_idx, 
                                IndVarCounts& ind_var_counts, MAXFUN::ParameterConstIterator p_iter)
{
  size_t group_idx = p_iter->getGroupIndex();

  string  categorical_trait_name = p_iter->getName();
  RPED::RefMPedInfo  rmpi = my_fped.info();
  size_t  trait_index = rmpi.trait_find(categorical_trait_name);
  
  assert(trait_index != (size_t)(-1));
  
  const RPED::RefTraitInfo&  trait_info = rmpi.trait_info(trait_index);
  
  assert(trait_info.type() == RPED::RefTraitInfo::categorical_trait);
  
  const vector<string>&  categories = trait_info.get_categories();
  size_t  category_count = categories.size();
  for(size_t c = 0; c < category_count; ++c)
  {
    string  category = categories[c];
    
    vector<FPED::MemberConstPointer>  valid_members;
    enumerate_category_members(valid_members, trait_index, c);
    populateEffectMatrix(categorical_trait_name, valid_members);
    
    // If there are fewer than two valid members, skip the category.
    size_t  member_count = valid_members.size();
    if(member_count < 2)
      continue;

    // There are at least two valid members; add the effect.
    PointDensity category_pd;
    
    setPointDensity(category_pd, valid_members, term_idx, group_idx, categorical_trait_name);
    generatePhiTerms(phenotype_pds, valid_members, ind_var_counts, term_idx, group_idx);
    ++term_idx;
    my_pds.push_back(category_pd);
  }   
}    
                            
                       
void  
Calculator::enumerate_category_members(vector<FPED::MemberConstPointer>& valid_members, 
                                                     size_t trait_index, size_t category_index)
{
  FPED::PedigreeConstIterator  ped_iter     = my_fped.pedigree_begin();
  FPED::PedigreeConstIterator  ped_end_iter = my_fped.pedigree_end();
  for(; ped_iter != ped_end_iter; ++ped_iter)
  {
    const FPED::FilteredPedigreeInfo&  fpi = ped_iter->info();
  
    FPED::MemberConstIterator  member_iter     = ped_iter->member_begin();
    FPED::MemberConstIterator  member_end_iter = ped_iter->member_end();
    for(; member_iter != member_end_iter; ++member_iter)
    {
      size_t  member_index = member_iter->mpindex();
      if(my_sampledata.isValid(member_index) &&
         fpi.trait(member_iter->index(), trait_index) == category_index)
      {
        valid_members.push_back(&(*member_iter));
      } 
    }
  }
}
                            

void
Calculator::addFamilyEffect(map<size_t, PointDensity>& phenotype_pds, size_t& term_idx, IndVarCounts& ind_var_counts)
{
  size_t family_idx = my_mgr.getParameter("Variance components", "Family").getGroupIndex();
  
  // Loop across pedigrees.
  FPED::PedigreeConstIterator  ped_iter     = my_fped.pedigree_begin();
  FPED::PedigreeConstIterator  ped_end_iter = my_fped.pedigree_end();
  for(; ped_iter != ped_end_iter; ++ped_iter)
  {
    // Loop across families.
    FPED::FamilyConstIterator fam_iter     = ped_iter->family_begin();
    FPED::FamilyConstIterator fam_end_iter = ped_iter->family_end();
    for(; fam_iter != fam_end_iter; ++fam_iter)
    {
      vector<FPED::MemberConstPointer>  valid_members;
      
      enumerateFamilyMembers(valid_members, fam_iter);
      populateEffectMatrix("Family", valid_members);
      
      // If there are fewer than two valid members, skip the family.
      size_t  member_count = valid_members.size();
      if(member_count < 2)
        continue;

      // There are at least two valid members; add the effect.
      PointDensity family_pd;
      
      setPointDensity(family_pd, valid_members, term_idx, family_idx, "Family");
      generatePhiTerms(phenotype_pds, valid_members, ind_var_counts, term_idx, family_idx);
      my_pds.push_back(family_pd);
      ++term_idx;
    }   
  }   
} 

void
Calculator::enumerateFamilyMembers(vector<FPED::MemberConstPointer>& valid_members, FPED::FamilyConstIterator fam_iter)
{
  if(my_sampledata.isValid(fam_iter->parent1()->mpindex()))
    valid_members.push_back(fam_iter->parent1());
  
  if(my_sampledata.isValid(fam_iter->parent2()->mpindex()))
    valid_members.push_back(fam_iter->parent2());
  
  FPED::OffspringConstIterator  offspr_iter     = fam_iter->offspring_begin();
  FPED::OffspringConstIterator  offspr_end_iter = fam_iter->offspring_end();
  for(; offspr_iter != offspr_end_iter; ++offspr_iter)
    if(my_sampledata.isValid(offspr_iter->mpindex()))
      valid_members.push_back(&*offspr_iter);
}

  
void  
Calculator::setPointDensity(PointDensity& pd, const vector<FPED::MemberConstPointer>& valid_members,
                            size_t term_idx, size_t group_idx, const string& generic_name)
{
  string member_names = "";
  
  size_t  member_count = valid_members.size();
  for(size_t i = 0; i < member_count; ++i)
  {
    member_names += valid_members[i]->pedigree()->name() + ":" + valid_members[i]->name() + " ";
    
    /*
    for(size_t j = 0; j < member_count; ++j)
    {
      if(i != j)
      {
        size_t  i_index = valid_members[i]->mpindex();
        size_t  j_index = valid_members[j]->mpindex();
        
        my_shared_effects[i_index][j_index][group_idx] = 1.0;
      }
    }
    */
    
  }
    
  my_term_names.push_back(generic_name + " effect shared by " + member_names);

  pd.name = my_term_names[term_idx];
  pd.coeff_pairs.push_back(CoeffPair(1.0, term_idx));
  pd.var_idxs.push_back(group_idx);
}


// Loop across valid individuals; for each one generate a new phi term, 
// and add the additive effect to the individual's phenotype phi term.
void  
Calculator::generatePhiTerms(map<size_t, PointDensity>& phenotype_pds, 
                             const vector<FPED::MemberConstPointer>& valid_members,
                             IndVarCounts& ind_var_counts, size_t term_idx, size_t group_idx)
{
  vector<FPED::MemberConstPointer>::const_iterator  ind_iter     = valid_members.begin();
  vector<FPED::MemberConstPointer>::const_iterator  ind_end_iter = valid_members.end();
  for(; ind_iter != ind_end_iter; ++ind_iter)
  {
    // Add the additive effect to the individual.
    phenotype_pds[(*ind_iter)->mpindex()].coeff_pairs.push_back(CoeffPair(1.0, term_idx));

    // Increment the count of additive family terms for this individual.
    ++(ind_var_counts[(*ind_iter)->mpindex()][group_idx]);
  } 
}

void
Calculator::addSiblingEffect(map<size_t, PointDensity>& phenotype_pds, size_t& term_idx, IndVarCounts& ind_var_counts)
{
  size_t sibling_idx = my_mgr.getParameter("Variance components", "Sibling").getGroupIndex();
  
  // Loop across pedigrees.
  FPED::PedigreeConstIterator  ped_iter     = my_fped.pedigree_begin();
  FPED::PedigreeConstIterator  ped_end_iter = my_fped.pedigree_end();
  for(; ped_iter != ped_end_iter; ++ped_iter)
  {
    // Loop across families.
    FPED::FamilyConstIterator  fam_iter     = ped_iter->family_begin();
    FPED::FamilyConstIterator  fam_end_iter = ped_iter->family_end();
    for(; fam_iter != fam_end_iter; ++fam_iter)
    {
      // Make a list of all the valid members:
      vector<FPED::MemberConstPointer>  valid_members;
    
      enumerateSibshipMembers(valid_members, fam_iter);
      populateEffectMatrix("Sibling", valid_members);    // For independent residuals.
      
      // Make sure there are at least two valid members.
      size_t  member_count = valid_members.size();
      if(member_count < 2)
        continue;

      // There are at least two valid members; add the effect:
      PointDensity  sibling_pd;
      
      setPointDensity(sibling_pd, valid_members, term_idx, sibling_idx, "Sibling");
      generatePhiTerms(phenotype_pds, valid_members, ind_var_counts, term_idx, sibling_idx);
      my_pds.push_back(sibling_pd);
      ++term_idx;      
    }    
  }    
}


void
Calculator::enumerateSibshipMembers(vector<FPED::MemberConstPointer>& valid_members, FPED::FamilyConstIterator fam_iter)
{
  FPED::OffspringConstIterator  offspr_iter = fam_iter->offspring_begin();
  FPED::OffspringConstIterator  offspr_end_iter = fam_iter->offspring_end();
  for(; offspr_iter != offspr_end_iter; ++offspr_iter)
    if(my_sampledata.isValid(offspr_iter->mpindex()))
      valid_members.push_back(&*offspr_iter);
}


void
Calculator::addMaritalEffect(map<size_t, PointDensity>& phenotype_pds, size_t& term_idx, IndVarCounts& ind_var_counts)
{
  size_t marital_idx = my_mgr.getParameter("Variance components", "Marital").getGroupIndex();
  
  // Loop across pedigrees.
  FPED::PedigreeConstIterator  ped_iter     = my_fped.pedigree_begin();
  FPED::PedigreeConstIterator  ped_end_iter = my_fped.pedigree_end();
  for(; ped_iter != ped_end_iter; ++ped_iter)
  {
    // Loop across families.
    FPED::FamilyConstIterator  fam_iter     = ped_iter->family_begin();
    FPED::FamilyConstIterator  fam_end_iter = ped_iter->family_end();
    for(; fam_iter != fam_end_iter; ++fam_iter)
    {
      // Make a list of all the valid parents.
      vector<FPED::MemberConstPointer>  valid_members;
      
      enumerateMaritalMembers(valid_members, fam_iter);
      populateEffectMatrix("Marital", valid_members);
      
      // Make sure there are two valid members.
      if(valid_members.size() != 2)
        continue;

      // There are two persons; add the effect.
      //size_t  marital_term = term_idx++;
      PointDensity  marital_pd;
      
      setPointDensity(marital_pd, valid_members, term_idx, marital_idx, "Marital");
      generatePhiTerms(phenotype_pds, valid_members, ind_var_counts, term_idx, marital_idx);
      my_pds.push_back(marital_pd);      
      ++term_idx;
    }    
  }   
} 


void
Calculator::enumerateMaritalMembers(vector<FPED::MemberConstPointer>& valid_members, FPED::FamilyConstIterator fam_iter)
{
  if(my_sampledata.isValid(fam_iter->parent1()->mpindex()))
    valid_members.push_back(fam_iter->parent1());
  
  if(my_sampledata.isValid(fam_iter->parent2()->mpindex()))
    valid_members.push_back(fam_iter->parent2());
}


void
Calculator::addPolygenicEffect(map<size_t, PointDensity>& phenotype_pds, size_t& term_idx)
{
  // Set up a structure to record the polygenic term-of-integration for all the founders:
  map<size_t, size_t> polygenic_terms;

  // Loop across all individuals and record terms-of-integration for each.
  size_t  individual_count = my_sampledata.getTotalIndividualCount();
  for(size_t i = 0; i < individual_count; ++i)
  {
    // Add a user-friendly name for this term, and set the term-of-integration for this individual:
    my_term_names.push_back("Polygenic effect for " + 
                            my_sampledata.getIndividual(i).pedigree()->name() + 
                            ":" + my_sampledata.getIndividual(i).name());
    polygenic_terms[my_sampledata.getIndividual(i).mpindex()] = term_idx++;
  }      
  
  // Loop across all individuals and add PD's:
  for(size_t i = 0; i < individual_count; ++i)
  {
    // Get the member and the polygenic term-of-integration:
    const FPED::Member&  ind = my_sampledata.getIndividual(i);
    size_t  polygenic_term = polygenic_terms.find(ind.mpindex())->second;
    
    // Set the name for this point density term as well as it's vector of coefficient pairs:
    PointDensity  pd;

    pd.name = "Polygenic effect for " + ind.pedigree()->name() + ":" + ind.name();
    pd.coeff_pairs.push_back(CoeffPair(1.0, polygenic_term));

    // If it's a nonfounder add the polygenic term for his/her parents to the coefficient pair list:
    if(ind.is_founder())
    {
      pd.var_idxs.push_back(my_mgr.getParameter("Variance components", "Polygenic").getGroupIndex());
    }
    else   // Nonfounder
    {
      pd.var_idxs.push_back(my_mgr.getParameter("Variance components", "0.5 * Polygenic").getGroupIndex());

      pd.coeff_pairs.push_back(CoeffPair(-0.5, polygenic_terms.find(ind.parent1()->mpindex())->second));
      pd.coeff_pairs.push_back(CoeffPair(-0.5, polygenic_terms.find(ind.parent2()->mpindex())->second));
    }
    
    // Add the PointDensity:
    my_pds.push_back(pd);
    
    // If the individual is valid, record the term as an additive adjustment to that individual's phenotype:
    // Note that we're NOT recording the count of additive polygenic terms for each
    // individual; this is because as long as polygenic is enabled, it must apply once and only once
    // to every valid individual.
    if(my_sampledata.isValid(i))
    {
      phenotype_pds[ind.mpindex()].coeff_pairs.push_back(CoeffPair(1.0, polygenic_term));
    }
  }
}


struct SortPairByValue 
{
  bool operator() (const pair<size_t, size_t>& x, const pair<size_t, size_t>& y) const 
  { 
    return x.second > y.second; 
  } 
};

//=======================================================================
//  reorganize_pds()
//=======================================================================
void
Calculator::reorganize_pds()
{
  // Ok, what we need to do at this point is look at the efficiency
  // of the organization of terms-of-integration. Specifically, we want
  // to see, for some term i, how many pds have non-zero coefficients of i.
  // Then, given that we still integrate from last term to first, we want
  // to first integrate out the terms-of-integration that have the FEWEST
  // corresponding pds.
  
  // A vector where for each index i:
  //   i is a term-of-integration
  //   the corresponding value is the number of pds that have a non-zero coefficient of i
  IdxVector pd_counts(my_term_count, 0);
  
  // Populate pd_counts.
  size_t  pd_count = my_pds.size();
  for(size_t pd = 0; pd < pd_count; ++pd)
  {
    CoeffPairVector::const_iterator  coeff_iter     = my_pds[pd].coeff_pairs.begin();
    CoeffPairVector::const_iterator  coeff_end_iter = my_pds[pd].coeff_pairs.end();
    for(; coeff_iter != coeff_end_iter; ++coeff_iter)
      ++(pd_counts[coeff_iter->term_idx]);
  }
  
  // Now sort it all.
  typedef multiset<pair<size_t, size_t>, SortPairByValue> SortedPairs;
  
  SortedPairs sorted_terms;

  // Populate the sorted_terms (notice that we SKIP term #0, because that's the constant term).
  for(size_t term = 1; term < pd_counts.size(); ++term)
    sorted_terms.insert(make_pair(term, my_config.reverse_sort ? term : pd_counts[term]));

  // Now create a lookup table, where for each index i:
  //   i is a term-of-integration
  //   the corresponding value is the new index for the term-of-integration
  IdxVector new_term_idxs(my_term_count, 0);
  
  // Populate it.
  size_t new_idx = 1;
  
  for(SortedPairs::const_iterator i = sorted_terms.begin(); i != sorted_terms.end(); ++i)
    new_term_idxs[i->first] = new_idx++;
  
  // Now let's fix up the TermNameVector:
  TermNameVector new_term_names(my_term_count);
  
  for(size_t i = 0; i < my_term_names.size(); ++i)
    new_term_names[new_term_idxs[i]] = my_term_names[i];
  
  my_term_names = new_term_names;

  // Now loop across all the pds and fix up the coeff pairs.
  for(size_t i = 0; i < pd_count; ++i)
  {
    CoeffPairVector::iterator  coeff_iter     = my_pds[i].coeff_pairs.begin();
    CoeffPairVector::iterator  coeff_end_iter = my_pds[i].coeff_pairs.end();
    for(; coeff_iter != coeff_end_iter; ++coeff_iter)
      coeff_iter->term_idx = new_term_idxs[coeff_iter->term_idx];
  }

  #ifdef DEBUG_CALC
    cout << "PointDensity's:" << endl; 
    for(size_t i = 0; i < my_pds.size(); ++i) 
    { 
      my_pds[i].dump(my_term_names); 
    } 
    
    cout << "END" << endl;
  #endif
}

//=======================================================================
//  build_pdms()
//=======================================================================
void
Calculator::build_pdms()
{
  #ifdef DEBUG_CALC
    cout << "Building PDMs..." << endl;
  #endif

  // Start by clearing it.
  my_pdms.clear();
  my_h_matrices.resize(my_term_count);
  my_r_matrices.resize(my_term_count,(size_t)-1);
  my_constant_matrices.clear();
  
  // Lookup table for the terms:
  my_lookups.resize(my_term_count);

  // Add matrices built from the PointDensity's.
  size_t  pd_count = my_pds.size();
  for(size_t i = 0; i < pd_count; ++i)
  {
    const PointDensity& pd = my_pds[i];
    
    // Insert a pre-sized available PDM into the PDM vector.
    my_pdms.push_back(generate_pdm(pd));
    
    // If this term represents a constant (has no coefficients for any terms-of-integration), then
    // we should record it as a single cell matrix.
    if(pd.coeff_pairs.size() == 1 && pd.coeff_pairs[0].term_idx == 0)
    {
      my_constant_matrices.insert(my_pdms.size() - 1);
    }
    else
    {
      // Loop across all the PD's CoeffPair's, and record in the lookup table the fact that every
      // term-of-integration in the CoeffPairVector has non-zero values in the newly added PDM.
      CoeffPairVector::const_iterator  coeff_iter     = pd.coeff_pairs.begin();
      CoeffPairVector::const_iterator  coeff_end_iter = pd.coeff_pairs.end();
      for(; coeff_iter != coeff_end_iter; ++coeff_iter)
        my_lookups[coeff_iter->term_idx].insert(my_pdms.size() - 1);
    }
  }

  // Loop through every term of integration (from the last to the first) and
  // perform the matrix magic.
  for(size_t term_idx = my_lookups.size() - 1; term_idx > 0; --term_idx)
  {
    #ifdef DEBUG_CALC
      cout << "  For term #" << term_idx << " (" << my_term_names[term_idx] << ")..." << endl;
    #endif

    // Figure out which terms are going into the H matrix.
    IdxSet idx_set;                                  // Will store the absolute indices that go into the H matrix
    size_t n       = 0;                              // The number of valid (available) A matrices.
  
    // Loop across A matrices to populate the reverse lookup.
    IdxSet::const_iterator idx_iter     = my_lookups[term_idx].begin();
    IdxSet::const_iterator idx_end_iter = my_lookups[term_idx].end();
    for(; idx_iter != idx_end_iter; ++idx_iter)
    {
      // Get the id of the A matrix:
      size_t a_matrix_idx = *idx_iter;
      
      // Add all the terms from this A matrix into the H matrix's idx set.
      IdxVector::const_iterator sub_iter     = my_pdms[a_matrix_idx].get_idxs().begin();
      IdxVector::const_iterator sub_end_iter = my_pdms[a_matrix_idx].get_idxs().end();
      for(; sub_iter != sub_end_iter; ++sub_iter)
        idx_set.insert(*sub_iter);
        
      // Increment the count of A matrices.
      ++n;
      
      // Now remove all references to this A matrix from all other term-of-integration lookups.
      for(size_t i = term_idx - 1; i > 0; --i)
        my_lookups[i].erase(a_matrix_idx);
    }

    // Skip this term if there are no available A matrices for it, or the index list is empty:
    if(n == 0 || idx_set.empty())
      continue;
      
    // Populate the lookup tables:
    PointDensityMatrix&  h_matrix       = my_h_matrices[term_idx].first;
    IdxVector&           h_reverse_idxs = my_h_matrices[term_idx].second;

    h_matrix.set_name("H matrix incorporating term " + my_term_names[term_idx]);
    h_matrix.get_matrix().resize_fill(idx_set.size(), idx_set.size(), 0.0);
    h_reverse_idxs.resize(my_term_count, (size_t)-1);
    h_matrix.get_idxs().clear();
    
    IdxSet::const_iterator  i     = idx_set.begin();
    IdxSet::const_iterator  i_end = idx_set.end();
    for(; i != i_end; ++i)
    {
      h_matrix.get_idxs().push_back(*i);
      h_reverse_idxs[*i] = h_matrix.get_idxs().size() - 1;
    }

    // Get the number of terms in the H matrix:
    size_t h_terms = h_matrix.get_matrix().rows();

    // If the H-matrix is at least 2-by-2, then generate an R matrix:
    if(h_terms > 1)
    {
      // Create the R matrix:
      PointDensityMatrix R_matrix("R matrix consolidating term " + my_term_names[term_idx]);
    
      // Resize it accordingly:
      R_matrix.get_matrix().resize_fill(h_terms - 1, h_terms - 1, 0.0);

      // Copy over the H matrix's idxs, leaving out the last one:
      for(size_t i = 0; i < h_matrix.get_idxs().size() - 1; ++i)
        R_matrix.get_idxs().push_back(h_matrix.get_idxs()[i]);

      // Add it to the vector:
      my_pdms.push_back(R_matrix);
      
      // Get the id:
      size_t r_matrix_id = my_pdms.size() - 1;
    
      // Record it in the lookup table:
      for(size_t i = 0; i < R_matrix.get_idxs().size(); ++i)
        my_lookups[R_matrix.get_idxs()[i]].insert(r_matrix_id);
        
      // Also record it in the R matrix lookup table:
      my_r_matrices[term_idx] = r_matrix_id;
      
      // If this is a 1-by-1 R matrix with index 0 corresponding to term-of-integration 0, then
      // we know it needs to be treated as a residual constant value. Record it!
      if(R_matrix.get_idxs().size() == 1 && R_matrix.get_idxs()[0] == 0)
        my_constant_matrices.insert(r_matrix_id);

    } // End if-h-matrix-big-enough
  } // End term-of-integration loop

  #ifdef DEBUG_CALC
    cout << "Finished building PDMs." << endl;
  #endif
}

//=======================================================================
//  dump_pdms()
//=======================================================================
void
Calculator::dump_pdms(const PDMVector& pdms) const
{
  // Debug output:
  for(PDMVector::const_iterator i = pdms.begin(); i != pdms.end(); ++i)
  {
    i->dump(my_term_names);
  }
}


void 
Calculator::dump_matrix(const string & name, const Matrix & matrix) const
{
  // dump matrix:
  OUTPUT::Table t(name);

    t << OUTPUT::TableColumn("");
    
    for(size_t i = 0; i < matrix.size(); ++i)
    {
      t << OUTPUT::TableColumn(my_term_names[i]);
    }
      
    for(size_t i = 0; i < matrix.size(); ++i)
    {
      OUTPUT::TableRow row;
      
      row << my_term_names[i];
      
      for(size_t j = 0; j < matrix.size(); ++j)
      {
        row << matrix[i][j];
      }
      
      t << row;
    }
    
    cout << t << endl;
}

//======================================================================
//  set_phenotype_coefficients()
//======================================================================
void 
Calculator::set_phenotype_coefficients()
{
  // Set the correct phenotypic coefficients in the PointDensity vector for 
  // every valid individual.
  for(size_t i = 0; i < my_phenotype_terms.size(); ++i)
    if(my_phenotype_terms[i] != (size_t)-1)
      my_pds[my_phenotype_terms[i]].coeff_pairs[0].coefficient = my_mcc.getZScore(i);
}

//=====================================================================
//  populate_matrices()
//=====================================================================
log_double 
Calculator::populate_matrices(MAXFUN::ParameterMgr & mgr)
{
#ifdef DEBUG_CALC
  cout << "Pdms:" << endl;
#endif
  // Get the variances:
  vector<double> variances(mgr.getParamCount("Variance components"), 0.0);
  
  for(size_t i = 0; i < variances.size(); ++i)
  {
    variances[i] = mgr.getParameter("Variance components", i).getCurrentEstimate();
#ifdef DEBUG_CALC
    cout << mgr.getParameter("Variance components", i).getName() << " = " << variances[i] << endl;
#endif
  }
    
  // Set up the product-of-standard-deviations:
  log_double stdev_product(1.0);
  
  // Now populate the matrices:
  for(size_t i = 0; i < my_pds.size(); ++i)
  {
    // Get the PointDensity and PointDensityMatrix:
    PointDensity       & pd  = my_pds  [i];
    PointDensityMatrix & pdm = my_pdms [i];

#ifdef DEBUG_CALC
    cout << "Using the following pd:" << endl;
    pd.dump(my_term_names);
#endif    
    // Calculate the variance:
    double variance = 0.0;

    for(size_t j = 0; j < pd.var_idxs.size(); ++j)
    {
#ifdef DEBUG_CALC    
      cout << "Adding var#" << pd.var_idxs[j] << " which = " << variances[pd.var_idxs[j]] << endl;
#endif
      variance += variances[pd.var_idxs[j]];
    }
    
    // Make sure to factor it into the product of stdev's:
    stdev_product *= sqrt(variance);

    // Loop across all the elements and populate the matrix:
    for(size_t row = 0; row < pd.coeff_pairs.size(); ++row)
    {
      for(size_t col = row; col < pd.coeff_pairs.size(); ++col)
      {
        pdm.get_matrix()(row, col) = (pd.coeff_pairs[row].coefficient * pd.coeff_pairs[col].coefficient) / variance;
        pdm.get_matrix()(col, row) = pdm.get_matrix()(row, col);
      }
    }

#ifdef DEBUG_CALC
    cout << "variance = " << variance << endl;
    pdm.dump(my_term_names);  
#endif
  }
  
#ifdef DEBUG_CALC
  cout << "Done showing pdms." << endl;
#endif

  return stdev_product;
}

void 
Calculator::populate_r_matrix(const PointDensityMatrix & H_matrix, PointDensityMatrix & R_matrix) const
{      
  // Fetch the size of the R matrix:
  size_t r_terms = H_matrix.get_matrix().rows() - 1;
  
  // Fetch the last row/col cell from the H matrix:
  double last_row_last_col_cell = H_matrix.get_matrix()(r_terms, r_terms);

  // Populate the upper-right half with values:
  for(size_t row = 0; row < r_terms; ++row)
  {
    // Note that because it's a square matrix, we only have to do one half:
    for(size_t col = row; col < r_terms; ++col)
    {
      double row_col_cell           = H_matrix.get_matrix()(row,     col),
             row_last_col_cell      = H_matrix.get_matrix()(row,     r_terms),
             last_row_col_cell      = H_matrix.get_matrix()(r_terms, col),
             new_cell               = row_col_cell - ((row_last_col_cell * last_row_col_cell) / last_row_last_col_cell);
             
      R_matrix.get_matrix()(row, col) = new_cell;
      R_matrix.get_matrix()(col, row) = new_cell;
    }
  }

#ifdef DEBUG_CALC
  R_matrix.dump(my_term_names);
#endif
}

void
Calculator::populate_h_matrix(PointDensityMatrix & H_matrix, const IdxVector & reverse_lookup, const IdxSet & applicable_pdms) const
{
  // Clear it out:
  H_matrix.get_matrix().fill(0.0);

#ifdef DEBUG_CALC
  H_matrix.dump(my_term_names);
#endif

  // Now loop across A matrices and add it all up (just the one half, since it's a symmetric matrix).
  // Loop across A matrices to populate the reverse lookup:
  for(IdxSet::iterator idx_iter = applicable_pdms.begin(); idx_iter != applicable_pdms.end(); ++idx_iter)
  {
    // Fetch the PDM:
    const PointDensityMatrix & pdm = my_pdms[*idx_iter];

#ifdef DEBUG_CALC
    pdm.dump(my_term_names);
#endif

    for(size_t a_matrix_row = 0; a_matrix_row < pdm.get_matrix().rows(); ++a_matrix_row)
    {
      size_t absolute_row = pdm.get_idxs()[a_matrix_row],
             h_matrix_row = reverse_lookup[absolute_row];
    
      for(size_t a_matrix_col = 0; a_matrix_col < pdm.get_matrix().rows(); ++a_matrix_col)
      {
        size_t absolute_col = pdm.get_idxs()[a_matrix_col],
               h_matrix_col = reverse_lookup[absolute_col];
        
        H_matrix.get_matrix()(h_matrix_row, h_matrix_col) += pdm.get_matrix()(a_matrix_row, a_matrix_col);
      }
    }
  }

#ifdef DEBUG_CALC
  H_matrix.dump(my_term_names);
#endif
}

//=======================================================================
//  calculateLh()
//=======================================================================
double
Calculator::calculateLh(MAXFUN::ParameterMgr & mgr)
{
#ifdef DEBUG_CALC
  cout << "Beginning lh calculation..." << endl;
#endif
  // Set the coefficients for the phenotypic terms:
  set_phenotype_coefficients();

  // Set up likelihood (at the same time populating the matrices):
  log_double lh = log_double(1.0) / populate_matrices(mgr);
  
#ifdef DEBUG_CALC
  cout << "  Starting lh = " << lh << endl;
#endif

  // Now the fun part: loop through every term of integration (from the last to the first) and
  // perform the matrix magic:
  for(size_t term_idx = my_lookups.size() - 1; term_idx > 0; --term_idx)
  {
    #ifdef DEBUG_CALC
      cout << "Integrating " << my_term_names[term_idx] << endl;
    #endif
    
    // Figure out which terms are going into the H matrix:
    const IdxSet& applicable_pdms = my_lookups[term_idx];
  
    // Skip this term if there are no available A matrices for it, or the index list is empty:
    if(applicable_pdms.empty())
      continue;

    // Get the H matrix:
    PointDensityMatrix & H_matrix = my_h_matrices[term_idx].first;

    // Populate it:
    populate_h_matrix(H_matrix, my_h_matrices[term_idx].second, applicable_pdms);

    // Multiply the likelihood:
    size_t h_terms       = H_matrix.get_matrix().rows();
    double coeff_of_term = 1.0 / sqrt(H_matrix.get_matrix()(h_terms - 1, h_terms - 1)),
           pi_factor     = pow(PI_2, -double(applicable_pdms.size()) - 2.0 / 2.0);

    lh *= (coeff_of_term * pi_factor);


#ifdef DEBUG_CALC
    cout << "  Lh becomes " << lh << endl;
#endif

    // If the H-matrix is at least 2-by-2, then populate the R matrix:
    if(h_terms > 1)
      populate_r_matrix(H_matrix, my_pdms[my_r_matrices[term_idx]]);
      
  } // End term-of-integration loop

  // Now we need to deal with the remaining 1-by-1 constant matrices:
  for(IdxSet::const_iterator pdm_iter = my_constant_matrices.begin(); pdm_iter != my_constant_matrices.end(); ++pdm_iter)
  {
    // Fetch the remaining cell from the pdm, square root it, and plug it into a normal pdf evaluation:
    lh *= NUMERICS::log_normal_pdf(sqrt(my_pdms[*pdm_iter].get_matrix()(0, 0)), 1.0);

#ifdef DEBUG_CALC
    cout << "  fetching single cell " 
              << my_pdms[*pdm_iter].get_matrix()(0, 0) 
              << ", multiplying lh by " 
              << NUMERICS::log_normal_pdf(sqrt(my_pdms[*pdm_iter].get_matrix()(0, 0)), 1.0)
              << ", lh now = " 
              << lh 
              << endl;
#endif
  }

#ifdef DEBUG_CALC
  cout << "Final lh = " << lh << endl;
  if(!finite(lh.get_log()) || isnan(lh.get_log())) { exit(0); }
exit(0);
#endif

  // Return the likelihood:
  return lh.get_log();
}
 
} // End namespace ASSOC
} // End namespace SAGE

