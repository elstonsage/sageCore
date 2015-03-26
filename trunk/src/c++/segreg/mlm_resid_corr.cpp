//===========================================================================
//  File:       mlm_resid_corr.cpp
//
//  Author:     Kai He                          
//
//  History:
//                                                 
//  Copyright (c) 2001 R. C. Elston
//  All rights reserved
//===========================================================================

#include "segreg/mlm_resid_corr.h"

#define FM resid_model::fm
#define FS resid_model::fs
#define MS resid_model::ms
#define FD resid_model::fd
#define MD resid_model::md
#define BB resid_model::bb
#define BS resid_model::bs
#define SS resid_model::ss

using namespace SAGE;
using namespace SEGREG;
using namespace std;

//===========================================================================
//
//  constructor - build data structures and run calculations
//
//===========================================================================
MlmResidCorrelationCalculator::MlmResidCorrelationCalculator()
{ }

MlmResidCorrelationCalculator::MlmResidCorrelationCalculator
    (const binary_member_calculator& b,
     const PedigreeDataSet::SubpedigreeCursor& spl, 
     const model& md)
{
   calculate(b, spl, md);
}

MlmResidCorrelationCalculator::MlmResidCorrelationCalculator
    (const MlmResidCorrelationCalculator& r)
  : my_means     (r.my_means),
    my_theta_i   (r.my_theta_i),
    my_theta_u_i (r.my_theta_u_i),
    my_tau_means (r.my_tau_means),
    my_corr_ij   (r.my_corr_ij),
    my_corr_ij_u (r.my_corr_ij_u)
{ }

MlmResidCorrelationCalculator& MlmResidCorrelationCalculator::operator=
(const MlmResidCorrelationCalculator& r)
{
  if(this != &r)
  {
    my_means     = r.my_means;
    my_theta_i   = r.my_theta_i;
    my_theta_u_i = r.my_theta_u_i;
    my_tau_means = r.my_tau_means;
    my_corr_ij   = r.my_corr_ij;
    my_corr_ij_u = r.my_corr_ij_u;
  }
  return *this;
}

void
MlmResidCorrelationCalculator::classify_members
   (const binary_member_calculator& b,
    const PedigreeDataSet::SubpedigreeCursor& spl)
{
  // First, clear our lists
  
  my_mother_list.clear();
  my_father_list.clear();
  my_parent_list.clear();
  my_child_list .clear();

  // Iterate through the subpedigrees.
  
  for(PedigreeDataSet::SubpedigreeIterator subped = spl.first;
          subped != spl.second; ++subped)
  {
    // Iterate through the subpedigree's families.
    for(FPED::FamilyConstIterator fam  = subped->family_begin(); 
                                  fam != subped->family_end(); 
                                ++fam)
    {
      // Classify the parents.  We assume that the parents are both sexed.
      // If this is not the case, then the parents are arbitrarily assigned
      // p1 as the mother and p2 as the father.  Note that the parents are both
      // sexed or neither are, or we wouldn't get this far.
      FPED::MemberConstPointer p1 = fam->parent1();
      FPED::MemberConstPointer p2 = fam->parent2();
      
      if(p1->is_male() && p2->is_female())	
      {
        if(b.is_member_valid(*p2)) my_mother_list.push_back(p2);
        if(b.is_member_valid(*p1)) my_father_list.push_back(p1);
      }
      else
      {
        if(b.is_member_valid(*p1)) my_mother_list.push_back(p1);
        if(b.is_member_valid(*p2)) my_father_list.push_back(p2);
      }
      
      // Push them back as parents too
      my_parent_list.push_back(p1);
      my_parent_list.push_back(p2);
      
      // Classify the children of this family.
      for(FPED::OffspringConstIterator sib  = fam->offspring_begin();
                                       sib != fam->offspring_end();
                                     ++sib)
      {
        if(b.is_member_valid(*sib)) my_child_list.push_back(&*sib);
      }
    } //end of fam iteration
  } //end of subped iteration
}

//===========================================================================
//
//  calculate(subplist,md) - calculate all 
//
//===========================================================================
void 
MlmResidCorrelationCalculator::calculate
    (const binary_member_calculator& b, 
     const PedigreeDataSet::SubpedigreeCursor& spl, 
     const model& md)
{
  // First, clear our information

  my_means.assign(0.0);
  
  double& my_fa_mean_sum = my_means[0];
  double& my_mo_mean_sum = my_means[1];
  double& my_pa_mean_sum = my_means[2];
  double& my_ch_mean_sum = my_means[3];
  
  my_tau_means[0].assign(0.0);
  my_tau_means[1].assign(0.0);
  my_tau_means[2].assign(0.0);
  my_tau_means[3].assign(0.0);

  // Classify the members of our subpedigrees into the father, mother, parent and child
  // subsets
  
  classify_members(b, spl);

  // Get our covariate info from the susceptibility sub model

  const vector<CovariateSubmodel::covariate>& cov = md.susc_cov_sub_model.covariates();

  // Create a vector to store information for each covariate

  vector<class_cov> cov_info(cov.size());

  // Iterate through our covariates:
  // 1. Gather the data sample for the covariate for the mother, father, parent and offspring
  // 2. Calculate the sum of covariate means * coefficients
  // 3. Calculate the interaction sums
  
  for(size_t i = 0; i < cov.size(); ++i)
  {
    // Get a reference to the current covariate
    
    class_cov& current = cov_info[i];

    // 1. GATHER SAMPLES
    // =================
    
    // Gather the mother sample
    
    current.mo_sinfo.clear();

    for(member_list::const_iterator j  = my_mother_list.begin();
                                    j != my_mother_list.end();
                                  ++j)
    {
      double cov_value = (*j)->pedigree()->info().trait((*j)->index(),cov[i].trait_index);
      
      current.mo_sinfo.add(cov_value);
    }

    // Gather the father sample
    
    current.fa_sinfo.clear();

    for(member_list::const_iterator j  = my_father_list.begin();
                                    j != my_father_list.end();
                                  ++j)
    {
      double cov_value = (*j)->pedigree()->info().trait((*j)->index(),cov[i].trait_index);
      
      current.fa_sinfo.add(cov_value);
    }

    // Gather the parent sample
    
    current.pa_sinfo.clear();

    for(member_list::const_iterator j  = my_parent_list.begin();
                                    j != my_parent_list.end();
                                  ++j)
    {
      double cov_value = (*j)->pedigree()->info().trait((*j)->index(),cov[i].trait_index);
      
      current.pa_sinfo.add(cov_value);
    }

    // Gather the child sample
    
    current.ch_sinfo.clear();

    for(member_list::const_iterator j  = my_child_list.begin();
                                    j != my_child_list.end();
                                  ++j)
    {
      double cov_value = (*j)->pedigree()->info().trait((*j)->index(),cov[i].trait_index);
      
      current.ch_sinfo.add(cov_value);
    }

    // 2. Calculate the sum of covariate means
    
    my_mo_mean_sum += cov[i].coefficient * current.mo_sinfo.mean();
    my_fa_mean_sum += cov[i].coefficient * current.fa_sinfo.mean();
    my_pa_mean_sum += cov[i].coefficient * current.pa_sinfo.mean();
    my_ch_mean_sum += cov[i].coefficient * current.ch_sinfo.mean();
    
    // 3. Calculate the interaction terms for each genotype, if applicable
    
    if(cov[i].has_interaction == true)
    {
      for(genotype_index gi = index_AA; gi != index_INVALID; ++gi)    
      {
        double Beta_u = md.susc_sub_model.parameter(gi);
        double tau_u  = cov[i].i_taus[gi];

        my_tau_means[0][gi] += Beta_u * tau_u * current.fa_sinfo.mean();
        my_tau_means[1][gi] += Beta_u * tau_u * current.mo_sinfo.mean();
        my_tau_means[2][gi] += Beta_u * tau_u * current.pa_sinfo.mean();
        my_tau_means[3][gi] += Beta_u * tau_u * current.ch_sinfo.mean();
      }
    }
  }

  // calculate theta_u_i given genotype, class status, and Beta_u

  calculate_theta_u_i(md);

  // calculate theta_i given q_AA, theta_u_i

  calculate_theta_i(md);

  // calculate corr_ij given Delta_corr and theta_i

  calculate_corrs(md);

  // testing outputs of data member
  ///dump_data(md);
}


//===========================================================================
//
//  calculate_corr()
//
//===========================================================================
double
MlmResidCorrelationCalculator::calculate_corr
    (double delta_ij,
     double theta_i,
     double theta_j) const
{
  double numerator   = delta_ij * exp( (theta_i + theta_j) / 2.0 );
  double denominator = (1.0 + exp(theta_i)) * (1.0 + exp(theta_j));
  double corr        = numerator / denominator;
  
  return corr;
}

//===========================================================================
//
//  Theta_u_i() - calculate theta_u_i with status i and genotype and fill
//  up the data structure
//
//===========================================================================
void 
MlmResidCorrelationCalculator::calculate_theta_u_i(const model& md)
{
  for(size_t i = 0; i < my_means.size(); ++i)
  {
    for(genotype_index gi = index_AA; gi != index_INVALID; ++gi)
    {
        double Beta_u = md.susc_sub_model.parameter(gi);
   
        my_theta_u_i[i][gi] = Beta_u + 
                              my_means[i] +
                              my_tau_means[i][gi];
    }
  }
}
//===========================================================================
//
//  Theta_i() - calculate theta_i with status i and fill up data structure
//
//===========================================================================
void 
MlmResidCorrelationCalculator::calculate_theta_i(const model& md)
{
  double A = md.freq_sub_model.freq_A();

  double aa =  A * A;
  double ab =  2.0 * A  * (1.0 - A);
  double bb = (1.0 - A) * (1.0 - A);

  for(size_t i = 0; i < my_means.size(); ++i)
  {
      // i = fa, mo, pa, ss
      my_theta_i[i] = aa * my_theta_u_i[i][index_AA] + 
                      ab * my_theta_u_i[i][index_AB] + 
                      bb * my_theta_u_i[i][index_BB]; 
  }
}

void
MlmResidCorrelationCalculator::calculate_ij_corrs
    (resid_model::corr ctype,
     double delta,
     size_t index1,
     size_t index2)
{
  double theta_i = my_theta_i[index1];
  double theta_j = my_theta_i[index2];

  my_corr_ij[ctype] = calculate_corr(delta, theta_i, theta_j);
  
  for(genotype_index gi = index_AA; gi != index_INVALID; ++gi)
  {
    double theta_u_i = my_theta_u_i[index1][gi];
    double theta_u_j = my_theta_u_i[index2][gi];

    my_corr_ij_u[ctype][gi]   = calculate_corr(delta, theta_u_i, theta_u_j);
  }
}

void
MlmResidCorrelationCalculator::calculate_corrs(const model& md)
{
  double Delta   = 0.0;
  
  size_t m = 2;
  size_t f = 2;

  // If we have arbitrary (f != m) correlations, use the sex-specific means
  if(md.resid_sub_model.option() == residual_correlation_sub_model::arb)
  {
    m = 1;
    f = 0;
  }
  
  Delta   = md.resid_sub_model.father_mother_correlation();

  calculate_ij_corrs(resid_model::fm, Delta, f, m);

  Delta   = md.resid_sub_model.father_son_correlation();

  calculate_ij_corrs(resid_model::fs, Delta, f, 3);
  calculate_ij_corrs(resid_model::fd, Delta, f, 3);

  Delta   = md.resid_sub_model.mother_son_correlation();

  calculate_ij_corrs(resid_model::ms, Delta, m, 3);
  calculate_ij_corrs(resid_model::md, Delta, m, 3);

  Delta   = md.resid_sub_model.sister_sister_correlation();

  calculate_ij_corrs(resid_model::bb, Delta, 3, 3);
  calculate_ij_corrs(resid_model::bs, Delta, 3, 3);
  calculate_ij_corrs(resid_model::ss, Delta, 3, 3);
}

/*
vector<double> 
MlmResidCorrelationCalculator::get_theta_i()
{
   return my_theta_i;
}

vector<vector<double> > 
MlmResidCorrelationCalculator::get_theta_u_i()
{
   return my_theta_u_i;
}
*/
double 
MlmResidCorrelationCalculator::get_corr_ij(resid_model::corr i)
{
   return my_corr_ij[i];
}

double 
MlmResidCorrelationCalculator::get_corr_ij_u(resid_model::corr i, size_t u)
{
   return my_corr_ij_u[i][u];
}
//===========================================================================
//
//  dump_data() - for testing and debugging
//
//===========================================================================
void 
MlmResidCorrelationCalculator::dump_data(const model& md)
{
  cout<<"Beta_susc.parameter(index_AA)="<<md.susc_sub_model.parameter(index_AA)<<endl;
  cout<<"Beta_susc.parameter(index_AB)="<<md.susc_sub_model.parameter(index_AB)<<endl;
  cout<<"Beta_susc.parameter(index_BB)="<<md.susc_sub_model.parameter(index_BB)<<endl;
  cout<<"qA.freq_A()="<<md.freq_sub_model.freq_A()<<endl;

  cout<<"Delta(fm)="<<md.resid_sub_model.father_mother_correlation()<<endl;
  cout<<"Delta(fs)="<<md.resid_sub_model.father_son_correlation()<<endl;
  cout<<"Delta(ms)="<<md.resid_sub_model.mother_son_correlation()<<endl;
  cout<<"Delta(ss)="<<md.resid_sub_model.sister_sister_correlation()<<endl;
  cout<<"Delta(fd)="<<md.resid_sub_model.father_daughter_correlation()<<endl;
  cout<<"Delta(md)="<<md.resid_sub_model.mother_daughter_correlation()<<endl;
  cout<<"Delta(bb)="<<md.resid_sub_model.brother_brother_correlation()<<endl;
  cout<<"Delta(bs)="<<md.resid_sub_model.brother_sister_correlation()<<endl;

  cout<<endl;  
  cout<<"my_means.size="<<my_means.size()<<endl;
  cout<<"my_theta_i.size="<<my_theta_i.size()<<endl;
  cout<<"my_theta_u_i.size="<<my_theta_u_i.size()<<endl;

  cout<<"my_means:"<<endl;
  for(size_t i = 0; i < my_means.size(); ++i)
      cout<<"my_means["<<i<<"]="
          <<my_means[i]<<endl;

  cout<<endl;
  cout<<"theta_i:"<<endl;
  for(size_t i = 0; i < my_theta_i.size(); ++i)
      cout<<"my_theta_i["<<i<<"]="
          <<my_theta_i[i]<<endl;

  cout<<endl;
  cout<<"theta_u_i:"<<endl;
  for(size_t i = 0; i < my_theta_u_i.size(); ++i)
  {
      for(size_t j = 0; j < 3; ++j)
          cout<<"my_theta_u_i["<<i
              <<"]["<<j<<"]="
              <<my_theta_u_i[i][j]<<endl;;
  }      

  cout<<endl;
  cout<<"corr_ij overall:"<<endl;
  for(int i = 0; i < NUM_CORR; ++i)
      cout<<"corr_ij["<<i<<"]="
          <<my_corr_ij[i]<<endl;

  cout<<endl;
  cout<<"corr_ij_u genotype:"<<endl;
  for(int i = 0; i < NUM_CORR; ++i)
  {
     for(int j = 0; j < NUM_GENOTYPES; ++j)
         cout<<"corr_ij_u["<<i<<"]["<<j<<"]="
             <<my_corr_ij_u[i][j]<<endl;
  }
}

