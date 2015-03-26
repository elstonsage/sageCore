#ifndef FAM_RESID_ADJ_H
#include "pedcalc/fam_resid_adj.h"
#endif

namespace SAGE
{
namespace PED_CALC
{

/// Basic constructor
///
/// \param resid The functor for accessing residual associations
/// \param aff   The functor for accessing individual affection status
/// \param pen   The functor for accessing the penetrance for individuals
template <typename GENOTYPE, typename MPTYPE>
  template <typename RESID_FUNC,
            typename AFF_FUNC,
            typename PEN_FUNC>
inline
FraBase<GENOTYPE,MPTYPE>::FraBase
  (const RESID_FUNC&  resid,
   const AFF_FUNC&    aff,
   const PEN_FUNC&    pen)
  : my_resid_func  (resid),
    my_aff_func    (aff),
    my_pen_func    (pen)
{ }

/// Copy constructor
///
/// \param rhs The object to be copied.
template <typename GENOTYPE, typename MPTYPE>
FraBase<GENOTYPE,MPTYPE>::FraBase(const FraBase& rhs)
  : my_resid_func  (rhs.my_resid_func),
    my_aff_func    (rhs.my_aff_func),
    my_pen_func    (rhs.my_pen_func)
{ }

/// Copy Operator
///
/// \param rhs The object to be copied.    
template <typename GENOTYPE, typename MPTYPE>
FraBase<GENOTYPE,MPTYPE>& 
  FraBase<GENOTYPE,MPTYPE>::operator=(const FraBase& rhs)
{
  if(this != &rhs)
  {
    my_resid_func = rhs.my_resid_func;
    my_aff_func   = rhs.my_aff_func;
    my_pen_func   = rhs.my_pen_func;
  }
  
  return *this;
}
    
/// Destructor
///
template <typename GENOTYPE, typename MPTYPE>
FraBase<GENOTYPE,MPTYPE>::~FraBase()
{ }

/// MemberData constructor
///
/// \param pen The member's penetrance
/// \param aff The member's affection                             
template <typename GENOTYPE, typename MPTYPE>
inline
FraBase<GENOTYPE,MPTYPE>::MemberData::MemberData(double pen, bool aff)
  : penetrance(pen),
    affection (aff)
{ }

/// Default Constructor
///
/// This constructor sets the penetrance to 0.0 and the affection to false.
template <typename GENOTYPE, typename MPTYPE>
inline
FraBase<GENOTYPE,MPTYPE>::MemberData::MemberData()
  : penetrance(0),
    affection (false)
{ }
      
/// Calculates the term
///
/// \f[
///     (-1)^{y_a+y_b}
/// \f]
///
/// where \f$y_a\f$ and \f$y_b\f$ are members' affection statuses.
///
/// This is equivalent to saying 1.0 if the members' affections are the same, 
/// and -1.0 if they are different.
///
/// \param aff1 The affection of the first member
/// \param aff2 The affection of the second member
template <typename GENOTYPE, typename MPTYPE>
double
FraBase<GENOTYPE,MPTYPE>::pow_sign(bool aff1, bool aff2) const
{
  return (aff1 == aff2) ? 1.0 : -1.0;
}
  
/// Calculate the familial adjustment based upon parental genotypes and
/// child information.
///
/// \param fam         The family whose adjustment we want.
/// \param mother_geno The genotype of the mother
/// \param father_geno The genotype of the father
/// \param child_data  The penetrance and affection of the children.  These
///                    are provided by the derived class since the exact and
///                    approximate do different things here.
template <typename GENOTYPE, typename MPTYPE>
inline double
FraBase<GENOTYPE, MPTYPE>::calc_adjustment
  (const FamilyType&              fam,
   const GENOTYPE&                mother_geno,
   const GENOTYPE&                father_geno,
   const std::vector<MemberData>& child_data) const
{
  // Get parents sorted
  const MemberType* mother = fam.get_mother();
  const MemberType* father = fam.get_father();

  // If we don't have parents, we can still calculate a residual, if
  // the correlations between parents and children are the same
  if(!mother)
  {
    if(my_resid_func(resid_MS) == my_resid_func(resid_FS) &&
       my_resid_func(resid_MD) == my_resid_func(resid_FD)    )
    {
      mother = fam.parent1();
      father = fam.parent2();
    }
  }
  
  assert(mother && father && "Parents unsexed and residuals don't match!");
  
  // Get mother and father data
  MemberData mother_data(my_pen_func(*mother, mother_geno),
                         my_aff_func(*mother));

  if(isnan(mother_data.penetrance)) mother_data.penetrance = 1.0;
                         
  MemberData father_data(my_pen_func(*father, father_geno),
                         my_aff_func(*father));
  
  if(isnan(father_data.penetrance)) father_data.penetrance = 1.0;

  // Calculate adjustment
  return calc_adjustment(mother_data, father_data, child_data);
}

/// Do the actual calculation
///
/// \param mother_data The mother's data
/// \param father_data The father's data
/// \param child_data  The data for the children
template <typename GENOTYPE, typename MPTYPE>
inline double
FraBase<GENOTYPE, MPTYPE>::calc_adjustment
  (MemberData                     mother_data,
   MemberData                     father_data,
   const std::vector<MemberData>& child_data) const
{
  double ms_term = 
      calc_parent_child_term
          (my_resid_func(resid_MS), mother_data, child_data);
           
  double fs_term = 
      calc_parent_child_term
          (my_resid_func(resid_FS), father_data, child_data);
           
  double ss_term =
      calc_sibling_term
          (my_resid_func(resid_BB), child_data);
          
  double fm_term =
      calc_parent_parent_term
          (my_resid_func(resid_FM), mother_data, father_data);
          
  return 1.0 + ms_term + fs_term + ss_term + fm_term;
}         

/// Calculate the parent to child term.  This is:
///
/// \f[
///   \delta_{ps}\left(1-Pen(p,u_{p})\right)\sum_{\lambda\in C_{mf}}(-1)^{y_{p}+y_{\lambda}}\left(1-Pen(\lambda,u_{\lambda})\right)
/// \f]
/// 
/// where \f$p\f$ is the parent.
///
/// \param resid      The residual
/// \param par_data   The parental data
/// \param child_data The data for the children
template <typename GENOTYPE, typename MPTYPE>
double
FraBase<GENOTYPE, MPTYPE>::calc_parent_child_term
    (double                         resid,
     MemberData                     par_data,
     const std::vector<MemberData>& child_data) const
{
  double sum = 0.0;
  
  for(size_t i = 0; i != child_data.size(); ++i)
  {
    sum +=   pow_sign(par_data.affection, child_data[i].affection)
           * (1.0 - child_data[i].penetrance);
  }
  
  double result = resid * (1.0 - par_data.penetrance) * sum;
  
  return result;
}

/// Calculate the sibling-sibling term
///
/// \f[
///    \delta_{ss}{\displaystyle \sum_{\begin{array}{c}
///                                    \lambda<j
///                                \\  j,\lambda\in C_{mf}
///                                    \end{array}}         (-1)^{y_{\lambda}+y_{j}}\left(1-Pen(\lambda,u_{\lambda})\right)\left(1-Pen(j,u_{j})\right)}
/// \f]
///
/// \param resid      The sib-sib residual
/// \param child_data The data for the siblings
template <typename GENOTYPE, typename MPTYPE>
double
FraBase<GENOTYPE, MPTYPE>::calc_sibling_term
    (double                         resid,
     const std::vector<MemberData>& child_data) const
{
  double sum = 0.0;
  
  for(size_t i = 1; i < child_data.size(); ++i)
  {
    for(size_t j = 0; j < i; ++j)
    {
      sum +=   pow_sign(child_data[i].affection, child_data[j].affection)
             * (1.0 - child_data[i].penetrance)
             * (1.0 - child_data[j].penetrance);
    }
  }
  
  return resid * sum;
}  

/// Calculate the parent-parent term
///
/// \f[
///    \delta_{fm}(-1)^{y_{m}+y_{f}}\left(1-Pen(f,u_{f})\right)\left(1-Pen(m,u_{m})\right)
/// \f]
///
/// \param mother_data The mother's data
/// \param father_data The father's data
template <typename GENOTYPE, typename MPTYPE>
inline double 
FraBase<GENOTYPE, MPTYPE>::calc_parent_parent_term
    (double resid, MemberData mother_data,
                   MemberData father_data) const
{
  return resid * pow_sign(mother_data.affection, father_data.affection)
               * (1.0 - mother_data.penetrance) 
               * (1.0 - father_data.penetrance);
}

/// Constructor
///
/// \param resid The functor for accessing residual associations
/// \param aff   The functor for accessing individual affection status
/// \param pen   The functor for accessing the penetrance for individuals
template <typename GENOTYPE, typename MPTYPE>
  template <typename RESID_FUNC,
            typename AFF_FUNC,
            typename PEN_FUNC>
inline
ExactFamResidAdj<GENOTYPE,MPTYPE>::ExactFamResidAdj
  (const RESID_FUNC&  resid,
   const AFF_FUNC&    aff,
   const PEN_FUNC&    pen)
  : BaseType(resid, aff, pen)
{ }

/// Copy Constructor
///
/// \param rhs The object to be copied.
template <typename GENOTYPE, typename MPTYPE>
ExactFamResidAdj<GENOTYPE,MPTYPE>::ExactFamResidAdj(const SelfType& rhs)
  : BaseType(rhs)  
{ }

/// Copy Operator
///
/// \param rhs The object to be copied.
template <typename GENOTYPE, typename MPTYPE>
ExactFamResidAdj<GENOTYPE,MPTYPE>&
ExactFamResidAdj<GENOTYPE,MPTYPE>::operator=(const SelfType& rhs)
{
  if(this != &rhs)
  {
    BaseType::operator=(rhs);
  }
  
  return *this;
}

/// Calculate the adjustment for the family given the family genotypes
///
/// \param fam             The family
/// \param mother_genotype The mother's genotype
/// \param father_genotype The father's genotype
/// \param child_genotypes The genotypes of the children.
template <typename GENOTYPE, typename MPTYPE>
double 
ExactFamResidAdj<GENOTYPE,MPTYPE>::calculate_adjustment
  (const typename MPTYPE::family_type& fam,
   const Genotype&                     mother_genotype,
   const Genotype&                     father_genotype,
   const std::vector<GENOTYPE>&        child_genotypes) const
{
  std::vector<typename SelfType::MemberData> child_data(fam.offspring_count());
  
  size_t i = 0;
  for(typename MPTYPE::offspring_const_iterator child = fam.offspring_begin();
      child != fam.offspring_end(); ++child, ++i)
  {
    child_data[i].penetrance  = my_pen_func  (*child, child_genotypes[i]);
    child_data[i].affection   = my_aff_func  (*child);
    
    if(isnan(child_data[i].penetrance)) child_data[i].penetrance = 1.0;
  }
  
  return calc_adjustment(fam, mother_genotype, father_genotype, child_data);
}

/// Constructor
///
/// \param resid  The functor for accessing residual associations
/// \param aff    The functor for accessing individual affection status
/// \param susc   The functor for accessing the susceptibility for individuals
/// \param transm The functor for accessing the genotypic transmission probabilities
/// \param gbegin Iterator to the beginning of the genotypes
/// \param gend   Iterator to the end of the genotypes
template <typename GENOTYPE, typename MPTYPE, typename GTYPE_ITER>
template <typename RESID_FUNC,
          typename AFF_FUNC,
          typename SUSC_FUNC,
          typename TRANSM_FUNC>
ApproximateFamResidAdj<GENOTYPE, MPTYPE, GTYPE_ITER>::ApproximateFamResidAdj
    (const RESID_FUNC&  resid,
     const AFF_FUNC&    aff,
     const SUSC_FUNC&   susc,
     const TRANSM_FUNC& transm,
     GenotypeIter       gbegin,
     GenotypeIter       gend)
  : BaseType(resid, aff, 
    BinaryPenetranceCalculator<MemberType, GENOTYPE>(susc, aff)),
    my_transm_func(transm),
    my_susc_func(susc),
    my_gbegin(gbegin),
    my_gend(gend)
{ }
                        
/// Constructor
///
/// \param resid  The functor for accessing residual associations
/// \param aff    The functor for accessing individual affection status
/// \param susc   The functor for accessing the susceptibility for individuals
/// \param transm The functor for accessing the genotypic transmission probabilities
/// \param pen    The functor for accessing the penetrance for individuals
/// \param gbegin Iterator to the beginning of the genotypes
/// \param gend   Iterator to the end of the genotypes
template <typename GENOTYPE, typename MPTYPE, typename GTYPE_ITER>
template <typename RESID_FUNC,
          typename AFF_FUNC,
          typename SUSC_FUNC,
          typename TRANSM_FUNC,
          typename PEN_FUNC>
ApproximateFamResidAdj<GENOTYPE, MPTYPE, GTYPE_ITER>::ApproximateFamResidAdj
    (const RESID_FUNC&  resid,
     const AFF_FUNC&    aff,
     const SUSC_FUNC&   susc,
     const TRANSM_FUNC& transm,
     GenotypeIter       gbegin,
     GenotypeIter       gend,
     const PEN_FUNC&    pen)
  : BaseType(resid, aff, susc, pen),
    my_transm_func(transm),
    my_susc_func(susc),
    my_gbegin(gbegin),
    my_gend(gend)
{ }

/// Copy Constructor
///
/// \param rhs The object to be copied.                       
template <typename GENOTYPE, typename MPTYPE, typename GTYPE_ITER>
ApproximateFamResidAdj<GENOTYPE, MPTYPE, GTYPE_ITER>::ApproximateFamResidAdj
  (const SelfType& rhs)
  : BaseType(rhs),
    my_transm_func(rhs.my_transm_func),
    my_gbegin(rhs.my_gbegin),
    my_gend(rhs.my_gend)
{ }

/// Copy Operator
///
/// \param rhs The object to be copied.
template <typename GENOTYPE, typename MPTYPE, typename GTYPE_ITER>
ApproximateFamResidAdj<GENOTYPE, MPTYPE, GTYPE_ITER>& 
  ApproximateFamResidAdj<GENOTYPE, MPTYPE, GTYPE_ITER>::operator=
    (const SelfType& rhs)
{
  if(this != &rhs)
  {
    BaseType::operator=(rhs);
    
    my_transm_func = rhs.my_transm_func;
    my_gbegin      = rhs.my_gbegin;
    my_gend        = rhs.my_gend;
  }
  
  return *this;
}
    
/// Calculate the adjustment
///
/// \param fam             The family whose adjustment we want.
/// \param mother_genotype The mother's genotype
/// \param father_genotype The father's genotype
template <typename GENOTYPE, typename MPTYPE, typename GTYPE_ITER>
double
ApproximateFamResidAdj<GENOTYPE, MPTYPE, GTYPE_ITER>::calculate_adjustment
  (const FamilyType& fam,
   const Genotype&   mother_genotype,
   const Genotype&   father_genotype) const
{
  std::vector<typename SelfType::MemberData> child_data(fam.offspring_count());
  
  size_t i = 0;
  for(typename MPTYPE::offspring_const_iterator child = fam.offspring_begin();
      child != fam.offspring_end(); ++child, ++i)
  {
    child_data[i].penetrance  = calc_est_penetrance(*child, mother_genotype, father_genotype);
    child_data[i].affection   = my_aff_func(*child);
    
    if(isnan(child_data[i].penetrance)) child_data[i].penetrance = 1.0;
  }
  
  return calc_adjustment(fam, mother_genotype, father_genotype, child_data);
}
                         
/// Calculate a member's estimated penetrance:
///
/// \f[
///     \frac{e^{\hat{\theta}_{\lambda}y_{\lambda}}}{1+e^{\hat{\theta}_{\lambda}}}
/// \f]
///
/// where \f$\hat{\theta}_{\lambda}\f$ is the estimated susceptibility for the
/// member
///
/// \param mem   The member we're calculating
/// \param mgeno The genotype of the mother (needed for estimated susceptibilty calc)
/// \param fgeno The genotype of the father (needed for estimated susceptibilty calc)
template <typename GENOTYPE, typename MPTYPE, typename GTYPE_ITER>
inline double
ApproximateFamResidAdj<GENOTYPE, MPTYPE, GTYPE_ITER>::calc_est_penetrance
  (const MemberType& mem,
   const Genotype&   mgeno,
   const Genotype&   fgeno) const
{
  double mean = calc_est_susceptibility(mem,mgeno,fgeno);
  
  bool affection = my_aff_func(mem);
  
  return exp(mean * affection) / (1.0 + exp(mean));
}
    
/// Calculate the member's estimated susceptibility
///
/// \f[
///     \hat{\theta}_{\lambda}=\frac{\sum_{u_{\lambda}}P\left(u_{\lambda}\mid u_{m},u_{f}\right)\frac{e^{\theta_{u}(\lambda)y_{\lambda}}}{1+e^{\theta_{u}(\lambda)}}\theta_{u}(\lambda)}
///                                 {\sum_{u_{\lambda}}P\left(u_{\lambda}\mid u_{m},u_{f}\right)\frac{e^{\theta_{u}(\lambda)y_{\lambda}}}{1+e^{\theta_{u}(\lambda)}}}
/// \f]
///
/// where \f$\theta_{u}(\lambda)\f$ is the susceptibility of the member
/// \f$\lambda\f$ for genotype \f$u\f$.
///
/// \param mem   The member we're calculating
/// \param mgeno The genotype of the mother (needed for estimated susceptibilty calc)
/// \param fgeno The genotype of the father (needed for estimated susceptibilty calc)
template <typename GENOTYPE, typename MPTYPE, typename GTYPE_ITER>
inline double
ApproximateFamResidAdj<GENOTYPE, MPTYPE, GTYPE_ITER>::calc_est_susceptibility
  (const MemberType& mem,
   const Genotype&   mgeno,
   const Genotype&   fgeno) const
{
  double numerator_sum    = 0.0;
  double denominator_sum  = 0.0;

  for(GenotypeIter gtype = my_gbegin; gtype != my_gend; ++gtype)
  { 
    double transmission   = my_transm_func (mem, mgeno, fgeno, (Genotype) *gtype);
    double penetrance     = my_pen_func    (mem, (Genotype) *gtype);
    double susceptibility = my_susc_func   (mem, (Genotype) *gtype);
    
    double denom_term = transmission * penetrance;
    double num_term   = transmission * penetrance * susceptibility;

    numerator_sum   += num_term;
    denominator_sum += denom_term;
  }

  return  numerator_sum / denominator_sum;
}
                                    
} // End Namespace PED_CALC
} // End Namespace SAGE
                                
