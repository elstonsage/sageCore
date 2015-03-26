#ifndef FAM_RESID_ADJ_H
#define FAM_RESID_ADJ_H

#include "pedcalc/binary_penetrance_calculator.h"
#include "functors/const_ref.h"
#include "boost/function.hpp"
#include "boost/bind.hpp"
#include "boost/iterator/counting_iterator.hpp"
#include <vector>
#include <iomanip>
#include <iostream>

namespace SAGE
{
namespace PED_CALC
{

/// \brief Enumerations of the residual associations used by the familial residual association objects
enum FraResidType
{
  resid_FM, ///< Father-mother   residual
  resid_FS, ///< Father-son      residual
  resid_MS, ///< Mother-son      residual
  resid_FD, ///< Father-daughter residual
  resid_MD, ///< Mother-daughter residual
  resid_BB, ///< Brother-brother residual
  resid_BS, ///< Brother-sister  residual
  resid_SS  ///< Sister-sister   residual
};

/// \internal
///
/// \brief Base class for calculating familial residual adjustments due to binary traits.
///
/// This class calculates the adjustment applied to a nuclear family's likelihood
/// due to residual associations of a binary trait modeled on some
/// underlying model.
///
/// The equation it calculates is:
///
/// \f[
///    \begin{array}{lll}
///    \rho=1 & + & \delta_{fs}\left(1-Pen(f,u_{f})\right){\displaystyle \sum_{\lambda\in C_{mf}}(-1)^{y_{f}+y_{\lambda}}\left(1-Pen(\lambda,u_{\lambda})\right)}
///    \\     & + & \delta_{ms}\left(1-Pen(m,u_{m})\right){\displaystyle \sum_{\lambda\in C_{mf}}(-1)^{y_{m}+y_{\lambda}}\left(1-Pen(\lambda,u_{\lambda})\right)}
///    \\     & + & \delta_{ss}{\displaystyle \sum_{\begin{array}{c}
///                                                 \lambda<j 
///    \\                                           j,\lambda\in C_{mf}
///                                                 \end{array}}         (-1)^{y_{\lambda}+y_{j}}\left(1-Pen(\lambda,u_{\lambda})\right)\left(1-Pen(j,u_{j})\right)}
///    \\     & + & \delta_{fm}(-1)^{y_{m}+y_{f}} 
///                 \left(1-Pen(f,u_{f})\right)
///                 \left(1-Pen(m,u_{m})\right)
///    \end{array}
/// \f]
///
/// where:
///
///  - \f$\delta_{xy}\f$  is the residual association between two relatives (f = father,
///                       m = mother, s = sibling).
///  - \f$Pen(x,u_{x})\f$ is the penetrance value for member \f$x\f$ given 
///                       genotype \f$u_{x}\f$. 
///  - \f$y_x\f$          is the affection status of the member.
///
/// Note that there is no explicit calculation of the penetrance using susceptibility
/// This is because firstly, the penetrance is often externally calculated, and
/// there is no reason to repeat work, and second because the suceptibilities
/// of the children may be approximate rather than exact.  By writing the
/// equation in this generic form, both the ExactFamResidAdj and the ApproximateFamResidAdj
/// can use this object as a base type.
///
/// \bold USEAGE
///
/// The class works by way of functors.  This is to make the class as generic as
/// possible.  Rather than basing the class on some originating types where
/// the data is stored, it retrieves values by use of three functions that
/// are passed to it at creation.  The functions have the following prototypes:
///
///  - Residual: double (FraResidType)
///  - Penetrance: double (const MPTYPE::member_type&, const GENOTYPE&)
///  - Affection: bool (const MPTYPE::member_type&)
///
/// Any function or functor which is convertible to these prototypes is
/// acceptable (see boost::function for details)
///
/// Because this class is a base class and should never be used directly,
/// the entire class interface is protected or private.  This makes it difficult
/// to use this class incorrectly.
template <typename GENOTYPE, typename MPTYPE>
class FraBase
{
  protected:
  
    typedef typename MPTYPE::family_type FamilyType;
    typedef typename MPTYPE::member_type MemberType;
    
  /// \name Object Management
  ///
  /// All object management is in the protected namespace.  This
  /// makes it difficult to use the object incorrectly.
  //@}
    template <typename RESID_FUNC,
              typename AFF_FUNC,
              typename PEN_FUNC>
    FraBase(const RESID_FUNC&  resid,
            const AFF_FUNC&    aff,
            const PEN_FUNC&    pen);
                            
    FraBase(const FraBase&);
    
    FraBase& operator=(const FraBase&);
    
    ~FraBase();
    
  //@}
        
                         
    /// \brief Stores a member's penetrance and affection
    ///
    /// This struct is a simple mechanism for storing penetrance and affection
    /// data for computation.
    struct MemberData
    {
      MemberData(double, bool);
      MemberData();
      
      double penetrance; ///< The penetrance of the member
      bool   affection;  ///< The affection status of the member
    };
    
    double pow_sign(bool, bool) const;
  
    double calc_adjustment         (const FamilyType&              fam,
                                    const GENOTYPE&                mother_data,
                                    const GENOTYPE&                father_data,
                                    const std::vector<MemberData>& child_data) const;
                                    
    double calc_adjustment         (MemberData                     mother_data,
                                    MemberData                     father_data,
                                    const std::vector<MemberData>& child_data) const;
                                
    double calc_parent_child_term  (double                         resid,
                                    MemberData                     par_data,
                                    const std::vector<MemberData>& child_data) const;
                                  
    double calc_sibling_term       (double                         resid,
                                    const std::vector<MemberData>& child_data) const;
    double calc_parent_parent_term (double                         resid,
                                    MemberData                     mother_data,
                                    MemberData                     father_data) const;

    typedef boost::function <double (FraResidType)>      ResidFuncType;
    typedef boost::function <bool   (const MemberType&)> AffFuncType;
    typedef boost::function <double (const MemberType&,
                                     const GENOTYPE&)>   PenFuncType;
    
    ResidFuncType  my_resid_func; ///< Residual Functor
    AffFuncType    my_aff_func;   ///< Affection Functor
    PenFuncType    my_pen_func;   ///< Penetrance Functor
};

/// \internal
///
/// \brief Class for calculating exact familial residual adjustments due to binary traits.
///
/// This class calculates the adjustment applied to a nuclear family's likelihood
/// due to residual associations of a binary trait modeled on some
/// underlying model.
///
/// The equation it calculates is:
///
/// \f[
///    \begin{array}{lll}
///    \rho=1 & + & \delta_{fs}\left(1-Pen(f,u_{f})\right){\displaystyle \sum_{\lambda\in C_{mf}}(-1)^{y_{f}+y_{\lambda}}\left(1-Pen(\lambda,u_{\lambda})\right)}
///    \\     & + & \delta_{ms}\left(1-Pen(m,u_{m})\right){\displaystyle \sum_{\lambda\in C_{mf}}(-1)^{y_{m}+y_{\lambda}}\left(1-Pen(\lambda,u_{\lambda})\right)}
///    \\     & + & \delta_{ss}{\displaystyle \sum_{\begin{array}{c}
///                                                 \lambda<j 
///    \\                                           j,\lambda\in C_{mf}
///                                                 \end{array}}         (-1)^{y_{\lambda}+y_{j}}\left(1-Pen(\lambda,u_{\lambda})\right)\left(1-Pen(j,u_{j})\right)}
///    \\     & + & \delta_{fm}(-1)^{y_{m}+y_{f}} 
///                 \left(1-Pen(f,u_{f})\right)
///                 \left(1-Pen(m,u_{m})\right)
///    \end{array}
/// \f]
///
/// where:
///
///  - \f$\delta_{xy}\f$  is the residual association between two relatives (f = father,
///                       m = mother, s = sibling).
///  - \f$Pen(x,u_{x})\f$ is the penetrance value for member \f$x\f$ given 
///                       genotype \f$u_{x}\f$. 
///  - \f$y_x\f$          is the affection status of the member.
///
/// Note that there is no explicit calculation of the penetrance using susceptibility
/// This is because the penetrance is often externally calculated, and
/// there is no reason to repeat work.
///
/// \bold USEAGE
///
/// The class works by way of functors.  This is to make the class as generic as
/// possible.  Rather than basing the class on some originating types where
/// the data is stored, it retrieves values by use of three functions that
/// are passed to it at creation.  The functions have the following prototypes:
///
///  - Residual: double (FraResidType)
///  - Penetrance: double (const MPTYPE::member_type&, const GENOTYPE&)
///  - Affection: bool (const MPTYPE::member_type&)
///
/// Any function or functor which is convertible to these prototypes is
/// acceptable (see boost::function for details)
///
/// Because this class is a base class and should never be used directly,
/// the entire class interface is protected or private.  This makes it difficult
/// to use this class incorrectly.
///
/// \internal
///
/// This class uses its base, FraBase<GENOTYPE, MPTYPE> to do the calculations.
/// The class itself is simply a mechanism for providing a more user-oriented
/// interface.  It looks up any values it needs and passes them on to the base class.
template <typename GENOTYPE,  typename MPTYPE>
class ExactFamResidAdj : protected FraBase<GENOTYPE, MPTYPE>
{
  public:
  
    typedef          MPTYPE                    MpedType; 
    typedef typename MPTYPE::family_type       FamilyType;
    typedef typename MPTYPE::member_type       MemberType;

    typedef FraBase<GENOTYPE, MPTYPE>          BaseType;
    typedef ExactFamResidAdj<GENOTYPE, MPTYPE> SelfType;

    typedef GENOTYPE Genotype;

  /// \name Object Management
  //@{
    template <typename RESID_FUNC,
              typename AFF_FUNC,
              typename PEN_FUNC>
    ExactFamResidAdj(const RESID_FUNC&  resid,
                     const AFF_FUNC&    aff,
                     const PEN_FUNC&    pen);
                            
    ExactFamResidAdj(const SelfType&);
    
    ExactFamResidAdj& operator=(const SelfType&);
  //@}
  
    double calculate_adjustment
              (const typename MPTYPE::family_type& fam,
               const Genotype&                     mother_genotype,
               const Genotype&                     father_genotype,
               const std::vector<Genotype>&        child_genotypes) const;
};

/// \internal
///
/// \brief Class for calculating exact familial residual adjustments due to binary traits.
///
/// This class calculates the adjustment applied to a nuclear family's likelihood
/// due to residual associations of a binary trait modeled on some
/// underlying model.
///
/// The equation it calculates is:
///
/// \f[
///    \begin{array}{lll}
///    \rho=1 & + & \delta_{fs}\left(1-Pen(f,u_{f})\right){\displaystyle \sum_{\lambda\in C_{mf}}(-1)^{y_{f}+y_{\lambda}}\left(1-Pen(\lambda)\right)}
///    \\     & + & \delta_{ms}\left(1-Pen(m,u_{m})\right){\displaystyle \sum_{\lambda\in C_{mf}}(-1)^{y_{m}+y_{\lambda}}\left(1-Pen(\lambda)\right)}
///    \\     & + & \delta_{ss}{\displaystyle \sum_{\begin{array}{c}
///                                                 \lambda<j 
///    \\                                           j,\lambda\in C_{mf}
///                                                 \end{array}}         (-1)^{y_{\lambda}+y_{j}}\left(1-Pen(\lambda)\right)\left(1-Pen(j)\right)}
///    \\     & + & \delta_{fm}(-1)^{y_{m}+y_{f}} 
///                 \left(1-Pen(f,u_{f})\right)
///                 \left(1-Pen(m,u_{m})\right)
///    \end{array}
/// \f]
///
/// where:
///
///  - \f$\delta_{xy}\f$  is the residual association between two relatives (f = father,
///                       m = mother, s = sibling).
///  - \f$Pen(x,u_{x})\f$ is the penetrance value for member \f$x\f$ given 
///                       genotype \f$u_{x}\f$. 
///  - \f$Pen(x)\f$       is a penetrance based upon the weighted sum of all
///                       genotypes for an individua \f$x\f$ (see below).
///  - \f$y_x\f$          is the affection status of the member.
///
/// In the approximate method, the child genotypes are not given explicitly.
/// Instead, the child suscptibility is calculated as:
///
/// \f[
///     \hat{\theta}_{\lambda}=\frac{\sum_{u_{\lambda}}P\left(u_{\lambda}\mid u_{m},u_{f}\right)\frac{e^{\theta_{u}(\lambda)y_{\lambda}}}{1+e^{\theta_{u}(\lambda)}}\theta_{u}(\lambda)}
///                                 {\sum_{u_{\lambda}}P\left(u_{\lambda}\mid u_{m},u_{f}\right)\frac{e^{\theta_{u}(\lambda)y_{\lambda}}}{1+e^{\theta_{u}(\lambda)}}}
/// \f]
///
/// where:
///  - \f$\theta_{u}(\lambda)\f$ is the susceptibility of the child \f$\lambda\f$ for genotype \f$u\f$.
///
/// From this it uses the standard form to calculate child penetrance (see the 
/// BinaryPenetranceCalculator).
///
/// \bold USEAGE
///
/// The class works by way of functors.  This is to make the class as generic as
/// possible.  Rather than basing the class on some originating types where
/// the data is stored, it retrieves values by use of three functions that
/// are passed to it at creation.  The functions have the following prototypes:
///
///  - Residual: double (FraResidType)
///  - Susceptiblity: double (const MPTYPE::member_type&, const GENOTYPE&)
///  - Affection: bool (const MPTYPE::member_type&)
///  - Transmission: double (const MPTYPE::member_type&, const GENOTYPE& child, const GENOTYPE& mother, const GENOTYPE& father)
///  - Penetrance: double (const MPTYPE::member_type&, const GENOTYPE&) (optional)
///
/// Any function or functor which is convertible to these prototypes is
/// acceptable (see boost::function for details)
///
/// Note thate while the Penetrance is optional, it is best to provide it if
/// it has already been calculated elsewhere in the program.
///
/// Additional to the functors, to do the approximation, the class needs genotype iterators
/// to allow it to calculate the approximate susceptibility for the children.
///
/// Because this class is a base class and should never be used directly,
/// the entire class interface is protected or private.  This makes it difficult
/// to use this class incorrectly.
///
/// \internal
///
/// This class uses its base, FraBase<GENOTYPE, MPTYPE> to do the calculations.
/// The class itself is simply a mechanism for providing a more user-oriented
/// interface and doing the approximation calculations which it uses.  The
/// actual calculation is performed by the base class.
template <typename GENOTYPE,
          typename MPTYPE,
          typename GTYPE_ITER = boost::counting_iterator<GENOTYPE> >
class ApproximateFamResidAdj : protected FraBase<GENOTYPE, MPTYPE>
{
  public:
  
    typedef          MPTYPE                                      MpedType; 
    typedef typename MPTYPE::family_type                         FamilyType;
    typedef typename MPTYPE::member_type                         MemberType;

    typedef FraBase<GENOTYPE, MPTYPE>                            BaseType;
    typedef ApproximateFamResidAdj<GENOTYPE, MPTYPE, GTYPE_ITER> SelfType;

    typedef GENOTYPE                                             Genotype;
    typedef GTYPE_ITER                                           GenotypeIter;

  /// \name Object Management
  //@{
    template <typename RESID_FUNC,
              typename AFF_FUNC,
              typename SUSC_FUNC,
              typename TRANSM_FUNC>
    ApproximateFamResidAdj(const RESID_FUNC&   resid,
                           const AFF_FUNC&     aff,
                           const SUSC_FUNC&    susc,
                           const TRANSM_FUNC&  transm,
                           GenotypeIter        gbegin,
                           GenotypeIter        gend);
                            
    template <typename RESID_FUNC,
              typename AFF_FUNC,
              typename SUSC_FUNC,
              typename TRANSM_FUNC,
              typename PEN_FUNC>
    ApproximateFamResidAdj(const RESID_FUNC&   resid,
                           const AFF_FUNC&     aff,
                           const SUSC_FUNC&    susc,
                           const TRANSM_FUNC&  transm,
                           GenotypeIter        gbegin,
                           GenotypeIter        gend,
                           const PEN_FUNC&     pen);
                            
    ApproximateFamResidAdj(const SelfType&);
    
    ApproximateFamResidAdj& operator=(const SelfType&);
  //@}
    
    double calculate_adjustment
              (const FamilyType& fam,
               const Genotype&   mother_genotype,
               const Genotype&   father_genotype) const;
                         
  private:
  
    double calc_est_penetrance     (const MemberType&,
                                    const Genotype&,
                                    const Genotype&) const;
                                    
    double calc_est_susceptibility (const MemberType&,
                                    const Genotype&,
                                    const Genotype&) const;
                                    
    typedef boost::function <double (const MemberType&,
                                     const GENOTYPE&)>   SuscFuncType;
    typedef boost::function <double (const MemberType&,
                                     const GENOTYPE&,
                                     const GENOTYPE&,
                                     const GENOTYPE&)>   TransmFuncType;
    
    TransmFuncType my_transm_func; ///< Transmission Functor
    SuscFuncType   my_susc_func;   ///< Susceptibility Functor
    
    GenotypeIter my_gbegin;        ///< Genotype begin
    GenotypeIter my_gend;          ///< Genotype end
};

}
}

#include "pedcalc/fam_resid_adj.ipp"

#endif
