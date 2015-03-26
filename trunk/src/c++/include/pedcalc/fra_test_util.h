#ifndef FRA_TEST_UTIL_H
#define FRA_TEST_UTIL_H

#include "pedcalc/fam_resid_adj.h"
#include "output/Output.h"
#include "mped/mp.h"
#include "boost/iterator/counting_iterator.hpp"
#include <iostream>
#include <sstream>


namespace SAGE
{
namespace PED_CALC
{

/// Generate an output of the familial residual adjustment for binary traits
/// for an entire multipedigree.
///
/// The output generated is a section with subsections for each pedigree
/// in the data set.
template <class FRA, class GENO_ITER>
OUTPUT::Section
    generate_fra_test_output
        (const FRA&                    fra,
         const typename FRA::MpedType& mp,
         const GENO_ITER&              gbegin,
         const GENO_ITER&              gend)
{
  OUTPUT::Section test_results(
      "Familial Binary Residual Adjustment Test on multipedigree");
      
  for(typename FRA::MpedType::pedigree_iterator p_iter = mp.pedigree_begin();
      p_iter != mp.pedigree_end();
      ++p_iter)
  {
    test_results << generate_fra_test_output(fra, *p_iter, gbegin, gend);
  }
  
  return test_results;
}

/// Generate an output of the familial residual adjustment for binary traits
/// for a pedigree.
///
/// The output generated is a section with subsections for each family
/// in the pedigree.
template <class FRA, class GENO_ITER>
OUTPUT::Section
    generate_fra_test_output
        (const FRA&                                   fra,
         const typename FRA::MpedType::pedigree_type& ped,
         const GENO_ITER&                             gbegin,
         const GENO_ITER&                             gend)
{
  OUTPUT::Section test_results(
      "Familial Binary Residuals for pedigree: " + ped.name());
      
  for(typename FRA::MpedType::family_const_iterator f_iter = ped.family_begin();
      f_iter != ped.family_end();
      ++f_iter)
  {
    test_results << generate_fra_test_output(fra, *f_iter, gbegin, gend);
  }
  
  return test_results;
}

template <typename T>
std::string convert_to_string(const T& t)
{
  std::stringstream s;
  
  s << t;
  
  return s.str();
}


template <class GENOTYPE, class MPTYPE, class GENO_ITER>
OUTPUT::Table
    generate_fra_test_output
        (const ExactFamResidAdj<GENOTYPE, MPTYPE>& fra,
         const typename MPTYPE::family_type&       fam,
         const GENO_ITER&                          gbegin,
         const GENO_ITER&                          gend)
{
  // Determine mother and father
  
  OUTPUT::Table test_results
      ("Family Adjustments for family: " + fam.name());

  const typename MPTYPE::member_type* mother = fam.parent1();
  const typename MPTYPE::member_type* father = fam.parent2();
  
  // I hate doing this, but there's no way to get the mother/father of a family!
  if(mother->is_male() || father->is_female())
    std::swap(mother, father);
  
  // Set up the columns in the output table
  
  test_results << OUTPUT::TableColumn("Mother: " + mother->name())
               << OUTPUT::TableColumn("Father: " + father->name());
               
  for(typename MPTYPE::offspring_const_iterator child = fam.offspring_begin();
      child != fam.offspring_end(); ++child)
  {
    test_results << OUTPUT::TableColumn("Child: " + child->name());
  }
  
  test_results << OUTPUT::TableColumn("Value");
  
  // Generate results to populate the table
  for(GENO_ITER mgeno = gbegin; mgeno != gend; ++mgeno)
  {
    for(GENO_ITER fgeno = gbegin; fgeno != gend; ++fgeno)
    {
      // We can't do *every* child model, so we'll try just a few:
      // 1.  All children the same, for each genotype
      // 2.  Each child with a genotype, the first child having the first
      //     genotype, second the second, etc.  If there aren't enough genotypes
      //     for each child, just repeat genotypes following the pattern.
      // 3.  Same as #2, but first having the *last* genotype, second last - 1,
      //     etc.
      
      // 1.
      for(GENO_ITER i = gbegin; i != gend; ++i)
      {
        std::vector<GENOTYPE> ch_genos(fam.offspring_count(), (GENOTYPE) *i);
        
        OUTPUT::TableRow r;
        
        r << convert_to_string(*mgeno) << convert_to_string(*fgeno);
        
        for(size_t ch = 0; ch != fam.offspring_count(); ++ch)
          r << convert_to_string(*i);
        
        r << fra.calculate_adjustment(fam,
                                      (GENOTYPE) *mgeno, 
                                      (GENOTYPE) *fgeno,
                                      ch_genos);
        
        test_results << r;
      }
      
     // 2.
      {
        std::vector<GENOTYPE> ch_genos(fam.offspring_count(), (GENOTYPE) *gbegin);
        
        GENO_ITER g = gbegin;
        for(size_t ch = 0; ch != fam.offspring_count(); ++ch, ++g)
        {
          if(g == gend) g = gbegin;
          
          ch_genos[ch] = (GENOTYPE) *g;
        }
        
        OUTPUT::TableRow r;
        
        r << convert_to_string(*mgeno) << convert_to_string(*fgeno);

        for(size_t ch = 0; ch != fam.offspring_count(); ++ch)
          r << convert_to_string(ch_genos[ch]);
        
        r << fra.calculate_adjustment(fam, (GENOTYPE) *mgeno, (GENOTYPE) *fgeno, ch_genos);
        
        test_results << r;
      }
      
      // 3.
      {
        std::vector<GENOTYPE> ch_genos(fam.offspring_count(), (GENOTYPE) *gbegin);
        
        GENO_ITER g = gend;
        for(size_t ch = 0; ch != fam.offspring_count(); ++ch)
        {
          if(g == gbegin) g = gend;
          
          --g;
          
          ch_genos[ch] = (GENOTYPE) *g;
        }
        
        OUTPUT::TableRow r;
        
        r << convert_to_string(*mgeno) << convert_to_string(*fgeno);
        
        for(size_t ch = 0; ch != fam.offspring_count(); ++ch)
          r << convert_to_string(ch_genos[ch]);
        
        r << fra.calculate_adjustment(fam, (GENOTYPE) *mgeno, (GENOTYPE) *fgeno, ch_genos);
        
        test_results << r;
      }
    }
  }
  
  return test_results;
}

template <class GENOTYPE, class MPTYPE, class GENO_ITER>
OUTPUT::Table
    generate_fra_test_output
        (const ApproximateFamResidAdj<GENOTYPE, MPTYPE, GENO_ITER>& fra,
         const typename MPTYPE::family_type&                        fam,
         const GENO_ITER&                                           gbegin,
         const GENO_ITER&                                           gend)
{
  // Determine mother and father
  
  OUTPUT::Table test_results
      ("Approximate Family Adjustments for family: " + fam.name());

  const typename MPTYPE::member_type* mother = fam.parent1();
  const typename MPTYPE::member_type* father = fam.parent2();
  
  // I hate doing this, but there's no way to get the mother/father of a family!
  if(mother->is_male() || father->is_female())
    std::swap(mother, father);
  
  // Set up the columns in the output table
  
  test_results << OUTPUT::TableColumn("Mother: " + mother->name())
               << OUTPUT::TableColumn("Father: " + father->name());
               
  test_results << OUTPUT::TableColumn("Value");
  
  // Generate results to populate the table
  for(GENO_ITER mgeno = gbegin; mgeno != gend; ++mgeno)
  {
    for(GENO_ITER fgeno = gbegin; fgeno != gend; ++fgeno)
    {
      OUTPUT::TableRow r;
      r << convert_to_string(*mgeno)
        << convert_to_string(*fgeno) 
        << fra.calculate_adjustment(fam, (GENOTYPE) *mgeno,
                                         (GENOTYPE) *fgeno);
                                         
      test_results << r;
    }
  }
  
  return test_results;
}

}
}

#endif
