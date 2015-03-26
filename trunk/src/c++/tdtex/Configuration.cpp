#include "output/Output.h"
#include "tdtex/Configuration.h"

namespace SAGE  {
namespace TDTEX {

Configuration::Configuration() :
  my_marker                ((size_t)-1),
  my_trait                 ((size_t)-1),
  my_parent_trait          ((size_t)-1),
  my_method                (ALLELES),
  my_max_children          (1),
  my_max_sib_pairs         (0),
  my_skip_permutation_test (true),
  my_skip_mc_test          (true),
  my_skip_mcmh_test        (true),
  my_sex_differential      (false),
  my_ofilename             ("tdtex.out")
{ }
    
Configuration::Configuration(const Configuration & other) :
  my_marker                (other.my_marker),
  my_trait                 (other.my_trait),
  my_parent_trait          (other.my_parent_trait),
  my_method                (other.my_method),
  my_max_children          (other.my_max_children),
  my_max_sib_pairs         (other.my_max_sib_pairs),
  my_skip_permutation_test (other.my_skip_permutation_test),
  my_skip_mc_test          (other.my_skip_mc_test),
  my_skip_mcmh_test        (other.my_skip_mcmh_test),
  my_sex_differential      (other.my_sex_differential),
  my_ofilename             (other.my_ofilename)
{
}

Configuration & 
Configuration::operator= (const Configuration & other)
{
  if(this != &other)
  {
    my_marker                = other.my_marker;
    my_trait                 = other.my_trait;
    my_parent_trait          = other.my_parent_trait;
    my_method                = other.my_method;
    my_max_children          = other.my_max_children;
    my_max_sib_pairs         = other.my_max_sib_pairs;
    my_skip_permutation_test = other.my_skip_permutation_test;
    my_skip_mc_test          = other.my_skip_mc_test;
    my_skip_mcmh_test        = other.my_skip_mcmh_test;
    my_sex_differential      = other.my_sex_differential;
    my_ofilename             = other.my_ofilename;
  }
  
  return *this;
}
  
void
Configuration::dump() const
{
  std::cout << (OUTPUT::Table("Configuration dump")
    << (OUTPUT::TableRow() << "Marker" << get_marker())
    << (OUTPUT::TableRow() << "Trait"  << get_trait())
    << (OUTPUT::TableRow() << "Parent" << get_parent_trait())
    << (OUTPUT::TableRow() << "method" << (get_method() == ALLELES ? "alleles" : "genotypes"))
    << (OUTPUT::TableRow() << "max children" << get_max_children())
    << (OUTPUT::TableRow() << "max sibpairs" << get_max_sib_pairs())
    << (OUTPUT::TableRow() << "sex differential" << get_sex_differential())
    << (OUTPUT::TableRow() << "skip mc test" << get_skip_mc_test())
    << (OUTPUT::TableRow() << "skip mcmh test" << get_skip_mcmh_test())
    << (OUTPUT::TableRow() << "skip perm test" << get_skip_permutation_test())
    << (OUTPUT::TableRow() << "ofilename" << get_ofilename()))
    << std::flush;
}


} // End namespace TDTEX
} // End namespace SAGE

