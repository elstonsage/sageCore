//============================================================================
//  File:       genotype.ipp
//
//  Author:     Bob Steagall
//
//  History:    Version 1.00.00
//              X, Y-linkage added   - yjs  Mar. 2002
//
//  Notes:    
//
//  Copyright (c) 1999 R.C. Elston
//  All Rights Reserved

#ifndef GENOTYPE_H
#include "mlocus/genotype.h"
#endif

namespace SAGE   {
namespace MLOCUS {

namespace PRIVATE {
//============================================================================
//  IMPLEMENTATION: allele_info
//============================================================================
//

inline
allele_info::allele_info()
  : name(), frequency(0.0), id(NPOS)
{}

inline
allele_info::allele_info(const string& n, double f, uint i)
  : name(n), frequency(f), id(i), my_sex_type(st_NORMAL)
{}

inline
allele_info::allele_info(SexType sex_type)
  : frequency(1.0), id(SEX_ALLELE_ID), my_sex_type(sex_type)
{
       if(my_sex_type == st_NULL_Y) name = "~Y";
  else if(my_sex_type == st_NULL_X) name = "~X";
  else assert(false);
}

inline
allele_info::allele_info(const allele_info& a)
  : name(a.name), frequency(a.frequency), id(a.id), my_sex_type(a.my_sex_type)
{}

inline
allele_info& allele_info::operator=(const allele_info& a)
{
    if(this != &a)
    {
      name      = a.name;
      frequency = a.frequency;
      id        = a.id;
      
      my_sex_type = a.my_sex_type;
    }

    return *this;
}

inline bool
allele_info::operator==(const allele_info& ai) const
{
    //lint --e{777}
    return name == ai.name  &&  frequency == ai.frequency  &&  id == ai.id &&
           my_sex_type == ai.my_sex_type;
}

inline bool
allele_info::operator!=(const allele_info& ai) const
{
  return !(*this == ai);
}

//============================================================================
//  IMPLEMENTATION: phased_genotype_info
//============================================================================
//
inline
genotype_info::genotype_info()
  : my_alleles(),
    my_phased_id(NPOS),
    my_flipped_id(NPOS),
    my_unphased_id(NPOS),
    my_ginfo(NULL)
{}

inline
genotype_info::genotype_info
    (const allele_info*         a1,
     const allele_info*         a2,
     uint                       phased_id,
     uint                       flipped_id,
     uint                       unphased_id,
     const genotype_model_info* info)
  : my_alleles(make_pair(a1,a2)),
    my_phased_id(phased_id),
    my_flipped_id(flipped_id),
    my_unphased_id(unphased_id),
    my_ginfo(info)
{ }

inline
genotype_info::genotype_info(const genotype_info& other)
  : my_alleles(other.my_alleles),
    my_phased_id(other.my_phased_id),
    my_flipped_id(other.my_flipped_id),
    my_unphased_id(other.my_unphased_id),
    my_ginfo(other.my_ginfo)
{ }

inline
genotype_info& genotype_info::operator=
    (const genotype_info& other)
{
  if(this != &other)
  {
    my_alleles     = other.my_alleles;
    my_phased_id   = other.my_phased_id;
    my_flipped_id  = other.my_flipped_id;
    my_unphased_id = other.my_unphased_id;
    my_ginfo       = other.my_ginfo;
  }
  return *this;
}

//============================================================================
//  IMPLEMENTATION: genotype_model_info
//============================================================================
//

inline
genotype_model_info::~genotype_model_info()
{ }

inline phased_genotype
genotype_model_info::get_phased_genotype(allele a1, allele a2) const
{
    if(!a1.is_valid() || !a2.is_valid()) return phased_genotype();

    const PRIVATE::PhasedGenotypeByAlleles& allele_pg =
            my_phased_genotypes.get<PRIVATE::alleles>();
            
    PRIVATE::PhasedGenotypeByAlleles::const_iterator i =
      allele_pg.find(make_pair(a1.my_info, a2.my_info));
  
    if(i != allele_pg.end()) return phased_genotype(&*i);

    return phased_genotype();
}

inline phased_genotype
genotype_model_info::get_phased_genotype(uint id) const
{
    if(id >= my_phased_genotypes.size()) return phased_genotype();
    
    return phased_genotype(&*my_phased_genotypes_quick_lookup[id]);
}

inline unphased_genotype
genotype_model_info::get_unphased_genotype(uint id) const
{
    if(id >= my_unphased_genotypes.size()) return unphased_genotype();
    
    return unphased_genotype(&*my_unphased_genotypes_quick_lookup[id]);
}

}

//============================================================================
//  IMPLEMENTATION: allele
//============================================================================
//
inline allele::allele()
  : my_info(&PRIVATE::invalid_allele_info)
{}


//lint -e{1554}
inline allele::allele(const allele& a)
  : my_info(a.my_info)
{}


inline allele::allele(const PRIVATE::allele_info* info)
  : my_info(info)
{}

inline allele::~allele()
{
  my_info = &PRIVATE::invalid_allele_info;
}

inline allele&
allele::operator=(const allele& a)
{
    if(&a != this)
    {
        my_info = a.my_info;
    }
    return *this;
}


inline bool
allele::operator==(const allele& rhs) const
{
    return my_info == rhs.my_info;
}

inline bool
allele::operator!=(const allele& rhs) const
{
    return !(*this == rhs);
}

inline bool
allele::operator<(const allele& rhs) const
{
    return is_valid() &&
           rhs.is_valid() && 
           my_info->id < rhs.my_info->id;
}

inline double
allele::frequency() const
{
    return my_info->frequency;
}

inline uint
allele::id() const
{
    return my_info->id;
}

inline const string&
allele::name() const
{
    return my_info->name;
}

inline bool
allele::is_valid() const
{
  return my_info != &PRIVATE::invalid_allele_info && my_info->id != NPOS;
}


/// Returns \c true if the allele is sex specific, ie, a null allele on the
/// opposing sex chromosome from the one of interest.  Returns \c false otherwise
inline bool
allele::is_sex_allele() const
{
  return my_info->my_sex_type != PRIVATE::allele_info::st_NORMAL;
}

/// Returns \c true if the allele represents the null state on the Y chromosome for
/// x-linked markers, \c false otherwise
inline bool
allele::is_null_y_allele() const
{
  return my_info->my_sex_type == PRIVATE::allele_info::st_NULL_Y;
}

/// Returns \c true if the allele represents the null state on the X chromosome for
/// Y-linked markers, \c false otherwise
inline bool
allele::is_null_x_allele() const
{
  return my_info->my_sex_type == PRIVATE::allele_info::st_NULL_X;
}

//============================================================================
//  IMPLEMENTATION: phased_genotype
//============================================================================
//
inline phased_genotype::phased_genotype()
  : my_info(&PRIVATE::invalid_phased_genotype_info)
{}


//lint -e{1554} <- pointer copy ok
inline phased_genotype::phased_genotype(const phased_genotype& g)
  : my_info(g.my_info)
{}


inline phased_genotype::phased_genotype
    (const PRIVATE::genotype_info* info)
  : my_info(info)
{}


inline phased_genotype::~phased_genotype()
{
  my_info = NULL;
}


inline phased_genotype&
phased_genotype::operator=(const phased_genotype& g)
{
    if (&g != this)
    {
        my_info = g.my_info;
    }
    return *this;
}


inline bool
phased_genotype::equivalent(const phased_genotype& g) const
{
    return my_info == g.my_info;
}


inline bool
phased_genotype::equivalent(const unphased_genotype& g) const
{
    return my_info->my_unphased_id == g.my_info->my_unphased_id;
}


inline string
phased_genotype::name(char c) const
{
    return allele1().name() + c + allele2().name();
}


inline bool
phased_genotype::operator==(const phased_genotype& rhs) const
{
    return my_info == rhs.my_info;
}

inline bool
phased_genotype::operator!=(const phased_genotype& rhs) const
{
    return !(*this == rhs);
}

inline bool
phased_genotype::operator<(const phased_genotype& rhs) const
{
    return my_info < rhs.my_info;
}

inline bool
phased_genotype::homozygous() const
{
    return allele1() == allele2();
}

inline allele
phased_genotype::allele1() const
{
    return allele(my_info->my_alleles.first);
}

inline allele
phased_genotype::allele2() const
{
    return allele(my_info->my_alleles.second);
}

inline const string
phased_genotype::name() const
{
    return allele1().name() + my_info->my_ginfo->separators[PhasedForward] +
           allele2().name();
}

inline double
phased_genotype::frequency() const
{
    return allele1().frequency() * allele2().frequency();
}

inline uint
phased_genotype::get_id() const
{
    return my_info->my_phased_id;
}

inline bool
phased_genotype::is_flippable() const
{
  return my_info->my_phased_id != my_info->my_flipped_id;
}

inline
phased_genotype phased_genotype::get_flipped_phased_genotype() const
{
  return my_info->my_ginfo->get_phased_genotype(my_info->my_flipped_id);
}

inline
unphased_genotype phased_genotype::get_equivalent_unphased_genotype() const
{
  return my_info->my_ginfo->get_unphased_genotype(my_info->my_unphased_id);
}

inline bool
phased_genotype::is_valid() const
{
  return my_info != &PRIVATE::invalid_phased_genotype_info;
}

inline bool
phased_genotype::is_sex_specific()      const
{
  return my_info->my_ginfo->my_type != AUTOSOMAL;
}

inline bool
phased_genotype::is_male_compatible()   const
{
  if(!is_valid()) return false;
  
  switch(my_info->my_ginfo->my_type)
  {
    case AUTOSOMAL : return true;
    case X_LINKED  : return  allele2().is_sex_allele();
    case Y_LINKED  : return !allele2().is_sex_allele();
  }
  return false;
}

inline bool
phased_genotype::is_female_compatible() const
{
  if(!is_valid()) return false;
  
  switch(my_info->my_ginfo->my_type)
  {
    case AUTOSOMAL : return true;
    case X_LINKED  : return !allele2().is_sex_allele();
    case Y_LINKED  : return  allele2().is_sex_allele();
  }
  return false;
}

//============================================================================
//  IMPLEMENTATION: unphased_genotype
//============================================================================
//
inline unphased_genotype::unphased_genotype()
  : my_info(&PRIVATE::invalid_unphased_genotype_info)
{}


//lint -e{1554}
inline unphased_genotype::unphased_genotype(const unphased_genotype& g)
  : my_info(g.my_info)
{}


inline unphased_genotype::unphased_genotype
    (const PRIVATE::genotype_info* info)
  : my_info(info)
{}


inline unphased_genotype::~unphased_genotype()
{
  my_info = NULL;
}


inline unphased_genotype&
unphased_genotype::operator=(const unphased_genotype& g)
{
    //lint -e{1555}

    if (&g != this)
    {
        my_info = g.my_info;
    }
    return *this;
}


inline bool
unphased_genotype::equivalent(const phased_genotype& g) const
{
    return g.equivalent(*this);
}


inline bool
unphased_genotype::equivalent(const unphased_genotype& g) const
{
    return my_info == g.my_info;
}


inline string
unphased_genotype::name(char c) const
{
    return allele1().name() + c + allele2().name();
}


inline bool
unphased_genotype::operator==(const unphased_genotype& rhs) const
{
    return my_info == rhs.my_info;
}

inline bool
unphased_genotype::operator!=(const unphased_genotype& rhs) const
{
    return !(*this == rhs);
}

inline bool
unphased_genotype::operator<(const unphased_genotype& rhs) const
{
    return my_info < rhs.my_info;
}

inline bool
unphased_genotype::homozygous() const
{
    return allele1() == allele2();
}

inline size_t
unphased_genotype::get_equivalent_phased_genotype_count() const
{
  return (my_info->my_phased_id == my_info->my_flipped_id) ? 1 : 2;
}

inline allele
unphased_genotype::allele1() const
{
    return allele(my_info->my_alleles.first);
}

inline allele
unphased_genotype::allele2() const
{
    return allele(my_info->my_alleles.second);
}

inline string
unphased_genotype::name() const
{
    return allele1().name() + my_info->my_ginfo->separators[0] + allele2().name();
}

inline uint
unphased_genotype::get_id() const
{
    return my_info->my_unphased_id;
}

inline double
unphased_genotype::frequency() const
{
    double d = allele1().frequency() * allele2().frequency();

    if(allele1() != allele2()) d *= 2.0;

    return d;
}

inline phased_genotype
unphased_genotype::get_equivalent_phased_genotype1() const
{
  return my_info->my_ginfo->get_phased_genotype(my_info->my_phased_id);
}

inline phased_genotype
unphased_genotype::get_equivalent_phased_genotype2() const
{
  return my_info->my_ginfo->get_phased_genotype(my_info->my_flipped_id);
}

inline bool
unphased_genotype::is_valid() const
{
  return my_info != &PRIVATE::invalid_unphased_genotype_info;
}
    
inline bool
unphased_genotype::is_sex_specific()      const
{
  return my_info->my_ginfo->my_type != AUTOSOMAL;
}

inline bool 
unphased_genotype::is_male_compatible()   const
{
  if(!is_valid()) return false;
  
  switch(my_info->my_ginfo->my_type)
  {
    case AUTOSOMAL : return true;
    case X_LINKED  : return  allele2().is_sex_allele();
    case Y_LINKED  : return !allele2().is_sex_allele();
  }
  return false;
}

inline bool
unphased_genotype::is_female_compatible() const
{
  if(!is_valid()) return false;
  
  switch(my_info->my_ginfo->my_type)
  {
    case AUTOSOMAL : return true;
    case X_LINKED  : return !allele2().is_sex_allele();
    case Y_LINKED  : return  allele2().is_sex_allele();
  }
  return false;
}

//============================================================================
//  IMPLEMENTATION: allele_iterator
//============================================================================
//

inline
allele_iterator::allele_iterator()
  : allele_iterator::iterator_adaptor_()
{}

inline
allele_iterator::allele_iterator(const allele_iterator& i)
  : allele_iterator::iterator_adaptor_(i)
{}

inline
allele_iterator::allele_iterator
  (const std::vector<PRIVATE::allele_info>::const_iterator& c)
  : allele_iterator::iterator_adaptor_(c)
{ }

inline
allele allele_iterator::dereference() const
{
  return allele(&*this->base_reference());
}

//============================================================================
//  IMPLEMENTATION: phased_genotype_iterator
//============================================================================
//

inline
phased_genotype_iterator::phased_genotype_iterator()
  : phased_genotype_iterator::iterator_adaptor_()
{}

inline
phased_genotype_iterator::phased_genotype_iterator
    (const phased_genotype_iterator& i)
  : phased_genotype_iterator::iterator_adaptor_(i)
{ }

inline
phased_genotype_iterator::phased_genotype_iterator
  (const PRIVATE::PhasedGenotypeBySequence::const_iterator& c)
  : phased_genotype_iterator::iterator_adaptor_(c)
{ }

inline
phased_genotype phased_genotype_iterator::dereference() const
{
  return phased_genotype(&*this->base_reference());
}
  
//============================================================================
//  IMPLEMENTATION: unphased_genotype_iterator
//============================================================================
//

inline
unphased_genotype_iterator::unphased_genotype_iterator()
  : unphased_genotype_iterator::iterator_adaptor_()
{}

inline
unphased_genotype_iterator::unphased_genotype_iterator
    (const unphased_genotype_iterator& i)
  : unphased_genotype_iterator::iterator_adaptor_(i)
{ }

inline
unphased_genotype_iterator::unphased_genotype_iterator
  (const PRIVATE::UnphasedGenotypeBySequence::const_iterator& c)
  : unphased_genotype_iterator::iterator_adaptor_(c)
{ }

inline
unphased_genotype unphased_genotype_iterator::dereference() const
{
  return unphased_genotype(&*this->base_reference());
}
  
//============================================================================
//  IMPLEMENTATION: genotype_model
//============================================================================
//

inline
genotype_model::genotype_model()
  : my_info(new PRIVATE::genotype_model_info())
{}


inline
genotype_model::genotype_model(const string& n)
  : my_info(new PRIVATE::genotype_model_info(n))
{}


inline
genotype_model::genotype_model(const genotype_model& gm)
  : my_info(gm.my_info)
{}


inline
genotype_model::~genotype_model()
{}


inline genotype_model&
genotype_model::operator=(const genotype_model& gm)
{
    if (&gm != this  &&  gm.my_info != my_info)
    {
        my_info = gm.my_info;
    }
    return *this;
}

inline allele
genotype_model::get_allele(uint i) const
{
  if(i == SEX_ALLELE_ID)   return get_sex_specific_allele();
  if(i >= allele_count())  return allele();
  
  return allele(&my_info->alleles[i]);
}

inline phased_genotype
genotype_model::get_phased_genotype(allele a1, allele a2) const
{
    return my_info->get_phased_genotype(a1, a2);
}

inline unphased_genotype
genotype_model::get_unphased_genotype(allele a1, allele a2) const
{
    if(!a1.is_valid() || !a2.is_valid()) return unphased_genotype();

    bool is_sex_based = a1.is_sex_allele() || a2.is_sex_allele();
    
    // This will put autosomal alleles in the right order.
    if(!is_sex_based && a2 < a1) std::swap(a1,a2);

    // Get our index based on alleles
    const PRIVATE::UnphasedGenotypeByAlleles& allele_pg =
            my_info->my_unphased_genotypes.get<PRIVATE::alleles>();

    // Look it up
    PRIVATE::UnphasedGenotypeByAlleles::const_iterator i =
      allele_pg.find(make_pair(a1.my_info, a2.my_info));
  
    if(i != allele_pg.end()) return unphased_genotype(&*i);

    if(is_sex_based)
    {
      // Attempt the opposite ordering for sex based genotyeps
      i = allele_pg.find(make_pair(a2.my_info, a1.my_info));
      
      if(i != allele_pg.end()) return unphased_genotype(&*i);
    }

    return unphased_genotype();
}

inline uint
genotype_model::allele_count() const
{
    return my_info->alleles.size();
}

inline uint
genotype_model::phased_genotype_count() const
{
    return my_info->my_phased_genotypes.size();
}

inline uint
genotype_model::unphased_genotype_count() const
{
    return my_info->my_unphased_genotypes.size();
}

inline allele_iterator
genotype_model::allele_begin() const
{
    return allele_iterator(my_info->alleles.begin());
}

inline allele_iterator
genotype_model::allele_end() const
{
    return allele_iterator(my_info->alleles.begin() + allele_count());
}

inline phased_genotype_iterator
genotype_model::phased_genotype_begin() const
{
    return phased_genotype_iterator(my_info->my_phased_genotypes.begin());
}

inline phased_genotype_iterator
genotype_model::phased_genotype_end() const
{
    return phased_genotype_iterator(my_info->my_phased_genotypes.end());
}

inline unphased_genotype_iterator
genotype_model::unphased_genotype_begin() const
{
    return unphased_genotype_iterator(my_info->my_unphased_genotypes.begin());
}

inline unphased_genotype_iterator
genotype_model::unphased_genotype_end() const
{
    return unphased_genotype_iterator(my_info->my_unphased_genotypes.end());
}

inline phased_genotype_iterator
genotype_model::phased_genotype_begin(const MPED::SexCode& s) const
{
  switch(s)
  {
    case MPED::SEX_MALE    : return phased_genotype_iterator(my_info->my_male_phased_genotype_begin);
    case MPED::SEX_FEMALE  : return phased_genotype_iterator(my_info->my_female_phased_genotype_begin);
    case MPED::SEX_MISSING :
    default                : return phased_genotype_begin();
  }
}

inline phased_genotype_iterator
genotype_model::phased_genotype_end(const MPED::SexCode& s) const
{
  switch(s)
  {
    case MPED::SEX_MALE    : return phased_genotype_iterator(my_info->my_male_phased_genotype_end);
    case MPED::SEX_FEMALE  : return phased_genotype_iterator(my_info->my_female_phased_genotype_end);
    case MPED::SEX_MISSING : 
    default                : return phased_genotype_end();
  }
}

inline unphased_genotype_iterator
genotype_model::unphased_genotype_begin(const MPED::SexCode& s) const
{
  switch(s)
  {
    case MPED::SEX_MALE    : return unphased_genotype_iterator(my_info->my_male_unphased_genotype_begin);
    case MPED::SEX_FEMALE  : return unphased_genotype_iterator(my_info->my_female_unphased_genotype_begin);
    case MPED::SEX_MISSING : 
    default                : return unphased_genotype_begin();
  }
}

inline unphased_genotype_iterator
genotype_model::unphased_genotype_end(const MPED::SexCode& s) const
{
  switch(s)
  {
    case MPED::SEX_MALE    : return unphased_genotype_iterator(my_info->my_male_unphased_genotype_end);
    case MPED::SEX_FEMALE  : return unphased_genotype_iterator(my_info->my_female_unphased_genotype_end);
    case MPED::SEX_MISSING : 
    default                : return unphased_genotype_end();
  }
}

inline phased_genotype
genotype_model::get_phased_genotype(uint id) const
{
  return my_info->get_phased_genotype(id);
}

inline unphased_genotype
genotype_model::get_unphased_genotype(uint id) const
{
  return my_info->get_unphased_genotype(id);
}

inline allele
genotype_model::get_sex_specific_allele() const
{
  if(my_info->my_type == AUTOSOMAL) return allele(&PRIVATE::invalid_allele_info);
  return allele(&my_info->my_sex_allele_info);
}

inline const string&
genotype_model::name() const
{
    return my_info->name;
}

inline const string&
genotype_model::missing_allele_name() const
{
    return my_info->missing_allele_name;
}

inline bool
genotype_model::dynamic_alleles() const
{
    return my_info->dynamic_alleles;
}

inline char
genotype_model::backward_separator() const
{
    return my_info->separators[2];
}

inline char
genotype_model::forward_separator() const
{
    return my_info->separators[1];
}

inline char
genotype_model::unphased_separator() const
{
    return my_info->separators[0];
}

inline string
genotype_model::separators() const
{
    return string(my_info->separators, 3);
}

inline void
genotype_model::set_name(const string& newname)
{
    uniquify();
    my_info->name = newname;
}

inline void
genotype_model::set_dynamic_alleles(bool d)
{
    uniquify();
    my_info->dynamic_alleles = d;
}

inline void
genotype_model::set_backward_separator(char c)
{
    uniquify();
    my_info->separators[2] = c;
}


inline void
genotype_model::set_forward_separator(char c)
{
    uniquify();
    my_info->separators[1] = c;
}


inline void
genotype_model::set_unphased_separator(char c)
{
    uniquify();
    my_info->separators[0] = c;
}

inline void
genotype_model::set_separators(const string& seps)
{
    uniquify();

    if (seps.size() >= 3)
    {
        set_backward_separator(seps[2]);
        set_forward_separator(seps[1]);
        set_unphased_separator(seps[0]);
    }
}

inline void
genotype_model::set_missing_allele_name(const string& n)
{
    string  nm = parse_allele_name(n);

    if (nm.size() > 0)
    {
        uniquify();
        my_info->missing_allele_name = nm;
    }
}



inline bool
genotype_model::is_x_linked() const
{
  return my_info->my_type == X_LINKED;
}

inline bool
genotype_model::is_y_linked() const
{
  return my_info->my_type == Y_LINKED;
}

inline bool
genotype_model::is_autosomal() const
{
  return my_info->my_type == AUTOSOMAL;
}


inline GenotypeModelType
genotype_model::get_model_type() const
{
  return my_info->my_type;
}

inline void
genotype_model::set_model_type(GenotypeModelType gt)
{
  if(gt == my_info->my_type) return;
  
  uniquify();

  my_info->configure_sex_type(gt);
}

// Added for allele_frequency adjustment - yjs Dec. 2002
//----------
//
inline void
genotype_model::modify_allele_frequency(allele al, double freq)
{
  my_info->alleles[al.id()].frequency = freq;
}

inline void
genotype_model::modify_allele_frequency(allele_iterator ai, double freq)
{
  modify_allele_frequency(*ai, freq);
}

inline void
genotype_model::modify_allele_frequency(const string& name, double freq)
{
  modify_allele_frequency(get_allele(name), freq);
}

//============================================================================
//  IMPLEMENTATION: child_genotype_set
//============================================================================
//
inline child_genotype_set::iterator
child_genotype_set::begin() const
{
    return iterator(my_genotypes);
}

inline child_genotype_set::iterator
child_genotype_set::end() const
{
    return iterator(my_genotypes + my_size);
}

inline const phased_genotype&
child_genotype_set::operator[](uint i) const
{
    return my_genotypes[i];
}

inline bool
child_genotype_set::contains(const phased_genotype& p) const
{
  for(iterator i = begin(); i != end(); ++i)
    if(i->equivalent(p)) return true;

  return false;
}

inline bool
child_genotype_set::contains(const unphased_genotype& p) const
{
  for(iterator i = begin(); i != end(); ++i)
    if(i->equivalent(p)) return true;

  return false;
}

inline uint
child_genotype_set::size() const
{
  return my_size;
}

} // End namespace MLOCUS
} // End namespace SAGE

