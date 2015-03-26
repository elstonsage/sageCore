//============================================================================
//  File:       penmodel.ipp
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

#ifndef PENMODEL_H
#include "mlocus/penmodel.h"
#endif

namespace SAGE   {
namespace MLOCUS {

//============================================================================
//  IMPLEMENTATION: penetrance_model_info
//============================================================================
//
inline penetrance_model_info::penetrance_model_info()
  : gmodel(), phmodel(), name(), codominant(0), strict_codominant(0)
{
  init();
}


inline penetrance_model_info::penetrance_model_info(const penetrance_model_info& p)
  : gmodel                (p.gmodel),
    phmodel               (p.phmodel),
    codominant_phenotypes (p.codominant_phenotypes),
    strict_phenotypes     (p.strict_phenotypes),
    name                  (p.name), 
    phased_penetrance     (p.phased_penetrance),
    unphased_penetrance   (p.unphased_penetrance),
    codominant            (p.codominant),
    strict_codominant     (p.strict_codominant),
    last_alleles          (p.last_alleles)
{ }


inline penetrance_model_info::penetrance_model_info(const string& n)
  : gmodel(), phmodel(), name(n), codominant(0), strict_codominant(0)
{
  init();
}


inline penetrance_model_info::penetrance_model_info(const genotype_model& gm)
  : gmodel(gm), phmodel(gmodel), name(), codominant(0), strict_codominant(0)
{
  init();
}


inline penetrance_model_info::penetrance_model_info
(const genotype_model& gm, const string& n)
  : gmodel(gm), phmodel(gmodel), name(n), codominant(0), strict_codominant(0)
{
  init();
}

inline penetrance_model_info::penetrance_model_info
(const genotype_model& gm, const phenotype_model& pm)
  : gmodel(gm), phmodel(pm), name(), codominant(0), strict_codominant(0)
{
  init();
}

inline penetrance_model_info::penetrance_model_info
(const genotype_model& gm, const phenotype_model& pm, const string& n)
  : gmodel(gm), phmodel(pm), name(n), codominant(0), strict_codominant(0)
{
  init();
}

inline void
penetrance_model_info::init()
{
  //lint --e{534} <- Don't care about set_default_value()s return value
  unphased_penetrance.set_default_value(0.0);
  phased_penetrance.set_default_value(0.0);
  last_alleles = 0;
}

//============================================================================
//  IMPLEMENTATION: unphased_penetrance_iterator
//============================================================================
//
inline bool
penetrance_model::unphased_penetrance_iterator::is_phased() const
{
  // If we're homozygous, there is nothing to do.
  
  if(unphased_geno().homozygous()) return false;

  // Determine which of the two penetrance values for the phased genotypes
  // is != 0
  
  bool pen1 = my_host->phased_penetrance(phenotype_id(),  geno_id()) != 0.0;
  bool pen2 = my_host->phased_penetrance(phenotype_id(), -geno_id()) != 0.0;

  // Return true if only one of pen1 and pen2 are true.
  
  return ((pen1 && !pen2) || (!pen1 && pen2));
}

//============================================================================
//  IMPLEMENTATION: penetrance_model
//============================================================================
//
inline uint
penetrance_model::allele_count() const
{
    return my_info->gmodel.allele_count();
}

inline uint
penetrance_model::phased_genotype_count() const
{
    return my_info->gmodel.phased_genotype_count();
}

inline uint
penetrance_model::unphased_genotype_count() const
{
    return my_info->gmodel.unphased_genotype_count();
}

inline uint
penetrance_model::phenotype_count() const
{
    return my_info->phmodel.phenotype_count();
}

//----------
//
inline allele_iterator
penetrance_model::allele_begin() const
{
    return my_info->gmodel.allele_begin();
}

inline allele_iterator
penetrance_model::allele_end() const
{
    return my_info->gmodel.allele_end();
}

//----------
//
inline phased_genotype_iterator
penetrance_model::phased_genotype_begin() const
{
    return my_info->gmodel.phased_genotype_begin();
}

inline phased_genotype_iterator
penetrance_model::phased_genotype_end() const
{
    return my_info->gmodel.phased_genotype_end();
}

inline unphased_genotype_iterator
penetrance_model::unphased_genotype_begin() const
{
    return my_info->gmodel.unphased_genotype_begin();
}

inline unphased_genotype_iterator
penetrance_model::unphased_genotype_end() const
{
    return my_info->gmodel.unphased_genotype_end();
}

//----------
//
inline penetrance_model::phenotype_iterator
penetrance_model::phenotype_begin() const
{
    return my_info->phmodel.phenotype_begin();
}

inline penetrance_model::phenotype_iterator
penetrance_model::phenotype_end() const
{
    return my_info->phmodel.phenotype_end();
}

//----------
//
inline penetrance_model::phased_penetrance_iterator
penetrance_model::phased_penetrance_begin(uint ptid) const
{
    return phased_penetrance_iterator(*this, ptid, false);
}

inline penetrance_model::phased_penetrance_iterator
penetrance_model::phased_penetrance_end(uint ptid) const
{
    return phased_penetrance_iterator(*this, ptid, true);
}

inline penetrance_model::phased_penetrance_iterator
penetrance_model::phased_penetrance_begin(const phenotype& p) const
{
    return phased_penetrance_iterator(*this, p.id(), false);
}

inline penetrance_model::phased_penetrance_iterator
penetrance_model::phased_penetrance_end(const phenotype& p) const
{
    return phased_penetrance_iterator(*this, p.id(), true);
}

inline penetrance_model::unphased_penetrance_iterator
penetrance_model::unphased_penetrance_begin(uint ptid) const
{
    return unphased_penetrance_iterator(*this, ptid, false);
}

inline penetrance_model::unphased_penetrance_iterator
penetrance_model::unphased_penetrance_end(uint ptid) const
{
    return unphased_penetrance_iterator(*this, ptid, true);
}

inline penetrance_model::unphased_penetrance_iterator
penetrance_model::unphased_penetrance_begin(const phenotype& p) const
{
    return unphased_penetrance_iterator(*this, p.id(), false);
}

inline penetrance_model::unphased_penetrance_iterator
penetrance_model::unphased_penetrance_end(const phenotype& p) const
{
    return unphased_penetrance_iterator(*this, p.id(), true);
}

inline uint
penetrance_model::phased_penetrance_count(uint id) const
{
  return (uint) my_info->phased_penetrance.row_elements((int) id);
}

inline uint
penetrance_model::phased_penetrance_count(const phenotype& p) const
{
    return phased_penetrance_count(p.id());
}

inline uint
penetrance_model::unphased_penetrance_count(uint id) const
{
  return (uint) my_info->unphased_penetrance.row_elements((int) id);
}

inline uint
penetrance_model::unphased_penetrance_count(const phenotype& p) const
{
    return unphased_penetrance_count(p.id());
}



//----------
//
inline allele 
penetrance_model::get_allele(uint id) const
{
  return gmodel().get_allele(id);
}

inline phased_genotype
penetrance_model::get_phased_genotype(int id) const
{
    return gmodel().get_phased_genotype(id);
}

inline unphased_genotype
penetrance_model::get_unphased_genotype(uint id) const
{
    return gmodel().get_unphased_genotype(id);
}

inline const phenotype&
penetrance_model::get_phenotype(uint id) const
{
    return my_info->phmodel.get_phenotype(id);
}

//----------
//
inline allele
penetrance_model::get_allele(const string& n) const
{
    return my_info->gmodel.get_allele(n);
}

inline phased_genotype
penetrance_model::get_phased_genotype(const string& n) const
{
    return my_info->gmodel.get_phased_genotype(n);
}

inline unphased_genotype
penetrance_model::get_unphased_genotype(const string& n) const
{
    return my_info->gmodel.get_unphased_genotype(n);
}

inline phenotype
penetrance_model::get_phenotype(const string& n) const
{
    return my_info->phmodel.get_phenotype(n);
}


inline phenotype
penetrance_model::get_missing_phenotype() const
{
    return my_info->phmodel.get_missing_phenotype();
}



//----------
//
inline uint
penetrance_model::get_phenotype_id(const string& n) const
{
    return my_info->phmodel.get_phenotype_id(n);
}

inline uint
penetrance_model::get_missing_phenotype_id() const
{
    return my_info->phmodel.get_missing_phenotype_id();
}

//----------
//
inline const genotype_model& penetrance_model::gmodel() const
{
    return my_info->gmodel;
}

//lint -e{1762} <- This is the non-const access to the gmodel
inline genotype_model& penetrance_model::gmodel()
{
    return my_info->gmodel;
}

inline const phenotype_model& penetrance_model::phmodel() const
{
    return my_info->phmodel;
}

//lint -e{1762} <- This is the non-const access to the phmodel
inline phenotype_model& penetrance_model::phmodel()
{
    return my_info->phmodel;
}

inline const string&
penetrance_model::name() const
{
    return my_info->name;
}

inline const string&
penetrance_model::set_name(const string& s)
{
  uniquify();

  my_info->name = s;

  return my_info->name;
}

inline const string&
penetrance_model::missing_allele_name() const
{
    return my_info->gmodel.missing_allele_name();
}

inline const string&
penetrance_model::missing_phenotype_name() const
{
    return my_info->phmodel.missing_phenotype_name();
}

inline bool
penetrance_model::genotype_informative_phenotype(uint id) const
{
  uint phased_count   = phased_penetrance_count(id);

  if(phased_count > 0 && phased_count < phased_genotype_count()) return true;

  uint unphased_count = unphased_penetrance_count(id);

  if(unphased_count > 0 && unphased_count < unphased_genotype_count()) return true;

  return false;
}

inline bool
penetrance_model::genotype_informative_phenotype(const phenotype& p) const
{
    return genotype_informative_phenotype(p.id());
}

inline bool
penetrance_model::penetrance_informative_phenotype(const phenotype& p) const
{
    return penetrance_informative_phenotype(p.id());
}

inline bool
penetrance_model::strict_phenotype(uint id) const
{
  //lint -e{56} -e{48} -e{734} <- Spurious messages
  return my_info->strict_phenotypes[id];
}

inline bool
penetrance_model::strict_phenotype(const phenotype& p) const
{
  //lint -e{56} -e{48} -e{734} <- Spurious messages
  return my_info->strict_phenotypes[p.id()];
}

inline bool
penetrance_model::codominant(bool strict) const
{
  if(strict)
    return my_info->strict_codominant == 0;
  else
    return my_info->codominant == 0;
}

inline bool
penetrance_model::codominant(uint ptid, bool strict) const
{
  //lint -e{56} -e{48} -e{734} <- Spurious messages
  if(strict || my_info->strict_phenotypes[ptid])
    return my_info->codominant_phenotypes[ptid];

  return true;
}

inline bool
penetrance_model::genotype_informative() const
{
  for(uint pid = 0; pid < phenotype_count()+1; ++pid)
  {
    if(genotype_informative_phenotype(pid)) return true;
  }

  return false;
}

inline bool
penetrance_model::penetrance_informative() const
{
  for(uint pid = 0; pid < phenotype_count()+1; ++pid)
  {
    if(penetrance_informative_phenotype(pid)) return true;
  }

  return false;
}

inline void
penetrance_model::alias_phenotype(const string& pname, const string& alias)
{
  uniquify();

  my_info->phmodel.alias_phenotype(pname, alias);
}

inline void
penetrance_model::alias_phenotype(uint pid, const string& alias)
{
  uniquify();

  my_info->phmodel.alias_phenotype(pid, alias);
}

inline void
penetrance_model::set_phenotype_strict(uint pid, bool strict)
{
  //lint --e{56} --e{48} --e{737} --e{63} <- Spurious messages

  uniquify();

  if(strict == my_info->strict_phenotypes[pid]) return;

  my_info->strict_phenotypes[pid] = strict;

  if(strict)
  {
    if(!codominant(pid, true)) ++my_info->codominant;
  }
  else
  {
    if(!codominant(pid, true)) --my_info->codominant;
  }  
}


//----------
//
inline double
penetrance_model::phased_penetrance(uint ptid, int gtid) const
{
    return my_info->phased_penetrance((int) ptid, gtid);
}

inline double
penetrance_model::unphased_penetrance(uint ptid, uint gtid) const
{
    return my_info->unphased_penetrance((int) ptid, (int) gtid);
}

inline double
penetrance_model::phased_penetrance(const phenotype& p, const phased_genotype& g) const
{
    return my_info->phased_penetrance((int) p.id(), g.get_id());
}

inline double
penetrance_model::unphased_penetrance(const phenotype& p, const phased_genotype& g) const
{
    return my_info->unphased_penetrance((int) p.id(), g.get_id());
}


inline bool
penetrance_model::is_x_linked() const
{
  return my_info->gmodel.is_x_linked();
}

inline bool
penetrance_model::is_y_linked() const
{
  return my_info->gmodel.is_y_linked();
}


inline bool
penetrance_model::is_autosomal() const
{
  return gmodel().is_autosomal();
}
inline GenotypeModelType
penetrance_model::get_model_type() const
{
  return gmodel().get_model_type();
}

inline void
penetrance_model::set_model_type(GenotypeModelType t)
{
  if(t == get_model_type()) return;
  
  bool gen_pheno = phmodel().has_generated_phenotypes();
  bool ext_pheno = phmodel().has_external_phenotypes();
  
  // External phenotypes point to genotypes which are likely to be
  // invalid now, and we can't create new phenotypes without likely screwing things up,
  // so call this a critical error
  if(ext_pheno)
    SAGE_internal_error();

  my_info->gmodel.set_model_type(t);

  // If we have generated phenotypes, we must create a new batch and reindex.
  if(gen_pheno)
  {
    penetrance_model m(gmodel(), name(), true, true);
    std::swap(*this,m);
  }
}

//----------
//
inline void
penetrance_model::mark_for_remap(const string& n)
{
  uniquify();
  my_info->gmodel.mark_for_remap(n);
}


inline void
penetrance_model::remap()
{
    uniquify();
    remap(*this, *this);
}


inline void
penetrance_model::remap(penetrance_model& dst) const
{
    remap(*this, dst);
}

//----------
//
//lint -e{1762} <- Member function not const
inline void
penetrance_model::reset_codominance(uint ptid, bool codom)
{
    //lint --e{56} --e{48} --e{737} --e{63} <- Spurious Messages

    if(codom != my_info->codominant_phenotypes[ptid])
    {
      my_info->codominant_phenotypes[ptid] = codom;

      if(codom) --my_info->strict_codominant;
      else      ++my_info->strict_codominant;

      if(my_info->strict_phenotypes[ptid])
        if(codom) --my_info->codominant;
        else      ++my_info->codominant;
    }
}

} // End namespace MLOCUS
} // End namespace SAGE

