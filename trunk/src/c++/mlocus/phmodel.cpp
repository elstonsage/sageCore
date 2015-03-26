//============================================================================
//  File:       phmodel.cpp
//
//  Author:     Bob Steagall
//
//  History:    Version 1.00.00
//
//  Notes:    
//
//  Copyright (c) 1999 R.C. Elston
//  All Rights Reserved
//============================================================================
//
#ifdef _MSC_VER
    #include <app/SAGEconfig.h>
    #pragma hdrstop
#endif

#include "mlocus/phmodel.h"
 
namespace SAGE   {
namespace MLOCUS {

#ifdef _DEBUG
    static  const char  REMAP[]   = "~";
    static  const char  MISSING[] = "*";
#else
    static  const char  MISSING[] = "*missing";
#endif

const string    nulls;

string  trim(const string& src);

//============================================================================
//  IMPLEMENTATION: phenotype_model_info
//============================================================================
//
phenotype_model_info::phenotype_model_info()
  : name(), missing_ptname(MISSING),
    missing_ptid(0),
    my_has_generated_phenotypes(false),
    my_has_external_phenotypes(false)
{}


phenotype_model_info::phenotype_model_info
(const string& n)
  : name(n), missing_ptname(MISSING),
    missing_ptid(0),
    my_has_generated_phenotypes(false),
    my_has_external_phenotypes(false)
{}



//============================================================================
//  IMPLEMENTATION: penetrance_model support classes
//============================================================================
//

/*
struct phenotype_sorter
{
    bool    operator ()(const phenotype& pa, const phenotype& pb) const
            {
                return pa.name() < pb.name();
            }
};
*/

//============================================================================
//  IMPLEMENTATION: phenotype_model
//============================================================================
//
phenotype_model::phenotype_model()
  : my_info(new phenotype_model_info())
{
  add_missing_phenotype();
}


phenotype_model::phenotype_model(const phenotype_model& pm)
  : my_info(pm.my_info)
{ }


phenotype_model::phenotype_model(const genotype_model& gm)
  : my_info(new phenotype_model_info())
{
    build_phenotypes_symmetric(gm);
}


phenotype_model::phenotype_model(const genotype_model& gm, const string& n)
  : my_info(new phenotype_model_info(n))
{
    build_phenotypes_symmetric(gm);
}


phenotype_model::phenotype_model(phenotype_model_info* info)
  : my_info(info)
{}


phenotype_model::~phenotype_model()
{}


phenotype_model&
phenotype_model::operator =(const phenotype_model& pm)
{
    if (&pm != this  &&  pm.my_info != my_info)
    {
        my_info = pm.my_info;
    }
    return *this;
}


//----------
//
phenotype
phenotype_model::get_phenotype(const string& nm) const
{
    if (nm == my_info->missing_ptname)
    {
        return my_info->phenotypes[my_info->missing_ptid];
    }
    else
    {
        phenotype_map::const_iterator i = my_info->phenotype_names.find(nm);

        if (i == my_info->phenotype_names.end()) return phenotype();
        
        //lint -e{115}
        return my_info->phenotypes[i->second];
    }
}


const phenotype&
phenotype_model::get_missing_phenotype() const
{
    return my_info->phenotypes[my_info->missing_ptid];
}

//----------
//
uint
phenotype_model::get_phenotype_id(const string& nm) const
{
    if (nm == my_info->missing_ptname)
    {
        return my_info->missing_ptid;
    }
    else
    {
        phenotype_map::const_iterator i = my_info->phenotype_names.find(nm);

        if (i == my_info->phenotype_names.end()) return NPOS;
        
        //lint -e{115}
        return i->second;
    }
}


//----------
//
phenotype_model
phenotype_model::clone() const
{
    phenotype_model    pm(*this);

    pm.uniquify();

    return pm;
}


//----------
//
void
phenotype_model::add_phenotype(const string& pstr)
{
    string  pname;

    pname = trim(pstr);

    if(pname.size() == 0) return;
    if(pname == my_info->missing_ptname) return;
    
    if(my_info->phenotype_names.find(pname) != my_info->phenotype_names.end()) return;

    uniquify();

    uint id = phenotype_count() + 1;

    my_info->phenotypes.push_back( phenotype(pname, id) );
    my_info->phenotype_names[pname] = id;
    my_info->my_has_external_phenotypes = true;
    
}

void
phenotype_model::alias_phenotype(const string& pstr, const string& alstr)
{
    uniquify();

    string  pname, alname;

    pname  = trim(pstr);
    alname = trim(alstr);

    phenotype_map::iterator p = my_info->phenotype_names.find(pname);
    phenotype_map::iterator a = my_info->phenotype_names.find(alname);
    phenotype_map::iterator alias_end = my_info->phenotype_names.end();

    if(p == alias_end || a != alias_end) return;
    
    my_info->phenotype_names[alname] = p->second;
}

void
phenotype_model::alias_phenotype(uint pid, const string& alias)
{
    uniquify();

    string  alname;

    alname = trim(alias);

    phenotype_map::iterator a = my_info->phenotype_names.find(alname);
    phenotype_map::iterator alias_end = my_info->phenotype_names.end();

    if(pid > phenotype_count() || a != alias_end) return;

    my_info->phenotype_names[alname] = pid;
}


void
phenotype_model::clear()
{
    uniquify();
    my_info->phenotypes.clear();
    my_info->phenotype_names.clear();
    my_info->my_has_generated_phenotypes  = false;
    my_info->my_has_external_phenotypes  = false;
}


void
phenotype_model::set_missing_phenotype_name(const string& nm)
{
    uniquify();

    // Remove the old mapping

    phenotype_map::iterator i = my_info->phenotype_names.find(my_info->missing_ptname);

    my_info->phenotype_names.erase(i);

    // Modify the new mappings.

    my_info->missing_ptname = nm;

    uint id = my_info->missing_ptid;

    my_info->phenotype_names[nm] = id;

    my_info->phenotypes[id] = phenotype(nm, id);
}

void
generate_unphased_sex_specific_aliases
    (phenotype_model& ph, const string& name,
     allele a1, allele a2, const genotype_model& gmodel)
{
  if(!a1.is_sex_allele() && !a2.is_sex_allele())
    return;
  
  if(a1 == a2) return;
  
  string missing = gmodel.missing_allele_name();
  
  if(a1.is_sex_allele())
  {
    ph.alias_phenotype(name, missing + gmodel.unphased_separator() + a2.name());
    ph.alias_phenotype(name, a2.name() + gmodel.unphased_separator() + missing);
    
    if(a2.is_null_y_allele())
      ph.alias_phenotype(name, a2.name() + gmodel.unphased_separator() + a2.name() + "(male)");
    else
      ph.alias_phenotype(name, a2.name() + gmodel.unphased_separator() + a2.name());      
  }
  if(a2.is_sex_allele())
  {
    ph.alias_phenotype(name, missing + gmodel.unphased_separator() + a1.name());
    ph.alias_phenotype(name, a1.name() + gmodel.unphased_separator() + missing);

    if(a2.is_null_y_allele())
      ph.alias_phenotype(name, a1.name() + gmodel.unphased_separator() + a1.name() + "(male)");
    else
      ph.alias_phenotype(name, a1.name() + gmodel.unphased_separator() + a1.name());
  }
}

void
generate_phased_sex_specific_aliases
    (phenotype_model& ph, const string& name,
     allele a1, allele a2, const genotype_model& gmodel)
{
  if(!a1.is_sex_allele() && !a2.is_sex_allele())
    return;
  
  if(a1 == a2) return;
  
  string missing = gmodel.missing_allele_name();
  
  if(a1.is_sex_allele())
  {
    ph.alias_phenotype(name, missing + gmodel.forward_separator() + a2.name());
    ph.alias_phenotype(name, a2.name() + gmodel.backward_separator() + missing);

    if(a1.is_null_y_allele())
    {
      ph.alias_phenotype(name, a2.name() + gmodel.forward_separator() + a2.name() + "(male)");
      ph.alias_phenotype(name, a2.name() + gmodel.backward_separator() + a2.name() + "(male)");
    }
    else
    {
      ph.alias_phenotype(name, a2.name() + gmodel.forward_separator() + a2.name());
      ph.alias_phenotype(name, a2.name() + gmodel.backward_separator() + a2.name());
    }
  }
  if(a2.is_sex_allele())
  {
    ph.alias_phenotype(name, missing + gmodel.backward_separator() + a1.name());
    ph.alias_phenotype(name, a1.name() + gmodel.forward_separator() + missing);
    
    if(a2.is_null_y_allele())
    {
      ph.alias_phenotype(name, a1.name() + gmodel.forward_separator() + a1.name() + "(male)");
      ph.alias_phenotype(name, a1.name() + gmodel.backward_separator() + a1.name() + "(male)");
    }
    else
    {
      ph.alias_phenotype(name, a1.name() + gmodel.forward_separator() + a1.name());
      ph.alias_phenotype(name, a1.name() + gmodel.backward_separator() + a1.name());
    }
  }
}

//lint -e{1762}
void
phenotype_model::push_phenotype(const string& nm, const genotype_model& gmodel)
{
    string name2, alname1, alname2;

    Ordering order;

    uint ptid = get_phenotype_id(nm);

    if(ptid == NPOS)
    {
        ptid = phenotype_count() + 1;

        my_info->phenotypes.push_back( phenotype(nm, ptid) );

        my_info->phenotype_names[nm] = ptid;

        // Create the alias's

        //lint -e{534}
        gmodel.parse_genotype_name(nm, alname1, alname2, order);

        if(order == Unphased)
        {
          name2 = alname2 + gmodel.unphased_separator() + alname1;

          my_info->phenotype_names[name2] = ptid;

          if(!gmodel.is_autosomal())
            generate_unphased_sex_specific_aliases(*this, name2,
                                                   gmodel.get_allele(alname1),
                                                   gmodel.get_allele(alname2), gmodel);
        }
        else
        {
          // Create the alternative phenotype name

          name2 = alname2 + gmodel.backward_separator() + alname1;

          my_info->phenotype_names[name2] = ptid;
          if(!gmodel.is_autosomal())
            generate_phased_sex_specific_aliases(*this, name2,
                                                 gmodel.get_allele(alname1),
                                                 gmodel.get_allele(alname2), gmodel);
        }
    }
}

void
phenotype_model::add_genotypes(const genotype_model& gmodel)
{
    // Now all the rest of them.

    string  pname;

    string  missing_allele_name    = gmodel.missing_allele_name();

    my_info->my_has_generated_phenotypes = true;

    //- Add a phenotype for each genotype name, both phased and unphased.
    //
    for(unphased_genotype_iterator upi = gmodel.unphased_genotype_begin();
        upi != gmodel.unphased_genotype_end(); ++upi)
    {
        pname = upi->name();

        push_phenotype(pname, gmodel);
    }

    phased_genotype_iterator    pgf = gmodel.phased_genotype_begin();       
    phased_genotype_iterator    pgl = gmodel.phased_genotype_end();     

    for ( ;  pgf != pgl;  ++pgf)
    {
        pname = (*pgf).name();

        push_phenotype(pname, gmodel);
    }
}

void
phenotype_model::add_autosomal_phenotypes_by_allele(const string& nm, const genotype_model& gmodel)
{
    string  pname1, pname2, pname3;

    string  missing_allele_name = gmodel.missing_allele_name();

    char    unphased_separator  = gmodel.unphased_separator();
    char    forward_separator   = gmodel.forward_separator();

    //- Add a phenotype for each genotype name, both phased and unphased.
    //
    for (allele_iterator  a = gmodel.allele_begin();  a != gmodel.allele_end();
         ++a)
    {
        pname1 = a->name() + unphased_separator + nm;
        push_phenotype(pname1, gmodel);
    }

    for (allele_iterator  a = gmodel.allele_begin();  a != gmodel.allele_end();
         ++a)
    {
        pname2 = a->name() + forward_separator + nm;
        pname3 = nm  + forward_separator + a->name();

        push_phenotype(pname2, gmodel);
        push_phenotype(pname3, gmodel);
    }
}

void
phenotype_model::add_female_x_linked_phenotypes_by_allele
    (const string& nm, const genotype_model& gmodel)
{
    add_autosomal_phenotypes_by_allele(nm, gmodel);
}

void
phenotype_model::add_male_x_linked_phenotypes_by_allele
    (const string& nm, const genotype_model& gmodel)
{
    char    unphased_separator  = gmodel.unphased_separator();
    char    forward_separator   = gmodel.forward_separator();

    string y_allele = gmodel.get_sex_specific_allele().name();
    string missing  = gmodel.missing_allele_name();
            
    string uph_sexed = nm + unphased_separator + y_allele;
    string pph_sexed = nm + forward_separator + y_allele;
            
    push_phenotype(uph_sexed, gmodel);
    push_phenotype(pph_sexed, gmodel);
}

void
phenotype_model::add_female_y_linked_phenotypes_by_allele
    (const string& nm, const genotype_model& gmodel)
{
  // Does nothing
}

void
phenotype_model::add_male_y_linked_phenotypes_by_allele
    (const string& nm, const genotype_model& gmodel)
{
    char    unphased_separator  = gmodel.unphased_separator();
    char    forward_separator   = gmodel.forward_separator();

    allele x_allele = gmodel.get_sex_specific_allele();
    
    string uph_sexed = x_allele.name() + unphased_separator + nm;
    string pph_sexed = x_allele.name() + forward_separator + nm;
    
    push_phenotype(uph_sexed, gmodel);
    push_phenotype(pph_sexed, gmodel);
}

void
phenotype_model::add_allele(const string& nm, const genotype_model& gmodel)
{
    switch(gmodel.get_model_type())
    {
      case AUTOSOMAL :
          add_autosomal_phenotypes_by_allele(nm, gmodel);
          break;
          
      case X_LINKED  :
          add_female_x_linked_phenotypes_by_allele(nm, gmodel);
          add_male_x_linked_phenotypes_by_allele(nm,gmodel);
          break;
      case Y_LINKED  :
          add_female_y_linked_phenotypes_by_allele(nm, gmodel);
          add_male_y_linked_phenotypes_by_allele(nm,gmodel);
          break;
    }

    my_info->my_has_generated_phenotypes = true;
}


//============================================================================
//============================================================================
//
void
phenotype_model::build_phenotypes_symmetric(const genotype_model& gmodel)
{
    add_missing_phenotype();
    
    add_missing_phenotype_aliases(gmodel);
  
    add_genotypes(gmodel);
}

void
phenotype_model::add_missing_phenotype()
{
    // First, add the missing phenotype

    my_info->phenotypes.push_back( phenotype(my_info->missing_ptname, 0) );
    my_info->phenotype_names[my_info->missing_ptname] = 0;

    my_info->missing_ptid = 0;
}

void 
phenotype_model::add_missing_phenotype_aliases(const genotype_model& gm)
{
  alias_phenotype(missing_phenotype_name(), gm.missing_allele_name() + gm.unphased_separator() + gm.missing_allele_name());
  alias_phenotype(missing_phenotype_name(), gm.missing_allele_name() + gm.forward_separator() + gm.missing_allele_name());
  alias_phenotype(missing_phenotype_name(), gm.missing_allele_name() + gm.backward_separator() + gm.missing_allele_name());
}

void
phenotype_model::uniquify()
{
    if (!my_info.unique())
    {
        my_info = boost::shared_ptr<phenotype_model_info>
                     (new phenotype_model_info(*my_info));
    }
}


} // End namespace MLOCUS
} // End namespace SAGE

