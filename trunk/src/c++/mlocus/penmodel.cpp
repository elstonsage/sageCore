//============================================================================
//  File:       penmodel.cpp
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

#include <iostream>
#include "mlocus/penmodel.h"
 
namespace SAGE   {
namespace MLOCUS {

#ifdef _DEBUG
    static  const char  REMAP[]   = "~";
    static  const char  MISSING[] = "*";
#else
//    static  const char  REMAP[]   = "~remap";
//    static  const char  MISSING[] = "*missing";
#endif

const string    nulls;
const double    PNVAL = -1.0;

string  trim(const string& src);

//============================================================================
//  IMPLEMENTATION: penetrance_info
//============================================================================
//
penetrance_info::penetrance_info()
  : pvalue(PNVAL), pname(), a1name(), a2name(), order(Unphased)
{}


penetrance_info::penetrance_info
(double pval, const string& p, const string& a1, const string& a2, Ordering o)
  : pvalue(pval), pname(p), a1name(a1), a2name(a2), order(o)
{}


//============================================================================
//  IMPLEMENTATION: penetrance_model
//============================================================================
//
penetrance_model::penetrance_model()
  : my_info(new penetrance_model_info())
{
  resize_matrices();
  create_missing_phenotype_penetrance(get_missing_phenotype_id());
}


penetrance_model::penetrance_model(const penetrance_model& pm)
  : my_info(pm.my_info)
{}


penetrance_model::penetrance_model
    (const genotype_model& gm, bool ph, bool pen)
  : my_info(new penetrance_model_info())
{
  my_info->gmodel = gm;

  if(ph)
  {
    my_info->phmodel = phenotype_model(gm);

    if(pen) build_symmetric_penetrances();
  }
  resize_matrices();

  create_missing_phenotype_penetrance(get_missing_phenotype_id());
}

penetrance_model::penetrance_model
    (const genotype_model& gm, const phenotype_model& pm, bool pen)
  : my_info(new penetrance_model_info(gm, pm))
{
  if(pen) build_symmetric_penetrances();

  resize_matrices();

  create_missing_phenotype_penetrance(get_missing_phenotype_id());
}


penetrance_model::penetrance_model
    (const genotype_model& gm, const string& n, bool ph, bool pen)
  : my_info(new penetrance_model_info())
{
  my_info->gmodel = gm;
  my_info->name   = n;

  if(ph)
  {
    my_info->phmodel = phenotype_model(gm);

    if(pen) build_symmetric_penetrances();
  }

  resize_matrices();

  create_missing_phenotype_penetrance(get_missing_phenotype_id());
}

penetrance_model::penetrance_model
    (const genotype_model& gm, const phenotype_model& pm, const string& n, bool pen)
  : my_info(new penetrance_model_info(gm, pm, n))
{
  resize_matrices();

  if(pen) build_symmetric_penetrances();

  create_missing_phenotype_penetrance(get_missing_phenotype_id());
}


penetrance_model::penetrance_model(penetrance_model_info* info)
  : my_info(info)
{}


penetrance_model::~penetrance_model()
{}


penetrance_model&
penetrance_model::operator =(const penetrance_model& pm)
{
    if (&pm != this  &&  pm.my_info != my_info)
    {
        my_info = pm.my_info;
    }
    return *this;
}


//----------
//
//----------
//
penetrance_model
penetrance_model::clone() const
{
    penetrance_model    pm(*this);

    pm.uniquify();

    return pm;
}

//----------
//
void
penetrance_model::add_allele(const string& n, double freq,
                             bool ph, bool pen)
{
  uint count = allele_count();

  my_info->gmodel.add_allele(n, freq);

  if(count != allele_count())
  {
      uniquify();

      //lint --e{534}
      resize_matrices();

      create_missing_phenotype_penetrance(get_missing_phenotype_id());

      if(ph)
      {
          // We can call the incremental version of add_allele only when no
          // other alleles were added without a previous re-build.
          if(count == my_info->last_alleles)
            my_info->phmodel.add_allele(n, my_info->gmodel);
          else
            my_info->phmodel.add_genotypes(my_info->gmodel);

          my_info->last_alleles = allele_count();

          resize_matrices();
          
          // NOTE:  this is expensive.  Should be fixed later.
          if(pen) build_symmetric_penetrances();
      }
  }
}

/// \internal
///
/// Internal Utility function used by penetrance model.  Not exported.
///
/// Given an allele string, adds the allele to the marker if it isn't already
/// present with an allele frequency of 0.0.
inline bool add_allele_to_model(penetrance_model& marker_info, const std::string& allele)
{
  const MLOCUS::genotype_model& gmodel = marker_info.gmodel();
  
  if( allele.size()                          &&
      allele != gmodel.missing_allele_name() && 
      !gmodel.get_allele(allele).is_valid()     )
  {
    marker_info.add_allele(allele, 0.0, true, true);
    return true;
  }
  return false;
}

/// \internal
///
/// Internal Utility function used by penetrance model.  Not exported.
///
/// Returns \c true if the marker allows expansion alleles, \c false otherwise.
///
inline bool marker_is_expandable(penetrance_model& marker_info)
{
  return marker_info.gmodel().dynamic_alleles() && marker_info.codominant();
}

/// If the marker is dynamic (allows new alleles), adds alleles which may be
/// missing.  Returns \c true if any new alleles were added, \c false otherwise.
///
/// Since this function will only operate if dynamic behavior is already specified,
/// phenotypes and penetrances are always added.
///
/// Sex alleles are dealt with appropriately, but ordering of alleles is unimportant
/// to this function.  Specifically, 'backwards' sex-allele based genotypes ("~Y/A" for
/// example) are still valid, and the non-sex-specific alleles are still added to the
/// model, but only valid phenotypes/penetrances are added.
bool penetrance_model::add_genotype_dynamically(const std::string& genotype)
{
    bool added1 = false;
    bool added2 = false;
  
    if(marker_is_expandable(*this))
    {
      string   allele1;
      string   allele2;
      MLOCUS::Ordering order;
 
      gmodel().parse_genotype_name(genotype, allele1, allele2, order);

      added1 = add_allele_to_model(*this, allele1);
      added2 = add_allele_to_model(*this, allele2);
    }
    
    return added1 || added2;
}

uint
penetrance_model::add_phenotype(const string& pname, bool strict)
{
  uint pid = get_phenotype_id(pname);

  // Phenotype already exists
  if( pid != NPOS )
  {
    // Only make unique and modify if strictness changes
    if( strict_phenotype(pid) != strict )
    {
      uniquify();
      set_phenotype_strict(pid, strict);
    }
    return pid;
  }

  // Otherwise add the new phenotype

  uniquify();

  my_info->phmodel.add_phenotype(pname);

  resize_matrices();

  pid = get_phenotype_id(pname);

  set_phenotype_strict(pid, strict);

  return pid;
}

void
penetrance_model::copy_penetrance(uint psource, uint pdest, bool override)
{
  copy_penetrance_sexed(psource, pdest, MPED::SEX_MISSING, override);
}

void
penetrance_model::copy_penetrance
    (const penetrance_model& model, uint psource, uint pdest, bool override)
{
  copy_penetrance_sexed(model, psource, pdest, MPED::SEX_MISSING, override);
}



void
penetrance_model::copy_penetrance_sexed(uint psource, uint pdest, MPED::SexCode sex, bool override)
{
    //lint --e{534}

    if(psource == pdest) return;

    uniquify();

    if(override)
    {
        my_info->unphased_penetrance.clear_row((int) pdest);
        my_info->phased_penetrance  .clear_row((int) pdest);
    }

    unphased_penetrance_iterator bu = unphased_penetrance_begin(psource);
    unphased_penetrance_iterator eu = unphased_penetrance_end  (psource);

    for( ; bu != eu; ++bu)
    {
        if(MPED::is_sex_unknown(sex) ||
           (MPED::is_male(sex)   && bu.unphased_geno().is_male_compatible()) ||
           (MPED::is_female(sex) && bu.unphased_geno().is_female_compatible())  )
          my_info->unphased_penetrance.set((int) pdest, (int) bu.geno_id(), *bu);
    }

    phased_penetrance_iterator bp = phased_penetrance_begin(psource);
    phased_penetrance_iterator ep = phased_penetrance_end  (psource);

    for( ; bp != ep; ++bp)
    {
        if(MPED::is_sex_unknown(sex) ||
           (MPED::is_male(sex)   && bp.phased_geno().is_male_compatible()) ||
           (MPED::is_female(sex) && bp.phased_geno().is_female_compatible())  )
          my_info->phased_penetrance.set((int) pdest, bp.geno_id(), *bp);
    }

    check_for_codominance(pdest);
}

void
penetrance_model::copy_penetrance_sexed
    (const penetrance_model& model, uint psource, uint pdest, MPED::SexCode sex, bool override)
{
    //lint --e{534}

    if(psource == pdest && model.my_info == my_info) return;

    uniquify();

    if(override)
    {
        my_info->unphased_penetrance.clear_row((int) pdest);
        my_info->phased_penetrance  .clear_row((int) pdest);
    }

    unphased_penetrance_iterator bu = model.unphased_penetrance_begin(psource);
    unphased_penetrance_iterator eu = model.unphased_penetrance_end  (psource);

    for( ; bu != eu; ++bu)
    {
        if(MPED::is_sex_unknown(sex) ||
           (MPED::is_male(sex)   && bu.unphased_geno().is_male_compatible()) ||
           (MPED::is_female(sex) && bu.unphased_geno().is_female_compatible())  )
          my_info->unphased_penetrance.set((int) pdest, (int) bu.geno_id(), (*bu));
    }

    phased_penetrance_iterator bp = model.phased_penetrance_begin(psource);
    phased_penetrance_iterator ep = model.phased_penetrance_end  (psource);

    for( ; bp != ep; ++bp)
    {
        if(MPED::is_sex_unknown(sex) ||
           (MPED::is_male(sex)   && bp.phased_geno().is_male_compatible()) ||
           (MPED::is_female(sex) && bp.phased_geno().is_female_compatible())  )
          my_info->phased_penetrance.set((int) pdest, bp.geno_id(), (*bp));
    }

    check_for_codominance(pdest);
}

void
penetrance_model::clear_penetrance(uint pdest)
{
    //lint --e{534}

    uniquify();

    my_info->unphased_penetrance.clear_row((int) pdest);
    my_info->phased_penetrance.clear_row  ((int) pdest);

    check_for_codominance(pdest);
}

//----------
//

void
penetrance_model::add_penetrance(double pval, const string& pname, const string& gname)
{
    uniquify();
    add_penetrance(pval, pname, gname, Unphased, '\0');
}


void
penetrance_model::add_penetrance
(double pval, const string& pname, const string& gname, Ordering order)
{
    uniquify();

    if (order == Unphased)
    {
        add_penetrance(pval, pname, gname, order, my_info->gmodel.unphased_separator());
    }
    else if (order == PhasedForward)
    {
        add_penetrance(pval, pname, gname, order, my_info->gmodel.forward_separator());
    }
    else if (order == PhasedBackward)
    {
        add_penetrance(pval, pname, gname, order, my_info->gmodel.backward_separator());
    }
}


void
penetrance_model::add_penetrance
(double pval, const string& pstr, const string& gstr, Ordering order, char sep)
{
    uniquify();

    string  aname1;
    string  aname2;
    string  gname;
    string  pname;

    uint ptid;

    gname = my_info->gmodel.parse_genotype_name(gstr, aname1, aname2, order, sep);
    
    pname = trim(pstr);
    
    ptid = get_phenotype_id(pname);

    if(ptid == get_missing_phenotype_id()) return;

    // Sort the alleles for potential problems
    if(aname1 == missing_allele_name())
    {
       aname1.swap(aname2);

       if(order != Unphased) order = (Ordering) (3 - order);
    }
    
    if(order == Unphased)
    {
      unphased_genotype gtype = get_unphased_genotype(gname);
      
      if(gtype.is_valid())
      {
          add_unphased_penetrance(pval, ptid, gtype.get_id());
          return;
      }

      if(aname1 == missing_allele_name())
      {
          create_missing_phenotype_penetrance(ptid, pval);
      }
      else
      {
          if(aname2 == missing_allele_name())
          {
              allele a = get_allele(aname1);

              add_penetrance_with_missing_allele(ptid, a, Unphased, pval);
          }
      }

      //lint -e{534}
      check_for_codominance(ptid);

    }
    else
    {
      unphased_genotype gtype = get_unphased_genotype(gname);
      
      if(gtype.is_valid())
      {
          add_phased_penetrance(pval, ptid, gtype.get_id());

          return;
      }

      if(aname1 == missing_allele_name())
      {
          create_missing_phenotype_penetrance(ptid, pval);
      }
      else
      {
          if(aname2 == missing_allele_name())
          {
              allele a = get_allele(aname1);

              add_penetrance_with_missing_allele(ptid, a, order, pval);
          }
      }

      //lint --e{534}
      check_for_codominance(ptid);

    }
}

void
penetrance_model::add_phased_penetrance(double pval, uint pid, int gid)
{
  //lint --e{534}

  uniquify();

  my_info->phased_penetrance.set((int) pid, gid, (pval));

  check_for_codominance(pid);
}

void
penetrance_model::add_unphased_penetrance(double pval, uint pid, uint gid)
{
  //lint --e{534}

  uniquify();

  my_info->unphased_penetrance.set((int) pid, (int) gid, (pval));

  check_for_codominance(pid);
}

void
penetrance_model::remove_phased_penetrance(uint pid, int gid, bool c)
{
  //lint --e{534}

  uniquify();
  
  my_info->phased_penetrance.remove((int) pid, gid);
  
  phased_genotype pg = get_phased_genotype(gid);
  
  uint flipped_gid = pg.get_flipped_phased_genotype().get_id();

  if(c && phased_penetrance(pid, flipped_gid) == 0)
    my_info->unphased_penetrance.remove((int) pid, pg.get_equivalent_unphased_genotype().get_id());

  check_for_codominance(pid);
}

void
penetrance_model::remove_unphased_penetrance(uint pid, uint gid, bool c)
{
  //lint --e{534}

  uniquify();

  my_info->unphased_penetrance.remove((int) pid, (int) gid);

  if(c)
  {
    unphased_genotype ug = get_unphased_genotype(gid);
    
    phased_genotype g1 = ug.get_equivalent_phased_genotype1();
    phased_genotype g2 = ug.get_equivalent_phased_genotype2();

    my_info->phased_penetrance.remove((int) pid, g1.get_id());
    my_info->phased_penetrance.remove((int) pid, g2.get_id());
  }

  check_for_codominance(pid);
}


void
penetrance_model::set_missing_phenotype_name(const string& n)
{
    uniquify();
    my_info->phmodel.set_missing_phenotype_name(n);
}


//============================================================================
//============================================================================
//
void
penetrance_model::build_symmetric_penetrances()
{
    //lint --e{534}

    if(!my_info->phmodel.has_generated_phenotypes()) return;

    int    pid;

    string  maname = missing_allele_name();

    resize_matrices();
    
    unphased_genotype_iterator    ugf = unphased_genotype_begin();       
    unphased_genotype_iterator    ugl = unphased_genotype_end();     

    for( ; ugf != ugl; ++ugf)
    {
        pid = (int) get_phenotype_id(ugf->name());

        // Set for unphased genotype
        my_info->unphased_penetrance.set(pid, ugf->get_id(), (1.0));

        // Set for phased genotypes
        phased_genotype g1 = ugf->get_equivalent_phased_genotype1();
        phased_genotype g2 = ugf->get_equivalent_phased_genotype2();
        
        my_info->phased_penetrance.set(pid, g1.get_id(), (1.0));
        my_info->phased_penetrance.set(pid, g2.get_id(), (1.0));
    }

    phased_genotype_iterator    pgf = phased_genotype_begin();       
    phased_genotype_iterator    pgl = phased_genotype_end();     

    for ( ;  pgf != pgl;  ++pgf)
    {
        pid = (int) get_phenotype_id(pgf->name());

        // Set for phased genotype
        my_info->phased_penetrance.set(pid, pgf->get_id(), (1.0));

        // Set for unphased genotype
        unphased_genotype upg = pgf->get_equivalent_unphased_genotype();
      
        my_info->unphased_penetrance.set(pid, upg.get_id(), (1.0));
    }
}

//lint -e{1762}
void penetrance_model::add_penetrance_with_missing_allele
    (uint pid, allele al, Ordering ordering, double val)
{
    if(!al.is_valid()) return;

    char    unphased_separator = my_info->gmodel.unphased_separator();
    char    forward_separator  = my_info->gmodel.forward_separator();

    string  aname, anameA, anameB;

    aname = al.name();

    for(allele_iterator al2 = allele_begin(); al2 != allele_end(); ++al2)
    {
        anameA = aname;
        anameB = al2->name();

        if(ordering == PhasedBackward) anameA.swap(anameB);

        if(ordering == Unphased)
        {
          unphased_genotype gtype = get_unphased_genotype(anameA + unphased_separator + anameB);

          //lint -e{534}

          my_info->unphased_penetrance.set((int) pid, (int) gtype.get_id(), val);
        }
        else
        {
          phased_genotype gtype = get_phased_genotype(anameA + forward_separator + anameB);

          //lint -e{534}

          my_info->phased_penetrance.set((int) pid, gtype.get_id(), val);
        }
    }
}

void
penetrance_model::make_consistent()
{
    uniquify();

    int    ptid;

    phased_genotype    g1, g2;

    double  pv, pv1, pv2;

    bool pvs, pvs1, pvs2;

    //- Save various sizes in local variables.
    //
    uint    pcnt  = phenotype_count() + 1;

    //- Enforce consistency of penetrance matrix data.
    //
    for (ptid = 1;  ptid < (int) pcnt;  ++ptid)
    {
        unphased_genotype_iterator    ugf = unphased_genotype_begin();       
        unphased_genotype_iterator    ugl = unphased_genotype_end();     

        for( ; ugf != ugl; ++ugf)
        {
            //- Get the ids of the two relevant phased genotypes.
            //
            g1 = ugf->get_equivalent_phased_genotype1();
            g2 = ugf->get_equivalent_phased_genotype2();

            //- Next, check to see if corresponding phased and unphased
            //  penetrance values are not the default
            //
            pvs  = my_info->unphased_penetrance.set(ptid,ugf->get_id());
            pvs1 = my_info->phased_penetrance  .set(ptid,g1.get_id());
            pvs2 = my_info->phased_penetrance  .set(ptid,g2.get_id());

            //- There are two cases to consider.  First, if the unphased
            //  value is set, but the phased values are not, then the
            //  phased values must be set equal to the unphased value.
            //
            if (pvs)
            {
              if(!pvs1 && !pvs2)
              {
                pv  = my_info->unphased_penetrance(ptid, ugf->get_id());

                //lint --e{534} 

                my_info->phased_penetrance.set(ptid, g1.get_id(), (pv));
                my_info->phased_penetrance.set(ptid, g2.get_id(), (pv));
              }
            }

            //- Otherwise, if the unphased value is unset, but at least one
            //  of the two phased values is set, then the unphased value
            //  must be assigned the average of the two phased values.
            //
            //  Note, this is probably a bad idea, but is the best we can
            //  do, given the situation.  It should not occur most of the
            //  time.
            else
            {
                if (pvs1 ||  pvs2)
                {

                    pv1 = my_info->phased_penetrance(ptid, g1.get_id());
                    pv2 = my_info->phased_penetrance(ptid, g2.get_id());

                    pv = (pv1 + pv2)/2.0;

                    //lint -e{534}

                    my_info->unphased_penetrance.set(ptid, ugf->get_id(), (pv));
                }
            }
        }

        //lint -e{534}

        check_for_codominance((uint) ptid);
    }
}

/// Returns true if the phenotype is codominant ignoring phased genotypes
bool is_unphased_codominant_at_phenotype
    (const penetrance_model_info& pm, uint phenotype_id)
{
  // Is codominant if there are no unphased genotypes with non-zero values or
  // exactly one
  return pm.unphased_penetrance.row_elements(phenotype_id) < 2;
}

/// Returns true if the phenotype is codominant ignoring unphased genotypes
bool is_phased_codominant_at_phenotype
    (const penetrance_model_info& pm, uint phenotype_id)
{
  // We know it's codominant if the number of phased penetrances are
  // 0 or 1, and it's not if there's more than 2...
  if(pm.phased_penetrance.row_elements(phenotype_id) < 2) return true;
  if(pm.phased_penetrance.row_elements(phenotype_id) > 2) return false;
  
  // ..., 2 is a special case, though, and is ok only if the two phased genotypes
  // are alternatives of one another.
  
  phased_genotype p1 = pm.gmodel.get_phased_genotype(pm.phased_penetrance.row_begin(phenotype_id)->first.col);
  phased_genotype p2 = pm.gmodel.get_phased_genotype((++pm.phased_penetrance.row_begin(phenotype_id))->first.col);

  return (p1.get_flipped_phased_genotype() == p2);
}

bool is_codominant_at_phenotype
    (const penetrance_model_info& pm, uint phenotype_id)
{
    if (phenotype_id == pm.phmodel.get_missing_phenotype_id()) return true;

    bool is_unphased_codominant = is_unphased_codominant_at_phenotype (pm, phenotype_id);
    bool is_phased_codominant   = is_phased_codominant_at_phenotype   (pm, phenotype_id);

    if(!is_unphased_codominant || !is_phased_codominant) return false;
    
    // If either doesn't have any elements, we know we're ok.
    if(pm.phased_penetrance.row_elements(phenotype_id) == 0 ||
       pm.unphased_penetrance.row_elements(phenotype_id) == 0) return true;
    
    phased_genotype    pg = pm.gmodel.get_phased_genotype(pm.phased_penetrance.row_begin(phenotype_id)->first.col);
    unphased_genotype upg = pm.gmodel.get_unphased_genotype(pm.unphased_penetrance.row_begin(phenotype_id)->first.col);
    
    return pg.equivalent(upg);
}

bool
penetrance_model::check_for_codominance(uint ptid)
{
    bool codom = is_codominant_at_phenotype(*my_info, ptid);

    reset_codominance(ptid, codom);

    return codom;
}

void
penetrance_model::create_missing_phenotype_penetrance(uint pid, double val)
{
    //lint --e{534}

    uint id;

    //- Add penetrance values for all unphased genotypes
    //
    for (id = 0;  id < unphased_genotype_count();  ++id)
    {
        my_info->unphased_penetrance.set((int) pid, (int) id, val);
    }

    //- Add penetrance values for all phased genotypes
    //
    phased_genotype_iterator    pgf = my_info->gmodel.phased_genotype_begin();       
    phased_genotype_iterator    pgl = my_info->gmodel.phased_genotype_end();     

    for ( ;  pgf != pgl;  ++pgf)
    {
        my_info->phased_penetrance.set((int) pid, (*pgf).get_id(), val);
    }

    if(pid != get_missing_phenotype_id())
    {
      reset_codominance(pid, (allele_count() < 2));
    }
}

void
penetrance_model::resize_matrices()
{
    int    ptcnt = (int) phenotype_count() + 1;

    my_info->phased_penetrance.resize(ptcnt);
    my_info->unphased_penetrance.resize(ptcnt);

    //lint --e{1013} <- Spurious

    my_info->codominant_phenotypes.resize(ptcnt, true);
    my_info->strict_phenotypes.resize(ptcnt, true);
}

void
penetrance_model::clear()
{
    //lint --e{534}

    if(!my_info.unique())
    {
      my_info = boost::shared_ptr<penetrance_model_info>
                      (new penetrance_model_info());
    }
    else
    {
      uniquify();

      my_info->gmodel.clear();
      my_info->phmodel.clear();
      my_info->phased_penetrance.clear();
      my_info->unphased_penetrance.clear();
      my_info->codominant = false;
    }
}


void
penetrance_model::uniquify()
{
    if (!my_info.unique())
    {
        my_info = boost::shared_ptr<penetrance_model_info>
                      (new penetrance_model_info(*my_info));
        my_info->gmodel.uniquify();
        my_info->phmodel.uniquify();
    }
}


//============================================================================
//  IMPLEMENTATION: phased_penetrance_iterator
//============================================================================
//
penetrance_model::phased_penetrance_iterator::phased_penetrance_iterator
(const host_type& m, uint row, bool end)
  : my_host(&m), my_row(row), my_element()
{
    if (end)
    {
        my_element = my_host->my_info->phased_penetrance.row_end((int) my_row);
    }
    else
    {
        my_element = my_host->my_info->phased_penetrance.row_begin((int) my_row);
    }
}


//============================================================================
//  IMPLEMENTATION: unphased_penetrance_iterator
//============================================================================
//
penetrance_model::unphased_penetrance_iterator::unphased_penetrance_iterator
(const host_type& m, uint row, bool end)
  : my_host(&m), my_row(row), my_element()
{
    if (end)
    {
        my_element = my_host->my_info->unphased_penetrance.row_end((int) my_row);
    }
    else
    {
        my_element = my_host->my_info->unphased_penetrance.row_begin((int) my_row);
    }
}

void
penetrance_model::remap(const penetrance_model& src, penetrance_model& dst)
{
    genotype_model gtmp;

    src.my_info->gmodel.remap(gtmp);

    penetrance_model tmp(gtmp, src.my_info->phmodel);

    std::map<allele,allele> allele_map;

    allele remap_allele = gtmp.get_allele("~remap");

    for(allele_iterator i = src.allele_begin(); i != src.allele_end(); ++i)
    {
        allele dst_allele = gtmp.get_allele(i->name());
        
        if(dst_allele.is_valid()) allele_map[*i] = dst_allele;
        else                      allele_map[*i] = remap_allele;
    }
    
    if(!src.gmodel().is_autosomal())
      allele_map[src.gmodel().get_sex_specific_allele()] = gtmp.get_sex_specific_allele();
    
    unphased_penetrance_matrix& upm = tmp.my_info->unphased_penetrance;

    unphased_penetrance_matrix::row_iterator upb = 
        src.my_info->unphased_penetrance.row_begin(1);
    unphased_penetrance_matrix::row_iterator upe = 
        src.my_info->unphased_penetrance.row_end((int) src.phenotype_count());

    for( ; upb != upe; ++upb)
    {
        unphased_genotype oldg = src.get_unphased_genotype((uint) upb->first.col);

        allele a1 = allele_map[oldg.allele1()];
        allele a2 = allele_map[oldg.allele2()];
        
        unphased_genotype newg = tmp.my_info->gmodel.get_unphased_genotype(a1, a2);

        //lint -e{534}
        upm.set(upb->first.row, (int) newg.get_id(), upb->second);
    }

    phased_penetrance_matrix& ppm = tmp.my_info->phased_penetrance;

    phased_penetrance_matrix::row_iterator ppb = 
        src.my_info->phased_penetrance.row_begin(1);
    phased_penetrance_matrix::row_iterator ppe =
        src.my_info->phased_penetrance.row_end((int) src.phenotype_count());

    for( ; ppb != ppe; ++ppb)
    {
        phased_genotype oldg = src.get_phased_genotype(ppb->first.col);

        allele a1 = allele_map[oldg.allele1()];
        allele a2 = allele_map[oldg.allele2()];

        phased_genotype newg = tmp.my_info->gmodel.get_phased_genotype(a1, a2);

        //lint -e{534}
        ppm.set(ppb->first.row, newg.get_id(), ppb->second);
    }

    // Get all the phenotype flags set for the new model

    phenotype_iterator phen = tmp.phenotype_begin();
    phenotype_iterator pend = tmp.phenotype_end();

    for( ; phen != pend; ++phen)
    {
      bool b = src.strict_phenotype((*phen).id());

      tmp.set_phenotype_strict((*phen).id(), b);

      //lint -e{534}
      tmp.check_for_codominance((*phen).id());
    }

    std::swap(tmp, dst);
}

bool
penetrance_model::penetrance_informative_phenotype(uint id) const
{
    //lint --e{777} <- Comparing floats for equality a bad idea, but only
    //                 thing we got for now.

    if(genotype_informative_phenotype(id)) return true;

    // Otherwise, we must search for a pair that are different.

    // Unphased check

    unphased_penetrance_iterator unph    = unphased_penetrance_begin(id);
    unphased_penetrance_iterator unphend = unphased_penetrance_end  (id);

    if(unph != unphend)
    {
      double unph_pen = *unph;

      ++unph;

      // Search for non-equal element
      for( ; unph != unphend && unph_pen == *unph; ++unph)
        ;

      if(unph != unphend) return true;
    }


    // Phased check.

    phased_penetrance_iterator ph    = phased_penetrance_begin(id);
    phased_penetrance_iterator phend = phased_penetrance_end  (id);

    if(unph != unphend)
    {
      double ph_pen = *ph;

      ++ph;

      // Search for non-equal element
      for( ; ph != phend && ph_pen == *ph; ++ph)
        ;

      if(ph != phend) return true;
    }

    // No differences?  Ok, then we're not informative

    return false;
}

} // End namespace MLOCUS
} // End namespace SAGE

