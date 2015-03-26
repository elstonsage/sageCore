//============================================================================
//  File:       genotype.cpp
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
//============================================================================
//
#include "mlocus/genotype.h"
 
namespace SAGE   {
namespace MLOCUS {

const char  SEPARATORS[4] = {'/', '<', '>', '\0'};

#ifdef _DEBUG
    static  const char  REMAP[]   = "~";
    static  const char  MISSING[] = "*";
#else
    static  const char  REMAP[]   = "~remap";
    static  const char  MISSING[] = "*missing";
#endif

//============================================================================
//  IMPLEMENTATION: genotype model support classes and functions
//============================================================================
//
bool
valid_allele_name(const string& s)
{
    // FIXME: comparison may need to be made case-insensitive
    return s.size() > 0; // && s != REMAP; 
}

string
trim(const string& src)
{
    string  dst(src);
    uint    i, j;

    for (i = j = 0;  i < src.size();  ++i)
    {
        if ( !isspace((int)(uchar) src[i]) )
        {
            dst[j] = src[i];
            ++j;
        }
    }
    dst.resize(j);

    return dst;
}

namespace PRIVATE
{
  allele_info   invalid_allele_info            = allele_info();
  genotype_info invalid_phased_genotype_info   = genotype_info();
  genotype_info invalid_unphased_genotype_info = genotype_info();

//============================================================================
//  IMPLEMENTATION: genotype_model_info
//============================================================================
//

genotype_model_info::genotype_model_info(GenotypeModelType t)
  : name(),
    missing_allele_name(MISSING),
    alleles(),
    allele_names(),
    remap_buffer(),
    dynamic_alleles(false),
    my_phased_genotypes(),
    my_unphased_genotypes(),
    my_phased_genotypes_quick_lookup(),
    my_unphased_genotypes_quick_lookup()
{
    memcpy(separators, SEPARATORS, sizeof(SEPARATORS));
    
    configure_sex_type(t);
}

genotype_model_info::genotype_model_info(const genotype_model_info& gm)
  : name(gm.name),
    missing_allele_name(gm.missing_allele_name),
    alleles(gm.alleles),
    allele_names(gm.allele_names),
    remap_buffer(gm.remap_buffer),
    dynamic_alleles(gm.dynamic_alleles),
    my_type(gm.my_type),
    my_sex_allele_info(gm.my_sex_allele_info),
    my_phased_genotypes(),
    my_unphased_genotypes(),
    my_phased_genotypes_quick_lookup(),
    my_unphased_genotypes_quick_lookup()
{
    memcpy(separators, gm.separators, sizeof(SEPARATORS));
    
    rebuild_genotypes();
}

genotype_model_info::genotype_model_info(const string& n, GenotypeModelType t)
  : name(n),
    missing_allele_name(MISSING), 
    alleles(),
    allele_names(),
    remap_buffer(),
    dynamic_alleles(false),
    my_phased_genotypes(),
    my_unphased_genotypes(),
    my_phased_genotypes_quick_lookup(),
    my_unphased_genotypes_quick_lookup()
{
    memcpy(separators, SEPARATORS, sizeof(SEPARATORS));
    
    configure_sex_type(t);
}

void genotype_model_info::rebuild_genotypes()
{
  my_phased_genotypes.clear();
  my_unphased_genotypes.clear();
  
  switch(my_type)
  {
    case AUTOSOMAL : add_autosomal_genotypes(); return;
    case X_LINKED  : add_x_linked_genotypes(); return;
    case Y_LINKED  : add_y_linked_genotypes(); return;
  }
}

void genotype_model_info::add_autosomal_genotypes()
{
  my_phased_genotypes_quick_lookup.resize(alleles.size() * alleles.size());
  my_unphased_genotypes_quick_lookup.resize(alleles.size() * (alleles.size() + 1) / 2);

  for(allele_iterator a1 = allele_iterator(alleles.begin());
                      a1 != allele_iterator(alleles.end()); ++a1)
  {
    add_autosomal_genotypes_for_allele(a1);
  }
  
  initialize_sex_iterators();
}

void genotype_model_info::add_x_linked_genotypes()
{
  my_phased_genotypes_quick_lookup.resize(alleles.size() * alleles.size() + alleles.size());
  my_unphased_genotypes_quick_lookup.resize(alleles.size() * (alleles.size() + 1) / 2 + alleles.size());
  
  initialize_sex_iterators();

  for(allele_iterator a1 = allele_iterator(alleles.begin());
                      a1 != allele_iterator(alleles.end()); ++a1)
  {
    add_autosomal_genotypes_for_allele(a1);
    add_x_linked_male_genotype_for_allele(a1);
  }
}

void genotype_model_info::add_y_linked_genotypes()
{
  my_phased_genotypes_quick_lookup.resize(alleles.size() + 1);
  my_unphased_genotypes_quick_lookup.resize(alleles.size() + 1);

  initialize_sex_iterators();

  add_y_linked_female_genotype();
  add_y_linked_male_genotypes();
}

void genotype_model_info::add_autosomal_genotypes_for_allele(allele_iterator a1)
{
    size_t phased_id = my_phased_genotypes.size();
    size_t unphased_id = my_unphased_genotypes.size();
    
    for(allele_iterator a2 = allele_iterator(alleles.begin());
                        a2 <= a1; ++a2)
    {
      int flipped_phased_id = (a1 == a2) ? phased_id : phased_id+1;
      
      add_unphased_genotype(PRIVATE::genotype_info(a2->my_info, a1->my_info,
                            phased_id, flipped_phased_id, unphased_id, this));
          
      add_phased_genotype(PRIVATE::genotype_info(a2->my_info, a1->my_info,
                          phased_id, flipped_phased_id, unphased_id, this));
          
      if(a1 != a2)
      {
        add_phased_genotype(PRIVATE::genotype_info(a1->my_info, a2->my_info,
                            flipped_phased_id, phased_id, unphased_id, this));
          
        ++phased_id;
      }
      ++phased_id;
      ++unphased_id;
    }
}

void genotype_model_info::add_x_linked_male_genotype_for_allele(allele_iterator a1)
{
    size_t phased_id = my_phased_genotypes.size();
    size_t unphased_id = my_unphased_genotypes.size();
    
    add_unphased_genotype(PRIVATE::genotype_info(a1->my_info, &my_sex_allele_info,
                          phased_id, phased_id, unphased_id, this));
        
    add_phased_genotype(PRIVATE::genotype_info(a1->my_info, &my_sex_allele_info,
                        phased_id, phased_id, unphased_id, this));
}

void genotype_model_info::add_y_linked_female_genotype()
{
    size_t phased_id   = my_phased_genotypes.size();
    size_t unphased_id = my_unphased_genotypes.size();
  
    add_unphased_genotype(PRIVATE::genotype_info(&my_sex_allele_info,  &my_sex_allele_info,
                          phased_id, phased_id, unphased_id, this));
        
    add_phased_genotype(PRIVATE::genotype_info(&my_sex_allele_info,  &my_sex_allele_info,
                        phased_id, phased_id, unphased_id, this));
}

void genotype_model_info::add_y_linked_male_genotypes()
{
    size_t phased_id   = my_phased_genotypes.size();
    size_t unphased_id = my_unphased_genotypes.size();
  
    for(allele_iterator a1 = allele_iterator(alleles.begin());
                        a1 != allele_iterator(alleles.end()); ++a1, ++phased_id, ++unphased_id)
    {
      add_unphased_genotype(PRIVATE::genotype_info(&my_sex_allele_info, a1->my_info,
                            phased_id, phased_id, unphased_id, this));
          
      add_phased_genotype(PRIVATE::genotype_info(&my_sex_allele_info, a1->my_info,
                          phased_id, phased_id, unphased_id, this));
    }
}

void genotype_model_info::add_phased_genotype(const genotype_info& info)
{
  if(my_type == AUTOSOMAL)
  {
    add_unsexed_phased_genotype(info);
  }
  else
  {
    if(phased_genotype(&info).is_male_compatible())
      add_male_phased_genotype(info);
    else
      add_female_phased_genotype(info);
  }
}
void genotype_model_info::add_unsexed_phased_genotype(const genotype_info& info)
{
  my_phased_genotypes.insert(my_phased_genotypes.end(), info);

  my_phased_genotypes_quick_lookup[info.my_phased_id] =
        &*my_phased_genotypes.get<PRIVATE::id>().find(info.my_phased_id);
}
void genotype_model_info::add_male_phased_genotype(const genotype_info& info)
{
  PhasedGenotypeBySequence::iterator i = 
      my_phased_genotypes.insert(my_male_phased_genotype_end, info).first;

  my_phased_genotypes_quick_lookup[info.my_phased_id] =
        &*my_phased_genotypes.get<PRIVATE::id>().find(info.my_phased_id);

  if(my_male_phased_genotype_begin == my_male_phased_genotype_end)
    my_male_phased_genotype_begin = my_female_phased_genotype_end = i;
}
void genotype_model_info::add_female_phased_genotype(const genotype_info& info)
{
  my_phased_genotypes.insert(my_female_phased_genotype_end, info);

  my_phased_genotypes_quick_lookup[info.my_phased_id] =
        &*my_phased_genotypes.get<PRIVATE::id>().find(info.my_phased_id);
        
  my_female_phased_genotype_begin = my_phased_genotypes.begin();
}

void genotype_model_info::add_unphased_genotype(const genotype_info& info)
{
  if(my_type == AUTOSOMAL)
  {
    add_unsexed_unphased_genotype(info);
  }
  else
  {
    if(unphased_genotype(&info).is_male_compatible())
      add_male_unphased_genotype(info);
    else
      add_female_unphased_genotype(info);
  }
}
void genotype_model_info::add_unsexed_unphased_genotype(const genotype_info& info)
{
  my_unphased_genotypes.insert(my_unphased_genotypes.end(), info);

  my_unphased_genotypes_quick_lookup[info.my_unphased_id] =
        &*my_unphased_genotypes.get<PRIVATE::id>().find(info.my_unphased_id);
}

void genotype_model_info::add_male_unphased_genotype(const genotype_info& info)
{
  UnphasedGenotypeBySequence::iterator i = 
      my_unphased_genotypes.insert(my_male_unphased_genotype_end, info).first;

  my_unphased_genotypes_quick_lookup[info.my_unphased_id] =
        &*my_unphased_genotypes.get<PRIVATE::id>().find(info.my_unphased_id);

  if(my_male_unphased_genotype_begin == my_male_unphased_genotype_end)
    my_male_unphased_genotype_begin = my_female_unphased_genotype_end = i;
}

void genotype_model_info::add_female_unphased_genotype(const genotype_info& info)
{
  my_unphased_genotypes.insert(my_female_unphased_genotype_end, info);

  my_unphased_genotypes_quick_lookup[info.my_unphased_id] =
        &*my_unphased_genotypes.get<PRIVATE::id>().find(info.my_unphased_id);
        
  my_female_unphased_genotype_begin = my_unphased_genotypes.begin();
}

void genotype_model_info::configure_sex_type(GenotypeModelType gt)
{
  my_type = gt;
  
  switch(gt)
  {
    case X_LINKED :
      my_sex_allele_info = allele_info(allele_info::st_NULL_Y);
      break;
    case Y_LINKED :
      my_sex_allele_info = allele_info(allele_info::st_NULL_X);
      break;
    case AUTOSOMAL :
      my_sex_allele_info = allele_info();
  }

  rebuild_genotypes();
}

void genotype_model_info::initialize_sex_iterators()
{
  my_male_phased_genotype_begin = my_phased_genotypes.begin();
  my_male_phased_genotype_end   = my_phased_genotypes.end();  

  my_female_phased_genotype_begin = my_phased_genotypes.begin();
  my_female_phased_genotype_end   = my_phased_genotypes.end();  

  my_male_unphased_genotype_begin = my_unphased_genotypes.begin();
  my_male_unphased_genotype_end   = my_unphased_genotypes.end();  

  my_female_unphased_genotype_begin = my_unphased_genotypes.begin();
  my_female_unphased_genotype_end   = my_unphased_genotypes.end();  
}

} // End PRIVATE namespace

//============================================================================
//  IMPLEMENTATION: genotype_model
//============================================================================
//
allele
genotype_model::get_allele(const string& n) const
{
    if(!is_autosomal() && n == my_info->my_sex_allele_info.name)
      return allele(&my_info->my_sex_allele_info);
  
    allele_map::const_iterator i = my_info->allele_names.find(n);

    //lint -e{58}
    bool test = (i != my_info->allele_names.end());

    return test ? allele(&my_info->alleles[i->second]) : allele();
}

phased_genotype
genotype_model::get_phased_genotype(const string& n) const
{
    string an1, an2;

    Ordering order;

    //lint -e{534}
    parse_genotype_name(n, an1, an2, order);

    if(!an1.size() || !an2.size() || order == Unphased)
      return phased_genotype();
    
    allele a1 = get_allele(an1);
    allele a2 = get_allele(an2);

    if(order == PhasedBackward) std::swap(a1, a2);

    return get_phased_genotype(a1, a2);
}

unphased_genotype
genotype_model::get_unphased_genotype(const string& n) const
{
    string an1, an2;

    Ordering order;

    //lint -e{534}
    parse_genotype_name(n, an1, an2, order);

    if(!an1.size() || !an2.size() || order != Unphased) return unphased_genotype();

    allele a1 = get_allele(an1);
    allele a2 = get_allele(an2);

    return get_unphased_genotype(a1, a2);
}

//----------
//
genotype_model
genotype_model::clone() const
{
    genotype_model  gm(*this);

    gm.uniquify();

    return gm;
}

//----------
//
void
genotype_model::add_allele(const string& n, double freq)
{
    string  nm = parse_allele_name(n);

    if (nm.size() == 0) return;

    allele_map::iterator i = my_info->allele_names.find(nm);

    //lint -e{81}
    if(i != my_info->allele_names.end()) return;

    uniquify();

    uint aid = my_info->alleles.size();

    my_info->allele_names[nm] = aid;

    my_info->alleles.push_back(allele_info(nm, freq, aid));

    my_info->rebuild_genotypes();
}

void
genotype_model::normalize()
{
  uniquify();

  double total = 0.0;

  for(allele_iterator i = allele_begin(); i != allele_end(); ++i)
    total += (*i).frequency();

  if(total != 0.0)
  {
    for(uint i = 0; i != allele_count(); ++i)
    {
      allele_info& a = my_info->alleles[i];
      a = allele_info(a.name, a.frequency/total, a.id);
    }
  }
  else
  {
    double d = 1.0 / allele_count();

    for(uint i = 0; i != allele_count(); ++i)
    {
      allele_info& a = my_info->alleles[i];
      a = allele_info(a.name, d, a.id);
    }
  }
}

void
genotype_model::clear()
{
    uniquify();

    my_info->alleles.clear();
    my_info->allele_names.clear();
    my_info->remap_buffer.clear();
    my_info->rebuild_genotypes();
}


//----------
//
void
genotype_model::mark_for_remap(const string& n)
{
    string  nm = parse_allele_name(n);
    
    allele a = get_allele(nm);

    if (nm.size() > 0  &&  a.is_valid() && !a.is_sex_allele())
    {
        uniquify();
        my_info->remap_buffer.push_back(nm);
    }
}


void
genotype_model::remap()
{
    uniquify();
    remap(*this, *this);
}


void
genotype_model::remap(genotype_model& dst)
{
    uniquify();
    remap(*this, dst);
}


//----------
//
string
genotype_model::parse_allele_name(const string& n) const
{
    string  nm(n);

    nm = trim(nm);

    if (valid_allele_name(nm))
    {
        return nm;
    }
    else
    {
        return string();
    }
}


string
genotype_model::parse_genotype_name(const string& gstr, char sepc) const
{
  string aname1 = "";
  string aname2 = "";
  Ordering order;

  return parse_genotype_name(gstr, aname1, aname2, order, sepc);
}

string
genotype_model::parse_genotype_name
(const string& gstr, string& aname1, string& aname2, Ordering& order, char sepc) const
{
    string  gname;
    string  seps;
    string::size_type    sloc;

    //- Before starting the parse, all names are set to the empty string and
    //  ordering set to the default of unordered.
    //

    aname1 = aname2 = gname;

    order = Unphased;

    // If gstr is empty, there's really not much to be done.

    if(!gstr.size())
      return gname;

    //- The separator character argument is used as both a separator and
    //  a flag.  The null character indicates that we should parse based
    //  on the default separators stored in this marker instance.
    //
    if (sepc == '\0')
    {
        //lint -e{534}
        seps.assign(my_info->separators, 3);
    }
    else
    {
        seps = sepc;
    }
    
    //- Look for the first occurrence of the separator charactor in the 
    //  genotype name.  Make sure that it is a valid location.
    //
    sloc = gstr.find_first_of(seps);

    if (sloc == string::npos  ||  sloc == 0  ||  sloc == (gstr.size()-1))
    {
        return gname;
    }

    //- Make sure that there is not a second occurence of any separator 
    //  character in the genotype name.
    //
    if (gstr.find_first_of(seps, sloc+1) != string::npos)
    {
        return gname;
    }

    //- Determine the ordering in the case that it is unspecified.
    //
    if (sepc == '\0')
    {
        sepc = gstr[sloc];

        if (sepc == my_info->separators[0])
        {
            order = Unphased;
        }
        else if (sepc == my_info->separators[1])
        {
            order = PhasedForward;
        }
        else if (sepc == my_info->separators[2])
        {
            order = PhasedBackward;
        }
        else
        {
            return gname;
        }
    }

    //- Determine the allele names.
    //
    if (order == PhasedBackward)
    {
        aname2 = gstr.substr(0, sloc);
        aname1 = gstr.substr(sloc + 1, gstr.size() - (sloc + 1));
    }
    else
    {
        aname1 = gstr.substr(0, sloc);
        aname2 = gstr.substr(sloc + 1, gstr.size() - (sloc + 1));
    }

    aname1 = parse_allele_name(aname1);
    aname2 = parse_allele_name(aname2);

    //- Build a genotype name.
    //
    if (aname1.size() > 0  &&  aname2.size() > 0)
    {
        sepc  = (order == Unphased) ? my_info->separators[0] : my_info->separators[1];
        gname = aname1 + sepc + aname2;
    }

    return gname;
}


void
genotype_model::uniquify()
{
    if (!my_info.unique())
    {
        my_info = boost::shared_ptr<PRIVATE::genotype_model_info>
                      (new PRIVATE::genotype_model_info(*my_info));

        my_info->rebuild_genotypes();
    }
}


//lint -e{1764}
void
genotype_model::remap(const genotype_model& src, genotype_model& dst)
{
    if (src.my_info->remap_buffer.size() > 0)
    {
        genotype_model  tmp;

        double remap_freq = 0.0;

        allele_vector::iterator     ax = src.my_info->alleles.begin();
        allele_vector::iterator     al = src.my_info->alleles.end();

        //- Scan through the allele vector, looking at the alleles' names.
        //
        for (;  ax != al;  ++ax)
        {
            allele_name_buffer::iterator    nx = src.my_info->remap_buffer.begin();
            allele_name_buffer::iterator    nl = src.my_info->remap_buffer.end();

            //- If the current allele's name is NOT in the remap buffer,
            //  then it can be added to the new allele vector directly.
            //
            if (std::find(nx, nl, ax->name) == nl)
            {
                tmp.add_allele(ax->name, ax->frequency);
            }
            else
            {
                remap_freq += ax->frequency;
            }
        }

        //- Make sure that an allele representing the REMAPPed alleles is 
        //  also added to the new allele vector.
        //
        tmp.add_allele(REMAP, remap_freq);

        //- Finally, copy the new data into the destination.
        //      
        if(&src != &dst)
        {
            dst.my_info->name                = src.my_info->name;
            dst.my_info->separators[0]       = src.my_info->separators[0];
            dst.my_info->separators[1]       = src.my_info->separators[1];
            dst.my_info->separators[2]       = src.my_info->separators[2];
            dst.my_info->separators[3]       = src.my_info->separators[3];
            dst.my_info->missing_allele_name = src.my_info->missing_allele_name;
            dst.my_info->dynamic_alleles     = src.my_info->dynamic_alleles;
            dst.my_info->my_type             = src.my_info->my_type;
            dst.my_info->my_sex_allele_info  = src.my_info->my_sex_allele_info;
            // Don't have to copy genotypes, they're being rebuilt anyway.
        }

        std::swap(dst.my_info->alleles, tmp.my_info->alleles);
        std::swap(dst.my_info->allele_names, tmp.my_info->allele_names);
        
        dst.my_info->rebuild_genotypes();
    }
}

//============================================================================
//  IMPLEMENTATION: child_genotype_set
//============================================================================
//
child_genotype_set::child_genotype_set()
  : my_size(0)
{}


child_genotype_set::child_genotype_set
    (const phased_genotype& a,
     const phased_genotype& b,
     bool                   r)
  : my_size(0)
{
    // Can do nothing if the genotypes come from different models.
    if(a.my_info->my_ginfo != b.my_info->my_ginfo) return;
    
    const PRIVATE::genotype_model_info& ginf = *a.my_info->my_ginfo;

    my_genotypes[0] = ginf.get_phased_genotype(a.allele1(), b.allele1());
    my_genotypes[1] = ginf.get_phased_genotype(a.allele1(), b.allele2());
    my_genotypes[2] = ginf.get_phased_genotype(a.allele2(), b.allele1());
    my_genotypes[3] = ginf.get_phased_genotype(a.allele2(), b.allele2());

    if(!r) my_size = 4;
    else   reduce();
}

child_genotype_set::child_genotype_set
    (const unphased_genotype& a,
     const unphased_genotype& b,
     bool                     r)
  : my_size(0)
{
    // Can do nothing if the genotypes come from different models.
    if(a.my_info->my_ginfo != b.my_info->my_ginfo) return;

    const PRIVATE::genotype_model_info& ginf = *a.my_info->my_ginfo;

    my_genotypes[0] = ginf.get_phased_genotype(a.allele1(), b.allele1());
    my_genotypes[1] = ginf.get_phased_genotype(a.allele1(), b.allele2());
    my_genotypes[2] = ginf.get_phased_genotype(a.allele2(), b.allele1());
    my_genotypes[3] = ginf.get_phased_genotype(a.allele2(), b.allele2());

    if(!r) my_size = 4;
    else   reduce();
}

void child_genotype_set::reduce()
{
  uint        i, j;
  bool        skip;

  my_size = 1;

  //- Iterate over the elements and remove duplicates.

  for (i = 1;  i < 4;  ++i)
  {
      skip = false;

      //- Check for equivalence with each of the previous elements of
      //  'my_genotypes'.

      for (j = 0;  j < my_size;  ++j)
      {
          if ( my_genotypes[i].equivalent(my_genotypes[j]) )
          {
              skip = true;
              break;
          }
      }
      
      //- If no skip is indicated, we make sure to keep this element

      if (skip == false)
      {
          my_genotypes[my_size++] = my_genotypes[i];
      }
  }
}


void do_allele_tests(ostream& o, const genotype_model& m)
{
    // Allele count tests
    size_t diff  = m.allele_end() - m.allele_begin();
    size_t cnt   = m.allele_count();

    if (cnt == diff)
        o << ".............allele count OK" << endl;
    else
        o << ".............allele count FAIL" << endl;

    // Allele Iterator test
    
    cnt = 0;
    
    allele_iterator     af = m.allele_begin();      
    allele_iterator     al = m.allele_end();        

    bool good = true;

    for ( ;  af != al;  ++af, ++cnt)
    {
        if(!af->is_valid())
        {
            o << ".............allele invalid" << endl;
            good = false;
        }

        if (af->name() != m.get_allele(af->name()).name())
        {
            o << "............allele lookup FAIL at: " << af->name() << endl;
            good = false;
        }

        if (*af != m.get_allele(af->name()))
        {
            o << "............allele lookup FAIL at: " << af->name() << endl;
            good = false;
        }
    }
    if (good)
    {
        o << "............allele lookup OK" << endl;
    }
    assert(cnt == m.allele_count());

    // Missing Allele Lookup tests
    assert(!m.get_allele(m.missing_allele_name()).is_valid());
    assert(!m.get_allele(m.allele_count()).is_valid());
    if(m.is_autosomal())
    {
      assert(!m.get_allele((uint)-1).is_valid());
      
      assert(m.get_allele((uint)-1) == m.get_allele(m.allele_count()));
      assert(m.get_allele((uint)-1) == m.get_allele(m.missing_allele_name()));
    }
    else
    {
      assert(m.get_allele((uint)-1).is_valid());
      
      assert(m.get_allele((uint)-1) == m.get_sex_specific_allele());
    }
}

void test_unphased_genotype_autosomal(unphased_genotype g)
{
  if(!g.homozygous())
      assert(g.get_equivalent_phased_genotype1() != g.get_equivalent_phased_genotype2());
    
  assert(g.is_male_compatible() && g.is_female_compatible());
}

void test_unphased_genotype_x_linked(unphased_genotype g)
{
  assert(g.is_sex_specific());
  
  assert(!g.allele1().is_sex_allele());
    
  if(!g.allele2().is_sex_allele())  assert( g.is_female_compatible() && !g.is_male_compatible());
  else                              assert(!g.is_female_compatible() &&  g.is_male_compatible());
}

void test_unphased_genotype_y_linked(unphased_genotype g)
{
  assert(g.is_sex_specific());

  assert(g.allele1().is_sex_allele());
    
  if(!g.allele2().is_sex_allele())   assert( g.is_male_compatible() && !g.is_female_compatible());
  else                               assert(!g.is_male_compatible() &&  g.is_female_compatible());
}

bool test_unphased_genotype(ostream& o, const genotype_model& m,
                            unphased_genotype g)
{
  bool good = true;

  if (g != m.get_unphased_genotype(g.get_id()))
  {
      o << "....phased genotype index FAIL at: " << g.get_id() << endl;
      good = false;
  }

  if (g.get_id() != m.get_unphased_genotype(g.name()).get_id())
  {
      o << ".unphased genotype lookup FAIL at: " << g.get_id() << endl;
      good = false;
  }
  assert(g == g.get_equivalent_phased_genotype1().get_equivalent_unphased_genotype());
  assert(g == g.get_equivalent_phased_genotype2().get_equivalent_unphased_genotype());
  
  if(m.is_autosomal())
  {
    test_unphased_genotype_autosomal(g);
  }
  else if(m.is_x_linked())
  {
    test_unphased_genotype_x_linked(g);
  }
  else
  {
    test_unphased_genotype_y_linked(g);
  }
  
  return good;
}

void do_unphased_genotype_tests(ostream& o, const genotype_model& m)
{
    uint diff  = m.unphased_genotype_count();
    uint  cnt;
    
         if(m.is_autosomal()) cnt = ((m.allele_count() + 1)*m.allele_count())/2;
    else if(m.is_x_linked())  cnt = ((m.allele_count() + 3)*m.allele_count())/2;
    else                      cnt = m.allele_count()+1;

    if (cnt == diff)
        o << "..unphased genotype count OK" << endl;
    else
        o << "..unphased genotype count FAIL" << endl;

    unphased_genotype_iterator  gf = m.unphased_genotype_begin();        
    unphased_genotype_iterator  gl = m.unphased_genotype_end();      

    bool good = true;

    for ( ;  gf != gl;  ++gf)
    {
      if(!test_unphased_genotype(o, m, *gf))
        good = false;
    }
    if (good)
    {
        o << ".unphased genotype lookup OK" << endl;
    }
}

bool test_phased_genotype(ostream& o, const genotype_model& m,
                          phased_genotype g)
{
  bool good = true;
  
  if (g != m.get_phased_genotype(g.get_id()))
  {
      o << "....phased genotype index FAIL at: " << g.get_id() << endl;
      good = false;
  }

  if (g.get_id() != m.get_phased_genotype(g.name()).get_id())
  {
      o << "...phased genotype lookup FAIL at: " << g.get_id() << endl;
      good = false;
  }
  
  if(m.is_autosomal())
  {
    assert(!g.is_sex_specific());
    
    if(g.allele1() < g.allele2())
      assert(g == g.get_equivalent_unphased_genotype().get_equivalent_phased_genotype1());
    else
      assert(g == g.get_equivalent_unphased_genotype().get_equivalent_phased_genotype2());
  }
  else if(m.is_x_linked())
  {
    assert(g.is_sex_specific());
    assert(!g.allele1().is_sex_allele());
    
    if(g.allele2().is_sex_allele())
    {
      assert(g.is_male_compatible() && !g.is_female_compatible());
      assert(g == g.get_equivalent_unphased_genotype().get_equivalent_phased_genotype1());
      assert(g == g.get_equivalent_unphased_genotype().get_equivalent_phased_genotype2());
    }
    else
    {
      assert(!g.is_male_compatible() && g.is_female_compatible());
      if(g.allele1() < g.allele2())
        assert(g == g.get_equivalent_unphased_genotype().get_equivalent_phased_genotype1());
      else
        assert(g == g.get_equivalent_unphased_genotype().get_equivalent_phased_genotype2());
    }
  }

  return good;
}

void do_phased_genotype_tests(ostream& o, const genotype_model& m)
{
    uint diff = m.phased_genotype_count();
    uint cnt  = m.allele_count()*m.allele_count();

         if(m.is_autosomal()) cnt = m.allele_count()*m.allele_count();
    else if(m.is_x_linked())  cnt = m.allele_count()*m.allele_count() + m.allele_count();
    else                      cnt = m.allele_count()+1;

    if (cnt == diff)
        o << "....phased genotype count OK" << endl;
    else
        o << "....phased genotype count FAIL" << endl;

    phased_genotype_iterator    ogf = m.phased_genotype_begin();       
    phased_genotype_iterator    ogl = m.phased_genotype_end();     

    bool good = true;

    for ( ;  ogf!=ogl;  ++ogf)
    {
      if(!test_phased_genotype(o,m,*ogf))
        good = false;
    }

    if (good)
    {
        o << "...phased genotype lookup OK" << endl;
    }
}

void test_sex_based_iterators(const genotype_model& m)
{
  if(m.is_autosomal())
  {
    assert(m.unphased_genotype_begin() == m.unphased_genotype_begin(MPED::SEX_MISSING));
    assert(m.unphased_genotype_begin() == m.unphased_genotype_begin(MPED::SEX_MALE));
    assert(m.unphased_genotype_begin() == m.unphased_genotype_begin(MPED::SEX_FEMALE));

    assert(m.unphased_genotype_end()   == m.unphased_genotype_end(MPED::SEX_MISSING));
    assert(m.unphased_genotype_end()   == m.unphased_genotype_end(MPED::SEX_MALE));
    assert(m.unphased_genotype_end()   == m.unphased_genotype_end(MPED::SEX_FEMALE));

    assert(m.phased_genotype_begin() == m.phased_genotype_begin(MPED::SEX_MISSING));
    assert(m.phased_genotype_begin() == m.phased_genotype_begin(MPED::SEX_MALE));
    assert(m.phased_genotype_begin() == m.phased_genotype_begin(MPED::SEX_FEMALE));

    assert(m.phased_genotype_end()   == m.phased_genotype_end(MPED::SEX_MISSING));
    assert(m.phased_genotype_end()   == m.phased_genotype_end(MPED::SEX_MALE));
    assert(m.phased_genotype_end()   == m.phased_genotype_end(MPED::SEX_FEMALE));
  }
  else
  {
    assert(m.unphased_genotype_begin() == m.unphased_genotype_begin(MPED::SEX_MISSING));
    assert(m.unphased_genotype_begin() != m.unphased_genotype_begin(MPED::SEX_MALE));
    assert(m.unphased_genotype_begin() == m.unphased_genotype_begin(MPED::SEX_FEMALE));

    assert(m.unphased_genotype_end()   == m.unphased_genotype_end(MPED::SEX_MISSING));
    assert(m.unphased_genotype_end()   == m.unphased_genotype_end(MPED::SEX_MALE));
    assert(m.unphased_genotype_end()   != m.unphased_genotype_end(MPED::SEX_FEMALE));

    assert(m.unphased_genotype_begin(MPED::SEX_MALE) == m.unphased_genotype_end(MPED::SEX_FEMALE));
    
    for(unphased_genotype_iterator i = m.unphased_genotype_begin(MPED::SEX_FEMALE);
        i != m.unphased_genotype_begin(MPED::SEX_FEMALE); ++i)
    {
      assert(i->is_sex_specific());
      assert(i->is_female_compatible() && !i->is_male_compatible());
    }
    for(unphased_genotype_iterator i = m.unphased_genotype_begin(MPED::SEX_MALE);
        i != m.unphased_genotype_begin(MPED::SEX_MALE); ++i)
    {
      assert(i->is_sex_specific());
      assert(!i->is_female_compatible() && i->is_male_compatible());
    }

    assert(m.phased_genotype_begin() == m.phased_genotype_begin(MPED::SEX_MISSING));
    assert(m.phased_genotype_begin() != m.phased_genotype_begin(MPED::SEX_MALE));
    assert(m.phased_genotype_begin() == m.phased_genotype_begin(MPED::SEX_FEMALE));

    assert(m.phased_genotype_end()   == m.phased_genotype_end(MPED::SEX_MISSING));
    assert(m.phased_genotype_end()   == m.phased_genotype_end(MPED::SEX_MALE));
    assert(m.phased_genotype_end()   != m.phased_genotype_end(MPED::SEX_FEMALE));

    assert(m.phased_genotype_begin(MPED::SEX_MALE) == m.phased_genotype_end(MPED::SEX_FEMALE));
    
    for(phased_genotype_iterator i = m.phased_genotype_begin(MPED::SEX_FEMALE);
        i != m.phased_genotype_begin(MPED::SEX_FEMALE); ++i)
    {
      assert(i->is_sex_specific());
      assert(i->is_female_compatible() && !i->is_male_compatible());
    }
    for(phased_genotype_iterator i = m.phased_genotype_begin(MPED::SEX_MALE);
        i != m.phased_genotype_begin(MPED::SEX_MALE); ++i)
    {
      assert(i->is_sex_specific());
      assert(!i->is_female_compatible() && i->is_male_compatible());
    }
  }
}

void genotype_model_test(ostream& o, const genotype_model& m)
{
    do_allele_tests(o, m);
    do_unphased_genotype_tests(o, m);
    do_phased_genotype_tests(o,m);
    
    test_sex_based_iterators(m);

    o << endl;
}

void genotype_model_print(ostream& o, const genotype_model& m)
{
    //lint --e{522}

    o << "    alleles:";
    allele_iterator     af = m.allele_begin();
    allele_iterator     al = m.allele_end();
    for (;  af != al;  ++af)
    {
        o << " " << af->name();
    }
    o << endl;

    o << "unphased genotypes:";
    unphased_genotype_iterator  ugf = m.unphased_genotype_begin();
    unphased_genotype_iterator  ugl = m.unphased_genotype_end();
    for (;  ugf != ugl;  ++ugf)
    {
        o << " " << (*ugf).name();
    }
    o << endl;

    o << "phased genotypes:";
    phased_genotype_iterator    pgf = m.phased_genotype_begin();
    phased_genotype_iterator    pgl = m.phased_genotype_end();
    for (;  pgf != pgl;  ++pgf)
    {
        o << " " << (*pgf).name();
    }
    o << endl;
}

} // End namespace MLOCUS
} // End namespace SAGE

