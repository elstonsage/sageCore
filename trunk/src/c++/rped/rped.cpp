// rped.cpp -- Implements the PedigreeSort function
#include <cmath>
#include "mped/sp.h"
#include "rped/rped.h"

namespace SAGE {
namespace RPED {

#define DEBUG_PEDSORT(x)

RefMPedInfo::RefMPedInfo()
{
  set_sex_code_male           ("M");
  set_sex_code_female         ("F");
  set_sex_code_unknown        ("");
  set_individual_missing_code ("");

  my_trait_info.clear();
}

size_t 
RefMPedInfo::trait_find(const std::string &name) const
{
  std::string tn = toUpper(name);

  for( size_t t = 0; t < trait_count(); ++t )
    if( toUpper(trait_info(t).name()) == tn || toUpper(trait_info(t).alias_name()) == tn )
      return t;

  return (size_t)-1;
}

size_t 
RefMPedInfo::add_trait(const string &trait_name, RefTraitInfo::trait_t trait_type, RefTraitInfo::trait_use use)
{
  if(trait_exists(trait_name))
  {
    return trait_find(trait_name);
  }
  else
  {
    my_trait_info.push_back(RefTraitInfo(trait_name, trait_type, use));

    return trait_count() - 1;
  }
}


void
RefMPedInfo::remove_trait_info(size_t t_id)
{
  assert(t_id < my_trait_info.size());
  my_trait_info.erase(my_trait_info.begin() + t_id);
}

void
RefMPedInfo::remove_marker_info(size_t m_id)
{
  if( my_marker_info.name_find(my_marker_info.name(m_id)) != my_marker_info.name_end() )
    my_marker_info.erase(my_marker_info.name(m_id));
}

void
RefMPedInfo::remove_marker_info(string m_name)
{
  if( my_marker_info.name_find(m_name) != my_marker_info.name_end() )
    my_marker_info.erase(m_name);
}

void
RefMPedInfo::read_pheno_reader_info(const LSFBase* params, cerrorstream &errors)
{
  AttrList::const_iterator a, a1;
  AttrVal v;

  LSFList::const_iterator i;
  for( i = params->List()->begin(); i != params->List()->end(); ++i )
  {
    if( !*i || !(*i)->name().size() )
      continue;

    const LSFBase* param = *i;

    string name = toUpper( param->name() );

    if( name == "ALLELE_DELIMITER" )
    {
      v = attr_value(param, "ALLELE_DELIMITER",0);
      if( v.has_value() )
      {
        if( v.String().size() == 1 )
          my_pheno_reader_info.set_allele_delimiter(v.String()[0]);
        else
          my_pheno_reader_info.set_allele_delimiter(strip_ws(v.String())[0]);

        //cout << "ALLELE_DELIMITER = '" << v.String()[0] << "'" << endl;
      }
    }

    if( name == "ALLELE_MISSING" )
    {
      v = attr_value(param, "ALLELE_MISSING",0);
      if( v.has_value() )
      {
        my_pheno_reader_info.set_allele_missing(v.String());
        //cout << "ALLELE_MISSING = '" << v.String() << "'" << endl;
      }
    }

    if( name == "ALLELE_FREQUENCY" && param->attrs() )
    {
      if( (a=param->attrs()->find("equal")) != param->attrs()->end() )
        my_pheno_reader_info.set_allele_adjustment(PhenotypeReaderInfo::equal);
      else if( (a=param->attrs()->find("complement")) != param->attrs()->end() )
        my_pheno_reader_info.set_allele_adjustment(PhenotypeReaderInfo::comp);
      else if(    (a =param->attrs()->find("minimum")) != param->attrs()->end()
               && (a1=param->attrs()->find("maximum")) != param->attrs()->end() )
      {
        if(    a ->second.Real() < 0. || a ->second.Real() > 1.
            && a1->second.Real() < 0. || a1->second.Real() > 1. )
          errors << priority(error) << "Invalid 'minimum or maximum' value for parameter 'allele_frequency'." << endl;
        else
        {
          my_pheno_reader_info.set_allele_adjustment(PhenotypeReaderInfo::min_max);
          my_pheno_reader_info.set_min_allele_freq(a ->second.Real());
          my_pheno_reader_info.set_max_allele_freq(a1->second.Real());
        }
      }
      else if( (a =param->attrs()->find("minimum")) != param->attrs()->end() )
      {
        if( a ->second.Real() < 0. || a ->second.Real() > 1. )
          errors << priority(error) << "Invalid 'minimum' value for parameter 'allele_frequency'." << endl;
        else
        {
          my_pheno_reader_info.set_allele_adjustment(PhenotypeReaderInfo::min);
          my_pheno_reader_info.set_min_allele_freq(a->second.Real());
        }
      }
      else if( (a1 =param->attrs()->find("maximum")) != param->attrs()->end() )
      {
        if( a1->second.Real() < 0. || a1->second.Real() > 1. )
          errors << priority(error) << "Invalid 'maximum' value for parameter 'allele_frequency'." << endl;
        else
        {
          my_pheno_reader_info.set_allele_adjustment(PhenotypeReaderInfo::max);
          my_pheno_reader_info.set_max_allele_freq(a1->second.Real());
        }
      }
    }

    if( name == "COVARIATE_MOI"      || name == "COV_MOI" ||
        name == "COVARIATE_FUNCTION" || name == "COV_FUNC" )
    {
      v = attr_value(param, 0);
      if( v.has_value() && v.String().size() )
      {
        string cov_moi = toUpper(v.String());
        if( cov_moi == "ADD" || cov_moi == "ADDITIVE" || cov_moi == "A" )
          my_pheno_reader_info.set_covariate_moi("_ADD");
        else if( cov_moi == "DOM" || cov_moi == "DOMINANT" || cov_moi == "D" )
          my_pheno_reader_info.set_covariate_moi("_DOM");
        else if( cov_moi == "REC" || cov_moi == "RECESSIVE" || cov_moi == "R" )
          my_pheno_reader_info.set_covariate_moi("_REC");
        else
          errors << priority(error)
                 << "Invalid value for parameter 'covariate_function' specified.  Skipping..." << endl;
      }
      else
        errors << priority(error)
               << "No value for parameter 'covariate_function' specified.  Skipping..." << endl;

      if( has_attr(param, "BASE_ALLELE") )
      {
        v = attr_value(param, "BASE_ALLELE");
        if( v.has_value() && v.String().size() )
          my_pheno_reader_info.set_covariate_allele(v.String());
        else
          errors << priority(error)
                 << "No value for 'base_allele' specified.  Skipping..." << endl;
      }

      a = param->attrs()->find("ALLOW_HEMIZYGOTE");
      if( a == param->attrs()->end() )
        a = param->attrs()->find("ALLOW_HEMI");

      if( a != param->attrs()->end() && a->second.has_value() && a->second.String().size() )
      {
        string hemi = toUpper(a->second.String());
        if( hemi == "YES" || hemi == "TRUE" )
          my_pheno_reader_info.set_allow_hemizygote(true);
        else if( hemi == "NO" || hemi == "FALSE" )
          my_pheno_reader_info.set_allow_hemizygote(false);
        else
          errors << priority(error)
                 << "Invalid value for 'allow_hemizygote' specified.  Setting to 'false'..." << endl;
      }
    }
  }

  if( !my_pheno_reader_info.get_covariate_allele().size() )
    my_pheno_reader_info.set_covariate_moi("");

#if 0
  cout << "allele_delimiter = " << my_pheno_reader_info.get_allele_delimiter() << endl
       << "allele_missing   = " << my_pheno_reader_info.get_allele_missing() << endl
       << "allele_frequency = " << my_pheno_reader_info.get_allele_adjustment() << ", "
       << my_pheno_reader_info.get_min_allele_freq() << ", "
       << my_pheno_reader_info.get_max_allele_freq() << endl
       << "covariate_moi    = " << my_pheno_reader_info.get_covariate_moi() << endl
       << "covariate_allele = " << my_pheno_reader_info.get_covariate_allele() << endl
       << "allow_hemizygote = " << my_pheno_reader_info.get_allow_hemizygote() << endl;
#endif
}

void
RefPedInfo::remove_trait(size_t t_id)
{
  assert(t_id < my_traits.size());
  my_traits.erase(my_traits.begin() + t_id);

  --my_trait_count;
}

void
RefPedInfo::remove_marker(size_t m_id)
{
  assert(m_id < my_markers.size());
  my_markers[m_id] = my_markers[my_markers.size()-1];

  my_markers.pop_back();

  --my_marker_count;
}


// Returns: 0 - trait set ok
//          1 - trait set ok, but missing
//          2 - bad trait value, assumed missing
//          3 - invalid invididual or trait id
//          4 - trait not set due to trait settings
int
RefPedInfo::set_trait(size_t i, size_t t, const std::string & value, RefTraitInfo & trait_info)
{
  // Set up return value:
  int code = 0;

  // Trait # or member # out of range:
  if(t >= trait_count() || i >= member_count())
  {
    code = 3;
  }

  // Invalid trait type:
  else if(trait_info.type() == RefTraitInfo::invalid_trait)
  {
    code = 4;
  }

  // Categorical trait:
  else if(trait_info.type() == RefTraitInfo::categorical_trait)
  {
    if(value == trait_info.string_missing_code())
    {
      code = set_trait(i, t, numeric_limits<double>::quiet_NaN()) ? 0 : 3;
    }
    else if(trait_info.get_lockout() && (trait_info.get_id(value) == (size_t)-1)) // Lockout turned on and it's not in the list!
    {
      code = 2;
    }
    else // No lockout, or lockout's on and it IS in the list:
    {
      code = set_trait(i, t, trait_info.add_category(value)) ? 0 : 3;
    }
  }

  // Continuous or binary trait
  else
  {
    double d = str2doub(value);

    if(!finite(d))
      code = 2;

    std::string smiss  = trait_info.string_missing_code  ();
    double      nmiss  = trait_info.numeric_missing_code ();

    if( value == smiss || (finite(nmiss) && d == nmiss))
    {
      code = 1;
      d    = numeric_limits<double>::quiet_NaN();
    }
    else if( trait_info.type() == RefTraitInfo::binary_trait )
    {
      double thresh = trait_info.threshold();
        
      if(value == trait_info.string_affected_code() )
      {    
        code = 0;
        d    = 1.0;
      }
      else if( value == trait_info.string_unaffected_code() )
      {
        code = 0;
        d    = 0.0;
      }
      else if( finite(d) && d == trait_info.numeric_affected_code() )
      {    
        code = 0;
        d = 1.0;
      }
      else if( finite(d) && d == trait_info.numeric_unaffected_code() )
      {
        code = 0;
        d = 0.0;
      }
      else if( finite(d) && finite(thresh) )
      {
        code = 0;
        d    = (d > thresh) ? 1.0 : 0.0;
      }
      else
      {
        code = 2;
        d = numeric_limits<double>::quiet_NaN();
      }
    }

    if(!set_trait(i, t, d))
    {
      code = 3;
    }
  }


  return code;
}

inline
std::string compose_genotype(const std::string& value1,
                             const std::string& value2,
                             RefMarkerInfo&     marker_info)
{
  std::string value = value1;

  if(value2.size())
    value  += marker_info.gmodel().unphased_separator() + value2;

  // Add any genotype information we need to parse the genotype.  Does
  // nothing if inappropriate
  marker_info.add_genotype_dynamically(value);
  
  // Convert to canonical form, stripping internal white space, etc
  string canonical_value = marker_info.gmodel().parse_genotype_name(value);
  
  return canonical_value.size() ? canonical_value : value;
}

inline bool phenotype_not_found(uint phenotype)
{
  return phenotype == MLOCUS::NPOS;
}
inline bool allele_invalid(MLOCUS::allele a1)
{
  return !a1.is_valid();
}

/// Given an allele string, adds the allele to the marker if it isn't already
/// present with an allele frequency of 0.0
inline void add_marker_allele(RefMarkerInfo& marker_info, std::string allele)
{
  const MLOCUS::genotype_model& gmodel = marker_info.gmodel();
  
  if( allele.size()                          &&
      allele != gmodel.missing_allele_name() && 
      allele_invalid(gmodel.get_allele(allele)) )
  {
    marker_info.add_allele(allele, 0.0, true, true);
  }
}

/// Returns \c true if the marker allows expansion alleles, \c false otherwise.
///
inline bool marker_is_expandable(RefMarkerInfo& marker_info)
{
  return marker_info.gmodel().dynamic_alleles() && marker_info.codominant();
}

/// If the marker is expandable (allows new alleles), adds alleles which may be
/// missing
inline void expand_marker(RefMarkerInfo& marker_info, std::string phenotype)
{
    uint pid = marker_info.get_phenotype_id(phenotype);
    
    if( phenotype_not_found(pid) && marker_is_expandable(marker_info))
    {
      string   allele1;
      string   allele2;
      MLOCUS::Ordering order;
 
      marker_info.gmodel().parse_genotype_name(phenotype, allele1, allele2, order);

      add_marker_allele(marker_info, allele1);
      add_marker_allele(marker_info, allele2);
    }
}
// Returns: 0 - marker set ok
//          1 - marker set ok, but missing
//          2 - bad marker value, assumed missing
//          3 - invalid invididual or marker id
int
RefPedInfo::set_autosomal_phenotype(size_t i, size_t m, const std::string& value1,
                        const std::string& value2, RefMarkerInfo& marker_info)
{
  int code = 0;

  std::string value = compose_genotype(value1, value2, marker_info);

  uint p = marker_info.get_phenotype_id(value);
  
  if( value.empty() || p == marker_info.get_missing_phenotype_id())                                  // If value empty, marker considered missing.
  {
    set_phenotype(i, m, marker_info.get_missing_phenotype_id());
    
    return 1;
  }

  if( phenotype_not_found(p) || !set_phenotype(i, m, p))
  {
    code = 2;
  }

  return code;
}

// Returns: 0 - marker set ok
//          1 - marker set ok, but missing
//          2 - bad marker value, assumed missing
//          3 - invalid invididual or marker id
//          4 - sex based marker, but individual missing sex
//          5 - phenotype incompatible with sex of individual
int
RefPedInfo::set_x_linked_phenotype
  (size_t i, MPED::SexCode s, size_t m, const std::string& value1,
   const std::string& value2, RefMarkerInfo& marker_info)
{
  std::string value = compose_genotype(value1, value2, marker_info);
  
  if(value.empty()) // If value empty, marker considered missing.
  {
    set_phenotype(i, m, marker_info.get_missing_phenotype_id());
    
    return 1;
  }
  
  // Look up the phenotype id, if we can
  uint p = MLOCUS::NPOS;
  
  if(MPED::is_male(s))
    p = marker_info.get_phenotype_id(value + "(male)");
  
  if(p == MLOCUS::NPOS)
    p = marker_info.get_phenotype_id(value);

  if(p == marker_info.get_missing_phenotype_id())
  {
    set_phenotype(i, m, marker_info.get_missing_phenotype_id());
    
    return 1;
  }

  if(MPED::is_sex_unknown(s))
    return 4;
  
  if( phenotype_not_found(p) || !set_phenotype(i, m, p))
  {
    return 2;
  }

  // Check validity of the phenotype given the sex
  bool is_valid = false;
  
  for(RefMarkerInfo::phased_penetrance_iterator it = marker_info.phased_penetrance_begin(p);
      !is_valid && it != marker_info.phased_penetrance_end(p); ++it)
  {
    if((MPED::is_male(s)   && it.phased_geno().is_male_compatible()) ||
       (MPED::is_female(s) && it.phased_geno().is_female_compatible()) )
      is_valid = true;
  }
  for(RefMarkerInfo::unphased_penetrance_iterator it = marker_info.unphased_penetrance_begin(p);
      !is_valid && it != marker_info.unphased_penetrance_end(p); ++it)
  {
    if((MPED::is_male(s)   && it.unphased_geno().is_male_compatible()) ||
       (MPED::is_female(s) && it.unphased_geno().is_female_compatible()) )
      is_valid = true;
  }
  
  if(!is_valid)
  {
    set_phenotype(i, m, marker_info.get_missing_phenotype_id());
    
    return 5;
  }
  
  return 0;
}

// Returns: 0 - marker set ok
//          1 - marker set ok, but missing
//          2 - bad marker value, assumed missing
//          3 - invalid invididual or marker id
//          4 - sex based marker, but individual missing sex
//          5 - phenotype incompatible with sex of individual

int
RefPedInfo::set_y_linked_phenotype
  (size_t i, MPED::SexCode s, size_t m, const std::string& value1,
   const std::string& value2, RefMarkerInfo& marker_info)
{
  std::string value = compose_genotype(value1, value2, marker_info);

  if(MPED::is_female(s))
  {
    uint female_id = marker_info.get_phenotype_id("~X/~X");

    set_phenotype(i, m, female_id);

    if(value.empty()) // If value empty, marker considered missing.
    {
      return 1;
    }
    
    uint p = marker_info.get_phenotype_id(value);

    if(p == marker_info.get_missing_phenotype_id())
    {
      return 1;
    }
    
    if(phenotype_not_found(p))
      return 2;
    
    if(p != female_id)
    {
      return 5;
    }
    return 0;
  }
  
  // Male
  uint p = marker_info.get_phenotype_id(value);
  
  if(value.empty()) // If value empty, marker considered missing.
  {
    return 1;
  }

  if(p == marker_info.get_missing_phenotype_id())
  {
    set_phenotype(i, m, marker_info.get_missing_phenotype_id());
    
    return 1;
  }

  if(MPED::is_sex_unknown(s))
    return 4;
  
  if( phenotype_not_found(p) || !set_phenotype(i, m, p))
  {
    return 2;
  }

  if(p == marker_info.get_phenotype_id("~X/~X"))
  {
    set_phenotype(i, m, marker_info.get_missing_phenotype_id());

    return 5;
  }
  
  return 0;
}


// Returns: 0 - marker set ok
//          1 - marker set ok, but missing
//          2 - bad marker value, assumed missing
//          3 - invalid invididual or marker id
//          4 - sex based marker, but individual missing sex
//          5 - phenotype incompatible with sex of individual
int
RefPedInfo::set_phenotype(size_t i, size_t m, const std::string& value1,
                          const std::string& value2, RefMarkerInfo& marker_info)
{
  if( m >= marker_count() || i >= member_count() )
    return 3;

  return set_autosomal_phenotype(i,m,value1,value2, marker_info);
}

// Added for X & Y-linkage  - yjs  Mar. 2002
// Adjust the data according to the sex for X-linked & Y-linked marker.
//
//   X-linked marker:
//     female - no adjustment needed, missing allele is missing allele.
//       male - if one is missing, the missing allele is replaced with pseudo-allele "~Y".
//
//   Y-linked marker:
//     female - if allele exist, should be reported as error.
//              if one missing, replaced with "~X" to be reported as error.
//       male - if one missing, the missing allele is replaced with pseudo-allele "~X".
//
// Returns: 0 - marker set ok
//          1 - marker set ok, but missing
//          2 - bad marker value, assumed missing
//          3 - invalid invididual or marker id
//          4 - sex based marker, but individual missing sex
//          5 - phenotype incompatible with sex of individual

int
RefPedInfo::set_phenotype(size_t i, MPED::SexCode s, size_t m, const std::string& value1,
                          const std::string& value2, RefMarkerInfo& marker_info)
{
  if( m >= marker_count() || i >= member_count() )
    return 3;
  
  switch(marker_info.get_model_type())
  {
    case MLOCUS::AUTOSOMAL : return set_autosomal_phenotype(i,m,value1,value2, marker_info);
    case MLOCUS::X_LINKED  : return set_x_linked_phenotype(i,s,m,value1,value2, marker_info);
    case MLOCUS::Y_LINKED  : return set_y_linked_phenotype(i,s,m,value1,value2, marker_info);
  }
  
  return 3;
}

// Sort a reference pedigree into topological order based on lineage
void PedigreeSort(RefPedigree& p)
{
  int p1, p2, sp1, sp2;
  const RefPedigree::member_type *mid;

  RefPedInfo& info = p.info();

  for(size_t count = 0; count < p.member_count(); ++count)
  {
    mid = &p.member_index(count);
    p1 = p2 = -1;
    if(mid->parent1())
      p1 = mid->parent1()->index();
    if(mid->parent2())
      p2 = mid->parent2()->index();

    /*
    if(p1 == mid->index() || p2 == mid->index() )
    {
      cerr << "Your pedigree is broken.  Please fix it." << endl;
      cerr << "Indiviual " << mid->name() << " in pedigree " 
           << mid->pedigree()->name() << " is his/her own "
           << "parent (" << mid->parent1()->name() << "," 
                         << mid->parent2()->name() << ")" << endl;
      exit(1);
    }
    */
    
    // - Generalization of above check to insure that a person is not his/her own
    //   ancestor.  -djb  8/1/3
    //
    set<size_t>  ancestors;
    const RefPedigree::member_type*  member = ancestor_error(ancestors, mid);
    if(member)
    {
      cerr << "Problem detected in pedigree '" << mid->pedigree()->name() << "'."
           << "  One or more members is his/her own ancestor.\n"
           << "Member '" << member->name() << "' is involved.  Please "
           << "fix this problem and rerun program." << endl;
      exit(1);
    }

    DEBUG_PEDSORT(std::cout << "i=" << mid->index() << ", p=(" << p1 << "," << p2 << ")" << std::endl;)

    if(p1 == -1 && p2 == -1) continue;

    size_t maxp = max(p1,p2);

    if(mid->index() < maxp)
    {
      info.swap_members(count, maxp);
      p.member_index_swap(count, maxp);
      --count;
    }
  }

  // Need to sort subpedigrees too.
  RefMultiPedigree::subpedigree_iterator si = p.subpedigree_begin();
  for( ; si != p.subpedigree_end(); ++si )
  {
    RefMultiPedigree::member_iterator mi = si->member_begin();
    for( ; mi != si->member_end(); ++mi )
    {
      sp1 = sp2 = -1;
      if( mi->parent1() )
        sp1 = mi->parent1()->subindex();
      if( mi->parent2() )
        sp2 = mi->parent2()->subindex();

      if( sp1 == -1 && sp2 == -1 ) continue;

      size_t maxsp = max(sp1,sp2);

      if( mi->subindex() < maxsp )
      {
        mi->subpedigree()->member_index_swap(mi->subindex(), maxsp);
        --mi;
      }
    }
  }
}


// - If individual is in ancestor set, return pointer to that individual 
//   otherwise, return 0.  Check ancestry via parent 1.  Check ancestry via parent 2.
//
const RefPedigree::member_type*
ancestor_error(std::set<size_t>& ancestors, const RefPedigree::member_type* ind)
{
  if(! ind)
  {
    return  0;
  }
  
  if(ancestors.find(ind->index()) != ancestors.end())
  {
    return  ind;
  }
  
  ancestors.insert(ind->index());
  const RefPedigree::member_type*  p1_error = ancestor_error(ancestors, ind->parent1());
  
  if(p1_error)
  {
    ancestors.erase(ind->index());
    return  p1_error;
  }
  else
  {
    const RefPedigree::member_type*  p2_error = ancestor_error(ancestors, ind->parent2());
    ancestors.erase(ind->index());  
    return  p2_error;
  }
}

class PedigreeSorter // BUGGY!
{
public:
  PedigreeSorter();
  PedigreeSorter(RefPedigree& p);
  void sort(RefPedigree& p);

private:
  void             clear();
  void            visit(RefPedigree::member_pointer m);
  void      emit_members(RefPedigree& p);

  std::vector<size_t> my_new_order;
  size_t              my_next_index;
};


PedigreeSorter::PedigreeSorter()
{
  clear();
}

PedigreeSorter::PedigreeSorter(RefPedigree& p)
{
  sort(p);
}

void
PedigreeSorter::clear()
{
  my_next_index = 0;
  my_new_order.resize(0);
}

void
PedigreeSorter::sort(RefPedigree& p)
{
  if( p.member_count() < 2 || !p.member_count() )
    return;

  size_t member_count = p.member_count(); 

  my_next_index = 0;

  my_new_order.resize(0);
  my_new_order.resize(member_count, (size_t)-1 );

  for( size_t i = 0 ; i < member_count; ++i )
    visit( &p.member_index(i) );  

  emit_members(p);
  clear();

  CheckPedigreeSort(p);
}

void PedigreeSorter::visit(RefPedigree::member_pointer m)
{
  if( !m )
    return;

  size_t index = m->index();

  if( my_new_order[index] != (size_t)-1 )
    return;

  --my_new_order[index];

  visit( m->parent1() );
  visit( m->parent2() );

  my_new_order[index] = my_next_index++;
}

void PedigreeSorter::emit_members(RefPedigree& p)
{
  for(size_t i = 0; i < p.member_count(); ++i)
  {
    size_t old_index = i;
    size_t new_index = my_new_order[i];

    if( new_index == (size_t)-1 )
      continue;

    cout << "emitting " << p.member_index(old_index).name() << " to " 
         << new_index   << endl;

    p.info().swap_members(new_index, old_index);
    p.member_index_swap(  new_index, old_index);

    my_new_order[i] = (size_t)-1;
  }
}

void BuggyPedigreeSort(RefPedigree& p)
{
  PedigreeSorter sorter;
  sorter.sort(p);
}

bool CheckPedigreeSort(RefPedigree& p)
{
  for(size_t i = 0; i < p.member_count(); ++i)
  {
    RefPedigree::member_pointer m = &p.member_index(i);
    size_t p1 = m->parent1()? m->parent1()->index() : 0;
    size_t p2 = m->parent2()? m->parent2()->index() : 0;
    if( i < p1 || i < p2 )
    {
      cout << "Bad order!" << endl;
      return false;
    }
  }
  return true;
}

size_t 
RefTraitInfo::get_id(const std::string & category) const
{
  for(size_t i = 0; i < my_categories.size(); ++i)
    if(my_categories[i] == category)
      return i;
      
  return (size_t)-1;
}
    

} // end namespace RPED
} // end namespace SAGE

