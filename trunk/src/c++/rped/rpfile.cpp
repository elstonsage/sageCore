#include <string>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <iomanip>
#include "mped/sp.h"
#include "rped/rpfile.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "error/bufferederrorstream.h"

#define DEBUG_RPEDFILE(x)

namespace SAGE {
namespace RPED {

//==================================================================
//
//  CONSTRUCTOR
//
//==================================================================
RefPedigreeFile::RefPedigreeFile(cerrorstream &err) : errors(err)
{
  set_treat_as_sibs          (false);
  set_verbose_output         (10);
  set_reject_partial_lineage ();
  set_require_record         (false);
  set_skip_traits            (false);
  set_skip_markers           (false);
  set_dynamic_markers        (false);
  set_sex_linked_exist       (false);
  set_no_sex_field           (false);
  set_no_sex_ok_option       (false);
  set_sex_code_trait         (false);
  set_pedigree_id_trait      (false);
  set_sex_field_name         ("");
  set_pedigree_id_name       ("");
  reset_counts               ();
  validate                   ();

  my_inds.clear();
}

//==================================================================   
//
//  DESTRUCTOR
//
//==================================================================
RefPedigreeFile::~RefPedigreeFile() 
{
  errors.flush();
}

//==================================================================   
//
//  reset_counts()
//
//==================================================================
void 
RefPedigreeFile::reset_counts()
{
  my_skip_count           = 0;
  my_study_id_count       = 0;
  my_pedigree_id_count    = 0;
  my_individual_id_count  = 0;
  my_parent_id_count      = 0;
  my_sex_count            = 0;
  my_marker_count         = 0;
  my_invalid_marker_count = 0;
  my_marker_cov_count         = 0;
  my_invalid_marker_cov_count = 0;
  my_trait_count          = 0;
  my_invalid_trait_count  = 0;
  my_string_count         = 0;
  my_invalid_string_count = 0;
}

//==================================================================   
//
//  input(...)
//
//==================================================================
bool 
RefPedigreeFile::input(RefMultiPedigree &p, const string &filename, ostream &messages)
{
  // Read in pedigree meta-information:
  if(!input_pedigree(p, filename, messages, false))
    return false;

  // Read in actual pedigree data:
  if(!input_data(p, filename, messages, true))
    return false;

  if(sex_linked_exist())
    update_sex_linked_marker_info(p);

  return do_no_sex_structural_test(p);
}

//==================================================================   
//
//  print_mped(...)
//
//==================================================================
void 
RefPedigreeFile::print_mped(
	const RefMultiPedigree & p, 
	const string           & filename, 
	      ostream          & messages,
	      bool               dump_trait, 
	      bool               dump_marker)
{
  const RefMPedInfo & mped_info = p.info();

  print_family_structure_header(messages, filename);

  // Loop through the first X individuals and print out their info:
  for(size_t i = 0; i < my_ind_list.size(); ++i)
  {
    std::string              ped_name = my_ind_list[i].first,
                             ind_name = my_ind_list[i].second;
    RPED::MemberConstPointer ind      = p.member_find(ped_name, ind_name);
    std::string              p1       = ind->parent1() ? ind->parent1()->name() : mped_info.individual_missing_code(),
                             p2       = ind->parent2() ? ind->parent2()->name() : mped_info.individual_missing_code(),
                             ind_sex  = "";

    if( ind->is_male() )
      ind_sex = mped_info.sex_code_male    ();
    else if( ind->is_female() )
      ind_sex = mped_info.sex_code_female  ();
    else
      ind_sex = mped_info.sex_code_unknown ();

    print_family_structure(messages, ped_name, ind_name, ind_sex, p1, p2);

  } // End loop across first X individuals

  print_family_structure_footer(messages);

  if( dump_trait && (mped_info.trait_count() || mped_info.string_count()) )
  {
    print_trait_header(messages, mped_info, filename);

    for( size_t i = 0; i < my_ind_list.size(); ++i )
    {
      std::string pn = my_ind_list[i].first,
                  id = my_ind_list[i].second;

      RefMultiPedigree::member_const_pointer mem = p.member_find(pn, id);

      const RefPedInfo &ped_info = mem->pedigree()->info();

      std::vector<std::pair<size_t,std::string> > trait_values;

      size_t sti = mped_info.trait_find("SEX_CODE");
      if( sti < mped_info.trait_count() )
      {
        std::string t_value = mped_info.trait_info(sti).string_missing_code();

        if( !ped_info.trait_missing(mem->index(), sti) )
          t_value = doub2str(ped_info.trait(mem->index(), sti));

        trait_values.push_back(make_pair(sti, t_value));
      }

      size_t t1 = mped_info.trait_find("FAMILIAL_INDICATOR");
      size_t t2 = mped_info.trait_find("FOUNDER_INDICATOR");
      size_t t3 = mped_info.trait_find("PEDIGREE_SIZE");

      if( t1 < mped_info.trait_count() )
      {
        std::string t_value = mped_info.trait_info(t1).string_missing_code();

        if( !ped_info.trait_missing(mem->index(), t1) )
          t_value = doub2str(ped_info.trait(mem->index(), t1));

        trait_values.push_back(make_pair(t1, t_value));
      }

      if( t2 < mped_info.trait_count() )
      {
        std::string t_value = mped_info.trait_info(t2).string_missing_code();

        if( !ped_info.trait_missing(mem->index(), t2) )
          t_value = doub2str(ped_info.trait(mem->index(), t2));

        trait_values.push_back(make_pair(t2, t_value));
      }

      if( t3 < mped_info.trait_count() )
      {
        std::string t_value = mped_info.trait_info(t3).string_missing_code();

        if( !ped_info.trait_missing(mem->index(), t3) )
          t_value = doub2str(ped_info.trait(mem->index(), t3));

        trait_values.push_back(make_pair(t3, t_value));
      }

      set<size_t> mcov_added;
      
      for(field_list_type::const_iterator field_info = my_fields.begin(); field_info != my_fields.end(); ++field_info )
      {
        if( !(field_info->type == trait || field_info->type == marker_cov || field_info->type == allele_cov) )
          continue;

        if( field_info->type != trait )
        {
          if( mcov_added.find(field_info->index) != mcov_added.end() )
            continue;
          else
            mcov_added.insert(field_info->index);
        }

        if( field_info->index < mped_info.trait_count() )
        {
          const RefTraitInfo & tinfo   = mped_info.trait_info(field_info->index);
                std::string    t_value = tinfo.string_missing_code();
                double         val     = ped_info.trait(mem->index(), field_info->index);

          if( !ped_info.trait_missing(mem->index(), field_info->index) )
          {
            t_value = tinfo.type() == RefTraitInfo::categorical_trait ? tinfo.get_categories()[(size_t)val] : doub2str(val);
          }

          trait_values.push_back(make_pair(field_info->index, t_value));
        }
      }

      std::vector<std::pair<size_t,std::string> > string_values;

      for( field_list_type::const_iterator field_info = my_fields.begin(); field_info != my_fields.end(); ++field_info )
      {
        if( field_info->type == string_field && field_info->index < mped_info.string_count() )
        {
          string_values.push_back(make_pair(field_info->index, ped_info.get_string(mem->index(), field_info->index)));
        }
      }

      print_trait(messages, mped_info, pn, id, trait_values, string_values);
    }

    print_trait_footer(messages);
  }

  if(dump_marker && mped_info.marker_count())
  {
    print_marker_header(messages, mped_info, filename);

    for( size_t i = 0; i < my_ind_list.size(); ++i )
    {
      std::string pn = my_ind_list[i].first;
      std::string id = my_ind_list[i].second;

      RefMultiPedigree::member_const_pointer mem = p.member_find(pn, id);

      const RefPedInfo &ped_info = mem->pedigree()->info();

      std::vector<std::pair<size_t, std::string> >  marker_values;
      set<size_t> marker_added;

      for(field_list_type::const_iterator field_info = my_fields.begin() ; field_info != my_fields.end(); ++field_info)
      {
        if( field_info->type != marker && field_info->type != allele )
          continue;

        if( marker_added.find(field_info->index) != marker_added.end() )
          continue;

        if( field_info->index < mped_info.marker_count())
        {
          marker_added.insert(field_info->index);
          size_t m_index = field_info->index;

          const RefMarkerInfo& minfo = mped_info.marker_info(m_index);

          uint pheno_id = minfo.get_missing_phenotype_id();

          if( !ped_info.phenotype_missing(mem->index(), m_index, minfo) )
            pheno_id = ped_info.phenotype(mem->index(), m_index);

          std::string pheno_name = minfo.get_phenotype(pheno_id).name();

          if( pheno_id == minfo.get_missing_phenotype_id() )
          {
            std::string a = minfo.missing_allele_name();
            if( a == "*missing" )
              a = "?";
            pheno_name = a + minfo.gmodel().unphased_separator() + a;
          }
          else if( minfo.codominant(pheno_id) )
          {
            std::string al1, al2;
            MLOCUS::Ordering order;

            std::string geno_name = minfo.unphased_penetrance_begin(pheno_id).unphased_geno().name();

            minfo.gmodel().parse_genotype_name(geno_name, al1, al2, order);

            if( !al1.size() || al1.substr(0,1) == "~" ) al1 = minfo.missing_allele_name();
            if( !al2.size() || al2.substr(0,1) == "~" ) al2 = minfo.missing_allele_name();

            if( al1 == "*missing" )
              al1 = "?";

            if( al2 == "*missing" )
              al2 = "?";

            pheno_name = al1 + minfo.gmodel().unphased_separator() + al2;
          }

          marker_values.push_back(make_pair(m_index, pheno_name));
        }
      }

      print_marker(messages, mped_info, pn, id, marker_values);
    }

    print_marker_footer(messages);
  }
}

//==================================================================   
//
//  print_family_structure_header(...)
//
//==================================================================
void 
RefPedigreeFile::print_family_structure_header(ostream &messages, const std::string &filename) const
{
  messages << std::endl 
           << "Family structure information on the first " 
           << verbose_output() 
           << " individuals read from file: " 
           << filename
           << std::endl 
           << std::endl
           << "     PED. ID       IND. ID       SEX       PARENT1       PARENT2     " << std::endl
           << "     ------------  ------------  --------  ------------  ------------" << std::endl;
}

//==================================================================   
//
//  print_family_structure(...)
//
//==================================================================
void 
RefPedigreeFile::print_family_structure(
	      ostream & messages,
        const std::string  & pn, 
	const std::string  & id,
        const std::string  & sex, 
        const std::string  & p1, 
	const std::string  & p2) const
{
  messages << left
           << "     " << setw(12) << pn
           << "  "    << setw(12) << id 
           << "  "    << setw(8)  << sex  
           << "  "    << setw(12) << p1  
           << "  "    << setw(12) << p2 
           << std::endl;
}

//==================================================================   
//
//  print_family_structure_footer(...)
//
//==================================================================
void 
RefPedigreeFile::print_family_structure_footer(ostream &messages) const
{
  messages << std::endl;
}

//==================================================================   
//
//  print_trait_header(...)
//
//==================================================================
void 
RefPedigreeFile::print_trait_header(ostream &messages, const RefMPedInfo &mped_info, const std::string &filename) const
{
  messages << std::endl
           << "Phenotypes for the first "
           << verbose_output()
           << " individuals read from file: " 
           << filename
           << std::endl
           << std::endl
           << "     PED. ID       IND. ID     " << left;

  size_t num_trait_printed  = 0;
  size_t num_string_printed = 0;

  size_t sti = mped_info.trait_find("SEX_CODE");
  if( sti < mped_info.trait_count() )
  {
    messages << "  " << setw(20) << mped_info.trait_info(sti).name();
    ++num_trait_printed;
  }

  size_t t1 = mped_info.trait_find("FAMILIAL_INDICATOR");
  size_t t2 = mped_info.trait_find("FOUNDER_INDICATOR");
  size_t t3 = mped_info.trait_find("PEDIGREE_SIZE");

  if( t1 < mped_info.trait_count() )
  {
    messages << "  " << setw(20) << mped_info.trait_info(t1).name();
    ++num_trait_printed;
  }

  if( t2 < mped_info.trait_count() )
  {
    messages << "  " << setw(20) << mped_info.trait_info(t2).name();
    ++num_trait_printed;
  }

  if( t3 < mped_info.trait_count() )
  {
    messages << "  " << setw(20) << mped_info.trait_info(t3).name();
    ++num_trait_printed;
  }

  set<size_t> mcov_added;

  for(field_list_type::const_iterator field_info = my_fields.begin() ; field_info != my_fields.end(); ++field_info)
  {
    if( !(field_info->type == trait || field_info->type == marker_cov || field_info->type == allele_cov) )
      continue;

    if( field_info->type != trait )
    {
      if( mcov_added.find(field_info->index) != mcov_added.end() )
        continue;
      else
        mcov_added.insert(field_info->index);
    }

    if( field_info->index < mped_info.trait_count() )
    {
      messages << "  " << setw(20) << mped_info.trait_info(field_info->index).name();
      ++num_trait_printed;
    }
  }

  for(field_list_type::const_iterator field_info = my_fields.begin(); field_info != my_fields.end() ; ++field_info)
  {
    if(field_info->type == string_field && field_info->index < mped_info.string_count())
    {
      messages << "  " << setw(20) << mped_info.string_info(field_info->index).name();
      ++num_string_printed;
    }
  }

  messages << std::endl << "     ------------  ------------";

  for( size_t i = 0; i < num_trait_printed + num_string_printed; ++i )
    messages << "  --------------------";

  messages << std::endl;
}

//==================================================================   
//
//  print_trait(...)
//
//==================================================================
void 
RefPedigreeFile::print_trait(
        ostream                                     & messages, 
  const RefMPedInfo                                 & mped_info,
  const std::string                                 & pn, 
  const std::string                                 & id,
  const std::vector<std::pair<size_t,std::string> > & trait_values,
  const std::vector<std::pair<size_t,std::string> > & string_values) const
{
  messages << "     " << setw(12) << pn << "  " << setw(12) << id;

  for(size_t i=0; i < trait_values.size(); ++i)
  {
    if( trait_values[i].first < mped_info.trait_count())
    {
      messages << "  " << setw(20) << trait_values[i].second;
    }
  }

  for( size_t i=0; i < string_values.size(); ++i)
  {
    if( string_values[i].first < mped_info.string_count())
    {
      messages << "  " << setw(20) << string_values[i].second;
    }
  }

  messages << std::endl;
}

//==================================================================   
//
//  print_trait_footer(...)
//
//==================================================================
void 
RefPedigreeFile::print_trait_footer(ostream &messages) const
{
  messages << std::endl;
}

//==================================================================   
//
//  print_marker_header(...)
//
//==================================================================
void 
RefPedigreeFile::print_marker_header(ostream &messages, const RefMPedInfo &mped_info, const std::string &filename) const
{
  messages << std::endl
           << "Markers for the first " 
           << verbose_output()
           << " individuals read from file: " 
           << filename
           << std::endl 
           << std::endl
           << "     PED. ID       IND. ID     " << left;

  size_t num_marker_printed  = 0;

  set<size_t> marker_added;

  for(field_list_type::const_iterator field_info = my_fields.begin() ; field_info != my_fields.end(); ++field_info)
  {
    if((field_info->type == marker || field_info->type == allele)   &&
       (marker_added.find(field_info->index) == marker_added.end()) &&
       (field_info->index < mped_info.marker_count()))
    {
      marker_added.insert(field_info->index);

      messages << "  " << setw(12) << mped_info.marker_info(field_info->index).name();

      ++num_marker_printed;
    }
  }

  messages << std::endl << "     ------------  ------------";

  for( size_t i = 0; i < num_marker_printed; ++i )
    messages << "  ------------";

  messages << std::endl;
}

//==================================================================   
//
//  print_marker(...)
//
//==================================================================
void 
RefPedigreeFile::print_marker(
        std::ostream                                 & messages, 
  const RefMPedInfo                                  & mped_info,
  const std::string                                  & pn, 
  const std::string                                  & id,
  const std::vector<std::pair<size_t, std::string> > & marker_values) const
{
  messages << "     " << setw(12) << pn << "  " << setw(12) << id;

  for(size_t m = 0; m < marker_values.size(); ++m)
  {
    if(marker_values[m].first < mped_info.marker_count())
    {
      messages << "  " << setw(12) << marker_values[m].second;
    }
  }

  messages << std::endl;
}

//==================================================================   
//
//  print_marker_footer(...)
//
//==================================================================
void 
RefPedigreeFile::print_marker_footer(ostream &messages) const
{
  messages << std::endl;
}

//==================================================================   
//
//  parse_dynamic_markers_missing(...)
//
//==================================================================
void 
RefPedigreeFile::parse_dynamic_markers_missing(RefMarkerInfo& marker, const std::string& missing)
{
  char        sep0 = marker.gmodel().separators()[0];
  std::string mv   = missing;

  if( strchr(mv.c_str(), sep0) )
  {
    std::string mv2 = strip_ws(mv.substr(mv.find_first_of(sep0), mv.size()), "/ \t");

    mv = strip_ws(mv.substr(0, mv.find_first_of(sep0)), " \t");

    if( mv != mv2 )
    {
      errors << priority(warning) 
             << "The Missing Value code must use the "
             << "same value for both alleles for '"
             << marker.name() 
             << "'.  Will use " 
             << mv 
             << "/" 
             << mv
             << " for the missing value." 
             << std::endl;
    }
  }

  if(mv != marker.missing_allele_name())
  {
    marker.gmodel().set_missing_allele_name(mv);
    marker.set_missing_phenotype_name(mv);

    std::string sep = marker.gmodel().separators();

    marker.alias_phenotype(marker.get_missing_phenotype_id(), mv + sep[0] + mv);
    marker.alias_phenotype(marker.get_missing_phenotype_id(), mv + sep[1] + mv);
    marker.alias_phenotype(marker.get_missing_phenotype_id(), mv + sep[2] + mv);
  }
}

//==================================================================   
//
//  validate_fields(...)
//
//==================================================================
bool 
RefPedigreeFile::validate_fields(bool data_only, bool quiet)
{
  // NOTE: We treat invalid markers and traits as non-fatal states

  reset_counts();

  typedef std::map<std::string, size_t> name_count_map;

  name_count_map trait_map,
                 string_map,
                 marker_map,
                 marker_cov_map;

  for(field_list_type::iterator i = my_fields.begin(); i != my_fields.end(); ++i)
  {
    switch( i->type )
    {
      case skip          : ++my_skip_count;                    break;
      case study_id      : ++my_study_id_count;                break;
      case pedigree_id   : ++my_pedigree_id_count;             break;
      case individual_id : ++my_individual_id_count;           break;
      case parent_id     : ++my_parent_id_count;               break;
      case sex_code      : ++my_sex_count;                     break;

      case trait         : trait_map  [toUpper(i->name)]++;    break;
      case string_field  : string_map [toUpper(i->name)]++;    break;
      case marker        : marker_map [toUpper(i->name)] += 2; break;
      case allele        : marker_map [toUpper(i->name)]++;    break;

      case marker_cov : marker_cov_map [toUpper(i->name)] += 2; break;
      case allele_cov : marker_cov_map [toUpper(i->name)]++;    break;
    }
  }

  // Set validity true to get started:
  bool valid_fields = true;

  // Check structural stuff:
  if(!data_only)
  {
    // Wrong count of any particular field:
    
    // Check study id count:
    if (my_study_id_count > 1 && !quiet)
    {
      errors << priority(warning) << "More than one study ID specified." << std::endl;
    }

    // Individual id:
    if(my_individual_id_count == 0)
    {
      valid_fields = false;
      if(!quiet)
        errors << priority(error) << "Individual ID field not found." << std::endl;
    }
    else if(my_individual_id_count > 1)
    {
      valid_fields = false;
      if(!quiet)
        errors << priority(error) << "Too many Individual ID fields found." << std::endl;
    }

    // Ped id:
    if(my_pedigree_id_count > 1)
    {
      valid_fields = false;
      if(!quiet)
        errors << priority(error) << "More than one pedigree id field indicated." << std::endl;
    }
  
    // Parent id:
    if(my_parent_id_count == 1)
    {
      valid_fields = false;
      if(!quiet)
        errors << priority(error) << "Only one parent id field indicated; either none or two are required." << std::endl;
    }
    else if(my_parent_id_count > 2)
    {
      valid_fields = false;
      if(!quiet)
        errors << priority(error) << "More than two parent id fields indicated." << std::endl;
    }
  
    // Parent id:
    if(my_sex_count > 1)
    {
      valid_fields = false;
      if(!quiet)
        errors << priority(error) << "More than one sex field indicated." << std::endl;
    }
  
    // Moving on to more complicated structural checking:  
    
    // Parents present, sex missing, but no_sex_field not set to true:
    if(my_parent_id_count == 2 && my_sex_count == 0 && my_no_sex_field == false)
    {
      valid_fields = false;
      if(!quiet)
        errors << priority(error) << "Parent ID's present, but sex field is missing and 'NO_SEX_FIELD' has not been specified." << std::endl;
    }
    
    // Parents present but treat_as_sibs specified:
    if(my_parent_id_count == 2 && get_treat_as_sibs() == true)
    {
      valid_fields = false;
      if(!quiet)
        errors << priority(error) << "User has specified 'TREAT_AS_SIBS', but still listed parent ids." << std::endl;
    }
    
    // Sex present but NO_SEX_FIELD also present:
    if(my_sex_count == 1 && my_no_sex_field == true)
    {
      valid_fields = false;
      if(!quiet)
        errors << priority(error) << "Both 'SEX_FIELD' and 'NO_SEX_FIELD' were specified in the parameter file." << std::endl;
    }
    
    // Sex missing, NO_SEX_FIELD specified, deliver an advance warning (but analysis still valid):
    if(my_sex_count == 0 && my_no_sex_field == true)
    {
      if(!quiet) 
        errors << priority(warning) << "Please note that no sex field is present in the data; "
                                    << "analyses that depend on sex-specific data may produce unpredictable results." 
                                    << std::endl;
    }
        
  } // End structural checking
  
  // Check traits:

  for(field_list_type::iterator i = my_fields.begin(); i != my_fields.end(); ++i)
  {
    if( i->type == trait )
    {
      size_t &count = trait_map[ toUpper(i->name) ];

      if( count == 1 )
      {
        ++my_trait_count;
      }
      else if( count != (size_t)-1 )  // Invalidate new bad trait
      {
        if(!quiet)
          errors << priority(warning) << "Trait assigned to more than one field.  Skipping trait '" << i->name << "'." << std::endl;

        // Count trait as bad and mark count to avoid further warnings
        ++my_invalid_trait_count;
        count = (size_t)-1;
        i->index = (size_t)-1;
      }
      else // invalidate known bad trait
      {
        i->index = (size_t)-1;
      }
    }
    else if( i->type == string_field )
    {
      size_t &count = string_map[ toUpper(i->name) ];

      if( count == 1 )
      {
        ++my_string_count;
      }
      else if( count != (size_t)-1 )  // Invalidate new bad string
      {
        if(!quiet)
          errors << priority(warning) << "String field assigned to more than one field.  Skipping field '" << i->name << "'." << std::endl;

        // Count field as bad and mark count to avoid further warnings
        ++my_invalid_string_count;
        count = (size_t)-1;
        i->index = (size_t)-1;
      }
      else // invalidate known bad field
      {
        i->index = (size_t)-1;
      }
    }
    else if( i->type == marker || i->type == allele )
    {
      size_t &count = marker_map[ toUpper(i->name) ];

      if( count == 2 )
      {
        ++my_marker_count;
      }
      else if( count != (size_t)-1 )  // Invalidate new bad marker
      {
        if(!quiet && count < 2)
        {
          errors << priority(warning) << "Marker assigned to too few allele fields.  Skipping marker '" << i->name << "'." << std::endl;
        }
        else if(!quiet && count > 2)
        {
          errors << priority(warning) << "Marker assigned to too many fields.  Skipping marker '" << i->name << "'." << std::endl;
        }

        // Count marker as bad and mark count to avoid further warnings
        ++my_invalid_marker_count;

        count    = (size_t)-1;
        i->index = (size_t)-1;
      }
      else // invalidate known bad marker
      {
        i->index = (size_t)-1;
      }
    }
    else if( i->type == marker_cov || i->type == allele_cov )
    {
      size_t &count = marker_cov_map[ toUpper(i->name) ];

      if( count == 2 )
      {
        ++my_marker_cov_count;
      }
      else if( count != (size_t)-1 )  // Invalidate new bad marker
      {
        if(!quiet && count < 2)
        {
          errors << priority(warning) << "Marker as covariate assigned to too few allele fields.  Skipping marker '" << i->name << "'." << std::endl;
        }
        else if(!quiet && count > 2)
        {
          errors << priority(warning) << "Marker as covariate assigned to too many fields.  Skipping marker '" << i->name << "'." << std::endl;
        }

        // Count marker as bad and mark count to avoid further warnings
        ++my_invalid_marker_cov_count;

        count    = (size_t)-1;
        i->index = (size_t)-1;
      }
      else // invalidate known bad marker
      {
        i->index = (size_t)-1;
      }
    }
  }

  // NOTE: We treat invalid markers and traits as non-fatal states

  return valid_fields;
}

//==================================================================   
//
//  build_pedigree(...)
//
//==================================================================
bool 
RefPedigreeFile::build_pedigree(RefMultiPedigree &p)
{
  const RefMPedInfo &mped_info = p.info();

  p.build();
  
  bool found_fatal_errors = false;

  // Construct the pedigrees and fill in index information
  // (much of this should/will be moved into the RefPedigree directly
  //  as a build action callback)
  RefMultiPedigree::pedigree_iterator i = p.pedigree_begin();
  for( int j=0; i != p.pedigree_end(); ++i, ++j)
  {
    if( i->error_count() )
    {
      if(!report_pedigree_build_errors(*i))
        found_fatal_errors = true;
    }
    
    if( i->warning_count() )
    {
      if(!report_pedigree_build_warnings(*i))
        found_fatal_errors = true;
    }

    // Create indices
    i->info().build(*i);
    i->info().resize_traits(  mped_info.trait_count() );
    i->info().resize_strings( mped_info.string_count() );
    i->info().resize_markers( mped_info.marker_count(), mped_info );
  }

  return !found_fatal_errors;
}

/// Reports structural errors detected when building the pedigree and a count
/// of the number of errors.  Returns a true/false based upon if any of the errors
/// are fatal (ie, unignorable).  Fatal errors currently include inconsistencies
/// in parental sex.
bool
RefPedigreeFile::report_pedigree_build_errors(const Pedigree &p) const
{
  bool found_fatal_error = false;
      
  // Create a buffer stream to hold the errors.  We want to report the error count
  // first, but some errors are suppressed in certain circumstances, so we need 
  // to count them first.  We can use a buffer to queue the messages.
  SAGE::bufferederrorstream<char> buf(errors);
  
  // Create a counter to count the reported errors
  size_t error_count = 0;
  
  RefMultiPedigree::pedigree_type::error_iterator err;
  for( err = p.error_begin(); err != p.error_end(); ++err )
  {
    switch( err->state )
    {
      case MPED::error_info::bad_sibship:
        buf << priority(warning)
            << "     Error in sibship due to individuals: '" 
            << err->name1 << "', '" << err->name2 << "'." << std::endl;
        ++error_count;
        break;

      case MPED::error_info::bad_marriage:
        buf << priority(warning)
            << "     Inconsistent mating due to individuals: '" 
            << err->name1 << "', '" << err->name2 << "'." << std::endl;
        ++error_count;
        break;

      case MPED::error_info::bad_lineage:
        buf << priority(warning)
            << "     Inconsistent lineage due to individuals: '" 
            << err->name2 << "', '" << err->name3 << "' -> '"
            << err->name1 << "'." << std::endl;
        ++error_count;
        break;

      case MPED::error_info::bad_gender:
        buf << priority(critical)
            << "     Individual '" << err->name1
            << "' is present as both male and female." << std::endl;
        found_fatal_error = true;
        ++error_count;
        break;

      case MPED::error_info::same_sex_marriage:
        buf << priority(critical)
            << "Both parents '" << err->name1 << "' and '" << err->name2
            << "' are labeled as or inferred to be " << err->name3 << "."
            << std::endl;
        found_fatal_error = true;
        ++error_count;
        break;
        
      case MPED::error_info::bad_marriage_loop:
        buf << priority(critical)
            << "Individual '" << err->name1 << "' belongs to a marriage loop "
            << "of odd length.  Such marriage loops cannot be assigned sex "
            << "consistently." << std::endl;
        found_fatal_error = true;
        ++error_count;
        break;
        
      case MPED::error_info::no_sex_parents:
        buf << priority(warning)
            << "Both parents '" << err->name1 << "' and '" << err->name2
            << "' have unknown sex and cannot be inferred." << std::endl;
        ++error_count;
        break;
        
      default: break;
    }
  }
  
  // If we have reported errors, report the number of errors first, then dump our
  // buffer
  if(error_count)
  {
    errors << priority(error) << "Error(s) building pedigree '" << p.name() << "': "
           << error_count << " errors." << std::endl;
    buf.flush_buffer();
  }
  
  return !found_fatal_error;
}

/// Reports structural warnings detected when building the pedigree.
/// Returns a true/false based upon if any of the errors
/// are fatal (ie, unignorable).  In the current version, no warning is fatal,
/// but the boolean is included for consistency.
bool
RefPedigreeFile::report_pedigree_build_warnings(const Pedigree &p) const
{
  RefMultiPedigree::pedigree_type::error_iterator err;
  for( err = p.warning_begin(); err != p.warning_end(); ++err )
  {
    switch( err->state )
    {
      case MPED::error_info::gender_inferred:
        errors << priority(warning)
               << "The sex of individual '" << err->name1
               << "' in pedigree '" << p.name()
               << "' is believed to be " << err->name2
               << ".  Sex code will be reset." << std::endl;
        break;

      default: break;
    }
  }
  
  return true;
}    

//==================================================================   
//
//  build_data(...)
//
//==================================================================
bool 
RefPedigreeFile::build_data(RefMultiPedigree &p)
{
  RefMPedInfo &mped_info = p.info();

  // fill in index information
  RefMultiPedigree::pedigree_iterator i = p.pedigree_begin();
  for( int j=0; i != p.pedigree_end(); ++i, ++j)
  {
    // Create indices
    i->info().build(*i);
    i->info().resize_traits(  mped_info.trait_count() );
    i->info().resize_strings( mped_info.string_count() );
    i->info().resize_markers( mped_info.marker_count(), mped_info );
  }

  // Build implicit traits
  if( sex_code_trait() )
  {
    size_t t = mped_info.trait_find("SEX_CODE");

    if( t < mped_info.trait_count() )
    {
      RefMultiPedigree::member_const_iterator j;

      for( i=p.pedigree_begin(); i != p.pedigree_end(); ++i )
      {
        RefMultiPedigree::pedinfo_type &info = i->info();

        for( j=i->member_begin(); j != i->member_end(); ++j )
        {
          int ind_num = j->index();

          if(j->is_male())
            info.set_trait(ind_num, t, mped_info.sex_code_male(), mped_info.trait_info(t));
          else if(j->is_female())
            info.set_trait(ind_num, t, mped_info.sex_code_female(), mped_info.trait_info(t));
        }
      }
    }
    ++my_trait_count;
  }

  if( pedigree_id_trait() )
  {
    size_t t1 = mped_info.trait_find("FAMILIAL_INDICATOR");
    size_t t2 = mped_info.trait_find("FOUNDER_INDICATOR");
    size_t t3 = mped_info.trait_find("PEDIGREE_SIZE");

    if(    t1 < mped_info.trait_count()
        && t2 < mped_info.trait_count()
        && t3 < mped_info.trait_count() )
    {
      RefMultiPedigree::member_const_iterator j;

      for( i=p.pedigree_begin(); i != p.pedigree_end(); ++i )
      {
        RefMultiPedigree::pedinfo_type &info = i->info();

        for( j=i->member_begin(); j != i->member_end(); ++j )
        {
          int ind_num = j->index();

          if( j->pedigree()->member_count() > 1 )
            info.set_trait(ind_num, t1, "1", mped_info.trait_info(t1));
          else
            info.set_trait(ind_num, t1, "0", mped_info.trait_info(t1));

          if( j->is_founder() )
            info.set_trait(ind_num, t2, "1", mped_info.trait_info(t2));
          else
            info.set_trait(ind_num, t2, "0", mped_info.trait_info(t2));

          info.set_trait(ind_num, t3, j->pedigree()->member_count());

        }
      }
    }

    my_trait_count += 3;
  }

  return true;
}

//==================================================================   
//
//  add_member(...)
//
//==================================================================
void 
RefPedigreeFile::add_member(
        RefMultiPedigree & multipedigree, 
  const std::string      & ped_name,
  const std::string      & ind_name, 
  const std::string      & sex,
  const std::string      & parent1, 
  const std::string      & parent2,
        size_t             line,
        size_t             count)
{
  if(count == (size_t)-1)
  {
     count = line;
  }
  
  if(!ped_name.size())
  {
    errors << priority(warning) << "[" << line << "] record is missing pedigree ID.  Skipping..." << std::endl;

    return;
  }

  if(ind_name == "" || ind_name == multipedigree.info().individual_missing_code())
  {
    errors << priority(warning) << "[" << line
           << "] Found an individual ID that is missing or"
           << " the same as the missing individual/parent code.  Skipping..."
           << std::endl;

    return;
  }
  
  if(reject_partial_lineage() && ((parent1 == multipedigree.info().individual_missing_code()) ^ (parent2 == multipedigree.info().individual_missing_code())))
  {
    errors << priority(warning) << "[" << line << "] Individual '" << ind_name << "' in pedigree '" << ped_name << "' has only one parent specified and will be skipped." << std::endl;

    return;
  }

  // Add entries for the given parents, assuming require_record is set to false (the default).
  if(!require_record())
  { 
    if(parent1 != multipedigree.info().individual_missing_code())
      multipedigree.add_member(ped_name, parent1, MPED::SEX_ARB);
      
    if(parent2 != multipedigree.info().individual_missing_code())
      multipedigree.add_member(ped_name, parent2, MPED::SEX_ARB);
  }

  if(sex == multipedigree.info().sex_code_male())
  {
    multipedigree.add_member(ped_name, ind_name, MPED::SEX_MALE);
  }
  else if(sex == multipedigree.info().sex_code_female())
  {
    multipedigree.add_member(ped_name, ind_name, MPED::SEX_FEMALE);
  }
  else
  {
    if(sex != multipedigree.info().sex_code_unknown())
    {
      errors << priority(warning) << "[" << line
             << "] Unable to read sex of individual "  << ind_name
             << " in pedigree " << ped_name
             << ".  Found '" << sex
             << "'.  Setting sex to unknown." << std::endl;
    }
    
    multipedigree.add_member(ped_name, ind_name, MPED::SEX_MISSING);
  }

  if(parent1 != multipedigree.info().individual_missing_code())
  {
    if(parent2 != multipedigree.info().individual_missing_code())
    {
      multipedigree.add_lineage(ped_name, ind_name, parent1, parent2);
    }
    else
    {
      multipedigree.add_lineage(ped_name, ind_name, parent1);
    }
  }
  else if(parent2 != multipedigree.info().individual_missing_code())
  {
    multipedigree.add_lineage(ped_name, ind_name, parent2);
  }

  // Add person to list of individuals:
  my_inds.insert(make_pair(ped_name, ind_name));

  // Check to see if there's a duplicate record:
  if(my_inds.count(make_pair(ped_name, ind_name)) > 1)
  {
    errors << priority(warning) 
           << "[" 
           << line
           << "] Duplicate record of individual "  
           << ind_name
           << " in pedigree " 
           << ped_name
           << "." 
           << std::endl;
  }
  else if(count <= verbose_output()) // Ok, there's no duplicate record, so stick 'em in the initial list:
  {
    my_ind_list.push_back(make_pair(ped_name, ind_name));
  }
}

//==================================================================   
//
//  normalize_allele_freq(...)
//
//==================================================================

// Added for allele frequency adjustment. - yjs Nov. 2002
//
void
RefPedigreeFile::normalize_allele_freq(RefMarkerInfo& marker)
{
  MLOCUS::genotype_model& gm = marker.gmodel();

  double f = 1.0 / gm.allele_count();

  MLOCUS::allele_iterator i = gm.allele_begin();

  for( ; i != gm.allele_end(); ++i )
  {
    gm.modify_allele_frequency(*i, f);
  }

  gm.normalize();
}

//==================================================================   
//
//  complement_allele_freq(...)
//
//==================================================================
void
RefPedigreeFile::complement_allele_freq(RefMarkerInfo& marker)
{
  MLOCUS::genotype_model& gm = marker.gmodel();

  MLOCUS::allele_iterator i = gm.allele_begin();

  for( ; i != gm.allele_end(); ++i )
  {
    double f = 1.0 - i->frequency();
    gm.modify_allele_frequency(*i, f);
  }

  gm.normalize();
}

//==================================================================   
//
//  adjust_min_allele_freq(...)
//
//==================================================================
void
RefPedigreeFile::adjust_min_allele_freq(RefMarkerInfo& marker, double f)
{
  MLOCUS::genotype_model& gm = marker.gmodel();

  MLOCUS::allele_iterator i = gm.allele_begin();

  for( ; i != gm.allele_end(); ++i )
  {
    if( i->frequency() < f )
    {
      gm.modify_allele_frequency(*i, f);
    }
  }

  gm.normalize();
}

//==================================================================   
//
//  adjust_max_allele_freq(...)
//
//==================================================================
void
RefPedigreeFile::adjust_max_allele_freq(RefMarkerInfo& marker, double f)
{
  MLOCUS::genotype_model& gm = marker.gmodel();

  MLOCUS::allele_iterator i = gm.allele_begin();

  for( ; i != gm.allele_end(); ++i )
  {
    if( i->frequency() > f )
    {
      gm.modify_allele_frequency(*i, f);
    }
  }

  gm.normalize();
}

//==================================================================   
//
//  parse_allele_frequency_adjustment(...)
//
//==================================================================
void
RefPedigreeFile::parse_allele_frequency_adjustment(RefMarkerInfo& marker, const LSFBase *param)
{
  if(!param)
    return;

  AttrList::const_iterator a, a1;
  AttrVal v;

  if(    (a=param->attrs()->find("equal_allele_freq")) != param->attrs()->end()
      || (a=param->attrs()->find("equal"))             != param->attrs()->end() )
    normalize_allele_freq(marker);
  else if(    (a=param->attrs()->find("complement_allele_freq")) != param->attrs()->end()
           || (a=param->attrs()->find("compl_allele_freq"))      != param->attrs()->end()
           || (a=param->attrs()->find("complement"))             != param->attrs()->end() )
    complement_allele_freq(marker);
  else
  {
    if(    (    (a =param->attrs()->find("minimum_allele_freq")) != param->attrs()->end()
             || (a =param->attrs()->find("minimum"))             != param->attrs()->end())
        && (    (a1=param->attrs()->find("maximum_allele_freq")) != param->attrs()->end()
             || (a1=param->attrs()->find("maximum"))             != param->attrs()->end()) )
    {
      if(    a ->second.Real() < 0. || a ->second.Real() > 1.
          && a1->second.Real() < 0. || a1->second.Real() > 1. )
        errors << priority(error) << "Invalid value for parameter 'm*_allele_freq'" << std::endl;
      else
      {
        adjust_min_allele_freq(marker, a ->second.Real());
        adjust_max_allele_freq(marker, a1->second.Real());
      }
    }
    else if(    (a =param->attrs()->find("minimum_allele_freq")) != param->attrs()->end()
             || (a =param->attrs()->find("minimum"))             != param->attrs()->end() )
    {
      if( a->second.Real() < 0. || a->second.Real() > 1. )
        errors << priority(error) << "Invalid value for parameter 'minimum_allele_freq'" << std::endl;
      else
      {
        adjust_min_allele_freq(marker, a->second.Real());
      }
    }
    else if(    (a =param->attrs()->find("maximum_allele_freq")) != param->attrs()->end()
             || (a =param->attrs()->find("maximum"))             != param->attrs()->end() )
    {
      if( a->second.Real() < 0. || a->second.Real() > 1. )
        errors << priority(error) << "Invalid value for parameter 'maximum_allele_freq'" << std::endl;
      else
      {
        adjust_max_allele_freq(marker, a->second.Real());
      }
    }
  }
}

//==================================================================   
//
//  update_marker_delimeter_info(...)
//
//==================================================================
void 
RefPedigreeFile::update_marker_delimiter_info(RefMarkerInfo& marker, char sep0)
{
  if( sep0 != marker.gmodel().separators()[0] )
  {
    marker.gmodel().set_unphased_separator(sep0);

    RefMarkerInfo::phenotype_iterator phi = marker.phenotype_begin();
    for( ; phi != marker.phenotype_end(); ++phi )
    {
      uint     id = phi->id();
      std::string   name = phi->name();
      std::string   allele1, allele2;
      MLOCUS::Ordering order;

      marker.gmodel().parse_genotype_name(name, allele1, allele2, order);

      std::string sep = marker.gmodel().separators();

      marker.alias_phenotype(id, allele1 + sep[0] + allele2);
      marker.alias_phenotype(id, allele1 + sep[1] + allele2);
      marker.alias_phenotype(id, allele1 + sep[2] + allele2);
    }
  }
}

//==================================================================   
//
//  update_marker_missing_info(...)
//
//==================================================================
void 
RefPedigreeFile::update_marker_missing_info(RefMarkerInfo& marker, const std::string& missing)
{
  parse_dynamic_markers_missing(marker, missing);
}

//==================================================================   
//
//  update_sex_linked_marker_info(...)
//
//==================================================================
void
RefPedigreeFile::update_sex_linked_marker_info(RefMultiPedigree &mp)
{
}

//
//--------------------------------------------------------------------------------------
//

bool
RefPedigreeFile::do_no_sex_structural_test(const RefMultiPedigree& mp)
{
  if(get_no_sex_ok_option()) return true;
  
  for(PedigreeConstIterator ped = mp.pedigree_begin(); ped != mp.pedigree_end(); ++ped)
  {
    // If the pedigree has no families, there is no structure, so we  skip it
    if(!ped->family_count()) continue;

    // Iterate through the subpedigrees
    for(SubpedigreeConstIterator sped = ped->subpedigree_begin();
        sped != ped->subpedigree_end(); ++sped)
    {
      // Iterate through the members, looking for members without sex information.
      for(MemberConstIterator mem = sped->member_begin(); mem != sped->member_end();
          ++mem)
      {
        // If we encounter an individual without sex, print the warning and exit.
        if(mem->is_sex_unknown())
        {
          errors << priority(critical)
                 << "Many S.A.G.E. algorithms rely on sex for pedigree structure "
                 << "information. Even when sex does not affect the outcome of an "
                 << "analysis, the lack of sex information may cause the program to "
                 << "behave improperly or even crash. Individuals without sex "
                 << "information, either direct or inferred, but with pedigree "
                 << "structure information have been detected in your data set. If "
                 << "your analyses involve pedigree structure, we recommend that "
                 << "this be corrected before continuing your analyses. If you would "
                 << "like to continue the analysis without correcting this issue, "
                 << "you must place a \"no_sex_ok\" parameter into your pedigree "
                 << "block." << endl;
                 
          return false;
        }
      }
    }
  }
  
  return true;
}

//
//--------------------------------------------------------------------------------------
//

//==================================================================   
//
//  CONSTRUCTOR
//
//==================================================================
RefLSFPedigreeFile::RefLSFPedigreeFile(cerrorstream &errors) : RefPedigreeFile(errors)
{}

//==================================================================   
//
//  process_parameters(...)
//
//==================================================================
bool
RefLSFPedigreeFile::process_parameters(RefMPedInfo &mped_info, const LSFBase* params)
{
  if(!params)
  {
    errors << priority(critical) << "RefPedigreeFile Error: No parameters specified." << std::endl;

    invalidate();
    
    return false;
  }
  else // params is valid
  {
    for(LSFList::const_iterator i=params->List()->begin(); i != params->List()->end(); ++i)
    {
      if(!process_parameter(mped_info, *i)) return false;
    }

    // Create sex_code covariate always as long as sex_field exist.
    //
    if( sex_field_name() != "" )
    {
      set_sex_code_trait(true);

      size_t t = mped_info.trait_find("SEX_CODE");

      // If not previously available, add it.
      if( t >= mped_info.trait_count() )
      {
        t = mped_info.add_binary_trait("SEX_CODE", RefTraitInfo::trait_covariate);

        mped_info.trait_info(t).set_string_missing_code     ("");
        mped_info.trait_info(t).set_numeric_missing_code    (-1);
        mped_info.trait_info(t).set_string_affected_code    (mped_info.sex_code_female());
        mped_info.trait_info(t).set_string_unaffected_code  (mped_info.sex_code_male());
        mped_info.trait_info(t).set_numeric_affected_code   (1);
        mped_info.trait_info(t).set_numeric_unaffected_code (0);
        mped_info.trait_info(t).set_alias_name(sex_field_name());
      }
    }

    // Create indicator covariates always as long as pedigree_id exist.
    //
    if( pedigree_id_name() != "" )
    {
      set_pedigree_id_trait(true);

      size_t t = mped_info.trait_find("FAMILIAL_INDICATOR");

      // If not previously available, add it.
      if( t >= mped_info.trait_count() )
      {
        t = mped_info.add_binary_trait("FAMILIAL_INDICATOR", RefTraitInfo::trait_covariate);

        mped_info.trait_info(t).set_string_missing_code     ("");
        mped_info.trait_info(t).set_numeric_missing_code    (-1);
        mped_info.trait_info(t).set_string_affected_code    ("1");
        mped_info.trait_info(t).set_string_unaffected_code  ("0");
        mped_info.trait_info(t).set_numeric_affected_code   (1);
        mped_info.trait_info(t).set_numeric_unaffected_code (0);
        mped_info.trait_info(t).set_alias_name("is_familiar");
      }

      t = mped_info.trait_find("FOUNDER_INDICATOR");

      // If not previously available, add it.
      if( t >= mped_info.trait_count() )
      {
        t = mped_info.add_binary_trait("FOUNDER_INDICATOR", RefTraitInfo::trait_covariate);

        mped_info.trait_info(t).set_string_missing_code     ("");
        mped_info.trait_info(t).set_numeric_missing_code    (-1);
        mped_info.trait_info(t).set_string_affected_code    ("1");
        mped_info.trait_info(t).set_string_unaffected_code  ("0");
        mped_info.trait_info(t).set_numeric_affected_code   (1);
        mped_info.trait_info(t).set_numeric_unaffected_code (0);
        mped_info.trait_info(t).set_alias_name("is_founder");
      }

      t = mped_info.trait_find("PEDIGREE_SIZE");

      // If not previously available, add it.
      if( t >= mped_info.trait_count() )
      {
        t = mped_info.add_continuous_trait("PEDIGREE_SIZE", RefTraitInfo::trait_covariate);

        mped_info.trait_info(t).set_string_missing_code     ("");
        mped_info.trait_info(t).set_numeric_missing_code    (std::numeric_limits<double>::quiet_NaN());
        mped_info.trait_info(t).set_alias_name("ped_size");
      }
    }
  }

  return true;
}

//==================================================================   
//
//  process_parameter(...)
//
//==================================================================
bool
RefLSFPedigreeFile::process_parameter(RefMPedInfo &mped_info, const LSFBase* param)
{
  if(!param)
    return true;

  // Set up local variables:
  
  std::string              name = toUpper(param->name());

       if(name == "FORMAT")                   return process_format                   (mped_info, param);
  else if(name == "VERBOSE")                  return process_verbose                  (mped_info, param);
  else if(name == "REQUIRE_RECORD")           return process_require_record           (mped_info, param);
  else if(name == "SKIP_TRAITS")              return process_skip_traits              (mped_info, param);
  else if(name == "SKIP_MARKERS")             return process_skip_markers             (mped_info, param);
  else if(name == "DYNAMIC_MARKERS")          return process_dynamic_markers          (mped_info, param);
  else if(name == "REJECT_PARTIAL_LINEAGE")   return process_reject_partial_lineage   (mped_info, param);
  else if(name == "INDIVIDUAL_MISSING_VALUE") return process_individual_missing_value (mped_info, param);
  else if(name == "SEX_CODE")                 return process_sex_code                 (mped_info, param);
  else if(name == "NO_SEX_OK")                return process_no_sex_ok                (mped_info, param);

  return true;
}

/// Processes the "FORMAT" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool
RefLSFPedigreeFile::process_format(RefMPedInfo &mped_info, const LSFBase* param)
{
  AttrVal val = attr_value(param, 0);

  if(val.has_value())
  {
    set_format(val.String());
  }

  return true;
}

/// Processes the "VERBOSE" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool
RefLSFPedigreeFile::process_verbose(RefMPedInfo &mped_info, const LSFBase* param)
{
  AttrVal val = attr_value(param, 0);

  if( val.has_value() && val.Int() > 0 )
  {
    set_verbose_output(val.Int());
  }

  return true;
}

/// Processes the "REQUIRE_RECORD" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool
RefLSFPedigreeFile::process_require_record(RefMPedInfo &mped_info, const LSFBase* param)
{
  AttrVal val = attr_value(param, 0);

  if(val.has_value())
  {
    if( toUpper(val.String()) == "TRUE")
    {
      set_require_record(true);
    }
    else if( toUpper(val.String()) == "FALSE")
    {
      set_require_record(false);
    }
    else
    {
      errors << priority(error) << "Unknown value for parameter 'require_record'" << std::endl;
    }
  }

  return true;
}

/// Processes the "SKIP_TRAITS" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool
RefLSFPedigreeFile::process_skip_traits(RefMPedInfo &mped_info, const LSFBase* param)
{
  AttrVal val = attr_value(param, 0);

  if(val.has_value() && !force_skip_traits())
  {
    if( toUpper(val.String()) == "TRUE")
    {
      set_skip_traits(true);
    }
    else if( toUpper(val.String()) == "FALSE")
    {
      set_skip_traits(false);
    }
    else
    {
      errors << priority(error) << "Unknown value for parameter 'skip_traits'" << std::endl;
    }
  }

  return true;
}

/// Processes the "SKIP_MARKERS" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool
RefLSFPedigreeFile::process_skip_markers(RefMPedInfo &mped_info, const LSFBase* param)
{
  AttrVal val = attr_value(param, 0);

  if(val.has_value() && !force_skip_markers())
  {
    if( toUpper(val.String()) == "TRUE")
    {
      set_skip_markers(true);
    }
    else if( toUpper(val.String()) == "FALSE")
    {
      set_skip_markers(false);
    }
    else
    {
      errors << priority(error) << "Unknown value for parameter 'skip_markers'" << std::endl;
    }
  }

  return true;
}

/// Processes the "DYNAMIC_MARKERS" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool
RefLSFPedigreeFile::process_dynamic_markers(RefMPedInfo &mped_info, const LSFBase* param)
{
  AttrVal val = attr_value(param, 0);

  if(val.has_value() && !force_dynamic_markers())
  {
    if( toUpper(val.String()) == "TRUE")
    {
      set_dynamic_markers(true);
    }
    else if( toUpper(val.String()) == "FALSE")
    {
      set_dynamic_markers(false);
    }
    else
    {
      errors << priority(error) << "Unknown value for parameter 'dynamic_markers'" << std::endl;
    }
  }

  return true;
}

/// Processes the "REJECT_PARTIAL_LINEAGE" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool
RefLSFPedigreeFile::process_reject_partial_lineage(RefMPedInfo &mped_info, const LSFBase* param)
{
  AttrVal val = attr_value(param, 0);

  if( val.has_value() )
  {
    if( toUpper(val.String()) == "TRUE")
    {
      set_reject_partial_lineage(true);
    }
    else if( toUpper(val.String()) == "FALSE")
    {
      set_reject_partial_lineage(false);
    }
    else
    {
      errors << priority(error) << "Unknown value for parameter 'reject_partial_lineage'" << std::endl;
    }
  }

  return true;
}

/// Processes the "INDIVIDUAL_MISSING_VALUE" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool
RefLSFPedigreeFile::process_individual_missing_value(RefMPedInfo &mped_info, const LSFBase* param)
{
  AttrVal val = attr_value(param, 0);

  if( val.has_value() )
  {
    mped_info.set_individual_missing_code(val.String());

    DEBUG_RPEDFILE(errors << priority(debug) << "Found individual missing value = " << mped_info.individual_missing_code() << std::endl;)
  }

  return true;
}

/// Processes the "SEX_CODE" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool
RefLSFPedigreeFile::process_sex_code(RefMPedInfo &mped_info, const LSFBase* param)
{
  AttrList::const_iterator a;
  
  if(param->attrs())
  {
    if( (a=param->attrs()->find("missing")) != param->attrs()->end() )
    {
      mped_info.set_sex_code_unknown(strip_ws(a->second.String()));

      DEBUG_RPEDFILE(errors << priority(debug) << "Found sex missing value code = " << mped_info.sex_code_unknown() << std::endl;)
    }
    if( (a=param->attrs()->find("male"))    != param->attrs()->end() )
    {
      mped_info.set_sex_code_male(strip_ws(a->second.String()));

      DEBUG_RPEDFILE(errors << priority(debug) << "Found sex male code = " << mped_info.sex_code_male() << std::endl;)
    }
    if( (a=param->attrs()->find("female"))  != param->attrs()->end() )
    {
      mped_info.set_sex_code_female(strip_ws(a->second.String()));

      DEBUG_RPEDFILE(errors << priority(debug) << "Found sex female code = " << mped_info.sex_code_female() << std::endl;)
    }
  }

  return true;
}

/// Processes the "NO_SEX_OK" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool
RefLSFPedigreeFile::process_no_sex_ok(RefMPedInfo &mped_info, const LSFBase* param)
{
  AttrVal val = attr_value(param, 0);

  if( val.has_value() )
  {
    if( toUpper(val.String()) == "TRUE")
    {
      set_no_sex_ok_option(true);
    }
    else if( toUpper(val.String()) == "FALSE")
    {
      set_no_sex_ok_option(false);
    }
    else
    {
      errors << priority(error) << "Unknown value for parameter 'no_sex_ok'" << std::endl;
    }
  }

  return true;
}

} // End namespace RPED
} // End namespace SAGE
