#include <string>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <iomanip>
#include "mped/sp.h"
#include "rped/rpfile.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "error/internal_error.h"
#include "util/StringUtils.h"

#define DEBUG_RPEDFILE(x)

namespace SAGE {
namespace RPED {

//=======================================================================
//=======================================================================
//=======================================================================
//
//  RefDelimitedPedigreeFile
//
//=======================================================================
//=======================================================================
//=======================================================================

//=======================================================================
//
//  Constructor
//
//=======================================================================
RefDelimitedPedigreeFile::RefDelimitedPedigreeFile(cerrorstream &err) : RefPedigreeFile(err)
{
  set_format_in_file              (true);
  set_whitespace                  (" \r\n");
  set_delimiters                  ("\t,");
  set_skip_consecutive_delimiters (false);
  set_skip_leading_delimiters     (false);
  set_skip_trailing_delimiters    (false);
  reset_counts                    ();
}

//=======================================================================
//
//  Destructor
//
//=======================================================================
RefDelimitedPedigreeFile::~RefDelimitedPedigreeFile()
{}

//=======================================================================
//
//  set_delimiters(...)
//
//=======================================================================
void 
RefDelimitedPedigreeFile::set_delimiters(const std::string &d)
{
  my_delimiters = d;

  if(my_delimiters.find("\n") == std::string::npos)
    my_delimiters += "\n";

  if(my_delimiters.find("\r") == std::string::npos)
    my_delimiters += "\r";
}

//=======================================================================
//
//  set_whitespace(...)
//
//=======================================================================
void 
RefDelimitedPedigreeFile::set_whitespace(const std::string &w)
{
  my_whitespace = w;

  if(my_whitespace.find("\n") == std::string::npos)
    my_whitespace += "\n";

  if(my_whitespace.find("\r") == std::string::npos)
    my_whitespace += "\r";
}


//=======================================================================
//
//  build_fields(...)
//
//=======================================================================
bool 
RefDelimitedPedigreeFile::build_fields(string_tokenizer&  header,
                                       const RefMPedInfo& mped_info,
                                       bool               quiet)
{
  my_fields.clear();

  for(string_tokenizer::iterator i = header.begin(); i != header.end(); ++i)
  {
    my_fields.push_back(field());

    if(i->empty())
       continue;

    field_map_type::const_iterator f = my_field_map.find( toUpper(*i) );
    
    if(f != my_field_map.end())
    {
      // If the field is in our map
      build_field(mped_info, f, quiet);
    }
    else
    {
      DEBUG_RPEDFILE(errors << priority(information) << "Field '" << *i << "' in pedigree file is skipped." << std::endl;)
    }
  }

  return true;
}

void 
RefDelimitedPedigreeFile::build_field(const RefMPedInfo&             mped_info,
                                      field_map_type::const_iterator f,
                                      bool                           quiet)
{
  // Verify that the field should not be skipped
  if(f->second.type == skip)
  {
    DEBUG_RPEDFILE(errors << priority(information) << "Field '" << *i << "' in pedigree file is skipped." << std::endl;)
    return;
  }

  field & current_field = my_fields.back();

  // This makes sense -- current field is a reference
  current_field = f->second;

  if( current_field.type == trait )
  {
    build_trait_field(mped_info, current_field);
  }
  else if( current_field.type == string_field )
  {
    build_string_field(mped_info, current_field);
  }
  else if( current_field.type == marker || current_field.type == allele )
  {
    build_marker_field(mped_info, current_field);
  }
  else if( current_field.type == marker_cov || current_field.type == allele_cov )
  {
    build_trait_field(mped_info, current_field);
  }

  DEBUG_RPEDFILE(errors << priority(information) << "Using field '" << *i << "'." << std::endl;)
}

void
RefDelimitedPedigreeFile::build_trait_field(const RefMPedInfo& mped_info,
                                            field&             current_field)
{
  size_t t = mped_info.trait_find(current_field.name);

  if(t < mped_info.trait_count())
  {
    current_field.index = t;
  }
  else
  {
    current_field.index = (size_t)-1;

    // Can never happen (yet) since we always add the trait to mped_info
    // in the LSF parser.  But we have to look out for the future.
    errors << priority(warning) << "Invalid trait '" << current_field.name << "'.  Skipping..." << std::endl;
  }
}

void
RefDelimitedPedigreeFile::build_marker_field(const RefMPedInfo& mped_info,
                                             field&             current_field)
{
  size_t m = mped_info.marker_find(current_field.name);

  if(m < mped_info.marker_count())
  {
    current_field.index = m;
  }
  else
  {
    current_field.index = (size_t)-1;

    if(mped_info.marker_count())
    {
      errors << priority(warning) << "Invalid marker '" << current_field.name << "'.  Skipping..." << std::endl;
    }
  }
}
void
RefDelimitedPedigreeFile::build_string_field(const RefMPedInfo& mped_info,
                                             field&             current_field)
{
  size_t s = mped_info.string_find(current_field.name);

  if(s < mped_info.string_count())
  {
    current_field.index = s;
  }
  else
  {
    current_field.index = (size_t)-1;

    // Can never happen (yet) since we always add the trait to mped_info
    // in the LSF parser.  But we have to look out for the future.
    errors << priority(warning) << "Invalid string field '" << current_field.name
           << "'.  Skipping..." << std::endl;
  }
}

//=================================================
//
//  check_header_vs_field_map(...)
//
//=================================================
bool 
RefDelimitedPedigreeFile::check_header_vs_field_map(string_tokenizer & header)
{
  bool match = true;

  for(field_map_type::const_iterator f = my_field_map.begin(); f != my_field_map.end(); ++f)
  {
    const field& current_field = f->second;

    if(current_field.type == trait)
    {
      bool exist_in_header = false;

      for(string_tokenizer::iterator h = header.begin(); h != header.end(); ++h)
      {
        if(toUpper(f->first) == toUpper(*h))
        {
          exist_in_header = true;

          break;
        }
      }

      if(!exist_in_header)
      {
        match = false;

        errors << priority(warning) << "Trait '" << current_field.name << "' doesn't exist in the pedigree file." << std::endl;
      }
    }
    else if(current_field.type == string_field)
    {
      bool exist_in_header = false;

      for(string_tokenizer::iterator h = header.begin() ; h != header.end(); ++h )
      {
        if(toUpper(f->first) == toUpper(*h))
        {
          exist_in_header = true;
          break;
        }
      }

      if( !exist_in_header )
      {
        match = false;

        errors << priority(warning) << "String '" << current_field.name << "' doesn't exist in the pedigree file." << std::endl;
      }
    }
  }

  return match;
}

/// Validates the format and file, then reads in the data.
///
bool
RefDelimitedPedigreeFile::input(RefMultiPedigree&  p,
                                const std::string& filename,
                                ostream&           messages)
{
  if(!validate_file(filename) || !validate_format(filename))
    return false;

  return RefPedigreeFile::input(p, filename, messages);
}

/// Verifies that the filename given exists and is valid.
///
bool
RefDelimitedPedigreeFile::validate_file(const std::string &filename)
{
  if(!filename.size())
  {
    errors << priority(fatal) << "No Family Data file specified." << std::endl;

    return false;
  }

  std::ifstream infile(filename.c_str());

  if(!infile.good())
  {
    errors << priority(fatal) << "Unable to open Family Data file '" << filename << "'. Please check your file." << std::endl;

    return false;
  }

  return true;
}

/// Checks to see if the format is valid, or, if it's in the file, reads it
/// from the file and checks that.  In either case, the format is then available
/// for use.
bool
RefDelimitedPedigreeFile::validate_format(const std::string &filename)
{
  // If the format should be in the file, we have to go and get it before we
  // can validate it.
  if(format_in_file())
  {
    // We assume that validate_file() has already been done
    std::ifstream infile(filename.c_str());

    std::string line;

    getline(infile, line);
    set_format(line);
  }
  
  // Check the format
  if(!format().size())
  {
    errors << priority(fatal) << "No format specified to read headerless Family Data file." << std::endl;

    return false;
  }

  return true;
}

//====================================================
//
//  input_pedigree(...)
//
//=====================================================
bool 
RefDelimitedPedigreeFile::input_pedigree(RefMultiPedigree& p,
                                         const string&     filename,
                                         ostream&          messages,
                                         bool              quiet)
{
  std::ifstream infile(filename.c_str());

  string_tokenizer tokenizer( format() );
  
  setup_tokenizer(tokenizer);

  std::string line;   // Current line

  // Skip the header if it's in the file.  validate_format() will have already
  // set it.
  if( format_in_file() )
  {
    getline(infile, line);
  }

  const RefMPedInfo &mped_info = p.info();

  if( !build_fields(tokenizer, mped_info, quiet) )
  {
    errors << priority(error) << "Cannot build list of fields to read.  Aborting." << std::endl;

    invalidate();
    return false;
  }

  if( !validate_fields(false, quiet) )
  {
    errors << priority(error) << "Invalid list of fields to read.  Aborting." << std::endl;

    invalidate();
    return false;
  }

  // Figure out if there's supposed to be a pedigree id field:
//  bool has_pedigree_id_field = false;
  
//  for(field_list_type::const_iterator itr = my_fields.begin(); itr != my_fields.end(); ++itr)
//  {
//    if(itr->type == pedigree_id)
//      has_pedigree_id_field = true;
//  }

  for( size_t count = 1; infile.good(); ++count )
  {
    getline(infile, line);

    if(!line.size())
      continue;

    std::string ped_name = "",                           // Pedigree id
                ind_name = "",                           // Individual id
                parent1  = "",                           // Parent1
                parent2  = "",                           // Parent2
                sex      = mped_info.sex_code_unknown(); // Sex code

    tokenizer.set_str(line);

    field_list_type::const_iterator  field_info     = my_fields.begin (),
                                     field_info_end = my_fields.end   ();
    string_tokenizer::const_iterator field          = tokenizer.begin (),
                                     field_end      = tokenizer.end   ();

    // Loop across tokens and field_info's at the same time, processing each:
    for( ; field_info != field_info_end && field != field_end; ++field_info, ++field)
    {
      switch( field_info->type )
      {
        case pedigree_id   : ped_name = *field;    break;
        case individual_id : ind_name = *field;    break;
        case parent_id     : if(!parent1.size())
                               parent1 = *field;
                             else
                               parent2 = *field;  break;
        case sex_code      : sex = *field;        break;
        default:                                  break;
      }

    } // End loop-across-tokens
    
    // Check for treat_ped_id options:
    if( !pedigree_id_count() )
      ped_name = "0";
    else if( get_treat_as_sibs() == true ) // There's a pedigree_id field and treat_as_sibs is enabled:
    {
      parent1 = ped_name + "_parent1";
      parent2 = ped_name + "_parent2";
      
      p.add_member(ped_name, parent1, MPED::SEX_MALE);
      p.add_member(ped_name, parent2, MPED::SEX_FEMALE);
    }

    // If parents are still empty for some reason, set them to be missing:
    parent1 = parent1.size() ? parent1 : mped_info.individual_missing_code();
    parent2 = parent2.size() ? parent2 : mped_info.individual_missing_code();
    
    // Skip this person if any essential bits of info are missing:
    if(!ped_name.size() && !ind_name.size() && !parent1.size() && !parent2.size() && !sex.size())
      continue;

    DEBUG_RPEDFILE(errors << priority(debug) << "Found (" << pn  << "," << id  << "," << sex << "," << p1 << "," << p2  << ")" << std::endl; )

    add_member(p, ped_name, ind_name, sex, parent1, parent2, count + format_in_file(), count);
  }

  infile.close();

  if( !build_pedigree(p) )
  {
    errors << priority(error)
           << "Fatal error building pedigree data structure.  Aborting." << std::endl;
    return false;    
  }

  return true;
}

bool 
RefDelimitedPedigreeFile::input_data(RefMultiPedigree& p,
                                     const string&     filename,
                                     ostream&          messages,
                                     bool              quiet)
{
  // Check to see if a second pass is necessary
  // (traits and covariates were defered in the first pass)
  //
  RefMPedInfo &mped_info = p.info();

  if( trait_count()            == 0 &&
      !sex_code_trait()             &&
      !pedigree_id_trait()          &&
      marker_count()           == 0 && 
      marker_cov_count()       == 0 &&
      string_count()           == 0 &&
      mped_info.trait_count()  == 0 &&
      mped_info.marker_count() == 0 &&
      mped_info.string_count() == 0 )
    return true;

  // We assume validate_file() has been called already
  std::ifstream infile( filename.c_str());

  string_tokenizer tokenizer( format() );
  
  setup_tokenizer(tokenizer);

  string line;   // Current line

  // Skip header if it is in the file.  We assume that it's already been 
  // set by validate_format()
  if( format_in_file() )
  {
    getline(infile, line);
  }

  if( !build_fields(tokenizer, mped_info, quiet) )
  {
    errors << priority(error)
           << "Cannot build list of fields to read.  Aborting." << std::endl;
    invalidate();
    return false;
  }

  if( !validate_fields(true, quiet) )
  {
    errors << priority(error)
           << "Invalid list of fields to read.  Aborting." << std::endl;
    invalidate();
    return false;
  }

  if( !check_header_vs_field_map(tokenizer) )
  {
    errors << priority(error)
           << "Pedigree block specifies misspelled or missing pedigree data fields."
           << " Please see inf file for details."<< std::endl;
    invalidate();
    return false;
  }

  if( !build_data(p) )
  {
    errors << priority(error)
           << "Error building data indices.  Aborting." << std::endl;
    invalidate();
    return false;
  }

  size_t verbose = verbose_output();

  verbose = verbose_output();

  field_list_type::const_iterator field_info      = my_fields.begin();
  field_list_type::const_iterator field_info_end  = my_fields.end();

  vector<pair<size_t,string> >  trait_values( trait_count() );
  vector<pair<size_t,string> > string_values( string_count() );

  // FIXME: This approach to reading markers will certainly have serious
  //        performance problems in the large scale and is here only for the
  //        simplicity.  A hash_map implementation would be much better.  A
  //        separate index vector is the optimal solution in the long term.
  typedef pair<string, string>          string_pair;
  typedef std::map<size_t, string_pair> marker_value_map;

  for( size_t count=1; infile.good(); ++count )
  {
    getline(infile, line);

    if(!line.size())
      continue;

    string pn;     // Pedigree id
    string id;     // Individual id
    string p1, p2; // Parents
    string sex;    // Sex code

    tokenizer.set_str(line);

    string_tokenizer::const_iterator field     = tokenizer.begin();
    string_tokenizer::const_iterator field_end = tokenizer.end();

    marker_value_map marker_values;
    marker_value_map marker_covariate_values;

    for( size_t t = 0; t < trait_count(); ++t )
      trait_values[t].first = (size_t)-1;

    for( size_t s = 0; s < string_count(); ++s )
      string_values[s].first = (size_t)-1;

    size_t tfound = 0;
    size_t sfound = 0;

    for( field_info = my_fields.begin() ;
         field_info != field_info_end && field != field_end ;
         ++field_info, ++field)
    {
      switch( field_info->type )
      {
        case   pedigree_id: pn = *field;    break;
        case individual_id: id = *field;    break;
        case         trait:
                            if( field_info->index >= mped_info.trait_count() )
                              break;
                            trait_values[tfound].first  = field_info->index;
                            trait_values[tfound].second = *field;
                            ++tfound;
                            break;

        case  string_field:
                            if( field_info->index >= mped_info.string_count() )
                              break;
                            string_values[sfound].first  = field_info->index;
                            string_values[sfound].second = strip_ws(*field);
                            ++sfound;
                            break;

        case        allele:
        case        marker:
        {
                            size_t index = field_info->index;
                            if( index >= mped_info.marker_count() )
                              break;

                            string_pair &values = marker_values[index];
                            if( !values.first.size())
                              values.first = strip_ws(*field);
                            else
                              values.second = strip_ws(*field);
                            break;
        }

        case        allele_cov:
        case        marker_cov:
        {
                            size_t index = field_info->index;

                            if( index >= mped_info.trait_count() )
                              break;

                            string_pair &values = marker_covariate_values[index];
                            if( !values.first.size())
                              values.first = strip_ws(*field);
                            else
                              values.second = strip_ws(*field);
                            break;
        }

        default:                            break;
      }
    }

    if( !pedigree_id_count() )
      pn = "0";

    if( !pn.size() || !id.size() )
      continue;

    if( !id.size() || id == mped_info.individual_missing_code() )
    {
      errors << priority(warning) << "[" << count + format_in_file()
             << "] Found an individual ID that is missing or"
             << " the same as the missing individual/parent code.  Skipping..."
             << std::endl;
      continue;
    }

    RefMultiPedigree::member_pointer mem = p.member_find(pn,id);

    // This should only happen when an error has already been reported
    if( !mem ) continue;

    RefMultiPedigree::pedinfo_type &info = mem->pedigree()->info();
    int ind_num = mem->index();

    // Make the traits
    for( size_t tt=0; tt < tfound; ++tt )
    {
      size_t t = trait_values[tt].first;
      const string &value = trait_values[tt].second;
      int code = info.set_trait(ind_num, t, value, mped_info.trait_info(t));

      switch(code)
      {
        case 0:                       // trait ok
        case 1:                       // trait ok, but missing
        case 4:                       // trait not set legitimately
                 break;
        case 3:                       // invalid ind. or trait id
                 errors << priority(warning) << "["<<count+format_in_file()
                        << "] Cannot set trait '"
                        << mped_info.trait_info(t).name()
                        << "' for individual '" << id << "' in pedigree '"
                        << pn << "'." << std::endl;
                 break;
        case 2:                       // bad trait value
                 errors << priority(warning) << "["<<count+format_in_file()
                        << "] Unrecognized value for trait '"
                        << mped_info.trait_info(t).name()
                        << "' for individual '" << id << "' in pedigree '"
                        << pn << "'.  Found '" << value << "'." << std::endl;

                 break;
       default:
                 errors << priority(error) << "[" <<count+format_in_file()
                        << "] Unexpected error setting trait '"
                        << mped_info.trait_info(t).name()
                        << "' for individual '" << id << "' in pedigree '"
                        << pn << "'." << std::endl;
                 break;
      }
    }

    for( size_t ss=0; ss < sfound; ++ss )
    {
      size_t s = string_values[ss].first;
      const string &value = string_values[ss].second;
      bool code = info.set_string(ind_num, s, value);

      if(!code)
        errors << priority(warning) << "["<<count+format_in_file() 
	       << "] Cannot set string field '"
               << mped_info.string_info(s).name()
               << "' for individual '" << id << "' in pedigree '"
               << pn << "'." << std::endl;
    }

    marker_value_map::const_iterator mm;
    for( mm = marker_values.begin(); mm != marker_values.end(); ++mm )
    {
      size_t m = mm->first;
      const string_pair &values = mm->second;

      int code;
      code = info.set_phenotype(ind_num, mem->get_effective_sex(), m, values.first, values.second,
                                mped_info.marker_info(m));

      switch(code)
      {
        case 0:                       // marker ok
        case 1:                       // marker ok, but missing
        case 3:                       // invalid ind. or marker id <-- Shouldn't ever happen
                 break;
        case 2:                       // bad marker value
                 errors << priority(warning) << "["<<count+format_in_file()
                        << "] Unrecognized value for marker '"
                        << mped_info.marker_info(m).name()
                        << "' of individual '" << id << "' in pedigree '"
                        << pn << "': Found '" << values.first << "'";
                 if(values.second.size())
                   errors << ", '" << values.second << "'";
                 errors << ". Marker will be set to missing for this individual." << std::endl;
                 break;
         case 4:
                 errors << priority(warning) << "[" << count+format_in_file()
                        << "] Marker '" << mped_info.marker_info(m).name()
                        << "' is sex-dependent, but the individual '" << id 
                        << "' in pedigree '" << pn 
                        << "' has unknown sex.  Marker will be set to missing for this individual." << std::endl;
                 break;
         case 5:
                 errors << priority(warning) << "[" << count+format_in_file()
                        << "] Phenotype '" << values.first << "'";
                 if(values.second.size())
                   errors << ", '" << values.second << "'";
                 errors << " at sex-dependent marker '" << mped_info.marker_info(m).name()
                        << "' is inconsistent with the sex of individual '" << id
                        << "' in pedigree '" << pn 
                        << "'.  Marker will be set to missing for this individual." << std::endl;
                 break;
      }
    }

    for( mm = marker_covariate_values.begin(); mm != marker_covariate_values.end(); ++mm )
    {
      size_t t = mm->first;
      const string_pair &values = mm->second;

      string mcov_name = mped_info.trait_info(t).name();

      const string &value = get_marker_covariate_value(mcov_name, values.first, values.second);

      int code = info.set_trait(ind_num, t, value, mped_info.trait_info(t));

      switch(code)
      {
        case 0:                       // trait ok
        case 1:                       // trait ok, but missing
        case 4:                       // trait not set legitimately
                 break;
        case 3:                       // invalid ind. or trait id
                 errors << priority(warning) << "["<<count+format_in_file()
                        << "] Cannot set marker covariate '"
                        << mped_info.trait_info(t).name()
                        << "' for individual '" << id << "' in pedigree '"
                        << pn << "'." << std::endl;
                 break;
        case 2:                       // bad trait value
                 errors << priority(warning) << "["<<count+format_in_file()
                        << "] Unrecognized value for marker covariate '"
                        << mped_info.trait_info(t).name()
                        << "' for individual '" << id << "' in pedigree '"
                        << pn << "'.  Found '" << value << "'." << std::endl;

                 break;
       default:
                 errors << priority(error) << "[" <<count+format_in_file()
                        << "] Unexpected error setting marker covariate '"
                        << mped_info.trait_info(t).name()
                        << "' for individual '" << id << "' in pedigree '"
                        << pn << "'." << std::endl;
                 break;
      }
    }
  }

  return true;
}

bool
RefDelimitedPedigreeFile::output(RefMultiPedigree& p,
                                 const string&     filename,
                                 ostream&          messages,
                                 bool              quiet)
{
  if( !filename.size() || !delimiters().size()) return false;

  const RefMPedInfo &mped_info = p.info();

  if( format().size() )
  {
    string_tokenizer tokenizer( format() );
    tokenizer.set_whitespace( whitespace() );
    tokenizer.set_delimiters( delimiters() );
    tokenizer.set_skip_consecutive_delimiters( skip_consecutive_delimiters() );
    tokenizer.set_skip_leading_delimiters( skip_leading_delimiters() );
    tokenizer.set_skip_trailing_delimiters( skip_trailing_delimiters() );

    if(!build_fields(tokenizer, mped_info, quiet))
    {
      errors << priority(error)
             << "Cannot build list of fields to write.  Aborting..." << std::endl;
      invalidate();
      return false;
    }

    if(!validate_fields(false, quiet))
    {
      errors << priority(error)
             << "Invalid list of fields to write.  Aborting." << std::endl;
      invalidate();
      return false;
    }
  }

  if( !my_fields.size() )
  {
    errors << priority(error)
           << "No fields to write.  Aborting..." << std::endl;
    invalidate();
    return false;
  }

  std::ofstream outfile(filename.c_str());

  if(!outfile.good())
  {
    errors << priority(fatal) << "Unable to open output Family Data file '" << filename
           << "'." << std::endl;
    return false;
  }

  string delimiter = delimiters().substr(0,1);

  field_list_type::const_iterator field_info;
  field_list_type::const_iterator field_info_end  = my_fields.end();

  bool first_field = true;
  for( field_info = my_fields.begin();
       field_info != field_info_end ;
     ++field_info )
  {
    if( !first_field )
      outfile << delimiter;
    first_field = false;
    outfile << field_info->field_name;
  }
  outfile << std::endl;

  RefPedigree::member_type* mid;
  RefPedigree::member_type* pid;
  size_t parent_count;

  vector<size_t> markers_written( mped_info.marker_count() );

  RefMultiPedigree::pedigree_iterator ped;
  for(ped = p.pedigree_begin(); ped != p.pedigree_end(); ++ped)
  {
    const RefPedInfo &ped_info = ped->info();

    for(unsigned int i = 0; i < ped->member_count(); ++i)
    {
      mid = &ped->member_index(i);
      parent_count = 0;
      markers_written.clear();

      first_field = true;

      for(field_info = my_fields.begin(); field_info != field_info_end; ++field_info)
      {
        if( !first_field )
          outfile << delimiter;
        first_field = false;

        switch( field_info->type )
        {
          case pedigree_id:   outfile << ped->name();
                              break;

          case individual_id: outfile << mid->name();
                              break;

          case parent_id:
          {
                              string value = mped_info.individual_missing_code();
                              pid = NULL;
                              if(parent_count == 0 && mid->parent1())
                              {
                                pid = mid->parent1();
                                parent_count = 1;
                              }
                              if(!pid && parent_count < 2 && mid->parent2())
                              {
                                pid = mid->parent2();
                                parent_count = 2;
                              }
                              if(pid)
                                value = pid->name();
                              else
                                parent_count = 2;
                              outfile << value;
                              break;
          }
          case      sex_code:
                              if( mid->is_male() )
                                outfile << mped_info.sex_code_male();
                              else if( mid->is_female() )
                                outfile << mped_info.sex_code_female();
                              else
                                outfile << mped_info.sex_code_unknown();
                              break;

          case         trait:
          {
                              if( field_info->index >= mped_info.trait_count() )
                                break;

                              if( ped_info.trait_missing(i, field_info->index) )
                              {
                                outfile << mped_info.trait_info(
                                      field_info->index).string_missing_code();
                                break;
                              }

                              outfile << ped_info.trait(i, field_info->index);
                              break;
          }

          case  string_field:
          {
                              if( field_info->index >= mped_info.string_count() )
                                break;
                              outfile << ped_info.get_string(i, field_info->index);
                              break;
          }

          case        allele:
          {
                              if( field_info->index >= mped_info.marker_count() )
                                break;

                              if(markers_written[field_info->index] < 2)
                              {
                                const RefMarkerInfo& minfo =
                                   mped_info.marker_info(field_info->index);

                                uint pheno_id = ped_info.phenotype(i, field_info->index);

                                string pheno_name = minfo.get_phenotype(pheno_id).name();

                                string al1, al2;
                                MLOCUS::Ordering order;

                                minfo.gmodel().parse_genotype_name(pheno_name, al1, al2, order);

                                if(!al1.size()) al1 = minfo.missing_allele_name();
                                if(!al2.size()) al2 = minfo.missing_allele_name();

                                if(markers_written[field_info->index] == 0)
                                  outfile << al1;
                                else
                                  outfile << al2;
                                markers_written[field_info->index]++;
                              }
                              break;
          }

          case        marker:
          {
                              if( field_info->index >= mped_info.marker_count() )
                                break;
                              if(markers_written[field_info->index] == 0)
                              {
                                const RefMarkerInfo& minfo =
                                   mped_info.marker_info(field_info->index);

                                uint pheno = ped_info.phenotype(i, field_info->index);

                                outfile << minfo.get_phenotype(pheno).name();
                                markers_written[field_info->index] = 2;
                              }
                              break;
          }

          default:
                              break;
        }
      }
      outfile << std::endl;
    }
  }
  outfile.close();
  return true;
}

void 
RefDelimitedPedigreeFile::setup_tokenizer(string_tokenizer& tokenizer)
{
  tokenizer.set_whitespace( whitespace() );
  tokenizer.set_delimiters( delimiters() );
  tokenizer.set_skip_consecutive_delimiters( skip_consecutive_delimiters() );
  tokenizer.set_skip_leading_delimiters( skip_leading_delimiters() );
  tokenizer.set_skip_trailing_delimiters( skip_trailing_delimiters() );
}

//============================================================
//============================================================
//============================================================
//
//  RefLSFDelimitedPedigreeFile
//
//============================================================
//============================================================
//============================================================

RefLSFDelimitedPedigreeFile::RefLSFDelimitedPedigreeFile(cerrorstream &errors)
{
  set_error_sink(errors);

  // Set defaults
  set_delimiters(",\t");
  set_whitespace(" ");
  set_skip_trailing_delimiters(false);
  set_skip_leading_delimiters(false);
  set_skip_consecutive_delimiters(false);
}

bool
RefLSFDelimitedPedigreeFile::process_parameters(RefMPedInfo &mped_info, const LSFBase *params)
{
  return RefLSFPedigreeFile::process_parameters(mped_info, params);
}
  
bool 
RefLSFDelimitedPedigreeFile::process_parameter(RefMPedInfo &mped_info, const LSFBase *param)
{
  if(!param)
    return true;

  // Allow generic super-class to handle parameters
  RefLSFPedigreeFile::process_parameter(mped_info, param);

  string name = toUpper( param->name() );

       if(name == "FORMAT")                      return process_format2                     (mped_info, param);
  else if(name == "WHITESPACE")                  return process_whitespace                  (mped_info, param);
  else if(name == "DELIMITER_MODE")              return process_delimiter_mode              (mped_info, param);
  else if(name == "DELIMITERS")                  return process_delimiters                  (mped_info, param);
  else if(name == "SKIP_LEADING_DELIMITERS")     return process_skip_leading_delimiters     (mped_info, param);
  else if(name == "SKIP_TRAILING_DELIMITERS")    return process_skip_trailing_delimiters    (mped_info, param);
  else if(name == "SKIP_CONSECUTIVE_DELIMITERS") return process_skip_consecutive_delimiters (mped_info, param);
  else if(name == "STUDY_ID")                    return process_study_id                    (mped_info, param);
  else if(name == "PEDIGREE_ID")                 return process_pedigree_id                 (mped_info, param);
  else if(name == "TREAT_AS_SIBS")               { set_treat_as_sibs(true); return true; }
  else if(name == "INDIVIDUAL_ID")               return process_individual_id               (mped_info, param);
  else if(name == "PARENT_ID")                   return process_parent_id                   (mped_info, param);
  else if(name == "SEX_FIELD")                   return process_sex_field                   (mped_info, param);
  else if(name == "NO_SEX_FIELD")                { set_no_sex_field(true); return true; }
  else if(name == "MARKER" || 
          name == "ALLELE" ||
          name == "TRAIT_MARKER" )               return process_marker                      (mped_info, param);
  else if(name == "PHENOTYPE" ||
          name == "TRAIT"     ||
          name == "COVARIATE")                   return process_phenotype                   (mped_info, param);
  else if(name == "STRING")                      return process_string                      (mped_info, param);
  else if(name == "MARKER_LIST")                 return process_marker_list                 (mped_info, param);
  else if(name == "COVARIATE_LIST")              return process_covariate_list              (mped_info, param);

  else return true;
}

/// Validates the format and file, then reads in the data.
///
bool
RefLSFDelimitedPedigreeFile::input(RefMultiPedigree&  p,
                                   const std::string& filename,
                                   ostream&           messages)
{
  if(!validate_file(filename) || !validate_format(filename))
    return false;

  if(!build_marker_list_parameters(p.info())) return false;
  if(!build_covariate_list_parameters(p.info())) return false;
  
  if( my_marker_covariates.size() )
  {
    if( !process_marker_covariates(p.info()) )
      return false;
  }

  return RefPedigreeFile::input(p, filename, messages);
}

/// Searches the list of MarkerListElement for one that begins with start.
///
/// \param start The start condition of the marker list
std::list<RefLSFDelimitedPedigreeFile::MarkerListElement>::iterator 
RefLSFDelimitedPedigreeFile::find_marker_list(string start)
{
  std::list<MarkerListElement>::iterator i = my_marker_lists.begin();
  
  for( ; i != my_marker_lists.end() && i->start_marker != start; ++i);
  
  return i;
}

bool
RefLSFDelimitedPedigreeFile::build_marker_list_parameters(RefMPedInfo &mped_info)
{
  string_tokenizer header(format());
  setup_tokenizer(header);
  
  for(string_tokenizer::iterator i = header.begin(); i != header.end(); ++i)
  {
    if(i->empty())
       continue;

    // Check to see if the field begins a marker list
    std::list<MarkerListElement>::iterator mlist_index = find_marker_list(toUpper(*i));
      
    if(mlist_index != my_marker_lists.end())
    {
      // If the field begins a marker_list, we can now add the fields to the map.
      if(!build_marker_list(header, i, mped_info, mlist_index))
        return false;
      
      my_marker_lists.erase(mlist_index);
    }
  }
  
  // Verify that all marker lists have been found.  If any remain, then there
  // is a problem!
  
  if(!my_marker_lists.empty())
  {
    for(std::list<MarkerListElement>::iterator i = my_marker_lists.begin();
        i != my_marker_lists.end(); ++i)
    {
      errors << priority(fatal) << "Marker list with start marker '"
             << i->start_marker << "' and end marker '"
             << i->end_marker << "' not found in pedigree file.   Please fix this and re-run S.A.G.E."
             << std::endl;
    }
    
    return false;
  }
  
  return true;
}

bool
RefLSFDelimitedPedigreeFile::build_marker_list(string_tokenizer&                      header,
                                               string_tokenizer::iterator             i, 
                                               RefMPedInfo&                           mped_info,
                                               std::list<MarkerListElement>::iterator mlist_index)
{
  MarkerListElement& mlist = *mlist_index;
  
  // Create the list of markers we must build by iterating through the 
  // header, looking for the end marker in the marker list.
  
  std::list<string> mlist_fields;
  
  for( ; i != header.end(); ++i)
  {
    // Check to make sure it's not also another parameter.  If it is, we've
    // got problems.
    field_map_type::const_iterator f = field_map().find( toUpper(*i) );
    
    if(f != field_map().end())
    {
      errors << priority(fatal) << "Overlap detected with marker_list with "
             << "field " << *i << ".  Marker list fields must be exclusive with "
             << "regards to other fields in the pedigree block.  Please fix this and re-run S.A.G.E."
             << std::endl;
      return false;
    }
    mlist_fields.push_back(toUpper(*i));
    
    if(toUpper(*i) == mlist.end_marker) break;
  }
  
  
  // If we didn't find the end, this is a fatal error
  if(i == header.end())
  {
    errors << priority(fatal) << "End marker '" << mlist.end_marker
           << "' of marker list starting with '" << mlist.start_marker 
           << "' not found in pedigree file fields following the start marker. "
           << "Please check your parameters."
           << std::endl;
           
    return false;
  }
  
  // Add the found markers to the marker list

  for(std::list<string>::const_iterator m = mlist_fields.begin(); m != mlist_fields.end(); ++m)
  {
    mlist.marker_params->attrs(true)->set(0, *m);
    if(!process_marker(mped_info, &*mlist.marker_params))
      return false;
  }
  return true;
}

std::list<RefLSFDelimitedPedigreeFile::MarkerListElement>::iterator 
RefLSFDelimitedPedigreeFile::find_covariate_list(string start)
{
  std::list<MarkerListElement>::iterator i = my_covariate_lists.begin();
  
  for( ; i != my_covariate_lists.end() && i->start_marker != start; ++i );
  
  return i;
}

bool
RefLSFDelimitedPedigreeFile::build_covariate_list_parameters(RefMPedInfo &mped_info)
{
  string_tokenizer header(format());
  setup_tokenizer(header);
  
  for( string_tokenizer::iterator i = header.begin(); i != header.end(); ++i )
  {
    if( i->empty() )
       continue;

    // Check to see if the field begins a covariate list
    std::list<MarkerListElement>::iterator clist_index = find_covariate_list(toUpper(*i));
      
    if( clist_index != my_covariate_lists.end() )
    {
      // If the field begins a covariate_list, we can now add the fields to the map.
      if( !build_covariate_list(header, i, mped_info, clist_index) )
        return false;
      
      my_covariate_lists.erase(clist_index);
    }
  }
  
  // Verify that all covariate lists have been found.  If any remain, then there
  // is a problem!
  
  if( !my_covariate_lists.empty() )
  {
    for( std::list<MarkerListElement>::iterator i = my_covariate_lists.begin();
        i != my_covariate_lists.end(); ++i )
    {
      errors << priority(fatal) << "Covariate list with start covariate '"
             << i->start_marker << "' and end covariate '"
             << i->end_marker << "' not found in pedigree file.   Please fix this and re-run S.A.G.E."
             << std::endl;
    }
    
    return false;
  }
  
  return true;
}

bool
RefLSFDelimitedPedigreeFile::build_covariate_list(string_tokenizer&                      header,
                                                  string_tokenizer::iterator             i, 
                                                  RefMPedInfo&                           mped_info,
                                                  std::list<MarkerListElement>::iterator clist_index)
{
  MarkerListElement& clist = *clist_index;
  
  // Create the list of covariates we must build by iterating through the 
  // header, looking for the end covariate in the covariate list.
  
  std::list<string> clist_fields;
  
  for( ; i != header.end(); ++i)
  {
    // Check to make sure it's not also another parameter.  If it is, we've
    // got problems.
    field_map_type::const_iterator f = field_map().find( toUpper(*i) );
    
    if( f != field_map().end() )
    {
      errors << priority(fatal) << "Overlap detected with covariate_list with "
             << "field " << *i << ".  Covariate list fields must be exclusive with "
             << "regards to other fields in the pedigree block.  Please fix this and re-run S.A.G.E."
             << std::endl;
      return false;
    }
    clist_fields.push_back(toUpper(*i));
    
    if( toUpper(*i) == clist.end_marker ) break;
  }
  
  
  // If we didn't find the end, this is a fatal error
  if( i == header.end() )
  {
    errors << priority(fatal) << "End covariate '" << clist.end_marker
           << "' of covariate list starting with '" << clist.start_marker 
           << "' not found in pedigree file fields following the start covariate. "
           << "Please check your parameters."
           << std::endl;
           
    return false;
  }
  
  // Add the found covariates to the covariate list

  for( std::list<string>::const_iterator m = clist_fields.begin(); m != clist_fields.end(); ++m )
  {
    clist.marker_params->attrs(true)->set(0, *m);
    if( !process_phenotype(mped_info, &*clist.marker_params) )
      return false;
  }
  return true;
}

/// Does LSFDelimited level processing of the "FORMAT" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_format2(RefMPedInfo &mped_info, const LSFBase *param)
{
  AttrVal v=attr_value(param, 0);
  if( v.has_value() )
    set_format_in_file(false);

  return true;
}

/// Processes the "WHITESPACE" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_whitespace(RefMPedInfo &mped_info, const LSFBase *param)
{
  AttrVal v=attr_value(param, 0);
  if( v.has_value() )
  {
    set_whitespace(v.String());
    DEBUG_RPEDFILE(errors << priority(debug) << "Found whitespace: '" << whitespace() << "'" << std::endl;)
  }

  return true;
}

/// Processes the "DELIMITER_MODE" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_delimiter_mode(RefMPedInfo &mped_info, const LSFBase *param)
{
  AttrVal v=attr_value(param, 0);

  if( v.has_value() )
  {
    if( toUpper(v.String()) == "SINGLE")
    {
      set_skip_trailing_delimiters    (false);
      set_skip_leading_delimiters     (false);
      set_skip_consecutive_delimiters (false);
    }
    else if( toUpper(v.String()) == "MULTIPLE")
    {
      set_skip_trailing_delimiters    (true);
      set_skip_leading_delimiters     (true);
      set_skip_consecutive_delimiters (true);
    }
    else
    {
      errors << priority(error) << "Unknown value for parameter 'delimiter_mode'" << std::endl;
    }

    DEBUG_RPEDFILE(errors << priority(debug) << "Found column delimited flag." << std::endl;)
  }

  return true;
}

/// Processes the "DELIMITERS" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_delimiters(RefMPedInfo &mped_info, const LSFBase *param)
{
  AttrVal v=attr_value(param, 0);
  if( v.has_value() )
  {
    stringstream all_del;
    all_del << mped_info.get_pheno_reader_info().get_allele_delimiter();
     
    if( all_del.str() == v.String() )
      errors << priority(error)
             << "The values for parameters 'delimiters' in pedigree block and 'allele_delimiter' in marker block should not be the same.  "
             << "The values for parameter 'marker' will be read incorrectly."
             << std::endl;

    set_delimiters(v.String());
    DEBUG_RPEDFILE(errors << priority(debug) << "Found delimiters: " << delimiters() << std::endl;)
  }

  return true;
}

/// Processes the "SKIP_LEADING_DELIMITERS" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_skip_leading_delimiters(RefMPedInfo &mped_info, const LSFBase *param)
{
  AttrVal v=attr_value(param, 0);
  if( v.has_value() )
  {
    if( toUpper(v.String()) == "TRUE")
      set_skip_leading_delimiters(true);
    else if( toUpper(v.String()) == "FALSE")
      set_skip_leading_delimiters(false);
    else
      errors << priority(error) << "Unknown value for parameter 'skip_leading_delimiters'" << std::endl;
  }

  return true;
}

/// Processes the "SKIP_TRAILING_DELIMITERS" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_skip_trailing_delimiters(RefMPedInfo &mped_info, const LSFBase *param)
{
  AttrVal v=attr_value(param, 0);
  if( v.has_value() )
  {
    if( toUpper(v.String()) == "TRUE")
      set_skip_trailing_delimiters(true);
    else if( toUpper(v.String()) == "FALSE")
      set_skip_trailing_delimiters(false);
    else
      errors << priority(error) << "Unknown value for parameter 'skip_trailing_delimiters'" << std::endl;
  }

  return true;
}

/// Processes the "SKIP_CONSECUTIVE_DELIMITERS" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_skip_consecutive_delimiters(RefMPedInfo &mped_info, const LSFBase *param)
{
  AttrVal v=attr_value(param, 0);
  if( v.has_value() )
  {
    if( toUpper(v.String()) == "TRUE")
      set_skip_consecutive_delimiters(true);
    else if( toUpper(v.String()) == "FALSE")
      set_skip_consecutive_delimiters(false);
    else
      errors << priority(error) << "Unknown value for parameter 'skip_consecutive_delimiters'" << std::endl;
  }

  return true;
}

/// Processes the "STUDY_ID" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_study_id(RefMPedInfo &mped_info, const LSFBase *param)
{
  AttrVal v=attr_value(param, 0);
  if( v.has_value() )
  {
    add_study_id_field(v.String());
    DEBUG_RPEDFILE(errors << priority(debug) << "Found study id field = " << v.String() << std::endl;)
  }

  return true;
}
/// Processes the "PEDIGREE_ID" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_pedigree_id(RefMPedInfo &mped_info, const LSFBase *param)
{
  AttrVal v=attr_value(param, 0);
  if( v.has_value() )
  {
    add_pedigree_id_field(v.String());
    set_pedigree_id_name(v.String());
    
    DEBUG_RPEDFILE(errors << priority(debug) << "Found pedigree id field = " << v.String() << std::endl;)
  }

  return true;
}

/// Processes the "INDIVIDUAL_ID" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_individual_id(RefMPedInfo &mped_info, const LSFBase *param)
{
  AttrVal v=attr_value(param, 0);
  if( v.has_value() )
  {
    add_individual_id_field(v.String());
    DEBUG_RPEDFILE(errors << priority(debug) << "Found individual id field = " << v.String() << std::endl;)
  }

  return true;
}

/// Processes the "PARENT_ID" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_parent_id(RefMPedInfo &mped_info, const LSFBase *param)
{
  AttrVal v=attr_value(param, 0);
  if( v.has_value() )
  {
    add_parent_id_field(v.String());
    DEBUG_RPEDFILE(errors << priority(debug) << "Found parent id field = " << v.String() << std::endl;)
  }

  return true;
}

/// Processes the "SEX_FIELD" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_sex_field(RefMPedInfo &mped_info, const LSFBase *param)
{
  AttrVal v=attr_value(param, 0);
  if( v.has_value() )
  {
    add_sex_field(v.String());
    set_sex_field_name(v.String());

    DEBUG_RPEDFILE(errors << priority(debug) << "Found sex field = " << v.String() << std::endl;)

    RefLSFPedigreeFile::process_sex_code(mped_info, param);
  }

  return true;
}

/// Processes the "MARKER", "TRAIT_MARKER", and "ALLELE" parameters
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_marker(RefMPedInfo &mped_info, const LSFBase *param)
{
  string name = toUpper(param->name());
  string field_name;

  if( param->attrs() )
    field_name = param->attrs()->StringAttr(0);

  if(!field_name.size())
  {
    errors << priority(warning) << "Marker with no field name specified.  Skipping." << std::endl;
    return true;
  }
  
  string marker_name = param->attrs()->StringAttr("name");

  if(!marker_name.size())
    marker_name = field_name;

  if( mped_info.get_pheno_reader_info().get_covariate_moi() != "" )
  {
    add_marker_covariate(mped_info, param, field_name, marker_name);

    return true;
  }

  // Test to see if the marker should be skipped, and return if so.
  if(test_skip_marker(param)) 
    return true;

  // Find the marker's model, if it exists
  size_t marker_id = mped_info.marker_find(marker_name);

  // Does it exist?
  bool has_model = marker_id < mped_info.marker_count();
  
  if(!has_model)
  {
    // Check if it should be read as covariate, not marker.

    if( use_as_covariate(param) )
    {
      add_marker_covariate(mped_info, param, field_name, marker_name);

      return true;
    }

    // We can only continue if dynamic parsing is allowed, so test that.
    if(!test_allow_dynamic(param))
    {
        errors << priority(warning) << "Marker '" << marker_name
               << "' not found.  Skipping." << std::endl;
        return true;
    }

    marker_id = setup_dynamic_marker(mped_info, param, marker_name);
  }
  else // This has a model, so we're not dynamic
  {
    // Check to see if there's any special allele frequency modifications
    if(param->attrs()->has_attr("equal_allele_freq")      ||
       param->attrs()->has_attr("complement_allele_freq") ||
       param->attrs()->has_attr("compl_allele_freq")      ||
       param->attrs()->has_attr("minimum_allele_freq")    ||
       param->attrs()->has_attr("maximum_allele_freq")    ||
       param->attrs()->has_attr("equal")                  ||
       param->attrs()->has_attr("complement")             ||
       param->attrs()->has_attr("minimum")                ||
       param->attrs()->has_attr("maximum") )
    {
      parse_allele_frequency_adjustment(mped_info.marker_info(marker_id), param);
    }
  }
  
  parse_marker_parameters(mped_info, param, marker_id);

  // Finally, we have set up everything we need to set up for adding the marker.
  if(name == "ALLELE")
    add_allele_field( field_name, marker_name );
  else
    add_marker_field( field_name, marker_name );

  return true;
}

/// Processes the "TRAIT", "COVARIATE", and "PHENOTYPE" parameters
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_phenotype(RefMPedInfo &mped_info, const LSFBase *param)
{
  AttrList::const_iterator a;

  std::string name       = toUpper(param->name()),
              field_name = "";

  if( param->attrs() )
  {
    field_name = param->attrs()->StringAttr(0);
  }

  if(!field_name.size())
  {
    errors << priority(warning) << "Trait with no field name specified.  Skipping." << std::endl;
    return true;
  }

  string trait_name;

  trait_name = param->attrs()->StringAttr("name");

  if(!trait_name.size())
    trait_name = field_name;

  bool skip = skip_traits();

  if( !force_skip_traits() && has_attr(param, "skip") )
  {
    skip = true;

    if( toUpper(attr_value(param, "skip").String()) == "FALSE" )
    {
      skip = false;
    }
  }

  if(skip)
  {
    return true;
  }

  add_trait_field( field_name, trait_name );

  size_t t         = mped_info.trait_find(trait_name);
  bool   has_model = t < mped_info.trait_count();

  if(has_model)
  {
    return true;
  }

  bool trvar = param->attrs()->has_attr("variate")   || name == "TRAIT",
       cov   = param->attrs()->has_attr("covariate") || name == "COVARIATE";

  RefTraitInfo::trait_use use = RefTraitInfo::unknown_use;

  if(trvar && cov)
  {
    errors << priority(warning) << "Trait '" << trait_name
           << "' has multiple types listed.  Setting to unknown." << std::endl;
  }
  else if(trvar)
  {
    use = RefTraitInfo::trait_variate;
  }
  else if(cov)
  {
    use = RefTraitInfo::trait_covariate;
  }

  // Binary!
  if( param->attrs()->has_attr("binary"))
  {
    std::string affected   = (a = param->attrs()->find("affected"))   == param->attrs()->end() ? "1" : strip_ws(a->second.String()),
                unaffected = (a = param->attrs()->find("unaffected")) == param->attrs()->end() ? "0" : strip_ws(a->second.String()),
                missing    = (a = param->attrs()->find("missing"))    == param->attrs()->end() ? ""  : strip_ws(a->second.String());
    double      threshold  = (a = param->attrs()->find("threshold"))  == param->attrs()->end() ? std::numeric_limits<double>::quiet_NaN() : a->second.Real();
    
    if( affected == unaffected )
    {
      errors << priority(error) << "Pedigree File: Affected and unaffected codes "
             << "must be specified and may not be the same." << std::endl;
      skip = true;
    }
    else if( affected == missing || unaffected == missing )
    {
      errors << priority(error) << "Pedigree File Error: Affected and unaffected codes "
             << "must not be the same as the missing value code." << std::endl;
      skip = true;
    }
    else
    {
      DEBUG_RPEDFILE(errors << priority(debug) << "Found binary trait = " << v.String() << std::endl;)

      size_t t = mped_info.add_binary_trait( trait_name, use );

      mped_info.trait_info(t).set_string_missing_code              (missing);
      mped_info.trait_info(t).set_numeric_missing_code    (str2doub(missing));
      mped_info.trait_info(t).set_string_affected_code             (affected);
      mped_info.trait_info(t).set_string_unaffected_code           (unaffected);
      mped_info.trait_info(t).set_numeric_affected_code   (str2doub(affected));
      mped_info.trait_info(t).set_numeric_unaffected_code (str2doub(unaffected));
      mped_info.trait_info(t).set_threshold                        (threshold);
    }
  }

  // Categorical!
  else if( param->attrs()->has_attr("categorical"))
  {
    std::string    missing    = (a = param->attrs()->find("missing")) == param->attrs()->end() ? ""  : strip_ws(a->second.String()),
                   values     = (a = param->attrs()->find("values"))  == param->attrs()->end() ? ""  : strip_ws(a->second.String());
    RefTraitInfo & tinfo      = mped_info.trait_info(mped_info.add_trait(trait_name, RefTraitInfo::categorical_trait, use));

    tinfo.set_string_missing_code  (missing);
    tinfo.set_numeric_missing_code (str2doub(missing));
    tinfo.set_lockout              (false);

    if(values != "")
    {
      UTIL::StringUtils::splitMultiDelimitedString(values, " \t\n,", tinfo.get_categories());
      tinfo.set_lockout(true);
    }
  }

  // Doesn't have binary or categorical attribute; must be continuous!
  else
  {
    string missing;

    if( (a=param->attrs()->find("missing"))  != param->attrs()->end() )
      missing = strip_ws(a->second.String());

    DEBUG_RPEDFILE(errors << priority(debug) << "Found continuous trait = " << trait_name << std::endl;)

    size_t t = mped_info.add_continuous_trait( trait_name, use );

    mped_info.trait_info(t).set_string_missing_code(missing);
    mped_info.trait_info(t).set_numeric_missing_code(str2doub(missing));
  }
  if(skip)
  {
    DEBUG_RPEDFILE(errors << priority(debug) << "Skipping trait = " << trait_name << std::endl;)
    mped_info.add_trait( trait_name, RefTraitInfo::invalid_trait );
  }

  return true;
}

/// Processes the "STRING" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_string(RefMPedInfo &mped_info, const LSFBase *param)
{
  string field_name;

  if( param->attrs() )
    field_name = param->attrs()->StringAttr(0);

  if(!field_name.size())
  {
    errors << priority(warning) << "String field with no field name specified.  Skipping." << std::endl;
    return true;
  }

  string string_name = param->attrs()->StringAttr("name");
  if(!string_name.size())
    string_name = field_name;

  bool skip = false;
  if( has_attr(param, "skip") )
  {
    skip = true;
    if( toUpper(attr_value(param, "skip").String()) == "FALSE" )
      skip = false;
  }

  if(skip)
    return true;

  add_string_field(field_name, string_name);
  mped_info.add_string_field(string_name);

  return true;
}

/// Processes the "MARKER_LIST" parameter
///
/// Note that because the format is generally in the file, the
/// markers cannot be added to the set of fields at this time.  This
/// isn't done until the input() method is called.
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_marker_list(RefMPedInfo& mped_info, const LSFBase* param)
{
  // Test for required attributes
  if(!param->attrs() || !has_attr(param,"start") || !has_attr(param,"end"))
  {
    errors << priority(fatal) << "marker_list parameter requires both start "
           << "and end attributes specified.  Please check your file and re-run S.A.G.E."
           << std::endl;
    return false;
  }
  
  // Make sure the required attributes have content
  string start_marker = toUpper(param->attrs()->StringAttr("start"));
  string end_marker = toUpper(param->attrs()->StringAttr("end"));

  if(start_marker.empty() || end_marker.empty())
  {
    errors << priority(fatal) << "marker_list parameter requires both start "
           << "and end attributes specified.  Please check your file and re-run S.A.G.E."
           << std::endl;
    return false;
  }

  // Make sure there isn't a marker list with this start already.
  //
  // Note:  This may not be necessary.  Marker lists which overlap are detected
  //        when their markers are inserted into the set of fields.  There
  //        may not be any reason to do this check.  But it's a better, and more
  //        informative error message, so it is left for that reason.
  if(find_marker_list(start_marker) != my_marker_lists.end())
  {
    errors << priority(fatal) << "Overlapping marker lists detected with duplicate first field "
           << start_marker << ".  Marker lists may not overlap.  Please correct this "
           << "and re-run S.A.G.E." << endl;
    return false;
  }
  
  // Create a new MarkerListElement for storing the list information until it can
  // be used.
  my_marker_lists.push_back(MarkerListElement());
  
  MarkerListElement& new_mlist = my_marker_lists.back();
  
  new_mlist.start_marker  = start_marker;
  new_mlist.end_marker    = end_marker;
  
  // We copy the attributes into a new LSFBase object so that they can be
  // used to set up the markers when it is time.
  new_mlist.marker_params = new LSFBase("MARKER_LIST");
  
  (*new_mlist.marker_params->attrs(true)) = (*param->attrs());

  return true;
}

bool 
RefLSFDelimitedPedigreeFile::process_covariate_list(RefMPedInfo& mped_info, const LSFBase* param)
{
  // Test for required attributes
  if( !param->attrs() || !has_attr(param,"start") || !has_attr(param,"end" ))
  {
    errors << priority(fatal) << "covariate_list parameter requires both start "
           << "and end attributes specified.  Please check your file and re-run S.A.G.E."
           << std::endl;
    return false;
  }
  
  // Make sure the required attributes have content
  string start_covariate = toUpper(param->attrs()->StringAttr("start"));
  string end_covariate = toUpper(param->attrs()->StringAttr("end"));

  if( start_covariate.empty() || end_covariate.empty() )
  {
    errors << priority(fatal) << "covariate_list parameter requires both start "
           << "and end attributes specified.  Please check your file and re-run S.A.G.E."
           << std::endl;
    return false;
  }

  // Make sure there isn't a covariate list with this start already.
  //
  // Note:  This may not be necessary.  Marker lists which overlap are detected
  //        when their markers are inserted into the set of fields.  There
  //        may not be any reason to do this check.  But it's a better, and more
  //        informative error message, so it is left for that reason.
  if( find_covariate_list(start_covariate) != my_covariate_lists.end() )
  {
    errors << priority(fatal) << "Overlapping covariate lists detected with duplicate first field "
           << start_covariate << ".  Covariate lists may not overlap.  Please correct this "
           << "and re-run S.A.G.E." << endl;
    return false;
  }
  
  // Create a new MarkerListElement for storing the list information until it can
  // be used.
  my_covariate_lists.push_back(MarkerListElement());
  
  MarkerListElement& new_clist = my_covariate_lists.back();
  
  new_clist.start_marker  = start_covariate;
  new_clist.end_marker    = end_covariate;
  
  // We copy the attributes into a new LSFBase object so that they can be
  // used to set up the covariates when it is time.
  new_clist.marker_params = new LSFBase("COVARIATE_LIST");
  
  (*new_clist.marker_params->attrs(true)) = (*param->attrs());

  return true;
}

/// Tests to see if a particular marker should be skipped.
///
/// The test is as follows:
///  -# If skip_markers() is \c true, we skip *unless* force_skip_markers() is
///     \c false and the parameter has a "skip=FALSE" attribute.
///  -# If skip_markers() is \c false, we skip *only* if force_skip_markers() is
///     \c false and the parameter has a "skip" attribute which is *not* "FALSE".
///
/// \param param The parameter containing the marker information
/// \returns \c true if the marker should be skipped, \c false otherwise
bool
RefLSFDelimitedPedigreeFile::test_skip_marker(const LSFBase* param) const
{
  // Get Default skip status
  bool skip = skip_markers();

  // Determine if skipping markers should be done
  if( !force_skip_markers() && has_attr(param, "skip") )
  {
    skip = true;
    if( toUpper(attr_value(param, "skip").String()) == "FALSE" )
      skip = false;
  }

  return skip;
}

/// Tests to see if a particular marker can be dynamic.
///
/// The test is as follows:
///  -# If dynamic_markers() is \c true, we allow dynamic *unless* 
///     force_dynamic_markers() is
///     \c false and the parameter has a "dynamic=FALSE" attribute.
///  -# If dynamic_markers() is \c false, we allow dynamic *only* if 
///     force_skip_markers() is
///     \c false and the parameter has a "dynamic" attribute which is *not* "FALSE".
///
/// \param param The parameter containing the marker information
/// \returns \c true if the marker is allowed to be dynamic, \c false otherwise
bool
RefLSFDelimitedPedigreeFile::test_allow_dynamic(const LSFBase* param) const
{
  // Get default dynamic status
  bool dynamic = dynamic_markers();

  if( !force_dynamic_markers() && has_attr(param, "dynamic") )
  {
    dynamic = true;
    if( toUpper(attr_value(param, "dynamic").String()) == "FALSE" )
      dynamic = false;
  }
  
  return dynamic;
}

/// Sets up a dynamic marker when one is parsed
///
/// \param param The parameter containing the dynamic marker
/// \returns The id of the new dynamic marker
size_t
RefLSFDelimitedPedigreeFile::setup_dynamic_marker
    (RefMPedInfo&   mped_info,
     const LSFBase* param,
     const string&  marker_name)
{
  size_t marker_id = mped_info.add_marker(marker_name);
  mped_info.marker_info(marker_id).gmodel().set_dynamic_alleles(true);

  const PhenotypeReaderInfo& pr = mped_info.get_pheno_reader_info();
  
  char   sep     = pr.get_allele_delimiter();
  string missing = pr.get_allele_missing();
  mped_info.marker_info(marker_id).gmodel().set_unphased_separator(sep);

  parse_dynamic_markers_missing(mped_info.marker_info(marker_id), missing);

  return marker_id;
}

/// This function parses the basic marker parameters which don't require a lot
/// of special setup.  This includes the x/y linked, delimiters, and missing
/// parameters.
void
RefLSFDelimitedPedigreeFile::parse_marker_parameters
    (RefMPedInfo &mped_info,
     const LSFBase* param,
     size_t marker_id)
{
  AttrList::const_iterator a;
  
  if(    (a=param->attrs()->find("allele_missing")) != param->attrs()->end()
      || (a=param->attrs()->find("missing"))        != param->attrs()->end() )
  {
    string global_missing = mped_info.get_pheno_reader_info().get_allele_missing();
    string local_missing  = strip_ws(a->second.String());

    if( global_missing != local_missing )
      update_marker_missing_info(mped_info.marker_info(marker_id), local_missing);
  }

  if( param->attrs()->has_attr("x_linked") )
  {
    mped_info.marker_info(marker_id).set_model_type(MLOCUS::X_LINKED);
    set_sex_linked_exist(true);
  }
  else if( param->attrs()->has_attr("y_linked") )
  {
    mped_info.marker_info(marker_id).set_model_type(MLOCUS::Y_LINKED);
    set_sex_linked_exist(true);
  }

  if(    (a=param->attrs()->find("allele_delimiter")) != param->attrs()->end()
      || (a=param->attrs()->find("delimiter"))        != param->attrs()->end() )
  {
    char global_sep = mped_info.get_pheno_reader_info().get_allele_delimiter();
    char local_sep  = global_sep;

    if( a->second.String().size() == 1 )
      local_sep = a->second.String()[0];
    else
      local_sep = strip_ws(a->second.String())[0];

    if( global_sep != local_sep )
      update_marker_delimiter_info(mped_info.marker_info(marker_id), local_sep);
  }
}

bool
RefLSFDelimitedPedigreeFile::use_as_covariate(const LSFBase* param)
{
  AttrList::const_iterator a;

  if(    (a=param->attrs()->find("covariate_moi")) == param->attrs()->end()
      && (a=param->attrs()->find("cov_moi"))       == param->attrs()->end()
      && (a=param->attrs()->find("cov_func"))      == param->attrs()->end()
      && (a=param->attrs()->find("covariate_function")) == param->attrs()->end() )
  {
    return false;
  }
  else
  {
    a = param->attrs()->find("covariate_moi");
    if( a == param->attrs()->end() )
      a = param->attrs()->find("cov_moi");
    if( a == param->attrs()->end() )
      a = param->attrs()->find("cov_func");
    if( a == param->attrs()->end() )
      a = param->attrs()->find("covariate_function");

    if( !a->second.String().size() )
      return false;
  }

  if(    (a=param->attrs()->find("covariate_allele")) == param->attrs()->end()
      && (a=param->attrs()->find("cov_allele"))       == param->attrs()->end()
      && (a=param->attrs()->find("base_allele"))      == param->attrs()->end() )
  {
    return false;
  }
  else
  {
    a = param->attrs()->find("covariate_allele");
    if( a == param->attrs()->end() )
      a = param->attrs()->find("cov_allele");
    if( a == param->attrs()->end() )
      a = param->attrs()->find("base_allele");

    if( !a->second.String().size() )
      return false;
  }

  return true;
}

void
RefLSFDelimitedPedigreeFile::add_marker_covariate(RefMPedInfo&   mped_info,
                                                  const LSFBase* param,
                                                  const string&  field_name,
                                                  const string&  marker_name)
{
  string cmoi       = mped_info.get_pheno_reader_info().get_covariate_moi();
  string callele    = mped_info.get_pheno_reader_info().get_covariate_allele();
  bool   allow_hemi = mped_info.get_pheno_reader_info().get_allow_hemizygote();
  string missing    = mped_info.get_pheno_reader_info().get_allele_missing();
  char   sep        = mped_info.get_pheno_reader_info().get_allele_delimiter();

  AttrList::const_iterator a;

  if(    (a=param->attrs()->find("covariate_moi")) != param->attrs()->end()
      || (a=param->attrs()->find("cov_moi"))       != param->attrs()->end()
      || (a=param->attrs()->find("covvariate_function")) != param->attrs()->end()
      || (a=param->attrs()->find("cov_func"))            != param->attrs()->end() )
  {
    string local_cmoi = toUpper(strip_ws(a->second.String()));

    if( local_cmoi == "ADD" || local_cmoi == "ADDITIVE" || local_cmoi == "A" )
      local_cmoi = "_ADD";
    else if( local_cmoi == "DOM" || local_cmoi == "DOMINANT" || local_cmoi == "D" )
      local_cmoi = "_DOM";
    else if( local_cmoi == "REC" || local_cmoi == "RECESSIVE" || local_cmoi == "R" )
      local_cmoi = "_REC";
    else
    {
      errors << priority(error)
             << "Invalid value for parameter 'covariate_moi'.  Using default..." << endl;
      local_cmoi = "_ADD";
    }

    if( cmoi != local_cmoi )
    {
      //cout << "different cov_moi for " << marker_name << endl;
      cmoi = local_cmoi;
    }
  }

  if(    (a=param->attrs()->find("covariate_allele")) != param->attrs()->end()
      || (a=param->attrs()->find("cov_allele"))       != param->attrs()->end()
      || (a=param->attrs()->find("base_allele"))      != param->attrs()->end() )
  {
    string local_callele = strip_ws(a->second.String());

    if( callele != local_callele )
    {
      //cout << "different cov_allele for " << marker_name << endl;
      callele = local_callele;
    }
  }

  if(    (a=param->attrs()->find("allow_hemizygote")) != param->attrs()->end()
      || (a=param->attrs()->find("allow_hemi"))       != param->attrs()->end() )
  {
    string hemi = toUpper(strip_ws(a->second.String()));

    bool local_hemi = false;
    if( hemi == "TRUE" || hemi == "YES" )
      local_hemi = true;

    if( allow_hemi != local_hemi )
    {
      //cout << "different allow_hemizygote for " << marker_name << endl;
      allow_hemi = local_hemi;
    }
  }

  if(    (a=param->attrs()->find("allele_missing")) != param->attrs()->end()
      || (a=param->attrs()->find("missing"))        != param->attrs()->end() )
  {
    string local_missing = strip_ws(a->second.String());

    if( missing != local_missing )
    {
      //cout << "different allele_missing for " << marker_name << endl;
      missing = local_missing;
    }
  }

  if(    (a=param->attrs()->find("allele_delimiter")) != param->attrs()->end()
      || (a=param->attrs()->find("delimiter"))        != param->attrs()->end() )
  {
    char local_sep = sep;

    if( a->second.String().size() == 1 )
      local_sep = a->second.String()[0];
    else
      local_sep = strip_ws(a->second.String())[0];

    if( sep != local_sep )
    {
      //cout << "different allele_delimiter for " << marker_name << endl;
      sep = local_sep;
    }
  }

  string mcov_name = marker_name + cmoi + "_" + callele;

#if 0
  cout << "field_name       = " << field_name << endl
       << "marker_name      = " << marker_name << endl
       << "mcov_name        = " << mcov_name << endl
       << "allele_delimiter = " << sep << endl
       << "allele_missing   = " << missing << endl
       << "covariate_moi    = " << cmoi << endl
       << "covariate_allele = " << callele << endl
       << "allow_hemizygote = " << allow_hemi << endl;
#endif

  if( toUpper(param->name()) == "ALLELE" )
    add_allele_cov_field(field_name, mcov_name);
  else
    add_marker_cov_field(field_name, mcov_name);

  if( my_marker_covariates.find(mcov_name) != my_marker_covariates.end() )
    return;

  marker_covariate_info new_mcov;

  new_mcov.type             = toUpper(param->name());
  new_mcov.field_name       = field_name;
  new_mcov.mcov_name        = mcov_name;
  new_mcov.allele_delimiter = sep;
  new_mcov.allele_missing   = missing;
  new_mcov.covariate_moi    = cmoi;
  new_mcov.covariate_allele = callele;
  new_mcov.allow_hemizygote = allow_hemi;
  new_mcov.trait_index      = (size_t)-1;

  my_marker_covariates[mcov_name] = new_mcov;

  return;
}

bool
RefLSFDelimitedPedigreeFile::process_marker_covariates(RefMPedInfo &mped_info)
{
  std::map<string, marker_covariate_info>::iterator mi     = my_marker_covariates.begin();
  std::map<string, marker_covariate_info>::iterator mi_end = my_marker_covariates.end();
  
  for( ; mi != mi_end; ++mi )
  {
    string mcov_name = mi->second.mcov_name;

    size_t t = mped_info.trait_find(mcov_name);

    if( t < mped_info.trait_count() )
      continue;

    t = mped_info.add_continuous_trait(mcov_name, RefTraitInfo::trait_covariate);

    mi->second.trait_index = t;

    mped_info.trait_info(t).set_string_missing_code("?");
    //mped_info.trait_info(t).set_string_missing_code(mi->second.allele_missing);
    mped_info.trait_info(t).set_numeric_missing_code(numeric_limits<double>::quiet_NaN());
  }

  return true;
}

string
RefLSFDelimitedPedigreeFile::get_marker_covariate_value(string mcov_name, const string& v1, const string& v2)
{
  const marker_covariate_info& mc_info = my_marker_covariates[mcov_name];

#if 0
  cout << "type             = " << mc_info.type << endl
       << "field_name       = " << mc_info.field_name << endl
       << "mcov_name        = " << mc_info.mcov_name << endl
       << "allele_delimiter = " << mc_info.allele_delimiter << endl
       << "allele_missing   = " << mc_info.allele_missing << endl
       << "covariate_moi    = " << mc_info.covariate_moi << endl
       << "covariate_allele = " << mc_info.covariate_allele << endl
       << "allow_hemizygote = " << mc_info.allow_hemizygote << endl
       << "trait_index      = " << mc_info.trait_index << endl
       << "v1 = " << v1 << ", v2 = " << v2 << endl;
#endif

  string allele1 = v1;
  string allele2 = v2;

  if( !v2.size() )
  {
    size_t sep = v1.find_first_of(mc_info.allele_delimiter);
    allele1 = v1.substr(0, sep);
    allele2 = v1.substr(sep+1, v1.size()-(sep+1));
  }

  if( allele1 == mc_info.allele_missing && allele2 == mc_info.allele_missing )
    return "?"; //mc_info.allele_missing;
  else if( allele1 == mc_info.allele_missing || allele2 == mc_info.allele_missing )
  {
    if( mc_info.allow_hemizygote )
    {
      if( allele1 == mc_info.allele_missing )
        allele1 = allele2;
      else
        allele2 = allele1;
    }
    else
      return "?"; //mc_info.allele_missing;
  }

  size_t allele_count = 0;

  if( allele1 == mc_info.covariate_allele )
    ++allele_count;

  if( allele2 == mc_info.covariate_allele )
    ++allele_count;

  if( allele_count == 1 )
  {
    if( mc_info.covariate_moi == "DOM" )
      allele_count = 2;
    else if( mc_info.covariate_moi == "REC" )
      allele_count = 0;
  }
#if 0
  cout << "a1 = " << allele1 << ", a2 = " << allele2 << ", a_count = " << allele_count << endl;
#endif

  return long2str(allele_count);
}

} // End namespace RPED
} // End namespace SAGE
