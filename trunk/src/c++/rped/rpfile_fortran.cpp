// Copyright(c) 1998 RC Elston
#include <string>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <iomanip>
#include "fortran/Tokenizer.h"
#include "fortran/Token_func.h"
#include "mped/sp.h"
#include "rped/rpfile.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "util/StringUtils.h"

#define DEBUG_RPEDFILE(x)

namespace SAGE {
namespace RPED {


RefFortranPedigreeFile::RefFortranPedigreeFile(cerrorstream &err)
       : RefPedigreeFile(err)
{}

RefFortranPedigreeFile::~RefFortranPedigreeFile()
{}

bool RefFortranPedigreeFile::
build_fields(const RefMPedInfo &mped_info, bool quiet)
{
  reset_counts();

  field_list_type::iterator i;
  for(i = field_list().begin(); i != field_list().end(); ++i)
  {
    if( i->type == trait )
    {
      size_t t = mped_info.trait_find(i->name);

      if(t < mped_info.trait_count() )
        i->index = t;
      else
      {
        i->index = (size_t)-1;

        errors << priority(warning) << "Invalid trait '" << i->name
               << "'.  Skipping..." << endl;
      }
    }
    else if( i->type == string_field )
    {
      size_t t = mped_info.string_find(i->name);

      if(t < mped_info.string_count() )
        i->index = t;
      else
      {
        i->index = (size_t)-1;

        errors << priority(warning) << "Invalid string field '" << i->name
               << "'.  Skipping..." << endl;
      }
    }
    else if( i->type == marker || i->type == allele )
    {
      size_t m = mped_info.marker_find(i->name);

      if(m < mped_info.marker_count() )
        i->index = m;
      else
      {
        i->index = (size_t)-1;

        if(mped_info.marker_count())
          errors << priority(warning) << "Invalid marker '" << i->name
                 << "'.  Skipping..." << endl;
      }
    }
  }

  return true;
}

bool RefFortranPedigreeFile::input_pedigree(RefMultiPedigree &p, const string &filename,
                                   ostream &messages, bool quiet)
{
  if( !filename.size() )
  {
    errors << priority(critical) << "No Family Data file specified." << endl;
    return false;
  }

  if( !format().size() )
  {
    errors << priority(critical)
           << "No FORTRAN Format statement given for Family Data File '"
           << filename << "'." << endl;
    return false;
  }

  std::ifstream infile( filename.c_str() );

  if(!infile.good())
  {
    errors << priority(critical) << "Unable to open Family Data file '" << filename
                    << "'. Please check your file." << endl;
    invalidate();
    return false;
  }

  const RefMPedInfo &mped_info = p.info();

  if(!build_fields(mped_info, quiet))
  {
    errors << priority(error)
           << "Cannot build list of fields to read.  Aborting..." << endl;
    invalidate();
    return false;
  }

  if(!validate_fields(false, quiet))
  {
    errors << priority(error)
           << "Invalid list of fields to read.  Aborting." << endl;
    invalidate();
    return false;
  }

  FortranFormatter ffortran;
  ffortran.grab(); // Make sure this doesn't run off
  Tokenizer tokens(&ffortran, infile);
  Tokenizer::token tok;

  ffortran.set_format( format() );

#if 0
  size_t verbose = verbose_output();
#endif

  field_list_type::const_iterator field_info      = my_fields.begin();
  field_list_type::const_iterator field_info_end  = my_fields.end();

#if 0
  if(verbose)
    print_family_structure_header(messages, filename);
#endif

  // Figure out if there's supposed to be a pedigree id field:
  
  bool has_pedigree_id_field = false;
  
  for(field_list_type::const_iterator itr = my_fields.begin(); itr != my_fields.end(); ++itr)
  {
    if(itr->type == pedigree_id)
      has_pedigree_id_field = true;
  }

  while(infile.good())
  {
    tokens.read_line();
    ffortran.reset();

    if(tokens.eof()) break;

    size_t record_starting_line = tokens.line_count();

    std::string ped_name = "",                           // Pedigree number
                ind_name = "",                           // Individual id
                parent1  = "",
                parent2  = "",                           // Parents
                sex      = mped_info.sex_code_unknown(); // Sex code

    for(field_info = my_fields.begin(); field_info != field_info_end ; ++field_info)
    {
      switch(field_info->type)
      {
        case pedigree_id   : ped_name = strip_ws(get_token(tok, &tokens).String());   break;
        case individual_id : ind_name = strip_ws(get_token(tok, &tokens).String());   break;
        case parent_id     : if(!parent1.size())
                               parent1 = strip_ws(get_token(tok, &tokens).String());
                             else
                               parent2 = strip_ws(get_token(tok, &tokens).String()); break;
        case sex_code      : sex = strip_ws(get_token(tok, &tokens).String());       break;
        default            : get_token(tok, &tokens);                                break;
      }
    }

    // Check for treat_ped_id options:
    
    if(has_pedigree_id_field == false) // If there isn't supposed to be a pedigree field, given them a dummy pedigree name:
    {
      ped_name = "0";
    }
    else if(has_pedigree_id_field && get_treat_as_sibs() == true) // There's a pedigree_id field and treat_as_sibs is enabled:
    {
      parent1 = ped_name + "_parent1";
      parent2 = ped_name + "_parent2";

      p.add_member(ped_name, parent1, MPED::SEX_MALE);
      p.add_member(ped_name, parent2, MPED::SEX_FEMALE);
    }

    // If parents are still empty for some reason, set them to be missing:
    parent1 = parent1.size() ? parent1 : mped_info.individual_missing_code();
    parent2 = parent2.size() ? parent2 : mped_info.individual_missing_code();

#if 0
    if(verbose)
    {
      --verbose;
      print_family_structure(messages, ped_name, ind_name, sex, parent1, parent2);
    }
#endif

    DEBUG_RPEDFILE(errors << priority(debug) << "Found (" << pn  << "," << id  << "," << sex << "," << p1 << "," << p2  << ")" << endl; )

//    if(record_starting_line % 100 == 0)
//      cout << "Found (" << pn  << "," << id  << "," << sex
//                            << "," << p1 << "," << p2  << ")" << endl;

    add_member(p, ped_name, ind_name, sex, parent1, parent2, record_starting_line);
  }

  infile.close();

#if 0
  if(verbose != verbose_output())
    print_family_structure_footer(messages);
#endif

  if(!build_pedigree(p))
  {
    errors << priority(error)
           << "Fatal error building pedigree data structure.  Aborting." << endl;
    return false;    
  }

  return true;
}

bool RefFortranPedigreeFile::input_data(RefMultiPedigree &p, const string &filename,
                                   ostream &messages, bool quiet)
{
  // Check to see if a second pass is necessary
  // (traits and covariates were defered in the first pass)
  if(trait_count() == 0 && !sex_code_trait() && !pedigree_id_trait() && marker_count() == 0)
    return true;

  if( !filename.size() )
  {
    errors << priority(critical) << "No Family Data file specified." << endl;
    return false;
  }

  if( !format().size() )
  {
    errors << priority(critical)
           << "No FORTRAN Format statement given for Family Data File '"
           << filename << "'." << endl;
    return false;
  }

  std::ifstream infile( filename.c_str() );

  if(!infile.good())
  {
    errors << priority(critical) << "Unable to open Family Data file '"
           << filename
           << "' to read traits and/or markers. Please check your file."
           << endl;
    invalidate();
    return false;
  }

  RefMPedInfo &mped_info = p.info();

  if(!build_fields(mped_info, quiet))
  {
    errors << priority(error)
           << "Cannot build list of fields to read.  Aborting..." << endl;
    invalidate();
    return false;
  }

  if(!validate_fields(true, quiet))
  {
    errors << priority(error)
           << "Invalid list of fields to read.  Aborting." << endl;
    invalidate();
    return false;
  }

  if(!build_data(p))
  {
    errors << priority(error)
           << "Error building data indices.  Aborting." << endl;
    invalidate();
    return false;
  }

  FortranFormatter ffortran;
  ffortran.grab(); // Make sure this doesn't run off
  Tokenizer tokens(&ffortran, infile);
  Tokenizer::token tok;

  ffortran.set_format( format() );

#if 0
  size_t verbose = verbose_output();
#endif

  field_list_type::const_iterator field_info      = my_fields.begin();
  field_list_type::const_iterator field_info_end  = my_fields.end();

  // PASS 2.  Should be seperate so that files of just traits+cov+markers
  //          can be read

  typedef vector<pair<size_t,string> > value_vector;

  // FIXME: This approach to reading markers will certainly have serious
  //        performance problems in the large scale and is here only for the
  //        simplicity.  A hash_map implementation would be much better.  A
  //        separate index vector is the optimal solution in the long term.
  typedef pair<string,string> string_pair;
  typedef map<size_t, string_pair> marker_value_map;

  value_vector trait_values( trait_count() );
  value_vector string_values( string_count() );

#if 0
  if(verbose)
    print_phenotype_header(messages, mped_info, filename);
#endif

  while( infile.good() )
  {
    tokens.read_line();
    ffortran.reset();

    if(tokens.eof())
      break;

    size_t record_starting_line = tokens.line_count();

    string pn;     // Pedigree number
    string id;     // Individual id
    string p1, p2; // Parents
    string sex;    // Sec code

    marker_value_map marker_values;

    for(size_t t = 0; t < trait_count(); ++t)
      trait_values[t].first = (size_t)-1;

    for(size_t s = 0; s < string_count(); ++s)
      string_values[s].first = (size_t)-1;

    size_t tfound = 0;
    size_t sfound = 0;

    for( field_info = my_fields.begin() ;
         field_info != field_info_end ;
         ++field_info)
    {
      switch( field_info->type )
      {
        case   pedigree_id: pn = strip_ws(get_token(tok, &tokens).String());
                            break;
        case individual_id: id = strip_ws(get_token(tok, &tokens).String());
                            break;
        case         trait:
                            if( field_info->index >= mped_info.trait_count() )
                            {
                              get_token(tok, &tokens);
                              break;
                            }
                            trait_values[tfound].first  = field_info->index;
                            trait_values[tfound].second = strip_ws(get_token(tok,&tokens).String());
                            ++tfound;
                            break;

        case  string_field:
                            if( field_info->index >= mped_info.string_count() )
                            {
                              get_token(tok, &tokens);
                              break;
                            }
                            string_values[sfound].first  = field_info->index;
                            string_values[sfound].second = strip_ws(get_token(tok,&tokens).String());
                            ++sfound;
                            break;

        case        allele:
        case        marker:
        {
                            size_t index = field_info->index;
                            if( index >= mped_info.marker_count() )
                            {
                              get_token(tok, &tokens);
                              break;
                            }
                            string_pair &values = marker_values[index];
                            if( !values.first.size())
                              values.first = strip_ws(get_token(tok,&tokens).String());
                            else
                              values.second = strip_ws(get_token(tok,&tokens).String());

                            break;
        }

        default:            get_token(tok, &tokens);
                            break;
      }
    }

#if 0
    if(verbose)
    {
      --verbose;
      print_phenotype(messages, mped_info, pn, id, trait_values, string_values);
    }
#endif

    bool has_pedigree_id_field = false;
  
    for(field_list_type::const_iterator itr = my_fields.begin(); itr != my_fields.end(); ++itr)
      if(itr->type == pedigree_id)
        has_pedigree_id_field = true;

    if(!has_pedigree_id_field)
      pn = "0";

    if(!pn.size() || !id.size() || id == mped_info.individual_missing_code())
      continue;

    RefMultiPedigree::member_pointer mem = p.member_find(pn,id);

    // This should only happen when an error has already been reported
    if( mem == NULL )
      continue;

    RefMultiPedigree::pedinfo_type &info = mem->pedigree()->info();
    int ind_num = mem->index();

    // Make the traits
    for(size_t tt=0; tt < tfound; ++tt)
    {
      size_t t = trait_values[tt].first;
      const string &value = trait_values[tt].second;

      int code = info.set_trait(ind_num, t, value, mped_info.trait_info(t));

      switch(code)
      {
        case 0:                       // trait ok
        case 1:                       // trait ok, but missing
        case 4:                       // trait not set legitamately (expected)
                 break;
        case 3:
                errors << priority(warning) << "[" << record_starting_line
                        << "] Cannot set trait '"
                        << mped_info.trait_info(t).name()
                        << "' for individual '" << id << "' in pedigree '"
                        << pn << "'." << endl;
                  break;
        case 2:                       // bad trait value
                 errors << priority(warning) << "[" << record_starting_line
                        << "] Unrecognized value for trait '"
                        << mped_info.trait_info(t).name()
                        << "' of individual '" << id << "' in pedigree '"
                        << pn << "'.  Found '" << value << "'." << endl;
       default:
                 errors << priority(error) << "[" << record_starting_line
                        << "] Unexpected error (" << code << ") setting trait '"
                        << mped_info.trait_info(t).name()
                        << "' for individual '" << id << "' in pedigree '"
                        << pn << "'." << endl;
                 break;
      }
    }

    for(size_t ss=0; ss < sfound; ++ss)
    {
      size_t s = string_values[ss].first;
      const string &value = string_values[ss].second;

      bool code = info.set_string(ind_num, s, value);

      if(!code)
        errors << priority(warning) << "[" << record_starting_line
               << "] Cannot set string field '"
               << mped_info.string_info(s).name()
               << "' for individual '" << id << "' in pedigree '"
               << pn << "'." << endl;
    }

    marker_value_map::const_iterator mm;
    for(mm = marker_values.begin(); mm != marker_values.end(); ++mm)
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
        case 3:                       // invalid ind. or marker id
                  break;
        case 2:                       // bad marker value
                 errors << priority(warning) << "[" << record_starting_line
                        << "] Unrecognized value for marker '"
                        << mped_info.marker_info(m).name()
                        << "' of individual '" << id << "' in pedigree '"
                        << pn << "': Found '" << values.first << "'";
                  if(values.second.size())
                    errors << ", '" << values.second << "'";
                  errors << ". Marker will be set to missing for this individual." << endl;

                 break;
         case 4:
                 errors << priority(warning) << "[" << record_starting_line
                        << "] Marker '" << mped_info.marker_info(m).name()
                        << "' is sex-dependent, but the individual '" << id 
                        << "' in pedigree '" << pn 
                        << "' has unknown sex.  Marker will be set to missing for this individual." << std::endl;
                 break;
         case 5:
                 errors << priority(warning) << "[" << record_starting_line
                        << "] Phenotype (found '" << values.first << "'";
                 if(values.second.size())
                   errors << ", '" << values.second << "'";
                 errors << " at sex-dependent marker '" << mped_info.marker_info(m).name()
                        << "' is inconsistent with the sex of individual '" << id
                        << "' in pedigree '" << pn 
                        << "'.  Marker will be set to missing for this individual." << std::endl;
                 break;
      }

    }
  }

#if 0
  if(verbose != verbose_output())
    print_phenotype_footer(messages);
#endif

  return true;
}

bool RefFortranPedigreeFile::output(RefMultiPedigree &p, const string &filename,
                                    ostream &messages, bool quiet)
{
  if( !filename.size() || !format().size()) return false;

  std::ofstream outfile(filename.c_str());

  if(!outfile.good())
  {
    errors << priority(critical) << "Unable to open output Family Data file '" << filename
           << "'." << endl;
    return false;
  }

  const RefMPedInfo &mped_info = p.info();

  if(!build_fields(mped_info, quiet))
  {
    errors << priority(error)
           << "Cannot build list of fields to write.  Aborting..." << endl;
    invalidate();
    return false;
  }

  if(!validate_fields(false, quiet))
  {
    errors << priority(error)
           << "Invalid list of fields to write.  Aborting." << endl;
    invalidate();
    return false;
  }

  FortranFormatter ffortran;
  ffortran.grab(); // Make sure this doesn't run off
  OutputTokenizer out(&ffortran, outfile);
  ffortran.set_format( format() );
  Tokenizer::token tok;

  field_list_type::const_iterator field_info;
  field_list_type::const_iterator field_info_end  = my_fields.end();

  RefPedigree::member_type* mid;
  RefPedigree::member_type* pid;
  size_t parent_count;

  vector<size_t> markers_written( mped_info.marker_count(), 0);

  RefMultiPedigree::pedigree_iterator ped;
  for(ped = p.pedigree_begin(); ped != p.pedigree_end(); ++ped)
  {
    const RefPedInfo &ped_info = ped->info();

    for(unsigned int i = 0; i < ped->member_count(); ++i)
    {
      ffortran.reset();
      mid = &ped->member_index(i);
      parent_count = 0;
      markers_written.resize(0);
      markers_written.resize(mped_info.marker_count(), 0);

      for( field_info = my_fields.begin();
           field_info != field_info_end ;
         ++field_info )
      {
        switch( field_info->type )
        {
          case   pedigree_id:
                              tok.first = Tokenizer_Action::String;
                              tok.second = ped->name();
                              out.put_token(tok);
                              break;

          case individual_id:
                              tok.first = Tokenizer_Action::String;
                              tok.second = mid->name();
                              out.put_token(tok);
                              break;

          case     parent_id:
          {
                              tok.first = Tokenizer_Action::String;
                              tok.second = mped_info.individual_missing_code();
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
                                tok.second = pid->name();
                              else
                                parent_count = 2;
                              out.put_token(tok);
                              break;
          }
          case      sex_code:
                              tok.first = Tokenizer_Action::String;
                              if( mid->is_male() )
                                tok.second = mped_info.sex_code_male();
                              else if( mid->is_female() )
                                tok.second = mped_info.sex_code_female();
                              else
                                tok.second = mped_info.sex_code_unknown();
                              out.put_token(tok);
                              break;

          case         trait:
          {
                              if( field_info->index >= mped_info.trait_count() )
                              {
                                tok.first = Tokenizer_Action::String;
                                tok.second = "";
                                out.put_token(tok);
                                break;
                              }

                              if( ped_info.trait_missing(i, field_info->index) )
                              {
                                tok.first = Tokenizer_Action::String;
                                tok.second = mped_info.trait_info(
                                      field_info->index).string_missing_code();
                                out.put_token(tok);
                                break;
                              }

                              tok.first  = Tokenizer_Action::Float;
                              tok.second = ped_info.trait(i, field_info->index);
                              out.put_token(tok);
                              break;
          }

          case    string_field:
          {
                              tok.first  = Tokenizer_Action::String;
                              if( field_info->index >= mped_info.string_count() )
                                tok.second = "";
                              else
                                tok.second = ped_info.get_string(i, field_info->index);
                              out.put_token(tok);
                              break;
          }

          case      allele:
          {
                              tok.first = Tokenizer_Action::String;
                              tok.second = "";
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
                                  tok.second = al1;
                                else
                                  tok.second = al2;
                                markers_written[field_info->index]++;
                              }
                              out.put_token(tok);
                              break;
          }

          case      marker:
          {
                              tok.first = Tokenizer_Action::String;
                              tok.second = "";
                              if( field_info->index >= mped_info.marker_count() )
                              {
                                out.put_token(tok);
                                break;
                              }
                              if(markers_written[field_info->index] == 0)
                              {
                                const RefMarkerInfo& minfo =
                                   mped_info.marker_info(field_info->index);

                                uint pheno = ped_info.phenotype(i, field_info->index);

                                tok.second = minfo.get_phenotype(pheno).name();
                                markers_written[field_info->index] = 2;
                              }
                              out.put_token(tok);
                              break;
          }

          default:
                              tok.first = Tokenizer_Action::String;
                              tok.second = "";
                              out.put_token(tok);
                              break;
        }
      }
      out.write_line();
    }
  }
  outfile.close();
  return true;
}

RefLSFFortranPedigreeFile::
RefLSFFortranPedigreeFile(cerrorstream &err)
{
  set_error_sink(err);
  nondefault_field_order = false;
}

bool
RefLSFFortranPedigreeFile::process_parameters(RefMPedInfo &mped_info, const LSFBase *params)
{
  RefLSFPedigreeFile::process_parameters(mped_info, params);

  if(!format().size())
  {
    errors << priority(error) << "RefPedigreeFile: No FORTRAN Format statement given.  Cannot read data...." << endl;
    invalidate();
  }

  // If they haven't specified their own order, we provide one.
  if(!nondefault_field_order)
  {
    field_list().push_front(parent_id);
    field_list().push_front(parent_id);
    field_list().push_front(sex_code);
    field_list().push_front(individual_id);
    field_list().push_front(pedigree_id);
  }
  return true;
}

bool
RefLSFFortranPedigreeFile::
process_parameter(RefMPedInfo &mped_info, const LSFBase *param)
{
  if(!param)
    return true;

  AttrList::const_iterator a;
  AttrVal v;

  // Allow generic super-class to handle parameters
  RefLSFPedigreeFile::process_parameter(mped_info, param);

  string name = toUpper( param->name() );

  if( name == "STUDY_ID" )
  {
    add_study_id();
  }

  if( name == "TREAT_AS_SIBS" )
  {
    set_treat_as_sibs(true);
  }

  if( name == "PEDIGREE_ID" )
  {
    nondefault_field_order = true;

    add_pedigree_id();
    set_pedigree_id_name("PEDIGREE_ID");
  }

  if( name == "INDIVIDUAL_ID" )
  {
    nondefault_field_order = true;

    add_individual_id();
  }

  if( name == "PARENT_ID" )
  {
    nondefault_field_order = true;

    add_parent_id();
  }

  if( name == "SEX_FIELD" )
  {
    nondefault_field_order = true;

    add_sex_field();
    set_sex_field_name("SEX_CODE");

    RefLSFPedigreeFile::process_sex_code(mped_info, param);
  }

  if( name == "NO_SEX_FIELD" )
  {
    set_no_sex_field(true);
  }

  if( name == "MARKER" || name == "ALLELE" || name == "TRAIT_MARKER" )
  {
    string marker_name;
    if( param->attrs() )
      marker_name = param->attrs()->StringAttr(0);

    if(!marker_name.size())
    {
      errors << priority(warning) << "Marker with no name specified.  Skipping." << endl;
      add_skip_field();
      return true;
    }

    bool skip = skip_markers();
    if( !force_skip_markers() && has_attr(param, "skip") )
    {
      skip = true;
      if( toUpper(attr_value(param, "skip").String()) == "FALSE" )
        skip = false;
    }

    if(skip)
    {
      add_skip_field();
      return true;
    }

    bool dynamic = dynamic_markers();
    if( !force_dynamic_markers() && has_attr(param, "dynamic") )
    {
      dynamic = true;
      if( toUpper(attr_value(param, "dynamic").String()) == "FALSE" )
        dynamic = false;
    }

    size_t m = mped_info.marker_find(marker_name);
    bool has_model = m < mped_info.marker_count();

    if( !has_model && !dynamic )
    {
        errors << priority(warning) << "Marker '" << marker_name
               << "' not found.  Skipping." << endl;
        add_skip_field();
        return true;
    }

    if(!has_model)
    {
      mped_info.add_marker(marker_name);
      m = mped_info.marker_find(marker_name);
      mped_info.marker_info(m).gmodel().set_dynamic_alleles(true);

      char   sep = mped_info.get_pheno_reader_info().get_allele_delimiter();
      mped_info.marker_info(m).gmodel().set_unphased_separator(sep);

      string missing = mped_info.get_pheno_reader_info().get_allele_missing();
      parse_dynamic_markers_missing(mped_info.marker_info(m), missing);
    }
    else if(    (a=param->attrs()->find("equal_allele_freq"))      != param->attrs()->end()
             || (a=param->attrs()->find("complement_allele_freq")) != param->attrs()->end()
             || (a=param->attrs()->find("compl"))                  != param->attrs()->end()
             || (a=param->attrs()->find("minimum_allele_freq"))    != param->attrs()->end()
             || (a=param->attrs()->find("maximum_allele_freq"))    != param->attrs()->end()
             || (a=param->attrs()->find("equal"))        != param->attrs()->end()
             || (a=param->attrs()->find("complement"))   != param->attrs()->end()
             || (a=param->attrs()->find("minimum"))      != param->attrs()->end()
             || (a=param->attrs()->find("maximum"))      != param->attrs()->end() )
    {
      parse_allele_frequency_adjustment(mped_info.marker_info(m), param);
    }

    if(name == "ALLELE")
      add_allele_field( marker_name );
    else
      add_marker_field( marker_name );

    if( param->attrs()->has_attr("x_linked"))
    {
       mped_info.marker_info(m).set_model_type(MLOCUS::X_LINKED);
       set_sex_linked_exist(true);
    }
    else if( param->attrs()->has_attr("y_linked"))
    {
       mped_info.marker_info(m).set_model_type(MLOCUS::Y_LINKED);
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
        update_marker_delimiter_info(mped_info.marker_info(m), local_sep);
    }

    if(    (a=param->attrs()->find("allele_missing")) != param->attrs()->end()
        || (a=param->attrs()->find("missing"))        != param->attrs()->end() )
    {
      string global_missing = mped_info.get_pheno_reader_info().get_allele_missing();
      string local_missing  = strip_ws(a->second.String());

      if( global_missing != local_missing )
        update_marker_missing_info(mped_info.marker_info(m), local_missing);
    }
  }

  if( name == "TRAIT" || name == "COVARIATE" || name == "PHENOTYPE")
  {
    string trait_name;
    if( param->attrs() )
      trait_name = param->attrs()->StringAttr(0);

    if(!trait_name.size())
    {
      errors << priority(warning) << "Trait with no name specified.  Skipping." << endl;
      add_skip_field();
      return true;
    }

    bool skip = skip_traits();
    if( !force_skip_traits() && has_attr(param, "skip") )
    {
      skip = true;
      if( toUpper(attr_value(param, "skip").String()) == "FALSE" )
        skip = false;
    }

    if(skip)
    {
      add_skip_field();
      return true;
    }

    add_trait_field(trait_name);

    size_t t = mped_info.trait_find(trait_name);
    bool has_model = t < mped_info.trait_count();

    if(has_model)
        return true;

    bool trvar = param->attrs()->has_attr("variate") || name == "TRAIT";
    bool cov   = param->attrs()->has_attr("covariate") || name == "COVARIATE";

    RefTraitInfo::trait_use use = RefTraitInfo::unknown_use;

    if(trvar && cov)
    {
      errors << priority(warning) << "Trait '" << trait_name
             << "' has multiple types listed.  Setting to unknown." << endl;
    }
    else if(trvar)
      use = RefTraitInfo::trait_variate;
    else if(cov)
      use = RefTraitInfo::trait_covariate;

    if( param->attrs()->has_attr("binary"))
    {
      string affected   = "1";
      string unaffected = "0";
      string missing    = "";
      double threshold  = std::numeric_limits<double>::quiet_NaN();

      if( (a=param->attrs()->find("affected"))    != param->attrs()->end() )
        affected = strip_ws(a->second.String());
      if( (a=param->attrs()->find("unaffected"))  != param->attrs()->end() )
        unaffected = strip_ws(a->second.String());
      if( (a=param->attrs()->find("missing"))     != param->attrs()->end() )
        missing = strip_ws(a->second.String());
      if( (a=param->attrs()->find("threshold"))   != param->attrs()->end() )
        threshold = a->second.Real();

      if( affected == unaffected )
      {
        errors << priority(error) << "Pedigree File: Affected and unaffected codes "
               << "must be specified and may not be the same." << endl;
        skip = true;
      }
      else if( affected == missing || unaffected == missing )
      {
        errors << priority(error) << "Pedigree File Error: Affected and unaffected codes "
               << "must not be the same as the missing value code." << endl;
        skip = true;
      }
      else
      {
        DEBUG_RPEDFILE(errors << priority(debug) << "Found binary trait = " << v.String() << endl;)
        size_t t = mped_info.add_binary_trait( trait_name, use );

        mped_info.trait_info(t).set_string_missing_code(missing);
        mped_info.trait_info(t).set_numeric_missing_code(str2doub(missing));
        mped_info.trait_info(t).set_string_affected_code(affected);
        mped_info.trait_info(t).set_string_unaffected_code(unaffected);
        mped_info.trait_info(t).set_numeric_affected_code(str2doub(affected));
        mped_info.trait_info(t).set_numeric_unaffected_code(str2doub(unaffected));
        mped_info.trait_info(t).set_threshold(threshold);
      }
    }

    // Categorical!
    else if( param->attrs()->has_attr("categorical"))
    {
      std::string    missing = (a = param->attrs()->find("missing")) == param->attrs()->end() ? ""  : strip_ws(a->second.String()),
                     values  = (a = param->attrs()->find("values"))  == param->attrs()->end() ? ""  : strip_ws(a->second.String());
      RefTraitInfo & tinfo   = mped_info.trait_info(mped_info.add_trait(trait_name, RefTraitInfo::categorical_trait, use));

      tinfo.set_string_missing_code  (missing);
      tinfo.set_numeric_missing_code (str2doub(missing));
      tinfo.set_lockout              (false);

      if(values != "")
      {
        UTIL::StringUtils::splitMultiDelimitedString(values, " \t\n,", tinfo.get_categories());
        tinfo.set_lockout(true);
      }
    }
    else
    {
      string missing;

      if( (a=param->attrs()->find("missing"))  != param->attrs()->end() )
        missing = strip_ws(a->second.String());

      DEBUG_RPEDFILE(errors << priority(debug) << "Found continuous trait = " << trait_name << endl;)
      size_t t = mped_info.add_continuous_trait( trait_name, use );

      mped_info.trait_info(t).set_string_missing_code(missing);
      mped_info.trait_info(t).set_numeric_missing_code(str2doub(missing));
    }

    if(skip)
    {
      DEBUG_RPEDFILE(errors << priority(debug) << "Skipping trait = " << trait_name << endl;)
      mped_info.add_trait( trait_name, RefTraitInfo::invalid_trait );
    }
  }

  if( name == "STRING" )
  {
    string string_name;
    if( param->attrs() )
      string_name = param->attrs()->StringAttr(0);

    if(!string_name.size())
    {
      errors << priority(warning) << "String field with no name specified.  Skipping." << endl;
      add_skip_field();
      return true;
    }

    bool skip = false;
    if( has_attr(param, "skip") )
    {
      skip = true;
      if( toUpper(attr_value(param, "skip").String()) == "FALSE" )
        skip = false;
    }

    if(skip)
    {
      add_skip_field();
      return true;
    }

    size_t s = mped_info.string_find(string_name);
    bool has_model = s < mped_info.string_count();

    if(!has_model)
      mped_info.add_string_field(string_name);

    add_string_field(string_name);
  }
  
  return true;
}

} // End namespace RPED
} // End namespace SAGE
