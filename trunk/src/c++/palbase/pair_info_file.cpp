//****************************************************************************
//* File:      pairinfo_file.cpp                                             *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0 Initial implementation                yjs Jan. 02 *
//*                                                                          *
//* Notes:     This file implements class for I/O of pair information file.  *
//*                                                                          *
//* Copyright (c) 2002 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include <fstream>
#include "palbase/pair_info_file.h"

namespace SAGE    {
namespace PALBASE {

RefPairInfoFile::RefPairInfoFile(cerrorstream &err)
               : errors(err)
{
  set_verbose_output(0);
  reset_counts();
  validate();
}

RefPairInfoFile::~RefPairInfoFile() 
{
  errors.flush();
}

void RefPairInfoFile::reset_counts()
{
  my_skip_count = 0;
  my_pedigree_id_count = 0;
  my_pair_id_count = 0;
  my_pair_covariate_count = 0;
  my_pair_weight_count = 0;
  my_invalid_covariate_count = 0;
  my_invalid_weight_count = 0;
}

void RefPairInfoFile::
print_pair_info_header(ostream &messages, const relative_pairs &rpairs, const string &filename) const
{
  messages << endl << "Phenotypes for the first " << verbose_output()
                   << " pairs read from file: " << filename
                   << endl << endl;

  messages << "     PED ID        IND. ID       IND. ID     ";

  field_list_type::const_iterator field_info      = my_fields.begin();
  field_list_type::const_iterator field_info_end  = my_fields.end();

  for( ; field_info != field_info_end ; ++field_info)
  {
    if(    field_info->type == pair_weight
        && field_info->index < rpairs.pair_weight_count() )
      messages << "  " << setw(20) << rpairs.pair_weight_info(field_info->index).get_name();

    if(    field_info->type == pair_covariate
        && field_info->index < rpairs.pair_covariate_count() )
      messages << "  " << setw(20) << rpairs.pair_covariate_info(field_info->index).get_name();
  }

  messages << endl << "     ------------  ------------  ------------";

  for( field_info = my_fields.begin(); 
       field_info != field_info_end ;
       ++field_info)
  {
    if(    field_info->type == pair_weight
        && field_info->index < rpairs.pair_weight_count())
      messages << "  --------------------";

    if(    field_info->type == pair_covariate
        && field_info->index < rpairs.pair_covariate_count())
      messages << "  --------------------";
  }

  messages << endl;
}

void RefPairInfoFile::
print_pair_info(ostream &messages, const relative_pairs &rpairs,
                const string& pn, const string& id1, const string& id2,
                const vector<pair<size_t,string> >& pair_weight_values,
                const vector<pair<size_t,string> >& pair_covariate_values) const
{
  messages << "     " << setw(12) << pn << "  " 
           << setw(12) << id1 << "  " << setw(12) << id2;

  for( size_t i = 0; i < pair_weight_values.size(); ++i )
    if( pair_weight_values[i].first < rpairs.pair_weight_info_count())
      messages << "  " << setw(20) << pair_weight_values[i].second;

  for( size_t i = 0; i < pair_covariate_values.size(); ++i )
    if( pair_covariate_values[i].first < rpairs.pair_covariate_info_count())
      messages << "  " << setw(20) << pair_covariate_values[i].second;
  messages << endl;
}

void RefPairInfoFile::print_pair_info_footer(ostream &messages) const
{
  messages << endl;
}

bool RefPairInfoFile::validate_fields(bool quiet)
{
  // NOTE: We treat invalid markers and traits as non-fatal states

  reset_counts();

  typedef std::map<string, size_t> name_count_map;

  name_count_map pair_covariate_map;
  name_count_map pair_weight_map;

  field_list_type::iterator i;
  for( i = my_fields.begin(); i != my_fields.end(); ++i )
    switch( i->type )
    {
      case skip          : ++my_skip_count; break;
      case pedigree_id   : ++my_pedigree_id_count; break;
      case pair_id       : ++my_pair_id_count; break;
      case pair_covariate: pair_covariate_map[ toUpper(i->name) ]++; break;
      case pair_weight   : pair_weight_map[ toUpper(i->name) ]++; break;
    }

  bool valid_fields = true;

  if(my_pedigree_id_count == 0)
  {
    valid_fields = false;
    if(!quiet)
      errors << priority(error)
             << "Pedigree ID field not found." << endl;
  }

  if(my_pedigree_id_count > 1)
  {
    valid_fields = false;
    if(!quiet)
      errors << priority(error)
             << "Too many Pedigree ID fields found." << endl;
  }

  if(my_pair_id_count == 0)
  {
    valid_fields = false;
    if(!quiet)
      errors << priority(error)
             << "Pair ID field not found." << endl;
  }

  if(my_pair_id_count > 2)
  {
    valid_fields = false;
    if(!quiet)
      errors << priority(error)
             << "Too many Pair ID fields found." << endl;
  }

  for(i = my_fields.begin(); i != my_fields.end(); ++i)
  {
    if(i->type == pair_covariate)
    {
      size_t &count = pair_covariate_map[ toUpper(i->name) ];

      if(count == 1)
        ++my_pair_covariate_count;
      else if( count != (size_t)-1 )  // Invalidate new bad covariate
      {
        if(!quiet)
          errors << priority(warning) 
                 << "Covariate assigned to more than one field.  Skipping covariate '" 
                 << i->name << "'." << endl;

        // Count covariate as bad and mark count to avoid further warnings
        ++my_invalid_covariate_count;
        count = (size_t)-1;
        i->index = (size_t)-1;
      }
      else // invalidate known bad trait
        i->index = (size_t)-1;
    }
    else if(i->type == pair_weight)
    {
      size_t &count = pair_weight_map[ toUpper(i->name) ];

      if(count == 1)
        ++my_pair_weight_count;
      else if( count != (size_t)-1 )  // Invalidate new bad weight
      {
        if(!quiet)
          errors << priority(warning) 
                 << "Weight assigned to more than one field.  Skipping weight '" 
                 << i->name << "'." << endl;

        // Count weight as bad and mark count to avoid further warnings
        ++my_invalid_weight_count;
        count = (size_t)-1;
        i->index = (size_t)-1;
      }
      else // invalidate known bad weight
        i->index = (size_t)-1;
    }
  }

  if(my_pair_weight_count > 1)
  {
    valid_fields = false;
    if(!quiet)
    {
      cout << endl;
      errors << priority(error)
             << "Too many pair_weight fields found." << endl;
    }
  }

  return valid_fields;
}

void RefPairInfoFile::build_pair_info(relative_pairs &rpairs)
{
  // Create indices
  if( rpairs.pair_covariate_info_count() )
    rpairs.resize_pair_covariates( rpairs.pair_covariate_info_count() );

  if( rpairs.pair_weight_info_count() )
    rpairs.resize_pair_weights( rpairs.pair_weight_info_count() );
}

size_t RefPairInfoFile::pair_find(relative_pairs &rpairs, const string &pn,
                                  const string &id1, const string &id2,
                                  size_t line)
{
  if( !pn.size() )
  {
    errors << priority(warning) << "[" << line
           << "] record is missing pedigree ID.  Skipping..."
           << endl;
    return (size_t)-1;
  }

  if( !id1.size() )
  {
    errors << priority(warning) << "[" << line
           << "] record is missing the first pair ID.  Skipping..."
           << endl;
    return (size_t)-1;
  }

  if( !id2.size() )
  {
    errors << priority(warning) << "[" << line
           << "] record is missing the second pair ID.  Skipping..."
           << endl;
    return (size_t)-1;
  }

  FPED::PedigreePointer p = rpairs.get_fped().pedigree_find(pn);

  if( !p )
  {
    errors << priority(warning) << "[" << line
           << "] Bad pedigree ID."
           << endl;
    return (size_t)-1;
  }
 
  FPED::MemberPointer m1 = p->member_find(id1);
  FPED::MemberPointer m2 = p->member_find(id2);

  if( !m1 )
  {
    errors << priority(warning) << "[" << line
           << "] Bad First Pair ID."
           << endl;
    return (size_t)-1;
  }
  
  if( !m2 )
  {
    errors << priority(warning) << "[" << line
           << "] Bad Second Pair ID."
           << endl;
    return (size_t)-1;
  }

  size_t pair_num = rpairs.find_pair(m1, m2);
  if( pair_num <= rpairs.pair_count() )
    return pair_num;

  return rpairs.find_pair(m2, m1);
}

//
//------------------------------------------------------------------
//

RefDelimitedPairInfoFile::RefDelimitedPairInfoFile(cerrorstream &err)
                        : RefPairInfoFile(err)
{
  set_whitespace(" \r\n");
  set_delimiters("\t, ");
  set_skip_consecutive_delimiters(true);
  set_skip_leading_delimiters(true);
  set_skip_trailing_delimiters(true);
  reset_counts();
}

RefDelimitedPairInfoFile::~RefDelimitedPairInfoFile() 
{
}

bool RefDelimitedPairInfoFile::
build_fields(string_tokenizer &header, relative_pairs &rpairs,
             bool quiet)
{
#if 0
  cout << endl << "before build_fields().." << endl;
  cout << "pair_c_info_count = " << rpairs.pair_covariate_info_count() << endl;
  cout << "pair_w_info_count = " << rpairs.pair_weight_info_count() << endl;
  cout << "field_list size = " << my_fields.size() << endl;
  cout << "field_map  size = " << my_field_map.size() << endl;
#endif

  // Set the field_map for the first three id fields.
  //
  string_tokenizer::iterator i;
  size_t field_index = 0;
  for( i = header.begin(); i != header.end() && field_index < 3; ++i, ++field_index )
  {
    if( !i->size() )
       continue;

    field_map_type::const_iterator f = my_field_map.find( toUpper(*i) );

    if( f == my_field_map.end() && field_index == 0 )
    {
      string field_name = *i;
      my_field_map[ toUpper(field_name) ] = field(pedigree_id, field_name);
    }

    if( f == my_field_map.end() && field_index == 1 )
    {
      string field_name = *i;
      my_field_map[ toUpper(field_name) ] = field(pair_id, field_name, field_name);
    }

    if( f == my_field_map.end() && field_index == 2 )
    {
      string field_name = *i;
      my_field_map[ toUpper(field_name) ] = field(pair_id, field_name, field_name);
    }
  }

  // If only ids in my_fields list, then the rest of fields becomes pair_covariates.
  //
  if(    !rpairs.pair_covariate_info_count()
      && !rpairs.pair_weight_info_count()
      && i != header.end() )
  {
#if 0
  cout << endl << "middle, no user-specified.." << endl;
  cout << "pair_c_info_count = " << rpairs.pair_covariate_info_count() << endl;
  cout << "pair_w_info_count = " << rpairs.pair_weight_info_count() << endl;
  cout << "field_list size = " << my_fields.size() << endl;
  cout << "field_map  size = " << my_field_map.size() << endl;
#endif

    for( ; i != header.end(); ++i )
    {
      if( !i->size() )
         continue;

      field_map_type::const_iterator f = my_field_map.find( toUpper(*i) );

      if( f == my_field_map.end() )
      {
        string field_name = *i;

        // Add into field_map.
        add_pair_covariate_field(field_name, field_name);

        if(!field_name.size())
        {
          errors << priority(warning) << "Pair Covariate with no name specified.  Skipping." << endl;
          continue;
        }

        // Add into rpairs.pair_covariate_info.
        size_t c = rpairs.pair_covariate_find( field_name );

        if( c < rpairs.pair_covariate_info_count() )  // Should never happen.
        {
          errors << priority(error) << "Should never happen!!  Skipping." << endl;
          continue;
        }

        string missing;
        pair_pheno_info::info_use usage = pair_pheno_info::mean;

        c = rpairs.add_pair_covariate(field_name, usage);

        // Update the my_fields.
        field_list().push_back( field(pair_covariate, field_name, field_name, c) );
      }
    }

    return validate_fields(quiet);
  }

  // If pair_covariate or pair_weight specified by user, then use only those fields.
  //
  assert( my_fields.size() == 3 );

#if 0
  cout << endl << "middle, with user-specified.." << endl;
  cout << "pair_c_info_count = " << rpairs.pair_covariate_info_count() << endl;
  cout << "pair_w_info_count = " << rpairs.pair_weight_info_count() << endl;
  cout << "field_list size = " << my_fields.size() << endl;
  cout << "field_map  size = " << my_field_map.size() << endl;
#endif

  for( ; i != header.end(); ++i )
  {
    my_fields.push_back( field() );

    if( !i->size() )
       continue;

    field_map_type::const_iterator f = my_field_map.find( toUpper(*i) );

    if( f == my_field_map.end() || f->second.type == skip)
    {
      continue;
    }

    field &current_field = my_fields.back();

    // This makes sense -- current field is a reference
    current_field = f->second;

    if( current_field.type == pair_covariate )
    {
      size_t c = rpairs.pair_covariate_find(current_field.name);

      if( c < rpairs.pair_covariate_info_count() )
        current_field.index = c;
      else
      {
        current_field.index = (size_t)-1;

        // Can never happen (yet) since we always add the pair_covariate to rpairs
        // in the LSF parser.  But we have to look out for the future.
        errors << priority(warning) << "Invalid pair_covariate '" << current_field.name 
               << "'.  Skipping..." << endl;
      }
    } 
    else if( current_field.type == pair_weight )
    {
      size_t w = rpairs.pair_weight_find(current_field.name);

      if( w < rpairs.pair_weight_info_count() )
        current_field.index = w;
      else
      {
        current_field.index = (size_t)-1;

        // Can never happen (yet) since we always add the pair_weight to rpairs
        // in the LSF parser.  But we have to look out for the future.
        errors << priority(warning) << "Invalid pair_weight '" << current_field.name 
               << "'.  Skipping..." << endl;
      }
    }
  }

#if 0
  cout << endl << "after build_fields().." << endl;
  cout << "pair_c_info_count = " << rpairs.pair_covariate_info_count() << endl;
  cout << "pair_w_info_count = " << rpairs.pair_weight_info_count() << endl;
  cout << "field_list size = " << my_fields.size() << endl;
  cout << "field_map  size = " << my_field_map.size() << endl;
#endif

  return validate_fields(quiet);
}

bool RefDelimitedPairInfoFile::input(relative_pairs &rpairs, const string &filename, 
                                     ostream &messages)
{
  if( !filename.size() )
  {
    errors << priority(critical) << "No Pair Information file specified." << endl;
    return false;
  }

  std::ifstream infile( filename.c_str());

  if(!infile.good())
  {
    errors << priority(critical) << "Unable to open Pair Information file '" << filename
                    << "'. Please check your file." << endl;
    return false;
  }

  string_tokenizer tokenizer( format() );
  tokenizer.set_whitespace( whitespace() );
  tokenizer.set_delimiters( delimiters() );
  tokenizer.set_skip_consecutive_delimiters( skip_consecutive_delimiters() );
  tokenizer.set_skip_leading_delimiters( skip_leading_delimiters() );
  tokenizer.set_skip_trailing_delimiters( skip_trailing_delimiters() );

  string line;   // Current line

  // Read header iff there is no format

  if( !format().size() )
  {
    getline(infile, line);
    set_format(line);
    tokenizer.set_str(line);
  }

  // Build up the my_fileds list to read the record.
  //
  if(!build_fields(tokenizer, rpairs, false))
  {
    invalidate();
    return false;
  }

  build_pair_info(rpairs);

  my_format = line;

  size_t verbose = verbose_output();

  if(verbose)
    print_pair_info_header(messages, rpairs, filename);

  for( size_t count = 1; infile.good(); ++count )
  {
    getline(infile, line);

    if(!line.size())
      continue;

    string pn;     // Pedigree id
    string id1;    // Pair id1
    string id2;    // Pair id2

    vector<pair<size_t,string> > pair_covariate_values( pair_covariate_count() );
    vector<pair<size_t,string> > pair_weight_values( pair_weight_count() );

    for( size_t c = 0; c < pair_covariate_count(); ++c )
      pair_covariate_values[c].first = (size_t)-1;

    for( size_t w = 0; w < pair_weight_count(); ++w )
      pair_weight_values[w].first = (size_t)-1;

    tokenizer.set_str(line);

    field_list_type::const_iterator field_info      = my_fields.begin();
    field_list_type::const_iterator field_info_end  = my_fields.end();
    string_tokenizer::const_iterator field     = tokenizer.begin();
    string_tokenizer::const_iterator field_end = tokenizer.end();

    size_t cfound = 0;
    size_t wfound = 0;

    for( ; 
         field_info != field_info_end && field != field_end ;
         ++field_info, ++field)
    {
      string field_value = *field;

#if 0
cout << endl
     << "field_info->index = " << field_info->index
     << ", field_value = " << field_value;
#endif

      switch( field_info->type )
      {
        case pedigree_id :
//cout << ", pedigree_id" << endl;
             pn = field_value;
             break;
        case pair_id :
//cout << ", pair_id" << endl;
             if(!id1.size())
               id1 = field_value;
             else
               id2 = field_value;
             break;
        case pair_covariate : 
//cout << ", pair_covariate" << endl;
             if( field_info->index >= rpairs.pair_covariate_info_count() )
               break;
             pair_covariate_values[cfound].first  = field_info->index;
             pair_covariate_values[cfound].second = field_value;
             ++cfound;
             break;
        case pair_weight :
//cout << ", pair_weight" << endl;
             if( field_info->index >= rpairs.pair_weight_info_count() )
               break;
             pair_weight_values[wfound].first  = field_info->index;
             pair_weight_values[wfound].second = field_value;
             ++wfound;
             break;
        case skip :// break;
//cout << ", skip" << endl; break;
        default   : break;
      }      
    }

    if( !pn.size() && !id1.size() && !id2.size() )
      continue;

    if(verbose)
    {
      --verbose;
      print_pair_info(messages, rpairs, pn, id1, id2,
                      pair_weight_values, pair_covariate_values);
    }

    // instead add, check.
    size_t pair_num = pair_find(rpairs, pn, id1, id2, count);
    if( pair_num > rpairs.pair_count() )
    { 
      errors << priority(warning) 
             << "Cannot find pair: " << pn << " (" << id1 << ", " << id2 << ").  Skipping..." << endl;
      continue;
    }

    for( size_t cc = 0; cc < cfound; ++cc )
    {
      size_t c = pair_covariate_values[cc].first;
      const string &value = pair_covariate_values[cc].second;
      int code = rpairs.set_pair_covariate(pair_num, c, value, rpairs.pair_covariate_info(c));

      switch(code)
      {
        case 0:                       // pair_covariate ok
        case 1:                       // pair_covariate ok, but missing
                 break;
        case 3:                       // invalid ind. or pair_covariate id
                 errors << priority(warning) << "[" << count 
                        << "] Cannot set pair_covariate '"  
                        << rpairs.pair_covariate_info(c).get_name() 
                        << "' for pair '" << id1 << ", " << id2 << "' in pedigree '"
                        << pn << "'." << endl;
                 break;
        case 2:                       // bad pair_covariate value
                 errors << priority(warning) << "[" << count
                        << "] Unrecognized value for pair_covariate '"
                        << rpairs.pair_covariate_info(c).get_name()
                        << "' for pair '" << id1 << ", " << id2 << "' in pedigree '"
                        << pn << "'.  Found '" << value << "'." << endl;

                 break;
       default:
                 errors << priority(error) << "[" << count 
                        << "] Unexpected error setting pair_covariate '"  
                        << rpairs.pair_covariate_info(c).get_name() 
                        << "' for pair '" << id1 << ", " << id2 << "' in pedigree '"
                        << pn << "'." << endl;
                 break;
      }
    }

    for( size_t ww = 0; ww < wfound; ++ww )
    {
      size_t w = pair_weight_values[ww].first;
      const string &value = pair_weight_values[ww].second;
      int code = rpairs.set_pair_weight(pair_num, w, value, rpairs.pair_weight_info(w));

      switch(code)
      {
        case 0:                       // pair_weight ok
        case 1:                       // pair_weight ok, but missing
                 break;
        case 3:                       // invalid ind. or pair_weight id
                 errors << priority(warning) << "[" << count 
                        << "] Cannot set pair_weight '"  
                        << rpairs.pair_weight_info(w).get_name() 
                        << "' for pair '" << id1 << ", " << id2 << "' in pedigree '"
                        << pn << "'." << endl;
                 break;
        case 2:                       // bad pair_weight value
                 errors << priority(warning) << "[" << count
                        << "] Unrecognized value for pair_weight '"
                        << rpairs.pair_weight_info(w).get_name()
                        << "' for pair '" << id1 << ", " << id2 << "' in pedigree '"
                        << pn << "'.  Found '" << value << "'." << endl;

                 break;
       default:
                 errors << priority(error) << "[" << count 
                        << "] Unexpected error setting pair_weight '"  
                        << rpairs.pair_weight_info(w).get_name() 
                        << "' for pair '" << id1 << ", " << id2 << "' in pedigree '"
                        << pn << "'." << endl;
                 break;
      }
    }
  }


  if(verbose != verbose_output() )
    print_pair_info_footer(messages);

  return true;
}

//
//------------------------------------------------------------------
//

RefLSFDelimitedPairInfoFile::
RefLSFDelimitedPairInfoFile(const LSFBase *params,
                            relative_pairs &rpairs, cerrorstream &err )
         : RefDelimitedPairInfoFile(err)
{
  if(!params)
  {
    errors << priority(critical) << "RefPairInfoFile Error: No parameters specified." << endl;
    invalidate();
    return;
  }

  // Set defaults
  set_verbose_output(10);
  set_skip_trailing_delimiters(true);
  set_skip_leading_delimiters(true);
  set_skip_consecutive_delimiters(true);

  // Set the order of the first three fields.
  add_pedigree_id_field();
  add_pair_id_field();
  add_pair_id_field();

  // Read the covariate & weight fields if there is any exist.
  if( !params->List() )
    return;

  LSFList::const_iterator i;
  AttrList::const_iterator a;
  AttrVal v;
  for( i = params->List()->begin(); i != params->List()->end(); ++i )
  {
    if( !*i ) continue;
    string name = toUpper( (*i)->name() );

    if( name == "PAIR_COVARIATE" )
    {
      string field_name;

      if( (*i)->attrs() )
        field_name = (*i)->attrs()->StringAttr(0);
      else
        continue;  // add error message?

      string pair_covariate_name;

      pair_covariate_name = (*i)->attrs()->StringAttr("name");

      if(!pair_covariate_name.size())
        pair_covariate_name = field_name;

      add_pair_covariate_field( field_name, pair_covariate_name );

      if(!field_name.size())
      {
        errors << priority(warning) << "Pair Covariate with no name specified.  Skipping." << endl;
        continue;
      }

      size_t c = rpairs.pair_covariate_find( pair_covariate_name );

      if( c < rpairs.pair_covariate_info_count() )
        continue;

      string missing;
      pair_pheno_info::info_use usage = pair_pheno_info::mean;
      double value = std::numeric_limits<double>::quiet_NaN();

      if( (a=(*i)->attrs()->find("missing"))  != (*i)->attrs()->end() )
        missing = strip_ws(a->second.String());

      if( (a=(*i)->attrs()->find("mean"))  != (*i)->attrs()->end() )
        value = (*i)->attrs()->FloatAttr("mean");
      else if( (a=(*i)->attrs()->find("minimum"))  != (*i)->attrs()->end() )
        usage = pair_pheno_info::minimum;

      c = rpairs.add_pair_covariate(pair_covariate_name, usage, value);

      rpairs.pair_covariate_info(c).set_string_missing_code(missing);
      rpairs.pair_covariate_info(c).set_numeric_missing_code(str2doub(missing));
    }

    if( name == "PAIR_WEIGHT" )
    {
      string field_name;

      if( (*i)->attrs() )
        field_name = (*i)->attrs()->StringAttr(0);
      else
        continue;  // add error message?

      string pair_weight_name;

      pair_weight_name = (*i)->attrs()->StringAttr("name");

      if(!pair_weight_name.size())
        pair_weight_name = field_name;

      add_pair_weight_field( field_name, pair_weight_name );

      if(!field_name.size())
      {
        errors << priority(warning) << "Pair Weight with no name specified.  Skipping." << endl;
        continue;
      }

      size_t w = rpairs.pair_weight_find( pair_weight_name );

      if( w < rpairs.pair_weight_info_count() )
        continue;

      string missing;

      if( (a=(*i)->attrs()->find("missing"))  != (*i)->attrs()->end() )
        missing = strip_ws(a->second.String());

      w = rpairs.add_pair_weight(pair_weight_name);

      rpairs.pair_weight_info(w).set_string_missing_code(missing);
      rpairs.pair_weight_info(w).set_numeric_missing_code(str2doub(missing));
    }
  }
}

void
read_pair_info_file(const vector<LSF_ptr<LSFBase> > f_list, relative_pairs &p,
                    cerrorstream &errors, ostream& info_file)
{
  if( !f_list.size() )
      return;

  cout << "Reading Pair Information File............." << endl;

  p.invalidate_pair_info();

  boost::shared_ptr<PALBASE::RefPairInfoFile> pair_reader;

  for( size_t i = 0; i < f_list.size(); ++i )
  {
    LSF_ptr<LSFBase> f = f_list[i];

    if( !f || !f->name().size() )
      return;

    if( toUpper(f->name()) == "PAIR_INFO_FILE" )
    {
      AttrVal a=attr_value(f,0);

      if( a.has_value() )
      {
        string filename = strip_ws(a.String());

        pair_reader = boost::shared_ptr<PALBASE::RefPairInfoFile>(
             new PALBASE::RefLSFDelimitedPairInfoFile(f, p, errors) );

        cout << "              from " << filename << flush;

        if( !pair_reader->input(p, filename, info_file) )
        {
          errors << priority(error)
                 << "Error reading pair information from file '" << filename
                 << "'.  Skipping..."
                 << endl;

          p.invalidate_pair_info();
          continue;
        }

        char old_fill = cout.fill('.');
        cout << setw(23-filename.size()) << "" << "done." << endl;
        cout.fill(old_fill);
      }
      else
      {
        errors << priority(error)
               << "No Pair Information File name specified.  Skipping..."
               << endl;
      }
    }
  }

  return;
}

} // end of namespace PALBASE
} // end of namespace SAGE

