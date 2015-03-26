//
//  IBD sharing file -- I/O of allele sharing IBD for General
//                      data structures
//
//  Author: Kevin Jacobs (jacobs@darwin.cwru.edu)
//
//  History:   0.1  kbj Initial implementation              Aug 25 98
//  Modified by:    Alexandra Borosova, Nov. 2003                    
//  Modification 1: "Zero marker information"
//  Last change: Alexandra Borosova, 02/05/04
//
//  Copyright (c) 1998  R.C. Elston
// -------------------------------------------------------------------------
// Introduction:
//   This class implements streaming of pairs to and from an ASCII data
//   file format.  This file stores only pairs and marker allele sharing
//   at a set of loci.
//

#include "ibd/ibdfile.h"

using namespace std;

namespace SAGE {

void
print_ibd_option(const ibd_option_type& opt)
{
  cout << endl << "#" << endl;
  cout << "ANALYSIS\n#------" << endl;
  cout << "title           = " << opt.title << endl;
  cout << "region          = " << opt.region;
  if( opt.x_linked )
    cout << ", x_linked";
  cout << endl;
  if( opt.exact )
    cout << "max_pedigree    = " << opt.max_pedigree << endl;
  cout << "scan_type       = " << opt.scan_type << endl;
  cout << "allow_loops     = " << opt.allow_loops << endl;
  cout << "ibd_mode        = " << opt.ibd_mode;
  if( opt.exact )
    cout << ", exact";
  cout << endl;
  cout << "split_pedigrees = " << opt.split_pedigrees << endl;
  cout << "use_simulation  = " << opt.use_simulation << endl;
  cout << endl
       << "#" << endl;
}

//
//----------------------------------------------------------------
//

bool
RefIBDReadFileStdIO::input(const std::string &fname, IBD *ibd)
{
  reset(fname);

  FILE *file = fopen(fname.c_str(), "r");

  if(!file || ferror(file) || feof(file) )
  {
    read_error("Cannot open filename '" + fname + "'.");
    return false;
  }
  bool ret = do_input(file, ibd);

  fclose(file);

  return ret;
}

bool
RefIBDReadFileStdIO::input_ibd_state(const std::string &fname, IBD *ibd)
{
  reset(fname);

  FILE *file = fopen(fname.c_str(), "r");

  if(!file || ferror(file) || feof(file) )
  {
    read_error("Cannot open filename '" + fname + "'.");
    return false;
  }
  bool ret = do_input_ibd_state(file, ibd);

  fclose(file);

  return ret;
}

bool
RefIBDReadFileStdIO::do_input(FILE *file, IBD *ibd)
{
  if(!ibd)
  {
    read_error("Invalid IBD storage");
    return false;
  }

  // Read Pair Data Header Information.
  if(!file || feof(file) || ferror(file))
  {
    read_error("Cannot read IBD file header information");
    return false;
  }

  // Reset position information
  reset();

  // Get the signature line and check it for correctness.
  static const string signature = "IBD File";

  string line;
  input_next_line(file,line);

  if( strip_ws(line).substr(0, signature.size()) != signature )
  {
    read_error("Invalid IBD file signature found: " + line);
    return false;
  }
  ibd_option_type ibd_option;

  vector<size_t>  markers;
  read_header(file, ibd, line, ibd_option, markers);

  string ped_name;
  string ind1_name;
  string ind2_name;

  string_tokenizer linebuf;

  while( !feof(file) && !ferror(file) )         // Read in until out of records
  {
    bool quiet = false;

    input_next_line(file,line);

    // Skip blank lines
    if( !line.size() )
      continue;

    linebuf.set_str(line);
    linebuf.set_whitespace("");
    linebuf.set_delimiters(",");

    string_tokenizer::iterator i     = linebuf.begin();
    string_tokenizer::iterator end   = linebuf.end();

    if( !read_pair_ids(i, end, ped_name, ind1_name, ind2_name) )
    {
      read_error("Pair skipped due to error...");
      continue;
    }

    linebuf.set_delimiters(", ");
    linebuf.set_skip_consecutive_delimiters(true);
    linebuf.set_skip_leading_delimiters(true);
    linebuf.set_skip_trailing_delimiters(true);

    if( i == end || ++i == end )
    {
      read_error("Pair contains no data");
      continue;
    }

    IBD::error_t err;
    id_pair p;

    if( ibd->has_pedigree(ped_name) )
    {
      p = ibd->get_pair(ped_name, ind1_name, ind2_name, err);

      if( err != IBD::no_error )
      {
        switch( err )
        {
          case IBD::bad_pedigree:
            read_error("cannot find pedigree '" + ped_name + "'");
            break;
          case IBD::bad_ind1:
            read_error("member id '" + ind1_name + "'"
                     + " in pedigree '" + ped_name + "' was not found");
            break;
          case IBD::bad_ind2:
            read_error("member id '" + ind2_name + "'"
                     + " in pedigree '" + ped_name + "' was not found");
            break;
          case IBD::bad_ind_both:
            read_error("member ids ('" + ind1_name + "','" + ind2_name
                      + "') in pedigree '" + ped_name + "' was not found");
            break;
          default:
            read_error("unknown error finding member ids ('" + ind1_name
                     + "','" + ind2_name + "') in pedigree '" + ped_name + "'");
            break;
        }

        continue;
      }

      if( require_pedigree() )
      {
        if( p.first == NULL && p.second == NULL )
        {
          read_error("Internal error.  Members ('" + ind1_name + "','"
              + ind2_name + "') in pedigree '" + ped_name + "' are invalid");
          continue;
        }

        if( p.first == NULL )
        {
           read_error("Internal error.  Member '" + ind1_name + "'"
                    + " in pedigree '" + ped_name + "' is invalid");
           continue;
        }

        if( p.second == NULL )
        {
           read_error("Internal error.  Member '" + ind2_name + "'"
                    + " in pedigree '" + ped_name + "' is invalid");
           continue;
        }
      }
    }
    else
    {
      read_error("Pair skipped: Invalid pedigree name '" + ped_name + "'.");
      continue;
    }

    // Skip pairs that will not be used for analysis
    if( skip_unused_pairs() && !ibd->use_pair(p.first, p.second) )
    {
      //read_error("Pair skipped since it won't be used.");
      continue;
    }

    size_t pair_num = ibd->add_pair( p.first, p.second );

    if( pair_num == (size_t)-1 )
    {
      if( warn_invalid_pairs() )
      {
        read_error("Warning: Members ('" + ind1_name + "','"
                 + ind2_name + "') in pedigree '" + ped_name
                 + "' cannot be used.");
      }
      continue;
    }

    size_t m;
    for( m = 0; i != end && m < markers.size(); ++m )
    {
      string field0 = *i++;
      string field1 = "-----------------";
      string field2;

      if( !ibd_option.old_ibd_format )
        field1 = *i++;

      if( i != end )
        field2 = *i++;

      if( !field0.size() || !field2.size() )
      {
        read_error("Invalid data format");
        break;
      }

      double f0 = QNAN;
      double f1 = QNAN;
      double f2 = QNAN;

      if( field0[0] != '-' || field2[0] != '-' )
      {
        f0 = str2doub(field0);
        f1 = str2doub(field1);
        f2 = str2doub(field2);
      }

      if( !ibd->set_ibd( pair_num, markers[m], f0, f1, f2 ) && !quiet )
      {
        read_error("Error setting IBD values for individuals ("
                    + ind1_name + "," + ind2_name + ") in pedigree '"
                    + ped_name + "'.");
        quiet = true;
      }
    }

    if( m < markers.size() )
    {
      read_error("Invalid IBD format for individuals ("
                      + ind1_name + "," + ind2_name + ") in pedigree '"
                      + ped_name + "'. Pair skipped");
      ibd->invalidate_pair( pair_num );
      continue;
    }
  }

  return true;
}

bool
RefIBDReadFileStdIO::do_input_ibd_state(FILE *file, IBD *ibd)
{
  if(!ibd)
  {
    read_error("Invalid IBD storage");
    return false;
  }

  // Read Pair Data Header Information.
  if(!file || feof(file) || ferror(file))
  {
    read_error("Cannot read IBD file header information");
    return false;
  }

  // Reset position information
  reset();

  // Get the signature line and check it for correctness.
  static const string signature = "IBD STATE File";

  string line;
  input_next_line(file,line);

  if( strip_ws(line).substr(0, signature.size()) != signature )
  {
    read_error("Invalid IBD state file signature found: " + line);
    return false;
  }

  ibd_option_type ibd_option;
  vector<size_t>  markers;
  read_header(file, ibd, line, ibd_option, markers);

  string_tokenizer linebuf;

  while( !feof(file) && !ferror(file) ) // Read in until out of records
  {
    markers_ibd_state  sped_ibd_state;

    vector< pair<string, string> > ind_names;
    vector< id_pair >              pair_ids;
    string                         ped_name;

    for( size_t m = 0; m < markers.size(); ++m )
    {
      input_next_line(file,line);
#if 0
      cout << "m = " << m << endl;
#endif
      linebuf.set_str(line);
      linebuf.set_whitespace("");
      linebuf.set_delimiters(" ");
      linebuf.set_skip_consecutive_delimiters(true);
      linebuf.set_skip_leading_delimiters(true);
      linebuf.set_skip_trailing_delimiters(true);

      string_tokenizer::iterator i   = linebuf.begin();
      string_tokenizer::iterator end = linebuf.end();

      string marker_name;
      ped_name = "";
      ind_names.resize(0);
      pair_ids.resize(0);

      if( !read_marker_pair_ids(i, end, marker_name, ped_name, ind_names) )
      {
        read_error("Marker skipped due to error...");
        continue;
      }

      if( m == 0 )
      {
        for( size_t pi = 0; pi < ind_names.size(); ++pi )
        {
          string ind1_name = ind_names[pi].first;
          string ind2_name = ind_names[pi].second;

          IBD::error_t err;
          id_pair p = ibd->get_pair(ped_name, ind1_name, ind2_name, err);
          size_t pair_num  = ibd->pair_index(p.first, p.second);

          if( pair_num >= ibd->pair_count() )
          {
            read_error("Error reading IBD STATE values for pairs ("
                        + ind1_name + "," + ind2_name + ") in pedigree '"
                        + ped_name + "'.");
          }

          pair_ids.push_back(p);
        }
      }

      input_next_line(file,line);

      map< vector<size_t>, double > current_ibd_map;

      while( line.size() > 0 )
      {
        linebuf.set_str(line);
        
        i   = linebuf.begin();
        end = linebuf.end();

        if( i == end )
        {
          read_error("Pair contains no data");
          continue;
        }

        string state_prob = *i++;

        vector<size_t> f_values;
        size_t pi;
        for( pi = 0; i != end && pi < ind_names.size(); ++pi )
        {
          string field0 = *i++;

          if( !field0.size() )
          {
            read_error("Invalid data format");
            break;
          }

          if( field0[0] != '-' )
          {
            if( field0.size() > 1 )
              field0 = field0[0];

            f_values.push_back(str2long(field0));
          }
        }
#if 0
        cout << "prop = " << state_prob << " ";
        for( size_t f = 0; f < f_values.size(); ++f )
          cout << f_values[f] << " ";
        cout << endl;
#endif
        if( pi < ind_names.size() || f_values.size() < ind_names.size() )
        {
          read_error("Invalid IBD state format for marker ("
                          + marker_name + ") in pedigree '"
                          + ped_name + "'. State skipped");
          continue;
        }

        //ibd_states.push_back(make_pair(str2doub(state_prob), f_values));

        if( current_ibd_map.find(f_values) != current_ibd_map.end() )
          current_ibd_map[f_values] = current_ibd_map[f_values] + str2doub(state_prob);
        else
          current_ibd_map[f_values] = str2doub(state_prob);

        input_next_line(file,line);
      }
#if 0
      cout << "current_ibd_map size = " << current_ibd_map.size() << endl;
#endif
      a_marker_ibd_state current_ibd_state;

      map< vector<size_t>, double >::const_iterator mi = current_ibd_map.begin();
      for( ; mi != current_ibd_map.end(); ++mi )
      {
        vector<size_t> i_state = mi->first;
        double         i_prob  = mi->second;

        current_ibd_state.push_back(make_pair(i_prob, i_state));
#if 0
        cout << i_prob << " ";
        for( size_t is = 0; is < i_state.size(); ++is )
          cout << i_state[is] << " ";
        cout << endl;
#endif
      }
#if 0
      cout << "current_ibd_states size = " << current_ibd_state.size() << endl;
#endif
      sped_ibd_state.push_back(current_ibd_state);
    }

    ibd_state_info sped_is_info;
    sped_is_info.pair_ids   = pair_ids;
    sped_is_info.ibd_states = sped_ibd_state;
    ibd->set_ibd_state(pair_ids[0].first->subpedigree(), sped_is_info);
#if 0
    cout << "ped_name = " << ped_name << endl;
    cout << "ind names size = " << ind_names.size() << endl;
    cout << "ped_ibd_states size = " << sped_ibd_state.size() << endl;
#endif
  }

  return true;
}

void
RefIBDReadFileStdIO::input_next_line(FILE *file, std::string &line)
{
  line.resize(0);
  if( !file || ferror(file) || feof(file) )
    return;

  do
  {
    line = strip_ws(getString(file, "\n\r", ""));

    if( !feof(file) && !ferror(file) )  // Eat end of line if not at EOF
    {
      getc(file);
      ++my_line;
    }
  }
  while( !ferror(file) && !feof(file) && line.size() && line[0] == '#' );

  if( line.size() && line[0] == '#' )
    line.resize(0);
}

bool
RefIBDReadFileStdIO::read_header(FILE *file, IBD *ibd,
                                 std::string&     line,
                                 ibd_option_type& ibd_option,
                                 vector<size_t>&  markers)
{
  input_next_line(file,line);

  ibd_option.old_ibd_format = true;

  if( line == "ANALYSIS" )
  {
    ibd_option.old_ibd_format = false;
    read_ibd_option(file, line, ibd_option);

    input_next_line(file,line);
  }

  ibd->set_ibd_option(ibd_option);

  if( line != "MARKERS" )
  {
    read_error("Can't find marker specification section xxx");
    return false;
  }

  std::vector<pair<string, bool> > marker_list;
  std::map<string, size_t>         marker_map;
  std::set<string>                 distance_set;

  input_next_line(file,line);

  while(line.size() > 0)
  {
    assert(ibd != NULL);

    // Need to get attribute from marker name for X-linkage.
    //
    string name = line;
    bool   is_x_linked = false;

    size_t pos  = line.find(',');

    if( pos > 0 && pos < line.size() )
    {
      string att = strip_ws(line.substr(pos+1, line.size() - (pos+1)), " ");
      if( toUpper(att) == "X_LINKED" )
      {
        name = line.substr(0, pos);
        is_x_linked = true;
      }
    }
    
    if( ibd_option.x_linked )
      is_x_linked = true;

    marker_list.push_back(make_pair(name, is_x_linked));

    size_t m_pos = name.find_last_of(' ', name.size());

    string region_name = name;
    string distance    = "";
    if( m_pos < name.size() )
    {
      region_name = name.substr(0, m_pos);
      distance = name.substr(m_pos+1, name.size());
    }

    marker_map[region_name] = marker_map[region_name] + 1;
    distance_set.insert(distance);

    input_next_line(file,line);
  }

  // For backward compatible with the old ibd files with/without cM
  //  after the marker name (new one has cM printed).
  // 
  for( size_t i = 0; i < marker_list.size(); ++i )
  {
    string name = marker_list[i].first;

    size_t m_pos = name.find_last_of(' ', name.size());

    string m_name = name.substr(0, m_pos);
    double m_dist = QNAN;

    if( m_pos < name.size() )
    {
      string dist = name.substr(m_pos+1, name.size());

      if(    dist.find_last_of('.', dist.size()) < dist.size()
          && !SAGE::isnan(str2doub(dist)) )
        m_dist = str2doub(dist);

      if( marker_map[m_name] > 1 )
        m_name = m_name + "_" + dist;
    }

    size_t m = ibd->marker_index(m_name);

    if( m >= ibd->marker_count() )
    {
      gmodel_type m_type = MLOCUS::AUTOSOMAL;
      if( marker_list[i].second )
        m_type = MLOCUS::X_LINKED;

      m = ibd->add_marker(m_name, m_dist, m_type);
    }

#if 0
  cout << "m = " << m << ", name = " << name
       << ", m_name = " << m_name << ", m_dist = " << m_dist;
  if( marker_list[i].second )
    cout << ", x_linked" << endl;
  else
    cout << ", autosomal" << endl; 
#endif

    markers.push_back(m);
  }

  return true;
}

bool
RefIBDReadFileStdIO::read_ibd_option(FILE *file, std::string &line, ibd_option_type& ibd_option)
{
  input_next_line(file,line);

  while(line.size() > 0)
  {
    // Need to get attribute from marker name for X-linkage.
    //
    size_t pos  = line.find('=');

    if( pos != size_t(-1) )
    {
      string opt = strip_ws(line.substr(0,     pos), " ");
      string val = strip_ws(line.substr(pos+1, line.size()-(pos+1)), " ");

      if( toUpper(opt) == "TITLE" )
        ibd_option.title = val;
      else if( toUpper(opt) == "REGION" )
      {
        size_t pos1  = val.find(',');

        if( pos1 != size_t(-1) )
        {
          string att = strip_ws(val.substr(pos1+1, val.size()-(pos1+1)), " ");

          if( toUpper(att) == "X_LINKED" )
            ibd_option.x_linked = true;

          ibd_option.region = strip_ws(val.substr(0, pos1), " ");
        }
        else
          ibd_option.region = val;
      }
      else if( toUpper(opt) == "MAX_PEDIGREE" )
        ibd_option.max_pedigree = val;
      else if( toUpper(opt) == "SCAN_TYPE" )
        ibd_option.scan_type = val;
      else if( toUpper(opt) == "ALLOW_LOOPS" )
        ibd_option.allow_loops = val;
      else if( toUpper(opt) == "IBD_MODE" )
      {
        ibd_option.ibd_mode = val;

        size_t pos1  = val.find(',');

        if( pos1 == size_t(-1) )
        {
          ibd_option.exact = false;
          ibd_option.ibd_mode = val;
        }
        else
          ibd_option.ibd_mode = strip_ws(val.substr(0, pos1), " ");
      }
      else if( toUpper(opt) == "SPLIT_PEDIGREES" )
        ibd_option.split_pedigrees = val;
      else if( toUpper(opt) == "USE_SIMULATION" )
        ibd_option.use_simulation = val;
    }

    input_next_line(file,line);
  }

  return true;
}

bool
RefIBDReadFileStdIO::read_pair_ids(string_tokenizer::iterator &tok, 
                                   string_tokenizer::iterator &end,
                                   std::string &ped_name,
                                   std::string &ind1_name,
                                   std::string &ind2_name)
{
  ped_name.resize(0);
  ind1_name.resize(0);
  ind2_name.resize(0);

  if( tok != end )
    ped_name = strip_ws(*tok++);
  if( tok != end )
    ind1_name = strip_ws(*tok++);
  if( tok != end )
    ind2_name = strip_ws(*tok);

  if( !ped_name.size() )
  {
    read_error("cannot read pedigree id");
    return false;
  }

  if( !ind1_name.size() )
  {
    read_error("cannot read first individual id");
    return false;
  }

  if( !ind2_name.size() )
  {
    read_error("cannot read second individual id");
    return false;
  }

  return true;
}

bool
RefIBDReadFileStdIO::read_marker_pair_ids(string_tokenizer::iterator     &tok, 
                                          string_tokenizer::iterator     &end,
                                          string                         &marker_name,
                                          string                         &ped_name,
                                          vector< pair<string, string> > &ind_names)
{
  marker_name.resize(0);
  ped_name.resize(0);
  ind_names.resize(0);

  if( tok != end )
    marker_name = strip_ws(*tok++);

  string pi_pair = strip_ws(*tok++);

  size_t p_pos = pi_pair.find_first_of('(');
  ped_name = pi_pair.substr(0, p_pos);

  size_t i_pos   = pi_pair.find_first_of(',');
  string i1_name = pi_pair.substr(p_pos+1, i_pos-p_pos-1);
  string i2_name = pi_pair.substr(i_pos+1, pi_pair.size()-1-i_pos-1);

  ind_names.push_back(make_pair(i1_name, i2_name));

  while( tok != end )
  {
    pi_pair = strip_ws(*tok++);

    p_pos   = pi_pair.find_first_of('(');
    i_pos   = pi_pair.find_first_of(',');
    i1_name = pi_pair.substr(p_pos+1, i_pos-p_pos-1);
    i2_name = pi_pair.substr(i_pos+1, pi_pair.size()-1-i_pos-1);

    ind_names.push_back(make_pair(i1_name, i2_name));
  }

#if 0
  cout << marker_name << " ";
  cout << ped_name << " ";
  for( size_t i = 0; i < ind_names.size(); ++i )
    cout << "(" << ind_names[i].first << "," << ind_names[i].second << ") ";
  cout << endl;
#endif

  if( !ped_name.size() )
  {
    read_error("cannot read pedigree id");
    return false;
  }

  if( !ind_names.size() )
  {
    read_error("cannot read first individual id");
    return false;
  }

  return true;
}

//
//----------------------------------------------------------------
//

RefIBDWriteFile::RefIBDWriteFile(const std::string &fname,
                                 std::ostream &output_messages)
               : messages(output_messages), file(fname.c_str()), filename(fname)
{
  reset(fname);

  if( !file )
  {
    write_error("Cannot write to filename '" + fname + "'.");
    return;
  }
}

RefIBDWriteFile::~RefIBDWriteFile()
{}

bool
RefIBDWriteFile::output_probability_header(const IBD *ibd)
{
  if( !file )
  {
    write_error("Cannot write to filename '" + filename + "'.");
    return false;
  }

  if( !is_valid_ibd(ibd) )
    return false;

  my_marker_count = ibd->marker_count();

  // DISPLAY Header
  file << "IBD File 2.0 : This File is automatically generated.  Do NOT edit!\n"
       << "#=================================================================\n";

  output_header_info(ibd);

  file << "#===============================================================================" << endl;

  // Write the Table of f0, f2 values
  size_t ped_col = 15;   // Initial Column width for pedigrees
  size_t ind_col = 10;   // Initial Column width for individuals
  size_t cw      = 17;   // Column width for Markers

  string dash = "-------------------------------------------------------";
  while( dash.size() < max(ind_col,ped_col) )
    dash += dash;

  // DISPLAY headers of columns
  file << "#"
       << right
       << setw(ped_col-1) << "Pedigree" << "  "       // Print Headers
       << setw(ind_col)   << "Ind 1"    << "  "
       << setw(ind_col)   << "Ind 2"    << "  ";

  for( size_t j = 0; j < ibd->marker_count(); ++j )
  {
    string m_name = ibd->marker_name(j);
    size_t m_pos  = m_name.find_last_of(' ', m_name.size());

    m_name = m_name.substr(0, m_pos);

    file << setw(cw - 4) << m_name << " f0,  "
         << setw(cw - 9) << m_name << " f1m-f1p,  "
         << setw(cw - 4) << m_name << " f2";

    if( j+1 < ibd->marker_count() )
      file << ",  ";
    else
      file << ".";
  }

  file << endl;

  file << "#"
       << dash.substr(0,ped_col-1) << "  "
       << dash.substr(0,ind_col) << "  "
       << dash.substr(0,ind_col);

  for(size_t j = 0; j < ibd->marker_count(); ++j )
    file << "  " << dash.substr(0,cw)
         << "  " << dash.substr(0,cw)
         << "  " << dash.substr(0,cw);

  file << endl;

  return true;
}
/*
bool
RefIBDWriteFile::output_covariance_header(const IBD *ibd)
{
  if( !file )
  {
    write_error("Cannot write to filename '" + filename + "'.");
    return false;
  }

  if( !is_valid_ibd(ibd) )
    return false;

  my_marker_count = ibd->marker_count();

  // DISPLAY Header
  file << "IBD COVARIANCE File 1.0 : This File is automatically generated.  Do NOT edit!\n"
       << "#============================================================================\n";
  //file << endl;

  output_header_info(ibd);

  file << "#===============================================================================" << endl;

  // Write the Table of f0, f2 values
  size_t ped_col = 15;   // Initial Column width for pedigrees
  size_t ind_col = 10;   // Initial Column width for individuals
  size_t cw      = 17;   // Column width for Markers

  string dash = "-------------------------------------------------------";
  while( dash.size() < max(ind_col,ped_col) )
    dash += dash;

  // DISPLAY headers of columns
  file << "#"
       << right
       << setw(ped_col-1) << "Pedigree" << "  "       // Print Headers
       << setw(ind_col)   << "P1 Ind 1" << "  "
       << setw(ind_col)   << "P1 Ind 2" << "  "
       << setw(ind_col)   << "P2 Ind 1" << "  "
       << setw(ind_col)   << "P2 Ind 2";

  for( size_t j = 0; j < ibd->marker_count(); ++j )
  {
    string m_name = ibd->marker_name(j);
    size_t m_pos  = m_name.find_last_of(' ', m_name.size());

    m_name = m_name.substr(0, m_pos);

    file << "  " << setw(cw) << m_name;
  }

  file << endl;

  file << "#"
       << dash.substr(0,ped_col-1) << "  "
       << dash.substr(0,ind_col) << "  "
       << dash.substr(0,ind_col) << "  "
       << dash.substr(0,ind_col) << "  "
       << dash.substr(0,ind_col);

  for(size_t j = 0; j < ibd->marker_count(); ++j )
    file << "  " << dash.substr(0,cw);

  file << endl;

  return true;
}
*/
bool
RefIBDWriteFile::output_state_header(const IBD *ibd)
{
  if( !file )
  {
    write_error("Cannot write to filename '" + filename + "'.");
    return false;
  }

  if( !is_valid_ibd(ibd) )
    return false;

  my_marker_count = ibd->marker_count();

  // DISPLAY Header
  file << "IBD STATE File 1.0 : This File is automatically generated.  Do NOT edit!\n"
       << "#=======================================================================\n";

  output_header_info(ibd);

  file << "#===============================================================================" << endl;

  // Find the maximum ped, ind name length.
  size_t ped_col = 0;   // Initial Column width for pedigrees
  size_t ind_col = 0;   // Initial Column width for individuals

  for( size_t i = 0; i < ibd->pair_count(); ++i )
  {
    std::string ped_name;
    std::string ind1_name;
    std::string ind2_name;

    if( !ibd->get_pair(i, ped_name, ind1_name, ind2_name) )
    {
      // Warn pair skipped?
      continue;
    }

    ped_col = std::max(ped_name.size(), ped_col);
    ind_col = std::max(max(ind1_name.size(), ind2_name.size()), ind_col);
  }

  // Find the maximum marker name length.
  size_t m_col = 10;

  for( size_t k = 0; k < ibd->marker_count(); ++k )
  {
    string m_name = ibd->marker_name(k);
    //size_t m_pos = m_name.find_last_of(' ', m_name.size());

    //if( m_pos < m_name.size() )
    //  m_name = m_name.substr(0, m_pos);

    m_col = std::max(m_name.size(), m_col);
  }

  string dash = "-------------------------------------------------------";
  while( dash.size() < max(ind_col,ped_col) )
    dash += dash;

  // DISPLAY headers of columns
  file << left << setw(m_col) << "#Marker"
       << "  "
       << "Ped(Ind1,Ind2)"
       << endl;

  file << "#"
       << dash.substr(0,m_col-1)
       << "  "
       << dash.substr(0,14);
       //<< endl;

  return true;
}

bool
RefIBDWriteFile::output_ibd_probability(const IBD *ibd)
{
  if( !file )
  {
    write_error("Cannot write to filename '" + filename + "'.");
    return false;
  }

  if( !is_valid_ibd(ibd) )
    return false;

  if( my_marker_count != ibd->marker_count() )
  {
    write_error("Cannot output different set of markers to IBD file");
    return false;
  }

  // Write the Table of f0, f2 values
  size_t ped_col = 15;   // Initial Column width for pedigrees
  size_t ind_col = 10;   // Initial Column width for individuals
  size_t cw      = 17;   // Column width for Markers

  string dash = "-------------------------------------------------------";
  while( dash.size() < max(ind_col,ped_col) )
    dash += dash;

  // Write out the Pairs
  for( size_t i = 0; i < ibd->pair_count(); ++i )
  {
    std::string ped_name;
    std::string ind1_name;
    std::string ind2_name;

    if( !ibd->get_pair(i, ped_name, ind1_name, ind2_name) )
    {
      // Warn pair skipped?
      continue;
    }

    file << setw(ped_col - 1) << ped_name   << ",  "
         << setw(ind_col - 1) << ind1_name  << ",  "
         << setw(ind_col - 1) << ind2_name  << ",  ";

    for( size_t k = 0; k < ibd->marker_count(); ++k )
    {
      if (k) file << "  ";

      double f0;
      double f2;
      double f1mp;

      ibd->get_ibd(i, k, f0, f1mp, f2);

      if( SAGE::isnan(f0) || SAGE::isnan(f2) )
      {
        file << setw(cw) << dash.substr(0,cw)
             << "  "
             << setw(cw) << dash.substr(0,cw)
             << "  "
             << setw(cw) << dash.substr(0,cw);
      }
      else
      {
        file << setw(cw) << doub2str(f0, cw, -1, ios::showpoint | ios::fixed);

        if( !SAGE::isnan(f1mp) )
        {
          file << "  "
               << setw(cw) << doub2str(f1mp, cw, -1, ios::showpoint | ios::fixed);
        }
        else
          file << "  "
               << setw(cw) << dash.substr(0,cw);

        file << "  "
             << setw(cw) << doub2str(f2, cw, -1, ios::showpoint | ios::fixed);
      }
    }

    file << endl;
  }

  return true;
}

bool
RefIBDWriteFile::output_ibd_state(const IBD *ibd)
{
  if( !file )
  {
    write_error("Cannot write to filename '" + filename + "'.");
    return false;
  }

  if( !is_valid_ibd(ibd) )
    return false;

  if( my_marker_count != ibd->marker_count() )
  {
    write_error("Cannot output different set of markers to IBD file");
    return false;
  }

  // Find the maximum ped, ind name length.
  size_t ped_col = 0;   // Initial Column width for pedigrees
  size_t ind_col = 0;   // Initial Column width for individuals

  vector<size_t> ind1_cols;
  vector<size_t> ind2_cols;

  for( size_t i = 0; i < ibd->pair_count(); ++i )
  {
    std::string ped_name;
    std::string ind1_name;
    std::string ind2_name;

    if( !ibd->get_pair(i, ped_name, ind1_name, ind2_name) )
    {
      // Warn pair skipped?
      continue;
    }

    ped_col = std::max(ped_name.size(), ped_col);
    ind_col = std::max(max(ind1_name.size(), ind2_name.size()), ind_col);
    ind1_cols.push_back(ind1_name.size());
    ind2_cols.push_back(ind2_name.size());
  }

  // Find the maximum marker name length.
  size_t m_col = 10;

  for( size_t k = 0; k < ibd->marker_count(); ++k )
  {
    string m_name = ibd->marker_name(k);
    //size_t m_pos = m_name.find_last_of(' ', m_name.size());

    //if( m_pos < m_name.size() )
    //  m_name = m_name.substr(0, m_pos);

    m_col = std::max(m_name.size(), m_col);
  }

  // Write out the ibd state for each marker.
  for( size_t k = 0; k < ibd->marker_count(); ++k )
  {
    string m_name = ibd->marker_name(k);
    //size_t m_pos = m_name.find_last_of(' ', m_name.size());

    //m_name = m_name.substr(0, m_pos);

    file << endl
         << left << setw(m_col) << m_name;

    for( size_t i = 0; i < ibd->pair_count(); ++i )
    {
      std::string ped_name;
      std::string ind1_name;
      std::string ind2_name;

      if( !ibd->get_pair(i, ped_name, ind1_name, ind2_name) )
      {
        // Warn pair skipped?
        continue;
      }

      file << "  "
           << setw(ped_col) << ped_name   << "("
           << setw(ind1_cols[i]) << ind1_name  << ","
           << setw(ind2_cols[i]) << ind2_name  << ")";
    }
    file << endl;

    vector< ibd_lvec_probability > ibd_state;

    ibd->get_ibd_state(k, ibd_state);

    for( size_t i = 0; i < ibd_state.size(); ++i )
    {
      file << setw(m_col) << doub2str(ibd_state[i].first, m_col, -1, ios::showpoint | ios::fixed);

      for( size_t j = 0; j < ibd_state[i].second.size(); ++j )
      {
        int i_val = ibd_state[i].second[j];

        if( i_val == 3 )
          file << "  " << right << setw(ped_col + ind1_cols[j] + 3) << "1m"
               << setw(ind2_cols[j]) << " ";
        else if( i_val == 4 )
          file << "  " << right << setw(ped_col + ind1_cols[j] + 3) << "1p"
               << setw(ind2_cols[j]) << " ";
        else
          file << "  " << right << setw(ped_col + ind1_cols[j] + 2) << i_val
               << setw(ind2_cols[j]+1) << " ";
      }
      file << endl;
    }
  }

  return true;
}

bool
RefIBDWriteFile::is_valid_ibd(const IBD *ibd)
{
  if( !ibd )
  {
    write_error("NULL IBD data structure");
    return false;
  }

  if( !ibd->built() )
  {
    write_error("Cannot create IBD file");
    return false;
  }

  if( !ibd->marker_count() )
  {
    write_error("No markers to output");
    return false;
  }

  if( !ibd->pair_count() )
  {
    write_error("No pairs to output");
    return false;
  }

  return true;
}

void
RefIBDWriteFile::output_header_info(const IBD *ibd)
{
  const ibd_option_type& opt = ibd->get_ibd_option();

  file << "#\nANALYSIS\n#-------\n";
  file << "title           = " << opt.title << endl;
  file << "region          = " << opt.region;
  if( opt.x_linked )
    file << ", x_linked";
  file << endl;
  if( opt.exact )
    file << "max_pedigree    = " << opt.max_pedigree << endl;
  file << "scan_type       = " << opt.scan_type << endl;
  file << "allow_loops     = " << opt.allow_loops << endl;
  file << "ibd_mode        = " << opt.ibd_mode;
  if( opt.exact )
    file << ", exact";
  file << endl;
  file << "split_pedigrees = " << opt.split_pedigrees << endl;
  file << "use_simulation  = " << opt.use_simulation << endl;
  file << endl;

  // Print marker header and marker names
  file << "#\nMARKERS\n#------\n";

  for( size_t k = 0; k < ibd->marker_count(); ++k )
  {
    file << ibd->marker_name(k);
    file << " " << doub2str(ibd->get_marker_info(k).distance, 0, 1, ios::showpoint | ios::fixed);
    if( opt.x_linked )
      file << ", x_linked";
    file << endl;
  }
  file << endl;

  return;
}

}
