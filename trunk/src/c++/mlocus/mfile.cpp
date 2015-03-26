#include <string>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include "LSF/parse_ops.h"
#include "error/errormanip.h"
#include "mlocus/mfile.h"

namespace SAGE   {
namespace MLOCUS {

InheritanceModelFile::InheritanceModelFile (const cerrorstream& err)
  : my_marker_verbose_output((uint)-1),
    my_genotype_verbose_output((uint)-1),
    my_valid(true),
    read_distances(false),
    errors(err)
{ }

InheritanceModelFile::~InheritanceModelFile() { }

bool InheritanceModelFile::input
  (inheritance_model_map& m, const string& fname)
{
  return input(m, fname, '/', cout);
}

bool InheritanceModelFile::input
  (inheritance_model_map& m, const string& fname, char separator)
{
  return input(m, fname, separator, cout);
}

bool InheritanceModelFile::input
  (inheritance_model_map& m, const string& fname, const RPED::PhenotypeReaderInfo& pi)
{
  return input(m, fname, pi, cout);
}

bool InheritanceModelFile::input
  (inheritance_model_map& m, const string& fname, char separator, ostream& messages)
{
  RPED::PhenotypeReaderInfo pri;

  pri.set_allele_delimiter(separator);

  return input(m, fname, pri, messages);
}

bool InheritanceModelFile::input
  (inheritance_model_map& m, const string& fname, const RPED::PhenotypeReaderInfo& pi, ostream& messages)
{
  //lint --e{534}

  if( !fname.size() )
  {
    errors << priority(critical) << "No Marker Locus Description file specified." << endl;
    return false;
  }

  std::ifstream infile( fname.c_str() );

  if(!infile.good())
  {
    errors << priority(critical) << "Unable to open Marker Locus Description file '" << fname
                    << "'. Please check your file." << endl;
    invalidate();
    return false;
  }

  if(!infile.good()) return false;

  kill_ws(infile);

  // Determine if the file is LSF or not

  if(infile.peek() == '#')
  {
    // Read the rest of the line to determine if it's of a type we can
    // handle
    
    string first_line = getString(infile, "\n\r", "");
    
    if(first_line == "#! Penetrance Probability File, version 1.0.")
    {
      errors << priority(critical) << "The current version of S.A.G.E. cannot "
             << "handle penetrance probabilities based on parental genotypes.  "
             << "Please use the Type Probability file instead." << endl;
             
      exit(1);
    }
    else if(first_line == "#! Type Probability File, version 1.0.")
    {
      errors << priority(critical) << "Type Probability File, version 1.0 "
             << "had serious errors, and is no longer supported.  Please "
             << "run your data with newer version of SEGREG." << endl;
             
      exit(1);
    }
    else if(first_line != "#! Type Probability File, version 1.1.")
    {
      errors << priority(critical) << "Unrecognized file format for "
             << "Trait or Locus Description file.  Please check your file."
             << endl;
             
      exit(1);
    }
    
    // Ok, so it's a type probability file.

    LSF_input load_state(infile, cout);
    assert(load_state.good());

    LSFBase* my_models = new LSFBase("MarkerModels");
    load_state.input_to(my_models, false);

    if( !my_models || !my_models->List() )
    {
      errors << priority(fatal)
             << "Error reading Trait or Locus Description file....  Exiting."
             << endl;

      exit(EXIT_FAILURE);
    }

    LSFInheritanceModelParser mparse(pi.get_allele_delimiter(), errors);
    
    return mparse.input_to(m, my_models);
  }

  // Ok, old-fashioned file

  char separator = pi.get_allele_delimiter();
  string mv  = getString( infile, "=\n\r", "");
  string mv2 = "";
  
  if (infile.peek() == '=' && 
      ( toUpper(mv = mv.substr(0,mv.find_last_not_of(" \t")+1)) == "MISSING"))
  {
    infile.get();  // Grab the =

    kill_ws(infile, "", " \t");

    mv = getString( infile, "\n\r", "");
    mv = mv.substr(0, std::min(mv.find_last_not_of(" \t")+1, (size_t) 9));

    if( strchr(mv.c_str(), separator) )
    {
      // Split the missing value into two strings.

      mv2 = strip_ws( mv.substr( mv.find_first_of(separator), mv.size()  ), "/ \t" );
      mv  = strip_ws( mv.substr( 0, mv.find_first_of(separator) ), " \t" );

      if(mv != mv2)
      {
        errors << priority(warning) << "The Missing Value code must use the "
               << "same value for both alleles.  Will use " << mv << separator << mv
               << " for the missing value." << endl;
               
        mv2 = mv;
      }
    }
  }
  else
  {
    infile.seekg(0);
    mv = "";
  }

  if( pi.get_allele_missing().size() )
    mv = pi.get_allele_missing();

  mv2 = mv;

  kill_ws(infile);

  while ( !infile.eof() )
  {
    string name = strip_ws(getString( infile, "\n\r", "" ));

//  name.resize( std::min(name.find_last_not_of(" \t")+1, (size_t) 10) );

    if(m.name_find(name) != m.name_end())
    {
      errors << priority(warning) << "Marker '" << name << "' is already "
             << "present in the set of markers.  It will be overwritten." << endl;
    }

    genotype_model gm(name);
    
    gm.set_unphased_separator(separator);

    if(mv.size())
      gm.set_missing_allele_name(mv);

    get_alleles(infile, gm, pi);
    
    if(test_codominant(infile))
    {
      penetrance_model& pm = m[name] = penetrance_model(gm, true, true);

      pm.set_name(name);

      if(mv.size())
      {
        pm.set_missing_phenotype_name(mv);
      
        if(mv2.size())
        {
          string sep = gm.separators();

          pm.alias_phenotype(pm.get_missing_phenotype_id(), mv + sep[0] + mv2);
          pm.alias_phenotype(pm.get_missing_phenotype_id(), mv + sep[1] + mv2);
          pm.alias_phenotype(pm.get_missing_phenotype_id(), mv + sep[2] + mv2);
        }
      }
    }
    else
    {
      penetrance_model& pm = m[name] = penetrance_model(gm, false, false);

      pm.set_name(name);

      if(mv.size())
      {
        pm.set_missing_phenotype_name(mv);
      
        if(mv2.size())
        {
          string sep = gm.separators();

          pm.alias_phenotype(pm.get_missing_phenotype_id(), mv + sep[0] + mv2);
          pm.alias_phenotype(pm.get_missing_phenotype_id(), mv + sep[1] + mv2);
          pm.alias_phenotype(pm.get_missing_phenotype_id(), mv + sep[2] + mv2);
        }
      }

      get_phenotypes(infile, pm);

      pm.make_consistent();
    }

    kill_ws(infile);
  }

//  output_markers(messages, m);

  return true;
}


void InheritanceModelFile::get_alleles(istream& infile, genotype_model& gm, const RPED::PhenotypeReaderInfo& pi)
{
  //lint --e{534}

#if 0
  cout << "Passed pi values : " << endl
       << "  pi.allele_delimiter = " << pi.get_allele_delimiter() << endl
       << "  pi.allele_missing   = " << pi.get_allele_missing() << endl
       << "  pi.allele_frequency = " << pi.get_allele_adjustment() << ", "
       << pi.get_min_allele_freq() << ", " << pi.get_max_allele_freq() << endl;
#endif

  const string& name = gm.name();

  string al, dtemp;
  double d = 0.0;
  double total = 0.0;

  bool first    = true;
  bool has_freq = true;

  do
  {
    kill_ws(infile);
    al = getString(infile, "=\t\n\r", "").substr(0,4);

    al = strip_ws(al," \t");

    kill_ws(infile, " \t");

    test_eof(infile);
    
    if(first)
    {
      if(infile.peek() != '=')
        has_freq = false;

      first = false;
    }

    if(has_freq)
    {
      if( infile.get() != '=')
      {
        errors << priority(error) << "Missing \'=\' for allele " << al
               << " in marker " << name << ". Unable to continue"
               << " reading Marker Locus Description File." << endl;
        exit(1);
      }

      kill_ws(infile," \t");

      test_eof(infile);

      dtemp = getString(infile, ";\n\r");
      d = atof( dtemp.c_str() );

      if( pi.get_allele_adjustment() == RPED::PhenotypeReaderInfo::comp )
        d = 1. - d;
      else if(    pi.get_allele_adjustment() != RPED::PhenotypeReaderInfo::none
               && pi.get_allele_adjustment() != RPED::PhenotypeReaderInfo::equal )
      {
        if(    (    pi.get_allele_adjustment() == RPED::PhenotypeReaderInfo::min
                 || pi.get_allele_adjustment() == RPED::PhenotypeReaderInfo::min_max )
            && ( d < pi.get_min_allele_freq() ) )
          d = pi.get_min_allele_freq();

        if(    (    pi.get_allele_adjustment() == RPED::PhenotypeReaderInfo::max
                 || pi.get_allele_adjustment() == RPED::PhenotypeReaderInfo::min_max )
            && ( d > pi.get_max_allele_freq() ) )
          d = pi.get_max_allele_freq();
      }

      if(d < 0.0)
      {
        errors << priority(warning) << "Negative Allele Frequency for "
               << "allele " << al << "in marker " << name << " invalid."
               << "  Setting to zero." << endl;
        d = 0.0;
      }
    }
    else
    {
      if( infile.peek() == '=')
      {
        errors << priority(error) 
               << "Frequencies must be present (or missing) for all alleles in marker "
               << name << ". Unable to continue reading Marker Locus Description File."
               << endl;
        exit(1);
      }
    }

    gm.add_allele(al, d);

    total += d;

    kill_ws(infile);

  } while (!infile.eof() && infile.peek() != ';');

  if( has_freq && pi.get_allele_adjustment() != RPED::PhenotypeReaderInfo::equal )
  {
    if(total == 0.0)
    {
      errors << priority(error) << "Allele frequencies for marker "
             << name << " sum to zero.  Please fix your file." << endl;

      exit(1);
    }

    // Only print message if far off.
    if(total <= 0.99 || total >= 1.01)
      errors << priority(warning) << "Allele frequencies for marker "
             << name << " do not sum to 1.  Allele frequencies will be"
             << " normalized to 1." << endl;

    //  Normalize if != 1.
    if(total != 1.0)
    {
      gm.normalize();
    }
  }
  else
  {
    // Set all frequencies to equiprobable
    gm.normalize();
  }

  test_eof(infile);

  infile.get();
}

void InheritanceModelFile::get_phenotypes
    (istream& infile, penetrance_model& pm)
{
  //lint --e{534}

  const string& name = pm.name();

  do
  {
    kill_ws(infile);
    string p = getString( infile, "=\n\r", "\t" );

    p = strip_ws(p, " \t");

    pm.add_phenotype(p);
    
    uint pid = pm.get_phenotype_id(p);

    kill_ws(infile, " \t");

    if(infile.get() != '=')
    {
      errors << priority(error) << "Missing \'=\' for phenotype " 
             << p << " in marker " << name << ". Unable to continue"
             << " reading Marker Locus Description File." << endl;

      exit(1);
    }

    kill_ws(infile, " \t");

    //lint -e{734}
    char open_char = infile.get();

    if(open_char != '{' && open_char != '<') // } <- to make bracket matching work

    {
      errors << priority(error) << "Missing opening character (\'{\' or \'<\' for phenotype " 
             << p << " in marker " << name << ". Unable to continue"
             << " reading Marker Locus Description File." << endl;

      exit(1);
    }

    do
    {
      kill_ws(infile);

      // { <- to make bracket matching work
      string genotype = getString( infile, "=}>", ",");

      genotype = strip_ws(genotype, " \t");

      unphased_genotype gtype = pm.get_unphased_genotype(genotype);

      if(!gtype.is_valid())           // We can't find the genotype.
      {
        errors << priority(warning) << "Unknown genotype " << genotype
               << " in marker " << name << ".  Genotype "
               << "will be ignored." << endl;

        if(infile.peek() == '=')             // If value attached, ignore it.
        {
          infile.get();
           // { <- to make bracket matching work
          getString(infile, "}>", ",");
        }
      }
      else
        if (infile.peek() != '=')
          pm.add_unphased_penetrance(1.0, pid, gtype.get_id());
        else
        {
          infile.get();
          kill_ws(infile);
          // { <- to make bracket matching work
          genotype = getString( infile, "}>", ",");
          pm.add_unphased_penetrance(atof( genotype.c_str() ), pid, gtype.get_id());
        }

      kill_ws(infile);
    } 
    while ( !infile.eof() && (infile.peek() != '}' && infile.peek() != '>') );

    test_eof(infile);

    infile.get();                    // Get rid of end bracket
    kill_ws(infile);
  }
  while ( !infile.eof() && infile.peek() != ';' );

  test_eof(infile);

  infile.get();
}

bool InheritanceModelFile::test_codominant(istream& infile)
{
  //lint --e{534}

  kill_ws(infile);

  test_eof(infile);

  if(infile.peek() == ';')
  {
    infile.get();

    return true;
  }

  return false;
}

void InheritanceModelFile::test_eof(const istream& infile)
{
  if(infile.eof())
  {
    errors << priority(error) << "Early end of file in Marker Locus "
           << "Description File. Please check your file." << endl;

    exit(1);
  }
}

void InheritanceModelFile::output_markers(ostream& o, const inheritance_model_map& m)
{
  o << "MARKER MODELS:\n============================" << endl << endl;

  uint v = my_marker_verbose_output;

  inheritance_model_map::index_const_iterator i = m.index_begin();

  for( ; v && i != m.index_end(); ++i, --v)
  {
    output_marker(o, i->second, my_genotype_verbose_output);
  }
}

void
InheritanceModelFile::output_marker
   (ostream& o, const penetrance_model& pm, uint verbosity) const 
{
  //lint --e{534}

  o << pm.name();
  if( pm.is_x_linked() )
    o << ", X-linked";
  else if( pm.is_y_linked() )
    o << ", Y-linked";
  o << endl;

  o << "----------------------------" << endl << endl << "Alleles:" << endl;;

  uint column = 0;

  o.precision(6);
  o << setiosflags(ios::fixed);

  for(allele_iterator i = pm.allele_begin(); i != pm.allele_end(); ++i)
  {
    if( i->name().find('~') < i->name().size() )
      continue;

    o << i->name() << "\t= " << i->frequency() << "\t";

    ++column;

    if(column == 4)
    {
      o << endl;
      column = 0;
    }
  }

  if(column) o << endl;

  o << endl << "Genotype    Genotypic Frequency" << endl;

  if( pm.is_y_linked() )
  {
    for(allele_iterator i = pm.allele_begin(); i != pm.allele_end(); ++i)
    {
      if( i->name().find('~') < i->name().size() )
        continue;

      o << setw(8) << right << i->name() + pm.gmodel().unphased_separator() + "X" << "    "
        << setw(9) << i->frequency() << endl;
    }
  }
  else
  {
    set<string> geno_names;

    unphased_genotype_iterator i = pm.unphased_genotype_begin();
    for( ; verbosity && i != pm.unphased_genotype_end(); --verbosity, ++i)
    {
      if( i->name().find('~') < i->name().size() )
        continue;

      geno_names.insert(i->name());
    }

    for( set<string>::iterator si = geno_names.begin(); si != geno_names.end(); ++si )
    {
      o << setw(8) << right << *si << "    "
        << setw(9) << pm.get_unphased_genotype(*si).frequency() << endl;
    }

    if(i != pm.unphased_genotype_end())
      o << " .\n .\n ." << endl;

    if( pm.is_x_linked() )
    {
      for(allele_iterator i = pm.allele_begin(); i != pm.allele_end(); ++i)
      {
        if( i->name().find('~') < i->name().size() )
          continue;

        o << setw(8) << right << i->name() + pm.gmodel().unphased_separator() + "Y" << "    "
          << setw(9) << i->frequency() << endl;
      }
    }
  }

  o.setf( (ios::fmtflags)0, ios::floatfield);

  o << endl << flush;
}

bool InheritanceModelFile::output
  (const inheritance_model_map& m, const string& filename, ostream& messages)
{
  ofstream genome_file;
  genome_file.open(filename.c_str());

  output_markers(genome_file, m);

  return true;
}

bool InheritanceModelFile::output
  (const inheritance_model_map& m, const string& filename)
{
  ofstream genome_file;
  genome_file.open(filename.c_str());

  output_markers(genome_file, m);

  return true;
}

bool InheritanceModelFile::output
  (const inheritance_model_map& m, ostream& o)
{
  output_markers(o, m);

  return true;
}

LSFInheritanceModelParser::LSFInheritanceModelParser
    (char sep, const cerrorstream& c)
  : separator(sep), errors(c)
{ }

bool LSFInheritanceModelParser::input_to(inheritance_model_map& m, LSFBase* b)
{
  LSFList::const_iterator i;
  for( i = b->List()->begin(); i != b->List()->end(); ++i )
  {
    if( !*i || !(*i)->name().size() || !(*i)->attrs())
      continue;

    string name = toUpper( (*i)->name() );

    if( name != "MODEL" )
    {
      errors << priority(warning) << "Unrecognized parameter in "
             << "Type Probability file." << endl;
      continue;
    }
    parse_model(m, *i);
  }

  return true;
}

void LSFInheritanceModelParser::parse_model(inheritance_model_map& m, LSFBase* model)
{
  string name = model->attrs()->StringAttr("VALUE");

  if(name.size() == 0) return;

  if(m.name_find(name) != m.name_end())
  {
    errors << priority(warning) << "Marker '" << name << "' is already "
           << "present in the set of markers.  It will be overwritten." << endl;
  }

  genotype_model gm(name + "~~~");
  
  gm.set_unphased_separator(separator);
  gm.set_missing_allele_name("MISSING");

  if(!get_alleles(model, gm)) return;

  penetrance_model& pm = m[name] = penetrance_model(gm, false, false);

  pm.set_name(name + "~~~");
  pm.set_missing_phenotype_name("MISSING");
  
  if(!get_phenotypes(model, pm)) return;

  pm.make_consistent();
}

bool LSFInheritanceModelParser::get_alleles(LSFBase* model, genotype_model& gm)
{
  LSFBase* alleles = model->List()->front();

  if(!alleles || !alleles->List() || toUpper(alleles->name()) != "ALLELES")
  {
    errors << priority(error) << "Format wrong in Type Probability file.  "
           << "The alleles are not found in the model.  Please check your file."
           << endl;
    return false;
  }

  double total = 0.0;
  
  LSFList::const_iterator i;
  for( i = alleles->List()->begin(); i != alleles->List()->end(); ++i )
  {
    if( !*i || !(*i)->name().size() )
      continue;

    string name = (*i)->name();
    
    AttrList::const_iterator iter = (*i)->attrs()->find("VALUE");
    if(iter == (*i)->attrs()->end())
    {
      errors << priority(warning) << "Malformed allele (no allele value).  Please check your file."
             << endl;
             
      continue;
    }
    
    double d = iter->second.Real();

    if(!finite(d) || d <= 0.0 || d > 1.0)
    {
      errors << priority(warning) << "Malformed allele (allele value out of range).  "
             << "Please check your file." << endl;
             
      continue;
    }
    total += d;
    
    gm.add_allele(name, d);
  }

  if(total == 0.0)
  {
    errors << priority(error) << "No valid alleles for model '"
           << gm.name() << "'. Model will be ignored." << endl;

    return false;
  }

  // Only print message if far off.
  else if(total <= 0.99 || total >= 1.01)

  //  Normalize if != 1.
  if(total != 1.0)
  {
    gm.normalize();
  }

  return true;
}

bool LSFInheritanceModelParser::get_phenotypes(LSFBase* model, penetrance_model& pm)
{
  const string& name = pm.name();

  // for each pedigree

  LSFList::const_iterator i = model->List()->begin();

  ++i;

  for( ; i != model->List()->end(); ++i)
  {
    if(!(*i) || !(*i)->List() || toUpper((*i)->name()) != "PEDIGREE")
    {
      errors << priority(warning) << "Misformed pedigree block in model '"
             << name << "' will be ignored."
             << endl;
      continue;
    }

    AttrList::const_iterator iter = (*i)->attrs()->find("VALUE");
    if(iter == (*i)->attrs()->end())
      continue;

    string pedigree = iter->second.String();

    // For each individual

    LSFList::const_iterator j = (*i)->List()->begin();
    for( ; j != (*i)->List()->end(); ++j)
    {
      if(!(*j) || !(*j)->List() || toUpper((*j)->name()) != "INDIVIDUAL")
      {
        errors << priority(warning) << "Misformed individual block in model '"
               << name << "' will be ignored."
               << endl;
        continue;
      }
      
      iter = (*j)->attrs()->find("VALUE");
      if(iter == (*j)->attrs()->end())
        continue;

      string ind = iter->second.String();
      
      string phen = pedigree + "~" + ind;

      pm.add_phenotype(phen);
      
      uint pid = pm.get_phenotype_id(phen);
      
      // For each genotype
      
      LSFList::const_iterator k = (*j)->List()->begin();
      for( ; k != (*j)->List()->end(); ++k)
      {
        if(!(*k) || toUpper((*k)->name()) != "GENOTYPE") continue;

        iter = (*k)->attrs()->find("VALUE");
        if(iter == (*k)->attrs()->end())
        {
          errors << priority(warning) << "Misformed genotype in model '"
                 << name << "' will be ignored."
                 << endl;
          continue;
        }

        string genotype = iter->second.String();

        size_t pos = genotype.find('/');
        
        if(pos >= genotype.size())
        {
          errors << priority(warning) << "Misformed genotype in model '"
                 << name << "' will be ignored."
                 << endl;
          continue;
        }
        
        genotype[pos] = separator;
        
        // Find the genotype in the model
        unphased_genotype gtype = pm.get_unphased_genotype(genotype);

        if(!gtype.is_valid())           // We can't find the genotype.
        {
          errors << priority(warning) << "Unknown genotype " << genotype
                 << " in model " << name << ".  Genotype "
                 << "will be ignored." << endl;
          continue;
        }

        AttrList::const_iterator iter = (*k)->attrs()->find("PENETRANCE");
        if(iter == (*k)->attrs()->end())
        {
          errors << priority(warning) << "Malformed genotype (no penetrance value)."
                 << "  Please check your file." << endl;
                 
          continue;
        }
        
        double d = iter->second.Real();

        if(!finite(d) || d < 0.0 || d > 1.0)
        {
          errors << priority(warning) << "Malformed penetrance (penetrance value out of range).  "
                 << "Please check your file." << endl;
                 
          continue;
        }

        if(d > 0.0) 
          pm.add_unphased_penetrance(d, pid, gtype.get_id());
      } 
    }
  }

  return true;
}

} // End namespace MLOCUS
} // End namespace SAGE

