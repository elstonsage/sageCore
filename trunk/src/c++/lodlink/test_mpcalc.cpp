//============================================================================
// File:      test_mpcalc.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   Created 12/3/2.                                                   
//                                                                          
// Notes:     Tests the likelihood calculator classes.
//                                                                          
// Copyright (c) 2000 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================
#include <string>
#include <cstdlib>
#include <cmath>
#include <functional>
#include <fstream>
#include <iomanip>
#include <cassert>
#include "LSF/LSFinit.h"
#include "LSF/LSFfile.h"
#include "app/SAGEapp.h"
#include "rped/rpfile.h"
#include "fped/fped.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "mlocus/mfile.h"
#include "mlocus/imodel.h"
#include "lodlink/likelihood.h"
#include "util/StringUtils.h"


using namespace std;
using namespace SAGE;
using namespace LODLINK;

int main(int argc, char* argv[])
{

  // ============  The preliminaries.  =============
  //
  if (argc < 6 || 7 < argc)
  {
    cerr << "usage: " << argv[0] << " <parameters> <pedigree> <locus> <male_theta> <female_theta>\n\n"
         << "Command line parameters:\n"     
         << "  parameters    - Parameter File\n"  
         << "  pedigree      - Pedigree Data File\n"  
         << "  locus         - Marker Locus Description File\n"
         << "  male_theta    - Male Recombination Fraction\n"
         << "  female_theta  - Female Recombination Fraction\n\n" 
         << "  alpha         - Proportion of families w. linkage\n\n" 
         << endl;
    exit(EXIT_FAILURE);
  }

  LSFInit();

  sage_cerr << prefix("%%TESTMPCALC-%P: ");

  // - Read parameter file.
  //
  LSFBase *params = loadLSFfile(argv[1], "TESTMPCALC Parameter file", sage_cerr, false);

  if(!params)
  {
    sage_cerr << priority(fatal) << "Error reading parameter file.... aborting." << endl;
    exit(EXIT_FAILURE);
  }

  
  ofstream info_file;
  info_file.open("testmpcalc.inf");

  if(!info_file)
  {
    sage_cerr << priority(fatal) 
              << "Cannot open output file: pedinfo.inf.  Exiting..." << endl;
    exit(EXIT_FAILURE);
  }

  ofstream out_file;
  out_file.open("testmpcalc.out");

  if(!out_file)
  {
    sage_cerr << priority(fatal) 
              << "Cannot open output file: testpeeler.out.  Exiting..." << endl;
    exit(EXIT_FAILURE);
  }

  SAGE::cerrormultistream errors;
  cerrorstream error_file(info_file);
  error_file.prefix("%%TESTMPCALC-%P:");
  errors.insert(sage_cerr);
  errors.restrict(r_ge, error);
  errors.insert(error_file);
  errors.restrict(r_ge, information);
  
  RPED::RefMultiPedigree p;

  // - Read marker locus description file.
  //
  MLOCUS::InheritanceModelFile   imf(errors);
  imf.input(p.info().markers(), argv[3]);
  

  // - Read pedigree data.
  //
  bool pedigree_loaded = false;

  LSFList::const_iterator i;
  AttrVal a;
  for(i=params->List()->begin(); i!=params->List()->end(); ++i)
  {
    if(! *i) continue;

    if( toUpper( (*i)->name() ) == "PEDIGREE" )
    {
      if( (*i)->attrs() && (*i)->attrs()->has_attr("column") )
      {
        RPED::RefLSFFortranPedigreeFile ped_reader(errors);
        ped_reader.set_force_skip_markers(false);
        ped_reader.set_force_skip_traits(false);
        ped_reader.set_force_dynamic_markers(true);
        ped_reader.process_parameters(p.info(), *i);
        
        if( !ped_reader.input(p, argv[2], info_file) )
        {
          errors << priority(fatal) << "Error reading pedigree data" << endl;
          exit(EXIT_FAILURE);
        }
        pedigree_loaded = true;
        
        // Need to be added to print 10 individuals. - yjs 
        //
        ped_reader.print_mped(p, argv[2], info_file);
        
        break;
      }
      else
      {
        RPED::RefLSFDelimitedPedigreeFile ped_reader(errors);
        ped_reader.set_force_skip_markers(false);
        ped_reader.set_force_skip_traits(false);
        ped_reader.set_force_dynamic_markers(true);
        ped_reader.process_parameters(p.info(), *i);

        if( !ped_reader.input(p,argv[2], info_file) )
        {
          errors << priority(fatal) << "Error reading pedigree data" << endl;
          exit(EXIT_FAILURE);
        }
        pedigree_loaded = true;
        
        // Need to be added to print 10 individuals. - yjs 
        //
        ped_reader.print_mped(p, argv[2], info_file);
        
        break;
      }
    }
  }  

  if( !pedigree_loaded )
  {
    errors << priority(fatal) << "Fatal Error: No pedigree specified!  Terminating..." << endl;
    exit(EXIT_FAILURE);
  }

  if( !p.pedigree_count() )
  {
    errors << priority(critical) << "No pedigrees to analyze... aborting." << endl;
    exit(EXIT_FAILURE);
  }

  size_t error_count = 0;
  RPED::RefMultiPedigree::pedigree_const_iterator j;
  for( j = p.pedigree_begin(); j != p.pedigree_end(); ++j)
    error_count += j->error_count();

  if( error_count )
    errors << priority(error) << "Errors appear in pedigree data.  " 
           << "Results may be incomplete." << endl;

  RPED::RefMultiPedigree::pedigree_iterator jj;
  for( jj = p.pedigree_begin(); jj != p.pedigree_end(); ++jj)
  {
    PedigreeSort( *jj );
  }
  
  FPED::FilteredMultipedigree  fp(p);
  FPED::MPFilterer::add_multipedigree_filtered_by_members(fp, p, FPED::always_keep());
  fp.construct();
  
  // ===============  Test of likelihood calculators  ================
  //
  typedef FPED::FilteredMultipedigree::pedigree_const_iterator     pedigree_const_iterator;
  typedef FPED::FilteredMultipedigree::subpedigree_const_iterator  subpedigree_const_iterator;
  typedef FPED::FilteredMultipedigree::member_const_iterator       member_const_iterator;
  
  const FPED::FilteredMultipedigreeInfo&  mp_info = fp.info();
  size_t  trait  = mp_info.marker_find("M1");
  size_t  marker = mp_info.marker_find("M2");
  assert(trait != (size_t)(-1));
  assert(marker != (size_t)(-1));
  
  double  male_theta   = atof(argv[4]);
  double  female_theta = atof(argv[5]);
  double  alpha = false;
  bool  sex_specific  = false;
  bool  use_alpha     = false;
  
  if(argc == 7)
  {
    use_alpha = true;
    alpha = atof(argv[6]);
  }
  
  if(male_theta != female_theta)
  {
    sex_specific = true;
  }
  
  mle_sub_model  mle(sex_specific, use_alpha);

  if(sex_specific)
  {
    mle.set_male_theta(male_theta);
    mle.set_female_theta(female_theta);
  }
  else
  {
    mle.set_average_theta(male_theta);
  }
  
  if(use_alpha)
  {
    mle.set_alpha(alpha);
  }
  
  mped_calculator  mp_calc(fp, mle, trait, marker);
          
  log_double  like = mp_calc.likelihood();
  log_double  unlinked_like = mp_calc.unlinked_likelihood();

  out_file 
    << (OUTPUT::Table()
    <<   (OUTPUT::TableRow() << "likelihood" << like.get_double())
    <<   (OUTPUT::TableRow() << "unlinked likelihood" << unlinked_like.get_double())
    <<   (OUTPUT::TableRow() << "log likelihood" << log10(like.get_double()) << "lod score" << (log10(like.get_double()) - log10(unlinked_like.get_double()))));
}



