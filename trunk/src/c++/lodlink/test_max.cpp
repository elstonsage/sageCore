//============================================================================
// File:      test_max.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   Created 12/9/2.                                                   
//                                                                          
// Notes:     Tests likelihood maximization.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
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
#include "maxfun/sub_model.h"
#include "lodlink/likelihood.h"
#include "lodlink/max_opt.h"

using namespace std;
using namespace SAGE;
using namespace LODLINK;

int main(int argc, char* argv[])
{

  // ============  The preliminaries.  =============
  //
  if (argc != 5)
  {
    cerr << "usage: " << argv[0] << " <parameters> <pedigree> <locus> <sex_specific>\n\n"
         << "Command line parameters:\n"     
         << "  parameters   - Parameter File\n"  
         << "  pedigree     - Pedigree Data File\n"  
         << "  locus        - Marker Locus Description File\n"
         << "  sex specific - Sex Specific Recombination Fractions (true | false).\n"
         << endl;
    exit(EXIT_FAILURE);
  }

  LSFInit();

  sage_cerr << prefix("%%TESTMAX-%P: ");

  // - Read parameter file.
  //
  LSFBase *params = loadLSFfile(argv[1], "TESTMAX Parameter file", sage_cerr, false);

  if(!params)
  {
    sage_cerr << priority(fatal) << "Error reading parameter file.... aborting." << endl;
    exit(EXIT_FAILURE);
  }

  
  ofstream info_file;
  info_file.open("testmax.inf");

  if(!info_file)
  {
    sage_cerr << priority(fatal) 
              << "Cannot open output file: pedinfo.inf.  Exiting..." << endl;
    exit(EXIT_FAILURE);
  }

  ofstream out_file;
  out_file.open("testmax.out");

  if(!out_file)
  {
    sage_cerr << priority(fatal) 
              << "Cannot open output file: testpeeler.out.  Exiting..." << endl;
    exit(EXIT_FAILURE);
  }

  SAGE::cerrormultistream errors;
  cerrorstream error_file(info_file);
  error_file.prefix("%%TESTMAX-%P:");
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
  
  // - Default is 'true'.
  //
  bool  sex_specific = (string(argv[4]) == "true");
  mle_sub_model  mle(sex_specific, false);
  
  mle_sub_model  unlinked_mle;
  unlinked_mle.set_average_theta(.5);
  
  mped_calculator  mp_calc(fp, mle, trait, marker);
  MAXFUN::ParameterMgr    param_mgr;
  MAXFUN::DebugCfg        debug_cfg;
  MAXFUN::SequenceCfg     sequence_cfg;
  
  int  failure = param_mgr.addSubModel(&mle);
  assert(! failure);
  
  set_maxfun_configuration(sequence_cfg);
  
  MAXFUN::Results  data = MAXFUN::maximize(param_mgr, mp_calc, sequence_cfg, debug_cfg);
  
  out_file << "last maxfun return code  " << data.getExitFlag() << endl;
  out_file << "final values of mle_sub_model:" << endl;
  out_file << mle << "\n" << endl;  
  
  double  result = data.getFinalFunctionValue();   // Natural log of likelihood
  
  mped_calculator  unlinked_mp_calc(fp, unlinked_mle, trait, marker);
  log_double  unlinked_like = unlinked_mp_calc.likelihood();  // Likelihood.
  
          
  out_file << "maximum log likelihood " << (result / log(10.0)) << ", " 
           << "maximum lod score  " << (result - unlinked_like.get_log()) / log(10.0)
           << endl;
}



