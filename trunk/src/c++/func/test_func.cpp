#include <string>
#include <functional>
#include <fstream>
#include <iomanip>
#include <cassert>
#include "LSF/LSFinit.h"
#include "LSF/LSFfile.h"
#include "app/SAGEapp.h"
#include "rped/rpfile.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "mlocus/mfile.h"
#include "mlocus/imodel.h"
#include "output/Output.h"
#include "func/evalfunc.h"

using namespace std;
using namespace SAGE;

int main(int argc, char* argv[])
{
  if (argc != 4)
  {
    cerr << "usage: " << argv[0] << " <parameters> <pedigree> <locus>" << endl 
         <<                                                               endl
         << "Command line parameters:"                                 << endl
         << "  parameters   - Parameter File"                          << endl
         << "  pedigree     - Pedigree Data File"                      << endl
         << "  locus        - Marker Locus Description File"           << endl
         <<                                                               endl 
         <<                                                               endl;

    exit(EXIT_FAILURE);
  }

  LSFInit();

  sage_cerr << prefix("%%TESTFUNC-%P: ");

  // - Read parameter file.
  //
  LSFBase *params = loadLSFfile(argv[1], "TESTFUNC Parameter file", sage_cerr, false);

  if(!params)
  {
    sage_cerr << priority(fatal) << "Error reading parameter file.... aborting." << endl;
    exit(EXIT_FAILURE);
  }

  ofstream info_file;
  info_file.open("testfunc.inf");

  if(!info_file)
  {
    sage_cerr << priority(fatal) 
              << "Cannot open output file: pedinfo.inf.  Exiting..." << endl;
    exit(EXIT_FAILURE);
  }

  ofstream out_file;
  out_file.open("testfunc.out");

  if(!out_file)
  {
    sage_cerr << priority(fatal) 
              << "Cannot open output file: testfunc.out.  Exiting..." << endl;
    exit(EXIT_FAILURE);
  }

  //print_title(info_file);
  //print_title(out_file);

  SAGE::cerrormultistream errors;
  cerrorstream error_file(info_file);

  error_file.prefix("%%TESTFUNC-%P:");

  errors.insert   (sage_cerr);
  errors.restrict (r_ge, error);
  errors.insert   (error_file);
  errors.restrict (r_ge, information);
  
  RPED::RefMultiPedigree p;

  // Read marker locus description file.

  MLOCUS::InheritanceModelFile   imf(errors);

  imf.input(p.info().markers(), argv[3]);

  bool pedigree_loaded = false;

  AttrVal a;

  for(LSFList::const_iterator i = params->List()->begin(); i != params->List()->end(); ++i)
  {
    if( !*i ) continue;

    if( toUpper( (*i)->name() ) == "PEDIGREE" )
    {
      if( (*i)->attrs() && (*i)->attrs()->has_attr("column") )
      {
        RPED::RefLSFFortranPedigreeFile ped_reader(errors);

        ped_reader.set_force_skip_markers    (false);
        ped_reader.set_force_skip_traits     (false);
        ped_reader.set_force_dynamic_markers (true);

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

        ped_reader.set_force_skip_markers    (false);
        ped_reader.set_force_skip_traits     (false);
        ped_reader.set_force_dynamic_markers (true);

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

  for(RPED::RefMultiPedigree::pedigree_const_iterator j = p.pedigree_begin(); j != p.pedigree_end(); ++j)
    error_count += j->error_count();

  if(error_count)
    errors << priority(error) << "Errors appear in pedigree data.  " << "Results may be incomplete." << endl;

  for(RPED::RefMultiPedigree::pedigree_iterator jj = p.pedigree_begin(); jj != p.pedigree_end(); ++jj)
    PedigreeSort(*jj);
  
  // Create traits specified by function blocks.

  SAGE::FUNC::FunctionEvaluator::evaluateFunction(p, errors, params);
  
  // Print trait values for each member in multipedigree.

  SAGE::OUTPUT::Table t;

  for(RPED::RefMultiPedigree::pedigree_iterator p_iter = p.pedigree_begin(); p_iter != p.pedigree_end(); ++p_iter)
  {
    for(RPED::RefPedigree::member_iterator m_iter = p_iter->member_begin(); m_iter != p_iter->member_end(); ++m_iter)
    {
      SAGE::OUTPUT::TableRow r;

      for(size_t i = 0; i < p_iter->info().trait_count(); ++i)
        r << p_iter->info().trait(m_iter->index(), i);
      
      for(size_t i = 0; i < p_iter->info().marker_count(); ++i)
        r << p_iter->info().phenotype(m_iter->index(), i);
        
      t << r;
    }
  }
  
  out_file << t;
}



