#include <string>
#include <functional>
#include <fstream>
#include <iomanip>
#include <cassert>
#include "LSF/LSFinit.h"
#include "LSF/LSFfile.h"
#include "app/SAGEapp.h"
#include "rped/rped.h"
#include "rped/rpfile.h"
#include "pedinfo/stats.h"
#include "pedinfo/stats_view.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "pairs/relmatrix.h"

using namespace std;

void print_title(ostream &o)
{
  o << "S.A.G.E. " << SAGE::APP::sage_release << " -- TESTRELMATRIX " << endl
    << endl;
}

int main(int argc, char* argv[])
{
  print_title(std::cout);
  //SAGEapp::expire();

  if (argc != 3)
  {
    cerr << "usage: " << argv[0] << " <parameters> <pedigree>"
         << endl << endl
         << "Command line parameters:"             << endl
         << "  parameters   - Parameter File"  << endl
         << "  pedigree     - Pedigree Data File"       << endl
         << endl << endl;
    exit(EXIT_FAILURE);
  }

  LSFInit();

  SAGE::sage_cerr << SAGE::prefix("%%TESTRELMATRIX-%P: ");

  std::cout << "Loading parameters................";
  LSFBase *params = loadLSFfile(argv[1], "TESTRELMATRIX Parameter file", SAGE::sage_cerr, false);
  std::cout << "done." << endl;

  if(!params)
  {
    SAGE::sage_cerr << SAGE::priority(SAGE::fatal) << "Error reading parameter file.... aborting." << endl;
    exit(EXIT_FAILURE);
  }

  ofstream info_file;
  info_file.open("test.inf");

  if(!info_file)
  {
    SAGE::sage_cerr << SAGE::priority(SAGE::fatal) 
              << "Cannot open output file: test.inf.  Exiting..." << endl;
    exit(EXIT_FAILURE);
  }

  ofstream out_file;
  out_file.open("test.out");

  if(!out_file)
  {
    SAGE::sage_cerr << SAGE::priority(SAGE::fatal) 
              << "Cannot open output file: test.out.  Exiting..." << endl;
    exit(EXIT_FAILURE);
  }

  //print_title(out_file);
  //print_title(info_file);

  SAGE::cerrormultistream errors;
  SAGE::cerrorstream error_file(info_file);
  error_file.prefix("%%TESTRELMATRIX-%P:");

  errors.insert   (SAGE::sage_cerr);
  errors.restrict (SAGE::r_ge, SAGE::error);
  errors.insert   (error_file);
  errors.restrict (SAGE::r_ge, SAGE::information);

  SAGE::RPED::RefMultiPedigree p;

  bool pedigree_loaded = false;

  LSFList::const_iterator i;
  AttrVal a;
  for(i=params->List()->begin(); i!=params->List()->end(); ++i)
  {
    if( !*i ) continue;

    if( toUpper( (*i)->name() ) == "PEDIGREE" )
    {
      if( (*i)->attrs() && (*i)->attrs()->has_attr("column") )
      {
        SAGE::RPED::RefLSFFortranPedigreeFile ped_reader;
        ped_reader.set_force_skip_markers(true);    
        ped_reader.set_force_skip_traits(false);    
        ped_reader.set_force_dynamic_markers(false);
        ped_reader.process_parameters(p.info(), *i);

        std::cout << "Reading pedigrees................." << flush;
        if( !ped_reader.input(p, argv[2], info_file) )
        {
          errors << SAGE::priority(SAGE::fatal) << "Error reading pedigree data" << endl;
          exit(EXIT_FAILURE);
        }
        pedigree_loaded = true;

        ped_reader.print_mped(p, argv[2], info_file);

        break;
      }
      else
      {
        SAGE::RPED::RefLSFDelimitedPedigreeFile ped_reader;
        ped_reader.set_force_skip_markers(true);    
        ped_reader.set_force_skip_traits(false);    
        ped_reader.set_force_dynamic_markers(false);
        ped_reader.process_parameters(p.info(), *i);

        std::cout << "Reading pedigrees................." << flush;
        if( !ped_reader.input(p,argv[2], info_file) )
        {
          errors << SAGE::priority(SAGE::fatal) << "Error reading pedigree data" << endl;
          exit(EXIT_FAILURE);
        }
        pedigree_loaded = true;

        ped_reader.print_mped(p, argv[2], info_file);

        break;
      }
    }
  }  

  if( !pedigree_loaded )
  {
    errors << SAGE::priority(SAGE::fatal) << "Fatal Error: No pedigree specified!  Terminating..." << endl;
    exit(EXIT_FAILURE);
  }

  std::cout << "done." << endl;

  std::cout << "Building pedigrees................" << flush;
  p.build();
  std::cout << "done." << endl;

  if( !p.pedigree_count() )
  {
    errors << SAGE::priority(SAGE::critical) << "No pedigrees to analyze... aborting." << endl;
    exit(EXIT_FAILURE);
  }

  size_t error_count = 0;
  SAGE::RPED::RefMultiPedigree::pedigree_const_iterator j;
  for( j = p.pedigree_begin(); j != p.pedigree_end(); ++j)
    error_count += j->error_count();

  if( error_count )
    errors << SAGE::priority(SAGE::error) << "Errors appear in pedigree data.  " 
           << "Results may be incomplete." << endl;

  // Sort each pedigree.
  std::cout << "Sorting pedigrees................." << flush;

  SAGE::RPED::RefMultiPedigree::pedigree_iterator                   ped;

  for( ped = p.pedigree_begin(); ped != p.pedigree_end(); ++ped )
    PedigreeSort(*ped);

  std::cout << "done." << endl;

  // Fill in relation type between two individuals into relation matrix.
  std::cout << "Generating relation matrix........" << flush;
  
  SAGE::RelationMatrix rel_matrix(&p);

  std::cout << "done." << endl;
 
  std::cout << "Writing relation table............" << flush;

  out_file  << std::endl << "Processed relation matrices:" << std::endl;
  rel_matrix.view(out_file, p, false);

  out_file  << std::endl << "Raw relation matrices:" << std::endl;
  rel_matrix.view(out_file, p, true);

  out_file  << std::endl << "Relative pairs:" << std::endl;
  rel_matrix.view_pairs(out_file, p);

  out_file  << std::endl << "Relative pairs w/ gender:" << std::endl;
  rel_matrix.view_pairs_gender(out_file, p);

  std::cout << "done." << endl;

  std::cout << endl << "Analysis complete!" << endl;  
  exit(EXIT_SUCCESS);

  return 0;
}
