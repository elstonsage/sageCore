//============================================================================
// File:      testrelpair_simple.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   Created 4/5/01
//                                                                          
// Notes:     Find and print pairs using pair_generator class.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include <string>
#include <functional>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <ios>
#include <stdio.h>
#include "LSF/LSFinit.h"
#include "LSF/LSFfile.h"
#include "app/SAGEapp.h"
#include "rped/rped.h"
#include "rped/rpfile.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "numerics/sinfo.h"
#include "pairs/relpair.h"


using namespace std;
using namespace SAGE;

string trait_t_to_string(RPED::RefTraitInfo::trait_t tt)
{
  switch(tt)
  { 
    case RPED::RefTraitInfo::invalid_trait:     return "invalid";
    case RPED::RefTraitInfo::continuous_trait:  return "continuous";
    case RPED::RefTraitInfo::binary_trait:      return "binary";
    case RPED::RefTraitInfo::categorical_trait: return "categorical";
    case RPED::RefTraitInfo::discrete_trait:    return "discrete";
    default:                                    return "unknown";
  }
}

void print_title(ostream &o)
{
  o << "S.A.G.E. " 
    << APP::sage_release 
    << " -- TESTRELPAIR_SIMPLE " 
    << endl
    << endl;
}

// - Print a relative pair.
//
void
print_pair(std::ostream& out, pair_generator::relative_pair pr) 
{
  out << std::setw(9) << left << pair_generator::pair_type_to_string(pr.type()) << "   ";
  out << "ped" << pr.member_one()->pedigree()->name() << " ";
  out << "member1    " << std::setw(3) << std::right << pr.member_one()->name() << "   ";
  out << "member2    " << std::setw(3) << std::right << pr.member_two()->name() << "   ";
  
  if(pr.connector_one())
  {
    out << "connector1 " << std::setw(3) << std::right << pr.connector_one()->name() << "   ";
  }
  
  if(pr.connector_two())
  {
    out << "connector2 " << std::setw(3) << std::right << pr.connector_two()->name() << "   ";
  }
  
  out << endl;
}

int main(int argc, char* argv[])
{
  // - HOUSEKEEPING
  //
  //print_title(cout);
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

  sage_cerr << prefix("%%TESTRELPAIR_SIMPLE-%P: ");

  cout << "Loading parameters..............";
  LSFBase *params = loadLSFfile(argv[1], "TESTRELPAIR_SIMPLE Parameter file", sage_cerr, false);
  cout << "done." << endl;

  if(!params)
  {
    sage_cerr << priority(fatal) << "Error reading parameter file.... aborting." << endl;
    exit(EXIT_FAILURE);
  }

  cout << "About to allocate multi-pedigree" << endl;
  RPED::RefMultiPedigree p;

  bool pedigree_loaded = false;

  cout << "About to read pedigree" << endl;

  LSFList::const_iterator i;
  AttrVal a;
  for(i=params->List()->begin(); i!=params->List()->end(); ++i)
  {
    if( !*i ) continue;

    if( toUpper( (*i)->name() ) == "PEDIGREE" )
    {
      if( (*i)->attrs() && (*i)->attrs()->has_attr("column") )
      {
        RPED::RefLSFFortranPedigreeFile ped_reader;
        ped_reader.set_force_skip_markers(true);
        ped_reader.set_force_skip_traits(false);
        ped_reader.set_force_dynamic_markers(false);
        ped_reader.process_parameters(p.info(), *i);

        cout << "Reading pedigrees..............." << flush;
        if( !ped_reader.input(p, argv[2], cout) )
        {
          cerr << "Error reading pedigree data" << endl;
          exit(EXIT_FAILURE);
        }
        pedigree_loaded = true;
        break;
      }
      else
      {
        RPED::RefLSFDelimitedPedigreeFile ped_reader;
        ped_reader.set_force_skip_markers(true);
        ped_reader.set_force_skip_traits(false);
        ped_reader.set_force_dynamic_markers(false);
        ped_reader.process_parameters(p.info(), *i);

        cout << "Reading pedigrees..............." << flush;
        if( !ped_reader.input(p,argv[2], cout) )
        {
          cerr << "Error reading pedigree data" << endl;
          exit(EXIT_FAILURE);
        }
        pedigree_loaded = true;
        break;
      }
    }
  }  

  if( !pedigree_loaded )
  {
    cerr << "Fatal Error: No pedigree specified!  Terminating..." << endl;
    exit(EXIT_FAILURE);
  }

  cout << "done." << endl;

  if( !p.pedigree_count() )
  {
    cerr << "No pedigrees to analyze... aborting." << endl;
    exit(EXIT_FAILURE);
  }

  size_t error_count = 0;
  RPED::RefMultiPedigree::pedigree_const_iterator j;
  for( j = p.pedigree_begin(); j != p.pedigree_end(); ++j)
    error_count += j->error_count();

  if( error_count )
    cerr << "Errors appear in pedigree data.  " 
           << "Results may be incomplete." << endl;

  RPED::RefMultiPedigree::pedigree_iterator jj;
  for( jj = p.pedigree_begin(); jj != p.pedigree_end(); ++jj)
  {
//  cout << "    Sorting pedigree: " << jj->name() << endl;
    PedigreeSort( *jj );
  }

  // - Create output file.
  //
  ofstream  out_file;
  out_file.open("testrelpair_simple.out");


  // - TEST PAIR GENERATOR. 
  //  
  // - Generate and print all pairs of each type.
  //  
  typedef pair_generator::pair_type   type;
  
  RPED::RefPedigree*  refped;    
  pair_generator  pg;      
  RPED::RefMultiPedigree::pedigree_iterator  ii;
  
  // - For each pedigree print pairs of all defined types.
  //
  for(ii = p.pedigree_begin(); ii != p.pedigree_end(); ++ii)
  {
    refped = &(*(ii));
    pg.set_pedigree(refped);
    out_file  << "\n#########  PEDIGREE: " << refped->name() << "  #########" << endl;
  
    // - For each pair type.
    //
    for(unsigned int i = 1; i < pair_generator::NULL_TYPE; ++i)
    {
      out_file << "\n";
    
      // - For each pair.
      //
      pg.set_type(static_cast<pair_generator::pair_type>(i));
      for(pair_generator::iterator iter = pg.begin(); iter != pg.end(); ++iter)
      {
        print_pair(out_file, *iter);
      }
    }
  }
  
  out_file << "\n\n" << endl;
  
  // - For each pedigree print all possible pairs w/o regard to type.
  //
  for(ii = p.pedigree_begin(); ii != p.pedigree_end(); ++ii)
  {
    refped = &(*(ii));
    pg.set_pedigree(refped);
    out_file  << "\n#########  PEDIGREE: " << refped->name() << "  #########" << endl;
  
    out_file << "\n";
  
    // - For each pair.
    //
    pg.set_type(pair_generator::EVERY);
    for(pair_generator::iterator iter = pg.begin(); iter != pg.end(); ++iter)
    {
      print_pair(out_file, *iter);
    }
  }
  
  
    
  exit(EXIT_SUCCESS);
}


