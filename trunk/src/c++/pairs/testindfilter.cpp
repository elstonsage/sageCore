//============================================================================
// File:      testindfilter.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                         
// History:   Created 10/20/00.
//                                                                          
// Notes:     Written to test SAGE::ind_filter_trait and SAGE::ind_filter classes.
//                                                                          
// Copyright (c) 2000 R.C. Elston                                           
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
    case RPED::RefTraitInfo::discrete_trait:    return "discrete";
    case RPED::RefTraitInfo::categorical_trait: return "categorical";
    default:                                    return "unknown";
  }
}

void print_title(ostream &o)
{
  o << "S.A.G.E. " << APP::sage_release << " -- TESTINDFILTER " << endl
    << endl;
}

int main(int argc, char* argv[])
{
  print_title(cout);
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

  sage_cerr << prefix("%%TESTINDFILTER-%P: ");

  cout << "Loading parameters..............";
  LSFBase *params = loadLSFfile(argv[1], "TESTINDFILTER Parameter file", sage_cerr, false);
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


  // - Read pedigree file.
  //
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



  // - Test ind_filter_trait and ind_filter. 
  //  
  RPED::RefTraitInfo                            rti;
  //RPED::RefTraitInfo::trait_t                   rt_type;
  RPED::RefMPedInfo                             rmpi = p.info();
  size_t                                  no_of_traits = rmpi.trait_count();
  
  RPED::RefMultiPedigree::pedigree_iterator     ii;
  RPED::RefPedigree*                            refped;
  RPED::RefPedigree::member_iterator            iter;
  ind_filter_trait                        filter_element;
  ind_filter                              filter;
  
  for(ii = p.pedigree_begin(); ii != p.pedigree_end(); ++ii)
  {
    refped = &(*(ii));
    cout << endl << "#########  PEDIGREE: " << refped->name() << "  #########" << endl;
    
    for(size_t t = 0; t < no_of_traits; ++t)
    {
      size_t total_count        = 0;
      size_t informative_count  = 0;
      size_t in_range_count     = 0;
      size_t affected_count     = 0;  
      size_t unaffected_count   = 0;
      
      filter_element.set_trait(t);
      if(t == 1)
      {
        filter_element.set_affected_min(40.0);
        // filter_element.set_min(30.0);
      }
      filter.clear_traits();
      filter.set_trait(filter_element);
      
      for(iter = refped->member_begin(); iter != refped->member_end(); ++iter)
      {
        ++total_count;
        if(filter.informative(&(*iter)))  ++informative_count;
        if(filter.in_range(&(*iter)))     ++in_range_count;
        if(filter.affected(&(*iter)))     ++affected_count;
        if(filter.unaffected(&(*iter)))   ++unaffected_count;
      }
      
      cout << endl;
      cout << "trait " << t << endl;
      cout << "informative "   << informative_count
           << "  in range "    << in_range_count
           << "  affected "    << affected_count
           << "  unaffected "  << unaffected_count << endl;
    } 
  }
  exit(EXIT_SUCCESS);
}


