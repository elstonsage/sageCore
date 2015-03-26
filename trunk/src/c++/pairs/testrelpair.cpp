//============================================================================
// File:      testrelpair.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:                                                      
//                                                                          
// Notes:     Rigorously tests pair_generator class for self consistency.
//            Also tests filtering_pair_generator.
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

namespace SAGE {

string trait_t_to_string(RPED::RefTraitInfo::trait_t tt)
{
  switch(tt)
  {
    case RPED::RefTraitInfo::invalid_trait:
      return "invalid";
    case RPED::RefTraitInfo::continuous_trait:
      return "continuous";
    case RPED::RefTraitInfo::binary_trait:
      return "binary";
    case RPED::RefTraitInfo::categorical_trait:
      return "categorical";
    case RPED::RefTraitInfo::discrete_trait:
      return "discrete";
    default:
      return "unknown";
  }
}

void print_title(ostream &o)
{
  o << "S.A.G.E. " << APP::sage_release << " -- TESTRELPAIR " << endl
    << endl;
}

} // End namespace SAGE

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

  SAGE::sage_cerr << SAGE::prefix("%%TESTRELPAIR-%P: ");

  cout << "Loading parameters..............";
  LSFBase *params = loadLSFfile(argv[1], "TESTRELPAIR Parameter file", SAGE::sage_cerr, false);
  cout << "done." << endl;

  if(!params)
  {
    SAGE::sage_cerr << SAGE::priority(SAGE::fatal) << "Error reading parameter file.... aborting." << endl;
    exit(EXIT_FAILURE);
  }

  cout << "About to allocate multi-pedigree" << endl;
  SAGE::RPED::RefMultiPedigree p;

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
        SAGE::RPED::RefLSFFortranPedigreeFile ped_reader;
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
        
        // Need to be added to print 10 individuals. - yjs 
        //
        ped_reader.print_mped(p, "");
        
        break;
      }
      else
      {
        SAGE::RPED::RefLSFDelimitedPedigreeFile ped_reader;
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
        
        // Need to be added to print 10 individuals. - yjs 
        //
        ped_reader.print_mped(p, "");

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
  SAGE::RPED::RefMultiPedigree::pedigree_const_iterator j;
  for( j = p.pedigree_begin(); j != p.pedigree_end(); ++j)
    error_count += j->error_count();

  if( error_count )
    cerr << "Errors appear in pedigree data.  " 
           << "Results may be incomplete." << endl;

  SAGE::RPED::RefMultiPedigree::pedigree_iterator jj;
  for( jj = p.pedigree_begin(); jj != p.pedigree_end(); ++jj)
  {
//  cout << "    Sorting pedigree: " << jj->name() << endl;
    PedigreeSort( *jj );
  }


  // - TEST PAIR GENERATOR. 
  //  
  // - Generate a count for each pair type.
  //  
  typedef SAGE::pair_generator::pair_type   type;
  
  size_t  totals[9];              // Pair type totals for a pedigree.
  size_t  multi_totals[9];        // Pair type totals for a multipedigree.
  
  for(int i = 0; i < 9; ++i)
  {
    multi_totals[i] = 0;
  }
  
  size_t  total;              // Single pedigree, single pair type.
  size_t  grand_total;        // Single pedigree, sum of all pair types calculated individually.
  size_t  check_total;        // Single pedigree, sume of all pair types calculated together.
  
  SAGE::RPED::RefPedigree*  refped;    
  SAGE::pair_generator  pg;      
  SAGE::RPED::RefMultiPedigree::pedigree_iterator  ii;
  
  // - For each pedigree.
  //
  for(ii = p.pedigree_begin(); ii != p.pedigree_end(); ++ii)
  {
    refped = &(*(ii));
    pg.set_pedigree(refped);
    cout << endl << "#########  PEDIGREE: " << refped->name() << "  #########" << endl;
    
  
    // - For each pair type.
    //
    for(unsigned int i = 1; i < SAGE::pair_generator::NULL_TYPE; ++i)
    {
    
      // - Count pairs of current type.
      //
      total = 0;
      pg.set_type(static_cast<SAGE::pair_generator::pair_type>(i));
      for(SAGE::pair_generator::iterator iter = pg.begin(); iter != pg.end(); ++iter)
      {
        ++total;
      }
      
      totals[i - 1] = total;
      multi_totals[i - 1] += total;
    }
    
    
    // - Calculate grand total.
    //  
    grand_total = 0;
    for(int i = SAGE::pair_generator::PARENTAL - 1; i < SAGE::pair_generator::NULL_TYPE - 1; ++i)
    {
      grand_total += totals[i];
    }
    
    // - Generate a total by iterating thru all pair types at once.
    //  
    check_total = 0;

    pg.set_types(pg.ALL_TYPES);
    for(SAGE::pair_generator::iterator iter = pg.begin(); iter != pg.end(); ++iter)
    {
      ++check_total;

//#define PG_VERBOSE
#ifdef PG_VERBOSE
        cout << *iter;
#endif
//#undef PG_VERBOSE
      
    }
    
    cout << endl << endl;
    
    cout << "        TOTALS FOR EACH PAIR TYPE" << endl << endl;
    cout << setw(15) << left << "Pair Type" << setw(25) << right << "Count" << endl;
    cout << setw(15) << left << "---------" << setw(25) << right << "--------" << endl;
    
    for(int i = SAGE::pair_generator::PARENTAL - 1; i < SAGE::pair_generator::NULL_TYPE - 1; ++i)
    {
      cout << setw(15) << left << SAGE::pair_generator::pair_type_to_string((type)(i + 1)) << setw(25) << right << totals[i] << endl;
    }
    cout << setw(40) << right << "========" << endl;
    cout << setw(15) << "GRAND TOTAL" << setw(25) << right << grand_total << endl << endl;
    cout << setw(15) << "CHECK TOTAL" << setw(25) << right << check_total << endl << endl;
    
    // - Check that ALL combinations of pair types produce consistent totals.
    //
    for(unsigned int i = SAGE::pair_generator::PARENTAL_MASK; i < SAGE::pair_generator::NULL_TYPE_MASK; ++i)  
    {
      total = 0;
      check_total = 0;
      pg.set_types(i);
      SAGE::pair_generator::pair_type current = pg.first_type();
      
#ifdef PG_VERBOSE
      cout << "Pair types:  " << SAGE::pair_generator::pair_type_to_string(current);
#endif

      check_total += totals[current - 1];
      while(pg.next_type(current) != pg.first_type())
      {
        current = pg.next_type(current);
        check_total += totals[current - 1];
        
#ifdef PG_VERBOSE
        cout << ", " << SAGE::pair_generator::pair_type_to_string(current);
#endif

      }
      
#ifdef PG_VERBOSE
      cout << endl;
#endif
      
      for(SAGE::pair_generator::iterator iter = pg.begin(); iter != pg.end(); ++iter)
      {
        ++total;
      }
      
#ifdef PG_VERBOSE
      cout << "Total number of pairs:  " << total << endl << endl << endl;
#endif
      
      assert(check_total == total);    // THE TEST !!!!
      
    }
    cout << endl << "Totals are consistent for all combinations." << endl << endl;
  
#define TEST_FILTERING
#ifdef TEST_FILTERING
  
    // - TEST FILTERING PAIR GENERATOR.
    //   For each trait, for each pair type, for each affection status, count.
    //   For continuous traits define affected as above the mean and unaffected as
    //   equal to or below the mean.
    //
    SAGE::RPED::RefTraitInfo              rti;
    SAGE::RPED::RefTraitInfo::trait_t     rt_type;
    SAGE::RPED::RefMPedInfo               rmpi = p.info();
    size_t                    no_of_traits = rmpi.trait_count();
    
    size_t                    ftotals[9][4];
    size_t                    across[9];
    size_t                    down[4];
    size_t                    all = 0;
    
    double                    trait_mean = 0.0;
        
    SAGE::filtering_pair_generator  fpg;
    SAGE::pair_filter_trait         trait;
    SAGE::pair_filter               pf;
    
    for(size_t i = 0; i < no_of_traits; ++i)  // ##### TRAITS ######
    {
      rti = rmpi.trait_info(i);
      rt_type = rti.type();
      cout << setw(25) << " ";
      cout << "Trait: " << rti.name() << "   Type: " << SAGE::trait_t_to_string(rt_type);
      
      if(rt_type == SAGE::RPED::RefTraitInfo::continuous_trait || rt_type == SAGE::RPED::RefTraitInfo::discrete_trait)
      {
        // Calculate trait mean.
        SAGE::SampleInfo sample;
        SAGE::RPED::RefPedigree::member_const_iterator iter;
        SAGE::RPED::RefPedInfo inf = refped->info();
        
        for(iter = refped->member_begin(); iter != refped->member_end(); ++iter)
        {
          sample += inf.trait(iter->index(), i);
        }
        trait_mean = sample.mean();
        cout << "   mean: " << trait_mean;
      }
      cout << endl;
      
      if(    rt_type == SAGE::RPED::RefTraitInfo::continuous_trait || rt_type == SAGE::RPED::RefTraitInfo::discrete_trait
          || rt_type == SAGE::RPED::RefTraitInfo::binary_trait)
      {
        trait.set_trait(i);
        if(rt_type != SAGE::RPED::RefTraitInfo::binary_trait)
        {
          trait.set_threshold(trait_mean);
        }
      
        // - Get counts for each pair type and affection status.
        //
        for(unsigned int j = SAGE::pair_generator::PARENTAL_MASK; j < SAGE::pair_generator::NULL_TYPE_MASK; j *= 2)  // ##### PAIR TYPES #### 
        {
          pg.set_types(j);
          fpg.set_pair_generator(pg);
          
          for(unsigned int k = SAGE::pair_filter_trait::CONCORD_UNAFF_MASK; k <= SAGE::pair_filter_trait::UNINFORM_MASK; k *= 2) // ### AFFECTION STATUSES ####
          {
            trait.set_statuses(k);
            pf.clear_traits();
            pf.set_trait(trait);
            fpg.set_filter(pf);
            total = 0;
            int row;
            int col;
          
            for(SAGE::filtering_pair_generator::iterator iter = fpg.begin(); iter != fpg.end(); ++iter)  // ##### PAIRS ######
            {
              ++total;
            }
            row = SAGE::pair_generator::mask_to_pair_type((SAGE::pair_generator::mask)j) - 1;
            col = SAGE::pair_filter_trait::mask_to_status((SAGE::pair_filter_trait::mask)k) - 1;
            ftotals[row][col] = total;
          }
        }
        
        // - Calculate totals.
        //
        for(int row = 0; row < 9; ++row)
        {
          across[row] = 0;
        
          for(int col = 0; col < 4; ++col)
          {
            across[row] += ftotals[row][col];
          }
        }
        
        for(int col = 0; col < 4; ++col)
        {
          down[col] = 0;
        
          for(int row = 0; row < 9; ++row)
          {
            down[col] += ftotals[row][col];
          }
        }
        
        all = 0;
        for(int row = 0; row < 9; ++row)
        {
          all += across[row];
        }
           
        // - Print the results.
        //
        cout << endl;                                     // header
        cout << setw(15) << left << "    ";
        cout << setw(15) << right << "CONCORD_UNAFF";    
        cout << setw(15) << right << "DISCORD";
        cout << setw(15) << right << "CONCORD_AFF";
        cout << setw(15) << right << "UNINFORM";
        cout << setw(15) << right << "Totals";
        cout << endl;
        cout << setw(15) << left << "    ";
        cout << setw(15) << right << "-------------";    
        cout << setw(15) << right << "-------";
        cout << setw(15) << right << "-----------";
        cout << setw(15) << right << "--------";
        cout << setw(15) << right << "------";
        cout << endl;
        
                                                          // pair types and data
        for(int row = SAGE::pair_generator::PARENTAL - 1; row < SAGE::pair_generator::NULL_TYPE - 1; ++row)
        {
          cout << setw(15) << left << SAGE::pair_generator::pair_type_to_string((type)(row + 1));
          
          for(int col = SAGE::pair_filter_trait::CONCORD_UNAFF - 1; col <= SAGE::pair_filter_trait::UNINFORM - 1; ++col)
          {
            cout << setw(15) << right << ftotals[row][col];
          }        
          cout << setw(15) << right << across[row];
          cout << endl;
        }
      }
      
      cout << setw(15) << left << "    ";               // footer
      cout << setw(15) << right << "=============";    
      cout << setw(15) << right << "=======";
      cout << setw(15) << right << "===========";
      cout << setw(15) << right << "========";
      cout << setw(15) << right << "======";
      cout << endl;
      
      cout << setw(15) << left << "Totals";
      for(int i = 0; i < 4; ++i)
      {
        cout << setw(15) << right << down[i];    
      }
      cout << setw(15) << right << all;
      cout << endl << endl << endl;
    }
  }
  
#endif
  
  // - Print totals over all pedigrees.
  //  
  cout << "        TOTALS OVER ALL PEDIGREES" << endl << endl;
  cout << setw(15) << left << "Pair Type" << setw(25) << right << "Count" << endl;
  cout << setw(15) << left << "---------" << setw(25) << right << "--------" << endl;
  
  for(int i = SAGE::pair_generator::PARENTAL - 1; i < SAGE::pair_generator::NULL_TYPE - 1; ++i)
  {
    cout << setw(15) << left << SAGE::pair_generator::pair_type_to_string((type)(i + 1)) << setw(25) << right << multi_totals[i] << endl;
  }

  exit(EXIT_SUCCESS);
}


