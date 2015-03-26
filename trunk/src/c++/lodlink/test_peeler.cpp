//============================================================================
// File:      test_peeler.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   Created 10/31/2.                                                   
//                                                                          
// Notes:     Tests the lodlink likelihood engine (peeler class).
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
#include "lodlink/peeler.h"


using namespace std;
using namespace SAGE;
using namespace LODLINK;

// - Equation (2) in Fernando, Stricker, Elston.  1993.
//
log_double  
likelihood(peeler& inst,
           const FPED::FilteredMultipedigree::subpedigree_const_iterator& subped_iter,
           size_t trait, size_t marker,
           const FPED::FilteredMultipedigree::member_const_iterator& m_iter)
{
  log_double  like(0);
          
  phenoset  ph_set(*subped_iter, trait, marker, *m_iter);
  
  phenoset::phenoset_iterator  ph_iter = ph_set.begin();
  for(; ph_iter != ph_set.end(); ++ph_iter)
  {
    log_double like_term(1);
    like_term *= inst.anterior(*m_iter, *ph_iter);
    like_term *= (*ph_iter).penetrance();
    like_term *= inst.posterior(*m_iter, *ph_iter);  
    
    like += like_term;
  }
  
  return like;
}

void  
print_diagnostic_info(ostream& out, peeler& inst,
                      const FPED::FilteredMultipedigree::subpedigree_const_iterator& subped_iter,
                      size_t trait, size_t marker,
                      const FPED::FilteredMultipedigree::member_const_iterator& m_iter)
{
  out << "\n\n";

  phenoset  ph_set(*subped_iter, trait, marker, *m_iter);
  
  phenoset::phenoset_iterator  ph_iter = ph_set.begin();
  for(; ph_iter != ph_set.end(); ++ph_iter)
  {
    out << "individual " << m_iter->name() << "  "
        << "genotype "   << joint_genotype(*ph_iter) << "  "
        << "anterior "   << inst.anterior(*m_iter, *ph_iter).get_double() << "       "
        << "posterior "  << inst.posterior(*m_iter, *ph_iter).get_double()
        << endl;
  }
}


int main(int argc, char* argv[])
{

  // ============  The preliminaries.  =============
  //
  if (argc != 6)
  {
    cerr << "usage: " << argv[0] << " <parameters> <pedigree> <locus> <male_theta> <female_theta>\n\n"
         << "Command line parameters:\n"     
         << "  parameters   - Parameter File\n"  
         << "  pedigree     - Pedigree Data File\n"  
         << "  locus        - Marker Locus Description File\n"
         << "  male_theta   - Male Recombination Fraction\n"
         << "  female_theta - Female Recombination Fraction\n\n" 
         << endl;
    exit(EXIT_FAILURE);
  }

  LSFInit();

  sage_cerr << prefix("%%TESTPEELER-%P: ");

  // - Read parameter file.
  //
  LSFBase *params = loadLSFfile(argv[1], "TESTPEELER Parameter file", sage_cerr, false);

  if(!params)
  {
    sage_cerr << priority(fatal) << "Error reading parameter file.... aborting." << endl;
    exit(EXIT_FAILURE);
  }

  
  ofstream info_file;
  info_file.open("testpeeler.inf");

  if(!info_file)
  {
    sage_cerr << priority(fatal) 
              << "Cannot open output file: pedinfo.inf.  Exiting..." << endl;
    exit(EXIT_FAILURE);
  }

  ofstream out_file;
  out_file.open("testpeeler.out");

  if(!out_file)
  {
    sage_cerr << priority(fatal) 
              << "Cannot open output file: testpeeler.out.  Exiting..." << endl;
    exit(EXIT_FAILURE);
  }

  SAGE::cerrormultistream errors;
  cerrorstream error_file(info_file);
  error_file.prefix("%%TESTPEELER-%P:");
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
  
  // ===============  Test of peeler  ================
  //
  typedef SAGE::FPED::FilteredMultipedigree::pedigree_const_iterator     pedigree_const_iterator;
  typedef SAGE::FPED::FilteredMultipedigree::subpedigree_const_iterator  subpedigree_const_iterator;
  typedef SAGE::FPED::FilteredMultipedigree::member_const_iterator       member_const_iterator;
  
  const FPED::FilteredMultipedigreeInfo&  mp_info = fp.info();
  size_t  trait  = mp_info.marker_find("M1");
  size_t  marker = mp_info.marker_find("M2");
  assert(trait != (size_t)(-1));
  assert(marker != (size_t)(-1));
  
  mle_sub_model  mle;
  double  male_theta   = atof(argv[4]);
  double  female_theta = atof(argv[5]);
  
  if(male_theta == female_theta)
  {
    mle.set_average_theta(male_theta);
  }
  else
  {
    mle.set(true, false);
    mle.set_male_theta(male_theta);
    mle.set_female_theta(female_theta);
    cout << "male " << male_theta << endl;
    cout << "female " << female_theta << endl;
  }
  
  mle_sub_model  unlinked_mle;
  unlinked_mle.set_average_theta(.5);
  
  pedigree_const_iterator  ped_iter = fp.pedigree_begin();
  for(; ped_iter != fp.pedigree_end(); ++ped_iter)
  {
    subpedigree_const_iterator  subped_iter = ped_iter->subpedigree_begin();
    {
      for(; subped_iter != ped_iter->subpedigree_end(); ++subped_iter)
      {
        peeler  inst(*subped_iter, mle, trait, marker);
        peeler  unlinked_inst(*subped_iter, unlinked_mle, trait, marker);
        
        // - This is a consistancy check.  The likelihood should be the same
        //   regardless of which is the pivotal member in the calculation.
        //        
        member_const_iterator  m_iter = subped_iter->member_begin();
        for(; m_iter != subped_iter->member_end(); ++m_iter)
        {
          // #define DEBUG_PEELER
          #ifdef DEBUG_PEELER
          
          print_diagnostic_info(out_file, inst, subped_iter, trait, marker, m_iter);       
          
          #undef DEBUG_PEELER
          #endif
        
          log_double  like = likelihood(inst, subped_iter, trait, marker, m_iter);
          log_double  unlinked_like = likelihood(unlinked_inst, subped_iter, trait, marker, m_iter);
          
          out_file << "pedigree " << ped_iter->name() << ", "
                   << "subpedigree " << subped_iter->name() << ", "
                   << "member " << m_iter->name() << ", "
                   << "log likelihood " << log10(like.get_double()) << ", " 
                   << "lod score  " << log10(like.get_double()) - log10(unlinked_like.get_double()) 
                   << endl;
        }
      }
    }
  }
}



