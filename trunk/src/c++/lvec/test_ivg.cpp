//==================================================================
//  File:     test_ivg.cpp                         
//
//  Author:   Yeunjoo Song                      
//
//  History:  Initial implementation.  Mar. 2002
//
//  Copyright (c) 2002 R. C. Elston
//==================================================================

#include <string>
#include <cassert>
#include <functional>
#include <vector>
#include <iomanip>
#include "LSF/LSF.h"
#include "LSF/LSFinit.h"
#include "LSF/LSFfactory.h"
#include "LSF/ErrorHandler.h"
#include "LSF/LSFsymbol.h"
#include "error/errormanip.h"
#include "error/errorstream.h"
#include "error/bufferederrorstream.h"
#include "gelim/geno_eliminate.h"
#include "gelim/pedigree_region.h"
#include "fped/fped.h"
#include "lvec/meiosis_map.h"
#include "lvec/test_ivg.h"
#include "lvec/test_ivg_input.h"

#include "test_iv.cpp"

using namespace SAGE;

TEST_IVG::TEST_IVG(int argc, char **argv)
          : APP::SAGEapp(APP::APP_MLOD, false, argc, argv), my_consistent_out(false), my_summary_file(NULL)
{
  if (arg_count < 2)
  {
    print_help(cerr);
    exit(EXIT_FAILURE);
  }

  LSFInit();
}

TEST_IVG::~TEST_IVG()
{}

void TEST_IVG::print_help(ostream &o)
{
  o << "usage: " << argv[0] << "  <parameters> <pedigree> <locus> "
    << endl << endl
    << "Command line parameters:" << endl
    << "  parameters   - Parameter File"         << endl
    << "  pedigree     - Pedigree Data File"     << endl
    << "  locus        - Locus Description File" << endl
    << endl << endl;
}

int TEST_IVG::main()
{
  // Create input data.
  //
  test_ivg_data mdata(name, debug());

  print_title(mdata.info());

  // Get the error stream and make it available locally
  //
  cerrorstream& errors = mdata.errors();

  mdata.input(argc, argv);

  // Start TESTIVG analysis.
  //
  cout << endl << "TEST_IVG analysis......." << endl << endl << flush;

  if( mdata.analysis().size() )
    parse_analyses(mdata.analysis()[0], mdata.pedigrees(), errors);

  boost::scoped_ptr<ofstream> s_file;
  string   sum_filename;
  ofstream sum_file;

  if( my_outfile_name.size() )
  {
    sum_filename = my_outfile_name + ".out";
    sum_file.open( sum_filename.c_str() );
    if(sum_file)
      print_title(sum_file);
    my_summary_file = &sum_file;
  }
  else
  {
    sum_filename = "test_ivg.out";
    if( s_file.get() == NULL )
    {
      s_file.reset( new ofstream(sum_filename.c_str()) );
      if(*s_file)
        print_title(*s_file);
    }      
    my_summary_file = s_file.get();
  }

  run_analyses(mdata);

  cout << endl << "Analysis complete!" << endl << endl;

  return EXIT_SUCCESS;
}

void TEST_IVG::parse_analyses(LSFBase* params, const RPED::RefMultiPedigree& mp, cerrorstream& errors)
{
  AttrVal v = attr_value(params, "out");
  if( !v.has_value() || !v.String().size() )
    v = attr_value(params, "output");
  if( v.has_value() && v.String().size() )
    my_outfile_name = v.String();

  if( !params->List() )
    return;

  LSFList::const_iterator i;
  AttrVal a;
  for( i = params->List()->begin(); i != params->List()->end(); ++i)
  {
    if( !*i || !((*i)->name().size()) ) continue;

    string name = toUpper( (*i)->name() );

    if( name == "SAMPLE_ID" )
    {
      a = attr_value(*i, "SAMPLE_ID", 0);
      if( a.has_value() )
      {
        size_t string_index = mp.info().string_find(a.String());

        if( string_index == (size_t) - 1 )
          errors << priority(information) << "Invalid sample_id name : "
                 << a.String() << "\n             Skipped... " << endl;
        else
          my_sample_ids.push_back(string_index);
      }
    }
    else if( name == "CONSISTENT_OUT" )
    {
      a = attr_value(*i, "CONSISTENT_OUT", 0);
      if( a.has_value() )
      {
        if(    toUpper(a.String()) == "TRUE"
            || toUpper(a.String()) == "YES" )
          my_consistent_out = true;
      }
    }
    else if( name == "PEDIGREE_SKIP" )
    {
      a = attr_value(*i, "PEDIGREE_SKIP", 0);
      if( a.has_value() )
      {
        if( mp.pedigree_find(a.String() ) )
        {
          size_t pedigree_index = mp.pedigree_find(a.String())->index();
          my_pedigree_skips.insert(pedigree_index);
        }
        else
          errors << priority(information) << "Invalid pedigree name : "
                 << a.String() << "\n             Skipped... " << endl;
      }
    }
    else
      errors << priority(information) << "Unknown option : "
             << name << "\n             Skipped... " << endl;
  }

  return;
}

void TEST_IVG::run_analyses(const test_ivg_data& mdata)
{
  genotype_eliminator gelim;

  const RPED::RefMultiPedigree& ref_mp = mdata.pedigrees();

  FPED::Multipedigree fped(ref_mp);

  FPED::MPFilterer::add_multipedigree(fped, ref_mp);

  fped.construct();

  FPED::PedigreeConstIterator ped     = fped.pedigree_begin();

  vector< string > ped_name_sorted;

  for( ; ped != fped.pedigree_end(); ++ped )
  {
    const FPED::Pedigree& pedigree = *ped;

    ped_name_sorted.push_back(pedigree.name());
  }

  sort(ped_name_sorted.begin(), ped_name_sorted.end());

  for( size_t ped_name = 0; ped_name < ped_name_sorted.size(); ++ped_name )
  {
    const FPED::Pedigree& pedigree = *fped.pedigree_find(ped_name_sorted[ped_name]);

    char old_fill = cout.fill('.');

    string s = "'" + pedigree.name() + "'";

    if( my_pedigree_skips.find(pedigree.index()) != my_pedigree_skips.end() )
    {
      cout << "Skipping pedigree " << s << "........." << endl << flush;
      continue;
    }

    cout << "Processing pedigree " << s << flush << endl;

#if 0
    cout << "member in mped : " << endl;
    RPED::filtered_multipedigree::pedigree_type::member_const_iterator m_i = pedigree.member_begin();
    for( ; m_i != pedigree.member_end(); ++m_i )
      cout << "name = " << m_i->name() << ", index = " << m_i->index() << endl;
    cout << endl;
#endif

    FPED::SubpedigreeConstIterator subped     = pedigree.subpedigree_begin();
    FPED::SubpedigreeConstIterator subped_end = pedigree.subpedigree_end();
    for( ; subped != subped_end; ++subped)
    {
        meiosis_map sect;

        sect.set_subpedigree(&(*subped));

        sect.build(); // Build the pedigree section.

        meiosis_map mm1(sect);

        test_meiosis_map(cout, mm1);

        inheritance_vector iv(mm1);

        test_inheritance_vector(cout, iv);

        test_descent_graph(cout, iv);

        SAGE::pedigree_imodel_generator gen;

        MLOCUS::inheritance_model model;

        for(size_t i = 0; i < fped.info().marker_count(); ++i)
        {
          gen.set_prior_remap(true);
          gen.set_genotype_elimination(true);
          gen.set_post_remap(true);

          model = gen(*subped, i);
          
          test_iv_generator(cout, mm1, model);
          
          if(model.codominant(false))
            test_codom_iv_generator(cout, mm1, model);
          else
            cout << "Skipping codominant test.  Model not codominant" << endl;
          
          test_fixed_bit_calculator(cout, *subped, model);
        }
    }

    cout << left << setw(22-s.size()) << "" << "done."<< endl;
    cout.fill(old_fill);
  }
}

int main(int argc, char **argv)
{
  free(malloc(1));

  TEST_IVG *test_ivg = new TEST_IVG(argc,argv);
  assert(test_ivg != NULL);

  int r = test_ivg->main();

  delete test_ivg;

  return r;
}
