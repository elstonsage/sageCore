#include "LSF/LSFfile.h"
#include "boost/smart_ptr.hpp"
#include "mlocus/mfile.h"
#include "LSF/LSFinit.h"
#include <string>
#include <cassert>
#include <functional>
#include <vector>
#include <iomanip>
#include "LSF/LSFsymbol.h"
#include "fped/fped.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "error/bufferederrorstream.h"
#include "gelim/test_gelim.h"
#include "mped/mp.h"
#include "gelim/ped_imodel_gen.h"
#include "rped/genome_description.h"
#include "gelim/pedigree_region.h"

#define DEBUG(x)

void test_marker_generation(const SAGE::RPED::RefMultiPedigree::subpedigree_type& sect, size_t marker_count);
void test_pedigree_region(const SAGE::RPED::RefMultiPedigree::subpedigree_type& sect, size_t marker_count);

namespace SAGE
{
  
int GELIMTEST::main()
{
  init();

  init_output();

  data.input(argc, argv);

  return 0;
}

void GELIMTEST::init()
{
  LSFInit();
}

void GELIMTEST::init_output()
{
// Initialize Error Handling.

  // Eventually add support for changing these names. - GCW
  ostream* o21 = new ofstream("genome.inf");

  data.messages.clear();
  data.messages.open("genibd.inf");

  if(!o21->good() || !data.messages.good())
  {
    sage_cerr << priority(fatal) << "Unable to initalize output files.  "
              << "Unable to continue.  Exiting." << endl;;
    exit(1);
  }

  print_title(*o21);
  print_title(data.messages);

  cerrorstream error_file(data.messages);
  cerrorstream info_file(cout);

  error_file.prefix("%%GELIMTEST-%p: ");
  info_file.prefix("     ");

  sage_cerr.prefix("%%GELIMTEST-%p: ");
  sage_cout.prefix("%%GELIMTEST-%p: ");

  data.information.insert(info_file);

  data.errors.insert(sage_cerr);
  data.errors.restrict(r_ge, error);

  data.errors.insert(error_file);
  data.errors.restrict(r_ge, information);
}

void GELIMTEST::print_help( ostream & )
{
}

} // SAGE_NS_END

using namespace SAGE;

int main(int argc, char **argv)
{
  free(malloc(1));

  GELIMTEST *genibd = new GELIMTEST(argc,argv);
  assert(genibd != NULL);
  genibd->print_title(cout);

  if( argc != 4 && argc != 5)
  {
    cerr << "Oops.  Wrong args." << endl;

    return -1;
  }

  int r = genibd->main();

  delete genibd;

  return r;
}

gelim_data::~gelim_data()
{ }

bool gelim_data::input(int _argc, char* _argv[])
{
  argc = _argc;
  argv = _argv;

#ifdef DEBUG_OUTPUT
  cout << "Reading Parameter File" << endl;
#endif

  read_par();

  // Read in our Marker Map.

  char c = '/';

  AttrVal v = attr_value(params->find("ALLELE_DELIMITER"),0);

  if( v.has_value() && v.String().size()) c = v.String()[0];
  else
  {
    // Included for backwards portability.
    AttrVal v = attr_value(params->find("MARKER_DELIMITER"),0);

    if( v.has_value() && v.String().size()) c = v.String()[0];
  }

#ifdef DEBUG_OUTPUT
  cout << "Reading Locus Description File" << endl;
#endif

  read_ldf(c);

#ifdef DEBUG_OUTPUT
  cout << "Reading Family Data File" << endl;
#endif

  read_fdf();

#ifdef DEBUG_OUTPUT
  cout << "Reading Genome Description File" << endl;
#endif

  return true;
}

void gelim_data::read_par()
{
  ifstream in_state(argv[1]);

  if( !in_state.good() )
  {
    sage_cerr << priority(fatal) << "Unable to load parameters.  Exiting"
              << endl;

    exit(3);
  }

  params = new SymbolTable();

  LSF_input load_state(in_state, cout);

  assert ( load_state.good() );

  LSF_ptr<LSFBase> s = new LSFBase("Parameters");

  load_state.input_to(s, false);

  if(!s || !s->List() )
  {
    sage_cerr << priority(fatal) 
              << "Unable to load parameters.  Please check your data file."
              << endl;

    exit(1);
  }
  for(LSFList::iterator i=s->List()->begin(); i!= s->List()->end(); i++)
  {
    if(!*i) continue;

    params->add( toUpper((*i)->name()), (*i) );
  }
}

void gelim_data::read_ldf(char sep)
{
  MLOCUS::inheritance_model_map& imap = pedigree_data.info().markers();

  MLOCUS::InheritanceModelFile ifile;

  ifile.set_marker_verbose_output(0);

  ifile.input(imap, argv[3]);
  ifile.output(imap, cout);
}

void gelim_data::read_fdf()
{
  cout << "Reading Pedigree File....................." << endl;

  LSFBase* pedigree = params->find("PEDIGREE");

  if( !pedigree )
  {
    sage_cerr << priority(fatal)
              << "No pedigree specified....  Exiting."
              << endl;

    exit(EXIT_FAILURE);
  }

  boost::scoped_ptr<RPED::RefLSFPedigreeFile> ped_reader;

  string pedigree_file;

  if( pedigree->attrs() )
    pedigree_file = pedigree->attrs()->StringAttr("file");

  if( !pedigree_file.size() )
    pedigree_file = argv[2];

  if( pedigree->attrs() && pedigree->attrs()->has_attr("column") )
    ped_reader.reset(new RPED::RefLSFFortranPedigreeFile());
  else
    ped_reader.reset(new RPED::RefLSFDelimitedPedigreeFile());

  cout << "              from " << pedigree_file << flush;

  ped_reader->set_force_skip_markers(false);
  ped_reader->set_force_skip_traits(false);
  ped_reader->set_force_dynamic_markers(false);

  ped_reader->process_parameters(pedigree_data.info(), pedigree);

  if( !ped_reader->input(pedigree_data, pedigree_file) )
  {
    cout << endl;
    sage_cerr << priority(fatal)
              << "Error reading pedigree data from file '" << pedigree_file
              << "'."
              << endl;
    
    exit(EXIT_FAILURE);
  }

  char old_fill = cout.fill('.');
  cout << setw(23-pedigree_file.size()) << "" << "done." << endl;
  cout.fill(old_fill);     

  ped_reader->print_mped(pedigree_data, pedigree_file, cout, true, true);

  if( !pedigree_data.pedigree_count() )
  {
    sage_cerr << priority(critical)
              << "No pedigree to analyze....  Exiting."
              << endl;
    exit(EXIT_FAILURE);
  }

  size_t error_count = 0;
  RPED::RefMultiPedigree::pedigree_const_iterator p;
  for( p = pedigree_data.pedigree_begin(); p != pedigree_data.pedigree_end(); ++p )
    error_count += p->error_count();

  if( error_count )
  {
    sage_cerr << priority(error)
              << "Errors appear in pedigree data. "
              << "Results may be incomplete."
              << endl;
  }
    
  cout << "Sorting pedigrees........................." << flush;
  RPED::RefMultiPedigree::pedigree_iterator j;
  for( j  = pedigree_data.pedigree_begin(); j != pedigree_data.pedigree_end(); ++j )
    PedigreeSort(*j);
  cout << "done." << endl;


  // Create pedigree subsets:

  RPED::RefMultiPedigree::pedigree_const_iterator ped = pedigree_data.pedigree_begin();

  for( ; ped != pedigree_data.pedigree_end(); ++ped)
  {
    RPED::RefMultiPedigree::subpedigree_const_iterator sect = ped->subpedigree_begin();

    for( ; sect != ped->subpedigree_end(); ++sect)    
    {
      test_marker_generation(*sect, pedigree_data.info().markers().size());
      test_pedigree_region(*sect, pedigree_data.info().markers().size());
    }
  }

}

void test_marker_generation(const SAGE::RPED::RefMultiPedigree::subpedigree_type& sect, size_t marker_count)
{
  SAGE::pedigree_imodel_generator gen;

  MLOCUS::inheritance_model model;

  // Build RPED::filtered_multipedigree.

  FPED::Multipedigree fped(*sect.multipedigree());
  
  typedef FPED::has_informative_loci<SAGE::RPED::Member>             has_inf_loci;
  
  FPED::MPFilterer::add_subpedigree_filtered_by_members(fped, sect, is_inf_within_sped(has_inf_loci(*sect.multipedigree())));
  
  fped.construct();

  const FPED::Subpedigree& fsect = fped.pedigree_index(0).subpedigree_index(0);

  for(size_t i = 0; i < marker_count; ++i)
  {
    cout << "Analyzing Marker: " << i << endl;

    // Initial test, only creating the pedigree based imodel.
/*
    gen.set_prior_remap(false);
    gen.set_genotype_elimination(false);
    gen.set_post_remap(false);

    model = gen(fsect, i);

    SAGE::inheritance_model_test(cout, model, "Initial Version");
    SAGE::MLOCUS::inheritance_model_print(cout, model, "Initial Version");

*/
    // Remap test1

    gen.set_prior_remap(true);
    gen.set_genotype_elimination(false);
    gen.set_post_remap(false);

    model = gen(fsect, i);

    SAGE::MLOCUS::inheritance_model_test(cout, model, "Remap 1");
    SAGE::MLOCUS::inheritance_model_print(cout, model, "Remap 1");

    cout << model.genotype_informative() << gen.informative() << gen.inconsistent()
         << endl;

    // Genotype_elimination test

    gen.set_prior_remap(true);
    gen.set_genotype_elimination(true);
    gen.set_post_remap(false);

    model = gen(fsect, i);

    SAGE::MLOCUS::inheritance_model_test(cout, model, "Genotype Elimination");
    SAGE::MLOCUS::inheritance_model_print(cout, model, "Genotype Elimination");

    cout << model.genotype_informative() << gen.informative() << gen.inconsistent()
         << endl;

    if(gen.inconsistent())
    {
      SAGE::inconsistency_printer printer;

      printer.print_table(gen.get_errors());

      gen.clear_errors();
    }

  }
}

void test_pedigree_region(const SAGE::RPED::RefMultiPedigree::subpedigree_type& sect, size_t marker_count)
{
  cout << "Testing Pedigree Region on " << sect.name() << endl;

  // First, create a genome_description to use

  RPED::genome_description genome(sect.multipedigree()->info());

  genome.add_region("All Markers");

  for(size_t i = 0; i < marker_count; ++i)
  {
    genome.add_locus(i, i);
  }

  genome.build();
  genome.freeze();

  // Build filtered_multipedigree.

  FPED::Multipedigree fped(*sect.multipedigree());

  typedef FPED::has_informative_loci<SAGE::RPED::Member>             has_inf_loci;
  
  FPED::MPFilterer::add_subpedigree_filtered_by_members(fped, sect, is_inf_within_sped(has_inf_loci(*sect.multipedigree())));

  fped.construct();

  const FPED::Subpedigree& fsect = fped.pedigree_index(0).subpedigree_index(0);

  // Now we build the pedigree_region

  pedigree_region ped_region(fsect, genome.region(0));

  // Output some stuff

  assert(genome.region(0).name() == ped_region.get_region().name());
  assert(fsect.pedigree()->name() == ped_region.get_subpedigree().pedigree()->name());

  cout << ped_region.get_region().name() << endl;
  cout << ped_region.get_subpedigree().pedigree()->name() << endl;

  assert(ped_region.is_built());
  assert(ped_region.inheritance_model_count() == marker_count);

  for(size_t i = 0; i < marker_count; ++i)
  {
    MLOCUS::inheritance_model model = ped_region[i];

    SAGE::MLOCUS::inheritance_model_test  (cout, model, "Pedigree Region");
    SAGE::MLOCUS::inheritance_model_print (cout, model, "Pedigree Region");

    cout << model.genotype_informative() << endl;
  }
  
}
