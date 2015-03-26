#include "LSF/LSFinit.h"
#include <string>
#include <cassert>
#include <functional>
#include <vector>
#include <iomanip>
#include "LSF/LSFsymbol.h"
#include "LSF/LSFfile.h"
#include "rped.new/rped.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "error/bufferederrorstream.h"
#include "rped.new/rpedtest.h"
#include "rped.new/rpedtest_input.h"
#include "mped.new/mp.h"
#include "boost/smart_ptr.hpp"
#include "rped.new/pedigree_section.h"

#define DEBUG(x)

namespace SAGE {
namespace RPED {

int RPEDTEST::main()
{
  init();

  init_output();

  data.input(argc, argv);

  return 0;
}

void RPEDTEST::init()
{
  LSFInit();
}

void RPEDTEST::init_output()
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

  error_file.prefix("%%RPEDTEST-%p: ");
  info_file.prefix("     ");

  sage_cerr.prefix("%%RPEDTEST-%p: ");
  sage_cout.prefix("%%RPEDTEST-%p: ");

  data.information.insert(info_file);

  data.errors.insert(sage_cerr);
  data.errors.restrict(r_ge, error);

  data.errors.insert(error_file);
  data.errors.restrict(r_ge, information);
}


void RPEDTEST::print_title( ostream &o )
{
  // This function non-existent.
  //print_program_title(o, "RPEDTEST" );
}

void RPEDTEST::print_help( ostream & )
{
}

rptest_data::~rptest_data()
{}

SAGE_NS_END

using namespace SAGE;

int main(int argc, char **argv)
{
  free(malloc(1));

  RPEDTEST *genibd = new RPEDTEST(argc,argv);
  assert(genibd != NULL);
  genibd->print_title(cout);
  genibd->expire();

  if( argc != 4 && argc != 5)
  {
    cerr << "\nusage: " << argv[0] << " <parameters> <pedigree> <locus> [map]" 
         << endl << endl
         << "  params   - Parameter File"  << endl
         << "  pedigree - Pedigree Data File"       << endl
         << "  locus    - Locus Descriptor File"    << endl
         << "  map      - Genome Map File (optional for single point analysis)"
         << endl << endl;

    return -1;
  }

  int r = genibd->main();

  delete genibd;

  return r;
}

bool rptest_data::input(int _argc, char* _argv[])
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

  read_gdf();

  return true;
}

// - Read parameter file.
//
void rptest_data::read_par()
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

// - Read genome description file.
//
void rptest_data::read_gdf()
{
  LSFBase* gdf_params;
  if(argc == 5)
  {
    gdf_params = loadLSFfile(argv[4], "genome description file", sage_cerr, false);
  }
  else
  {
    return;
  }

  if(! gdf_params)
  {
    sage_cerr << priority(fatal) << "Error reading genome description file.... aborting." << endl;
    exit(EXIT_FAILURE);
  }
  
  LSFBase*  regions = NULL;
  for(LSFList::iterator i = gdf_params->List()->begin(); i != gdf_params->List()->end(); ++i)
  {
    if(toUpper((*i)->name()) == "GENOME")
    {
      regions = (*i); 
      break;
    }
  }

  if(!regions || !regions->List())
  {
    sage_cerr << priority(error)
              << "Error reading genome file. Please check your genome file." 
              << endl;
    exit(3);
  }
  
  LSFDump  tester(cout);
  tester.dump(regions);
  
  // - Test genome_description class.
  //
  LSFgenome_description  gd_tester(pedigree_data.info(), errors, regions);
  
  std::ofstream  outfile("rpedtest_gdf.out");
  
  genome_description::iterator iter;
  for(iter = gd_tester.begin(); iter != gd_tester.end(); ++iter)
  {
    outfile << "marker " << (*iter).name() << "  info ptr. " << (*iter).locus() 
            << " location " << (*iter).location() 
            << " number of alleles " << (*iter).locus()->allele_count() << endl;
    cout    << "marker " << (*iter).name() << "  info ptr. " << (*iter).locus() 
            << " location " << (*iter).location() 
            << " number of alleles " << (*iter).locus()->allele_count() << endl;
  }
}

// - Read locus descriptor file.
//
void rptest_data::read_ldf(char sep)
{
  MLOCUS::inheritance_model_map& imap = pedigree_data.info().markers();

  InheritanceModelFile ifile;

  ifile.set_marker_verbose_output(0);

  ifile.input(imap, argv[3], sep);
}

void rptest_data::read_fdf()
{
  cout << "Reading pedigree data..........." << endl;

  LSFBase* pedigree = params->find("PEDIGREE");

  if(!pedigree)
  {
    errors << priority(fatal)
           << "OOPS!  No pedigree.  Bad test."
           << endl;
  }

  string pedigree_file = argv[2];

  boost::scoped_ptr<RefPedigreeFile> ped_reader;

  if( pedigree->attrs() && pedigree->attrs()->has_attr("column") )
    ped_reader = boost::scoped_ptr<RefPedigreeFile>(
             new RefLSFFortranPedigreeFile(pedigree, pedigree_data.info(), errors) );
  else
    ped_reader = boost::scoped_ptr<RefPedigreeFile>(
             new RefLSFDelimitedPedigreeFile(pedigree, pedigree_data.info(), errors) );

  cout << "       from " << pedigree_file << flush;

  if( !ped_reader->input(pedigree_data, pedigree_file) )
  {
    cout << endl;
    errors << priority(fatal)
           << "Error reading pedigree data from file '" << pedigree_file
           << "'." << endl;
  }
  char old_fill = cout.fill('.');
  cout << setw(20-pedigree_file.size()) << "" << "done." << endl;
  cout.fill(old_fill);

  cout << "Printing Pedigree" << endl;

  ped_reader->output(pedigree_data, "rpedtest.out");

  size_t t = pedigree_data.info().trait_find("T1");

  if(t != (size_t) -1)
  {
    // Create pedigree subsets:

    RefMultiPedigree::pedigree_const_iterator ped = pedigree_data.pedigree_begin();

    for( ; ped != pedigree_data.pedigree_end(); ++ped)
    {
      pedigree_section sect(pedigree_data);

      sect.set_ped(&*ped);

      RefMultiPedigree::member_const_iterator ind = ped->member_begin();

      for( ; ind != ped->member_end(); ++ind)
      {
        if(!ped->info().trait_missing(ind->index(), t))
          sect.add_ind(&*ind, true);
        else
          sect.add_ind(&*ind, false);
      }

      sect.build(); // Build the pedigree section.

      cout << "Individual Count: " << sect.individual_count() << endl
           << "Unused Count:     " << sect.unused_count() << endl
           << "total Count:      " << sect.total_count() << endl
           << "Founders:         " << sect.founder_count() << endl
           << "Nonfounders:      " << sect.nonfounder_count() << endl
           << "bits:             " << sect.bit_count() << endl
           << "Family Count:     " << sect.family_count() << endl;
      
      for(size_t i = 0; i < sect.individual_count(); ++i)
      {
        cout << i << "\t" << sect.member(i)->name() << " "
             << sect.mother_location(i) << ' ' << sect.father_location(i) << " "
             << sect.child_count(i) << ' ' << sect.mate_count(i) << ' '
             << sect.sib_count(i) << endl;
        cout << "\t" << sect.founder(i) << ' ' << sect.nonfounder(i) << ' '
             << sect.parent(i) << ' ' << sect.male(i) << endl;
        cout << "\t" << sect.includes(i) << ' ' << sect.unused(i) <<  ' '
             << sect.informative(i) << endl;
      }
    }
  }
}

} // End namespace RPED
} // End namespace SAGE

