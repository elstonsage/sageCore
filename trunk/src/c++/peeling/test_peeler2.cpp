#include "LSF/LSFfile.h"
#include "mlocus/mfile.h"
#include "LSF/LSFinit.h"
#include <string>
#include <cassert>
#include <functional>
#include <vector>
#include <iomanip>
#include "LSF/LSFsymbol.h"
#include "rped/rped.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "error/bufferederrorstream.h"
#include "peeling/test_peeler.h"
#include "peeling/peeler2.h"
#include "mped/mp.h"
#include "boost/smart_ptr.hpp"

class no_info1 {};

namespace SAGE {
namespace peeling {

template <>
class individual_cache<no_info1, int>
{
public:

  typedef no_info1   data_type;
  typedef int result_type;

  typedef RPED::RefMultiPedigree::member_type member_type;

  individual_cache() { t = 0; }

  // Functions to be implemented by Algorithm:

    // Information retrieval functions

  bool anterior_cached              (const data_type&) const;
  bool posterior_cached             (const data_type&) const { return false; }
  bool posterior_with_mate_cached   (const member_type& mate_index, const data_type&) const { return false; }
  bool posterior_except_mate_cached (const member_type& mate_index, const data_type&) const { return false; }

  result_type& anterior              (const data_type&)                    { return t; }
  result_type& posterior             (const data_type&)                    { return t; }
  result_type& posterior_with_mate   (const member_type& mate_index, const data_type&) { return t; }
  result_type& posterior_except_mate (const member_type& mate_index, const data_type&) { return t; }

  const result_type& anterior              (const data_type&)   const { return t; }
  const result_type& posterior             (const data_type&)   const { return t; }
  const result_type& posterior_with_mate   (const member_type& mate_index, const data_type&) const { return t; }
  const result_type& posterior_except_mate (const member_type& mate_index, const data_type&) const { return t; }

protected:

  // Actual data to store results are filled in by individual specializations

  result_type t;

};

bool individual_cache<no_info1, int>::anterior_cached(const data_type&) const
{ return false; }

}
} // End namespaces

class individual_counter : public SAGE::peeling::peeler<no_info1, int>
{
public:

  individual_counter(const subped_type& subped)
    : SAGE::peeling::peeler<no_info1, int>(subped) { }

protected:

  // Functions needing to be implemented by the algorithm developer

  const result_type& internal_anterior    (const member_type& ind, const data_type& g, result_type& r) 
  {
    // Initially count self
    
    int r2 = 1;

    // Ancestors through parents

    const member_type* mother = ind.parent1();
    const member_type* father = ind.parent2();

    // - No check for missing sex information is OK here because individual count 
    //   is not dependent on sex information. -djb
    //
    if(father->is_female() || mother->is_male())
      std::swap(mother, father);

    r2 += anterior(*mother, no_info1()) + anterior(*father, no_info1());

    // Add parents other mates

    r2 += posterior_except_mate(*mother, *father, no_info1());
    r2 += posterior_except_mate(*father, *mother, no_info1());

    // siblings and their posteriors

    SAGE::RPED::RefMultiPedigree::sibling_const_iterator sib     = ind.sibling_begin();
    SAGE::RPED::RefMultiPedigree::sibling_const_iterator sib_end = ind.sibling_end();

    for( ; sib != sib_end; ++sib)
      r2 += 1 + posterior(*sib, no_info1());

    return r=r2;
  }
  const result_type& internal_posterior   (const member_type& ind, const data_type& g, result_type& r)
  { 
    int r2 = 0;

    // Add the mate, the children of that mate, and the anterior of that mate for
    // all mates.

    SAGE::RPED::RefMultiPedigree::mate_const_iterator mate = ind.mate_begin();

    for( ; mate != ind.mate_end(); ++mate)
      r2 += posterior_with_mate(ind, mate->mate(), no_info1())
         +  anterior(mate->mate(), no_info1());
    
    return r = r2;
  }

  const result_type& internal_posterior_with_mate
    (const member_type& ind, const member_type& mate, const data_type& g, result_type& r)
  {
    int r2 = 0;

    SAGE::RPED::RefMultiPedigree::offspring_const_iterator child = ind.offspring_begin(mate);

    for( ; child != ind.offspring_end(); ++child)
    {
      r2 += 1 + posterior(*child, no_info1());
    }    

    return r=r2;
  
  }

  const result_type& internal_posterior_except_mate
    (const member_type& ind, const member_type& mate, const data_type& g, result_type& r)
  {
    return r = 0;

  }

  const result_type& internal_anterior_terminal(const member_type& ind, const data_type& g, result_type& r)
  {
    return r = 1;
  }

  const result_type& internal_posterior_terminal(const member_type& ind, const data_type& g, result_type& r)
  {
    return r = 0; 
  }
  

};

#define DEBUG(x)

void test_peeling(const SAGE::RPED::RefMultiPedigree::subpedigree_type& sp);

namespace SAGE
{

string fp(double d, size_t w, size_t p=4)
{
  string n = doub2str(d, w, p, ios::showpoint | ios::fixed);
  return n.substr(0,w);
}

int TESTPEELER::main()
{
  init();

  init_output();

  data.input(argc, argv);

  return 0;
}

void TESTPEELER::init()
{
  LSFInit();
}

void TESTPEELER::init_output()
{
// Initialize Error Handling.

  // Eventually add support for changing these names. - GCW
  ostream* o21 = new ofstream("genome.inf");

  data.messages.clear();
  data.messages.open("test_peeler.inf");

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

  error_file.prefix("%%TESTPEELER-%p: ");
  info_file.prefix("     ");

  sage_cerr.prefix("%%TESTPEELER-%p: ");
  sage_cout.prefix("%%TESTPEELER-%p: ");

  data.information.insert(info_file);

  data.errors.insert(sage_cerr);
  data.errors.restrict(r_ge, error);

  data.errors.insert(error_file);
  data.errors.restrict(r_ge, information);
}


void TESTPEELER::print_help( ostream & )
{
}

}

using namespace SAGE;

int main(int argc, char **argv)
{
  free(malloc(1));

  TESTPEELER *genibd = new TESTPEELER(argc,argv);
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

test_peeler_data::~test_peeler_data()
{ }

bool test_peeler_data::input(int _argc, char* _argv[])
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

void test_peeler_data::read_par()
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

void test_peeler_data::read_ldf(char sep)
{
  MLOCUS::inheritance_model_map& imap = pedigree_data.info().markers();

  MLOCUS::InheritanceModelFile ifile;

  ifile.set_marker_verbose_output(0);

  ifile.input(imap, argv[3]);
}

void test_peeler_data::read_fdf()
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

  boost::scoped_ptr<RPED::RefLSFPedigreeFile> ped_reader;

  if( pedigree->attrs() && pedigree->attrs()->has_attr("column") )
  {
    ped_reader.reset(new RPED::RefLSFFortranPedigreeFile(errors));
  }
  else
  {
    ped_reader.reset(new RPED::RefLSFDelimitedPedigreeFile(errors));
  }

  ped_reader->process_parameters(pedigree_data.info(), pedigree);

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

  ped_reader->output(pedigree_data, "test_peel.out");

//  size_t t = pedigree_data.info().trait_find("T1");

//  if(t != (size_t) -1)
  {
    // Create pedigree subsets:

    RPED::RefMultiPedigree::pedigree_const_iterator ped = pedigree_data.pedigree_begin();

    for( ; ped != pedigree_data.pedigree_end(); ++ped)
    {
      RPED::RefMultiPedigree::subpedigree_const_iterator subped = ped->subpedigree_begin();
      for( ; subped != ped->subpedigree_end(); ++subped)
      {
        std::cout << subped->name() << std::endl;

        test_peeling(*subped);
      }
    }
  }
}

void test_peeling(const SAGE::RPED::RefMultiPedigree::subpedigree_type& sp)
{
  individual_counter c(sp);

  for(int i = 0; i < (int)sp.member_count(); ++i)
  {
    cout << sp.member_index(i).name() << " : " << flush;
    cout << c.anterior (sp.member_index(i), no_info1()) << ' ' << flush; 
    cout << c.posterior(sp.member_index(i), no_info1()) << endl;
  }
}

