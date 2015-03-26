//==========================================================================
//  File:       test_mcmc_analysis.cpp
//
//  Author:     Yeunjoo Song
//
//  History:    Initial implementation.                              May. 04
//
//  Notes:      This class is the main driver performing analysis.
//
//  Copyright (c) 2004 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "mcmc/test_mcmc_analysis.h"
#include "mcmc/founder_allele_graph.h"

namespace SAGE
{

namespace MCMC
{

test_mcmc_analysis::test_mcmc_analysis(const RPED::RefMultiPedigree&  mp,
                                       const test_mcmc_parameters&    param,
                                       RPED::genome_description*      genome,
                                       ostream&                       inf,
                                       cerrorstream&                  err)
                  : info(inf), errors(err)
{
  if( !(&mp) || !genome )
    return;

  my_multipedigree = &mp;
  my_parameters    = &param;
  my_genome        = genome;
}

test_mcmc_analysis::~test_mcmc_analysis()
{}

bool
test_mcmc_analysis::run_analysis()
{
  if( !my_multipedigree || !my_genome )
    return false;

  cout << "Processing Analysis...." << endl;

  my_parameters->dump_parameters(info);

  if( !do_analysis() )
    return false;

  return true;
}

//
//-----------------------------------------------------------
//

bool
test_mcmc_analysis::do_analysis()
{
  my_ped_region.set_errorstream(errors);

//  string title = my_parameters->get_title();

  test_mcmc_region_iterator ri = my_parameters->region_begin();

  for( ; ri != my_parameters->region_end(); ++ri )
  {
    string output_name = ri->output;

    // Get rid of funky characters.  Replace with _
    //
    for( size_t j = 0; j < output_name.size(); ++j )
      if( strchr(" \n\r\t;:'\"", output_name[j]) )
        output_name[j]='_';

    cout << "Processing Region: " << ri->name << endl
         << "=======================================" << endl;

    info << endl
         << "Processing Region: " << ri->name << endl
         << "=======================================" << endl;

    // If we've removed markers

    const RPED::genome_description::region_type& r  = my_genome->region(ri->name);

    // Create a pedigree iteration
    for( size_t mp = 0; mp < my_multipedigree->pedigree_count(); ++mp )
    {
      const RPED::RefPedigree& fp = my_multipedigree->pedigree_index(mp);

      process_pedigree(output_name, fp, r);
    }

    cout << endl;
  }

  return true;
}

void
test_mcmc_analysis::process_pedigree(const string& output,
                                     const RPED::RefPedigree& rped,
                                     const RPED::genome_description::region_type& r)
{
  if( !rped.subpedigree_count() ) return;

  cout << endl << "Pedigree " << rped.name() << endl;
  info << endl << "Processing Pedigree " << rped.name() << ".............."
       << endl << endl;

  if( rped.subpedigree_count() > 1 )
    errors << priority(information)
           << "Pedigree '" << rped.name()
           << "' has " << rped.subpedigree_count()
           << " subpedigrees."
           << "  Each will be processed individually."
           << endl;

  for( size_t sp = 0; sp < rped.subpedigree_count(); ++sp )
  {
    const RPED::RefSubpedigree& subped = rped.subpedigree_index(sp);

    string name = "Pedigree '" + rped.name() + "'";
    string founders = "";

    FPED::Multipedigree fmp(*rped.multipedigree());

    typedef FPED::has_informative_loci<SAGE::RPED::Member>             has_inf_loci;
  

    FPED::FilterResults fresults =
        FPED::MPFilterer::add_subpedigree_filtered_by_members
          (fmp, subped, FPED::is_inf_within_sped(has_inf_loci(*subped.multipedigree())));

    fmp.construct();

    if( !fmp.pedigree_count() || !fmp.pedigree_index(0).subpedigree_count() )
      continue;

    const FPED::Subpedigree& fsubped = fmp.pedigree_index(0).subpedigree_index(0);

    McmcMeiosisMap mm(fsubped, r.is_x_linked());

    if( rped.subpedigree_count() > 1 )
    {
      name = "Subpedigree " + long2str(sp + 1) + " of pedigree '";
      name += (rped.name() + "'");

      founders = "founders (" + mm.get_subpedigree().member_index(0).name();

      for( size_t f = 1; f < fsubped.member_count(); ++f )
      {
        if( mm.get_subpedigree().member_index(f).is_founder() )
          founders += (", " + mm.get_subpedigree().member_index(f).name());
      }

      founders += ")";
    }

    if( subped.member_count() != fsubped.member_count() )
    {
      string removed_member = "";

      errors << priority(information)
             << "The following individual(s) were removed as uninformative "
             << "before processing " << name
             << "': (";

      FPED::FilterResults::MemberPtrIterator rem_mem = fresults.get_excluded_member_begin();

      errors << (*rem_mem)->name();

      for(++rem_mem; rem_mem != fresults.get_excluded_member_end(); ++rem_mem )
        errors << ", " << (*rem_mem)->name();

      errors << ")." << endl;
    }

    if( mm.get_nonfounder_count() < 1 )
    {
      errors << priority(information) << name
             << " has no valid pairs after removing uninformative "
             << "individuals.  It will be skipped." << endl;

      continue;
    }

    if( rped.subpedigree_count() > 1 )
    {
      cout << endl;

      cout << "    " << name << " : " << founders << endl << endl;

      errors << priority(information)
             << name <<" includes " << founders << "." << endl;
    }

    process_subpedigree(output, name, r, mm);
  }
}

void
test_mcmc_analysis::process_subpedigree(const string&           output,
                                        const string&           name,
                                        const RPED::genome_description::region_type&      r,
                                        const McmcMeiosisMap&   mmap)
{
  test_mcmc_engine_components(r, mmap);
}


void
test_mcmc_analysis::test_mcmc_engine_components(const RPED::genome_description::region_type&      r,
                                                const McmcMeiosisMap& mmap)
{
  my_ped_region.build(mmap.get_subpedigree(), r, true);

  for( size_t marker = 0; marker < r.locus_count(); ++marker )
    cout << "marker informative? " << boolalpha << my_ped_region[marker].genotype_informative() << endl;

  mmap.dump_map(cout);

  mcmc_data_accessor mda(mmap, r.locus_count());
  mda.dump_accessor(cout);

  mda.set_valid_locus(0);
  mda.set_valid_locus(1);

  marker_likelihood_calculator mlc(my_ped_region, mmap, mda);

  for( size_t pattern = 0; pattern < 16; ++pattern )
  {
    bit_field  b(6, (uint) pattern);

    cout << "inheritance vector  " << b.to_string() << endl;

    for(size_t i = 0; i < r.locus_count(); ++i)
    {
      // Don't really care about the log likelihood here, but this is how
      // we set the patterns
      mlc.log_likelihood(i, b);
    }

    mlc.dump_graphs(cout);
  }
}

} // end of namespace MCMC

} // end of namespace SAGE

