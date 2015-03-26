//==========================================================================
//  File:       analysis.cpp
//
//  Author:     Geoff Wedig & Yeunjoo Song
//
//  History:    Version 1.0  Initial implementation.
//                      2.0  Updated to new libraries            yjs Nov. 03
//
//  Notes:      This class is the main driver performing analysis.
//
//  Copyright (c) 2003 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "genibd/analysis.h"

namespace SAGE
{

namespace GENIBD
{

genibd_analysis::genibd_analysis(const RefMultiPedigree&   mp,
                                 const genibd_parameters&  param,
                                 genome_description*       genome,
                                 ostream&                  inf,
                                 cerrorstream&             err)
               : info(inf), errors(err)
{
  if( !(&mp) || !genome )
    return;

  my_multipedigree = &mp;
  my_parameters    = &param;
  my_genome        = genome;

  my_exact_ibd_analysis = NULL;
  my_pair_ibd_analysis  = NULL;
  my_sim_ibd_analysis   = NULL;
}

genibd_analysis::~genibd_analysis()
{
  if( my_exact_ibd_analysis )
    delete my_exact_ibd_analysis;

  if( my_pair_ibd_analysis )
    delete my_pair_ibd_analysis;

  if( my_sim_ibd_analysis )
    delete my_sim_ibd_analysis;
}

bool
genibd_analysis::run_analysis()
{
  if( !my_multipedigree || !my_genome )
    return false;

//  cout << "Processing Analysis...." << endl;

  // Write out the analysis so they remember what they're doing.
  //
  cout << "==================================================="
          "========================" << endl;

  my_parameters->dump_parameters(cout);

  cout << "==================================================="
          "========================" << endl;

  info << endl;
  info << "==================================================="
          "========================" << endl;

  my_parameters->dump_parameters(info);

  info << "==================================================="
          "========================" << endl;

  build();

  if( !do_analysis() )
    return false;

  return true;
}

//
//-----------------------------------------------------------
//

bool
genibd_analysis::build()
{
  init_likelihood_data();

  return true;
}

void
genibd_analysis::init_likelihood_data()
{
  analysis_data adata;

  bool header_printed = false;
  size_t total_length = 0;

  genibd_region_iterator ri = my_parameters->region_begin();

  for( ; ri != my_parameters->region_end(); ++ri )
  {
    region_type r  = my_genome->region(ri->name);

    // Generate the subpedigrees for this region

    my_bad_fam.resize(0);
    my_bad_ped.resize(0);

    my_bad_fam_name_size = 0;
    my_bad_ped_name_size = 0;

    string r_name = r.name();

    for( size_t mp = 0; mp < my_multipedigree->pedigree_count(); ++mp )  
    {
      const RefPedigree& ped = my_multipedigree->pedigree_index(mp);

      filtered_multipedigree fmp(*my_multipedigree);

      typedef FPED::has_informative_loci<SAGE::RPED::Member>             has_inf_loci;

      // Do local filtering for this region.
      //
      has_inf_loci hil(*my_multipedigree, false);

      for( size_t m = 0; m < r.locus_count(); ++m )
      {
        size_t m_index = r.locus(m).marker_index();
        hil.set_check_status_for_locus(m_index, true);
      }

      // Filter based upon being informative either locally (hil) or within
      // pedigrees.
      FPED::MPFilterer::add_pedigree_filtered_by_members(fmp, ped, is_inf_within_sped(hil));

      fmp.construct();

      if( !fmp.pedigree_count() || !fmp.pedigree_index(0).subpedigree_count() )
        continue;

      for( size_t fsp = 0; fsp < fmp.pedigree_index(0).subpedigree_count(); ++fsp )
      {
        const subped_type& fsubped = fmp.pedigree_index(0).subpedigree_index(fsp);

        meiosis_map mm(&fsubped, r.is_x_linked());

        if( split_pedigree(mm) )
        {
          for( size_t fam = 0; fam < fsubped.family_count(); ++fam )
          {
            const family_type& ffamily = fsubped.family_index(fam);

            filtered_multipedigree fmp2(*my_multipedigree);

            FPED::MPFilterer::add_family(fmp2, ffamily);

            fmp2.construct();

            if( !fmp2.pedigree_count() || !fmp2.pedigree_index(0).subpedigree_count() )
              continue;

            const subped_type& fsubfam = fmp2.pedigree_index(0).subpedigree_index(0);

            meiosis_map mm2(&fsubfam, r.is_x_linked());

            if( !check_pedigree(mm2, r, adata) )
            {
              my_bad_fam.push_back(mm2.get_pedigree()->name());
              my_bad_fam_name_size += (mm2.get_pedigree()->name()).size();
            }
          }
        }
        else if( !check_pedigree(mm, r, adata) )
        {
          my_bad_ped.push_back(mm.get_pedigree()->name());
          my_bad_ped_name_size += (mm.get_pedigree()->name()).size();
        }
      }
    }

    if( my_bad_ped.size() || my_bad_fam.size() )
    {
      size_t r_title = max(r_name.size(), size_t(6));

      if( !header_printed )
      {
        total_length += dump_header(r_title);
        header_printed = true;
      }

      info << " " << left << setw(r_title) << r.name() << " |";

      dump_bad_pedigrees(r_title);
    }
  }

  if( header_printed )
  {
    for( size_t a = 0; a < total_length; ++a )
      info << "-";
    info << endl << endl;
  }
}

size_t
genibd_analysis::dump_header(size_t r_title)
{
  size_t total_length = 0;

  info << endl
       << "There are ";

  if( my_bad_ped.size() )
    info << my_bad_ped.size() << " pedigrees (including subpedigrees)";

  if( my_bad_ped.size() && my_bad_fam.size() )
    info << " and ";
  
  if( my_bad_fam.size() )
    info << my_bad_fam.size() << " nuclear families";

  info << " which may not be processed." << endl
       << "Please check the information from the individual pedigree processing!" << endl;

  for( size_t r = 0; r < r_title + 3; ++r )
    info << "-";
  info << "-------------------";

  if( my_bad_ped_name_size + my_bad_fam_name_size < 7 )
    ;
  else if( my_bad_ped_name_size + my_bad_fam_name_size < 40 )
  {
    size_t total_size =   (my_bad_ped_name_size + my_bad_fam_name_size)
                        + (my_bad_ped.size() + my_bad_fam.size() - 1) * 2
                        - 7;
    for( size_t a = 0; a < total_size ; ++a )
      info << "-";
  }
  else
    info << "---------------------------------";

  info << endl;

  total_length = r_title + 22;
  if( my_bad_ped_name_size + my_bad_fam_name_size < 7 )
    ;
  else if( my_bad_ped_name_size + my_bad_fam_name_size < 40 )
    total_length += (  (my_bad_ped_name_size + my_bad_fam_name_size)
                     + (my_bad_ped.size() + my_bad_fam.size() - 1) * 2
                     - 7 );
  else
    total_length += 33;

  info << setw(r_title + 2)     << " Region " << "|";
  info << " Count | Pedigrees";
  info << endl;

  for( size_t r = 0; r < r_title + 3; ++r )
    info << "-";
  info << "-------------------";

  if( my_bad_ped_name_size + my_bad_fam_name_size < 7 )
    ;
  else if( my_bad_ped_name_size + my_bad_fam_name_size < 40 )
  {
    size_t total_size =   (my_bad_ped_name_size + my_bad_fam_name_size)
                        + (my_bad_ped.size() + my_bad_fam.size() - 1) * 2
                        - 7;
    for( size_t a = 0; a < total_size ; ++a )
      info << "-";
  }
  else
    info << "---------------------------------";

  info << endl;

  return total_length;
}

void
genibd_analysis::dump_bad_pedigrees(size_t r_title)
{
  if( my_bad_ped.size() )
  {
    info << " " << left << setw(6) << my_bad_ped.size()
                  << "| (";
    size_t a = 0;
    for( ; a < my_bad_ped.size()-1; ++a )
    {
      if(    (!(my_bad_ped_name_size < (40 - my_bad_fam_name_size)))
          && (a % 6 == 0)
          && a )
        info << endl
             << setw(r_title + 13) << " ";
      info << my_bad_ped[a] << ", ";
    }
    info << my_bad_ped[a];
    info << ")" << endl;
  }
  
  if( my_bad_ped.size() && my_bad_fam.size() )
    info << setw(r_title) << " " << " and ";
  
  if( my_bad_fam.size() )
  {
    if( my_bad_ped.size() )
      info << setw(r_title + 6) << " ";

    info << " " << setw(6) << my_bad_fam.size()
                  << "| (";
    size_t a = 0;
    for( ; a < my_bad_fam.size()-1; ++a )
    {
      if(    (!(my_bad_fam_name_size < (40 - my_bad_ped_name_size)))
          && (a % 6 == 0)
          && a )
        info << endl
             << setw(r_title + 13) << " ";
      info << my_bad_fam[a] << ", ";
    }
    info << my_bad_fam[a];
    info << ")" << endl;
  }
}

bool
genibd_analysis::do_analysis()
{
  my_ped_region.set_errorstream(errors);

  string title = my_parameters->title();

  genibd_region_iterator ri = my_parameters->region_begin();

  for( ; ri != my_parameters->region_end(); ++ri )
  {
    string output_name;

    if( ri->output.size() ) 
      output_name = ri->output;
    else
      output_name = my_parameters->output() + '.' + ri->name; // + ".ibd";

    // Get rid of funky characters.  Replace with _
    //
    for( size_t j = 0; j < output_name.size(); ++j )
      if( strchr(" \n\r\t;:'\"", output_name[j]) )
        output_name[j]='_';

    cout << endl
         << "  Processing Region: " << ri->name << endl
         << "  =======================================" << endl;

    info << endl
         << "Processing Region: " << ri->name << endl
         << "=======================================" << endl;
         
    // If we've removed markers

    const region_type& r  = my_genome->region(ri->name);

    // Create an output file

    my_ibd_prob_file  = NULL;
    my_ibd_state_file = NULL;

    // Create a pedigree iteration
    for( size_t mp = 0; mp < my_multipedigree->pedigree_count(); ++mp )  
    {
      const RefPedigree& fp = my_multipedigree->pedigree_index(mp);

      process_pedigree(title, output_name, fp, r);
    }

    cout << endl;
  }

  return true;
}

void
genibd_analysis::process_pedigree(const string&      title,
                                  const string&      output,
                                  const RefPedigree& rped,
                                  const region_type& r)
{
  if( !rped.subpedigree_count() ) return;
  
  cout << endl << "  " << title << ": Pedigree " << rped.name() << endl;
  info << endl << "Processing Pedigree " << rped.name() << ".............."
       << endl << endl;

  string name = "Pedigree '" + rped.name() + "'";
  string founders = "";

  filtered_multipedigree fmp0(*my_multipedigree);

  // Do local filtering for this region
  //
  typedef FPED::has_informative_loci<SAGE::RPED::Member> has_inf_loci0;

  has_inf_loci0 hil0(*my_multipedigree, false);

  for( size_t m = 0; m < r.locus_count(); ++m )
  {
    size_t m_index = r.locus(m).marker_index();
    hil0.set_check_status_for_locus(m_index, true);
  }

  // Filter based upon being informative either locally (hil) or within
  // subpedigrees.
  FPED::FilterResults fresults0 = 
      FPED::MPFilterer::add_pedigree_filtered_by_members(fmp0, rped, is_inf_within_sped(hil0));

  fmp0.construct();

  if( !fmp0.pedigree_count() || !fmp0.pedigree_index(0).subpedigree_count() )
  {
    errors << priority(information) << name
           << " has no valid members after removing uninformative "
           << "individuals.  It will be skipped." << endl;

    return;
  }

  const FPED::Pedigree& fped0 = fmp0.pedigree_index(0);

  if( rped.member_count() != fped0.member_count() )
  {
    errors << priority(information)
           << "The following individual(s) were removed as uninformative "
           << "before processing pedigree '" << name
           << "': (";

    FPED::FilterResults::MemberPtrIterator rem_mem0 = fresults0.get_excluded_member_begin();

    errors << (*rem_mem0)->name();

    for( ++rem_mem0; rem_mem0 != fresults0.get_excluded_member_end(); ++rem_mem0 )
      errors << ", " << (*rem_mem0)->name();

    errors << ")." << endl;
  }

  if( fped0.subpedigree_count() > 1 )
    errors << priority(information)
           << "Pedigree '" << rped.name()
           << "' has " << fped0.subpedigree_count()
           << " subpedigrees."
           << "  Each will be processed individually." // except for simulation."
           << endl;

  size_t total_sp = 0;

  for( size_t sp = 0; sp < rped.subpedigree_count(); ++sp )
  {
    const RefSubpedigree& subped = rped.subpedigree_index(sp);

    filtered_multipedigree fmp(*my_multipedigree);
    
    // Do local filtering for this region
    //
    typedef FPED::has_informative_loci<SAGE::RPED::Member> has_inf_loci;

    has_inf_loci hil(*my_multipedigree, false);

    for( size_t m = 0; m < r.locus_count(); ++m )
    {
       size_t m_index = r.locus(m).marker_index();
       hil.set_check_status_for_locus(m_index, true);
    }

    // Filter based upon being informative either locally (hil) or within
    // subpedigrees.
    FPED::FilterResults fresults =
        FPED::MPFilterer::add_subpedigree_filtered_by_members(fmp, subped, is_inf_within_sped(hil));

    fmp.construct();

    if( !fmp.pedigree_count() || !fmp.pedigree_index(0).subpedigree_count() )
      continue;

    const FPED::Pedigree& fped = fmp.pedigree_index(0);

    for( size_t fsp = 0; fsp < fped.subpedigree_count(); ++fsp, ++total_sp )
    {
      const subped_type& fsubped = fped.subpedigree_index(fsp);

      meiosis_map mm(&fsubped, r.is_x_linked());
      
      if( fped0.subpedigree_count() > 1 )
      {
        name = "Subpedigree " + long2str(total_sp + 1) + " of pedigree '";
        name += (rped.name() + "'");

        founders = "founders (" + mm.member(0)->name(); 

        for( size_type f = 1; f < fsubped.member_count(); ++f )
        {
          if( mm.founder(f) )
            founders += (", " + mm.member(f)->name());
        }

        founders += ")";
      }
      
      if( mm.bit_count() < 1 )
      {
        errors << priority(information) << name
               << " has no valid pairs after removing uninformative "
               << "individuals.  It will be skipped." << endl;

        continue;
      }

      if( split_pedigree(mm) )
      {
        errors << priority(information)
               << name << " with " << founders
               << " will be split into nuclear families for this analysis."
               << endl;

        for( size_t fam = 0; fam < fsubped.family_count(); ++fam )
        {
          const family_type& ffamily = fsubped.family_index(fam);

          name = "Nuclear Family (parents " + ffamily.parent1()->name()
                + " and " + ffamily.parent2()->name() + ")";

          cout << endl << "    " << name << ":" << endl << endl;

          filtered_multipedigree fmp2(*my_multipedigree);

          FPED::MPFilterer::add_family(fmp2, ffamily);

          fmp2.construct();

          if( !fmp2.pedigree_count() || !fmp2.pedigree_index(0).subpedigree_count() )
            continue;

          const subped_type& fsubfam = fmp2.pedigree_index(0).subpedigree_index(0);

          meiosis_map mm_fam(&fsubfam, r.is_x_linked());

          analysis_data adata;

          if( check_pedigree(mm_fam, r, adata) )
          {
            process_subpedigree(title, output, name, mm_fam, r, adata);
          }
          else
          {
            errors << priority(information)
                   << name << " cannot be processed on analysis '"
                   << title << "' on region '" << r.name() << "'.  It will be "
                   << "skipped for this analysis." << endl;
          }
        }
      }
      else
      {
        if( fped0.subpedigree_count() > 1 )
        {
          cout << endl;

          cout << "    " << name << " : " << founders << endl << endl;

          errors << priority(information)
                 << name <<" includes " << founders << "." << endl;
        }

        analysis_data adata;

        if( check_pedigree(mm, r, adata) )
        {
          process_subpedigree(title, output, name, mm, r, adata);
        }
        else
        {
          errors << priority(information)
                 << name << " cannot be processed on analysis '"
                 << title << "' on region '" << r.name() << "'.  It will be "
                 << "skipped for this analysis." << endl;
        }
      }
    }
  }

}

void
genibd_analysis::process_subpedigree(const string&         title,
                                     const string&         output,
                                     const string&         name,
                                     const meiosis_map&    mmap,
                                     const region_type&    r,
                                     const analysis_data&  adata)
{
#if 0
  cout << "meiosis_map dump :" << endl;
  mmap.dump_map();
#endif

  ibd_option_type opt;

  opt.title  = title;
  opt.region = r.name();

  opt.max_pedigree = long2str(my_parameters->max_exact_size());

  if( my_parameters->allow_family_splitting() == YES )
    opt.split_pedigrees = "yes";
  else if( my_parameters->allow_family_splitting() == ALWAYS )
    opt.split_pedigrees = "always";

  if( r.is_x_linked() )
    opt.x_linked = true;

  //bool expanded_ibds = false;

  IBD* ibd = NULL;

  if( adata.allow_single )
  {
    if( !do_single_genibd_analysis(title, mmap, r) )
      return;

    ibd = my_pair_ibd_analysis->ibd_adaptor();

    opt.exact = false;
    opt.ibd_mode = "singlepoint";
  }
  else if( adata.allow_exact_multi || adata.allow_exact_single )
  {
    if( !do_exact_genibd_analysis(title, mmap, r) )
      return;

    ibd = my_exact_ibd_analysis->ibd_adaptor();

    opt.exact = true;

    if( my_parameters->scan_interval() )
      opt.scan_type = "intervals";

    if( !my_parameters->is_multipoint() )
    {
      opt.ibd_mode = "singlepoint";
      opt.allow_loops = "on";
    }

    opt.ibd_state_out = my_parameters->output_ibd_state();
  }
  else if( adata.allow_sim )
  {
    if( !do_simulation_genibd_analysis(title, mmap, r) )
      return;

    ibd = my_sim_ibd_analysis->ibd_adaptor();

    opt.exact = false;

    if( my_parameters->allow_simulation() == YES )
      opt.use_simulation = "yes";
    else if( my_parameters->allow_simulation() == ALWAYS )
      opt.use_simulation = "always";

    if( my_parameters->scan_interval() )
      opt.scan_type = "intervals";
    
    if( !my_parameters->is_multipoint() )
    {
      opt.ibd_mode = "singlepoint";
      opt.allow_loops = "on";
    }

    if( my_parameters->scan_interval() )
    {
      ibd = expand_intervals(ibd, mmap);
    }
  }
  else
  {
    errors << priority(information)
           << name << " cannot be processed on analysis '"
           << title << "' on region '" << r.name() << "'.  It will be "
           << "skipped for this analysis." << endl;

    return;
  }

  if( !ibd || !ibd->built() || !ibd->pair_count() )
    return;

  ibd->set_ibd_option(opt);

  if( !my_ibd_prob_file )
  {
    my_ibd_prob_file = new RefIBDWriteFile(output + ".ibd", std::cerr);
    my_ibd_prob_file->output_probability_header(ibd);
  }

  bool valid_ibd_file = my_ibd_prob_file->output_ibd_probability(ibd);

  if(    valid_ibd_file
      && my_parameters->output_ibd_state()
      && ( adata.allow_exact_multi || adata.allow_exact_single) )
  {
    if( !my_ibd_state_file )
    {
      my_ibd_state_file = new RefIBDWriteFile(output + ".state", std::cerr);
      my_ibd_state_file->output_state_header(ibd);
    }

    my_ibd_state_file->output_ibd_state(ibd);
  }

  return;
}


bool
genibd_analysis::do_single_genibd_analysis(const string&      title,
                                           const meiosis_map& mmap,
                                           const region_type& r)
{
  if( my_pair_ibd_analysis )
    delete my_pair_ibd_analysis;

  my_pair_ibd_analysis = new pair_ibd_analysis();

  my_pair_ibd_analysis->set_errors(errors);

  if( !my_pair_ibd_analysis->build() )
    return false;

  my_ped_region.build(*mmap.get_subpedigree(), r, true);

  my_pair_ibd_analysis->set_pedigree(mmap.get_subpedigree(), my_ped_region);

  my_pair_ibd_analysis->build_ibds();

  if( !my_pair_ibd_analysis->built() )
    return false;

  add_pairs_single(*(mmap.get_subpedigree()));

  my_pair_ibd_analysis->compute(title);

#if 0
  vector<double>  f0, f1, f2;
  IBD* ibd = my_pair_ibd_analysis->ibd_adaptor();

  for( size_t i = 0; i < ibd->pair_count(); ++i )
  {
    cout << i << " ";

    ibd->get_ibd(i, f0, f1, f2);
    for( size_t m = 0; m < f0.size(); ++m )
    {
      cout << f0[m] << ", " << f1[m] << ", " << f2[m] << "	";
    }
    cout << endl;
  }
#endif

  // Just checking...
  if( !my_pair_ibd_analysis->valid() )
    return false;

  return true;
}

bool
genibd_analysis::do_exact_genibd_analysis (const string&      title,
                                           const meiosis_map& mmap,
                                           const region_type& r)
{
  size_t max_loci = r.locus_count();
  size_t max_bits = mmap.bit_count();

  if( my_exact_ibd_analysis )
    delete my_exact_ibd_analysis;

  my_exact_ibd_analysis = new exact_ibd_analysis(sage_cerr, cout, true);

  my_exact_ibd_analysis->set_errors(errors);

  if( !my_exact_ibd_analysis->build(max_loci, max_bits,
                                    my_parameters->is_multipoint(),
                                    !my_parameters->is_multipoint(),
                                    my_parameters->scan_interval()) )
    return false;

  my_ped_region.build(*mmap.get_subpedigree(), r, true);

  my_exact_ibd_analysis->set_pedigree(mmap, my_ped_region);

  bool use_intervals =    my_parameters->is_multipoint()
                       && my_parameters->scan_interval();

  my_exact_ibd_analysis->build_ibds(use_intervals);

  if( !my_exact_ibd_analysis->built() )
    return false;

  add_pairs_exact(*(mmap.get_subpedigree()));

  my_exact_ibd_analysis->compute(title, !my_parameters->is_multipoint(),
                                 my_parameters->scan_interval(),
                                 my_parameters->output_ibd_state());

#if 0
  vector<double>  f0, f1, f2;
  IBD* ibd = my_exact_ibd_analysis->ibd_adaptor();

  for( size_t i = 0; i < ibd->pair_count(); ++i )
  {
    cout << i << " ";

    ibd->get_ibd(i, f0, f1, f2);
    for( size_t m = 0; m < f0.size(); ++m )
    {
      cout << f0[m] << ", " << f1[m] << ", " << f2[m] << "	";
    }
    cout << endl;
  }
#endif

  // Just checking...
  if( !my_exact_ibd_analysis->valid() )
    return false;

  return true;
}

bool
genibd_analysis::do_simulation_genibd_analysis (const string&         title,
                                                const meiosis_map& mmap,
                                                const region_type&    r)
{
  bool codominant = true;

  my_ped_region.build(*mmap.get_subpedigree(), r, true);

  // Check codominance of the pedigree_region.
  for( size_t i = 0; i != r.locus_count(); ++i )
  {
    if(    !my_ped_region.model_informative(i) // same as good() in old code???
        && !my_ped_region[i].codominant() )
    {
      codominant = false;

      errors << priority(information) << "Marker '"
             << my_ped_region[i].name() << "' is not codominant for pedigree '"
             << mmap.get_pedigree()->name() << "' on region '"
             << my_ped_region.get_region().name() << "'.  It cannot be "
             << "simulated." << endl;
    }
  }

  if( !codominant )
    return false;

  if( my_sim_ibd_analysis )
    delete my_sim_ibd_analysis;

  my_sim_ibd_analysis = new sim_ibd_analysis(my_parameters->get_sim_parameters());

  my_sim_ibd_analysis->set_errors(errors);

  if( !my_sim_ibd_analysis->build() )
    return false;

  my_sim_ibd_analysis->set_pedigree(mmap, my_ped_region);

  my_sim_ibd_analysis->build_ibds();

  assert(my_sim_ibd_analysis->built());

  add_pairs_simulation(*(mmap.get_subpedigree()));

  return my_sim_ibd_analysis->compute(info);
}

void
genibd_analysis::add_pairs_exact(const subped_type& subped)
{
  generate_pairs(subped);

  //IBD* ibd = my_exact_ibd_analysis->ibd_adaptor();
  
  relpair_set_type::const_iterator t = my_relpair_set.begin();
  for( ; t != my_relpair_set.end(); ++t )
  {
    const pair_generator::relative_pair& rel_pair = *t;

    fmember_const_pointer member_one = subped.pedigree()->member_find(rel_pair.member_one()->name());
    fmember_const_pointer member_two = subped.pedigree()->member_find(rel_pair.member_two()->name());

    if( !member_one || !member_two )
      continue;

    if( both_members_exist(rel_pair, subped) )
    {
      if(    rel_pair.type() == pair_generator::SIBSIB
          || rel_pair.type() == pair_generator::HALFSIB
          || rel_pair.type() == pair_generator::COUSIN )
      {
        if( member_one->name() > member_two->name() )
          my_exact_ibd_analysis->add_pair(member_two, member_one, rel_pair.type());
        else
          my_exact_ibd_analysis->add_pair(member_one, member_two, rel_pair.type());
      }
      else
        my_exact_ibd_analysis->add_pair(member_one, member_two, rel_pair.type());

#if 0
  cout << pair_count << " "
       << rel_pair.member_one()->pedigree()->name() << ":("
       << rel_pair.member_one()->name() << ", "
       << rel_pair.member_two()->name() << ") "
       << endl;
#endif
    }
  }
}

void
genibd_analysis::add_pairs_single(const subped_type& subped)
{
  generate_pairs(subped);

  relpair_set_type::const_iterator t = my_relpair_set.begin();
  for( ; t != my_relpair_set.end(); ++t )
  {
    const pair_generator::relative_pair& rel_pair = *t;

    fmember_const_pointer member_one = subped.pedigree()->member_find(rel_pair.member_one()->name());
    fmember_const_pointer member_two = subped.pedigree()->member_find(rel_pair.member_two()->name());
    fmember_const_pointer connector_one = NULL;
    fmember_const_pointer connector_two = NULL;

    if( rel_pair.connector_one() != NULL )
      connector_one = subped.pedigree()->member_find(rel_pair.connector_one()->name());
    if( rel_pair.connector_two() != NULL )
      connector_two = subped.pedigree()->member_find(rel_pair.connector_two()->name());

    if( !member_one || !member_two )
      continue;

    if( both_members_exist(rel_pair, subped) )
    {
      if(    rel_pair.type() == pair_generator::SIBSIB
          || rel_pair.type() == pair_generator::HALFSIB )
      {
        if( member_one->name() > member_two->name() )
          my_pair_ibd_analysis->add_pair(member_two, member_one, connector_one, connector_two, rel_pair.type());
        else
          my_pair_ibd_analysis->add_pair(member_one, member_two, connector_one, connector_two, rel_pair.type());
      }
      else if( rel_pair.type() == pair_generator::COUSIN )
      {
        if( member_one->name() > member_two->name() )
          my_pair_ibd_analysis->add_pair(member_two, member_one, connector_two, connector_one, rel_pair.type());
        else
          my_pair_ibd_analysis->add_pair(member_one, member_two, connector_one, connector_two, rel_pair.type());
      }
      else
        my_pair_ibd_analysis->add_pair(member_one, member_two, connector_one, connector_two, rel_pair.type());

#if 0
  cout  << pair_gen.pair_type_to_string(rel_pair.type()) << ": "
       << rel_pair.member_one()->pedigree()->name() << ":("
       << rel_pair.member_one()->name() << ", "
       << rel_pair.member_two()->name() << ") (";
  if( rel_pair.connector_one() != NULL )
    cout << rel_pair.connector_one()->name() << ", ";
  else
    cout << "?, ";
  if( rel_pair.connector_two() != NULL )
    cout << rel_pair.connector_two()->name() << ") ";
  else
    cout << "?) ";
  cout << endl;
#endif
    }
  }
}

void
genibd_analysis::add_pairs_simulation(const subped_type& subped)
{
  generate_pairs(subped);

  relpair_set_type::const_iterator t = my_relpair_set.begin();
  for( ; t != my_relpair_set.end(); ++t )
  {
    const pair_generator::relative_pair& rel_pair = *t;

    fmember_const_pointer member_one = subped.pedigree()->member_find(rel_pair.member_one()->name());
    fmember_const_pointer member_two = subped.pedigree()->member_find(rel_pair.member_two()->name());

    if( !member_one || !member_two )
      continue;

    if( both_members_exist(rel_pair, subped) )
    {
      if(    rel_pair.type() == pair_generator::SIBSIB
          || rel_pair.type() == pair_generator::HALFSIB
          || rel_pair.type() == pair_generator::COUSIN )
      {
        if( member_one->name() > member_two->name() )
          my_sim_ibd_analysis->add_pair(member_two, member_one, rel_pair.type());
        else
          my_sim_ibd_analysis->add_pair(member_one, member_two, rel_pair.type());
      }
      else
        my_sim_ibd_analysis->add_pair(member_one, member_two, rel_pair.type());

#if 0
  cout  << pair_gen.pair_type_to_string(rel_pair.type()) << ": "
       << rel_pair.member_one()->pedigree()->name() << ":("
       << rel_pair.member_one()->name() << ", "
       << rel_pair.member_two()->name() << ") ";
  cout << endl;
#endif
    }
  }
}

void
genibd_analysis::generate_pairs(const subped_type& subped)
{
  unsigned int type_wanted = 0;
  
  switch(my_parameters->pair_category())
  {
    case SIB         : type_wanted = pair_generator::SIBSIB_MASK; 
                       break;
    case ALL_SIB     : type_wanted = pair_generator::SIBSIB_MASK |
                                     pair_generator::HALFSIB_MASK;
                       break;
    case RELATIVE    : type_wanted = pair_generator::SIBSIB_MASK |
                                     pair_generator::GRANDP_MASK |
                                     pair_generator::AVUNC_MASK |
                                     pair_generator::HALFSIB_MASK |
                                     pair_generator::COUSIN_MASK;
                       break;
    case ALL         : type_wanted = pair_generator::EVERY_MASK;
                       break;
  }

  const RefPedigree* rp = my_multipedigree->pedigree_find(subped.pedigree()->name());

  pair_generator pair_gen(const_cast<RefPedigree*>(rp), type_wanted);

  my_relpair_set.clear();

  pair_generator::iterator pi = pair_gen.begin();
  for( ; pi != pair_gen.end(); ++pi )
  {
    pair_generator::relative_pair& rel_pair = *pi;

    my_relpair_set.insert(rel_pair);
  }
}

IBD*
genibd_analysis::expand_intervals
    (IBD* old_ibd,
     const meiosis_map& mmap)
{
  region_type r = my_ped_region.get_region();

  IBD* new_ibd = new basic_storage_ibd(mmap, r, true);

  // Add the pairs

  for(size_t i = 0; i < old_ibd->pair_count(); ++i)
  {
    id_pair p = old_ibd->get_pair(i);

    new_ibd->add_pair(p.first, p.second);
  }

  new_ibd->build();

  // Copy the IBDs we calculated

  for(size_t i = 0; i < old_ibd->pair_count(); ++i)
  {
    for(size_t m = 0; m < old_ibd->marker_count(); ++m)
    {
      double f0, f2;
      
      old_ibd->get_ibd(i, m, f0, f2); 

      new_ibd->set_ibd(i, r.locus(m).point_prev_count(), f0, f2);
    }
  }

  return new_ibd;
}

} // end of namespace GENIBD

} // end of namespace SAGE

