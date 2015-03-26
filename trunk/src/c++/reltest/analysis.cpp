//==========================================================================
//  File:       analysis.cpp
//
//  Author:     Qing Sun & Yeunjoo Song
//
//  History:    Version 1.0  Initial implementation.
//                      2.0  Updated to new libraries            yjs Jul. 03
//
//  Notes:      This class is the main driver performing analysis.
//
//  Copyright (c) 2003 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "reltest/output.h"
#include "maxfunapi/maxfunapi.h"

namespace SAGE
{

namespace RELTEST
{

reltest_analysis::reltest_analysis(reltest_data& rtfile, const reltest_parser& rp, cerrorstream& err)
{
  if( !(&rtfile) )
    return;

  my_exact_ibd_analysis = NULL;
  my_input  = &rtfile;
  my_parser = &rp;

  my_multipedigree = &(my_input->pedigrees());

  my_total_map_points = 0.;
  my_total_genome_length = 0.;

  my_current_pairtype = SIB;

  my_adjusted_cutpoints[Cu] = my_cutpoints[Cu] = 0.;
  my_adjusted_cutpoints[Ch] = my_cutpoints[Ch] = 0.;
  my_adjusted_cutpoints[Cm] = my_cutpoints[Cm] = 3.27;
  my_adjusted_cutpoints[Cp] = my_cutpoints[Cp] = 0.;

  my_AMIC = 0.;

  my_unexpected_error = false; 
  errors = err;
}

reltest_analysis::~reltest_analysis()
{
  if( my_exact_ibd_analysis )
    delete my_exact_ibd_analysis;
}

bool
reltest_analysis::run_analysis(ostream& sum_out, ostream& nuc_out, ostream& det_out)
{
  if( !my_input )
    return false;

  build();

  if( my_unexpected_error || !do_analysis() )
    return false;

  reltest_output out_writer(this, sum_out, nuc_out, det_out, errors);

  return true;
}

void reltest_analysis::dump_pairs(ostream& o)  const
{
  for( size_t i = 0; i < my_ptt_pairs.size(); ++i )
    dump_current_pairtype(putative_type(i), o);
}

void reltest_analysis::dump_current_pairtype(putative_type p, ostream& o) const
{
  o << endl 
    << "Putative Type : " << get_pair_name(p) << endl
    << "         Size : " << my_ptt_pairs[p].size();
  
  for( size_t i = 0; i < my_ptt_pairs[p].size(); ++i )
  {
    if( !(i % 5) )
      o << endl << "               ";

    putative_pair pair = my_ptt_pairs[p][i];

    o << setw(2) << pair.first()->pedigree()->name() << " : "
      << setw(2) << pair.first()->name() << " - "
      << setw(2) << pair.second()->name() << "  ";
  }
  o << endl;
}

//
//-----------------------------------------------------------
//

bool
reltest_analysis::build()
{
  LSF_ptr<RPED::genome_description> mg = my_input->genome();

  my_ptt_pairs.resize(5);
  my_misc_pairs.resize(5);

  my_chrom_length.resize(mg->region_count());

  for( int c = 0; c < mg->region_count(); ++c )
    my_chrom_length[c] = mg->region(c).length();

  const vector<size_t>& analysis_regions = my_parser->get_analysis_regions();

  for( size_t k = 0; k < analysis_regions.size(); ++k )
  {
    my_total_genome_length += my_chrom_length[analysis_regions[k]];

    int invalid_markers = 0;  // number of markers not on 1 cM boundaries

    const RPED::genome_description::region_type& r = mg->region(analysis_regions[k]);

    for( size_t l = 0; l < r.locus_count(); ++l )
    {
      if( !test_location(r.locus_location(l)) )
      {
        ++invalid_markers;
      }
    }

    my_total_map_points += mg->region(analysis_regions[k]).point_count() - invalid_markers;
  }

#if 0
  cout << "my_total_genome_length = " << my_total_genome_length << endl
       << "my_total_map_points    = " << my_total_map_points << endl;
#endif

  if( my_total_genome_length < 1200 )
    errors << priority(warning) 
           << "Total length of genotyped genome is less than 1200cM."
           << " Pair types may not be easily distinguished." << endl;


  return true;
}

//*********************************************************************
//**
//** For each type of putative pairs, we do the following procdure:
//**   1. generate the pairs
//**   2. perform our analysis on the pairs over the whole genome
//**
//** Note: the statitics used to calculate cut points are calculated
//**       once during the phase of the analysis of putative sib pairs.
//**
//**********************************************************************   
bool
reltest_analysis::do_analysis()
{
  // generate pairs
  generate_pairs();

  for( size_t i = 0; i < 5; ++i ) //i=0:sib, 1:hsib, 2:MZtwin, 3:parent_child, 4:parents
  {
    // If putative sib pair type is not required by user and any of cutpoints 
    // is not preset, then we need to do our analysis on putative sib pair 
    // in order to obtain complete set of cutpoints regardless it is set or not.
    if(    i == 0 
        && !my_parser->get_analysis_pairtypes()[i]
        && !my_parser->calculate_cutpoints() )
      continue;

    // We only do our analyses on the types of putative pairs required, except
    // putative sib pair which is handled above.
    if( i && !my_parser->get_analysis_pairtypes()[i] )
      continue;

    my_current_pairtype = (putative_type) i;

    // If any individual in a pair have no data over whole genome,
    // then the pair is considered as uninformative pair and is removed.
    remove_uninformative_pairs();

    // No sib pairs available
    if(    i == 0
        && my_ptt_pairs[my_current_pairtype].size() == 0 )
    {
      errors << priority(information) 
             << "Dataset contains no informative sib pairs.  Analysis will not run." 
             << endl; 
      my_unexpected_error = true;

      return false;
    }

    // No pair of type i available
    if( my_ptt_pairs[my_current_pairtype].size() == 0 )
    {
      errors << priority(information)
             << "Dataset contains no informative "
             << get_pair_name(my_current_pairtype)
             << "." << endl;
      continue;
    }

    // Screen out
    string p_name = toLower(get_pair_name(my_current_pairtype));
    string dot = ".........................";

    string sub_dot = dot.substr(0, 34 - p_name.size());
    p_name += sub_dot;

    cout << "Testing " << setw(22) << p_name << flush;

    // Now, we perform our analysis on the current putative pair.
    pair_analysis();

    // For each analysis run, we need to compute our statistics based on IBDs of
    // putative sib pairs and we only need to do this once.
    if( i == 0 )
    {
      stats();

      // Before reclassfy the pairs, apply a nonparametric estimation procedure
      // to obtain the mean and variance of the sib-pair distributions of
      // our statistics. - new addition, 09-2001, yjs
      //
      nonparametric_estimation();
    }

    // Reclassify our current putative pairs based on their Yj and Yj*
    reclassify_pairs();

    // Find out the pairs that are misclassified.
    get_misclassified_pairs();

    // Test the deviation of mean from zero. - new addition 09-2001 yjs
    if( i == 0 )
      test_deviation();

    cout << "done." << flush << endl;
  }

  if( my_unexpected_error )
    return false;
  else
    return true;
}

string reltest_analysis::get_pair_name(putative_type type) const
{
  if( type==SIB )                     return "SIB PAIRS";
  else if( type==HSIB )               return "HALF SIB PAIRS";
  else if( type==MZTWIN )             return "MZTWINS";
  else if( type==MARITAL )            return "MOTHER_FATHER PAIRS"; 
  else if( type==PARENT_OFFSPRING )   return "PARENT_OFFSPRING PAIRS";
  else                                return "UNRELATED PAIRS";
}

void
reltest_analysis::generate_pairs()
{
  if( !my_multipedigree )
  {
    my_unexpected_error = true;
    return;
  }

  RPED::RefMultiPedigree::pedigree_const_iterator p = my_multipedigree->pedigree_begin();

  for( ; p != my_multipedigree->pedigree_end(); ++p )  
  {
    RPED::RefMultiPedigree::subpedigree_const_iterator sp = p->subpedigree_begin();
    for( ; sp != p->subpedigree_end(); ++sp )
    {
      if( !sp->family_count() )
        continue;

      RPED::RefMultiPedigree::family_const_iterator f = sp->family_begin();

      for( ; f != sp->family_end(); ++f )
      {
        // marital pair
        putative_pair parent_pair;

        parent_pair.first(f->parent1());
        parent_pair.second(f->parent2());

        my_ptt_pairs[MARITAL].push_back(parent_pair);

        RPED::RefMultiPedigree::offspring_const_iterator o1 = f->offspring_begin();

        for( ; o1 != f->offspring_end(); ++o1 )
        {
          // parent-offspring pair
          putative_pair po_pair1;

          po_pair1.first(f->parent1());
          po_pair1.second(&(*o1));

          my_ptt_pairs[PARENT_OFFSPRING].push_back(po_pair1);

          putative_pair po_pair2;

          po_pair2.first(f->parent2());
          po_pair2.second(&(*o1));

          my_ptt_pairs[PARENT_OFFSPRING].push_back(po_pair2);

          // sib pair
          RPED::RefMultiPedigree::offspring_const_iterator o2 = o1;
          ++o2;

          for( ; o2 != f->offspring_end(); ++o2 )
          {
            putative_pair sib_pair;

            sib_pair.first(&(*o1));
            sib_pair.second(&(*o2));

            my_ptt_pairs[SIB].push_back(sib_pair);
          }
        }
      }

      // half sib pair
      size_t member_count = p->member_count();

      for( size_t i1 = 0; i1 < member_count; ++i1 )
      {
        for( size_t i2 = 0; i2 <= i1; ++i2 )
        {
          ind_id s1 = &p->member_index(i1);
          ind_id s2 = &p->member_index(i2);

          if( is_hsib(s1, s2) )
          {
            putative_pair hsib_pair;

            hsib_pair.first(s1);
            hsib_pair.second(s2);

            my_ptt_pairs[HSIB].push_back(hsib_pair);
          }
        }
      }
    }
  }
}

//uninformative pair is the pair that contains at least one individual
//having no data at all markers over the whole genome.
// - Not whole genome, but regions included in analysis
//
void
reltest_analysis::remove_uninformative_pairs()
{
  //calculate total number of markers over the whole genome
  // ?? not whole genome, but just analysis regions.

//  double total_markers = my_input->genome()->locus_count();

  double total_markers = 0.;

  const vector<size_t>& analysis_regions = my_parser->get_analysis_regions();

  for( size_t k = 0; k < analysis_regions.size(); ++k )
  {
    const region& r = my_input->genome()->region(analysis_regions[k]);

    total_markers += r.locus_count();
  }

  double m1,m2; 
  vector<putative_pair>::iterator i = my_ptt_pairs[my_current_pairtype].begin();

  for( ; i != my_ptt_pairs[my_current_pairtype].end(); )
  {
    m1 = uninformative_markers(i->first());

    if(    m1 == total_markers
       || (m1 / total_markers) * 100. >= 99. )
    {
      i = my_ptt_pairs[my_current_pairtype].erase(i);
      continue;
    }

    m2 = uninformative_markers(i->second());

    if(    m2 == total_markers
       || (m2 / total_markers) * 100. >= 99. )
    {
      i = my_ptt_pairs[my_current_pairtype].erase(i);
      continue;
    }
   
    i->set_s1_missing_rate(m1/total_markers);
    i->set_s2_missing_rate(m2/total_markers);
    ++i;
  }
}

int
reltest_analysis::uninformative_markers(ind_id id)  const
{
  const RPED::RefMPedInfo& mped_info = my_multipedigree->info();
  const RPED::RefPedInfo&  ped_info  = id->pedigree()->info();

  int missing_count = 0;

  const vector<size_t>& analysis_regions = my_parser->get_analysis_regions();

  for( size_t k = 0; k < analysis_regions.size(); ++k )
  {
    const region& r = my_input->genome()->region(analysis_regions[k]);

    for( size_t l = 0; l < r.locus_count(); ++l )
    {
      size_t m = r.locus(l).marker_index();

      if( ped_info.phenotype_missing(id->index(), m, mped_info.marker_info(m)) )
        ++missing_count;
    }
  }

/* Not whole genome, but just the regions included in analysis.
  for( size_t c = 0; c < my_input->genome()->region_count(); ++c )
  {
    RPED::genome_description::region_type rt = my_input->genome()->region(c);

    int markers = rt.locus_count();
   
    for( size_t m = 0; m < markers; ++m )
    {
      size_t mi = rt.locus(m).marker_index();

      if( ped_info.phenotype_missing(id->index(), mi, mped_info.marker_info(mi)) )
        ++missing_count;
    }
  }
*/
  return missing_count;
}

void
reltest_analysis::pair_analysis()
{
  IBD* ibd = NULL;

  //for each chromosome(region), generate IBD probabilities for each putative pair
  //
  const vector<size_t>& analysis_regions = my_parser->get_analysis_regions();

  for( size_t c = 0; c < analysis_regions.size(); ++c )
  {
#if 0
    cout << "Processiong region " << c+1 << " of " << analysis_regions.size()
         << " with " << my_ptt_pairs[my_current_pairtype].size() << " pairs..." << endl;
#endif
    const region& r = my_input->genome()->region(analysis_regions[c]);

    for( size_t i = 0; i < my_ptt_pairs[my_current_pairtype].size(); ++ i )
    {
      const putative_pair& pair = my_ptt_pairs[my_current_pairtype][i];

      const subpedigree& subped = *(pair.first()->subpedigree());

#if 0
      // build temp multipedigree
      //
      RPED::RefMultiPedigree temp_mp;

      for( size_t c = 0; c < my_input->genome()->region_count(); ++c )
      {
        const region& r = my_input->genome()->region(c);

        for( size_t m = 0; m < r.locus_count(); ++m )
        {
          size_t mi = r.locus(m).marker_index();

          RPED::RefMarkerInfo rm(subped.multipedigree()->info().marker_info(mi));

          temp_mp.info().add_marker(rm.name(), rm);
        }
      }

      temp_mp.info().set_individual_missing_code(subped.multipedigree()->info().individual_missing_code());
      temp_mp.info().set_sex_code_male(subped.multipedigree()->info().sex_code_male());
      temp_mp.info().set_sex_code_female(subped.multipedigree()->info().sex_code_female());
      temp_mp.info().set_sex_code_unknown(subped.multipedigree()->info().sex_code_unknown());

      string ped_name = subped.pedigree()->name();

      temp_mp.add_member(ped_name, "dummy2", SAGE::male);
      temp_mp.add_member(ped_name, "dummy1", SAGE::female);

      temp_mp.add_member(ped_name, pair.first()->name(), pair.first()->gender());
      temp_mp.add_lineage(ped_name, pair.first()->name(), "dummy1", "dummy2");

      temp_mp.add_member(ped_name, pair.second()->name(), pair.second()->gender());
      temp_mp.add_lineage(ped_name, pair.second()->name(), "dummy1", "dummy2");

      temp_mp.build();
      
      RPED::RefMultiPedigree::pedigree_iterator p = temp_mp.pedigree_begin();
      for( ; p != temp_mp.pedigree_end(); ++p )
      {
        p->info().build(*p);
        p->info().resize_markers( temp_mp.info().marker_count() );
        p->info().resize_traits(0);
        p->info().resize_strings(0);
      }

      RPED::RefMultiPedigree::member_const_pointer sib1 = temp_mp.member_find(ped_name, pair.first()->name());
      RPED::RefMultiPedigree::member_const_pointer sib2 = temp_mp.member_find(ped_name, pair.second()->name());

      RPED::RefMPedInfo&                    mped_info = temp_mp.info();
      RPED::RefMultiPedigree::pedinfo_type& info      = temp_mp.pedigree_find(ped_name)->info();

      for( size_t c = 0; c < my_input->genome()->region_count(); ++c )
      {
        const region& r = my_input->genome()->region(c);

        for( size_t m = 0; m < r.locus_count(); ++m )
        {
          size_t mi = r.locus(m).marker_index();

          const RPED::RefMarkerInfo rm = subped.multipedigree()->info().marker_info(mi);

          size_t new_mi = mped_info.marker_find(rm.name());

          info.set_phenotype(sib1->index(), new_mi, subped.pedigree()->info().phenotype(pair.first()->index(), mi));
          info.set_phenotype(sib2->index(), new_mi, subped.pedigree()->info().phenotype(pair.second()->index(), mi));
        }
      }

      // Sort pedigree
      for( p  = temp_mp.pedigree_begin(); p != temp_mp.pedigree_end(); ++p )
        PedigreeSort(*p);

      const subpedigree& new_subped = temp_mp.pedigree_find(ped_name)->subpedigree_index(0);
      putative_pair      new_pair;
      new_pair.first(sib1);
      new_pair.second(sib2);

      ibd = get_pair_ibds(r, new_subped, new_pair);

      cout << endl;
      cout << "new RPED::RefMultipedigree :" << endl;
      cout << "                info :" << endl;
      const RPED::RefMPedInfo& m_info = temp_mp.info();
      cout << "  trait count  = " << m_info.trait_count() << endl;
      cout << "  marker count = " << m_info.marker_count() << endl;
      for( size_t mi = 0; mi < m_info.marker_count(); ++mi )
      {
        cout << "          " << mi << " : "
             << m_info.marker_info(mi).name() << "	"
             << m_info.marker_info(mi).allele_count() << "	"
             << m_info.marker_info(mi).phased_genotype_count() << "	"
             << m_info.marker_info(mi).unphased_genotype_count() << "	"
             << m_info.marker_info(mi).phenotype_count() << endl;
      }

      RPED::RefMultiPedigree::pedigree_const_iterator pi = temp_mp.pedigree_begin();
      for( ; pi != temp_mp.pedigree_end(); ++pi )
      {
        cout << pi->index() << ": " << pi->name() << endl;
        cout << "   mem count : " << pi->member_count() << endl;

        RPED::RefMultiPedigree::family_const_iterator fi = pi->family_begin();
        for( ; fi != pi->family_end(); ++fi )
        {
          cout << fi->index() << ": " << fi->name() << endl;
          cout << "p1 = " << fi->parent1()->name() << endl;
          cout << "p2 = " << fi->parent2()->name() << endl;
          RPED::RefMultiPedigree::offspring_const_iterator oi = fi->offspring_begin();
          for( ; oi != fi->offspring_end(); ++oi )
          {
            cout << "offspring " << oi->index() << " : " << oi->name() << endl;
          }
        }
      }
#endif

      FPED::Multipedigree fped(*subped.multipedigree());

      FPED::filter_to_sib_pair(fped, *(pair.first()), *(pair.second()));
      
      fped.construct();

#if 0
  cout << endl;
  cout << "filtered_multipedigree :" << endl;
  cout << "                  info :" << endl;
  const RPED::RefMPedInfo& f_info = fped.info();
  cout << "  trait count  = " << f_info.trait_count() << endl;
  cout << "  marker count = " << f_info.marker_count() << endl;

  for( size_t mi = 0; mi < f_info.marker_count(); ++mi )
  {
    cout << "          " << mi << " : "
         << f_info.marker_info(mi).name() << "	"
         << f_info.marker_info(mi).allele_count() << "	"
         << f_info.marker_info(mi).phased_genotype_count() << "	"
         << f_info.marker_info(mi).unphased_genotype_count() << "	"
         << f_info.marker_info(mi).phenotype_count() << endl;
  }

  fpedigree_const_iterator fpi = fped.pedigree_begin();
  for( ; fpi != fped.pedigree_end(); ++fpi )
  {
    cout << fpi->index() << ": " << fpi->name() << endl;
    cout << "   mem count : " << fpi->member_count() << endl;

    ffamily_const_iterator ffi = fpi->family_begin();
    for( ; ffi != fpi->family_end(); ++ffi )
    {
      cout << ffi->index() << ": " << ffi->name() << endl;
      cout << "p1 = " << ffi->parent1()->name() << endl;
      cout << "p2 = " << ffi->parent2()->name() << endl;
      foffspring_const_iterator foi = ffi->offspring_begin();
      for( ; foi != ffi->offspring_end(); ++foi )
      {
        cout << "offspring " << foi->index() << " : " << foi->name() << endl;
      }
    }
  }
#endif

      const FPED::Subpedigree& new_subped = fped.pedigree_index(0).subpedigree_index(0);
      
      ibd = get_pair_ibds(r, new_subped, pair);
      
      if( !ibd ) 
      {
        my_unexpected_error = true;
        return;
      }

      if( !update_pair_stats(i, ibd, r) )  
      { 
        my_unexpected_error = true;
        return; 
      }
    }
  }
}

IBD*
reltest_analysis::get_pair_ibds(const region& r, const FPED::Subpedigree& subped, const putative_pair& pair)
{
  size_t max_loci = r.locus_count();

  meiosis_map mm(&subped);

  size_t max_bits = mm.bit_count();

  pedigree_region  pr(subped, r, errors, true, true);

  if( my_exact_ibd_analysis )
    delete my_exact_ibd_analysis;
  
  my_exact_ibd_analysis = new exact_ibd_analysis();

  my_exact_ibd_analysis->set_errors(errors);

  if( !my_exact_ibd_analysis->build(max_loci, max_bits) )
    return NULL;

  my_exact_ibd_analysis->set_pedigree(mm, pr);

  my_exact_ibd_analysis->build_ibds(true);

  assert(my_exact_ibd_analysis->built());

  FPED::MemberConstPointer sib1 = subped.pedigree()->member_find(pair.first()->name());
  FPED::MemberConstPointer sib2 = subped.pedigree()->member_find(pair.second()->name());

  my_exact_ibd_analysis->add_pair(sib1, sib2, pair_generator::SIBSIB);
  
  my_exact_ibd_analysis->compute("test");

  assert(my_exact_ibd_analysis->valid());

  return my_exact_ibd_analysis->ibd_adaptor();
}

bool
reltest_analysis::update_pair_stats(int pair_index, IBD* ibd, const region& r)
{
  // ibd values for a chromosome(region) k.
  //
  vector<double>  f0, f2;

  putative_pair& p = my_ptt_pairs[my_current_pairtype][pair_index];

  string  n1,n2;
  string  pn;

  n1 = p.first()->name();
  n2 = p.second()->name();
  pn = p.first()->pedigree()->name();

//  size_t p_index = ibd->pair_index(p.first(), p.second());

//  ibd->get_ibd(p_index, f0, f2);
  ibd->get_ibd(0, f0, f2);
  
#if 0
  cout << "pair_index = " << pair_index << endl;
  cout << "point      = " << f0.size() << endl;
#endif

  double points = f0.size();

  double Yj  = 0.0;
  double Yjp = 0.0;
  double MIC = 0.0;

  for( size_t i = 0; i < f0.size(); ++i )
  {
    double d = r.point_location(i);

#if 0
  cout << "i = " << i << ", ";
  cout << f0[i] << ", " << f2[i] << ", d = " << d;
#endif

    if( !test_location(d) )
    {
      --points;
      continue;
    }

    //calculate mean allele sharing Yj and parents-offspring statistic Yjp

    double f_diff = f2[i] - f0[i];

    Yj  += f_diff;                        // Xs - 1 = f1+f2+f2-1 = (1-f0-f2)+f2+f2-1 = f2-f0 
    Yjp += 2.0 * (f2[i] + f0[i]) - 1.0;   // Xsp    = f2+f0-f1   = f2+f0-(1-f0-f2)   = 2(f0+f2)-1
    MIC += f2[i] * 3.0 - f0[i] + 1.0 - (f_diff + 1.0)*(f_diff + 1.0); //MIC = ?

#if 0
  cout << ", Yj = " << Yj << ", Yjp = " << Yjp << ", MIC = " << MIC << endl;
#endif
  }

  Yj  = Yj * sqrt(2.0) / points; 
  Yjp = Yjp            / points;

  p.add_Yj (Yj);
  p.add_Yjp(Yjp);
  p.add_MIC(MIC);

#if 0
  cout << "final Yj = " << Yj << ", Yjp = " << Yjp << ", MIC = " << MIC << endl;
#endif

  return true;
}

void reltest_analysis::stats()
{
  double sum_Var_Yj  = 0.0;
  double sum_Var_Yjp = 0.0;

  const vector<size_t>& analysis_regions = my_parser->get_analysis_regions();

  //calculate total genome length over all chromosomes and Var_Yj, Var_Yjp
  //
  for( size_t k = 0; k < analysis_regions.size(); ++k )
  {
    // Lk is the length of chromosome k in cM.
    // beta = 0.04 for Yj.
    double beta_Lk  = 0.04 * my_chrom_length[analysis_regions[k]];

    sum_Var_Yj  += (2.0 / beta_Lk) - (2.0 * (1.0 - exp(-beta_Lk)) / (beta_Lk*beta_Lk));

    // beta = 0.08 for Yj*.
    beta_Lk *= 2.0;
    sum_Var_Yjp += (2.0 / beta_Lk) - (2.0 * (1.0 - exp(-beta_Lk)) / (beta_Lk*beta_Lk));
  }

  my_Var_Yj  = sqrt(sum_Var_Yj);
  my_Var_Yjp = sqrt(sum_Var_Yjp);

  double temp = 0.0;
  for( size_t i = 0; i < my_ptt_pairs[my_current_pairtype].size(); ++i )
  {
    temp += my_ptt_pairs[my_current_pairtype][i].get_MIC(); 
  }

  my_AMIC = 1.0 - (temp*2.0)/(my_total_map_points*my_ptt_pairs[my_current_pairtype].size());

#if 0
  cout << "my_Var_Yj = " << my_Var_Yj << ", my_Var_Yjp = " << my_Var_Yjp << ", my_AMIC = " << my_AMIC << ", temp = " << temp << endl;
#endif

  double log_T       = log10(my_total_genome_length/150.0);
  double log_AMIC    = log10(my_AMIC);
  double sq_log_AMIC = log_AMIC*log_AMIC;

#if 0
  cout << "temp2 = " << log_T << ", temp3 = " << log_AMIC << ", temp4 = " << sq_log_AMIC << endl;
#endif

  const vector<double>& preset_cutpoints = my_parser->get_preset_cutpoints();

  if( SAGE::isnan(preset_cutpoints[Cu]) )
  {
    my_cutpoints[Cu] = 0.421 + 0.506*log_T + 1.162*log_AMIC + 0.472*sq_log_AMIC;
    my_cutpoints[Cu] = -pow(10.0, my_cutpoints[Cu]);
    my_adjusted_cutpoints[Cu] = my_cutpoints[Cu];

#if 0
  cout << "Cu cutpoint = " << my_cutpoints[Cu] << endl;
#endif
  }

  if( SAGE::isnan(preset_cutpoints[Ch]) )
  {
    my_cutpoints[Ch] = (-0.141) + 0.524*log_T + 0.237*log_AMIC - 0.861*sq_log_AMIC;
    my_cutpoints[Ch] = -pow(10.0, my_cutpoints[Ch]);
    my_adjusted_cutpoints[Ch] = my_cutpoints[Ch];

#if 0
  cout << "Ch cutpoint = " << my_cutpoints[Ch] << endl;
#endif
  }

  if( SAGE::isnan(preset_cutpoints[Cm]) )
  {
    my_cutpoints[Cm]=3.27;
    my_adjusted_cutpoints[Cm] = my_cutpoints[Cm];

#if 0
  cout << "Cm cutpoint = " << my_cutpoints[Cm] << endl;
#endif
  }

  if( SAGE::isnan(preset_cutpoints[Cp]) )
  {
//    my_cutpoints[Cp] = 0.475 + 0.518*log_T + 2.220*log_AMIC;
//    my_cutpoints[Cp] = 0.2 + 0.518*log_T + 2.220*log_AMIC;
    my_cutpoints[Cp] = 0.3 + 0.518*log_T + 2.220*log_AMIC;
    my_cutpoints[Cp] = -pow(10.0, my_cutpoints[Cp]);
    my_adjusted_cutpoints[Cp] = my_cutpoints[Cp];

#if 0
  cout << "Cp cutpoint = " << my_cutpoints[Cp] << endl;
#endif
  }
}

void reltest_analysis::reclassify_pairs()
{
  //classify pairs
  double temp_Yj;
  bool   sib;

  for( size_t p = 0; p < my_ptt_pairs[my_current_pairtype].size(); ++p )
  {
    my_ptt_pairs[my_current_pairtype][p].set_Yj(my_ptt_pairs[my_current_pairtype][p].get_Yj()/my_Var_Yj);
    my_ptt_pairs[my_current_pairtype][p].set_Yjp(my_ptt_pairs[my_current_pairtype][p].get_Yjp()/my_Var_Yjp);

    temp_Yj = my_ptt_pairs[my_current_pairtype][p].get_Yj();

    sib = false;

    if( temp_Yj <= my_adjusted_cutpoints[Cu] ) 
      my_ptt_pairs[my_current_pairtype][p].set_type(UNRELATED); 
    else if( temp_Yj > my_adjusted_cutpoints[Cu] && temp_Yj <= my_adjusted_cutpoints[Ch] )
      my_ptt_pairs[my_current_pairtype][p].set_type(HSIB); 
    else if( temp_Yj > my_adjusted_cutpoints[Ch] && temp_Yj <= my_adjusted_cutpoints[Cm] )
    {
      my_ptt_pairs[my_current_pairtype][p].set_type(SIB);
      sib = true;
    }
    else if( temp_Yj > my_adjusted_cutpoints[Cm] ) 
      my_ptt_pairs[my_current_pairtype][p].set_type(MZTWIN);

    if(    sib
        && my_ptt_pairs[my_current_pairtype][p].get_Yjp() <= my_adjusted_cutpoints[Cp] )
      my_ptt_pairs[my_current_pairtype][p].set_type(PARENT_OFFSPRING);
  }
}

void reltest_analysis::get_misclassified_pairs()
{
  putative_type type = my_current_pairtype;

  my_misc_pairs[my_current_pairtype].resize(0);

  //assume that parents pairs are of UNRELATED class type
  if( type == MARITAL )
    type = UNRELATED;

  for( size_t i = 0; i < my_ptt_pairs[my_current_pairtype].size(); ++i )
  {
    if( my_ptt_pairs[my_current_pairtype][i].get_pair_type() != type )
      my_misc_pairs[my_current_pairtype].push_back(my_ptt_pairs[my_current_pairtype][i]);		
  }
}

// Compute the new cutpoint here.
void
reltest_analysis::nonparametric_estimation()
{
  // Set the Yj & Yjp, but need to recover old one to use later
  // at the reclassify().
  for( size_t p = 0; p < my_ptt_pairs[my_current_pairtype].size(); ++p )
  {
    //cout << "processing pair " << p+1 << " out of "
    //     << my_ptt_pairs[my_current_pairtype].size();

    double Yj  = my_ptt_pairs[my_current_pairtype][p].get_Yj();
    double Yjp = my_ptt_pairs[my_current_pairtype][p].get_Yjp();

    //cout << ", Yj = " << Yj << ", Yjp = " << Yjp << endl;

    my_ptt_pairs[my_current_pairtype][p].set_Yj(Yj/my_Var_Yj);
    my_ptt_pairs[my_current_pairtype][p].set_Yjp(Yjp/my_Var_Yjp);
  }

  pair<double, size_t> re = estimate_shift(true, my_solution_Yj);

  my_mu_Yj     = re.first;
  my_picked_Yj = re.second;

  //cout << "my_mu_Yj = " << my_mu_Yj << ", my_picked_Yj = " << my_picked_Yj << endl;

  adjust_cutpoints(true);

  re = estimate_shift(false, my_solution_Yjp);

  my_mu_Yjp     = re.first;
  my_picked_Yjp = re.second;

  adjust_cutpoints(false);

  // Reset the Yj & Yjp.
  //
  for( size_t p = 0; p < my_ptt_pairs[my_current_pairtype].size(); ++p )
  {
    double Yj  = my_ptt_pairs[my_current_pairtype][p].get_Yj();
    double Yjp = my_ptt_pairs[my_current_pairtype][p].get_Yjp();

    my_ptt_pairs[my_current_pairtype][p].set_Yj(Yj*my_Var_Yj);
    my_ptt_pairs[my_current_pairtype][p].set_Yjp(Yjp*my_Var_Yjp);
  }
}

pair<double, size_t>
reltest_analysis::estimate_shift(bool is_Yj, vector<solution_type>& solution)
{
  solution.resize(0);
  solution.resize(3);

  // Method 1 - old maxfunapi
  //
/*
  vector< std::pair<double, double> > init_theta;
  encode_params_trial_vector(init_theta);

  L2_error_procedure shift_function(is_Yj, my_ptt_pairs[my_current_pairtype]);

  Maxfun maxfun(shift_function);

  initialize_maxfun(maxfun);
  
  for( size_t ti = 0; ti < init_theta.size(); ++ti )
  {
    maxfun.thin(0) = init_theta[ti].first;
    maxfun.thin(1) = init_theta[ti].second;

    run_optimum_maxfun(maxfun);

    solution[ti].function_value = maxfun.value();
    solution[ti].mean           = maxfun.param(0);
    solution[ti].variance       = maxfun.param(1);

#if 0
  cout << "function_value = " << solution[ti].function_value << endl;
  cout << "mean           = " << solution[ti].mean           << endl;
  cout << "variance       = " << solution[ti].variance       << endl;
  cout << "nfe : " << maxfun.evaluations()
       << ", lfl : " << maxfun.last_error() << endl << endl;
#endif
  }
*/

  // Method 2 - using new maxfunapi
  //
  do_new_L2_procedure(solution, is_Yj);

  const double qNaN = std::numeric_limits<double>::quiet_NaN();
  double mu     = qNaN;
  size_t picked = 3;

  double maximum_value = -numeric_limits<double>::infinity();
  double smallest_mean =  numeric_limits<double>::infinity();

  size_t max_pos = 3;

  for( size_t ti = 0; ti < 3; ++ti )
  {
    if(    solution[ti].function_value >= maximum_value
        && fabs(solution[ti].mean) <= smallest_mean )
    {
      maximum_value = solution[ti].function_value;
      smallest_mean = fabs(solution[ti].mean);
      max_pos = ti;
    }
  }

#if 0
  cout << "max_pos = " << max_pos << endl;
#endif

  bool valid_pos = check_mean_validity(max_pos, is_Yj);

  if( valid_pos )
  {
    picked = max_pos;
    mu = solution[max_pos].mean;

    return make_pair(mu, picked);
  }
  else
  {
    maximum_value = -numeric_limits<double>::infinity();
    smallest_mean =  numeric_limits<double>::infinity();

    size_t max_pos2 = 3;

    for( size_t ti = 0; ti < 3; ++ti )
    {
      if( ti == max_pos )
        continue;

      if(    solution[ti].function_value >= maximum_value
          && fabs(solution[ti].mean) <= smallest_mean )
      {
        maximum_value = solution[ti].function_value;
        smallest_mean = fabs(solution[ti].mean);
        max_pos2 = ti;
      }
    }
    
    valid_pos = check_mean_validity(max_pos2, is_Yj);

    if( valid_pos )
    {
      picked = max_pos2;
      mu = solution[max_pos2].mean;

      return make_pair(mu, picked);
    }
    else
    {
      maximum_value = -numeric_limits<double>::infinity();
      smallest_mean =  numeric_limits<double>::infinity();

      size_t max_pos3 = 3;

      for( size_t ti = 0; ti < 3; ++ti )
      {
        if( ti == max_pos || ti == max_pos2 )
          continue;

        if(    solution[ti].function_value >= maximum_value
            && fabs(solution[ti].mean) <= smallest_mean )
        {
          maximum_value = solution[ti].function_value;
          smallest_mean = fabs(solution[ti].mean);
          max_pos3 = ti;
        }
      }
      
      valid_pos = check_mean_validity(max_pos3, is_Yj);

      if( valid_pos )
      {
        picked = max_pos3;
        mu = solution[max_pos3].mean;

        return make_pair(mu, picked);
      }
    }
  }

#if 0
  cout << "??? ever ???" << endl;
#endif

  max_pos = 3;
  smallest_mean = numeric_limits<double>::infinity();

  for( size_t ti = 0; ti < 3; ++ti )
  {
    if( fabs(solution[ti].mean) <= smallest_mean )
    {
      smallest_mean = fabs(solution[ti].mean);
      max_pos = ti;
    }
  }

  picked = max_pos;
  mu = qNaN;

  return make_pair(mu, picked);
}

bool
reltest_analysis::check_mean_validity(size_t max_pos, bool is_Yj)
{
  if( max_pos >= 3 )
    return false;

  double Y_min =  numeric_limits<double>::infinity();
  double Y_max = -numeric_limits<double>::infinity();

  for( size_t p = 0; p < my_ptt_pairs[my_current_pairtype].size(); ++p )  
  {
    if( is_Yj )
    {
      Y_min = min(Y_min, my_ptt_pairs[my_current_pairtype][p].get_Yj());
      Y_max = max(Y_max, my_ptt_pairs[my_current_pairtype][p].get_Yj());
    }
    else
    {
      Y_min = min(Y_min, my_ptt_pairs[my_current_pairtype][p].get_Yjp());
      Y_max = max(Y_max, my_ptt_pairs[my_current_pairtype][p].get_Yjp());
    }
  }

  int max_bin = (int) (ceil( ((Y_max+0.1) - (Y_min-0.01)) / 0.25 ));

  Histogram hist_Y(max_bin, Y_min-0.01, (Y_min-0.01)+0.25*max_bin);

  for( size_t p = 0; p < my_ptt_pairs[my_current_pairtype].size(); ++p )  
  {
    if( is_Yj )
      hist_Y.add(my_ptt_pairs[my_current_pairtype][p].get_Yj());
    else
      hist_Y.add(my_ptt_pairs[my_current_pairtype][p].get_Yjp());      
  }

  size_t max_index = hist_Y.maxFrequency().first;

  double lower     = hist_Y.lower_bound();
  double bin_lower = hist_Y.bin_range(max_index).first  + lower;
  double bin_upper = hist_Y.bin_range(max_index).second + lower;
  double bin_mean  = (bin_lower + bin_upper) / 2.0;

  double mean = my_solution_Yj[max_pos].mean;
  if( !is_Yj )
    mean = my_solution_Yjp[max_pos].mean;

#if 0
  size_t max_count = hist_Y.maxFrequency().second;
  cout << "max_index = " << max_index << endl
       << "max_count = " << max_count << endl
       << "lower_bnd = " << lower     << endl
       << "total_pnt = " << hist_Y.total_point_count() << endl
       << "bin_range = (" << bin_lower << " to " << bin_upper << ")" << endl
       << "bin_ mean = " << bin_mean << endl
       << "my_mu     = " << mean << endl << endl;
#endif

  return ( mean >= (bin_mean - 1.5) && mean <= (bin_mean + 1.5) );
}

void
reltest_analysis::adjust_cutpoints(bool is_Yj)
{
  double mu;
  if( is_Yj )
  {
    mu = my_mu_Yj;
    if( SAGE::isnan(mu) )
      mu = my_solution_Yj[my_picked_Yj].mean;
  }
  else
  {
    mu = my_mu_Yjp;
    if( SAGE::isnan(mu) )
      mu = my_solution_Yjp[my_picked_Yjp].mean;
  }

  if( SAGE::isnan(mu) )
  {
    errors << priority(error)
           << "Mean is not a number. Adjusting the cutpoints failed."
           << "Old cutpoints will be used." << endl;
    return;
  }

  const vector<double>& preset_cutpoints = my_parser->get_preset_cutpoints();

  if( is_Yj )
  {
    if( SAGE::isnan(preset_cutpoints[Cu]) )
      my_adjusted_cutpoints[Cu] = my_cutpoints[Cu] + mu;

    if( SAGE::isnan(preset_cutpoints[Ch]) )
      my_adjusted_cutpoints[Ch] = my_cutpoints[Ch] + mu;

    if( SAGE::isnan(preset_cutpoints[Cm]) )
      my_adjusted_cutpoints[Cm] = my_cutpoints[Cm] + mu;
  }
  else
  {
    if( SAGE::isnan(preset_cutpoints[Cp]) )
      my_adjusted_cutpoints[Cp] = my_cutpoints[Cp] + mu;
  }
}

void
reltest_analysis::encode_params_trial_vector(vector< std::pair<double, double> >& init_theta)
{
  double init_mean[]     = { 0.0, -1.0, 1.0 };
  double init_variance[] = { 0.5,  0.5, 0.5 };
  
  init_theta.resize(3);

  for( size_t i = 0; i < init_theta.size(); ++i )
  {
    init_theta[i].first  = init_mean[i];
    init_theta[i].second = init_variance[i];
  }
}

void
reltest_analysis::initialize_maxfun(Maxfun& maxfun)
{
  const double inf    = std::numeric_limits<double>::infinity();
  const double ne_inf = -std::numeric_limits<double>::infinity();

  maxfun.nt() = 2;

  maxfun.thl(0)     = ne_inf;     //   lower bound
  maxfun.thu(0)     = inf;        //   upper bound
  maxfun.istin(0)   = 1;          //   independent, but used in depar

  maxfun.thl(1)     = 0.;         //   lower bound for variance
  maxfun.thu(1)     = inf;        //   upper bound
  maxfun.istin(1)   = 1;          //   independent, but used in depar
}

void
reltest_analysis::run_optimum_maxfun(Maxfun& maxfun)
{
  int    max_mthds[] = { 2,     5,     2,     6 };
  double max_epsc1[] = { 1e-3,  1e-3,  1e-4,  1e-12 };
  double max_epsc2[] = { 1e-15, 1e-15, 1e-15, 1e-15 };
  int    max_maxit[] = { 1,     20,    50,    20 };

  for( size_t m = 0; m < 4; ++m )
  {
    // Do not do maximization 2 or 3 if we've converged on the previous step.
    if( m > 0 && m < 3 && maxfun.last_error() < 4)
      continue;

    maxfun.method()  = max_mthds[m];
    maxfun.epsc1()   = max_epsc1[m];
    maxfun.epsc2()   = max_epsc2[m];
    maxfun.maxit()   = maxfun.nt() * max_maxit[m];

    if( m )
      for(int j = 0; j < maxfun.nt(); ++j)
        maxfun.thin(j) = maxfun.param(j);

    maxfun.run();
  }
}

void
reltest_analysis::do_new_L2_procedure(vector<solution_type>& solution, bool is_Yj)
{
  L2_error_procedure shift_function(is_Yj, my_ptt_pairs[my_current_pairtype]);

  const double inf    = std::numeric_limits<double>::infinity();
  const double ne_inf = -std::numeric_limits<double>::infinity();

  double init_mean[]     = { 0.0, -1.0, 1.0 };
  double init_variance[] = { 0.5,  0.5, 0.5 };

  for( size_t i = 0; i < 3; ++i )
  {
    MAXFUN::ParameterMgr pm;

    pm.addParameter("l2", "mean", MAXFUN::Parameter::INDEPENDENT, init_mean[i], ne_inf, inf);
    pm.addParameter("l2", "variance", MAXFUN::Parameter::INDEPENDENT, init_variance[i], 0.0, inf);

    MAXFUN::DebugCfg    db;
    //MAXFUN::DebugCfg    db(MAXFUN::DebugCfg::COMPLETE);
    MAXFUN::SequenceCfg sc(MAXFUN::SequenceCfg::USER_DEFINED, "A");

    sc.addRunCfg(MAXFUN::RunCfg::DIRECT_WITHOUT, 1);
    sc.getLatestRunCfg().epsilon1 = 1e-3;
    sc.getLatestRunCfg().epsilon2 = 1e-15;
    sc.getLatestRunCfg().var_cov = MAXFUN::RunCfg::FINAL;

    sc.addRunCfg(MAXFUN::RunCfg::VAR_METRIC_IDENTITY, 20);
    sc.getLatestRunCfg().epsilon1 = 1e-3;
    sc.getLatestRunCfg().epsilon2 = 1e-15;
    sc.getLatestRunCfg().var_cov = MAXFUN::RunCfg::FINAL;

    sc.addRunCfg(MAXFUN::RunCfg::DIRECT_WITHOUT, 50);
    sc.getLatestRunCfg().epsilon1 = 1e-4;
    sc.getLatestRunCfg().epsilon2 = 1e-15;
    sc.getLatestRunCfg().var_cov = MAXFUN::RunCfg::FINAL;

    sc.addRunCfg(MAXFUN::RunCfg::VAR_METRIC_ESTIMATE, 20);
    sc.getLatestRunCfg().epsilon1 = 1e-12;
    sc.getLatestRunCfg().epsilon2 = 1e-15;
    sc.getLatestRunCfg().var_cov = MAXFUN::RunCfg::FINAL;

    MAXFUN::Results maxfun_info = MAXFUN::Maximizer::Maximize(sc, pm, shift_function, db);

    solution[i].function_value = maxfun_info.getFinalFunctionValue();
    solution[i].mean           = maxfun_info.getParameterMgr().getParameter(0).getFinalEstimate();
    solution[i].variance       = maxfun_info.getParameterMgr().getParameter(1).getFinalEstimate();

#if 0
  std::cout << MAXFUN::OutputFormatter::convertEstimates (maxfun_info)
            << MAXFUN::OutputFormatter::convertMatrix    (maxfun_info)
            << MAXFUN::JointTest(maxfun_info, maxfun_info).summarizeAsTable();

  cout << "function_value = " << solution[i].function_value << endl;
  cout << "mean           = " << solution[i].mean           << endl;
  cout << "variance       = " << solution[i].variance       << endl;
  cout << "nfe : " << maxfun_info.getIterations()
       << ", lfl : " << maxfun_info.getExitFlag() << endl << endl;
#endif
  }

  return;
}

void
reltest_analysis::test_deviation()
{
  double sumYj   = 0.;
  double sumYjYj = 0.;
  double count   = 0.;

  for( size_t p = 0; p < my_ptt_pairs[my_current_pairtype].size(); ++p )
    if( my_ptt_pairs[my_current_pairtype][p].get_pair_type() == SIB )
    {
      double Yj = my_ptt_pairs[my_current_pairtype][p].get_Yj();

      sumYj   += Yj;
      sumYjYj += (Yj*Yj);
      count++;
    }

  my_Yj_mean        = sumYj / count;
  my_standard_error = sqrt( (sumYjYj - (sumYj * sumYj / count)) / (count * count) );

#if 0
  cout << "count      = " << count << endl;
  cout << "sumYj      = " << sumYj << endl;
  cout << "sumYjYj    = " << sumYjYj << endl;
  cout << "my_std_err = " << my_standard_error << endl;
  cout << "my_Yj_mean = " << my_Yj_mean << " +- " <<  2.0 * my_standard_error <<endl;
#endif
}

} // end of namespace RELTEST

} // end of naemspace SAGE
