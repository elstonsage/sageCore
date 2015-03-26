//====================================================================
//  Pairwise IBD Analysis - Single point Relative Pair IBD generator. 
//                                                                    
//  Author: Geoff Wedig (wedig@darwin.cwru.edu)                       
//                                                                    
//  History:   0.1  gcw Relpal's Generator          until July 1998   
//             0.2  gcw Adapted for GENIBD          Aug/Sept   1998   
//             1.0  Updated to new libraries.       yjs   May  2004   
//                                                                    
//  Copyright (c) 1998  R.C. Elston                                   
//====================================================================

#include "genibd/pair_ibd_analysis.h"

#undef PAIR_GEN_DEBUG

namespace SAGE
{

namespace GENIBD
{

pair_ibd_analysis::pair_ibd_analysis(cerrorstream& e)
                 : errors(e)
{
  my_peeler = NULL;
  my_ibds   = NULL;
  my_dots   = NULL;

  my_built = false;
  my_valid = false;
  my_verbose = false;

  my_current_marker = 0;
}

pair_ibd_analysis::~pair_ibd_analysis()
{
  if( my_ibds )
    delete my_ibds;

  if( my_dots )
    delete my_dots;
}

bool
pair_ibd_analysis::build()
{
  return my_built = true;  
}

bool
pair_ibd_analysis::set_pedigree(const meiosis_map& mmap, const pedigree_region& pr)
{
  my_meiosis_map = mmap;

  my_ped_region = pr;

  my_region = my_ped_region.get_region();

  my_pairs.resize(0);

  return true;
}

bool
pair_ibd_analysis::build_ibds()
{
  if( my_ibds )
    delete my_ibds;

  my_ibds = new basic_storage_ibd(my_meiosis_map, my_region, false);

  return true;
}

bool
pair_ibd_analysis::add_pair(fmember_const_pointer m1, fmember_const_pointer m2,
                            fmember_const_pointer c1, fmember_const_pointer c2,
                            pair_generator::pair_type t)
{
  if( !my_ibds )
    return false;

  filtered_relative_pair pair(m1, m2, c1, c2, t);

  my_ibds->add_pair(m1, m2);
  my_pairs.push_back(pair);

  return true;
}

bool
pair_ibd_analysis::compute(const string& title)
{
  if( my_dots )
    delete my_dots;

  my_dots = new text_dot_formatter(cout);

  my_dots->set_prefix_width(43);
  my_dots->set_trigger_count(my_region.locus_count());

  if( !built() || !my_meiosis_map.built() || !my_region.valid() )
    return my_valid = false;

#ifdef PAIR_GEN_DEBUG
  cout << "Beginning Generation of Pi-Hats" << endl;
#endif

  if( !my_ibds->built() )
  {
    my_ibds->build();

    assert(my_ibds->built());
  }

  my_dots->set_prefix("        Generating Single Point Likelihoods");

#ifdef PAIR_GEN_DEBUG
    cout << "Doing pedigree " << my_meiosis_map.get_pedigree()->name() << endl;
#endif

  for( my_current_marker = 0; my_current_marker < my_region.locus_count(); ++my_current_marker )
  {
    my_dots->trigger();

    const inheritance_model& model = my_ped_region[my_current_marker];

    if( !model.genotype_informative() && !model.penetrance_informative() )
      continue;

    if( my_peeler )
      delete my_peeler;

    my_peeler = new peeler(*(my_meiosis_map.get_subpedigree()), model);

    compute_subpedigree_likelihood(model); // calculate the ped likelihood

    if( my_subped_likelihood == 0.0 )
    {
      cout << endl << endl;

      errors << SAGE::priority(SAGE::error) << "Pedigree '"
             << my_meiosis_map.get_pedigree()->name() << "' has zero likelihood at marker '"
             << model.name()
             << "'.  Marker will be skipped for this pedigree."
             << endl;

      cout << endl;

      my_dots->spammed();

      continue;
    }


#ifdef PAIR_GEN_DEBUG
    cout << "Starting work on Pi Hats" << endl;
#endif

    generate_pi_hats(model);

#ifdef PAIR_GEN_DEBUG
    cout << "Finished marker '" << model.name() << "'" << endl;
#endif
  }

  return my_valid = true;
}

// Calculates the likelihood of a filtered_subpedigree at the current marker.
//
void
pair_ibd_analysis::compute_subpedigree_likelihood(const inheritance_model& model)
{
#if 0
  cout << "pair_ibd_analysis::compute_subpedigree_likelihood()..." << endl;

  fmember_const_iterator  mem2 = my_meiosis_map.get_subpedigree()->member_begin();
  for( ; mem2 != my_meiosis_map.get_subpedigree()->member_end(); ++mem2 )
  {
    cout << "ind = " << mem2->name() << endl;

    size_t m_index2 = mem2->subindex() + 1;

    phased_pen_iter m_iter2 = model.phased_penetrance_begin(m_index2);

    for( ; m_iter2 != model.phased_penetrance_end(m_index2); ++m_iter2 )
    {
      phased_pen_iter mi_iter2(m_iter2);

      cout << "anterior = " << my_peeler->anterior(*mem2, mi_iter2)
           << ", posterior = " << my_peeler->posterior(*mem2, mi_iter2)
           << ", penetrance = " << (*mi_iter2)
           << endl;
    }
  }
  cout << endl;
#endif

  fmember_const_iterator  mem = my_meiosis_map.get_subpedigree()->member_begin();  // Find the first founder of the pedigree
  size_t                  m_index = mem->subindex() + 1;

  log_double  like(0);

  phased_pen_iter m_iter = model.phased_penetrance_begin(m_index);

  for( ; m_iter != model.phased_penetrance_end(m_index); ++m_iter )
  {
    log_double like_term(1);

    phased_pen_iter mi_iter(m_iter);

    like_term *= my_peeler->anterior(*mem, mi_iter); 
    like_term *= (*mi_iter);
    like_term *= my_peeler->posterior(*mem, mi_iter);

    like += like_term;
  }

  my_subped_likelihood = like.get_double();

#if 0
  cout << "my_subped_likelihood = " << my_subped_likelihood << endl;
#endif
}

// generate_pi_hats iterates through the pairs for a pedigree and calculates
//   the pi^ values for the current marker.
//
void
pair_ibd_analysis::generate_pi_hats(const inheritance_model& model)
{

  for( size_t current_pair = 0; current_pair < my_pairs.size(); ++current_pair )
  {

#ifdef PAIR_GEN_DEBUG
    cout << "Doing:  Marker: " << model.name()
         << "  Pedigree: " << my_meiosis_map.get_pedigree()->name()
         << " Pair: " 
         << my_pairs[current_pair].member_one->name() << " and "
         << my_pairs[current_pair].member_two->name() << "..." << flush;
#endif

    switch( my_pairs[current_pair].type )
    {
      case pair_generator::SIBSIB  : generate_sib(current_pair, model);
                                     break;

      case pair_generator::HALFSIB : generate_hsib(current_pair, model);
                                     break;

      case pair_generator::AVUNC   : generate_avunc(current_pair, model);
                                     break;

      case pair_generator::GRANDP  : generate_gpar(current_pair, model);
                                     break;

      case pair_generator::COUSIN  : generate_cous(current_pair, model);

      default                      : break;
    }
#ifdef PAIR_GEN_DEBUG
    cout << "Done." << endl << flush;
#endif

  }
}

// Calculates the sibling pi^
//
void
pair_ibd_analysis::generate_sib(size_t current_pair, const inheritance_model& model)
{
  //Variables:
  //   ant_mother,  ant_father - anteriors of the parents.
  //  post_mother, post_father - posteriors of parents without current family.
  //                  ind_v[8] - ind_ac(1), ind_ad(1), ind_bc(1), .. ind_bd(2)
  //           v1,          v2 - The v value for the current set of genotypes for 1,2 ibd.
  //          lp1,         lp2 - Summation of second parent for 1,2 allele ibd.
  //          lt1,         lt2 - total summation for 1 and 2 allele ibd.
  //                      sibs - Posterior probablility of the additional sibs given parents.

  // New added variables for maternal & paternal bit split (m:mother,f:father) - yjs Jun. 2002
  //   vm,       vf
  //  lpm,      lpf
  //  ltm,      ltf

  filtered_relative_pair& sib_pair = my_pairs[current_pair];

//  if(    sib_pair.member_one->info().phenotype_missing(my_current_marker, model)
//      || sib_pair.member_two->info().phenotype_missing(my_current_marker, model) )
//    return;

  double lp1 = 0., lp2 = 0., lpm = 0., lpf = 0.;
  double lt1 = 0., lt2 = 0., ltm = 0., ltf = 0.;

  fmember_const_pointer sib1   = sib_pair.member_one;
  fmember_const_pointer sib2   = sib_pair.member_two;
  fmember_const_pointer mother = my_meiosis_map.mother(sib1);
  fmember_const_pointer father = my_meiosis_map.father(sib1);

  // Phenotype Indices
  size_t sib1_index   = sib1->subindex() + 1;
  size_t sib2_index   = sib2->subindex() + 1;
  size_t mother_index = mother->subindex() + 1;
  size_t father_index = father->subindex() + 1;

  phased_pen_iter moth_pen_iter = model.phased_penetrance_begin(mother_index);
  for( ; moth_pen_iter != model.phased_penetrance_end(mother_index); ++moth_pen_iter )
  {
    lp1 = lp2 = lpm = lpf = 0.;

    if( *moth_pen_iter == 0. )//|| moth_pen_iter.geno_id() < 0 )
       continue;

    phased_pen_iter fath_pen_iter = model.phased_penetrance_begin(father_index);
    for( ; fath_pen_iter != model.phased_penetrance_end(father_index); ++fath_pen_iter )
    {
      // If x-linked, skip over non-Y genotypes for the father.
      if( model.is_x_linked() && !is_Y_genotype(fath_pen_iter.phased_geno()) )
        continue;

      if( *fath_pen_iter == 0. )//|| fath_pen_iter.geno_id() < 0 )
         continue;

      vector<double> ind_v(8, 0.);

      // Get possible genotypes of sibs
 
      // Get the parental alleles
      ind_genotype mi_geno(moth_pen_iter);
      ind_genotype fi_geno(fath_pen_iter);
      
      MLOCUS::allele m_allele1 = mi_geno.mg.allele1();
      MLOCUS::allele m_allele2 = mi_geno.mg.allele2();
      MLOCUS::allele f_allele1 = fi_geno.mg.allele1();
      MLOCUS::allele f_allele2 = fi_geno.mg.allele2();

      // First Sibling
      phased_pen_iter s1_iter = model.phased_penetrance_begin(sib1_index);
      for( ; s1_iter != model.phased_penetrance_end(sib1_index); ++s1_iter )
      {
        if(    model.is_x_linked()
            && sib1->is_male()
            && !is_Y_genotype(s1_iter.phased_geno()) )
          continue;

        // Get sib's alleles
        MLOCUS::allele s_allele1 = s1_iter.phased_geno().allele1();
        MLOCUS::allele s_allele2 = s1_iter.phased_geno().allele2();
        
        // If first allele could come from mom's first and second from dad's first
        if( s_allele1 == m_allele1 && s_allele2 == f_allele1 )
        {
          ind_v[0] += (*s1_iter) * my_peeler->posterior(*sib1, s1_iter).get_double();
        }

        // If first allele could come from mom's first and second from dad's second
        if( s_allele1 == m_allele1 && s_allele2 == f_allele2 )
        {
          ind_v[1] += (*s1_iter) * my_peeler->posterior(*sib1, s1_iter).get_double();
        }

        if( s_allele1 == m_allele2  && s_allele2 == f_allele1 )
        {
          ind_v[2] += (*s1_iter) * my_peeler->posterior(*sib1, s1_iter).get_double();
        }

        if( s_allele1 == m_allele2 && s_allele2 == f_allele2 )
        {
          ind_v[3] += (*s1_iter) * my_peeler->posterior(*sib1, s1_iter).get_double();
        }
      }

      phased_pen_iter s2_iter = model.phased_penetrance_begin(sib2_index);
      for( ; s2_iter != model.phased_penetrance_end(sib2_index); ++s2_iter )
      {
        if(    model.is_x_linked()
            && sib2->is_male()
            && !is_Y_genotype(s2_iter.phased_geno()) )
          continue;

        // Get sib's alleles
        MLOCUS::allele s_allele1 = s2_iter.phased_geno().allele1();
        MLOCUS::allele s_allele2 = s2_iter.phased_geno().allele2();

        if( s_allele1 == m_allele1 && s_allele2 == f_allele1 )
        {
          ind_v[4] += (*s2_iter) * my_peeler->posterior(*sib2, s2_iter).get_double();
        }

        if( s_allele1 == m_allele1 && s_allele2 == f_allele2 )
        {
          ind_v[5] += (*s2_iter) * my_peeler->posterior(*sib2, s2_iter).get_double();
        }

        if( s_allele1 == m_allele2 && s_allele2 == f_allele1 )
        {
          ind_v[6] += (*s2_iter) * my_peeler->posterior(*sib2, s2_iter).get_double();
        }

        if( s_allele1 == m_allele2 && s_allele2 == f_allele2 )
        {
          ind_v[7] += (*s2_iter) * my_peeler->posterior(*sib2, s2_iter).get_double();
        }
      }

      // Calculate v1 and v2
      double both_sharing   =   ind_v[0] * ind_v[4] + ind_v[1] * ind_v[5]
                              + ind_v[2] * ind_v[6] + ind_v[3] * ind_v[7];
      double mother_sharing =   ind_v[0] * ind_v[5] + ind_v[1] * ind_v[4]
                              + ind_v[2] * ind_v[7] + ind_v[3] * ind_v[6];
      double father_sharing =   ind_v[0] * ind_v[6] + ind_v[1] * ind_v[7]
                              + ind_v[2] * ind_v[4] + ind_v[3] * ind_v[5];
#if 0
  cout << endl
       << "current_pair = " << current_pair
       << ", S1 = " << sib1->name() 
       << ", S2 = " << sib2->name() 
       << ", marker = " << my_current_marker
       << ", ma_geno = " << mi_geno.mg.name()
       << ", fa_geno = " << fi_geno.mg.name() << endl;

  for( size_t a = 0; a < ind_v.size(); ++a )
    cout << "	" << ind_v[a];
  cout << endl;

  cout << "  mother_sharing = " << mother_sharing << endl;
  cout << "  father_sharing = " << father_sharing << endl;
  cout << "  both_sharing   = " << both_sharing << endl;
#endif

      double v1 = 0.125 * ( mother_sharing + father_sharing );
      double v2 = 0.25  * ( both_sharing );

      double vm = 0.25 * ( mother_sharing);
      double vf = 0.25 * ( father_sharing);

      if( v1 || v2 )
      {
        double sibs        = my_peeler->likelihood_of_sibs_except_sib(mi_geno, fi_geno, *sib1, *sib2).get_double();
        double ant_father  = my_peeler->anterior(*father, fath_pen_iter).get_double();
        double post_father = my_peeler->posterior_except_mate(*father, *mother, fath_pen_iter).get_double();

        lp1 += (*fath_pen_iter) * ant_father * post_father * v1 * sibs;
        lp2 += (*fath_pen_iter) * ant_father * post_father * v2 * sibs;

#if 0
  cout << "  sibs         = " << sibs << endl
       << "  ant_father   = " << ant_father << endl
       << "  post_father  = " << post_father << endl
       << "  f_penetrance = " << (*fath_pen_iter) << endl;
#endif
      }

      if( vm || vf )
      {
        double sibs        = my_peeler->likelihood_of_sibs_except_sib(mi_geno, fi_geno, *sib1, *sib2).get_double();
        double ant_father  = my_peeler->anterior(*father, fath_pen_iter).get_double();
        double post_father = my_peeler->posterior_except_mate(*father, *mother, fath_pen_iter).get_double();

        lpm += (*fath_pen_iter) * ant_father * post_father * vm * sibs;
        lpf += (*fath_pen_iter) * ant_father * post_father * vf * sibs;
      }
    }  // End Second Parent (Father) Calculations

    double ant_mother  = my_peeler->anterior(*mother, moth_pen_iter).get_double();
    double post_mother = my_peeler->posterior_except_mate(*mother, *father, moth_pen_iter).get_double();

#if 0
  cout << "  ant_mother   = " << ant_mother << endl
       << "  post_mother  = " << post_mother << endl
       << "  m_penetrance = " << (*moth_pen_iter) << endl;
#endif

    if( lp1 ) lt1 += (*moth_pen_iter) * ant_mother * post_mother * lp1;
    if( lp2 ) lt2 += (*moth_pen_iter) * ant_mother * post_mother * lp2;

    if( lpm ) ltm += (*moth_pen_iter) * ant_mother * post_mother * lpm;
    if( lpf ) ltf += (*moth_pen_iter) * ant_mother * post_mother * lpf;
  } // End likelihood Calculations

  double pr1 = lt1 * 0.5  / my_subped_likelihood;
  double pr2 = lt2 * 0.25 / my_subped_likelihood;

  double prm = ltm * 0.25 / my_subped_likelihood;
  double prf = ltf * 0.25 / my_subped_likelihood;
  double prm_prf = numeric_limits<double>::quiet_NaN();

  if( finite(prm) && finite(prf) )
  {
    prm_prf = prm - prf;

    if( prm_prf < 0. && prm_prf > -std::numeric_limits<double>::epsilon() )
      prm_prf = 0.;
  }

#if 0
  cout << " pr1 = " << pr1 << ", pr2 = " << pr2 << ", prm = " << prm << ", prf = " << prf << endl;
#endif

  double pr0 = 1.0 - pr1 - pr2;

  if( 1.0 - pr1 - pr2 < 0.0 )
    pr0 = 0.0;
 
  if( model.is_x_linked() && is_brother_brother(sib1, sib2) )
  {
    pr0     = pr1;
    prm_prf = pr2;
    pr2     = 0.0;
  }
  
//  cout << "XX" << pr0 << ' ' << prm_prf << ' ' << pr2 << endl;

  my_ibds->set_ibd(current_pair, my_current_marker, pr0, prm_prf, pr2);
}

// half sib pi^s

void
pair_ibd_analysis::generate_hsib(size_t current_pair, const inheritance_model& model)
{
  //Variables:
  //    mutual_parent - Shared Parent's id.
  //   hsib1,   hsib2 - Half sib ids.
  // spouse1, spouse2 - Parent individuals (second parent for half sib n)
  //       ant_parent - Shared parent's anterior.
  //      post_parent - Posterior of the shared parent.
  //        post_mate - Posteriors with only the current mates.
  //            al[2] - Possible shared alleles.
  //                a - current allele for conditioning.
  //    fam_side1 & 2 - Probability of family of that side, given conditioning.
  //            total - Total likelihood.

  filtered_relative_pair& hsib_pair = my_pairs[current_pair];

//  if(    hsib_pair.member_one->info().phenotype_missing(my_current_marker, model)
//      || hsib_pair.member_two->info().phenotype_missing(my_current_marker, model) )
//    return;

  allele al[2];

  double total = 0.0;

  fmember_const_pointer hsib1         = hsib_pair.member_one;
  fmember_const_pointer hsib2         = hsib_pair.member_two;
  fmember_const_pointer mutual_parent = hsib_pair.connector_one;

  fmember_const_pointer spouse1 = NULL;
  fmember_const_pointer spouse2 = NULL;

  if( my_meiosis_map.mother(hsib1)->name() == my_meiosis_map.mother(hsib2)->name() )
  {
    spouse1 = my_meiosis_map.father(hsib1);
    spouse2 = my_meiosis_map.father(hsib2);
  }
  else if( my_meiosis_map.father(hsib1)->name() == my_meiosis_map.father(hsib2)->name() )
  {
    spouse1 = my_meiosis_map.mother(hsib1);
    spouse2 = my_meiosis_map.mother(hsib2);
  }

#if 0
  cout << endl << "current_pair = " << current_pair << ": "
       << hsib1->name() << "," << hsib2->name()
       << " -- mutual parent = " << mutual_parent->name()
       << ", sp1 = " << spouse1->name()
       << ", sp2 = " << spouse2->name() << endl;
#endif

  size_t hsib1_index         = hsib1->subindex() + 1;
  size_t hsib2_index         = hsib2->subindex() + 1;
  size_t mutual_parent_index = mutual_parent->subindex() + 1;
  //size_t spouse1_index       = spouse1->subindex() + 1;
  //size_t spouse2_index       = spouse2->subindex() + 1;

  phased_pen_iter h1_iter = model.phased_penetrance_begin(hsib1_index);
  for( ; h1_iter != model.phased_penetrance_end(hsib1_index); ++h1_iter )
  {
    if(    model.is_x_linked()
        && hsib1->is_male()
        && !is_Y_genotype(h1_iter.phased_geno()) )
      continue;

    phased_pen_iter h1i_iter(h1_iter);

    phased_pen_iter h2_iter = model.phased_penetrance_begin(hsib2_index);
    for( ; h2_iter != model.phased_penetrance_end(hsib2_index); ++h2_iter )
    {
      if(    model.is_x_linked()
          && hsib2->is_male()
          && !is_Y_genotype(h2_iter.phased_geno()) )
        continue;

      if( !condition(h1_iter.phased_geno(), h2_iter.phased_geno(), al) )
        continue;

      phased_pen_iter h2i_iter(h2_iter);

      // genotypes for hsibs ok (they share alleles ibs)

      for( size_t a = 0; a < 2 && al[a] != allele(); ++a )
      {
        phased_pen_iter mp_iter = model.phased_penetrance_begin(mutual_parent_index);
        for( ; mp_iter != model.phased_penetrance_end(mutual_parent_index); ++mp_iter )
        {
          if(    model.is_x_linked()
              && mutual_parent->is_male()
              && !is_Y_genotype(mp_iter.phased_geno()) )
            continue;

          if( !has_allele(mp_iter.phased_geno(), al[a]) )
            continue;

          phased_pen_iter mpi_iter(mp_iter);

          // Shared parent shares the allele, thus could be the source.

          double post_mate =   my_peeler->posterior_with_mate(*mutual_parent, *spouse1, mpi_iter).get_double()
                             * my_peeler->posterior_with_mate(*mutual_parent, *spouse2, mpi_iter).get_double();

          double fam_side1 = cond_post(mutual_parent, spouse1, hsib1, mpi_iter, h1i_iter, al[a], model);
          double fam_side2 = cond_post(mutual_parent, spouse2, hsib2, mpi_iter, h2i_iter, al[a], model);
#if 0
  cout << "  post_mate = " << post_mate << endl
       << "  fam_side1 = " << fam_side1 << endl
       << "  fam_side2 = " << fam_side2 << endl;
#endif
          if( post_mate && fam_side1 && fam_side2 )
          {
            double ant_parent  = my_peeler->anterior(*mutual_parent, mpi_iter).get_double();
            double post_parent = my_peeler->posterior(*mutual_parent, mpi_iter).get_double() / post_mate;

            total +=   ant_parent * post_parent
                     * (*mpi_iter)
                     * (*h1i_iter)
                     * (*h2i_iter)
                     * fam_side1 * fam_side2 
                     * has_allele(mp_iter.phased_geno(), al[a]);

#if 0
  cout << "  ant_parent  = " << ant_parent << endl
       << "  post_parent = " << post_parent << endl
       << "  total       = " << total << endl;
#endif

          }
        } // end shared parent 
      } //
    } // end half sib 2
  } // end half sib 1

  double f1mp = 0.0;
  if( mutual_parent->is_female() )
    f1mp = total / my_subped_likelihood;
  else
    f1mp = -(total / my_subped_likelihood);

#if 0
  cout << " total       = " << total << endl
       << " subped_like = " << my_subped_likelihood << endl
       << " f1mp        = " << f1mp << endl;
#endif

  double pr0 = 1.0 - total / my_subped_likelihood;

  if( model.is_x_linked() && is_brother_brother(hsib1, hsib2) )
  {
    pr0  = total / my_subped_likelihood;
    f1mp = 0.0;
  }

  my_ibds->set_ibd(current_pair, my_current_marker, pr0, f1mp, 0.0);
}

// avuncular pi^s

void pair_ibd_analysis::generate_avunc(size_t current_pair, const inheritance_model& model)
{
  //Variables:
  //  ant_grandp1,  ant_grandp2 - Grandparental anteriors.
  //                 post_uncle - avunc. posterior.
  //                post_parent - Parent's posterior without child's family.
  // post_grandp1, post_grandp2 - Grandparent's posteriors without current family.
  //                      al[2] - Possible shared alleles.
  //                          a - current allele_id.
  //                   cond_fam - Family of parent and child conditioned for allele.
  //                       sibs - Posterior Probability of sibs of parent and avunc.
  //                      total - Likelihood after conditioning.
  //                      sub_t - Partial multiplication of some terms to avoid extra work
  //               ha1,     ha2 - Granparent n has the shared allele.

  filtered_relative_pair& avunc_pair = my_pairs[current_pair];

//  if(    avunc_pair.member_one->info().phenotype_missing(my_current_marker, model)
//      || avunc_pair.member_two->info().phenotype_missing(my_current_marker, model) )
//    return;

  allele al[2];

  double total = 0.0;

  fmember_const_pointer uncle   = avunc_pair.member_one;
  fmember_const_pointer nephew  = avunc_pair.member_two;
  fmember_const_pointer parent1 = avunc_pair.connector_one;
  fmember_const_pointer grandp1 = my_meiosis_map.mother(parent1);
  fmember_const_pointer grandp2 = my_meiosis_map.father(parent1);

  fmember_const_pointer parent2;
  if( parent1 == my_meiosis_map.father(nephew) )
    parent2 = my_meiosis_map.mother(nephew);
  else
    parent2 = my_meiosis_map.father(nephew);

#if 0
  cout << endl << "current_pair = " << current_pair << ": "
       << uncle->name() << "," << nephew->name()
       << " -- parent = " << parent1->name()
       << ", gp1 = " << grandp1->name()
       << ", gp2 = " << grandp2->name() << endl;
#endif

  size_t uncle_index   = uncle->subindex() + 1;
  size_t nephew_index  = nephew->subindex() + 1;
  size_t parent1_index = parent1->subindex() + 1;
  //size_t parent2_index = parent2->subindex() + 1;
  size_t grandp1_index = grandp1->subindex() + 1;
  size_t grandp2_index = grandp2->subindex() + 1;

  phased_pen_iter n_iter = model.phased_penetrance_begin(nephew_index);
  for( ; n_iter != model.phased_penetrance_end(nephew_index); ++n_iter )
  {
    if(    model.is_x_linked()
        && nephew->is_male()
        && !is_Y_genotype(n_iter.phased_geno()) )
      continue;

    phased_pen_iter ni_iter(n_iter);

    phased_pen_iter u_iter = model.phased_penetrance_begin(uncle_index);
    for( ; u_iter != model.phased_penetrance_end(uncle_index); ++u_iter )
    {
      if(    model.is_x_linked()
          && uncle->is_male()
          && !is_Y_genotype(u_iter.phased_geno()) )
        continue;

      if( !condition(n_iter.phased_geno(), u_iter.phased_geno(), al) )
        continue;

#if 0
  cout << "allele shared : ";
  if( al[0] != allele() )
    cout << al[0].name();
  else
    cout << ".";
  cout << "/";
  if( al[1] != allele() )
    cout << al[1].name();
  else
    cout << ".";
  cout << endl;
#endif

      phased_pen_iter ui_iter(u_iter);

      // Child and avunc. share an allele

      double post_uncle = my_peeler->posterior(*uncle, ui_iter).get_double();

      for( size_t a = 0; a < 2 && al[a] != allele(); ++a )
      {
        phased_pen_iter p1_iter = model.phased_penetrance_begin(parent1_index);
        for( ; p1_iter != model.phased_penetrance_end(parent1_index); ++p1_iter )
        {
          if(    model.is_x_linked()
              && parent1->is_male()
              && !is_Y_genotype(p1_iter.phased_geno()) )
            continue;

          phased_pen_iter p1i_iter(p1_iter);

          double cond_fam = cond_post(parent1, parent2, nephew, p1i_iter, ni_iter, al[a], model);

          if( !has_allele(p1_iter.phased_geno(), al[a]) || !cond_fam )
            continue;

          // Avunc, child and parent all have the allele (may still be ibs)
          // and the conditional posterior of the parent is not 0.

          double post_parent = my_peeler->posterior_except_mate(*parent1, *parent2, p1i_iter).get_double();

#if 0
  cout << "  cond_fam    = " << cond_fam    << endl
       << "  parent_post = " << post_parent << endl
       << "  uncle_post  = " << post_uncle  << endl;

  size_t gp_i = 0;
#endif

          phased_pen_iter gp1_iter = model.phased_penetrance_begin(grandp1_index);
          for( ; gp1_iter != model.phased_penetrance_end(grandp1_index); ++gp1_iter )
          {
            if(    model.is_x_linked()
                && grandp1->is_male()
                && !is_Y_genotype(gp1_iter.phased_geno()) )
              continue;

            phased_pen_iter gp1i_iter(gp1_iter);

            double ant_grandp1  = my_peeler->anterior(*grandp1, gp1i_iter).get_double();
            double post_grandp1 = my_peeler->posterior_except_mate(*grandp1, *grandp2, gp1i_iter).get_double();

#if 0
  cout << endl << "  " << gp_i << endl;
  cout << "    gpant1  = " << ant_grandp1  << endl
       << "    gppost1 = " << post_grandp1 << endl;
#endif

            phased_pen_iter gp2_iter = model.phased_penetrance_begin(grandp2_index);
            for( ; gp2_iter != model.phased_penetrance_end(grandp2_index); ++gp2_iter )
            {
              if(    model.is_x_linked()
                  && grandp2->is_male()
                  && !is_Y_genotype(gp2_iter.phased_geno()) )
                continue;

              size_t ha1 = has_allele(gp1_iter.phased_geno(), al[a]);
              size_t ha2 = has_allele(gp2_iter.phased_geno(), al[a]);

              if( !(ha1 || ha2) ) continue;

#if 0
  cout << "    gp allele shared : " << al[a].name()
       << " - " << gp1_iter.phased_geno().allele1().name() << "/" << gp1_iter.phased_geno().allele2().name()
       << " - " << gp2_iter.phased_geno().allele1().name() << "/" << gp2_iter.phased_geno().allele2().name()
       << endl;
#endif
              phased_pen_iter gp2i_iter(gp2_iter);

              // At least one Grandparent has the shared allele

              double ant_grandp2  = my_peeler->anterior(*grandp2, gp2i_iter).get_double();
              double post_grandp2 = my_peeler->posterior_except_mate(*grandp2, *grandp1, gp2i_iter).get_double();

              double sub_t =   ant_grandp1 * post_grandp1 * ant_grandp2 * post_grandp2
                             * (*gp2i_iter)
                             * (*gp1i_iter)
                             * (*p1i_iter)
                             * (*ui_iter)
                             * (*ni_iter)
                             * cond_fam * post_parent * post_uncle;

#if 0
  cout << "      gpant2  = " << ant_grandp2  << endl
       << "      gppost2 = " << post_grandp2 << endl
       << "      sub_t   = " << sub_t        << endl;
  gp_i++;
#endif
              conditional_genotype con_geno1;

              bool ch_g1   = ch_geno1(gp1_iter.phased_geno(),
                                      gp2_iter.phased_geno(),
                                      grandp1->is_female(),
                                      p1_iter.phased_geno(),
                                      al[a], con_geno1);
              double cont1 = cond_tran(u_iter.phased_geno(), con_geno1);
              double sibs1 = my_peeler->likelihood_of_sibs_except_sib(gp1_iter.phased_geno(), gp2_iter.phased_geno(), *parent1, *uncle).get_double();

              if( ha1 && ch_g1 && cont1 && sibs1 )
              {
                total +=   sub_t * sibs1 * ha1
                         * cond_tran(p1_iter.phased_geno(), con_geno1)
                         * cond_tran(u_iter.phased_geno(), con_geno1);
#if 0
  cout << "      ch_g1 = " << ch_g1 << endl
       << "      cont1 = " << cont1 << endl
       << "      sibs1 = " << sibs1 << endl;
#endif
              }

              conditional_genotype con_geno2;

              bool ch_g2   = ch_geno1(gp2_iter.phased_geno(),
                                      gp1_iter.phased_geno(),
                                      !grandp1->is_female(),
                                      p1_iter.phased_geno(),
                                      al[a], con_geno2);
              double cont2 = cond_tran(u_iter.phased_geno(), con_geno2);
              double sibs2 = my_peeler->likelihood_of_sibs_except_sib(gp1_iter.phased_geno(), gp2_iter.phased_geno(), *parent1, *uncle).get_double();

              if( ha2 && ch_g2 && cont2 && sibs2 )
              {
                total +=   sub_t * sibs2 * ha2
                         * cond_tran(p1_iter.phased_geno(), con_geno2)
                         * cond_tran(u_iter.phased_geno(), con_geno2);
#if 0
  cout << "      ch_g2 = " << ch_g2 << endl
       << "      cont2 = " << cont2 << endl
       << "      sibs2 = " << sibs2 << endl;
#endif
              }
            } // end GP2
#if 0
  cout << "     gp2 total = " << total << endl;
#endif

          } // end GP1
#if 0
  cout << "    gp1 total = " << total << endl;
#endif

        } // end Parent
      }
    } // end Avunc
  } // end Child

#if 0
  cout << " * total = " << total << endl;
#endif

  my_ibds->set_ibd(current_pair, my_current_marker, 1.0 - total / my_subped_likelihood, QNAN, 0.0);
}

// grandparental pi^s

void pair_ibd_analysis::generate_gpar(size_t current_pair, const inheritance_model& model)
{
  //Variables:
  //  ant_grandp1,  ant_grandp2 - Grandparental anteriors.
  //                post_parent - Parent's postior without child's family.
  // post_grandp1, post_grandp2 - Grandparental posteriors without each other.
  //                      al[2] - Possible shared alleles.
  //                          a - allele counter.
  //                   cond_fam - Family of parent and child, conditioned on shared allele.
  //                       sibs - Prob. of sibs of parent given gparents genotypes.
  //                      total - Total likelihood.

  // Only process if the child phenotype isn't missing or error and at least
  //   one grandparent isn't missing or in error.  It only sets those
  //   grandparents whose phenotypes are not in error.

  filtered_relative_pair& grand_pair = my_pairs[current_pair];

  allele al[2];

  double total = 0.0;

  FPED::MemberConstPointer grandp1 = grand_pair.member_one;
  FPED::MemberConstPointer child   = grand_pair.member_two;
  FPED::MemberConstPointer parent1 = grand_pair.connector_one;

  // Figure out the second grandparent
  FPED::MemberConstPointer grandp2;

  if( grandp1 == my_meiosis_map.father(parent1) )
    grandp2 = my_meiosis_map.mother(parent1);
  else
    grandp2 = my_meiosis_map.father(parent1);

  // If the child is uninformative or both grandparents are uninformative,
  // skip this pair.
  //
  // XXXX ====> I'm not sure this is correct.  What if they are 'uninformative'
  //            but genotype elimination gives them a unique genotype?  Large
  //            pedigrees have this issue!  GCW 2006.08.17
//  if(    child->info().phenotype_missing(my_current_marker, model)
//      || (    grandp1->info().phenotype_missing(my_current_marker, model)
//           && grandp2->info().phenotype_missing(my_current_marker, model)) )
//    return;

  // Find the other pair.  We can calculate both of them at once.
  size_t other_pair = my_ibds->pair_index(grandp2, child);

  // If other pair is smaller, then this pair has already done so just exit
  if( other_pair < current_pair )
    return;

  // Get the other parent of the child, and determine who is the mother and who
  // the father
  FPED::MemberConstPointer parent2;
  if( parent1 == my_meiosis_map.father(child) )
    parent2 = my_meiosis_map.mother(child);
  else
    parent2 = my_meiosis_map.father(child);

#if 0
  cout << endl << "current_pair = " << current_pair << ": "
       << grandp1->name() << "," << child->name()
       << " -- parent = " << parent1->name()
       << ", " << parent2->name()
       << ", gp2 = " << grandp2->name()
       << ", other_pair = " << other_pair
       << endl;

  size_t p1_id = grandp1->info().phenotype(my_current_marker);
  size_t p2_id = grandp2->info().phenotype(my_current_marker);

  const inheritance_model& org_model = my_meiosis_map.get_multipedigree().info().marker_info(my_current_marker);
  cout << "gp1 pheno = " << org_model.get_phenotype(p1_id).name() << endl;
  cout << "gp2 pheno = " << org_model.get_phenotype(p2_id).name() << endl;
#endif

  // Get the phenotype indices of every important individual involved
  size_t grandp1_index = grandp1->subindex() + 1;
  size_t child_index   = child->subindex() + 1;
  size_t parent1_index = parent1->subindex() + 1;
  size_t grandp2_index = grandp2->subindex() + 1;

  phased_pen_iter c_iter = model.phased_penetrance_begin(child_index);
  for( ; c_iter != model.phased_penetrance_end(child_index); ++c_iter )
  {
    // If we're x-linked and the child is a male, skip genotypes that don't have
    // a Y-Chromosome
    if(    model.is_x_linked()
        && child->is_male()
        && !is_Y_genotype(c_iter.phased_geno()) )
      continue;

    phased_pen_iter gp1_iter = model.phased_penetrance_begin(grandp1_index);
    for( ; gp1_iter != model.phased_penetrance_end(grandp1_index); ++gp1_iter )
    {
      // If we're x-linked and the gp1 is a male, skip genotypes that don't have
      // a Y-Chromosome
      if(    model.is_x_linked()
          && grandp1->is_male()
          && !is_Y_genotype(gp1_iter.phased_geno()) )
        continue;

      // If there are no shared alleles between GP and child, we don't have to
      // proceed.
      if( !condition(c_iter.phased_geno(), gp1_iter.phased_geno(), al) )
        continue;

#if 0
  cout << "Child Genotype: " << c_iter.phased_geno().name() << endl;
  cout << "GP    Genotype: " << gp1_iter.phased_geno().name() << endl;

  cout << "allele shared : ";
  if( al[0] != allele() )
    cout << al[0].name();
  else
    cout << ".";
  cout << "/";
  if( al[1] != allele() )
    cout << al[1].name();
  else
    cout << ".";
  cout << endl;
#endif

      // Child and Granp potentially share an allele

      // Calculate GP's posterior and anterior
      double post_grandp1 = my_peeler->posterior_except_mate(*grandp1, *grandp2, gp1_iter).get_double();
      double ant_grandp1  = my_peeler->anterior(*grandp1, gp1_iter).get_double();

      // For each allele the potentially share
      for( size_t a = 0; a < 2 && al[a] != allele(); ++a )
      {
        phased_pen_iter p1_iter = model.phased_penetrance_begin(parent1_index);
        for( ; p1_iter != model.phased_penetrance_end(parent1_index); ++p1_iter )
        {
          // If xlinked and parent is the male, skip genotypes without a Y chromosome.
          // Note that there's some efficiency issues here that can be used, since
          // gp and ch can only share an allele if both are same sex, and if both
          // are male, can only share Y when parent is male too, and all this
          // calculation is unnecessary.  This should be incorporated.
          // XXX ==> GCW 2006.08.17
          if(    model.is_x_linked()
              && parent1->is_male()
              && !is_Y_genotype(p1_iter.phased_geno()) )
            continue;

          // Calculate the posterior of parent1 with parent 2, conditional on child
          // having allele al[a]
          double cond_fam = cond_post(parent1, parent2, child, p1_iter, c_iter, al[a], model);

          if( !has_allele(p1_iter.phased_geno(), al[a]) || !cond_fam )
            continue;

          // Parent shares the allele too and parent's conditional posterior
          // is != 0

          double post_parent = my_peeler->posterior_except_mate(*parent1, *parent2, p1_iter).get_double();
#if 0
  cout << "  cond_fam     = " << cond_fam    << endl
       << "  parent_post  = " << post_parent << endl
       << "  grandp1_post = " << post_grandp1  << endl
       << "  grandp1_ant  = " << ant_grandp1  << endl;
#endif
          phased_pen_iter gp2_iter = model.phased_penetrance_begin(grandp2_index);
          for( ; gp2_iter != model.phased_penetrance_end(grandp2_index); ++gp2_iter )
          {
            // If we're x-linked and the gp2 is a male, skip genotypes that don't have
            // a Y-Chromosome
            if(    model.is_x_linked()
                && grandp2->is_male()
                && !is_Y_genotype(gp2_iter.phased_geno()) )
              continue;

#if 0
  cout << "Par   Genotype: " << p1_iter.phased_geno().name() << endl;
  cout << "GP1   Genotype: " << gp1_iter.phased_geno().name() << endl;
  cout << "GP2   Genotype: " << gp2_iter.phased_geno().name() << endl;
#endif
            conditional_genotype con_geno;

            bool ch_g  = ch_geno1(gp1_iter.phased_geno(),
                                  gp2_iter.phased_geno(),
                                  grandp1->is_female(),
                                  p1_iter.phased_geno(),
                                  al[a], con_geno);
                                  
            if( ch_g )
            {
              double ant_grandp2  = my_peeler->anterior(*grandp2, gp2_iter).get_double();
              double post_grandp2 = my_peeler->posterior_except_mate(*grandp2, *grandp1, gp2_iter).get_double();

              double cont = cond_tran(p1_iter.phased_geno(), con_geno);
              double sibs;
              if( grandp1 == my_meiosis_map.mother(parent1) )
                sibs = my_peeler->likelihood_of_sibs(gp1_iter.phased_geno(), gp2_iter.phased_geno(), *parent1).get_double();
              else
                sibs = my_peeler->likelihood_of_sibs(gp2_iter.phased_geno(), gp1_iter.phased_geno(), *parent1).get_double();
#if 0
  cout << endl;
  
  cout << "    gpant2  = " << ant_grandp2  << endl
       << "    gppost2 = " << post_grandp2 << endl;
  cout << "    sibs    = " << sibs << endl;
  cout << "    cont    = " << cont << endl;
#endif

                double tot = 
                         ant_grandp1 * post_grandp1 * ant_grandp2 * post_grandp2
                       * (*gp1_iter)
                       * (*gp2_iter)
                       * (*p1_iter)
                       * (*c_iter)
                       * has_allele(gp1_iter.phased_geno(), al[a])
                       * cont * cond_fam * post_parent * sibs;
                total += tot;
            }
          } // end GP2
        } // end Parent
      }
    } // end GP1
  } // end Grandchild

#if 0
  cout << " * total = " << total << endl
       << "   pr0   = " << 1.0 - total / my_subped_likelihood << endl;
#endif

  // This shouldn't happen!!!!!!!!   Numerical issues?? XXXXX
  if( total > my_subped_likelihood )
    total = my_subped_likelihood;

//  if( !grandp1->info().phenotype_missing(my_current_marker, model) )
    my_ibds->set_ibd(current_pair, my_current_marker, 1.0 - total / my_subped_likelihood, QNAN, 0.0);

//  if( !grandp2->info().phenotype_missing(my_current_marker, model) )
    my_ibds->set_ibd(other_pair, my_current_marker, total / my_subped_likelihood, QNAN, 0.0);
}

// Cousin pi^s

void pair_ibd_analysis::generate_cous(size_t current_pair, const inheritance_model& model)
{
  //Variables:
  //  ant_grandp1,  ant_grandp2 - Grandparental anteriors.
  // post_parent1, post_parent2 - Parent's postior without children's family.
  // post_grandp1, post_grandp2 - Grandparent's posteriors without current family.
  //                      al[2] - Possible shared alleles.
  //                          a - current allele_id.
  //                   cond_fam - Family of parent and child conditioned for allele.
  //                       sibs - Posterior Probability of sibs of parent and avunc.
  //                      total - Likelihood after conditioning.
  //                      sub_t - Partial multiplication of some terms to avoid extra work
  //          ha1,          ha2 - Granparent n has the shared allele.

  filtered_relative_pair& cousin_pair = my_pairs[current_pair];

//  if(    cousin_pair.member_one->info().phenotype_missing(my_current_marker, model)
//      || cousin_pair.member_two->info().phenotype_missing(my_current_marker, model) )
//    return;

  allele al[2];

  double total = 0.0;

  fmember_const_pointer cousin1 = cousin_pair.member_one;
  fmember_const_pointer cousin2 = cousin_pair.member_two;

  fmember_const_pointer parent1 = cousin_pair.connector_one;
  fmember_const_pointer parent2 = cousin_pair.connector_two;

  fmember_const_pointer grandp1 = my_meiosis_map.mother(parent1);
  fmember_const_pointer grandp2 = my_meiosis_map.father(parent1);

  fmember_const_pointer spouse1;
  if( parent1 == my_meiosis_map.father(cousin1) )
    spouse1 = my_meiosis_map.mother(cousin1);
  else
    spouse1 = my_meiosis_map.father(cousin1);

  fmember_const_pointer spouse2;
  if( parent2 == my_meiosis_map.father(cousin2) )
    spouse2 = my_meiosis_map.mother(cousin2);
  else
    spouse2 = my_meiosis_map.father(cousin2);

#if 0
  cout << endl << "current_pair = " << current_pair << ": "
       << cousin1->name() << "," << cousin2->name()
       << " -- parent of 1 = " << parent1->name() << ", " << spouse1->name()
       << " -- parent of 2 = " << parent2->name() << ", " << spouse2->name()
       << ", gp1 = " << grandp1->name()
       << ", gp2 = " << grandp2->name() << endl;
#endif

  size_t cousin1_index = cousin1->subindex() + 1;
  size_t cousin2_index = cousin2->subindex() + 1;
  size_t parent1_index = parent1->subindex() + 1;
  size_t parent2_index = parent2->subindex() + 1;
  size_t grandp1_index = grandp1->subindex() + 1;
  size_t grandp2_index = grandp2->subindex() + 1;
  //size_t spouse1_index = spouse1->subindex() + 1;
  //size_t spouse2_index = spouse2->subindex() + 1;

  phased_pen_iter c1_iter = model.phased_penetrance_begin(cousin1_index);
  for( ; c1_iter != model.phased_penetrance_end(cousin1_index); ++c1_iter )
  {
    if(    model.is_x_linked()
        && cousin1->is_male()
        && !is_Y_genotype(c1_iter.phased_geno()) )
      continue;

    phased_pen_iter c1i_iter(c1_iter);

    phased_pen_iter c2_iter = model.phased_penetrance_begin(cousin2_index);
    for( ; c2_iter != model.phased_penetrance_end(cousin2_index); ++c2_iter )
    {
      if(    model.is_x_linked()
          && cousin2->is_male()
          && !is_Y_genotype(c2_iter.phased_geno()) )
        continue;

      if( !condition(c1_iter.phased_geno(), c2_iter.phased_geno(), al) )
        continue;

      phased_pen_iter c2i_iter(c2_iter);
#if 0
  cout << "allele shared : ";
  if( al[0] != allele() )
    cout << al[0].name();
  else
    cout << ".";
  cout << "/";
  if( al[1] != allele() )
    cout << al[1].name();
  else
    cout << ".";
  cout << endl;
#endif

      // Cousins share an allele

      for( size_t a = 0; a < 2 && al[a] != allele(); ++a )
      {
        // First parent
        phased_pen_iter p1_iter = model.phased_penetrance_begin(parent1_index);
        for( ; p1_iter != model.phased_penetrance_end(parent1_index); ++p1_iter )
        {
          if(    model.is_x_linked()
              && parent1->is_male()
              && !is_Y_genotype(p1_iter.phased_geno()) )
            continue;

          phased_pen_iter p1i_iter(p1_iter);

          double cond_fam1 = cond_post(parent1, spouse1, cousin1, p1i_iter, c1i_iter, al[a], model);

          if( !has_allele(p1_iter.phased_geno(), al[a]) || !cond_fam1 )
            continue;

          // First parent has the allele and shares it with its child

          double post_parent1 = my_peeler->posterior_except_mate(*parent1, *spouse1, p1i_iter).get_double();

          // Second parent
          phased_pen_iter p2_iter = model.phased_penetrance_begin(parent2_index);
          for( ; p2_iter != model.phased_penetrance_end(parent2_index); ++p2_iter )
          {
            if(    model.is_x_linked()
                && parent2->is_male()
                && !is_Y_genotype(p2_iter.phased_geno()) )
              continue;

            phased_pen_iter p2i_iter(p2_iter);

            double cond_fam2 = cond_post(parent2, spouse2, cousin2, p2i_iter, c2i_iter, al[a], model);

            if( !has_allele(p2_iter.phased_geno(), al[a]) || !cond_fam2 )
              continue;

            // Second parent has the allele and shares it with its child

            double post_parent2 = my_peeler->posterior_except_mate(*parent2, *spouse2, p2i_iter).get_double();

            phased_pen_iter gp1_iter = model.phased_penetrance_begin(grandp1_index);
            for( ; gp1_iter != model.phased_penetrance_end(grandp1_index); ++gp1_iter )
            {
              if(    model.is_x_linked()
                  && grandp1->is_male()
                  && !is_Y_genotype(gp1_iter.phased_geno()) )
                continue;

              phased_pen_iter gp1i_iter(gp1_iter);

              double ant_grandp1  = my_peeler->anterior(*grandp1, gp1i_iter).get_double();
              double post_grandp1 = my_peeler->posterior_except_mate(*grandp1, *grandp2, gp1i_iter).get_double();

              phased_pen_iter gp2_iter = model.phased_penetrance_begin(grandp2_index);
              for( ; gp2_iter != model.phased_penetrance_end(grandp2_index); ++gp2_iter )
              {
                if(    model.is_x_linked()
                    && grandp2->is_male()
                    && !is_Y_genotype(gp2_iter.phased_geno()) )
                  continue;

                size_t ha1 = has_allele(gp1_iter.phased_geno(), al[a]);
                size_t ha2 = has_allele(gp2_iter.phased_geno(), al[a]);

                if( !(ha1 || ha2) ) continue;

                phased_pen_iter gp2i_iter(gp2_iter);

                // At least one Grandparent has the shared allele

                double ant_grandp2  = my_peeler->anterior(*grandp2, gp2i_iter).get_double();
                double post_grandp2 = my_peeler->posterior_except_mate(*grandp2, *grandp1, gp2i_iter).get_double();

                double sub_t =   ant_grandp1 * post_grandp1 * ant_grandp2 * post_grandp2
                               * (*gp2i_iter)
                               * (*gp1i_iter)
                               * (*p2i_iter)
                               * (*p1i_iter)
                               * (*c2i_iter)
                               * (*c1i_iter)
                               * cond_fam2 * cond_fam1 * post_parent2 * post_parent1;

                conditional_genotype con_geno1;

                bool ch_g1   = ch_geno1(gp1_iter.phased_geno(),
                                        gp2_iter.phased_geno(),
                                        grandp1->is_female(),
                                        p1_iter.phased_geno(),
                                        al[a], con_geno1);
                double cont1 = cond_tran(p2_iter.phased_geno(), con_geno1);
                double sibs1 = my_peeler->likelihood_of_sibs_except_sib(gp1_iter.phased_geno(), gp2_iter.phased_geno(), *parent1, *parent2).get_double();

                if( ha1 && ch_g1 && cont1 && sibs1 )
                {
                  total +=   sub_t * sibs1 * ha1
                           * cond_tran(p1_iter.phased_geno(), con_geno1)
                           * cond_tran(p2_iter.phased_geno(), con_geno1);
                }

                conditional_genotype con_geno2;

                bool ch_g2   = ch_geno1(gp2_iter.phased_geno(),
                                        gp1_iter.phased_geno(),
                                        !grandp1->is_female(),
                                        p1_iter.phased_geno(),
                                        al[a], con_geno2);
                double cont2 = cond_tran(p2_iter.phased_geno(), con_geno2);
                double sibs2 = my_peeler->likelihood_of_sibs_except_sib(gp1_iter.phased_geno(), gp2_iter.phased_geno(), *parent1, *parent2).get_double();

                if( ha2 && ch_g2 && cont2 && sibs2 )
                {
                  total +=   sub_t * sibs2 * ha2
                           * cond_tran(p1_iter.phased_geno(), con_geno2)
                           * cond_tran(p2_iter.phased_geno(), con_geno2);
                }
              } // end GP2
            } // end GP1
          } // end Parent 2
        } // end Parent 1
      }
    } // end Child 2
  } // end Child 1

#if 0
  cout << " * total = " << total << endl;
#endif

  my_ibds->set_ibd(current_pair, my_current_marker, 1.0 - total / my_subped_likelihood, QNAN, 0.0);
}

// cond_post calculates the posterior probability of a person for a given
//   genotype, conditional upon a child having a certain genotype, and a
//   certain allele is shared between them (ie, the posterior of a person
//   with a child with a certain genotype that shares an allele with the
//   parent).  Used for all pair types except siblings.

// Note:  If the parent is homozygous, this calculates wrong, but there are
//        adjustment factors in the generate functions to account for this.
//
double
pair_ibd_analysis::cond_post(fmember_const_pointer    mp,
                             fmember_const_pointer    spouse,
                             fmember_const_pointer    child,
                             const phased_pen_iter&   mpi_iter,
                             const phased_pen_iter&   child_iter,
                             const allele&            al,
                             const inheritance_model& model) const
{
#if 0
  cout << "cond_post(" << mp->name() << ", " << spouse->name()
       << ", " << child->name() << ", ";
  if( al != allele() )
    cout << al.name();
  else
    cout << ".";
  cout << ")" << endl;
#endif

  double               result = 0.0;
  conditional_genotype con_geno;
 
  const MLOCUS::phased_genotype& mp_geno = mpi_iter.phased_geno();
  const MLOCUS::phased_genotype& ch_geno = child_iter.phased_geno();

  size_t spouse_index = spouse->subindex() + 1;

  phased_pen_iter sp_iter = model.phased_penetrance_begin(spouse_index);
  for( ; sp_iter != model.phased_penetrance_end(spouse_index); ++sp_iter )
  {
    if( !ch_geno1(mp_geno,
                  sp_iter.phased_geno(),
                  (spouse != my_meiosis_map.mother(child)),
                  ch_geno, al, con_geno) )
      continue;

    double p2ant     = my_peeler->anterior(*spouse, sp_iter).get_double();
    double p2postexc = my_peeler->posterior_except_mate(*spouse, *mp, sp_iter).get_double();

    double sibprob;
    if( spouse == my_meiosis_map.mother(child) )
      sibprob = my_peeler->likelihood_of_sibs(ind_genotype(sp_iter), ind_genotype(mpi_iter), *child).get_double();
    else
      sibprob = my_peeler->likelihood_of_sibs(ind_genotype(mpi_iter), ind_genotype(sp_iter), *child).get_double();

    double penet     = (*sp_iter);
    double condtran  = cond_tran(ch_geno, con_geno);
    double cpost     = my_peeler->posterior(*child, child_iter).get_double();

    result +=   p2ant * p2postexc * sibprob * penet * condtran * cpost;

#if 0
  cout << " p2ant     = " << p2ant << endl
       << " p2postexc = " << p2postexc << endl
       << " sibprob   = " << sibprob << endl
       << " penet     = " << penet << endl
       << " condtran  = " << condtran << endl
       << " cpost     = " << cpost << endl;
#endif
  }

#if 0
  cout << " * result = " << result << endl;
#endif

  return result;
}

// ch_geno (children's genotypes) calculates the places the allele a passes
// to from the first parent.  The returned bools have true for those that
// have the Allele a from the first parent, and false for those who don't.
// If both alleles from the first parent include a, then it only keeps track
// of the first one, and the value is adjusted accordingly later.

bool
pair_ibd_analysis::ch_geno1
    (const phased_genotype& pg1,
     const phased_genotype& pg2,
     bool                   first_is_mother,
     const phased_genotype& cg,
     const allele&          al,
     conditional_genotype&  con_geno) const
{
  if(first_is_mother)
  {
    con_geno.child_geno = child_genotype_set(pg1, pg2);

    if( al == pg1.allele1() )
    {
      con_geno.b[0] = con_geno.b[1] = true;
      con_geno.b[2] = con_geno.b[3] = false;
  
      return (    con_geno.child_geno[0] == cg
               || con_geno.child_geno[1] == cg );
    }
    else
    {
      con_geno.b[2] = con_geno.b[3] = true;
      con_geno.b[0] = con_geno.b[1] = false;
  
      return (    con_geno.child_geno[2] == cg
               || con_geno.child_geno[3] == cg );
    }
  }
  else
  {
    con_geno.child_geno = child_genotype_set(pg2, pg1);

    if( al == pg1.allele1() )
    {
      con_geno.b[0] = con_geno.b[2] = true;
      con_geno.b[1] = con_geno.b[3] = false;
  
      return (    con_geno.child_geno[0] == cg
               || con_geno.child_geno[2] == cg );
    }
    else
    {
      con_geno.b[1] = con_geno.b[3] = true;
      con_geno.b[0] = con_geno.b[2] = false;
  
      return (    con_geno.child_geno[1] == cg
               || con_geno.child_geno[3] == cg );
    }
  }
}

/// condition calculates the alleles we need to condition on.  It returns a
///   bool indicating wether there are any, and sets al to the one or two
///   alleles shared by g1 and g2.

bool
pair_ibd_analysis::condition(const phased_genotype& g1, const phased_genotype& g2, allele* al) const
{
#if 0
  cout << "condition("
       << g1.allele1().name() << "/" << g1.allele2().name() << ", "
       << g2.allele1().name() << "/" << g2.allele2().name() << ", "
       << ")" << endl;
#endif

  al[0] = allele();
  al[1] = allele();

  // If the first allele of both people is the same...
  if( g1.allele1() == g2.allele1() )
  {
    al[0] = g1.allele1();

    // If the second alleles are the same, but different from the first...
    if( g1.allele2() == g2.allele2() && g1.allele2() != g1.allele1() )
      al[1] = g1.allele2();

#if 0
  cout << "return true " << al[0].name() << "/" << al[1].name() << endl;
#endif

    return true;
  }
  else if( g1.allele1() == g2.allele2() )
  {
    al[0] = g1.allele1();
 
    if( g1.allele2() == g2.allele1() )
      al[1] = g1.allele2();

#if 0
  cout << "return true " << al[0].name() << "/" << al[1].name() << endl;
#endif

    return true;
  }

  if( g1.allele2() == g2.allele1() ) al[0] = g1.allele2();
  if( g1.allele2() == g2.allele2() ) al[0] = g1.allele2();

#if 0
  cout << "return " << al[0].name() << "/" << al[1].name() << endl;
#endif

  if( al[0] != allele() ) return true;

  return false;
}

} // end of namespace GENIBD

} // end of namespace SAGE
