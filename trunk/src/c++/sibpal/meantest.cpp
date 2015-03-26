#include "sibpal/meantest.h"

using namespace std;

#define DEBUG_MEAN(x)

namespace SAGE   {
namespace SIBPAL {

SibMeanTest::SibMeanTest(relative_pairs& p, cerrorstream& err)
           : pairs(p), errors(err)
{}

void
SibMeanTest::regress()
{
  if( !parameters().get_use_full_sibs() && !parameters().get_use_half_sibs() )
    return;
  
  // FIXME: Make the check for a valid trait better encapsulated
  if(   parameters().trait().trait != (size_t)-1 
     && parameters().trait().trait >= pairs.trait_count() )
    return; 

  int a = parameters().trait().affected_types();

  // FIXME: We can only perform this test when there is no trait
  //        or when we are looking at only one affected pair type.
  if( a != 3 && a != 1)
    return;

  GLS3 gls(2);

  for(size_t i=0; i < parameters().marker_count(); ++i)
  {
    gls.reset();

    size_t m = parameters()[i].marker;

    parameters()[i].affecteds = parameters().trait().affected_count();
    parameters()[i].estimate.clear();

    if( parameters().get_use_full_sibs() && !parameters().get_use_half_sibs() )
      parameters()[i].estimate.set_w(parameters().get_w());

    regress_marker(gls, m);

    parameters()[i].pair_count = gls.observation_count;

    if(!gls.beta)
      continue;

    parameters()[i].estimate.set_pi(gls.beta(0,0));
    parameters()[i].estimate.set_f1(gls.beta(1,0));

    if(!gls.Variance)
      continue;

#if 0
  cout << "gls beta : " << endl;
  print_matrix(gls.beta, cout);
  cout << "gls Var : " << endl;
  print_matrix(gls.Variance, cout);
#endif

    parameters()[i].estimate.set_cov(gls.Variance);
  }
  parameters().validate();
}

void
SibMeanTest::regress_marker(GLS3& gls, size_t marker)
{
  gls.reset(1,2);

  matrix A;                    // design matrix
  matrix y;                    // trait vector

  pair_filter pf;
  pf.add_marker(marker);

  if( parameters().trait().trait != (size_t)-1 )
  {
    pf.add_trait( parameters().trait().trait,
                  parameters().trait().affection );
  }

  for(size_t i = 0; i < parameters().subset_count(); ++i)
  {
    pf.add_trait( parameters().subsets(i).trait,
                  parameters().subsets(i).affection,
                  std::numeric_limits<double>::epsilon() );
  }

  vector<sib_cluster> my_sib_clusters;
  bool use_fsib = get_use_pairs().first;
  bool use_hsib = get_use_pairs().second;

  // Iterate over all sib_clusters
  sibship_cluster_const_iterator c_iter = pairs.sibship_cluster_begin();
  sibship_cluster_const_iterator c_end  = pairs.sibship_cluster_end();

  for( ; c_iter != c_end; ++c_iter )
  {
    if(    use_fsib && !use_hsib
        && c_iter->full_sibship_map.size() )
    {
      const vector<size_t>& hsib = c_iter->hsib_pairs;

      map< id_pair, vector<size_t> >::const_iterator mi = c_iter->full_sibship_map.begin();
      for( ; mi != c_iter->full_sibship_map.end(); ++mi )
      {
        const vector<size_t>& fsib = mi->second;

        sib_cluster sc(fsib, hsib, pairs, use_fsib, use_hsib, pf);

        if( !sc.valid_pair_count() ) continue;

        my_sib_clusters.push_back(sc);
      }
    }
    else
    {
      sib_cluster sc(c_iter->fsib_pairs, c_iter->hsib_pairs, pairs, use_fsib, use_hsib, pf);

      if( !sc.valid_pair_count() ) continue;

      my_sib_clusters.push_back(sc);
    }
  }

#if 0
  for( size_t i = 0; i < my_sib_clusters.size(); ++i )
  {
    const sib_cluster& sc = my_sib_clusters[i];

    cout << endl << "sib cluster " << i << " ";
    sc.dump();

    cout << " All valid pairs:" << endl;
    for( size_t j = 0; j < sc.valid_pair_count(); ++j )
    {
      const sib_pair this_pair = sc[j];
      cout << "  " << i << " " << this_pair.pair_number()
           << "(" << this_pair.sibs().first->name() 
           << "," << this_pair.sibs().second->name()
           << ")";
    }
    cout << endl;
  }
  cout << endl;
#endif
                                                                                     

  for( size_t c = 0; c < my_sib_clusters.size(); ++c )
  {
    const sib_cluster& sc = my_sib_clusters[c];

    size_t n = sc.valid_pair_count();

    // Compute basis matrices
    A.resize_fill(n,1,1.0);

    if( parameters().get_use_full_sibs() && !parameters().get_use_half_sibs() )
      marker_matrix(sc, y, marker, 0.0, parameters().get_w(), 1.0);
    else
      marker_matrix(sc, y, marker, 0.0, 0.5, 1.0);

#if 0
  cout << "y : " << endl;
  print_matrix(y, cout);
  cout << "A : " << endl;
  print_matrix(A, cout);
#endif

    gls.add_block(y, A);
  }

  if(gls.observation_count - 1 <= 0)
    errors << priority(warning) << "Not enough valid pairs." << endl;

  gls.compute();
}

} //end of namespace SIBPAL
} //end of namespace SAGE
