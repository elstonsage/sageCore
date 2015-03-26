#include "sibpal/sib_cluster.h"

namespace SAGE   {
namespace SIBPAL {

sib_cluster::sib_cluster()
{}

sib_cluster::sib_cluster(const vector<size_t>& fsib_pairs,
                         const vector<size_t>& hsib_pairs,
                         relative_pairs& rp,
                         bool fs, bool hs,
                         const pair_filter& p,
                         bool x_linked, bool mm, bool mf, bool ff)
           : my_filter(&p), my_data(&rp)
{
  my_valid_fsib_pairs.resize(0);
  my_valid_hsib_pairs.resize(0);

  my_sibs.clear();

  size_t f_start = 0;
  size_t f_stop  = 0;

  if( fs && fsib_pairs.size() )
  {
    f_start = fsib_pairs[0];
    f_stop  = fsib_pairs[fsib_pairs.size()-1] + 1;

    //cout << "f_start = " << f_start << endl
    //     << "f_stop  = " << f_stop  << endl;

    for( size_t i = 0; i < fsib_pairs.size(); ++i )
    {
      sib_pair this_pair(fsib_pairs[i], &rp);
    
      if( my_filter && my_filter->valid(this_pair) )
      {
        if( mm && this_pair.is_mm_pair() )
        {
          my_valid_fsib_mm_pairs.push_back(fsib_pairs[i]);

          my_valid_fsib_pairs.push_back(fsib_pairs[i]);

          my_sibs.insert( this_pair.rels().pair.first  );
          my_sibs.insert( this_pair.rels().pair.second );
        }
        else if( mf && this_pair.is_mf_pair() )
        {
          my_valid_fsib_mf_pairs.push_back(fsib_pairs[i]);

          my_valid_fsib_pairs.push_back(fsib_pairs[i]);

          my_sibs.insert( this_pair.rels().pair.first  );
          my_sibs.insert( this_pair.rels().pair.second );
        }
        else if( ff && this_pair.is_ff_pair() )
        {
          my_valid_fsib_ff_pairs.push_back(fsib_pairs[i]);

          my_valid_fsib_pairs.push_back(fsib_pairs[i]);

          my_sibs.insert( this_pair.rels().pair.first  );
          my_sibs.insert( this_pair.rels().pair.second );
        }
        else if( !x_linked )
        {
          my_valid_fsib_pairs.push_back(fsib_pairs[i]);

          my_sibs.insert( this_pair.rels().pair.first  );
          my_sibs.insert( this_pair.rels().pair.second );
        }
      }
    }
  }

  size_t h_start = 0;
  size_t h_stop  = 0;

  if( hs && hsib_pairs.size() )
  {
    h_start = hsib_pairs[0];
    h_stop  = hsib_pairs[hsib_pairs.size()-1] + 1;

    //cout << "h_start = " << h_start << endl
    //     << "h_stop  = " << h_stop  << endl;

    for( size_t i = 0; i < hsib_pairs.size(); ++i )
    {
      sib_pair this_pair(hsib_pairs[i], &rp);

      if( my_filter && my_filter->valid(this_pair) )
      {
        my_valid_hsib_pairs.push_back(hsib_pairs[i]);

        my_sibs.insert( this_pair.rels().pair.first  );
        my_sibs.insert( this_pair.rels().pair.second );
      }
    }
  }
}

void
sib_cluster::dump() const
{
  cout << endl << "=== sib_cluster dump:" << endl;
  
  cout << endl
       << "pair count, valid = " << valid_pair_count() << endl
       << "valid   pair count = " << valid_fsib_pair_count() << endl
       << "valid h pair count = " << valid_hsib_pair_count() << endl
       << "valid sib count = " << valid_sib_count() << endl;

  cout << "valid_fsib_pairs size = " << my_valid_fsib_pairs.size() << endl;

  for( size_t j = 0; j < my_valid_fsib_pairs.size(); ++j )
    cout << "  " << j << ":" << my_valid_fsib_pairs[j];
  cout << endl;

  cout << "brother-brother pairs size = " << my_valid_fsib_mm_pairs.size() << endl;

  for( size_t j = 0; j < my_valid_fsib_mm_pairs.size(); ++j )
    cout << "  " << j << ":" << my_valid_fsib_mm_pairs[j];
  cout << endl;

  cout << "brother-sister pairs size = " << my_valid_fsib_mf_pairs.size() << endl;

  for( size_t j = 0; j < my_valid_fsib_mf_pairs.size(); ++j )
    cout << "  " << j << ":" << my_valid_fsib_mf_pairs[j];
  cout << endl;

  cout << "sister-sister pairs size = " << my_valid_fsib_ff_pairs.size() << endl;

  for( size_t j = 0; j < my_valid_fsib_ff_pairs.size(); ++j )
    cout << "  " << j << ":" << my_valid_fsib_ff_pairs[j];
  cout << endl;

  cout << "valid_hsib_pairs size = " << my_valid_hsib_pairs.size() << endl;

  for( size_t j = 0; j < my_valid_hsib_pairs.size(); ++j )
    cout << "  " << j << ":" << my_valid_hsib_pairs[j];
  cout << endl;

  cout << "h.brother-h.brother pairs size = " << my_valid_hsib_mm_pairs.size() << endl;

  for( size_t j = 0; j < my_valid_hsib_mm_pairs.size(); ++j )
    cout << "  " << j << ":" << my_valid_hsib_mm_pairs[j];
  cout << endl;

  cout << "h.brother-h.sister pairs size = " << my_valid_hsib_mf_pairs.size() << endl;

  for( size_t j = 0; j < my_valid_hsib_mf_pairs.size(); ++j )
    cout << "  " << j << ":" << my_valid_hsib_mf_pairs[j];
  cout << endl;

  cout << "h.sister-h.sister pairs size = " << my_valid_hsib_ff_pairs.size() << endl;

  for( size_t j = 0; j < my_valid_hsib_ff_pairs.size(); ++j )
    cout << "  " << j << ":" << my_valid_hsib_ff_pairs[j];
  cout << endl;
}

} // end namespace SIBPAL
} // end namespace SAGE
