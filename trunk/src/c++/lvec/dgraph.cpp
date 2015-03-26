#include "lvec/dgraph.h"

namespace SAGE
{

void descent_graph::move_to(equivalence_class e, const meiosis_map& m)
{
  if(&m != mm)
  {
    mm = &m;

    av.resize(mm->get_subpedigree()->member_count());
  }

  allele_type allele = 0;

  for(size_t i = 0; i < mm->get_subpedigree()->member_count(); ++i)
  {
    meiosis_map::size_type mloc = mm->mother_index(i);

    if(mloc == (index) -1)
    {
      av[i].first  = allele++;
      av[i].second = allele++;
    }
    else
    {
      meiosis_map::size_type floc = mm->father_index(i);

      meiosis_map::index mmei = mm->mother_meiosis(i);
      meiosis_map::index fmei = mm->father_meiosis(i);

      int p1 = (mmei < SAGE::meiosis_map::meiosis_bits) ?
               0 : e & (1 << (mmei - SAGE::meiosis_map::meiosis_bits));
      int p2 = (fmei < SAGE::meiosis_map::meiosis_bits) ?
               0 : e & (1 << (fmei - SAGE::meiosis_map::meiosis_bits));

      av[i].first  = p1 ? av[mloc].second : av[mloc].first;
      av[i].second = p2 ? av[floc].second : av[floc].first;
    }
  }
}

void descent_graph::print_allele_vector() const
{
  cout << "dgraph::allele_vector.." << endl;
  for( size_t i = 0; i < av.size(); ++i )
  {
    cout << i << "	" << av[i].first << "	" << av[i].second << endl;
  }
}

}
