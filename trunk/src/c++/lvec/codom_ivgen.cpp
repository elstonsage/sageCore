#include <algorithm>
#include <cassert>
#include "mlocus/imodel.h"
#include "lvec/inheritance_vector.h"
#include "lvec/codom_ivgen.h"
#include "lvec/node.h"

namespace SAGE
{

vector<edge> edges;
vector<node> nodes;

bool codominant_iv_generator::build
    (const SAGE::meiosis_map& mm, const SAGE::MLOCUS::inheritance_model& im,
     bool use_pop_freq) const
{
  use_pf = use_pop_freq;

  my_mmap   = mm;
  my_imodel = im;

  if(   !my_mmap.get_subpedigree()->member_count()
      || my_mmap.get_subpedigree()->member_count() != my_imodel.phenotype_count() )
    return false;

  ped_id pid = my_mmap.get_subpedigree();

  if(!pid || !my_imodel.codominant(false))
    return false;

  my_inds.resize(my_mmap.get_subpedigree()->member_count());

  edges.resize(0);
  nodes.resize(0);

  edges.reserve(my_mmap.get_subpedigree()->member_count());
  nodes.reserve(my_mmap.founder_count() * 2);

  my_penetrance = 1.0;
  
  for(size_t j = 0; j < my_mmap.get_subpedigree()->member_count(); ++j)
  {
    uint pheno_id = my_mmap.get_subpedigree()->member_index(j).subindex() + 1;
//    uint pheno_id = j + 1;

    my_inds[j].process     = false;
    my_inds[j].edge_number = -1;

    size_t upc = my_imodel.unphased_penetrance_count(pheno_id);
    size_t ppc = my_imodel.phased_penetrance_count(pheno_id);

#if 0
  cout << "j = " << j << ", pheno_id = " << pheno_id
       << ", name = " << my_mmap.get_subpedigree()->member_index(j).name()
       << ", upc = " << upc  << ", ppc = " << ppc
       << ", pheno = " << my_imodel.get_phenotype(pheno_id).name() << endl;
#endif

    if(upc < 2)
    {
      my_inds[j].edge_number = edges.size();

      SAGE::MLOCUS::inheritance_model::phased_penetrance_iterator ppi =
          my_imodel.phased_penetrance_begin(pheno_id);

      SAGE::MLOCUS::phased_genotype g = ppi.phased_geno();
    
      my_penetrance *= *ppi;

      bool ph = (ppc == 1);

      edges.push_back(
         edge(my_mmap[j], make_pair(g.allele1(), g.allele2()), ph));

      back(j);
    }

    if(my_mmap.founder(j))
    {
      my_inds[j].node1 = nodes.size();
      my_inds[j].node2 = nodes.size() + 1;

      nodes.resize(nodes.size() + 2);

      nodes[nodes.size() - 2] = node(nodes.size() - 2);
      nodes[nodes.size() - 1] = node(nodes.size() - 1);
      
      if(my_imodel.get_model_type() != MLOCUS::AUTOSOMAL)
      {
        SAGE::MLOCUS::inheritance_model::phased_penetrance_iterator ppi =
            my_imodel.phased_penetrance_begin(pheno_id);

        SAGE::MLOCUS::phased_genotype g = ppi.phased_geno();
        
        nodes[nodes.size() - 2].set_sex_null_state(g.allele1().is_sex_allele());
        nodes[nodes.size() - 1].set_sex_null_state(g.allele2().is_sex_allele());
      }

      if(my_inds[j].edge_number != -1)
      {
        nodes[my_inds[j].node1].add_edge(edges[my_inds[j].edge_number], 0);
        nodes[my_inds[j].node2].add_edge(edges[my_inds[j].edge_number], 1);
      }
    }
  }

  // This has to wait until all individuals are processed
  for(size_t j = 0; j < my_mmap.get_subpedigree()->member_count(); ++j)
  {
    if( !my_inds[j].process && !my_mmap.founder(j) )
    {
      my_hflags.push_back(my_mmap.mother_meiosis(j));
      my_hflags.push_back(my_mmap.father_meiosis(j));
    }
  }

#if 0

  cout << "Number of individuals: " << my_mmap.get_subpedigree()->member_count() << endl;
  cout << "name\tmeiosis\tprocess\tedge_no\tnode1\tnode2" << endl;

  for(size_t j = 0; j < my_mmap.get_subpedigree()->member_count(); ++j)
  {
    cout << my_mmap.member(j)->name() << '\t';
    if(!my_mmap.founder(j))
      cout << my_mmap.mother_meiosis(j) << ' ' << my_mmap.father_meiosis(j);

    cout << '\t' << my_inds[j].process
         << '\t' << my_inds[j].edge_number
         << '\t' << my_inds[j].node1
         << '\t' << my_inds[j].node2
         << endl;
  }

  cout << "my_penetrance = " << my_penetrance << endl;
  for( list<size_t>::const_iterator h = my_hflags.begin(); h != my_hflags.end(); ++h )
    cout << *h << "	";
  cout << endl;

  cout << "1. my_ivector :" << endl;
  cout << "     nonf = " << my_ivector.get_nonfounders() << endl
       << "        f = " << my_ivector.get_founders() << endl;

  cout << "NODES: " << nodes.size() << endl;
  cout << "i\tstate\tfixed\tch\ta[0]/a[1]/f_a\tconnection edges_name.." << endl;
  for(size_t a = 0; a < nodes.size(); ++a)
    nodes[a].dump(cout, a);

  cout << "EDGES: " << edges.size() << endl;
  cout << "id()\tstate\ta[0]/a[1]\tch\tn0+n1\tnode_ids.." << endl;
  for(size_t a = 0; a < edges.size(); ++a)
    edges[a].dump(cout);

#endif

  my_ivector(my_mmap);
  
#if 0
  cout << "2. my_ivector :" << endl;
  cout << "     nonf = " << my_ivector.get_nonfounders() << endl
       << "        f = " << my_ivector.get_founders() << endl;
#endif

  build_ind(0);

  for(size_t j = 0; j < my_mmap.get_subpedigree()->member_count(); ++j)
  {
    if(my_mmap.founder(j) && my_inds[j].edge_number != -1)
    {
      nodes[my_inds[j].node2].remove_edge();
      nodes[my_inds[j].node1].remove_edge();
    }
  }

#if 0
  cout << "3. my_ivector :" << endl;
  cout << "     nonf = " << my_ivector.get_nonfounders() << endl
       << "        f = " << my_ivector.get_founders() << endl;
#endif

  return true;
}

void codominant_iv_generator::back(int j) const
{
  if(my_inds[j].process || my_mmap.founder(j)) return;

  my_inds[j].process = true;

  back(my_mmap.mother_index(j));
  back(my_mmap.father_index(j));
}

void codominant_iv_generator::build_ind(size_t i) const
{
  if(i == my_mmap.get_subpedigree()->member_count()) { eval_ind(0, 1.0); return; }

#if 0
  cout << "Building Individual " << i << ' ' << my_mmap.member(i)->name() << endl;
  cout << &edges[0] << ' ' << endl;
#endif

  // Determine if person is processed or not.
  if(my_inds[i].process)
  {
    long ml = my_mmap.mother_index(i);
    long fl = my_mmap.father_index(i);

    long mn[2]; mn[0] = my_inds[ml].node1; mn[1] = my_inds[ml].node2;
    long fn[2]; fn[0] = my_inds[fl].node1; fn[1] = my_inds[fl].node2;

    size_t mmei = my_mmap.mother_meiosis(i), fmei = my_mmap.father_meiosis(i);

//    for(size_t j = 0; j < 1 + (mmei >= SAGE::meiosis_map::meiosis_bits); ++j)
//      for(size_t k = 0; k < 1 + (fmei >= SAGE::meiosis_map::meiosis_bits); ++k)
    for(size_t j = 0; j < 2 ; ++j)
      for(size_t k = 0; k < 2; ++k)
      {
#if 0
        cout << "Setting " << my_mmap.member(i)->name() << " to " << j << "," << k << endl;
        cout << "       mother_index = " << ml << ", father_index = " << fl << endl;
        cout << "       mn = " << mn[0] << ", " << mn[1] << endl;
        cout << "       fn = " << fn[0] << ", " << fn[1] << endl;
        cout << "       mmei = " << mmei << ", fmei = " << fmei << endl;
#endif

        my_inds[i].node1 = mn[j];
        my_inds[i].node2 = fn[k];

        my_ivector[mmei] = j; my_ivector[fmei] = k;

#if 0
  cout << "name\tmeiosis\tprocess\tedge_no\tnode1\tnode2" << endl;

  for(int m = 0; m < my_mmap.get_subpedigree()->member_count(); ++m)
  {
    cout << my_mmap.member(m)->name() << '\t';
    if(!my_mmap.founder(m))
      cout << my_mmap.mother_meiosis(m) << ' ' << my_mmap.father_meiosis(m);

    cout << '\t' << my_inds[m].process
         << '\t' << my_inds[m].edge_number
         << '\t' << my_inds[m].node1
         << '\t' << my_inds[m].node2
         << endl;
  }

  cout << "4. my_ivector :" << endl;
  cout << "     nonf = " << my_ivector.get_nonfounders() << endl
       << "        f = " << my_ivector.get_founders() << endl;
  cout << "     mmei = " << mmei << ", j = " << j << endl
       << "     fmei = " << fmei << ", k = " << k << endl;
#endif

        if(my_inds[i].edge_number != -1)
          add_edge(i, mn[j], fn[k]);
        else
          build_ind(i+1);
      }

  }
  else  // if the person is not processed
  {
    build_ind(i+1);
  }
} 

void codominant_iv_generator::add_edge(int i, int j, int k) const
{
#if 0
  cout << "codominant_iv_generator::add_edge("
       << i << "," << j << "," << k << ")" << endl;
  cout << "edge_no\tid()\tstate\ta[0]/a[1]\tch\tn0+n1\tnode_ids.." << endl;
  cout << my_inds[i].edge_number << "\t"; edges[my_inds[i].edge_number].dump(cout);

  cout << "from\ti\tstate\tfixed\tch\ta[0]/a[1]/f_a\tconnection edges_name.." << endl;
  cout << "mat\t"; nodes[j].dump(cout, j);
  cout << "pat\t"; nodes[k].dump(cout, k);
#endif

  if(!nodes[j].add_edge(edges[my_inds[i].edge_number], 0)) return;

  if(!nodes[k].add_edge(edges[my_inds[i].edge_number], 1))
  {
    nodes[j].remove_edge();
    return;
  }

#if 0
  cout << "NODES: " << nodes.size() << endl;
  cout << "i\tstate\tfixed\tch\ta[0]/a[1]/f_a\tconnection edges_name.." << endl;
  for(int a = 0; a < nodes.size(); ++a)
    nodes[a].dump(cout, a);

  cout << "EDGES: " << edges.size() << endl;
  cout << "id()\tstate\ta[0]/a[1]\tch\tn0+n1\tnode_ids.." << endl;
  for(int a = 0; a < edges.size(); ++a)
    edges[a].dump(cout);
#endif

  build_ind(i+1);

  nodes[k].remove_edge();
  nodes[j].remove_edge();
}

void codominant_iv_generator::eval_ind(size_t i, double prob) const
{
#if 0
  if(i == 0)
  {
    cout << "codominant_iv_generator::eval_ind(" << i << ", " << prob << ")" << endl;
    cout << "NODES: " << nodes.size() << endl;
    cout << "i\tstate\tfixed\tch\ta[0]/a[1]/f_a\tconnection edges_name.." << endl;
    for(int a = 0; a < nodes.size(); ++a)
      nodes[a].dump(cout, a);
  
    cout << "EDGES: " << edges.size() << endl;
    cout << "id()\tstate\ta[0]/a[1]\tch\tn0+n1\tnode_ids.." << endl;
    for(int a = 0; a < edges.size(); ++a)
      edges[a].dump(cout);
  }
#endif

  if(i == nodes.size()) { evaluate(prob); return; }

  double d = nodes[i].probability(my_imodel, use_pf);

#if 0
  cout << "i = " << i << " d = " << d << endl;
#endif

  if(d == 0.0) return;

  if(d != -1)
  {
    eval_ind(i+1, prob * d);
    return;
  }


  if(nodes[i].fix(0))
  {
#if 0
    cout << "Fixing " << nodes[i].id() << " to first" << endl;

    cout << "NODES: " << nodes.size() << endl;
    for(int a = 0; a < nodes.size(); ++a)
      nodes[a].dump(cout, a);

    cout << "EDGES: " << edges.size() << endl;
    for(int a = 0; a < edges.size(); ++a)
      edges[a].dump(cout);
#endif

    d = nodes[i].probability(my_imodel, use_pf);
    if(d == 0.0) return;

    eval_ind(i+1, prob * d);

    nodes[i].unfix();

#if 0
    cout << "After fix" << endl;
    for(int a = 0; a < nodes.size(); ++a)
      nodes[a].dump(cout, a);

    for(int a = 0; a < edges.size(); ++a)
      edges[a].dump(cout);
#endif
  }
  if(nodes[i].fix(1))
  {
#if 0
    cout << "Fixing " << nodes[i].id() << " to second" << endl;

    for(int a = 0; a < nodes.size(); ++a)
      nodes[a].dump(cout, a);

    for(int a = 0; a < edges.size(); ++a)
      edges[a].dump(cout);
#endif

    d = nodes[i].probability(my_imodel, use_pf);
    if(d == 0.0) return;

    eval_ind(i+1, prob * d);

    nodes[i].unfix();

#if 0
    cout << "After fix" << endl;
    for(int a = 0; a < nodes.size(); ++a)
      nodes[a].dump(cout, a);

    for(int a = 0; a < edges.size(); ++a)
      edges[a].dump(cout);
#endif
  }
}

void codominant_iv_generator::evaluate(double prob) const
{
  list<size_t>::const_iterator i = my_hflags.begin();

  set_bit(i, prob);
}

void codominant_iv_generator::set_bit(list<size_t>::const_iterator& i, double prob)  const
{
//cout << "codominant_iv_generator::set_bit(" << prob << ")" << endl;

  boost::shared_ptr<iv_acceptor> iva2 = my_iva;

  // If we're at the end of the half bit flags, we're done.
  if(i == my_hflags.end())
  {
/*
    cout << "NODES: " << nodes.size() << endl;
    for(int a = 0; a < nodes.size(); ++a)
      nodes[a].dump(cout, a);

    cout << "EDGES: " << edges.size() << endl;
    for(int a = 0; a < edges.size(); ++a)
      edges[a].dump(cout);
*/

    iva2->accept(my_ivector.get_equivalence_class(), prob * my_penetrance);
    

    return;
  }

  // Get rid of repeats due to first founder meiosis.
  if(*i < my_mmap.founder_count())
  {
    my_ivector[*i] = 0;

    set_bit(++i, prob);

    --i;

    return;
  }


  int j = *i;
  ++i;

  my_ivector[j] = 0;

  set_bit(i, prob);

  my_ivector[j] = 1;

  set_bit(i, prob);

  --i;

}

}
