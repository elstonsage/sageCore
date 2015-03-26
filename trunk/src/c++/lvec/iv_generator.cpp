#include <algorithm>
#include "lvec/inheritance_vector.h"
#include "lvec/iv_generator.h"

namespace SAGE
{

// The following structure (parent_vector and helpers) is used to keep track
// of the child genotypes of each parent as we procede through the pedigree.

struct parent_pair
{
  typedef meiosis_map::member_const_pointer ind_id;

  ind_id          moth, fath;

  SAGE::MLOCUS::child_genotype_set cg;
  bool               b;

#ifdef HFLAGS
  int             valid;  //< Which child genotypes are valid? May change XXX
#endif
  
  bool operator()(ind_id m, ind_id f) { return m == moth && f == fath; }

  bool operator<(const parent_pair& rhs) const
  {
    if(moth < rhs.moth) return true;
    if(moth == rhs.moth && fath < rhs.fath) return true;

    return false;
  }

  parent_pair(ind_id m = NULL, ind_id f = NULL) : moth(m), fath(f) { }
};

class parent_vector : public vector<parent_pair>
{
public:

  typedef vector<parent_pair>                    vect;
  typedef meiosis_map::subpedigree_const_pointer sped_id;
  typedef meiosis_map::member_const_pointer      ind_id;

  parent_vector() : spid(NULL) { }
  parent_vector(const SAGE::meiosis_map& m, const SAGE::MLOCUS::inheritance_model& i) : spid(NULL)
  {
    set_family(m, i);
  }

  ~parent_vector() { }

  void set_family(const SAGE::meiosis_map& m, const SAGE::MLOCUS::inheritance_model& imodel)
  {
    spid = m.get_subpedigree();

    if(!spid) return;
      
    resize(m.get_subpedigree()->family_count());
    
    size_t parents = 0;

    // Iterate through each member of the subpedigree.  For each non-founder, check
    // to see if its family has been shown yet.  If not, add it to the vector.  Note
    // that this can probably be done using the Family Iterator much more efficiently.
    for(size_t i = 0; i < m.get_subpedigree()->member_count(); ++i)
    {
      if( !m.founder(i) )
      {
        bool b = false;

        // Determine if the family has already been shown.
        for(size_t j = 0; !b && j < parents; ++j)
          if((*this)[j].moth == m.mother(i) &&
             (*this)[j].fath == m.father(i)) b = true;

        // If we've seen this family, we don't have to do anything.
        if(b) continue;
        
        (*this)[parents].moth = m.mother(i);
        (*this)[parents].fath = m.father(i);
        (*this)[parents].b    = false;

        ++parents;
      }
    }

    std::sort(begin(), end());
  }

  int find(ind_id moth, ind_id fath)
  {
    iterator i = std::lower_bound(begin(), end(), parent_pair(moth, fath));

    int x = i - begin();

    return x;
  }
    
protected:

  sped_id spid;
};

bool iv_generator::build
   (const SAGE::meiosis_map& _m, const SAGE::MLOCUS::inheritance_model& imod,
    bool use_pop_freq) const
{
  use_pf = use_pop_freq;

  mm = _m;

  if(   !mm.get_subpedigree()->member_count()
      || mm.get_subpedigree()->member_count() != imod.phenotype_count() )
    return false;

  subpedigree_const_pointer pid = mm.get_subpedigree();

  if(!pid) return false;

  genotypes.resize(mm.get_subpedigree()->member_count());

  parent_vector pv(mm, imod);

  iv(mm);
  
  build_ind(0, 1, imod, pv);

  return true;
}

void iv_generator::build_ind(size_t i, double prob, const SAGE::MLOCUS::inheritance_model& imod,
                             parent_vector& pv) const
{
  if( i == mm.get_subpedigree()->member_count() )
  {
    evaluate(prob);
    return;
  }

  // Determine if founder or not.
  if( !mm.founder(i) )
  {
    // Find the parents and create their genotypes if not already done
    int par = pv.find(mm.mother(i), mm.father(i));

    bool b = pv[par].b;
    
    if( !b )
    {
      size_t m = mm.mother_index(i);
      size_t f = mm.father_index(i);

      pv[par].b = true;

      pv[par].cg = 
          SAGE::MLOCUS::child_genotype_set(genotypes[m], genotypes[f]);

#ifdef HFLAGS
      // Determine which genotypes we need to check.  Genotypes that
      // are the same due to homozygosity we can avoid evaluating
      // twice.

      pv[par].valid = 0;
      if(genotypes[m].homozygous()) pv[par].valid += 2;
      if(genotypes[f].homozygous()) pv[par].valid += 1;
#endif
    }
    

#ifdef HFLAGS
    // If the parents are homozygous, then all children from them are
    // unknown.  We keep a list and only generate for one of the two
    // options

    if(pv[par].valid & 1) hflags.push_back(mm.father_meiosis(i));
    if(pv[par].valid & 2) hflags.push_back(mm.mother_meiosis(i));
#endif

    int count;

    // For each possible genotype from those parents, see if plausible for
    // current individual
    for(count = 0; count < 4; ++count)
    {
#ifdef HFLAGS
      // If the genotype isn't valid due to homozygosity, skip it.
      if(pv[par].valid & count) continue;
#endif

      // If the child's genotype is valid

      genotypes[i] = pv[par].cg[count];

      double penetrance = imod.phased_penetrance(i+1, genotypes[i].get_id());

      if(penetrance > 0)
      {
        // Set the inheritance vector bits
        iv[mm.mother_meiosis(i)] = count & 2;
        iv[mm.father_meiosis(i)] = count & 1;
        
        // Go to the next individual
        build_ind(i+1, prob * penetrance, imod, pv);
      }
    }

    // Remove the homozygous parents
#ifdef HFLAGS
    if(pv[par].valid & 1) hflags.pop_back();
    if(pv[par].valid & 2) hflags.pop_back();
#endif

    // reset the parents generated flag.
    pv[par].b = b;
  }
  else  // if the person is a founder
  {
    SAGE::MLOCUS::inheritance_model::phased_penetrance_iterator p =
        imod.phased_penetrance_begin(i+1);

    SAGE::MLOCUS::inheritance_model::phased_penetrance_iterator pend =
        imod.phased_penetrance_end(i+1);

    for( ; p != pend; ++p)
    {
      genotypes[i] = p.phased_geno();

      double frequency = 1.0;

      if(use_pf) frequency = genotypes[i].frequency();

      double penetrance = *p;

//      cout << "   " << genotypes[i].name() << ' ' << frequency << ' ' << penetrance << endl;

      build_ind(i+1, prob * frequency * penetrance, imod, pv);
    }
  }

} 

void iv_generator::evaluate(double prob) const
{
#ifdef HFLAGS
  list<size_t>::iterator i = hflags.begin();

  set_bit(i, prob);
#else
  iva->accept(iv, prob);
#endif
}

#ifdef HFLAGS
void iv_generator::set_bit(list<size_t>::iterator& i, double prob)  const
{
  // This is required to make the object non-const.
  boost::shared_ptr<iv_acceptor> iva2 = iva;

  if(i == hflags.end())
  {
/*
    for(size_t i = 0; i < genotypes.size(); ++i)
    {
      cout << genotypes[i].name() << endl;
    }
*/

    iva2->accept(iv.get_equivalence_class(), prob);
    return;
  }

  // Get rid of repeats due to first founder meiosis.
  if(*i < mm.founder_meiosis_count())
  {
    iv[*i] = 0;

    set_bit(++i, prob*2);

    --i;

    return;
  }

  int j = *i;
  ++i;

  iv[j] = 0;

  set_bit(i, prob);

  iv[j] = 1;

  set_bit(i, prob);

  --i;

}

#endif

}
