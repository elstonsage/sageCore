// ================
// Inline Functions
// ================

inline
descent_graph::descent_graph(const inheritance_vector& i, bool a)
    : mm(NULL), autosomal(a)
{
  move_to(i);
}

inline
descent_graph::descent_graph(equivalence_class e, const meiosis_map& m, bool a)
    : mm(NULL), autosomal(a)
{
  move_to(e, m);
}

inline descent_graph::~descent_graph() { }

inline void descent_graph::move_to(const inheritance_vector& iv)
{
  move_to(iv.get_equivalence_class(), iv.get_meiosis_map());
}

inline
descent_graph::allele_pair descent_graph::alleles(member_const_pointer id) const
{ return av[id->subindex()]; }

inline
descent_graph::allele_pair descent_graph::alleles(index i)  const
{ return av[i]; }

inline
int descent_graph::sharing(member_const_pointer i1, member_const_pointer i2) const
{ return sharing(i1->subindex(), i2->subindex()); }

inline
int descent_graph::sharing(index i1,  index i2)  const
{
  int p;

  p = (av[i1].first  == av[i2].first)  + 
      (av[i1].second == av[i2].second) +
      (av[i1].second == av[i2].first)  +
      (av[i1].first  == av[i2].second);


//
// Move autosomal testing into separate class.
/*
  if(autosomal)
  {
    p = (av[i1].first  == av[i2].first)  + 
        (av[i1].second == av[i2].second) +
        (av[i1].second == av[i2].first)  +
        (av[i1].first  == av[i2].second);
  }
  else
  {
    MP_Base::sex s1 = mm->get_mps()->get_sex(mm->id(i1)),
                 s2 = mm->get_mps()->get_sex(mm->id(i1));

    if( (s1 == MP_Base::Male || s1 == MP_Base::BMale) &&
        (s2 == MP_Base::Male || s2 == MP_Base::BMale)    )
    {
      return (av[i1].first  == av[i2].first);
    }
    else
    {
      p = (av[i1].first  == av[i2].first)  + 
          (av[i1].second == av[i2].second) +
          (av[i1].second == av[i2].first)  +
          (av[i1].first  == av[i2].second);
    }
  }
*/

  // If the pair is hyper-correlated (they are A/A & A/A, all IBD) we return
  // 2 alleles IBD.  This is not quite correct, but otherwise our model
  // fails.
  if(p > 4) return 2;

  return p;
}

inline
int descent_graph::mother_sharing(member_const_pointer i1, member_const_pointer i2) const
{ return mother_sharing(i1->subindex(), i2->subindex()); }

inline
int descent_graph::mother_sharing(index i1,  index i2)  const
{
  int p;

  p = (av[i1].first == av[i2].first);

  return p;
}

inline
int descent_graph::father_sharing(member_const_pointer i1, member_const_pointer i2) const
{ return father_sharing(i1->subindex(), i2->subindex()); }

inline
int descent_graph::father_sharing(index i1,  index i2)  const
{
  int p;

  p = (av[i1].second == av[i2].second);

  return p;
}
