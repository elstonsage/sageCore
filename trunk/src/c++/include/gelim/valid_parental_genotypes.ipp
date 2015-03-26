inline
  valid_parental_genotypes::valid_parental_genotypes()
  : my_genotype_pairs(),
    my_mother_valid(),
    my_father_valid(),
    my_mother_invalid(),
    my_father_invalid()
{ }

inline
  valid_parental_genotypes::valid_parental_genotypes
  (const family_type& fam,
   const imodel&      model,
   bool               phased)
{
  //lint --e{534}
  generate_valid_parental_genotypes(fam, model, phased);
}

inline
  valid_parental_genotypes::valid_parental_genotypes
  (const valid_parental_genotypes& v)
  : my_genotype_pairs(v.my_genotype_pairs),
    my_mother_valid(v.my_mother_valid),
    my_father_valid(v.my_father_valid),
    my_mother_invalid(v.my_mother_invalid),
    my_father_invalid(v.my_father_invalid)
{ }

inline
valid_parental_genotypes&
  valid_parental_genotypes::operator=(const valid_parental_genotypes& v)
{
  if(&v == this) return *this;

  my_genotype_pairs = v.my_genotype_pairs;

  my_mother_valid   = v.my_mother_valid;
  my_father_valid   = v.my_father_valid;
  my_mother_invalid = v.my_mother_invalid;
  my_father_invalid = v.my_father_invalid;

  return *this;
}


inline
valid_parental_genotypes::genotype_pair_iterator
  valid_parental_genotypes::genotype_pair_begin() const
{
  return my_genotype_pairs.begin();
}

inline
valid_parental_genotypes::genotype_pair_iterator
  valid_parental_genotypes::genotype_pair_end()   const
{
  return my_genotype_pairs.end();
}

inline
valid_parental_genotypes::parental_genotype_pair
  valid_parental_genotypes::genotype_pair(size_t i) const
{
  return my_genotype_pairs[i];
}
  
inline size_t valid_parental_genotypes::genotype_pair_count() const
{
  return my_genotype_pairs.size();
}


inline
valid_parental_genotypes::genotype_iterator
    valid_parental_genotypes::mother_valid_genotype_begin() const
{
  return my_mother_valid.begin();
}

inline
valid_parental_genotypes::genotype_iterator
    valid_parental_genotypes::mother_valid_genotype_end()   const
{
  return my_mother_valid.end();
}

  
inline
uint valid_parental_genotypes::mother_valid_genotype(size_t i) const
{
  return my_mother_valid[i];
}
  
inline size_t valid_parental_genotypes::mother_valid_genotype_count() const
{
  return my_mother_valid.size();
}


inline
valid_parental_genotypes::genotype_iterator
    valid_parental_genotypes::father_valid_genotype_begin() const
{
  return my_father_valid.begin();
}

inline
valid_parental_genotypes::genotype_iterator
    valid_parental_genotypes::father_valid_genotype_end()   const
{
  return my_father_valid.end();
}

  
inline
uint valid_parental_genotypes::father_valid_genotype(size_t i) const
{
  return my_father_valid[i];
}
  
inline size_t valid_parental_genotypes::father_valid_genotype_count() const
{
  return my_father_valid.size();
}


inline
valid_parental_genotypes::genotype_iterator
    valid_parental_genotypes::mother_invalid_genotype_begin() const
{
  return my_mother_invalid.begin();
}

inline
valid_parental_genotypes::genotype_iterator
    valid_parental_genotypes::mother_invalid_genotype_end()   const
{
  return my_mother_invalid.end();
}

  
inline
uint valid_parental_genotypes::mother_invalid_genotype(size_t i) const
{
  return my_mother_invalid[i];
}
  
inline size_t valid_parental_genotypes::mother_invalid_genotype_count() const
{
  return my_mother_invalid.size();
}


inline
valid_parental_genotypes::genotype_iterator
    valid_parental_genotypes::father_invalid_genotype_begin() const
{
  return my_father_invalid.begin();
}

inline
valid_parental_genotypes::genotype_iterator
    valid_parental_genotypes::father_invalid_genotype_end()   const
{
  return my_father_invalid.end();
}

  
inline
uint valid_parental_genotypes::father_invalid_genotype(size_t i) const
{
  return my_father_invalid[i];
}
  
inline size_t valid_parental_genotypes::father_invalid_genotype_count() const
{
  return my_father_invalid.size();
}

inline void valid_parental_genotypes::clean_internals()
{
  my_genotype_pairs = parental_genotype_pair_vector();

  my_mother_valid   = parental_genotype_vector();
  my_father_valid   = parental_genotype_vector();
  my_mother_invalid = parental_genotype_vector();
  my_father_invalid = parental_genotype_vector();
}
