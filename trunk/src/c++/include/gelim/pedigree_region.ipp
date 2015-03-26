inline
pedigree_region::pedigree_region
  (SAGE::cerrorstream& err, bool q)
    : my_errors(err),
      my_is_quiet(q),
      my_is_built(false),
      my_subpedigree()
{ }

inline 
pedigree_region::pedigree_region
  (const subpedigree&  sped,
   const region&       r,
   SAGE::cerrorstream& err,
   bool                quiet,
   bool                eliminate)
    : my_errors  (err),
      my_is_quiet(quiet),
      my_is_built(false)
{
  build(sped, r, eliminate);
}

inline 
pedigree_region::pedigree_region
  (const subpedigree&     sped,
   const pedigree_region& pr,
   const vector<uint>&    pr_ids,
   SAGE::cerrorstream&    err,
   bool                   quiet,
   bool                   eliminate)
    : my_errors  (err),
      my_is_quiet(quiet),
      my_is_built(false)
{
  build(sped, pr, pr_ids, eliminate);
}

inline
pedigree_region::pedigree_region
  (const pedigree_region& p)
    : my_errors        (p.my_errors),
      my_is_quiet      (p.my_is_quiet),
      my_is_built      (p.my_is_built),
      my_subpedigree(p.my_subpedigree),
      my_region     (p.my_region)
{
  my_markers = p.my_markers;

  my_model_consistencies = p.my_model_consistencies;
  my_model_informatives  = p.my_model_informatives;
}

inline
pedigree_region::~pedigree_region()
{ }

inline
pedigree_region& pedigree_region::operator=
  (const pedigree_region& p)
{
  if(&p == this) return *this;
  
  my_errors   = p.my_errors;
  my_is_quiet = p.my_is_quiet;
  my_is_built = p.my_is_built;

  my_subpedigree = p.my_subpedigree;
  my_region      = p.my_region;
  my_markers     = p.my_markers;

  my_model_consistencies = p.my_model_consistencies;
  my_model_informatives  = p.my_model_informatives;

  return *this;
}

inline
SAGE::cerrorstream& pedigree_region::get_errorstream()
{
  return my_errors;
}

inline
void pedigree_region::set_errorstream(SAGE::cerrorstream& err)
{
  my_errors = err;
}

inline
bool
pedigree_region::is_quiet() const
{
  return my_is_quiet;
}

inline
void
pedigree_region::set_quiet(bool q)
{
  my_is_quiet = q;
}

inline
const pedigree_region::subpedigree& pedigree_region::get_subpedigree() const
{
  return *my_subpedigree;
}

inline
const pedigree_region::region& pedigree_region::get_region() const
{
  return my_region;
}

inline bool pedigree_region::is_built() const { return my_is_built; }

inline
MLOCUS::inheritance_model&
    pedigree_region::operator[] (size_type t)
{
  return my_markers[t];
}

inline
const MLOCUS::inheritance_model&
    pedigree_region::operator[] (size_type t) const
{
  return my_markers[t];
}

inline bool
pedigree_region::model_consistent(size_type t) const
{
  return my_model_consistencies[t];
}
inline bool
pedigree_region::model_informative(size_type t) const
{
  return my_model_informatives[t];
}

inline
pedigree_region::size_type pedigree_region::inheritance_model_count() const
{
  return my_markers.size();
}

inline
void pedigree_region::reset()
{
  my_is_built = false;

  my_markers.clear();
}
