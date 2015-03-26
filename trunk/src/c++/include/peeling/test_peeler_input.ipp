// ================
// Inline Functions
// ================

inline SymbolTable* test_peeler_data::parameters() const
{ return params; }

inline AttrVal test_peeler_data::parameter (string s) const
{ return params->query(s); }

inline const RPED::RefMultiPedigree* test_peeler_data::pedigrees()     const
{ return &pedigree_data; }

