// ================
// Inline Functions
// ================

inline SymbolTable* segreg_data::parameters() const
{ return params; }

inline AttrVal segreg_data::parameter (string s) const
{ return params->query(s); }

inline const RPED::RefMultiPedigree* segreg_data::pedigrees()     const
{ return &pedigree_data; }

