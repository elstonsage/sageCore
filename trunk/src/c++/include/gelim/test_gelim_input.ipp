// ================
// Inline Functions
// ================

inline SymbolTable* gelim_data::parameters() const
{ return params; }

inline AttrVal gelim_data::parameter (string s) const
{ return params->query(s); }

inline const RPED::RefMultiPedigree* gelim_data::pedigrees()     const
{ return &pedigree_data; }

