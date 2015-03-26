namespace SAGE {
namespace RPED {

//============================================
//
//        PedData
//
//============================================

template <class GI, class FI, class SI, class PI, class MI>
inline PedData<GI, FI, SI, PI, MI>::PedData() 
{}

template <class GI, class FI, class SI, class PI, class MI>
inline
PedData<GI, FI, SI, PI, MI>::PedData(const MPED::multipedigree_base& m)
{
  update_storage(m);
}

template <class GI, class FI, class SI, class PI, class MI>
inline size_t PedData<GI, FI, SI, PI, MI>::pedigree_info_count    ()         const { return my_storage.ped_data.size();             }

template <class GI, class FI, class SI, class PI, class MI>
inline size_t PedData<GI, FI, SI, PI, MI>::subpedigree_info_count (size_t p) const { return my_storage.ped_data[p].sub_data.size(); }

template <class GI, class FI, class SI, class PI, class MI>
inline size_t PedData<GI, FI, SI, PI, MI>::family_info_count      (size_t p) const { return my_storage.ped_data[p].fam_data.size(); }

template <class GI, class FI, class SI, class PI, class MI>
inline size_t PedData<GI, FI, SI, PI, MI>::individual_info_count  (size_t p) const { return my_storage.ped_data[p].ind_data.size(); }

template <class GI, class FI, class SI, class PI, class MI>
inline size_t PedData<GI, FI, SI, PI, MI>::subpedigree_info_count (const MPED::pedigree_base& p) const { return my_storage.ped_data[p.index()].sub_data.size(); }

template <class GI, class FI, class SI, class PI, class MI>
inline size_t PedData<GI, FI, SI, PI, MI>::family_info_count      (const MPED::pedigree_base& p) const { return my_storage.ped_data[p.index()].fam_data.size(); }

template <class GI, class FI, class SI, class PI, class MI>
inline size_t PedData<GI, FI, SI, PI, MI>::individual_info_count  (const MPED::pedigree_base& p) const { return my_storage.ped_data[p.index()].ind_data.size(); }

template <class GI, class FI, class SI, class PI, class MI>
inline const typename PedData<GI, FI, SI, PI, MI>::mpinfo_type & PedData<GI, FI, SI, PI, MI>::multipedigree_info ()                   const { return my_storage.mped_info;               }

template <class GI, class FI, class SI, class PI, class MI>
inline const typename PedData<GI, FI, SI, PI, MI>::pedinfo_type& PedData<GI, FI, SI, PI, MI>::pedigree_info      (size_t p)           const { return my_storage.ped_data[p].ped_info;    }

template <class GI, class FI, class SI, class PI, class MI>
inline const typename PedData<GI, FI, SI, PI, MI>::subinfo_type& PedData<GI, FI, SI, PI, MI>::subpedigree_info   (size_t p, size_t s) const { return my_storage.ped_data[p].sub_data[s]; }

template <class GI, class FI, class SI, class PI, class MI>
inline const typename PedData<GI, FI, SI, PI, MI>::subinfo_type& PedData<GI, FI, SI, PI, MI>::family_info        (size_t p, size_t f) const { return my_storage.ped_data[p].fam_data[f]; }

template <class GI, class FI, class SI, class PI, class MI>
inline const typename PedData<GI, FI, SI, PI, MI>::indinfo_type& PedData<GI, FI, SI, PI, MI>::individual_info    (size_t p, size_t i) const { return my_storage.ped_data[p].ind_data[i]; }

template <class GI, class FI, class SI, class PI, class MI>
inline typename PedData<GI, FI, SI, PI, MI>::mpinfo_type & PedData<GI, FI, SI, PI, MI>::multipedigree_info ()                          { return my_storage.mped_info;               }

template <class GI, class FI, class SI, class PI, class MI>
inline typename PedData<GI, FI, SI, PI, MI>::pedinfo_type& PedData<GI, FI, SI, PI, MI>::pedigree_info      (size_t p)                  { return my_storage.ped_data[p].ped_info;    }

template <class GI, class FI, class SI, class PI, class MI>
inline typename PedData<GI, FI, SI, PI, MI>::subinfo_type& PedData<GI, FI, SI, PI, MI>::subpedigree_info   (size_t p, size_t s)        { return my_storage.ped_data[p].sub_data[s]; }

template <class GI, class FI, class SI, class PI, class MI>
inline typename PedData<GI, FI, SI, PI, MI>::subinfo_type& PedData<GI, FI, SI, PI, MI>::family_info        (size_t p, size_t f)        { return my_storage.ped_data[p].fam_data[f]; }

template <class GI, class FI, class SI, class PI, class MI>
inline typename PedData<GI, FI, SI, PI, MI>::indinfo_type& PedData<GI, FI, SI, PI, MI>::individual_info    (size_t p, size_t i)        { return my_storage.ped_data[p].ind_data[i]; }

template <class GI, class FI, class SI, class PI, class MI>
inline const typename PedData<GI, FI, SI, PI, MI>::pedinfo_type& PedData<GI, FI, SI, PI, MI>::pedigree_info    (const MPED::pedigree_base    & p)            const  { return my_storage.ped_data[p.index()].ped_info;                        }

template <class GI, class FI, class SI, class PI, class MI>
inline const typename PedData<GI, FI, SI, PI, MI>::subinfo_type& PedData<GI, FI, SI, PI, MI>::subpedigree_info (const MPED::pedigree_base    & p, size_t s)  const  { return my_storage.ped_data[p.index()].sub_data[s];                     }

template <class GI, class FI, class SI, class PI, class MI>
inline const typename PedData<GI, FI, SI, PI, MI>::subinfo_type& PedData<GI, FI, SI, PI, MI>::subpedigree_info (const MPED::subpedigree_base & s)            const  { return my_storage.ped_data[s.pedigree()->index()].sub_data[s.index()]; }

template <class GI, class FI, class SI, class PI, class MI>
inline const typename PedData<GI, FI, SI, PI, MI>::faminfo_type& PedData<GI, FI, SI, PI, MI>::family_info      (const MPED::pedigree_base    & p, size_t f)  const  { return my_storage.ped_data[p.index()].fam_data[f];                     }

template <class GI, class FI, class SI, class PI, class MI>
inline const typename PedData<GI, FI, SI, PI, MI>::faminfo_type& PedData<GI, FI, SI, PI, MI>::family_info      (const MPED::subpedigree_base & s, size_t f)  const  { return my_storage.ped_data[s.pedigree()->index()].fam_data[f];         }

template <class GI, class FI, class SI, class PI, class MI>
inline const typename PedData<GI, FI, SI, PI, MI>::faminfo_type& PedData<GI, FI, SI, PI, MI>::family_info      (const MPED::family_base      & f)            const  { return my_storage.ped_data[f.pedigree()->index()].fam_data[f.index()]; }

template <class GI, class FI, class SI, class PI, class MI>
inline const typename PedData<GI, FI, SI, PI, MI>::faminfo_type& PedData<GI, FI, SI, PI, MI>::family_info      (const MPED::member_base      & m)            const  { return my_storage.ped_data[m.pedigree()->index()].fam_data[m.index()]; }

template <class GI, class FI, class SI, class PI, class MI>
inline const typename PedData<GI, FI, SI, PI, MI>::indinfo_type& PedData<GI, FI, SI, PI, MI>::individual_info  (const MPED::pedigree_base    & p, size_t i)  const  { return my_storage.ped_data[p.index()].ind_data[i];                     }

template <class GI, class FI, class SI, class PI, class MI>
inline const typename PedData<GI, FI, SI, PI, MI>::indinfo_type& PedData<GI, FI, SI, PI, MI>::individual_info  (const MPED::subpedigree_base & s, size_t i)  const  { return my_storage.ped_data[s.pedigree()->index()].ind_data[i];         }

template <class GI, class FI, class SI, class PI, class MI>
inline const typename PedData<GI, FI, SI, PI, MI>::indinfo_type& PedData<GI, FI, SI, PI, MI>::individual_info  (const MPED::member_base      & m)            const  { return my_storage.ped_data[m.pedigree()->index()].ind_data[m.index()]; }

template <class GI, class FI, class SI, class PI, class MI>
inline typename PedData<GI, FI, SI, PI, MI>::pedinfo_type& PedData<GI, FI, SI, PI, MI>::pedigree_info    (const MPED::pedigree_base    & p)           { return my_storage.ped_data[p.index()].ped_info;                        }

template <class GI, class FI, class SI, class PI, class MI>
inline typename PedData<GI, FI, SI, PI, MI>::subinfo_type& PedData<GI, FI, SI, PI, MI>::subpedigree_info (const MPED::pedigree_base    & p, size_t s) { return my_storage.ped_data[p.index()].sub_data[s];                     }

template <class GI, class FI, class SI, class PI, class MI>
inline typename PedData<GI, FI, SI, PI, MI>::subinfo_type& PedData<GI, FI, SI, PI, MI>::subpedigree_info (const MPED::subpedigree_base & s)           { return my_storage.ped_data[s.pedigree()->index()].sub_data[s.index()]; }

template <class GI, class FI, class SI, class PI, class MI>
inline typename PedData<GI, FI, SI, PI, MI>::faminfo_type& PedData<GI, FI, SI, PI, MI>::family_info      (const MPED::pedigree_base    & p, size_t f) { return my_storage.ped_data[p.index()].fam_data[f];                     }

template <class GI, class FI, class SI, class PI, class MI>
inline typename PedData<GI, FI, SI, PI, MI>::faminfo_type& PedData<GI, FI, SI, PI, MI>::family_info      (const MPED::subpedigree_base & s, size_t f) { return my_storage.ped_data[s.pedigree()->index()].fam_data[f];         }

template <class GI, class FI, class SI, class PI, class MI>
inline typename PedData<GI, FI, SI, PI, MI>::faminfo_type& PedData<GI, FI, SI, PI, MI>::family_info      (const MPED::family_base      & f)           { return my_storage.ped_data[f.pedigree()->index()].fam_data[f.index()]; }

template <class GI, class FI, class SI, class PI, class MI>
inline typename PedData<GI, FI, SI, PI, MI>::faminfo_type& PedData<GI, FI, SI, PI, MI>::family_info      (const MPED::member_base      & m)           { return my_storage.ped_data[m.pedigree()->index()].fam_data[m.index()]; }

template <class GI, class FI, class SI, class PI, class MI>
inline typename PedData<GI, FI, SI, PI, MI>::indinfo_type& PedData<GI, FI, SI, PI, MI>::individual_info  (const MPED::pedigree_base    & p, size_t i) { return my_storage.ped_data[p.index()].ind_data[i];                     }

template <class GI, class FI, class SI, class PI, class MI>
inline typename PedData<GI, FI, SI, PI, MI>::indinfo_type& PedData<GI, FI, SI, PI, MI>::individual_info  (const MPED::subpedigree_base & s, size_t i) { return my_storage.ped_data[s.pedigree()->index()].ind_data[i];         }

template <class GI, class FI, class SI, class PI, class MI>
inline typename PedData<GI, FI, SI, PI, MI>::indinfo_type& PedData<GI, FI, SI, PI, MI>::individual_info  (const MPED::member_base      & m)           { return my_storage.ped_data[m.pedigree()->index()].ind_data[m.index()]; }

template <class GI, class FI, class SI, class PI, class MI>
inline typename PedData<GI, FI, SI, PI, MI>::mpinfo_type & PedData<GI, FI, SI, PI, MI>::info()                                 { return my_storage.mped_info;                                           }

template <class GI, class FI, class SI, class PI, class MI>
inline typename PedData<GI, FI, SI, PI, MI>::pedinfo_type& PedData<GI, FI, SI, PI, MI>::info(const MPED::pedigree_base    & p) { return my_storage.ped_data[p.index()].ped_info;                        }

template <class GI, class FI, class SI, class PI, class MI>
inline typename PedData<GI, FI, SI, PI, MI>::subinfo_type& PedData<GI, FI, SI, PI, MI>::info(const MPED::subpedigree_base & s) { return my_storage.ped_data[s.pedigree()->index()].sub_data[s.index()]; }

template <class GI, class FI, class SI, class PI, class MI>
inline typename PedData<GI, FI, SI, PI, MI>::faminfo_type& PedData<GI, FI, SI, PI, MI>::info(const MPED::family_base      & f) { return my_storage.ped_data[f.pedigree()->index()].fam_data[f.index()]; }

template <class GI, class FI, class SI, class PI, class MI>
inline typename PedData<GI, FI, SI, PI, MI>::indinfo_type& PedData<GI, FI, SI, PI, MI>::info(const MPED::member_base      & m) { return my_storage.ped_data[m.pedigree()->index()].ind_data[m.index()]; }

template <class GI, class FI, class SI, class PI, class MI>
inline const typename PedData<GI, FI, SI, PI, MI>::mpinfo_type & PedData<GI, FI, SI, PI, MI>::info()                                 const { return my_storage.mped_info;                                           }

template <class GI, class FI, class SI, class PI, class MI>
inline const typename PedData<GI, FI, SI, PI, MI>::pedinfo_type& PedData<GI, FI, SI, PI, MI>::info(const MPED::pedigree_base    & p) const { return my_storage.ped_data[p.index()].ped_info;                        }

template <class GI, class FI, class SI, class PI, class MI>
inline const typename PedData<GI, FI, SI, PI, MI>::subinfo_type& PedData<GI, FI, SI, PI, MI>::info(const MPED::subpedigree_base & s) const { return my_storage.ped_data[s.pedigree()->index()].sub_data[s.index()]; }

template <class GI, class FI, class SI, class PI, class MI>
inline const typename PedData<GI, FI, SI, PI, MI>::faminfo_type& PedData<GI, FI, SI, PI, MI>::info(const MPED::family_base      & f) const { return my_storage.ped_data[f.pedigree()->index()].fam_data[f.index()]; }

template <class GI, class FI, class SI, class PI, class MI>
inline const typename PedData<GI, FI, SI, PI, MI>::indinfo_type& PedData<GI, FI, SI, PI, MI>::info(const MPED::member_base      & m) const { return my_storage.ped_data[m.pedigree()->index()].ind_data[m.index()]; }

template <class GI, class FI, class SI, class PI, class MI>
inline typename PedData<GI, FI, SI, PI, MI>::pedinfo_type& PedData<GI, FI, SI, PI, MI>::operator[](const MPED::pedigree_base    & p) { return my_storage.ped_data[p.index()].ped_info;                        }

template <class GI, class FI, class SI, class PI, class MI>
inline typename PedData<GI, FI, SI, PI, MI>::subinfo_type& PedData<GI, FI, SI, PI, MI>::operator[](const MPED::subpedigree_base & s) { return my_storage.ped_data[s.pedigree()->index()].sub_data[s.index()]; }

template <class GI, class FI, class SI, class PI, class MI>
inline typename PedData<GI, FI, SI, PI, MI>::faminfo_type& PedData<GI, FI, SI, PI, MI>::operator[](const MPED::family_base      & f) { return my_storage.ped_data[f.pedigree()->index()].fam_data[f.index()]; }

template <class GI, class FI, class SI, class PI, class MI>
inline typename PedData<GI, FI, SI, PI, MI>::indinfo_type& PedData<GI, FI, SI, PI, MI>::operator[](const MPED::member_base      & m) { return my_storage.ped_data[m.pedigree()->index()].ind_data[m.index()]; }

template <class GI, class FI, class SI, class PI, class MI>
inline const typename PedData<GI, FI, SI, PI, MI>::pedinfo_type& PedData<GI, FI, SI, PI, MI>::operator[](const MPED::pedigree_base    & p) const { return my_storage.ped_data[p.index()].ped_info;                        }

template <class GI, class FI, class SI, class PI, class MI>
inline const typename PedData<GI, FI, SI, PI, MI>::subinfo_type& PedData<GI, FI, SI, PI, MI>::operator[](const MPED::subpedigree_base & s) const { return my_storage.ped_data[s.pedigree()->index()].sub_data[s.index()]; }

template <class GI, class FI, class SI, class PI, class MI>
inline const typename PedData<GI, FI, SI, PI, MI>::faminfo_type& PedData<GI, FI, SI, PI, MI>::operator[](const MPED::family_base      & f) const { return my_storage.ped_data[f.pedigree()->index()].fam_data[f.index()]; }

template <class GI, class FI, class SI, class PI, class MI>
inline const typename PedData<GI, FI, SI, PI, MI>::indinfo_type& PedData<GI, FI, SI, PI, MI>::operator[](const MPED::member_base      & m) const { return my_storage.ped_data[m.pedigree()->index()].ind_data[m.index()]; }

template <class GI, class FI, class SI, class PI, class MI>
inline
void PedData<GI, FI, SI, PI, MI>::update_storage(const MPED::multipedigree_base& m)
{
  size_t n = m.pedigree_count();

  my_storage.ped_data.resize(n);

  for(size_t i = 0; i < n; ++i)
    update_storage( m.pedigree_index(i) );
}

template <class GI, class FI, class SI, class PI, class MI>
inline
void PedData<GI, FI, SI, PI, MI>::update_storage(const MPED::pedigree_base& p)
{
  if( my_storage.ped_data.size() <= p.index())
    my_storage.ped_data.resize(p.index() + 1);

  my_storage.ped_data[p.index()].sub_data.resize( p.subpedigree_count() );
  my_storage.ped_data[p.index()].fam_data.resize( p.family_count()      );
  my_storage.ped_data[p.index()].ind_data.resize( p.member_count()      );
}

template <class GI, class FI, class SI, class PI, class MI>
inline
void PedData<GI, FI, SI, PI, MI>::update_storage(const MPED::subpedigree_base& s)
{
  update_storage( *s.pedigree() );
}

template <class GI, class FI, class SI, class PI, class MI>
inline
void PedData<GI, FI, SI, PI, MI>::update_storage(const MPED::family_base& f)
{
  update_storage( *f.pedigree() );
}

template <class GI, class FI, class SI, class PI, class MI>
inline
void PedData<GI, FI, SI, PI, MI>::update_storage(const MPED::member_base& m)
{
  update_storage( *m.pedigree() );
}

//============================================
//
//        RefTraitInfo
//
//============================================

inline
RefTraitInfo::RefTraitInfo()
{ 
  init(); 
  set_type(unknown_trait);
  set_usage(unknown_use);
}

inline
RefTraitInfo::RefTraitInfo(string trait_name, trait_t trait_type, trait_use usage)
{
  init();
  set_name(trait_name);
  set_alias_name(trait_name);
  set_type(trait_type);
  set_usage(usage);
}

inline const string & RefTraitInfo::name()       const { return my_name;       }
inline const string & RefTraitInfo::alias_name() const { return my_alias_name; }

inline RefTraitInfo::trait_t   RefTraitInfo::type()  const { return my_type; }
inline RefTraitInfo::trait_use RefTraitInfo::usage() const { return my_use;  }

inline double RefTraitInfo::threshold()               const { return my_threshold; }
inline double RefTraitInfo::numeric_missing_code()    const { return my_numeric_missing_code;    }
inline double RefTraitInfo::numeric_affected_code()   const { return my_numeric_affected_code;   }
inline double RefTraitInfo::numeric_unaffected_code() const { return my_numeric_unaffected_code; }

inline const string & RefTraitInfo::string_missing_code()    const { return my_string_missing_code;    }
inline const string & RefTraitInfo::string_affected_code()   const { return my_string_affected_code;   }
inline const string & RefTraitInfo::string_unaffected_code() const { return my_string_unaffected_code; }

inline void RefTraitInfo::set_name                   (const string    & n)   { my_name = n;                  }
inline void RefTraitInfo::set_alias_name             (const string    & n)   { my_alias_name = n;            }
inline void RefTraitInfo::set_type                   (      trait_t     typ) { my_type = typ;                }
inline void RefTraitInfo::set_usage                  (      trait_use   tru) { my_use = tru;                 }
inline void RefTraitInfo::set_threshold              (      double      th)  { my_threshold = th;            }
inline void RefTraitInfo::set_numeric_missing_code   (      double      nm)  { my_numeric_missing_code = nm; }
inline void RefTraitInfo::set_string_missing_code    (const string    & sm)  { my_string_missing_code = sm;  }
inline void RefTraitInfo::set_string_affected_code   (const string    & ac)  { my_string_affected_code    = ac; }
inline void RefTraitInfo::set_string_unaffected_code (const string    & uc)  { my_string_unaffected_code  = uc; }
inline void RefTraitInfo::set_numeric_affected_code  (      double      ac)  { my_numeric_affected_code   = ac; }
inline void RefTraitInfo::set_numeric_unaffected_code(      double      uc)  { my_numeric_unaffected_code = uc; }

inline 
void RefTraitInfo::init()
{
  double qNaN              = numeric_limits<double>::quiet_NaN();

  my_threshold             = my_numeric_missing_code    = qNaN;
  my_numeric_affected_code = my_numeric_unaffected_code = qNaN;
  my_string_affected_code  = my_string_unaffected_code  = "";
  my_string_missing_code   = "";
}

//==============================================
//
//      RefStringInfo
//
//==============================================

inline
RefStringInfo::RefStringInfo()
{ 
  init(); 
}

inline
RefStringInfo::RefStringInfo(string name)
{
  init();
  set_name(name);
}

inline const string & RefStringInfo::name                () const { return my_name; }
inline const string & RefStringInfo::string_missing_code () const { return my_string_missing_code; }

inline void RefStringInfo::set_name                (const string & n)  { my_name = n;                  }
inline void RefStringInfo::set_string_missing_code (const string & sm) { my_string_missing_code = sm;  }

inline
void RefStringInfo::init()
{
  my_string_missing_code   = "";
}

//===============================================
//
//        PhenotypeReaderInfo
//
//===============================================

inline
PhenotypeReaderInfo::PhenotypeReaderInfo()
{
  set_allele_adjustment(none);
  set_min_allele_freq(0.);
  set_max_allele_freq(1.);
  set_allele_delimiter('/');
  set_allele_missing("");
  set_covariate_moi("");
  set_covariate_allele("");
  set_allow_hemizygote(false);
}

inline void PhenotypeReaderInfo::set_allele_adjustment (allele_adj a)   { my_allele_adjustment   = a; }
inline void PhenotypeReaderInfo::set_min_allele_freq   (double m)       { my_minimum_allele_freq = m; }
inline void PhenotypeReaderInfo::set_max_allele_freq   (double m)       { my_maximum_allele_freq = m; }
inline void PhenotypeReaderInfo::set_allele_delimiter  (char s)         { my_allele_delimiter    = s; }
inline void PhenotypeReaderInfo::set_allele_missing    (const string s) { my_allele_missing      = s; }

inline void PhenotypeReaderInfo::set_covariate_moi     (const string c) { my_cov_moi = c; }
inline void PhenotypeReaderInfo::set_covariate_allele  (const string s) { my_cov_allele = s; }
inline void PhenotypeReaderInfo::set_allow_hemizygote  (bool b)         { my_allow_hemizygote = b; }

inline PhenotypeReaderInfo::allele_adj PhenotypeReaderInfo::get_allele_adjustment() const { return my_allele_adjustment; }

inline double PhenotypeReaderInfo::get_min_allele_freq() const { return my_minimum_allele_freq; }
inline double PhenotypeReaderInfo::get_max_allele_freq() const { return my_maximum_allele_freq; }

inline char         PhenotypeReaderInfo::get_allele_delimiter() const { return my_allele_delimiter; }
inline const string PhenotypeReaderInfo::get_allele_missing  () const { return my_allele_missing; }

inline const string PhenotypeReaderInfo::get_covariate_moi()    const { return my_cov_moi; }
inline const string PhenotypeReaderInfo::get_covariate_allele() const { return my_cov_allele; }
inline bool         PhenotypeReaderInfo::get_allow_hemizygote() const { return my_allow_hemizygote; }

//=============================================
//
//           RefMPedInfo
//
//=============================================

inline size_t RefMPedInfo::trait_count  () const  { return my_trait_info.size();    }
inline size_t RefMPedInfo::string_count () const  { return my_string_info.size();   }
inline size_t RefMPedInfo::marker_count () const  { return my_marker_info.size();   }

inline       MLOCUS::inheritance_model_map & RefMPedInfo::markers()       { return my_marker_info; }
inline const MLOCUS::inheritance_model_map & RefMPedInfo::markers() const { return my_marker_info; }

inline       RefTraitInfo & RefMPedInfo::trait_info(size_t t)       { return my_trait_info[t]; }
inline const RefTraitInfo & RefMPedInfo::trait_info(size_t t) const { return my_trait_info[t]; }

inline       RefStringInfo & RefMPedInfo::string_info(size_t s)       { return my_string_info[s]; }
inline const RefStringInfo & RefMPedInfo::string_info(size_t s) const { return my_string_info[s]; }

inline       RefMarkerInfo & RefMPedInfo::marker_info(size_t m)        { return my_marker_info[m]; }
inline const RefMarkerInfo & RefMPedInfo::marker_info(size_t m) const  { return my_marker_info[m]; }
inline       RefMarkerInfo & RefMPedInfo::marker_info(const string& m) { return my_marker_info[m]; }

inline
size_t RefMPedInfo::string_find(const std::string &name) const
{
  string sn = toUpper(name);

  for(size_t s = 0; s < string_count(); ++s)
    if( toUpper(string_info(s).name()) == sn )
      return s;

  return (size_t)-1;
}

inline
size_t RefMPedInfo::marker_find(const std::string &name) const
{
  return my_marker_info.index(name);
}

inline
bool RefMPedInfo::trait_exists(const string& trait_name) const
{
  return trait_find(trait_name) != (size_t)(-1);
}

inline
bool RefMPedInfo::marker_exists(const string& marker_name) const
{
  return marker_find(marker_name) != (size_t)(-1);
}

inline 
RefTraitInfo::trait_t RefMPedInfo::get_trait_type(const string& trait_name) const
{
  size_t  trait_index = trait_find(trait_name);

  if(trait_index != (size_t)(-1))
  {
    return trait_info(trait_index).type();
  }
  else
  {
    return RefTraitInfo::invalid_trait;
  }
}

inline
size_t RefMPedInfo::add_continuous_trait(const string &trait_name, 
                 RefTraitInfo::trait_use use)
{
  return add_trait(trait_name, RefTraitInfo::continuous_trait, use);
}

inline
size_t RefMPedInfo::add_binary_trait(const string &trait_name, RefTraitInfo::trait_use use)
{
  return add_trait(trait_name, RefTraitInfo::binary_trait, use);
}

inline
size_t RefMPedInfo::add_continuous_covariate(const string &covariate_name)
{
  return add_trait(covariate_name, RefTraitInfo::continuous_trait, RefTraitInfo::trait_covariate);
}

inline
size_t RefMPedInfo::add_binary_covariate(const string &covariate_name)
{
  return add_trait(covariate_name, RefTraitInfo::binary_trait, RefTraitInfo::trait_covariate);
}

inline
size_t RefMPedInfo::add_string_field(const string &name)
{
  size_t s = string_find(name);
  if( s < string_count() )
    return s;

  my_string_info.push_back( RefStringInfo(name) );
  return string_count();
}

inline
size_t RefMPedInfo::add_marker(const string &marker_name)
{
  size_t m = marker_find(marker_name);

  if( m < marker_count() )
    return m;

  RefMarkerInfo& info = my_marker_info[toUpper(marker_name)] = RefMarkerInfo();
  info.set_name(marker_name);
  info.gmodel().set_name(marker_name);
  info.phmodel().set_name(marker_name);        
  return marker_count()-1;
}

inline
size_t RefMPedInfo::add_marker(const string &marker_name, const MLOCUS::inheritance_model& model)
{
  size_t m = marker_find(marker_name);

  if( m < marker_count() )
    return m;

  my_marker_info[toUpper(marker_name)] = model;
  return marker_count()-1;
}
  
inline
void RefMPedInfo::remove_last_trait() 
{
  my_trait_info.pop_back();
}

inline const string & RefMPedInfo::individual_missing_code () const { return my_individual_missing_code; }
inline const string & RefMPedInfo::sex_code_male           () const { return my_sex_code_male;    }
inline const string & RefMPedInfo::sex_code_female         () const { return my_sex_code_female;  }
inline const string & RefMPedInfo::sex_code_unknown        () const { return my_sex_code_unknown; }

inline void RefMPedInfo::set_individual_missing_code (const string &imc) { my_individual_missing_code = strip_ws(imc); }
inline void RefMPedInfo::set_sex_code_male           (const string &c)   { my_sex_code_male = c; }
inline void RefMPedInfo::set_sex_code_female         (const string &c)   { my_sex_code_female = c; }
inline void RefMPedInfo::set_sex_code_unknown        (const string &c)   { my_sex_code_unknown = c; }

inline       PhenotypeReaderInfo & RefMPedInfo::get_pheno_reader_info()       { return my_pheno_reader_info; }
inline const PhenotypeReaderInfo & RefMPedInfo::get_pheno_reader_info() const { return my_pheno_reader_info; }
  
//======================================================
//
//          RefPedInfo
//
//======================================================

inline
RefPedInfo::RefPedInfo() : my_trait_count(0) 
{ }

inline size_t RefPedInfo::trait_count  ()    const { return my_trait_count;   }
inline size_t RefPedInfo::string_count ()    const { return my_string_count;  }
inline size_t RefPedInfo::marker_count ()    const { return my_marker_count;  }
inline size_t RefPedInfo::member_count ()    const { return my_member_count;  }

inline
double RefPedInfo::trait(size_t i, size_t t) const 
{ 
  if(i >= member_count() || t >= trait_count())
     return std::numeric_limits<double>::quiet_NaN();
  return my_traits[t][i]; 
}

inline
bool RefPedInfo::trait_missing(size_t i, size_t t) const 
{ 
  if(i >= member_count() || t >= trait_count())
     return true;
  return SAGE::isnan(trait(i,t)); 
}

inline
bool RefPedInfo::set_trait(size_t i, size_t t, double d)
{ 
  if(i >= member_count() || t >= trait_count())
     return false;
  my_traits[t][i] = d; 
  return true;
}

inline
string RefPedInfo::get_string(size_t i, size_t s) const 
{ 
  if(i >= member_count() || s >= string_count())
     return "";
  return my_strings[s][i]; 
}

inline
bool RefPedInfo::set_string(size_t i, size_t s, const string& val)
{ 
  if(i >= member_count() || s >= string_count())
     return false;
  my_strings[s][i] = val; 
  return true;
}

inline
uint RefPedInfo::phenotype(size_t i, size_t m) const 
{ 
  if(i >= member_count() || m >= marker_count())
     return MLOCUS::NPOS;
  return my_markers[m][i]; 
}

inline
bool RefPedInfo::phenotype_missing(size_t i, size_t m, const RefMarkerInfo& mi) const 
{ 
//    return my_markers[m][i] == mi.get_missing_phenotype_id(); 
  return    my_markers[m][i] == mi.get_missing_phenotype_id()
         || my_markers[m][i] == MLOCUS::NPOS;
}

inline
bool RefPedInfo::set_phenotype(size_t i, size_t m, uint p)
{ 
  if(i >= member_count() || m >= marker_count())
     return false;
  my_markers[m][i] = p; 
  return true;
}

inline
void RefPedInfo::resize_traits(size_t traits)
{
  my_trait_count = traits;
  my_traits.resize(traits);

  double qnan = numeric_limits<double>::quiet_NaN();
  for(size_t i=0; i < traits; ++i)
  {
// WHY???   my_traits[i].resize(0);
    my_traits[i].resize(member_count(), qnan );
  }
}

inline
void RefPedInfo::resize_strings(size_t s)
{
  my_string_count = s;
  my_strings.resize(s);

  string null("");
  for(size_t i=0; i < s; ++i)
    my_strings[i].resize(member_count(), null);
}

inline
void RefPedInfo::resize_markers(size_t markers, const RefMPedInfo& mped_info)
{
  my_marker_count = markers;
  my_markers.resize(markers);

  for(size_t i=0; i < markers; ++i)
  {
    my_markers[i].resize(member_count(), mped_info.marker_info(i).get_missing_phenotype_id());
  }
}

inline
void RefPedInfo::swap_members(size_t m1, size_t m2)
{
  if(m1 >= member_count() || m2 >=  member_count() )
    return;
      
  for(size_t t=0; t < trait_count(); ++t)
    std::swap( my_traits[t][m1], my_traits[t][m2] );

  for(size_t s=0; s < string_count(); ++s)
    std::swap( my_strings[s][m1], my_strings[s][m2] );

  for(size_t m=0; m < marker_count(); ++m)
    std::swap( my_markers[m][m1], my_markers[m][m2] );
}
      
inline void RefPedInfo::build(MPED::pedigree_base &pedbase)
{
  my_member_count = pedbase.member_count();
}

} // End namespace RPED
} // End namespace SAGE
