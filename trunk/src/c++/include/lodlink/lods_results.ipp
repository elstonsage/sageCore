//============================================================================
// File:      lods_results.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   2/12/3 created        -djb
//                                                                          
// Notes:     Inline implementation of lod results classes.
//                                                                          
// Copyright (c) 2003 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


//============================================================================
// IMPLEMENTATION:  lods_result
//============================================================================
//
inline
lods_result::~lods_result()
{}

inline double
lods_result::get_subpedigree_lod_score(const multipedigree_lod_score& mpls, 
                                          const string& pedigree, const string& member)
{
  vector<pedigree_lod_score>::const_iterator  ped_iter;
  ped_iter = find(mpls.ped_lod_scores.begin(), mpls.ped_lod_scores.end(),
                  pedigree_lod_score(pedigree));
  
  if(ped_iter != mpls.ped_lod_scores.end())
  {
    vector<subpedigree_lod_score>::const_iterator  subped_iter;  
    subped_iter = find(ped_iter->sub_lod_scores.begin(),
                       ped_iter->sub_lod_scores.end(),
                       subpedigree_lod_score(member));
                       
    if(subped_iter != ped_iter->sub_lod_scores.end())
    {
      return  subped_iter->lod_score;
    }
  }
  
  return  QNAN;
}


//============================================================================
// IMPLEMENTATION:  subpedigree_lod_score
//============================================================================
//
inline
subpedigree_lod_score::subpedigree_lod_score(const string& mn, const double& lod)
      : member_name(mn), lod_score(lod)
{}

inline bool
subpedigree_lod_score::operator==(const subpedigree_lod_score& right) const
{
  return  member_name == right.member_name;
}

//============================================================================
// IMPLEMENTATION:  pedigree_lod_score
//============================================================================
//
inline
pedigree_lod_score::pedigree_lod_score(const string& pn, const double& lod)
      : pedigree_name(pn), lod_score(lod)
{}

inline bool
pedigree_lod_score::operator==(const pedigree_lod_score& right) const
{
  return  pedigree_name == right.pedigree_name;
}

inline void
pedigree_lod_score::write(ostream& out)
{
  out << "\n\n"
      << pedigree_name << "   " << endl;
  
  for(size_t i = 0; i < sub_lod_scores.size(); ++i)
  {
    out << "   " << sub_lod_scores[i].member_name << "   " 
                 << sub_lod_scores[i].lod_score << endl;
  }
}


//============================================================================
// IMPLEMENTATION:  multipedigree_lod_score
//============================================================================
//
inline
multipedigree_lod_score::multipedigree_lod_score(const double& lod)
      : lod_score(lod)
{}

inline void
multipedigree_lod_score::write(ostream& out)
{
  
}


//============================================================================
// IMPLEMENTATION:  non_ss_lods_result
//============================================================================
//
inline
non_ss_lods_result::non_ss_lods_result()
{}

inline void
non_ss_lods_result::write_summary(ostream& out) const
{
  ios::fmtflags old_flags = out.flags();

  out << left << setw(_LOCUS.lw()) << marker 
      << setfill(' ') << setw(_LOCUS.sw()) << "";
      
  out << right << resetiosflags(ios::floatfield)  << setprecision(PRC2);
  for(size_t i = 0; i < lod_scores.size(); ++i)
  {
    out << setw(RECOM_SZ);
    write_double(out, lod_scores[i].second.lod_score);
    out << setw(SPACE_SZ) << ""; 
  }
  
  out << endl;
  
  out.flags(old_flags);
}

inline void
non_ss_lods_result::write_detail(ostream& out) const
{}

inline void
non_ss_lods_result::write_family_detail(ostream& out, const string& pedigree,
                                                      const string& member   ) const
{
  ios::fmtflags old_flags = out.flags();

  out << left << setw(_LOCUS.lw()) << marker 
      << setfill(' ') << setw(_LOCUS.sw()) << "";
      
  out << right << resetiosflags(ios::floatfield)  << setprecision(PRC2);
  for(size_t i = 0; i < lod_scores.size(); ++i)
  {
    out << setw(RECOM_SZ);
    write_double(out, get_subpedigree_lod_score(lod_scores[i].second, pedigree, member));
    out << setw(SPACE_SZ) << ""; 
  }
  
  out << endl;
  
  out.flags(old_flags);
}

inline void
non_ss_lods_result::write_vc_matrix(ostream& out) const
{

}


//============================================================================
// IMPLEMENTATION:  ss_lods_result
//============================================================================
//
inline
ss_lods_result::ss_lods_result()
{}

inline void
ss_lods_result::write_summary(ostream& out) const
{
  ios::fmtflags old_flags = out.flags();

  out << left << setw(_LOCUS.lw()) << marker 
      << setfill(' ') << setw(_LOCUS.sw()) << "";
      
  out << right << resetiosflags(ios::floatfield)  << setprecision(PRC2);
  for(size_t i = 0; i < lod_scores.size(); ++i)
  {
    out << setw(RECOM_SZ);
    write_double(out, lod_scores[i].second.lod_score);
    out << setw(SPACE_SZ) << ""; 
  }
  
  out << endl;
  
  out.flags(old_flags);    
}

inline void
ss_lods_result::write_detail(ostream& out) const
{}

inline void
ss_lods_result::write_family_detail(ostream& out, const string& pedigree,
                                                      const string& member   ) const
{
  ios::fmtflags old_flags = out.flags();

  out << left << setw(_LOCUS.lw()) << marker 
      << setfill(' ') << setw(_LOCUS.sw()) << "";
      
  out << right << resetiosflags(ios::floatfield)  << setprecision(PRC2);
  for(size_t i = 0; i < lod_scores.size(); ++i)
  {
    out << setw(RECOM_SZ);
    write_double(out, get_subpedigree_lod_score(lod_scores[i].second, pedigree, member));
    out << setw(SPACE_SZ) << ""; 
  }
  
  out << endl;
  
  out.flags(old_flags);
}

inline void
ss_lods_result::write_vc_matrix(ostream& out) const
{

}
