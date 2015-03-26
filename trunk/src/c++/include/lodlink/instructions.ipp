//============================================================================
// File:      instructions.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   9/5/2 created        -djb
//                                                                          
// Notes:     Inline implementation of struct, instructions.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


inline bool
operator<(const theta_pair& left, const theta_pair& right)
{
  if(left.male_theta == right.male_theta)
  {
    return  left.female_theta < right.female_theta;
  }
  else
  {
    return  left.male_theta < right.male_theta;
  }
}

inline bool
operator==(const theta_pair& left, const theta_pair& right)
{
  return left.male_theta   == right.male_theta   &&
         left.female_theta == right.female_theta   ;
}

inline ostream& 
operator<<(ostream& out, const theta_pair& p)
{
  out << "Male theta:    " << p.male_theta    << "\n"
      << "Female theta:  " << p.female_theta  << "\n"
      << endl;
      
  return out;
}

inline ostream&
operator<<(ostream& out, const instructions& m)
{
  out << std::boolalpha
      << "Valid:                            " << m.valid                << "\n" 
      << "Output root name:                 " << m.file_name_root       << "\n"
      << "Linkage option:                   " << m.link_option()        << "\n"
      << "Title:                            " << m.title                << "\n"
      << "Trait:                            " << (m.trait.empty() ? "none" : m.trait) << "\n"
      << "Linkage test:                     " << m.linkage_test         << "\n"
      << "Linkage test sex-specific:        " << m.linkage_sex_specific << "\n"
      << "Linkage test assumes homogeneity: " << m.linkage_homog        << "\n"
      << "Smith's test:                     " << m.smiths_test          << "\n"
      << "Smith's test sex-specific:        " << m.smiths_sex_specific  << "\n"
      << "Morton's test:                    " << m.mortons_test         << "\n"
      << "Morton's test sex-specific:       " << m.mortons_sex_specific << "\n";
      
  std::map<string, group>::const_iterator  g_iter;
  for(g_iter = m.groups.begin(); g_iter != m.groups.end(); ++g_iter)
  {
    out << "group " << g_iter->first << " pedigrees:" << endl;
    group::const_iterator  p_iter;
    for(p_iter = g_iter->second.begin(); p_iter != g_iter->second.end(); ++p_iter)
    {
      out << *p_iter << endl;
    }
  }
      
  // - Print in order established by parser.
  //
  for(size_t i = 0; i < m.male_female_thetas.size(); ++i)
  {
    out << m.male_female_thetas[i];
  }
  
  for(size_t i = 0; i < m.average_thetas.size(); ++i)
  {
    out << "Average theta:   " << m.average_thetas[i] << endl;
  }
  
  out << "\nGenotypes:                      " << m.genotypes << "\n"
      << "Genotypes sex-specific:           " << m.genotypes_sex_specific << endl;
  
  out << std::noboolalpha << endl;
  
  return out;
}

//============================================================================
// IMPLEMENTATION:  theta_pair
//============================================================================
//
inline
theta_pair::theta_pair(double ml_theta, double fml_theta)
      : male_theta(ml_theta), female_theta(fml_theta)
{}


//============================================================================
// IMPLEMENTATION:  instructions
//============================================================================
//
inline
instructions::instructions(cerrorstream& errors)
{
  reset();
}

inline void
instructions::reset()
{
  file_name_root = "";
  title          = "";
  
  linkage = TRAIT;
  trait   = "";
  
  linkage_test          = true;
  linkage_sex_specific  = false;
  linkage_homog         = true;
  smiths_test           = false;
  smiths_sex_specific   = false;
  mortons_test          = false;
  mortons_sex_specific  = false;
  groups.clear();

  reset_thetas();
  
  genotypes = false;
  genotypes_sex_specific = false;

  // - Note:  not *really* valid until parser::init_parse() is called and user supplied
  //   trait is set.
  //  
  valid = true;
}

inline void
instructions::reset_thetas()
{
  male_female_thetas.resize(0);
  average_thetas = std::vector<double>(FIRST_THETA, LAST_THETA);
}
    
inline string
instructions::linkage_option_2_string(linkage_option l)
{
  string  option = "";
  
  switch(l)
  {
    case MARKER:
      option = "marker";
      break;
    case TRAIT:
      option = "trait";
      break;
    default:
      assert(false);
  }
  
  return option;
}
    
inline string
instructions::link_option() const
{
  return  instructions::linkage_option_2_string(linkage);
}

// - Does any part of the analysis use sex-specific recombination fractions?
//
inline bool
instructions::sex_specific() const
{
  return  (linkage_test && linkage_sex_specific)   ||
          (smiths_test && smiths_sex_specific)     ||
          (mortons_test && mortons_sex_specific)   ||
          (genotypes && genotypes_sex_specific)    ||
           male_female_thetas.size();
}


