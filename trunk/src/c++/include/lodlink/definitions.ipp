//============================================================================
// File:      definitions.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   10/18/2 - created.                                djb
//                                                                          
// Notes:     
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

inline std::ostream&
operator<<(std::ostream& out, const joint_genotype& jg)
{
  out << jg.tg.allele1().name() << jg.mg.allele1().name() << "/"
      << jg.tg.allele2().name() << jg.mg.allele2().name();
      
  return out;
}

inline std::ostream&
operator<<(std::ostream& out, const haplotype& h)
{
  out << h.ta.name() << "/" << h.ma.name();
  
  return out;
}

//============================================================================
// IMPLEMENTATION:  joint_genotype
//============================================================================
//
inline  
joint_genotype::joint_genotype(MLOCUS::phased_genotype t, MLOCUS::phased_genotype m)
    : tg(t), mg(m)
{}

inline
joint_genotype::joint_genotype(const joint_pen_iter& jpi)
    : tg(jpi.trait_iter.phased_geno()), mg(jpi.marker_iter.phased_geno())
{}

inline double
joint_genotype::frequency() const
{
  return tg.frequency() * mg.frequency();
}

inline void
joint_genotype::print() const
{
  cout << "trait " << tg.name() << "  "
       << "marker " << mg.name();
}

//============================================================================
// IMPLEMENTATION:  haplotype
//============================================================================
//
inline  
haplotype::haplotype(MLOCUS::allele t, MLOCUS::allele m)
    : ta(t), ma(m)
{}

//============================================================================
// IMPLEMENTATION:  joint_pen_iter
//============================================================================
//
inline  
joint_pen_iter::joint_pen_iter(pen_iter ti, pen_iter mi)
    : trait_iter(ti), marker_iter(mi)
{}

inline log_double
joint_pen_iter::penetrance() const
{
  return log_double((*trait_iter) * (*marker_iter));
}


