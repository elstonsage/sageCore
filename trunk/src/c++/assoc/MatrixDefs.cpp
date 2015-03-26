#include "assoc/MatrixDefs.h"

namespace SAGE  {
namespace ASSOC {

double calculate_variance(const PointDensity & pd, const MAXFUN::ParameterMgr & mgr)
{
  // Calculate variance:
  double variance = 0.0;
  
  for(size_t var_idx = 0; var_idx < pd.var_idxs.size(); ++var_idx)
  {
    variance += mgr.getParameter("Variance components", pd.var_idxs[var_idx]).getCurrentEstimate();
  }

  return variance;
}

PointDensityMatrix generate_pdm(const PointDensity & pd)
{
  // Create the PointDensityMatrix:
  PointDensityMatrix pdm;
  
  // Set the name:
  pdm.set_name(pd.name);

  // Resize it accordingly:
  pdm.get_matrix().resize_fill(pd.coeff_pairs.size(), pd.coeff_pairs.size(), 0.0);

  // Figure out the ids:
  pdm.get_idxs().resize(pd.coeff_pairs.size());
  
  for(size_t coeff_idx = 0; coeff_idx < pd.coeff_pairs.size(); ++coeff_idx)
  {
    pdm.get_idxs()[coeff_idx] = pd.coeff_pairs[coeff_idx].term_idx;
  }

  // Return the PointDensityMatrix:
  return pdm;
}

void PointDensity::dump(const TermNameVector & term_names) const
{
  std::cout << "PointDensity dump" << std::endl;
  
  std::cout << "  Name: " << name << std::endl;
      
  std::cout << "  CoeffPairs: ";
      
  for(size_t i = 0; i < coeff_pairs.size(); ++i)
    std::cout << (i ? ", " : "") << coeff_pairs[i].coefficient << " * " << term_names[coeff_pairs[i].term_idx];
        
  std::cout << std::endl;
      
  std::cout << "  Variance components: ";
      
  for(size_t i = 0; i < var_idxs.size(); ++i)
    std::cout << var_idxs[i] << " ";
        
  std::cout << std::endl << std::endl;
}

void PointDensityMatrix::dump() const
{
  size_t term_count = 0;

  for(size_t i = 0; i < get_idxs().size(); ++i)
    term_count = get_idxs()[i] + 1 > term_count ? get_idxs()[i] + 1 : term_count;

  TermNameVector t(term_count);
  
  for(size_t i = 0; i <  t.size(); ++i)
  {
    std::ostringstream s; s << "Term #" << i; t[i] = s.str();
  }
  
  dump(t);
}

void PointDensityMatrix::dump(const TermNameVector & term_names) const
{
  OUTPUT::Table t("PointDensityMatrix: " + my_name);
      
  t << OUTPUT::TableColumn("");
      
  for(size_t i = 0; i < my_idxs.size(); ++i)
  {
    std::ostringstream title;
       
    title << term_names[my_idxs[i]];

    t << OUTPUT::TableColumn(title.str());
  }
      
  for(size_t i = 0; i < my_matrix.rows(); ++i)
  {
    OUTPUT::TableRow row;
        
    std::ostringstream title;
        
    title << term_names[my_idxs[i]];

    row << title.str();
        
    for(size_t j = 0; j < my_matrix.rows(); ++j)
      row << my_matrix(i, j);
          
    t << row;
  }

  std::cout << t << std::endl;
}

} // End namespace ASSOC
} // End namespace SAGE

