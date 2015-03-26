#include "segreg/mean_split.h"
#include "segreg/model.h"

// changes made so that values vector only includes 
// members not in conditional subset (due to JA)

namespace SAGE
{
namespace SEGREG
{

bool 
mean_split::calculate_trait_statistics(const FPED::Multipedigree& ped_data, const string& trait_name)
{
  clear();

  // Get our trait values

  size_t trait_number = ped_data.info().trait_find(trait_name);

  vector<double> values;

  for(FPED::PedigreeConstIterator
      ped  = ped_data.pedigree_begin();
      ped != ped_data.pedigree_end();
    ++ped)

    for(size_t
        member_idx = 0;
        member_idx < ped->info().member_count();
      ++member_idx)

// this loop skips over members on the conditioned subset
// model::cond_mem_set is empty unless there is ascertainment
    {
      FPED::MemberConstPointer mem = &ped->member_index(member_idx);
      if (SAGE::SEGREG::model::cond_mem_set.size() > 0) {
      set<FPED::MemberConstPointer>::iterator memit;
      memit = SAGE::SEGREG::model::cond_mem_set.find(mem);
        if (memit == SAGE::SEGREG::model::cond_mem_set.end()) {
        double value = ped->info().trait(member_idx,trait_number);
        if(finite(value)) values.push_back(value);
        }
      } else { // no ascertainment 
         double value = ped->info().trait(member_idx,trait_number);
         if(finite(value)) values.push_back(value);
      }
    }

  if(values.size() < 3)
  {
    return my_status = false;
  }
  
  // Sort the vector:
  sort(values.begin(), values.end());

  // If the vector is all the same value
  if(values.front() == values.back())
  {
    y1bar = values[0];
    y2bar = values[0];
    S1    = 0;
    S2    = 0;
    n     = values.size();
    n1    = n / 2;
    n2    = n - n1;
    p     = (1.0 * n1) / n;

    mean    = values[0];
    var     = n;           // This is not true.  Var is actually 0, but we can
                           // always use a greater variance than the true one.
    min_val = values[0];
    max_val = values[0];
    
    return my_status = true;
  }

  // These values won't change, so we only calculate them once.
  n = values.size();

  min_val = values.front();
  max_val = values.back();

  // Find the other values and return if we could find a valid one
  my_status = find_best_stats(values);

  calculate_mean  (values);
  calculate_var   (values);

  return my_status;
}

//===================================================================================
inline bool 
mean_split::find_best_stats(vector<double> & values)
{
  // Calculate the max of the equation.

  double calc_max = -std::numeric_limits<double>::infinity();

  uint bestn = (uint) -1;

  for(uint tn1 = 0; tn1 <= n; ++tn1)
  {
    calculate_values(tn1, values);

    double temp_val = 1.0 * n1 * n2 * (y1bar - y2bar) * (y1bar - y2bar) / (S1 + S2);

//    cout << tn1 << '\t' << n1 << ' ' << n2 << endl;
//    cout << "  " << y1bar << '\t' << y2bar << endl;
//    cout << "  " << S1 << '\t' << S2 << endl;
//    cout << "  " << temp_min << endl;

    if(!SAGE::isnan(temp_val) && temp_val > calc_max)
    {
      calc_max = temp_val;
      bestn    = tn1;
    }      
  }

  if(bestn == (uint) -1)
  {
    clear();
    return false;
  }    

  calculate_values(bestn, values);

  p = (1.0 * n1) / n;

  return true;
}

}
}
