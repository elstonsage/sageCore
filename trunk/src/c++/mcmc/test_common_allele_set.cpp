#include "mcmc/common_allele_set.h"
#include "mcmc/founder_allele_graph.h"

#include <iostream>
#include <iomanip>

void print_csa_status(const std::string& state_info,
                      const SAGE::MCMC::CommonAlleleSet& csa)
{
  std::cout << "CSA in state:    " << state_info << std::endl;
  std::cout << "Max Alleles:     " << csa.get_max_allele_count() << std::endl;
  std::cout << "Is Dirty:        " << csa.is_dirty() << std::endl;
  std::cout << "Current State:   " << csa.get_status() << std::endl;
  
  if(csa.get_allele_id(0) != SAGE::MLOCUS::NPOS)
    std::cout << "Allele 0:        " << csa.get_allele_id(0) << std::endl;
  else
    std::cout << "Allele 0:        MLOCUS::NPOS" << std::endl;
  if(csa.get_allele_id(1) != SAGE::MLOCUS::NPOS)
    std::cout << "Allele 1:        " << csa.get_allele_id(1) << std::endl;
  else
    std::cout << "Allele 1:        MLOCUS::NPOS" << std::endl;

  std::cout << "Number of pairs: " << csa.get_allele_pair_count() << std::endl
            << std::endl;

  for(size_t i = 0; i < csa.get_max_allele_count(); ++i)
  {
    std::cout << "Num of pairs for Allele #" << std::setw(2) << i << ": "
              << csa.get_allele_pair_count_for_allele(i) << std::endl;
  }
  std::cout << std::endl;

  std::cout << "Is Dirty:        " << csa.is_dirty() << std::endl;

  std::cout << std::endl;
}

// This program tests the CommonAlleleSet object.  Since that is a fairly
// simple class, we don't need the app or lots of data, so this test can be
// kept simple.

main()
{
  
  // Create a common Allele set
  SAGE::MCMC::CommonAlleleSet csa(4);

  // Print out set status
  print_csa_status("Initialization", csa);
  
  // Stuff a single pair into it
  std::cout << "Potential Change Flag: " << csa.add_allele_pair(2,1) << std::endl << std::endl;
  
  // Print out set status
  print_csa_status("One Pair", csa);
  
  // Stuff a few more (identical) pairs into it
  std::cout << "Potential Change Flag: " << csa.add_allele_pair(1,2) << std::endl << std::endl;
  std::cout << "Potential Change Flag: " << csa.add_allele_pair(2,1) << std::endl << std::endl;
  std::cout << "Potential Change Flag: " << csa.add_allele_pair(1,2) << std::endl << std::endl;
  
  // Print out set status
  print_csa_status("Four Identical Pairs", csa);
  
  // Remove a pair
  std::cout << "Potential Change Flag: " << csa.remove_allele_pair(1,2) << std::endl << std::endl;
  
  // Print out set status
  print_csa_status("Removed one pair", csa);
  
  // Stuff a (non-identical) pair into it
  std::cout << "Potential Change Flag: " << csa.add_allele_pair(1,3) << std::endl << std::endl;

  // Print out set status
  print_csa_status("Non-Identical Pairs (one allele still valid)", csa);
  
  // Remove one of the identical pairs
  std::cout << "Potential Change Flag: " << csa.remove_allele_pair(1,2) << std::endl << std::endl;
  
  // Print out set status
  print_csa_status("Removed (1,2) pair", csa);
  
  // Remove the non-identical pair
  std::cout << "Potential Change Flag: " << csa.remove_allele_pair(1,3) << std::endl << std::endl;
  
  // Print out set status
  print_csa_status("Remove (1,3) pair", csa);
  
  // Stuff a pair that makes it invalid
  std::cout << "Potential Change Flag: " << csa.add_allele_pair(3,0) << std::endl << std::endl;
  
  // Print out set status
  print_csa_status("Invalid set", csa);
  
  // Remove pair that makes it invalid
  std::cout << "Potential Change Flag: " << csa.remove_allele_pair(3,0) << std::endl << std::endl;

  // Print out set status
  print_csa_status("Invalid set valid again", csa);

  // Call test complete
  std::cout << "Testing Complete!" << std::endl << std::endl;
}
