#include "data/ArgumentRuleset.h"

using namespace SAGE::APP;

char * name = "foo";

void test_single_display(const SAGE::APP::ArgumentRuleset & ruleset)
{
  std::cout << "===========================================================================================" << std::endl
            << "Testing display usage for ruleset:" << std::endl;

  ruleset.dump();
  
  SAGE::APP::ArgumentParser::display_usage(name, ruleset, std::cout);
}

void test_display()
{
  SAGE::APP::ArgumentRuleset ruleset;
  
  ruleset.add_rule(SAGE::APP::ArgumentRuleset::Rule(SAGE::APP::PARAMETER_FILE,   SAGE::APP::ArgumentRuleset::ZERO_OR_ONE));
  ruleset.add_rule(SAGE::APP::ArgumentRuleset::Rule(SAGE::APP::PEDIGREE_FILE,    SAGE::APP::ArgumentRuleset::ZERO_OR_ONE));
  ruleset.add_rule(SAGE::APP::ArgumentRuleset::Rule(SAGE::APP::LOCUS_FILE,       SAGE::APP::ArgumentRuleset::ZERO_OR_ONE));
  ruleset.add_rule(SAGE::APP::ArgumentRuleset::Rule(SAGE::APP::GENOME_FILE,      SAGE::APP::ArgumentRuleset::ZERO_OR_ONE));
  ruleset.add_rule(SAGE::APP::ArgumentRuleset::Rule(SAGE::APP::TRAIT_MODEL_FILE,   SAGE::APP::ArgumentRuleset::ZERO_OR_ONE));
  ruleset.add_rule(SAGE::APP::ArgumentRuleset::Rule(SAGE::APP::IBD_FILE,         SAGE::APP::ArgumentRuleset::ZERO_OR_ONE));
  test_single_display(ruleset);
    
  ruleset.clear();
  ruleset.add_rule(SAGE::APP::ArgumentRuleset::Rule(SAGE::APP::PEDIGREE_FILE,    SAGE::APP::ArgumentRuleset::ONE));
  ruleset.add_rule(SAGE::APP::ArgumentRuleset::Rule(SAGE::APP::PARAMETER_FILE,   SAGE::APP::ArgumentRuleset::ZERO_OR_ONE));
  ruleset.add_rule(SAGE::APP::ArgumentRuleset::Rule(SAGE::APP::LOCUS_FILE,       SAGE::APP::ArgumentRuleset::ZERO_OR_ONE));
  ruleset.add_rule(SAGE::APP::ArgumentRuleset::Rule(SAGE::APP::GENOME_FILE,      SAGE::APP::ArgumentRuleset::ZERO_OR_ONE));
  ruleset.add_rule(SAGE::APP::ArgumentRuleset::Rule(SAGE::APP::TRAIT_MODEL_FILE,   SAGE::APP::ArgumentRuleset::ZERO_OR_ONE));
  ruleset.add_rule(SAGE::APP::ArgumentRuleset::Rule(SAGE::APP::IBD_FILE,         SAGE::APP::ArgumentRuleset::ZERO_OR_ONE));  
  test_single_display(ruleset);
  
  ruleset.clear();
  ruleset.add_rule(SAGE::APP::ArgumentRuleset::Rule(SAGE::APP::PARAMETER_FILE,   SAGE::APP::ArgumentRuleset::ONE));
  ruleset.add_rule(SAGE::APP::ArgumentRuleset::Rule(SAGE::APP::PEDIGREE_FILE,    SAGE::APP::ArgumentRuleset::ONE));
  ruleset.add_rule(SAGE::APP::ArgumentRuleset::Rule(SAGE::APP::TRAIT_MODEL_FILE, SAGE::APP::ArgumentRuleset::ONE));
  ruleset.add_rule(SAGE::APP::ArgumentRuleset::Rule(SAGE::APP::GENOME_FILE,      SAGE::APP::ArgumentRuleset::ZERO_OR_ONE));
  
  test_single_display(ruleset);
}

void test_ex(int argc, char ** argv, const SAGE::APP::ArgumentRuleset & ruleset)
{
  std::cout << "==============================================================" << std::endl
            << "Execution test:" << std::endl;
  
  ruleset.dump();

  std::cout << "Run with parameters: ";

  for(int i = 0; i < argc; ++i)
    std::cout << " " << argv[(size_t)i];
    
  std::cout << std::endl;
  
  SAGE::APP::process_commandline(argc, argv, name, ruleset, std::cout, false).dump();
}

void test_simple()
{
  ArgumentRuleset  ruleset;
  
  ruleset.add_rule(SAGE::APP::ArgumentRuleset::Rule(SAGE::APP::PARAMETER_FILE,   SAGE::APP::ArgumentRuleset::ONE));
  ruleset.add_rule(SAGE::APP::ArgumentRuleset::Rule(SAGE::APP::PEDIGREE_FILE,    SAGE::APP::ArgumentRuleset::ONE));

  int    argc1   = 3;  char * argv1[] = { "freq", "param", "ped"             }; test_ex(argc1, argv1, ruleset); // succeed
  int    argc2   = 5;  char * argv2[] = { "freq", "-p", "param", "-d", "ped" }; test_ex(argc2, argv2, ruleset); // succeed
  int    argc3   = 5;  char * argv3[] = { "freq", "-d", "ped", "-p", "param" }; test_ex(argc3, argv3, ruleset); // succeed
  int    argc4   = 2;  char * argv4[] = { "freq", "param"                    }; test_ex(argc4, argv4, ruleset); // fail
  int    argc5   = 3;  char * argv5[] = { "freq", "-p", "param"              }; test_ex(argc5, argv5, ruleset); // fail
  int    argc6   = 4;  char * argv6[] = { "freq", "param", "ped", "foo"      }; test_ex(argc6, argv6, ruleset); // succeed
  int    argc7   = 1;  char * argv7[] = { "freq"                             }; test_ex(argc7, argv7, ruleset); // fail
}

void test_medium()
{
  ArgumentRuleset ruleset;
  ruleset.add_rule(SAGE::APP::ArgumentRuleset::Rule(SAGE::APP::PARAMETER_FILE,   SAGE::APP::ArgumentRuleset::ONE));
  ruleset.add_rule(SAGE::APP::ArgumentRuleset::Rule(SAGE::APP::PEDIGREE_FILE,    SAGE::APP::ArgumentRuleset::ONE));
  ruleset.add_rule(SAGE::APP::ArgumentRuleset::Rule(SAGE::APP::LOCUS_FILE,       SAGE::APP::ArgumentRuleset::ZERO_OR_ONE));
  ruleset.add_rule(SAGE::APP::ArgumentRuleset::Rule(SAGE::APP::TRAIT_MODEL_FILE, SAGE::APP::ArgumentRuleset::ZERO_OR_ONE));

  int    argc1   = 3;  char * argv1[] = { "freq", "param", "ped"                                         }; test_ex(argc1, argv1, ruleset); // succeed
  int    argc2   = 4;  char * argv2[] = { "freq", "param", "ped", "loc"                                  }; test_ex(argc2, argv2, ruleset); // succeed    
  int    argc3   = 5;  char * argv3[] = { "freq", "param", "ped", "loc", "tld"                           }; test_ex(argc3, argv3, ruleset); // succeed    
  int    argc4   = 5;  char * argv4[] = { "freq", "-p", "param", "-d", "ped"                             }; test_ex(argc4, argv4, ruleset); // succeed
  int    argc5   = 7;  char * argv5[] = { "freq", "-p", "param", "-d", "ped", "-l", "locus"              }; test_ex(argc5, argv5, ruleset); // succeed
  int    argc6   = 7;  char * argv6[] = { "freq", "-p", "param", "-d", "ped", "-m", "tld"                }; test_ex(argc6, argv6, ruleset); // succeed
  int    argc7   = 9;  char * argv7[] = { "freq", "-p", "param", "-d", "ped", "-l", "locus", "-t", "tld" }; test_ex(argc7, argv7, ruleset); // succeed
}

void test_complex()
{
  ArgumentRuleset ruleset;
  ruleset.add_rule(SAGE::APP::ArgumentRuleset::Rule(SAGE::APP::PARAMETER_FILE,   SAGE::APP::ArgumentRuleset::ONE));
  ruleset.add_rule(SAGE::APP::ArgumentRuleset::Rule(SAGE::APP::PEDIGREE_FILE,    SAGE::APP::ArgumentRuleset::ONE));
  ruleset.add_rule(SAGE::APP::ArgumentRuleset::Rule(SAGE::APP::IBD_FILE,         SAGE::APP::ArgumentRuleset::ONE_OR_MORE));
  ruleset.add_rule(SAGE::APP::ArgumentRuleset::Rule(SAGE::APP::LOCUS_FILE,       SAGE::APP::ArgumentRuleset::ZERO_OR_ONE));
  ruleset.add_rule(SAGE::APP::ArgumentRuleset::Rule(SAGE::APP::TRAIT_MODEL_FILE, SAGE::APP::ArgumentRuleset::ZERO_OR_ONE));

  int    argc1   = 3;  char * argv1[] = { "freq", "param", "ped"                                         }; test_ex(argc1, argv1, ruleset); // fail
  int    argc2   = 5;  char * argv2[] = { "freq", "-p", "param", "-d", "ped"                             }; test_ex(argc2, argv2, ruleset); // fail
  int    argc3   = 7;  char * argv3[] = { "freq", "-p", "param", "-d", "ped", "-l", "loc"                }; test_ex(argc3, argv3, ruleset); // fail
  int    argc4   = 7;  char * argv4[] = { "freq", "-p", "param", "-d", "ped", "-i", "ibd1"               }; test_ex(argc4, argv4, ruleset); // succeed
  int    argc5   = 9;  char * argv5[] = { "freq", "-p", "param", "-d", "ped", "-i", "ibd1", "-i", "ibd2" }; test_ex(argc5, argv5, ruleset); // succeed
}


int main()
{
  test_display();

  test_simple();
  
  test_medium();
  
  test_complex();

  return 0;
}

