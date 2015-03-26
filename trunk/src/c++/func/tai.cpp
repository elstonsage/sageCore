//=============================================================================
// File:      tai.cpp (Transmitted Allele indicator)
//
// Author:    Kai He
//
//
// History:   10/2003 Initial version
//
// Notes:
//
// Copyright (c) 2002 R.C. Elston
// All Rights Reserved
//============================================================================

#include "func/tai.h"
#include <iostream>

namespace SAGE {
namespace FUNC {

//----------------------------------------------------------------------------
//
// Constructor
//
//---------------------------------------------------------------------------- 
transmitted_allele_indicator::transmitted_allele_indicator(cerrorstream & err) : my_errors(err)
{
  my_marker_name = "";
  my_allele_name = "";
}

//============================================================================
//
//  Copy constructor
//
//============================================================================
transmitted_allele_indicator::transmitted_allele_indicator(const transmitted_allele_indicator & other)
{
  my_errors      = other.my_errors;
  my_marker_name = other.my_marker_name;
  my_allele_name = other.my_allele_name;
  my_type        = other.my_type;
}

//============================================================================
//
//  operator=(...)
//
//---------------------------------------------------------------------------- 
transmitted_allele_indicator & 
transmitted_allele_indicator::operator=(const transmitted_allele_indicator & other)
{
  if(this != &other)
  {
    my_errors      = other.my_errors;
    my_marker_name = other.my_marker_name;
    my_allele_name = other.my_allele_name;
    my_type        = other.my_type;
  }

  return *this;
}

//=============================================================================
//
//  Destructor
//
//=============================================================================
transmitted_allele_indicator::~transmitted_allele_indicator()
{}

//----------------------------------------------------------------------------
//
//  createTaiStatusTrait(...)
//
//---------------------------------------------------------------------------- 
void
transmitted_allele_indicator::createTaiStatusTrait(RPED::RefMultiPedigree & mp, const FunctionParser::TraitData& pd)
{
  // Parses the given function block (through pd); assigns marker name and allele name.

  parseTaiExpression(pd.expr);

  // Verify that the marker/allele names are good:

  if(verifyParameters(mp) == false)
    return;

  // Populate the RefMultiPedigree with TAI status variables:

  populateRefMultiPedigree(mp, pd);
}

//==============================================================================
//
//  parseTaiExpression(...)
//
//==============================================================================
void 
transmitted_allele_indicator::parseTaiExpression (string expr)
{
  if(check_parenthese(expr) != 1)
  {
     cout << endl
          << "One or more parenthese(s) is(are) missing in the expression of "
          << "the function block(s), check and try again!"
          << endl;

     exit(1);
  }

  // get my_process

  string process_name;
  size_t first = expr.find('(',0);
  for(size_t i=0; i<first;++i)
      process_name+=expr[i];

  my_type = process_name == "tai" ? TAI : UTAI;

  // remove parentheses

  string str = deparenthese(expr);
  parameter_list plst;
  plst.clear();

  if(str.size() > 0)
  {
    char* ch  = &str[0];
    char* pch = strtok (ch,",");

    while (pch != NULL)
    {
      plst.push_back(pch);
      pch = strtok (NULL, " ,.");
    }
  }

  my_marker_name = *plst.begin();

  plst.pop_front();

  my_allele_name = *plst.begin();
}

//----------------------------------------------------------------------------
//
// check_parenthese(...)
//
//---------------------------------------------------------------------------- 
size_t 
transmitted_allele_indicator::check_parenthese(string seq)
{
  size_t c1,c2; c1=c2=0;
  for(size_t i=0; i< seq.size(); ++i)
  {
    if(seq[i]=='(') c1++;
    if(seq[i]==')') c2++;
  }
  if(c1==c2)
     return 1;
  else if(c1>c2)
     return 2;
  else if(c1<c2)
     return 3;
  else
     return 0;
}

//============================================================================
//
//  get_allele_1(...)
//
//============================================================================
string
transmitted_allele_indicator::get_allele_1(string pheno_name)
{
  string a;
  size_t pos = pheno_name.find('/',0);

  for(size_t k = 0; k < pos; ++k)
    a += pheno_name[k];

  return a;
}

//=============================================================================
//
//  get_allele_2(...)
//
//=============================================================================
string
transmitted_allele_indicator::get_allele_2(string pheno_name)
{
  string a;
  size_t pos = pheno_name.find('/',0);
  for(size_t k = pos + 1; k < pheno_name.size(); ++k)
    a += pheno_name[k];

  return a;
}  

//----------------------------------------------------------------------------
//
//
//
//---------------------------------------------------------------------------- 

bool
transmitted_allele_indicator::is_homozygous(string a1, string a2, string ale)
{
  return a1 == ale && a2 == ale;
}

//----------------------------------------------------------------------------
//
//
//
//---------------------------------------------------------------------------- 

bool
transmitted_allele_indicator::is_heterozygous(string a1, string a2, string ale)
{
  return (a1 == ale && a2 != ale) || (a2 == ale && a1 != ale);
}

//----------------------------------------------------------------------------
//
//
//
//---------------------------------------------------------------------------- 

bool 
transmitted_allele_indicator::allele_in_genotype(string c1, string p1, string p2)
{
  return (c1 == p1 || c1 == p2);
}


//----------------------------------------------------------------------------
//
//
//
//---------------------------------------------------------------------------- 
string 
transmitted_allele_indicator::deparenthese(string str)
{
  string str2;
  size_t first = str.find('(',0),
         last  = str.find(')',0);

  for(size_t i = first + 1; i < last; ++i)
      str2 += str[i];

  return str2;
}

//=============================================================================
//
//  verifyParameters(...)
//
//=============================================================================
bool 
transmitted_allele_indicator::verifyParameters(const RPED::RefMultiPedigree & mp)
{
  const RPED::RefMPedInfo & mp_info = mp.info();

  if(mp_info.marker_exists(my_marker_name))
  {
    const RPED::RefMarkerInfo & marker_info = mp_info.marker_info(mp_info.marker_find(my_marker_name));

    if(marker_info.gmodel().get_allele(my_allele_name).is_valid())
    {
      return true;
    }
    else
    {
      my_errors << priority(error)   
                << "Allele '"    
                << my_allele_name
                << "' does not exist for marker '"
                << my_marker_name
                << "' (used in a tai(...) function block)."
                << endl;

      return false;
    }
  }
  else
  {
    my_errors << priority(error)
              << "Marker '"
              << my_marker_name
              << "' does not exist (used in a tai(...) function block)."
              << endl;

    return false;
  }
}

//----------------------------------------------------------------------------
//
//  getTaiStatus(...)
//
//---------------------------------------------------------------------------- 
double
transmitted_allele_indicator::getTaiStatus(const RPED::RefMember & mem)
{
  // First, make sure that this member is a non-founder:

  if(MPED::mp_utilities::is_founder(mem) == true)
  {
    return std::numeric_limits<double>::quiet_NaN();
  }

  // Now that we know the individual is a non-founder, we need to make
  // sure that the three people are fully informative:

        size_t                marker_id   = mem.multipedigree()->info().marker_find (my_marker_name);
  const RPED::RefMarkerInfo & marker_info = mem.multipedigree()->info().marker_info (marker_id);
        size_t                parent1_id  = mem.parent1()->index();
        size_t                parent2_id  = mem.parent2()->index();
        size_t                child_id    = mem.index(); 
  
  if(mem.pedigree()->info().phenotype_missing(parent1_id, marker_id, marker_info) ||
     mem.pedigree()->info().phenotype_missing(parent2_id, marker_id, marker_info) ||
     mem.pedigree()->info().phenotype_missing(child_id,   marker_id, marker_info))
  {
    return std::numeric_limits<double>::quiet_NaN();
  }

  // Now that we know that everyone is fully informative, we need to make
  // sure that the homozygous/heterozygous mating combinations are valid for
  // calculating a TAI status (see interface documentation for more details).

  enum genotype { HOMOZYGOUS, HETEROZYGOUS };

  // Pull out parent 1 information:

  uint     phenotype_id_p1   = mem.pedigree()->info().phenotype(parent1_id, marker_id);
  string   phenotype_name_p1 = marker_info.get_phenotype(phenotype_id_p1).name(),
           allele1_p1        = get_allele_1(phenotype_name_p1),
           allele2_p1        = get_allele_2(phenotype_name_p1);
  genotype genotype_p1       = allele1_p1 == allele2_p1 ? HOMOZYGOUS : HETEROZYGOUS;

  // Pull out parent 2 information:

  uint     phenotype_id_p2   = mem.pedigree()->info().phenotype(parent2_id, marker_id);
  string   phenotype_name_p2 = marker_info.get_phenotype(phenotype_id_p2).name(),
           allele1_p2        = get_allele_1(phenotype_name_p2),
           allele2_p2        = get_allele_2(phenotype_name_p2);
  genotype genotype_p2       = allele1_p2 == allele2_p2 ? HOMOZYGOUS : HETEROZYGOUS;

  // Pull out child information:

  uint     phenotype_id_child   = mem.pedigree()->info().phenotype(child_id, marker_id);
  string   phenotype_name_child = marker_info.get_phenotype(phenotype_id_child).name(),
           allele1_child        = get_allele_1(phenotype_name_child),
           allele2_child        = get_allele_2(phenotype_name_child);
  genotype genotype_child       = allele1_child == allele2_child ? HOMOZYGOUS : HETEROZYGOUS;

  // Check the various homozygous/heterozygous combinations:

  if(((genotype_p1 == HOMOZYGOUS)   && (genotype_p2 == HOMOZYGOUS)) ||
     ((genotype_p1 == HETEROZYGOUS) && (genotype_p2 == HETEROZYGOUS) && (genotype_child == HETEROZYGOUS)))
  {
    return std::numeric_limits<double>::quiet_NaN();
  }

  // Ok, now we now that the member is a nonfounder, everyone is
  // informative, and we have valid homozygous/heterozygous combinations. Now we
  // calculate the TAI status variable:

  double tai_status;

  if(is_homozygous(allele1_p1, allele2_p1, my_allele_name) || is_homozygous(allele1_p2, allele2_p2, my_allele_name))
  {
    if(allele1_child == my_allele_name && allele2_child == my_allele_name)
      tai_status = (my_type == TAI ? 1.0 : 0.0);
    else
      tai_status = (my_type == TAI ? 0.0 : 1.0);
  }
  else
  {
    if(allele1_child == my_allele_name || allele2_child == my_allele_name)
      tai_status = (my_type == TAI ? 1.0 : 0.0);
    else
      tai_status = (my_type == TAI ? 0.0 : 1.0);
  }

  return tai_status;
}

//----------------------------------------------------------------------------
//
//  populateRefMultiPedigree(...)
//
//---------------------------------------------------------------------------- 
void 
transmitted_allele_indicator::populateRefMultiPedigree(RPED::RefMultiPedigree& mp, const FunctionParser::TraitData& pd)
{
  // 1. Create the trait entry in both the RefMPedInfo:

  RPED::RefMPedInfo & mp_info   = mp.info();
  size_t              trait_num = mp_info.add_continuous_trait(pd.trait_name, pd.usage);

  mp_info.trait_info(trait_num).set_string_missing_code("");
  mp_info.trait_info(trait_num).set_numeric_missing_code(-1);

  // 2. Populate the RefMultiPedigree with values:

  for(RPED::RefMultiPedigree::pedigree_iterator ped = mp.pedigree_begin(); ped != mp.pedigree_end(); ++ped)
  {
    ped->info().resize_traits(ped->info().trait_count() + 1);

    for(RPED::RefMultiPedigree::member_const_iterator mem = ped->member_begin(); mem != ped->member_end(); ++mem)
      ped->info().set_trait(mem->index(), trait_num, getTaiStatus(*mem));
  }
}

} // End namespace FUNC
} // End namespace SAGE
