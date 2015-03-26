//============================================================================
// File:      pedparse.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   Created 11/00                                       
//                                                                          
// Notes:     Implementation of the following classes -
//              pedinfo_parser
//                                                                          
// Copyright (c) 2000 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "pedinfo/pedparse.h"

using namespace SAGE;

//============================================================================
// IMPLEMENTATION:  pedinfo_parser
//============================================================================
//
pedinfo_parser::pedinfo_parser()
{
  clear();
}

// - Check attributes of pedinfo_analysis parameters in parameter file.
//   Set data members according to attributes.  Generate error messages 
//   as necessary.
//
void
pedinfo_parser::parse(const LSFBase* analysis, cerrormultistream& errors,
                      const RPED::RefMultiPedigree* mp)
{
  clear();

  if(!analysis || !mp)
  {
    return;
  }
  
  // - Output.
  //
  AttrList* a_list = analysis->attrs();
  if(a_list)
  {
    string  out_value = a_list->StringAttr("OUT");
    string  output_value = a_list->StringAttr("OUTPUT");
    if(out_value.size())
    {
      my_file_name = out_value + ".out";
    }
    else if(output_value.size())
    {
      my_file_name = output_value + ".out";
    }
  }
  
  const LSFList* p_list = analysis->List();
  if(p_list == 0)
  {
    errors << priority(warning) << "No parameters read for Pedinfo analysis block." << endl;
    
    return;
  }
  
  LSFList::const_iterator   a_iter;
  AttrVal                   attribute;
  
  // - Look at each attribute of the analysis parameter.
  //
  for(a_iter = p_list->begin(); a_iter != p_list->end(); ++a_iter)
  {
    if(! *a_iter)
    {
      continue;
    }
    
    string attribute_name = toUpper((*a_iter)->name());
    
    // - each_pedigree.
    //
    if(attribute_name == "EACH_PEDIGREE")
    {
      attribute = attr_value(*a_iter, "EACH_PEDIGREE", 0);
      if(attribute.has_value())
      {
        if(toUpper(attribute.String()) == "TRUE")
        {
          my_show_each_pedigree = true;
        }
        else if(toUpper(attribute.String()) == "FALSE")
        {
          my_show_each_pedigree = false;
        }
        else
        {
          cout << endl;
          errors << priority(error) << "Bad value, '" + attribute.String() + "', for 'each_pedigree'"
                                    << " sub-parameter of 'pedinfo_analysis' parameter.\n"
                                    << "Default value, 'false', was used." << endl; 
        }
      }
      else
      {
        cout << endl;
        errors << priority(error) << "Missing value for 'each_pedigree' sub-parameter of "
                                  << "'pedinfo_analysis' parameter.\n"
                                  << "Default value, 'false', was used." << endl; 
      }
    }
    
    // - suppress_general.
    //
    else if(attribute_name == "SUPPRESS_GENERAL")
    {
      attribute = attr_value(*a_iter, "SUPPRESS_GENERAL", 0);
      if(attribute.has_value())
      {
        if(toUpper(attribute.String()) == "TRUE")
        {
          my_suppress_general = true;
        }
        else if(toUpper(attribute.String()) == "FALSE")
        {
          my_suppress_general = false;
        }
        else
        {
          cout << endl;
          errors << priority(error) << "Bad value, '" + attribute.String() + "', for 'suppress_general'"
                                    << " sub-parameter of 'pedinfo_analysis' parameter.\n"
                                    << "Default value, 'false', was used." << endl; 
        }
      }
      else
      {
        cout << endl;
        errors << priority(error) << "Missing value for 'suppress_general' sub-parameter of "
                                  << "'pedinfo_analysis' parameter.\n"
                                  << "Default value, 'false', was used." << endl; 
      }
    }    
    
    // - trait.
    //
    else if(attribute_name == "TRAIT")
    {
      attribute = attr_value(*a_iter, "TRAIT", 0);
      if(attribute.has_value())
      {
        size_t trait = mp->info().trait_find(attribute.String());
        if(trait == (size_t)(-1))
        {
          cout << endl;
          errors << priority(error) << "Bad value, '" + attribute.String() + "', for 'trait'"
                                    << " sub-parameter of 'pedinfo_analysis' parameter.\n"
                                    << "Sub-parameter ignored." << endl; 
        }
        else
        {
          my_traits.push_back(trait); 
        }
      }
      else
      {
        cout << endl;
        errors << priority(error) << "Missing value for 'trait' sub-parameter of "
                                  << "'pedinfo_analysis' parameter.\n"
                                  << "Sub-parameter ignored." << endl;
      }
    }
    
    // - covariate.
    //
    else if(attribute_name == "COVARIATE")
    {
      attribute = attr_value(*a_iter, "COVARIATE", 0);
      if(attribute.has_value())
      {
        size_t trait = mp->info().trait_find(attribute.String());
        if(trait == (size_t)(-1))
        {
          cout << endl;
          errors << priority(error) << "Bad value, '" + attribute.String() + "', for 'covariate'"
                                    << " sub-parameter of 'pedinfo_analysis' parameter.\n"
                                    << "Sub-parameter ignored." << endl; 
        }
        else
        {
          my_traits.push_back(trait); 
        }
      }
      else
      {
        cout << endl;
        errors << priority(error) << "Missing value for 'covariate' sub-parameter of "
                                  << "'pedinfo_analysis' parameter.\n"
                                  << "Sub-parameter ignored." << endl;
      }
    }
    
    // - phenotype.
    //
    else if(attribute_name == "PHENOTYPE")
    {
      attribute = attr_value(*a_iter, "PHENOTYPE", 0);
      if(attribute.has_value())
      {
        size_t trait = mp->info().trait_find(attribute.String());
        if(trait == (size_t)(-1))
        {
          cout << endl;
          errors << priority(error) << "Bad value, '" + attribute.String() + "', for 'phenotype'"
                                    << " sub-parameter of 'pedinfo_analysis' parameter.\n"
                                    << "Sub-parameter ignored." << endl; 
        }
        else
        {
          my_traits.push_back(trait); 
        }
      }
      else
      {
        cout << endl;
        errors << priority(error) << "Missing value for 'phenotype' sub-parameter of "
                                  << "'pedinfo_analysis' parameter.\n"
                                  << "Sub-parameter ignored." << endl;
      }
    }
    else
    {
      cout << endl;
      errors << priority(error) << "Bad sub-parameter, '" + (*a_iter)->name() + "' for "
                                << "'pedinfo_analysis' parameter.\n"
                                << "Sub-parameter ignored." << endl;
    }
  }
}

void
pedinfo_parser::clear()
{
  my_traits.clear();
  my_show_each_pedigree = false;
  my_suppress_general = false;
  my_file_name = "pedinfo.out";
}

