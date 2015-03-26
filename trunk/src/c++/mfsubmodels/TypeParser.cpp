#include "mfsubmodels/TypeParser.h"

namespace SAGE        {
namespace MFSUBMODELS {

//=============================================================
//
//  parse(...)
//
//=============================================================
bool
TypeSpecificParser::parse(TypeSpecificSubmodel & submodel, const LSFBase * param, cerrorstream & errors)
{
  if(!param->List())
    return false;

  TypeSpecificSubmodel temp_submodel(submodel.getCategory());

  for(LSFList::const_iterator iter = param->List()->begin(); iter != param->List()->end(); ++iter)
  {
    if(!(*iter))
      continue;

    string param_name = toUpper((*iter)->name());

    if(param_name == "OPTION") 
    {
      parseOption(temp_submodel, *iter, errors);
    }
    else if(param_name == toUpper(temp_submodel.getBriefName()))
    {
      parseGenotype(temp_submodel, *iter, errors);
    }
    else
    {
      errors << priority(error) << "Unrecognized parameter '" << param_name << "' in " << getSubblockName(submodel) << " sub-block." << endl;
    }
  }

  int warning_code = 0;

  int err_code = temp_submodel.isInitialSetupConsistent(warning_code);

  if(warning_code != 0)
  {
    string message = "";

    if(warning_code == 1)
    {
      message = "OPTION=TWO, AA and BB are given but not equal; interpreting as " + TypeConstraint(submodel.getDefaultTwoTreatment()).toString();
    }
    else if(warning_code == 2)
    {
      message = "OPTION=TWO, AB and BB are given but not equal; interpreting as TWO_DOM";
    }

    errors << priority(warning)
           << "Unclear statement in " 
           << getSubblockName(submodel)
           << " subblock:"
           << endl
           << message
           << endl;
  }

  if(err_code == 0)
  {
    submodel = temp_submodel;

    return true;
  }
  else
  {
    string message = "";

    if(err_code == 1)
    {
      message = "OPTION=ONE, but at least two different values given";
    }
    else if(err_code == 2)
    {
      message = "OPTION=TWO_DOM, but different values given for AA and AB";
    }
    else if(err_code == 3)
    {
      message = "OPTION=TWO_DOM, valid value given for dominant type (AA/AB), but no value given for recessive type (BB)";
    }
    else if(err_code == 4)
    {
      message = "OPTION=TWO_REC, but different values given for AB and BB";
    }
    else if(err_code == 5)
    {
      message = "OPTION=TWO_REC, valid value given for dominant type (AB/BB), but no value given for recessive type (AA)";
    }
    else if(err_code == 6)
    {
      message = "OPTION=TWO/TWO_DOM/TWO_REC, but one one value specified";
    }
    else if(err_code == 7)
    {
      message = "OPTION=THREE_ADD, three values given, but do not meet relational requirement";
    }
    else if(err_code == 8)
    {
      message = "OPTION=THREE_DEC, three values given, but do not meet relational requirement";
    }
    else if(err_code == 9)
    {
      message = "OPTION=THREE_INC, three values given, but do not meet relational requirement";
    }
    else if(err_code == 10)
    {
      message = "OPTION=THREE_ADD/THREE_DEC/THREE_INC, at least one value given but not all three";
    }
    else if(err_code == 11)
    {
      message = "OPTION=TWO, but three different values given";
    }
    else if(err_code == 12)
    {
      message = "OPTION=TWO, AA equals AB, but BB not given; not sure how to interpret this";
    }
    else if(err_code == 13)
    {
      message = "OPTION=TWO, AB equals BB, but AA not given; not sure how to interpret this";
    }
    else if(err_code == 14)
    {
      message = "OPTION=TWO, AA equals BB, and AB != AA or BB; not sure how to interpret this";
    }
    else if(err_code == 100)
    {
      message = "OPTION=ONE, at least two different fixedness statuses given";
    }
    else if(err_code == 101)
    {
      message = "OPTION=TWO (interpreted as TWO_DOM) or TWO_DOM, conflicting fixed status given for AA and AB";
    }
    else if(err_code == 102)
    {
      message = "OPTION=TWO (interpreted as TWO_REC) or TWO_REC, conflicting fixed status given for AB and BB";
    }
    else if(err_code == 200)
    {
      message = "Variance sub-block indicated, value less than zero given";
    }

    errors << priority(error)
           << "Inconsistency in " 
           << getSubblockName(submodel)
           << " subblock:"
           << endl
           << message
           << endl
           << "Ignoring subblock..."
           << endl;
  }

  return false;
}

//============================================================
//
//  areMeanAndVarianceCompatible(...)
//
//============================================================
bool 
TypeSpecificParser::areMeanAndVarianceCompatible(
  const TypeSpecificSubmodel & mean_submodel, 
  const TypeSpecificSubmodel & var_submodel,
  cerrorstream & errors)
{
  int error_code = TypeSpecificSubmodel::areMeanAndVarianceCompatible(mean_submodel, var_submodel);

  if(error_code != 0)
  {
    string message = "";

    if(error_code == 1)
    {
      message = "Number of variances ("       + var_submodel  .getTypeConstraint().getTypeCountAsString() +
                ") exceeds number of means (" + mean_submodel .getTypeConstraint().getTypeCountAsString() + ")";
    }

    errors << priority(error)
           << "Inconsistency in " 
           << getSubblockName(mean_submodel)
           << " / "
           << getSubblockName(var_submodel)
           << " subblocks:"
           << endl
           << message
           << endl;

    return false;
  }

  return true;
}


//============================================================
//
//  getSubblockName(...)
//
//============================================================
string 
TypeSpecificParser::getSubblockName(const TypeSpecificSubmodel & submodel)
{
  return "type_" + submodel.getBriefName();
}

//============================================================
//
//  parseOption(...)
//
//============================================================
void 
TypeSpecificParser::parseOption(TypeSpecificSubmodel & submodel, const LSFBase * param, cerrorstream & errors)
{
  if(!param)
    return;

  if(toUpper(param->name()) != "OPTION")
    return;

  // If this is an "OPTION" parameter, pull out the option value:

  AttrVal a = attr_value(param, 0);

  if(a.has_value() == true)
  {
    // Fetch the ConstraintType, check it, and set it:

    TypeConstraint t;

    if(t.setOption(toUpper(a.String())) ? !submodel.setTypeConstraint(t) : true)
    {
      errors << priority(error) 
             << "Value '" 
             << a.String()
             << "'for"
             << " option parameter of type_mean or type_susceptibility" 
             << " sub-block not recognized or not allowed. "
             << " Ignoring ..." 
             << endl;
    }
  }
  else
  {
    errors << priority(error)
           << "There is no value given for OPTION. Ignoring..."
           << endl;
  }
}

//============================================================
//
//  parseGenotype(...)
//
//============================================================
void 
TypeSpecificParser::parseGenotype(TypeSpecificSubmodel & submodel, const LSFBase * param, cerrorstream & errors)
{
  if(!param)
    return;

  if((toUpper(param->name()) != toUpper(submodel.getBriefName())) || !param->attrs())
    return;

  // Set up variables that will/may be populated with data:

  vector<size_t> genotypes(0);
  double         initial_val     = std::numeric_limits<double>::quiet_NaN();
  bool           fixedness       = false,
                 genotype_set    = false,
                 initial_val_set = false,
                 fixedness_set   = false;

  // Loop through attribute list:

  for(AttrList::const_iterator iter = param->attrs()->begin(); iter != param->attrs()->end(); ++iter)
  {
    string attr_name = toUpper(AttrNameMgr.query(iter->first));

         if(attr_name == "VALUE") parseGenotypeValue (genotypes,   genotype_set,    iter, errors);
    else if(attr_name == "VAL")   parseGenotypeVal   (initial_val, initial_val_set, iter, errors);
    else if(attr_name == "FIXED") parseGenotypeFixed (fixedness,   fixedness_set,   iter, errors);
    else
    {
      errors << priority(error) << "Invalid value '" << attr_name << "' given as an attribute " << param->name() << endl;
    }
  }
   
  // Having looped through all the attributes, make sure we got the ones we needed:

  if(genotype_set == false)
  {
    errors << priority(error) << "No genotype given for '" << param->name() << "'" << endl;
  }
  else if(fixedness_set && !initial_val_set)
  {
    errors << priority(error) << "Fixedness set for parameter '" << param->name() << "', but val not indicated." << endl;
  }
  else if(!initial_val_set && !fixedness_set)
  {
    errors << priority(error) << "Parameter '" << param->name() << "' indicated, but no attributes given." << endl;
  }
  else
  {
    for(vector<size_t>::iterator i = genotypes.begin(); i != genotypes.end(); ++i)
    {
      submodel.setTypeValues(*i, initial_val, fixedness_set, fixedness);
    }
  }
}

//============================================================
//
//  parseGenotypeValue(...)
//
//============================================================
void 
TypeSpecificParser::parseGenotypeValue(
  vector<size_t>           & genotypes,
  bool                     & genotype_set,
  AttrList::const_iterator   iter, 
  cerrorstream             & errors)
{
  genotype_set = false;

  string attr_name = toUpper(AttrNameMgr.query(iter->first));

  if(attr_name != "VALUE")
    return;

  AttrVal a = iter->second;

  if(a.has_value())
  {
    string genotype_string = toUpper(a.String());

    genotype_set = true;

         if(genotype_string == "AA")   genotypes.push_back(0);
    else if(genotype_string == "AB")   genotypes.push_back(1);
    else if(genotype_string == "BB")   genotypes.push_back(2);
    else if(genotype_string == "A*") { genotypes.push_back(0); 
                                       genotypes.push_back(1); }
    else if(genotype_string == "B*") { genotypes.push_back(1); 
                                       genotypes.push_back(2); }
    else if(genotype_string == "**") { genotypes.push_back(0); 
                                       genotypes.push_back(1); 
                                       genotypes.push_back(2); }
    else
    {
      genotype_set = false;

      errors << priority(error) << "Invalid value '" << genotype_string << "' given, only 'AA', 'AB', 'BB', 'A*', 'B*', and '**' allowed." << endl;
    }
  }
  else
  {
    errors << priority(error) << "No value specified for VALUE attribute. Ignoring..." << endl;
  }
}

//============================================================
//
//  parseGenotypeVal(...)
//
//============================================================
void 
TypeSpecificParser::parseGenotypeVal(
  double                   & initial_val, 
  bool                     & initial_val_set, 
  AttrList::const_iterator   iter, 
  cerrorstream             & errors)
{
  initial_val_set = false;

  string attr_name = toUpper(AttrNameMgr.query(iter->first));

  if(attr_name != "VAL")
    return;

  AttrVal a = iter->second;

  if(a.has_value())
  {
    if(finite(a.Real()))
    {
      initial_val     = a.Real();
      initial_val_set = true;
    }   
    else
    {
      errors << priority(error) << "Value of " << a.String() << " not understood.  Ignoring ..." << endl;
    }
  }
  else
  {
    errors << priority(error) << "No value specified for VAL attribute. Ignoring..." << endl;
  }
}

//============================================================
//
//  parseGenotypeFixed(...)
//
//============================================================
void 
TypeSpecificParser::parseGenotypeFixed(
  bool                     & fixedness,
  bool                     & fixedness_set, 
  AttrList::const_iterator   iter, 
  cerrorstream             & errors)
{
  fixedness_set = false;

  string attr_name = toUpper(AttrNameMgr.query(iter->first));

  if(attr_name != "FIXED")
    return;

  AttrVal a = iter->second;

  if(a.has_value())
  {
    string fixedness_string = toUpper(a.String());

    if(fixedness_string == "TRUE")
    {
      fixedness_set = true;
      fixedness     = true;
    }
    else if(fixedness_string == "FALSE")
    {
      fixedness_set = true;
      fixedness     = false;
    }
    else
    {
      errors << priority(error) << "Invalid value '" << fixedness_string << "' given for FIXED attribute." << endl;
    }
  }
  else
  {
    errors << priority(error) << "No value specified for FIXED attribute. Ignoring..." << endl;
  }
}








//=============================================================
//
//  parse(...)
//
//=============================================================
bool
NewTypeSpecificParser::parse(NewTypeSpecificSubmodel & submodel, const LSFBase * param, cerrorstream & errors)
{
  // Instantiate temporary submodel:

  MAXFUN::NewSubmodelShPtr   temp_submodel_shptr = submodel.clone();
  NewTypeSpecificSubmodel  & temp_submodel       = *(static_cast<NewTypeSpecificSubmodel*>(temp_submodel_shptr.get()));

  // Set category:

  string category_name = toUpper(param->name().substr(5));

       if(category_name == "MEAN")    temp_submodel.setCategory(NewTypeSpecificSubmodel::MEAN);
  else if(category_name == "VAR")     temp_submodel.setCategory(NewTypeSpecificSubmodel::VARIANCE);
  else if(category_name == "SUSCEPT") temp_submodel.setCategory(NewTypeSpecificSubmodel::SUSCEPTIBILITY);

  // Loop through parameters:

  if(!param->List())
    return false;

  for(LSFList::const_iterator iter = param->List()->begin(); iter != param->List()->end(); ++iter)
  {
    if(!(*iter))
      continue;

    string param_name = toUpper((*iter)->name());

    if(param_name == "OPTION") 
    {
      parseOption(temp_submodel, *iter, errors);
    }
    else if(param_name == toUpper(temp_submodel.getBriefName()))
    {
      parseGenotype(temp_submodel, *iter, errors);
    }
    else
    {
      errors << priority(error) << "Unrecognized parameter '" << param_name << "' in " << getSubblockName(temp_submodel) << " sub-block." << endl;
    }
  }

  int warning_code = 0;

  int err_code = temp_submodel.isInitialSetupConsistent(warning_code);

  if(warning_code != 0)
  {
    string message = "";

    if(warning_code == 1)
    {
      message = "OPTION=TWO, AA and BB are given but not equal; interpreting as " + TypeConstraint(temp_submodel.getDefaultTwoTreatment()).toString();
    }
    else if(warning_code == 2)
    {
      message = "OPTION=TWO, AB and BB are given but not equal; interpreting as TWO_DOM";
    }

    errors << priority(warning)
           << "Unclear statement in " 
           << getSubblockName(temp_submodel)
           << " subblock:"
           << endl
           << message
           << endl;
  }

  if(err_code == 0)
  {
    submodel = temp_submodel;

    return true;
  }
  else
  {
    string message = "";

    if(err_code == 1)
    {
      message = "OPTION=ONE, but at least two different values given";
    }
    else if(err_code == 2)
    {
      message = "OPTION=TWO_DOM, but different values given for AA and AB";
    }
    else if(err_code == 3)
    {
      message = "OPTION=TWO_DOM, valid value given for dominant type (AA/AB), but no value given for recessive type (BB)";
    }
    else if(err_code == 4)
    {
      message = "OPTION=TWO_REC, but different values given for AB and BB";
    }
    else if(err_code == 5)
    {
      message = "OPTION=TWO_REC, valid value given for dominant type (AB/BB), but no value given for recessive type (AA)";
    }
    else if(err_code == 6)
    {
      message = "OPTION=TWO/TWO_DOM/TWO_REC, but one one value specified";
    }
    else if(err_code == 7)
    {
      message = "OPTION=THREE_ADD, three values given, but do not meet relational requirement";
    }
    else if(err_code == 8)
    {
      message = "OPTION=THREE_DEC, three values given, but do not meet relational requirement";
    }
    else if(err_code == 9)
    {
      message = "OPTION=THREE_INC, three values given, but do not meet relational requirement";
    }
    else if(err_code == 10)
    {
      message = "OPTION=THREE_ADD/THREE_DEC/THREE_INC, at least one value given but not all three";
    }
    else if(err_code == 11)
    {
      message = "OPTION=TWO, but three different values given";
    }
    else if(err_code == 12)
    {
      message = "OPTION=TWO, AA equals AB, but BB not given; not sure how to interpret this";
    }
    else if(err_code == 13)
    {
      message = "OPTION=TWO, AB equals BB, but AA not given; not sure how to interpret this";
    }
    else if(err_code == 14)
    {
      message = "OPTION=TWO, AA equals BB, and AB != AA or BB; not sure how to interpret this";
    }
    else if(err_code == 100)
    {
      message = "OPTION=ONE, at least two different fixedness statuses given";
    }
    else if(err_code == 101)
    {
      message = "OPTION=TWO (interpreted as TWO_DOM) or TWO_DOM, conflicting fixed status given for AA and AB";
    }
    else if(err_code == 102)
    {
      message = "OPTION=TWO (interpreted as TWO_REC) or TWO_REC, conflicting fixed status given for AB and BB";
    }
    else if(err_code == 200)
    {
      message = "Variance sub-block indicated, value less than zero given";
    }

    errors << priority(error)
           << "Inconsistency in " 
           << getSubblockName(temp_submodel)
           << " subblock:"
           << endl
           << message
           << endl
           << "Ignoring subblock..."
           << endl;
  }

  return false;
}

//============================================================
//
//  areMeanAndVarianceCompatible(...)
//
//============================================================
bool 
NewTypeSpecificParser::areMeanAndVarianceCompatible(
  const NewTypeSpecificSubmodel & mean_submodel, 
  const NewTypeSpecificSubmodel & var_submodel,
  cerrorstream & errors)
{
  int error_code = NewTypeSpecificSubmodel::areMeanAndVarianceCompatible(mean_submodel, var_submodel);

  if(error_code != 0)
  {
    string message = "";

    if(error_code == 1)
    {
      message = "Number of variances ("       + var_submodel  .getTypeConstraint().getTypeCountAsString() +
                ") exceeds number of means (" + mean_submodel .getTypeConstraint().getTypeCountAsString() + ")";
    }

    errors << priority(error)
           << "Inconsistency in " 
           << getSubblockName(mean_submodel)
           << " / "
           << getSubblockName(var_submodel)
           << " subblocks:"
           << endl
           << message
           << endl;

    return false;
  }

  return true;
}


//============================================================
//
//  getSubblockName(...)
//
//============================================================
string 
NewTypeSpecificParser::getSubblockName(const NewTypeSpecificSubmodel & submodel)
{
  return "type_" + submodel.getBriefName();
}

//============================================================
//
//  parseOption(...)
//
//============================================================
void 
NewTypeSpecificParser::parseOption(NewTypeSpecificSubmodel & submodel, const LSFBase * param, cerrorstream & errors)
{
  if(!param)
    return;

  if(toUpper(param->name()) != "OPTION")
    return;

  // If this is an "OPTION" parameter, pull out the option value:

  AttrVal a = attr_value(param, 0);

  if(a.has_value() == true)
  {
    // Fetch the ConstraintType, check it, and set it:

    TypeConstraint t;

    if(t.setOption(toUpper(a.String())) ? !submodel.setTypeConstraint(t) : true)
    {
      errors << priority(error) 
             << "Value '" 
             << a.String()
             << "'for"
             << " option parameter of type_mean or type_susceptibility" 
             << " sub-block not recognized or not allowed. "
             << " Ignoring ..." 
             << endl;
    }
  }
  else
  {
    errors << priority(error)
           << "There is no value given for OPTION. Ignoring..."
           << endl;
  }
}

//============================================================
//
//  parseGenotype(...)
//
//============================================================
void 
NewTypeSpecificParser::parseGenotype(NewTypeSpecificSubmodel & submodel, const LSFBase * param, cerrorstream & errors)
{
  if(!param)
    return;

  if((toUpper(param->name()) != toUpper(submodel.getBriefName())) || !param->attrs())
    return;

  // Set up variables that will/may be populated with data:

  vector<size_t> genotypes(0);
  double         initial_val     = std::numeric_limits<double>::quiet_NaN();
  bool           fixedness       = false,
                 genotype_set    = false,
                 initial_val_set = false,
                 fixedness_set   = false;

  // Loop through attribute list:

  for(AttrList::const_iterator iter = param->attrs()->begin(); iter != param->attrs()->end(); ++iter)
  {
    string attr_name = toUpper(AttrNameMgr.query(iter->first));

         if(attr_name == "VALUE") parseGenotypeValue (genotypes,   genotype_set,    iter, errors);
    else if(attr_name == "VAL")   parseGenotypeVal   (initial_val, initial_val_set, iter, errors);
    else if(attr_name == "FIXED") parseGenotypeFixed (fixedness,   fixedness_set,   iter, errors);
    else
    {
      errors << priority(error) << "Invalid value '" << attr_name << "' given as an attribute " << param->name() << endl;
    }
  }
   
  // Having looped through all the attributes, make sure we got the ones we needed:

  if(genotype_set == false)
  {
    errors << priority(error) << "No genotype given for '" << param->name() << "'" << endl;
  }
  else if(fixedness_set && !initial_val_set)
  {
    errors << priority(error) << "Fixedness set for parameter '" << param->name() << "', but val not indicated." << endl;
  }
  else if(!initial_val_set && !fixedness_set)
  {
    errors << priority(error) << "Parameter '" << param->name() << "' indicated, but no attributes given." << endl;
  }
  else
  {
    for(vector<size_t>::iterator i = genotypes.begin(); i != genotypes.end(); ++i)
    {
      submodel.setTypeValues(*i, initial_val, fixedness_set, fixedness);
    }
  }
}

//============================================================
//
//  parseGenotypeValue(...)
//
//============================================================
void 
NewTypeSpecificParser::parseGenotypeValue(
  vector<size_t>           & genotypes,
  bool                     & genotype_set,
  AttrList::const_iterator   iter, 
  cerrorstream             & errors)
{
  genotype_set = false;

  string attr_name = toUpper(AttrNameMgr.query(iter->first));

  if(attr_name != "VALUE")
    return;

  AttrVal a = iter->second;

  if(a.has_value())
  {
    string genotype_string = toUpper(a.String());

    genotype_set = true;

         if(genotype_string == "AA")   genotypes.push_back(0);
    else if(genotype_string == "AB")   genotypes.push_back(1);
    else if(genotype_string == "BB")   genotypes.push_back(2);
    else if(genotype_string == "A*") { genotypes.push_back(0); 
                                       genotypes.push_back(1); }
    else if(genotype_string == "B*") { genotypes.push_back(1); 
                                       genotypes.push_back(2); }
    else if(genotype_string == "**") { genotypes.push_back(0); 
                                       genotypes.push_back(1); 
                                       genotypes.push_back(2); }
    else
    {
      genotype_set = false;

      errors << priority(error) << "Invalid value '" << genotype_string << "' given, only 'AA', 'AB', 'BB', 'A*', 'B*', and '**' allowed." << endl;
    }
  }
  else
  {
    errors << priority(error) << "No value specified for VALUE attribute. Ignoring..." << endl;
  }
}

//============================================================
//
//  parseGenotypeVal(...)
//
//============================================================
void 
NewTypeSpecificParser::parseGenotypeVal(
  double                   & initial_val, 
  bool                     & initial_val_set, 
  AttrList::const_iterator   iter, 
  cerrorstream             & errors)
{
  initial_val_set = false;

  string attr_name = toUpper(AttrNameMgr.query(iter->first));

  if(attr_name != "VAL")
    return;

  AttrVal a = iter->second;

  if(a.has_value())
  {
    if(finite(a.Real()))
    {
      initial_val     = a.Real();
      initial_val_set = true;
    }   
    else
    {
      errors << priority(error) << "Value of " << a.String() << " not understood.  Ignoring ..." << endl;
    }
  }
  else
  {
    errors << priority(error) << "No value specified for VAL attribute. Ignoring..." << endl;
  }
}

//============================================================
//
//  parseGenotypeFixed(...)
//
//============================================================
void 
NewTypeSpecificParser::parseGenotypeFixed(
  bool                     & fixedness,
  bool                     & fixedness_set, 
  AttrList::const_iterator   iter, 
  cerrorstream             & errors)
{
  fixedness_set = false;

  string attr_name = toUpper(AttrNameMgr.query(iter->first));

  if(attr_name != "FIXED")
    return;

  AttrVal a = iter->second;

  if(a.has_value())
  {
    string fixedness_string = toUpper(a.String());

    if(fixedness_string == "TRUE")
    {
      fixedness_set = true;
      fixedness     = true;
    }
    else if(fixedness_string == "FALSE")
    {
      fixedness_set = true;
      fixedness     = false;
    }
    else
    {
      errors << priority(error) << "Invalid value '" << fixedness_string << "' given for FIXED attribute." << endl;
    }
  }
  else
  {
    errors << priority(error) << "No value specified for FIXED attribute. Ignoring..." << endl;
  }
}









} // End namespace MFSUBMODELS
} // End namespace SAGE
