#ifndef MFSUBMODELS_TYPEPARSER_H
#define MFSUBMODELS_TYPEPARSER_H

#include "LSF/LSF.h"
#include "mfsubmodels/TypeSubmodel.h"

namespace SAGE        {
namespace MFSUBMODELS {

/** \brief Parses a TypeSpecificSubmodel
  *
  */
class TypeSpecificParser 
{
public:

  ///
  /// Parses the given submodel.
  /// \retval true Parsing was successful
  /// \retval false Parsing was not successful
  static bool parse(TypeSpecificSubmodel & submodel, const LSFBase * param, cerrorstream & errors);

  ///
  /// If your analysis uses both MEAN and VARIANCE submodels, you must make sure that the settings assigned
  /// in both submodels are compatible with each other. In particular, you need to make sure that
  /// the number of variances does not exceed the number of means. This function carries out such an evaluation,
  /// and reports any inconsistencies.
  /// \param mean_submodel The mean submodel
  /// \param var_submodel The variance submodel
  /// \param errors The errorstream to use
  /// \retval true The submodels are compatible
  /// \retval false The submodels are incompatible
  static bool areMeanAndVarianceCompatible(const TypeSpecificSubmodel & mean_submodel, const TypeSpecificSubmodel & var_submodel, cerrorstream & errors);

private:

  ///
  /// Returns the name of the parameter sub-block for the given submodel.
  /// (Ie: "type_suscept")
  static string getSubblockName(const TypeSpecificSubmodel & submodel);

  ///
  /// Parses the "OPTION" parameter for a submodel.
  static void parseOption(TypeSpecificSubmodel & submodel, const LSFBase * param, cerrorstream & errors);

  ///
  /// Parses the "MEAN", "VAR", or "SUSCEPT" parameter for a submodel.
  static void parseGenotype(TypeSpecificSubmodel & submodel, const LSFBase * param, cerrorstream & errors);

  ///
  /// Parses the "VALUE" attribute of a MEAN/VAR/SUSCEPT parameter.
  static void parseGenotypeValue(vector<size_t> & genotypes, bool & genotype_set, AttrList::const_iterator iter, cerrorstream & errors);

  ///
  /// Parses the "VAL" attribute of a MEAN/VAR/SUSCEPT parameter.
  static void parseGenotypeVal(double & initial_val, bool & initial_val_set, AttrList::const_iterator iter, cerrorstream & errors);

  ///
  /// Parses the "FIXED" attribute of a MEAN/VAR/SUSCEPT parameter.
  static void parseGenotypeFixed(bool & fixedness, bool & fixedness_set, AttrList::const_iterator iter, cerrorstream & errors);
};




/** \brief Parses a NewTypeSpecificSubmodel
  *
  */
class NewTypeSpecificParser 
{
public:

  ///
  /// Parses the given submodel.
  /// \retval true Parsing was successful
  /// \retval false Parsing was not successful
  static bool parse(NewTypeSpecificSubmodel & submodel, const LSFBase * param, cerrorstream & errors);

  ///
  /// If your analysis uses both MEAN and VARIANCE submodels, you must make sure that the settings assigned
  /// in both submodels are compatible with each other. In particular, you need to make sure that
  /// the number of variances does not exceed the number of means. This function carries out such an evaluation,
  /// and reports any inconsistencies.
  /// \param mean_submodel The mean submodel
  /// \param var_submodel The variance submodel
  /// \param errors The errorstream to use
  /// \retval true The submodels are compatible
  /// \retval false The submodels are incompatible
  static bool areMeanAndVarianceCompatible(const NewTypeSpecificSubmodel & mean_submodel, const NewTypeSpecificSubmodel & var_submodel, cerrorstream & errors);

private:

  ///
  /// Returns the name of the parameter sub-block for the given submodel.
  /// (Ie: "type_suscept")
  static string getSubblockName(const NewTypeSpecificSubmodel & submodel);

  ///
  /// Parses the "OPTION" parameter for a submodel.
  static void parseOption(NewTypeSpecificSubmodel & submodel, const LSFBase * param, cerrorstream & errors);

  ///
  /// Parses the "MEAN", "VAR", or "SUSCEPT" parameter for a submodel.
  static void parseGenotype(NewTypeSpecificSubmodel & submodel, const LSFBase * param, cerrorstream & errors);

  ///
  /// Parses the "VALUE" attribute of a MEAN/VAR/SUSCEPT parameter.
  static void parseGenotypeValue(vector<size_t> & genotypes, bool & genotype_set, AttrList::const_iterator iter, cerrorstream & errors);

  ///
  /// Parses the "VAL" attribute of a MEAN/VAR/SUSCEPT parameter.
  static void parseGenotypeVal(double & initial_val, bool & initial_val_set, AttrList::const_iterator iter, cerrorstream & errors);

  ///
  /// Parses the "FIXED" attribute of a MEAN/VAR/SUSCEPT parameter.
  static void parseGenotypeFixed(bool & fixedness, bool & fixedness_set, AttrList::const_iterator iter, cerrorstream & errors);
};




//=====================
//  INLINE FUNCTIONS
//=====================


} // End namespace MFSUBMODELS
} // End namespace SAGE

#endif
