#ifndef PARSER_H
#define PARSER_H

#include <limits>
#include <iostream>
#include <vector>
#include <string>
#include "error/errorstream.h"
#include "error/errorbuf.h"
#include "error/errormanip.h"
#include "LSF/LSFsymbol.h"

namespace SAGE {
namespace APP  {

class BasicParser;

/** \brief Converts LSF attribute data into basic data types.
  *
  * This class provides functionality for converting LSF data into specific 
  * data types, doing
  * validity checking as they do.  They return an error flag to indicate
  * problems.  For the most part, they are self explanatory, following a
  * common format.
  *
  * There are two versions of the functions.  The first version is quiet,
  * while the latter takes an errorstream as an argument and produces messages
  * for some (though not all) of the potential errors.
  */ 
class LSFConvert
{
public:

  ///
  /// Errors that the conversion functions may return
  enum error_t
  {
    GOOD,             ///< Parse was successful
    LSF_NULL,         ///< LSF object empty                                (debug msg)
    NO_LIST,          ///< LSF object has no attributes                    (quiet)  
    NO_ATTR,          ///< LSF object does not have XXX attribute          (quiet)
    ATTR_EMPTY,       ///< The attribute exists but has no value           (warning)
    ATTR_INVALID,     ///< The attribute is not understood (type specific) (error)
    INVALID           ///< Invalid error, should not be used               (no msg)
  };

  /// @name Default attribute extraction functions
  //@{
    
    ///
    /// Extracts a boolean from an LSFBase pointer's default attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param value The boolean reference into which the value will be placed.
    /// \returns An LSFConvert::error_t indicating the success/failure of the extraction.
    static error_t to_boolean (const LSFBase* param, bool&   value);

    ///
    /// Extracts a boolean from an LSFBase pointer's default attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param value The boolean reference into which the value will be placed.
    /// \param errors The errorstream to which error message should be directed.
    /// \returns An LSFConvert::error_t indicating the success/failure of the extraction.
    static error_t to_boolean (const LSFBase* param, bool&   value, cerrorstream& errors);

    ///
    /// Extracts an int from an LSFBase pointer's default attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param value The int reference into which the value will be placed.
    /// \returns An LSFConvert::error_t indicating the success/failure of the extraction.
    static error_t to_integer (const LSFBase* param, int&    value);

    ///
    /// Extracts an int from an LSFBase pointer's default attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param value The int reference into which the value will be placed.
    /// \param errors The errorstream to which error message should be directed.
    /// \returns An LSFConvert::error_t indicating the success/failure of the extraction.
    static error_t to_integer (const LSFBase* param, int&    value, cerrorstream& errors);

    ///
    /// Extracts a string from an LSFBase pointer's default attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param value The string reference into which the value will be placed.
    /// \returns An LSFConvert::error_t indicating the success/failure of the extraction.
    static error_t to_string  (const LSFBase* param, string& value);

    ///
    /// Extracts a string from an LSFBase pointer's default attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param value The string reference into which the value will be placed.
    /// \param errors The errorstream to which error message should be directed.
    /// \returns An LSFConvert::error_t indicating the success/failure of the extraction.
    static error_t to_string  (const LSFBase* param, string& value, cerrorstream& errors);

    ///
    /// Extracts a double from an LSFBase pointer's default attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param value The double reference into which the value will be placed.
    /// \returns An LSFConvert::error_t indicating the success/failure of the extraction.
    static error_t to_real    (const LSFBase* param, double& value);

    ///
    /// Extracts a double from an LSFBase pointer's default attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param value The double reference into which the value will be placed.
    /// \param errors The errorstream to which error message should be directed.
    /// \returns An LSFConvert::error_t indicating the success/failure of the extraction.
    static error_t to_real    (const LSFBase* param, double& value, cerrorstream& errors);

  //@}

  /// @name Named attribute extraction functions
  //@{

    ///
    /// Extracts a boolean from an LSFBase pointer's named attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param attr The name of the attribute whose value will be extracted.
    /// \param value The boolean reference into which the value will be placed.
    /// \returns An LSFConvert::error_t indicating the success/failure of the extraction.
    static error_t to_boolean (const LSFBase* param, const string& attr, bool&   value);

    ///
    /// Extracts a boolean from an LSFBase pointer's named attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param attr The name of the attribute whose value will be extracted.
    /// \param value The boolean reference into which the value will be placed.
    /// \param errors The errorstream to which error message should be directed.
    /// \returns An LSFConvert::error_t indicating the success/failure of the extraction.
    static error_t to_boolean (const LSFBase* param, const string& attr, bool&   value, cerrorstream& errors);

    ///
    /// Extracts an int from an LSFBase pointer's named attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param attr The name of the attribute whose value will be extracted.
    /// \param value The int reference into which the value will be placed.
    /// \returns An LSFConvert::error_t indicating the success/failure of the extraction.
    static error_t to_integer (const LSFBase* param, const string& attr, int&    value);

    ///
    /// Extracts an int from an LSFBase pointer's named attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param attr The name of the attribute whose value will be extracted.
    /// \param value The int reference into which the value will be placed.
    /// \param errors The errorstream to which error message should be directed.
    /// \returns An LSFConvert::error_t indicating the success/failure of the extraction.
    static error_t to_integer (const LSFBase* param, const string& attr, int&    value, cerrorstream& errors);

    ///
    /// Extracts a string from an LSFBase pointer's named attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param attr The name of the attribute whose value will be extracted.
    /// \param value The string reference into which the value will be placed.
    /// \returns An LSFConvert::error_t indicating the success/failure of the extraction.
    static error_t to_string  (const LSFBase* param, const string& attr, string& value);

    ///
    /// Extracts a string from an LSFBase pointer's named attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param attr The name of the attribute whose value will be extracted.
    /// \param value The string reference into which the value will be placed.
    /// \param errors The errorstream to which error message should be directed.
    /// \returns An LSFConvert::error_t indicating the success/failure of the extraction.
    static error_t to_string  (const LSFBase* param, const string& attr, string& value, cerrorstream& errors);

    ///
    /// Extracts a float from an LSFBase pointer's named attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param attr The name of the attribute whose value will be extracted.
    /// \param value The double reference into which the value will be placed.
    /// \returns An LSFConvert::error_t indicating the success/failure of the extraction.
    static error_t to_real    (const LSFBase* param, const string& attr, double& value);

    ///
    /// Extracts a float from an LSFBase pointer's named attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param attr The name of the attribute whose value will be extracted.
    /// \param value The double reference into which the value will be placed.
    /// \param errors The errorstream to which error message should be directed.
    /// \returns An LSFConvert::error_t indicating the success/failure of the extraction.
    static error_t to_real    (const LSFBase* param, const string& attr, double& value, cerrorstream& errors);

  //@}

  /// @name Numbered attribute extraction functions
  //@{

    ///
    /// Extracts a boolean from an LSFBase pointer's numbered attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param id The attr_id of the attribute whose value will be extracted.
    /// \param value The boolean reference into which the value will be placed.
    /// \returns An LSFConvert::error_t indicating the success/failure of the extraction.
    static error_t to_boolean (const LSFBase* param, attr_id id, bool&   value);

    ///
    /// Extracts a boolean from an LSFBase pointer's numbered attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param id The attr_id of the attribute whose value will be extracted.
    /// \param value The boolean reference into which the value will be placed.
    /// \param errors The errorstream to which error message should be directed.
    /// \returns An LSFConvert::error_t indicating the success/failure of the extraction.
    static error_t to_boolean (const LSFBase* param, attr_id id, bool&   value, cerrorstream& errors);

    ///
    /// Extracts an int from an LSFBase pointer's numbered attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param id The attr_id of the attribute whose value will be extracted.
    /// \param value The int reference into which the value will be placed.
    /// \returns An LSFConvert::error_t indicating the success/failure of the extraction.
    static error_t to_integer (const LSFBase* param, attr_id id, int&    value);

    ///
    /// Extracts an int from an LSFBase pointer's numbered attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param id The attr_id of the attribute whose value will be extracted.
    /// \param value The int reference into which the value will be placed.
    /// \param errors The errorstream to which error message should be directed.
    /// \returns An LSFConvert::error_t indicating the success/failure of the extraction.
    static error_t to_integer (const LSFBase* param, attr_id id, int&    value, cerrorstream& errors);

    ///
    /// Extracts a string from an LSFBase pointer's numbered attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param id The attr_id of the attribute whose value will be extracted.
    /// \param value The string reference into which the value will be placed.
    /// \returns An LSFConvert::error_t indicating the success/failure of the extraction.
    static error_t to_string  (const LSFBase* param, attr_id id, string& value);

    ///
    /// Extracts a string from an LSFBase pointer's numbered attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param id The attr_id of the attribute whose value will be extracted.
    /// \param value The string reference into which the value will be placed.
    /// \param errors The errorstream to which error message should be directed.
    /// \returns An LSFConvert::error_t indicating the success/failure of the extraction.
    static error_t to_string  (const LSFBase* param, attr_id id, string& value, cerrorstream& errors);

    ///
    /// Extracts a float from an LSFBase pointer's numbered attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param id The attr_id of the attribute whose value will be extracted.
    /// \param value The double reference into which the value will be placed.
    /// \returns An LSFConvert::error_t indicating the success/failure of the extraction.
    static error_t to_real    (const LSFBase* param, attr_id id, double& value);

    ///
    /// Extracts a float from an LSFBase pointer's numbered attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param id The attr_id of the attribute whose value will be extracted.
    /// \param value The double reference into which the value will be placed.
    /// \param errors The errorstream to which error message should be directed.
    /// \returns An LSFConvert::error_t indicating the success/failure of the extraction.
    static error_t to_real    (const LSFBase* param, attr_id id, double& value, cerrorstream& errors);

  //@}
};

///
/// Returns the parameter/attribute pair as a string for error messages.
/// Note that if id is 0, then only the parameter's name will be returned.
/// \param param The LSFBase pointer whose name will be extracted
/// \param id The attr_id identifying which attribute's name will be extracted
/// \param cap If set to \c true, the first letter of the returned string will be capitalized.
string parameter_title(const LSFBase* param, attr_id id, bool cap = false);

///
/// Checks for symbols that don't match the expected symbol list.
/// Often, we have a list of attributes to check.  We want to know if
/// there's any attributes we don't recognize.  This function does this,
/// producing error messages for each unrecognized symbol.  It returns a
/// count of the number of symbols it didn't recognize.
///
/// This function is templatized on T, the list of symbols to check.  This
/// structure must support the [] operator.
/// \param param The LSFBase pointer whose data will be examined
/// \param exp_sym_list An indexed structure of strings that represents all valid attribute names
/// \param exp_count The size of the exp_sym_list parameter
/// \param errors The errorstream to which error messages should be directed
/// \returns The number of errors detected in the attribute list
template <class T>
size_t check_for_incorrect_attributes
    (LSFBase*      param,
     T             exp_sym_list,
     size_t        exp_count,
     cerrorstream& errors);

/** \brief Virtual base class for all programs' parsers
  *
  * For every SAGE program that implements a parser, said parser should be
  * derived from this base class.
  */
class BasicParser
{
public:
  /// @name Constructor
  //@{

    ///
    /// Constructor
    /// \param err The errorstream to which error messages will be directed
    inline BasicParser(cerrorstream& err = SAGE::sage_cerr);
  
  //@}

  /// @name Required virtual interface
  //@{

    ///
    /// Find and parse symbols that are in a symbol table.
    /// \param syms The symbol table
    virtual void parse_symbols(const SymbolTable* syms) = 0;

    ///
    /// Parse a parameter *external* to the analysis block.  These should be used
    /// sparingly.
    virtual void parse_parameter(const LSFBase* param)  = 0;

    ///
    /// Parse a set of parameters contained in the LSFlist contained in
    /// params.
    virtual void parse_test_parameter_section(const LSFBase* params) = 0;

    ///
    /// Parse a single parameter.  This is generally called from
    /// parse_test_parameter_section.
    virtual void parse_test_parameter(const LSFBase* param) = 0;

  //@}

  /// @name Optional virtual interface
  //@{

    ///
    /// Destructor
    virtual inline ~BasicParser();

    ///
    /// Clear
    virtual inline void clear();

  //@}

  /// @name Validity funtions
  //@{

    ///
    /// Returns a bool indicating whether or not the thingy is valid.
    virtual inline bool valid() const;

    ///
    /// Sets the thingy to invalid.
    virtual inline void invalidate();

  //@}

protected:

  ///
  /// Error code from LSFConvert
  typedef LSFConvert::error_t error_t;

  /// @name Basic Parsing Functions
  /// These are useful overrides to the conversion functions which take a
  /// boolean verbosity flag instead of the ostream.  If verbose (by
  /// default), they use the BasicParser's errorstream.
  //@{
  //@}

  /// @name Default attribute parsing functions
  //@{

    ///
    /// Extracts a boolean from an LSFBase pointer's default attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param value The boolean reference into which the value will be placed.
    /// \param verbose Indicates whether or not error messages should be generated.
    /// \returns An BasicParser::error_t indicating the success/failure of the extraction.
    error_t parse_boolean (const LSFBase* param, bool&   value, bool verbose = true);

    ///
    /// Extracts a int from an LSFBase pointer's default attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param value The int reference into which the value will be placed.
    /// \param verbose Indicates whether or not error messages should be generated.
    /// \returns An BasicParser::error_t indicating the success/failure of the extraction.
    error_t parse_integer (const LSFBase* param, int&    value, bool verbose = true);

    ///
    /// Extracts a string from an LSFBase pointer's default attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param value The string reference into which the value will be placed.
    /// \param verbose Indicates whether or not error messages should be generated.
    /// \returns An BasicParser::error_t indicating the success/failure of the extraction.
    error_t parse_string  (const LSFBase* param, string& value, bool verbose = true);

    ///
    /// Extracts a double from an LSFBase pointer's default attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param value The double reference into which the value will be placed.
    /// \param verbose Indicates whether or not error messages should be generated.
    /// \returns An BasicParser::error_t indicating the success/failure of the extraction.
    error_t parse_real    (const LSFBase* param, double& value, bool verbose = true);

  //@}

  /// @name Named attribute parsing functions
  //@{

    ///
    /// Extracts a boolean from an LSFBase pointer's named attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param attr The name of the attribute whose value will be extracted.
    /// \param value The boolean reference into which the value will be placed.
    /// \param verbose Indicates whether or not error messages should be generated.
    /// \returns An BasicParser::error_t indicating the success/failure of the extraction.
    error_t parse_boolean (const LSFBase* param, const string& attr, bool&   value, bool verbose = true);

    ///
    /// Extracts a int from an LSFBase pointer's named attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param attr The name of the attribute whose value will be extracted.
    /// \param value The int reference into which the value will be placed.
    /// \param verbose Indicates whether or not error messages should be generated.
    /// \returns An BasicParser::error_t indicating the success/failure of the extraction.
    error_t parse_integer (const LSFBase* param, const string& attr, int&    value, bool verbose = true);

    ///
    /// Extracts a string from an LSFBase pointer's named attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param attr The name of the attribute whose value will be extracted.
    /// \param value The string reference into which the value will be placed.
    /// \param verbose Indicates whether or not error messages should be generated.
    /// \returns An BasicParser::error_t indicating the success/failure of the extraction.
    error_t parse_string  (const LSFBase* param, const string& attr, string& value, bool verbose = true);

    ///
    /// Extracts a double from an LSFBase pointer's named attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param attr The name of the attribute whose value will be extracted.
    /// \param value The double reference into which the value will be placed.
    /// \param verbose Indicates whether or not error messages should be generated.
    /// \returns An BasicParser::error_t indicating the success/failure of the extraction.
    error_t parse_real    (const LSFBase* param, const string& attr, double& value, bool verbose = true);

  //@}

  /// @name Numbered attribute parsing functions
  //@{

    ///
    /// Extracts a boolean from an LSFBase pointer's numbered attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param id The id of the attribute whose value will be extracted.
    /// \param value The boolean reference into which the value will be placed.
    /// \param verbose Indicates whether or not error messages should be generated.
    /// \returns An BasicParser::error_t indicating the success/failure of the extraction.
    error_t parse_boolean(const LSFBase* param, attr_id id, bool&   value, bool verbose = true);

    ///
    /// Extracts a int from an LSFBase pointer's numbered attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param id The id of the attribute whose value will be extracted.
    /// \param value The int reference into which the value will be placed.
    /// \param verbose Indicates whether or not error messages should be generated.
    /// \returns An BasicParser::error_t indicating the success/failure of the extraction.
    error_t parse_integer(const LSFBase* param, attr_id id, int&    value, bool verbose = true);

    ///
    /// Extracts a string from an LSFBase pointer's numbered attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param id The id of the attribute whose value will be extracted.
    /// \param value The string reference into which the value will be placed.
    /// \param verbose Indicates whether or not error messages should be generated.
    /// \returns An BasicParser::error_t indicating the success/failure of the extraction.
    error_t  parse_string(const LSFBase* param, attr_id id, string& value, bool verbose = true);

    ///
    /// Extracts a double from an LSFBase pointer's numbered attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param id The id of the attribute whose value will be extracted.
    /// \param value The double reference into which the value will be placed.
    /// \param verbose Indicates whether or not error messages should be generated.
    /// \returns An BasicParser::error_t indicating the success/failure of the extraction.
    error_t    parse_real(const LSFBase* param, attr_id id, double& value, bool verbose = true);

  //@}

  ///
  /// Errorstream to which error messages should be directed.
  cerrorstream& errors;

private:
  bool              my_valid; 
};

//==================================================================================
//  INLINE FUNCTIONS
//==================================================================================

inline BasicParser::BasicParser(cerrorstream& err) : errors(err)
{
  clear();
}                              
  
inline BasicParser::~BasicParser() 
{ 
  clear(); 
}

inline void 
BasicParser::clear()
{
  invalidate();
}
   
inline bool 
BasicParser::valid() const 
{ 
  return my_valid;  
}

inline void 
BasicParser::invalidate()
{ 
  my_valid = false; 
}

} // End namespace APP
} // End namespace SAGE

#include "app/aparser.ipp"

#endif
