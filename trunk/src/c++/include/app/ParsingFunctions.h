#ifndef APP_PARSING_FUNCTIONS
#define APP_PARSING_FUNCTIONS

#include "app/aparser.h"

namespace SAGE {
namespace APP  {

class ParsingFunctions
{
  private:
  
    class SimpleParser : public APP::BasicParser
    {
      friend class ParsingFunctions;
      
      void parse_symbols(const SymbolTable* syms) {}
      void parse_parameter(const LSFBase* param) {}
      void parse_test_parameter_section(const LSFBase* params) {}
      void parse_test_parameter(const LSFBase* param) {}
    };

  public:
  
  /// @name Default attribute parsing functions
  //@{

    ///
    /// Extracts a boolean from an LSFBase pointer's default attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param value The boolean reference into which the value will be placed.
    /// \param verbose Indicates whether or not error messages should be generated.
    /// \returns An BasicParser::error_t indicating the success/failure of the extraction.
    static LSFConvert::error_t parse_boolean (const LSFBase* param, bool&   value, bool verbose = true)
    {
      return SimpleParser().parse_boolean(param, value, verbose);
    }

    ///
    /// Extracts a int from an LSFBase pointer's default attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param value The int reference into which the value will be placed.
    /// \param verbose Indicates whether or not error messages should be generated.
    /// \returns An BasicParser::error_t indicating the success/failure of the extraction.
    static LSFConvert::error_t parse_integer (const LSFBase* param, int&    value, bool verbose = true)
    {
      return SimpleParser().parse_integer(param, value, verbose);
    }

    ///
    /// Extracts a string from an LSFBase pointer's default attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param value The string reference into which the value will be placed.
    /// \param verbose Indicates whether or not error messages should be generated.
    /// \returns An BasicParser::error_t indicating the success/failure of the extraction.
    static LSFConvert::error_t parse_string  (const LSFBase* param, string& value, bool verbose = true)
    {
      return SimpleParser().parse_string(param, value, verbose);
    }

    ///
    /// Extracts a double from an LSFBase pointer's default attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param value The double reference into which the value will be placed.
    /// \param verbose Indicates whether or not error messages should be generated.
    /// \returns An BasicParser::error_t indicating the success/failure of the extraction.
    static LSFConvert::error_t parse_real    (const LSFBase* param, double& value, bool verbose = true)
    {
      return SimpleParser().parse_real(param, value, verbose);
    }

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
    static LSFConvert::error_t parse_boolean (const LSFBase* param, const string& attr, bool&   value, bool verbose = true)
    {
      return SimpleParser().parse_boolean(param, attr, value, verbose);
    }

    ///
    /// Extracts a int from an LSFBase pointer's named attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param attr The name of the attribute whose value will be extracted.
    /// \param value The int reference into which the value will be placed.
    /// \param verbose Indicates whether or not error messages should be generated.
    /// \returns An BasicParser::error_t indicating the success/failure of the extraction.
    static LSFConvert::error_t parse_integer (const LSFBase* param, const string& attr, int&    value, bool verbose = true)
    {
      return SimpleParser().parse_integer(param, attr, value, verbose);
    }
              
    ///
    /// Extracts a string from an LSFBase pointer's named attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param attr The name of the attribute whose value will be extracted.
    /// \param value The string reference into which the value will be placed.
    /// \param verbose Indicates whether or not error messages should be generated.
    /// \returns An BasicParser::error_t indicating the success/failure of the extraction.
    static LSFConvert::error_t parse_string  (const LSFBase* param, const string& attr, string& value, bool verbose = true)
    {
      return SimpleParser().parse_string(param, attr, value, verbose);
    }
              
    ///
    /// Extracts a double from an LSFBase pointer's named attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param attr The name of the attribute whose value will be extracted.
    /// \param value The double reference into which the value will be placed.
    /// \param verbose Indicates whether or not error messages should be generated.
    /// \returns An BasicParser::error_t indicating the success/failure of the extraction.
    static LSFConvert::error_t parse_real    (const LSFBase* param, const string& attr, double& value, bool verbose = true)
    {
      return SimpleParser().parse_real(param, attr, value, verbose);
    }

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
    static LSFConvert::error_t parse_boolean(const LSFBase* param, attr_id id, bool&   value, bool verbose = true)
    {
      return SimpleParser().parse_boolean(param, id, value, verbose);
    }

    ///
    /// Extracts a int from an LSFBase pointer's numbered attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param id The id of the attribute whose value will be extracted.
    /// \param value The int reference into which the value will be placed.
    /// \param verbose Indicates whether or not error messages should be generated.
    /// \returns An BasicParser::error_t indicating the success/failure of the extraction.
    static LSFConvert::error_t parse_integer(const LSFBase* param, attr_id id, int&    value, bool verbose = true)
    {
      return SimpleParser().parse_integer(param, id, value, verbose);
    }

    ///
    /// Extracts a string from an LSFBase pointer's numbered attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param id The id of the attribute whose value will be extracted.
    /// \param value The string reference into which the value will be placed.
    /// \param verbose Indicates whether or not error messages should be generated.
    /// \returns An BasicParser::error_t indicating the success/failure of the extraction.
    static LSFConvert::error_t  parse_string(const LSFBase* param, attr_id id, string& value, bool verbose = true)
    {
      return SimpleParser().parse_string(param, id, value, verbose);
    }

    ///
    /// Extracts a double from an LSFBase pointer's numbered attribute.
    /// \param param The LSFBase pointer from which the value will be extracted.
    /// \param id The id of the attribute whose value will be extracted.
    /// \param value The double reference into which the value will be placed.
    /// \param verbose Indicates whether or not error messages should be generated.
    /// \returns An BasicParser::error_t indicating the success/failure of the extraction.
    static LSFConvert::error_t    parse_real(const LSFBase* param, attr_id id, double& value, bool verbose = true)
    {
      return SimpleParser().parse_real(param, id, value, verbose);
    }

  //@}
};

} // End namespace APP
} // End namespace SAGE

#endif
