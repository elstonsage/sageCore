#ifndef TDTEX_PARSER_H
#define TDTEX_PARSER_H

#include "LSF/LSF.h"
#include "app/aparser.h"
#include "error/errorstream.h"
#include "LSF/LSFsymbol.h"
#include "tdtex/TransmissionTables.h"
#include "tdtex/Configuration.h"

namespace SAGE  {
namespace TDTEX {

/** \class Parser
  * \brief tdtex-specialized parser derived from APP::BasicParser
  */
class Parser : public APP::BasicParser
{
public:

  ///
  /// Parses the given lsf params into a Configuration instance.
  /// Throws a std::exception if the Configuration is invalid for any reason.
  static Configuration parse_parameters(const RPED::RefMPedInfo & mped_info, const LSFBase * params, cerrorstream & err = SAGE::sage_cerr);

private:

  /// @name Constructors
  //@{

    ///
    /// Constructor
    Parser(const RPED::RefMPedInfo & mped_info, cerrorstream & err = SAGE::sage_cerr);

  //@}

    ///
    /// Returns the parsed Configuration.
    const Configuration & get_configuration() const { return my_config; }

    ///
    /// Returns the parsed Configuration.
    Configuration & get_configuration() { return my_config; }

  /// @name Parsing functions
  //@{

    // Required by virtual table, but not used:
    virtual void parse_symbols(const SymbolTable* syms) { }
    virtual void parse_parameter(const LSFBase* param)  { }

    // Loops through parameters and parses them:
    virtual void parse_test_parameter_section(const LSFBase* params);
    virtual void parse_test_parameter(const LSFBase* param);

  //@}

  ///@name Parsing functions
  //@{

    virtual void parse_limit(const LSFBase* param, attr_id id, size_t& value, const std::string& name);

    void parse_limit(const LSFBase* param, size_t& value, const std::string& name);

    void parse_limit(const LSFBase* param, const std::string& attr, size_t& value, const std::string& name);

    void parse_marker(const LSFBase* param);

    void parse_trait(const LSFBase* param);

    void parse_parent_trait(const LSFBase* param);

    void parse_sample(const LSFBase* param);

    void parse_parent_sex(const LSFBase* param);

  //@}

    void set_sample_alleles();

    void set_sample_genotypes();

  // Data members

    Configuration       my_config;

    const RPED::RefMPedInfo & my_mped_info;
};


} // End namespace TDTEX
} // End namespace SAGE


#endif
