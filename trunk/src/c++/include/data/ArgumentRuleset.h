#ifndef DATA_ARGUMENT_RULESET
#define DATA_ARGUMENT_RULESET

// #define DEPRECATION_MSG

#include <iostream>
#include <vector>
#include <string>
#include <map>

namespace SAGE {
namespace APP  {

/// Indicates what type of file an argument is.
enum FileType { PARAMETER_FILE, PEDIGREE_FILE, LOCUS_FILE, GENOME_FILE, TRAIT_MODEL_FILE, IBD_FILE };
    
/// The number of file types.
static const size_t file_type_count = 6;

static const std::string long_names  [file_type_count] = { "Parameter",   "Pedigree Data", "Locus Description", "Genome Description", "Trait Locus Description OR Type Probability", "IBD Sharing FIle" };
static const std::string short_names [file_type_count] = { "parameters",  "pedigree",      "locus",             "genome",             "trait model",                                 "ibd" };
static const char        flags       [file_type_count] = { 'p',           'd',             'l',                 'g',                  'm',                                           'i'   };

/// \brief Specifies which commandline arguments are optional, required, or forbidden.
/// Take a look at test_cmdline.cpp in the data directory. This is a good all-around test of the arguments system.
///
class ArgumentRuleset
{
  public:
  
  /// Indicates how many times a given argument can occur on the commandline
  enum Count { ZERO_OR_ONE, ONE, ONE_OR_MORE, ZERO_OR_MORE };

  /// A pair, where the first element is a file type, and the second is it's requirement status.
  typedef std::pair<FileType, Count> Rule;

  // A vector of Argument's.
  typedef std::vector<Rule> RuleVector;

  /// @name Constructors
  //@{
  
    ///
    /// Constructor
    ArgumentRuleset();
    
    ///
    /// Copy Constructor
    ArgumentRuleset(const ArgumentRuleset & other);
    

    ///
    /// operator=
    ArgumentRuleset & operator= (const ArgumentRuleset & other);
  
  //@}

  /// @name Setting/getting configuration
  //@{
  
    ///
    /// Returns the vector of allowable arguments.
    const RuleVector & get_rules() const { return my_rules; }

    ///
    /// Returns the vector of allowable arguments.
    RuleVector & get_rules() { return my_rules; }
    void  add_rule(const Rule& rule);
    void  clear();
    
  //@}

  /// Debugging:
  void dump() const;
  
  private:
  
    RuleVector my_rules;
    bool  first_optional_added;   
};
  

// Forward declaration:
class ArgumentParser;
  
/// \brief Specifies which commandline arguments were found, based on an ArgumentRuleset instance.
///
class ArgumentsFound
{
  public:
  
    /// Only the ArgumentParser is allowed to configure the ArgumentsFound instance
    friend class ArgumentParser;

    /// A vector of filenames.
    typedef std::vector<std::string> FilenameVector;

  /// @name Constructors
  //@{
  
    ///
    /// Constructor
    ArgumentsFound();
    
    ///
    /// Copy constructor
    ArgumentsFound(const ArgumentsFound & other);
    
    ///
    /// operator=
    ArgumentsFound& operator= (const ArgumentsFound & other);
  
  //@}

  /// @name Getting information
  //@{
  
    ///
    /// Returns whether or not at least one argument for the given FileType was specified.
    bool argument_specified(FileType f) const { return my_fvectors.find(f)->second.size(); }
    
    ///
    /// Returns the list of names given for this FileType.
    const FilenameVector & get_arguments(FileType f) const { return my_fvectors.find(f)->second; }
    
  //@}
  
  /// Debugging:
  void dump() const;

  private:
  
  /// @name Setting information
  //@{
  
    ///
    /// Sets the filename for the given argument.
    void add_argument(FileType f, const std::string & s) { my_fvectors[f].push_back(s); }
  
  //@}

  typedef std::map<FileType, FilenameVector> FilenameVectors;
    
  FilenameVectors my_fvectors;
};

/// \brief Parses a commandline, based on an ArgumentRuleset, and returns an ArgumentsFound instance.
///
class ArgumentParser
{
  public:
  
  /// @name Parsing features
  //@{
  
    typedef std::vector<std::string> ArgVector;

    ///
    /// Converts the argc/argv variables into a straightforward vector of strings.
    static ArgVector convert_args(int argc, char ** argv);

    ///
    /// Displays information about how to correctly invoke the program from the commandline.
    /// \param program_name The name of the program
    /// \param ruleset The rules for which arguments can occur on the commandline
    /// \param out The iostream to which the display will be sent
    static void display_usage(const std::string & program_name, const ArgumentRuleset & ruleset, std::ostream & out);

    ///
    /// Parses the given commandline (via argc/argv), using the rules given in the ArgumentRuleset instance.
    /// Returns an ArgumentsFound instance describing which arguments were found.
    /// NOTE: Throws an exception is the commandline is incompatible in any way with the given ruleset.
    static ArgumentsFound parse_commandline(const ArgVector & args, const ArgumentRuleset & ruleset);
  
  //@}
};

/// Processes a commandline, and returns an ArgumentsFound instance.
///
/// If the commandline is not processable for any reason, this function invokes ArgumentParser::display_usage().
/// If the exit_on_fail argument is set to true, the program aborts.
///
ArgumentsFound process_commandline(
        int                argc,
        char            ** argv,
  const std::string     &  program_name,
  const ArgumentRuleset &  ruleset,
        std::ostream    &  out,
        bool               exit_on_fail);
  


} // End namespace APP
} // End namespace SAGE

#endif

