#include "data/ArgumentRuleset.h"

#include <sstream>
#include <string>
#include <iostream>
#include <iomanip>

#include "output/Output.h"


namespace SAGE {
namespace APP  {


//=============================================================================
//=============================================================================
//
//                    ArgumentRuleset
//
//=============================================================================
//=============================================================================

ArgumentRuleset::ArgumentRuleset()
      : first_optional_added(false)
{}

ArgumentRuleset::ArgumentRuleset(const ArgumentRuleset & other) : my_rules(other.my_rules) { }
    
ArgumentRuleset & 
ArgumentRuleset::operator= (const ArgumentRuleset & other)
{
  if(this != &other)
  {
    my_rules = other.my_rules;
  }
  
  return *this;
}

void
ArgumentRuleset::clear()
{
  my_rules.clear();
  first_optional_added = false;
}

void
ArgumentRuleset::add_rule(const Rule& rule)
{
  Count  c = rule.second;
  
  /*                             subsequent
   *                        < >     ...    [ ]
   *             only  < >   A       A      A
   *   previous  any   ...   NA      NA     AF 
   *             any   [ ]   NA      NA     AF
   *
   *     < > ONE
   *     ... ONE_OR_MORE
   *     [ ] ZERO_OR_ONE
   *
   *     NA - Not allowed.  To insure backward compatability, a required element may not follow
   *          an optional one.  ONE_OR_MORE consists of a required element followed by optional
   *          ones.
   *     AF - Allowed, but only parsable with flags.  Not NA, but containing two or more optional
   *          elements.
   *     A  - Everything else.
   *                               - djb
   */
  
  assert(! (first_optional_added && (c == ONE || c == ONE_OR_MORE)));
  
  my_rules.push_back(rule);
  
  if(c == ZERO_OR_ONE || c == ONE_OR_MORE || c == ZERO_OR_MORE)
  {
    first_optional_added = true;
  }
}

void
ArgumentRuleset::dump() const
{
  OUTPUT::Table t("ArgumentRuleset");
  
  t << OUTPUT::TableColumn("File type") << OUTPUT::TableColumn("Count");
  
  for(RuleVector::const_iterator rule_itr = my_rules.begin(); rule_itr != my_rules.end(); ++rule_itr)
  {
    OUTPUT::TableRow r = OUTPUT::TableRow() << long_names[(size_t)rule_itr->first];
    
    switch(rule_itr->second)
    {
      case ZERO_OR_ONE  : r << "Zero or one"; break;
      case ONE          : r << "One";         break;
      case ONE_OR_MORE  : r << "One or more"; break;
      case ZERO_OR_MORE : r << "Zero or more"; break;
    }
    
    t << r;
  }
  
  std::cout << t << std::flush;
}

//=============================================================================
//=============================================================================
//
//                    ArgumentsFound
//
//=============================================================================
//=============================================================================

ArgumentsFound::ArgumentsFound()
{
  for(size_t t = 0; t < file_type_count; ++t)
  {
    my_fvectors[static_cast<FileType>(t)] = FilenameVector(0);
  }
}
       
ArgumentsFound::ArgumentsFound(const ArgumentsFound & other) : my_fvectors(other.my_fvectors) { }
    
ArgumentsFound & 
ArgumentsFound::operator= (const ArgumentsFound & other)
{
  if(this != &other)
  {
    my_fvectors = other.my_fvectors;
  }
  
  return *this;
}
  
void
ArgumentsFound::dump() const
{
  OUTPUT::Table t("ArgumentsFound");
  
  t << OUTPUT::TableColumn("File type");
  
  for(FilenameVectors::const_iterator f = my_fvectors.begin(); f != my_fvectors.end(); ++f)
  {
    if(f->second.size())
    {
      OUTPUT::TableRow r;
    
      r << long_names[(size_t)f->first];
    
      for(FilenameVector::const_iterator p = f->second.begin(); p != f->second.end(); ++p)
        r << *p;
      
      t << r;
    }
  }
  
  std::cout << t << std::flush;
}

//=============================================================================
//=============================================================================
//
//                    ArgumentParser
//
//=============================================================================
//=============================================================================

ArgumentParser::ArgVector
ArgumentParser::convert_args(int argc, char ** argv)
{
  ArgVector args;
  
  for(size_t i = 1; i < (size_t)argc; ++i)
  {
    args.push_back(std::string(argv[i]));
  }
    
  return args;
}

void 
ArgumentParser::display_usage(const std::string & program_name, const ArgumentRuleset & ruleset, std::ostream & out)
{
  // Figure out the widest short name:
  size_t max_width = 0;

  for(size_t i = 0; i < file_type_count; ++i)
    max_width = short_names[i].length() > max_width ? short_names[i].length() : max_width;  

  // Set up both parts of the output:
  std::ostringstream noflag_str;
  std::ostringstream flagged_str;
  std::ostringstream detailed_str;

  // Figure out the no_flag string:
  noflag_str << "usage (without flags): "  << program_name;

  bool first_optional_found = false;

  for(ArgumentRuleset::RuleVector::const_iterator arg = ruleset.get_rules().begin(); arg != ruleset.get_rules().end(); ++arg)
  {
    // Fetch the info:
    FileType               ftype = arg->first;
    ArgumentRuleset::Count count = arg->second;
    switch(count)
    {
      case ArgumentRuleset::ZERO_OR_ONE:
        if(! first_optional_found)  // Must use flag method for multiple optionals.
        {
          noflag_str << " [" << short_names[(size_t)ftype] << "]";
        }
        
        first_optional_found = true;   // No way to parse an optional after this w/o flags.
        break;
        
      case ArgumentRuleset::ONE:
        noflag_str << " <" << short_names[(size_t)ftype] << ">";
        break;

      case ArgumentRuleset::ONE_OR_MORE:
        noflag_str << " <" << short_names[(size_t)ftype] << "> ...";
        first_optional_found = true;   // No way to parse an optional after this w/o flags.
        break;

      case ArgumentRuleset::ZERO_OR_MORE:
        noflag_str << " [" << short_names[(size_t)ftype] << "] ...";
        first_optional_found = true;   // No way to parse an optional after this w/o flags.
        break;
        
      default:
        assert(false);      
    }
  }

  // Figure out the flagged_str string:
  flagged_str  << "usage (with flags): "     << program_name;

  // First print requireds:
  for(ArgumentRuleset::RuleVector::const_iterator arg = ruleset.get_rules().begin(); arg != ruleset.get_rules().end(); ++arg)
    if(arg->second == ArgumentRuleset::ONE)
      flagged_str << " <-"   << flags[(size_t)arg->first] << " " << short_names[(size_t)arg->first] << ">";

  // Secondly, print the one-or-more's:
  for(ArgumentRuleset::RuleVector::const_iterator arg = ruleset.get_rules().begin(); arg != ruleset.get_rules().end(); ++arg)
    if(arg->second == ArgumentRuleset::ONE_OR_MORE)
      flagged_str << " <-"   << flags[(size_t)arg->first] << " " << short_names[(size_t)arg->first] << "> ..."; 

  // Lastly, print the zero-or-one's:
  for(ArgumentRuleset::RuleVector::const_iterator arg = ruleset.get_rules().begin(); arg != ruleset.get_rules().end(); ++arg)
    if(arg->second == ArgumentRuleset::ZERO_OR_ONE)
      flagged_str << " [-"   << flags[(size_t)arg->first] << " " << short_names[(size_t)arg->first] << "]";

  for(ArgumentRuleset::RuleVector::const_iterator arg = ruleset.get_rules().begin(); arg != ruleset.get_rules().end(); ++arg)
    if(arg->second == ArgumentRuleset::ZERO_OR_MORE)
      flagged_str << " [-"   << flags[(size_t)arg->first] << " " << short_names[(size_t)arg->first] << "] ...";

  // Figure out the detailed string:
  detailed_str << "Command line parameters:" << std::endl;
  
  for(ArgumentRuleset::RuleVector::const_iterator arg = ruleset.get_rules().begin(); arg != ruleset.get_rules().end(); ++arg)
  {
    // Fetch the info:
    FileType               ftype = arg->first;
    ArgumentRuleset::Count count = arg->second;

    // Add the info to the detailed string:    
    detailed_str << "  " << std::left << std::setw(max_width) << short_names[(size_t)ftype] << " - ";

    switch(count)
    {
      case ArgumentRuleset::ZERO_OR_ONE  : detailed_str << long_names[(size_t)ftype] << " File (optional)" << std::endl; break;
      case ArgumentRuleset::ONE          : detailed_str << long_names[(size_t)ftype] << " File"            << std::endl; break;
      case ArgumentRuleset::ONE_OR_MORE  : detailed_str << long_names[(size_t)ftype] << " File(s)"         << std::endl; break;
      case ArgumentRuleset::ZERO_OR_MORE : detailed_str << long_names[(size_t)ftype] << " (File(s))"       << std::endl; break;
    }
  }
  
  // Add the trailing info:
  noflag_str   << std::endl;
  flagged_str  << std::endl;
  detailed_str << std::endl;
  
  // Send it to the output stream:
  out << noflag_str.str() << flagged_str.str() << std::endl << detailed_str.str();
}

ArgumentsFound 
ArgumentParser::parse_commandline(const ArgVector & args, const ArgumentRuleset & ruleset)
{
  // Create the arguments found instance:
  ArgumentsFound found;
  
  // Make sure there's at least one valid argument:
  if(args.size() == 0 || (args.size() && args[0].empty()))
    throw std::exception();
    
  // Figure out if this commandline has flags or doesn't have flags:
  enum { FLAGS, NO_FLAGS } flag_status = args[0][0] == '-' ? FLAGS : NO_FLAGS;
  
  // If it's flagged, process it one way:
  if(flag_status == FLAGS)
  {
    for(size_t cur_arg = 0; cur_arg < args.size(); )
    {
      // Make sure the flag has a filename following it:
      if(cur_arg + 1 >= args.size())
        throw std::exception();
        
      // Fetch the flag and filename:
      const std::string & flag     = args[cur_arg];
      const std::string & filename = args[cur_arg + 1];
    
      // Identify the flag:
      bool flag_found = false;
      
      for(size_t i = 0; i < file_type_count; ++i)
      {
        // Does the flag match?
        if((flag.size() == 2) && (flag[0] == '-') && (flag[1] == flags[i]))
        {
          // Fetch the corresponding FileType:
          FileType ftype = (FileType)i;
          
          // Make sure the filetype is allowed:
          bool filetype_allowed = false;
          
          for(ArgumentRuleset::RuleVector::const_iterator rule_itr = ruleset.get_rules().begin(); rule_itr != ruleset.get_rules().end(); ++rule_itr)
            filetype_allowed |= (rule_itr->first == ftype);
          
          if(filetype_allowed == false)
            throw std::exception();
            
          // It's allowable; set it in the found object; set the found flag; and abort the loop:
          found.add_argument(ftype, filename);
          
          flag_found = true;
          
          break;
        }

      } // End of loop-across-flags
      
      // If the flag wasn't found, abort!
      if(flag_found == false)
        throw std::exception();

      // Increment cur_arg past the pair:
      cur_arg += 2;

    } // End loop-across-arguments
  }
  
  // Otherwise it's NOT flagged, so process it that way:
  else
  {
  
#ifdef DEPRECATION_MSG
  std::cout << "\n  *** Note: use of command line file names without file type flags is deprecated. ***\n" << std::endl;
#endif
  
    bool  first_optional_found = false;
    ArgumentRuleset::RuleVector::const_iterator rule_itr = ruleset.get_rules().begin();
    for(ArgVector::const_iterator arg_itr = args.begin(); arg_itr != args.end() && rule_itr != ruleset.get_rules().end(); ++arg_itr)
    {
      // Fetch the rule info:
      std::string            filename = *arg_itr;
      FileType               ftype    =  rule_itr->first;
      ArgumentRuleset::Count count    =  rule_itr->second;
      
      // - Should not see any switches.
      //
      if(filename[0] == '-')
      {
        throw std::exception();
      }
      
      // Based on the ruletype, figure out what to do:
      switch(count)
      {
        case ArgumentRuleset::ZERO_OR_ONE:
          if(! first_optional_found)
          { 
            found.add_argument(ftype, filename);
            first_optional_found = true; 
            ++rule_itr;
          }
          else
          {
            throw std::exception();
          } 
          break;
          
        case ArgumentRuleset::ONE:
          found.add_argument(ftype, filename); 
          ++rule_itr; 
          break;
          
        case ArgumentRuleset::ONE_OR_MORE: 
          found.add_argument(ftype, filename);
          first_optional_found = true;             
          break;
        
        case ArgumentRuleset::ZERO_OR_MORE:
          if(! first_optional_found)
          { 
            found.add_argument(ftype, filename);
            first_optional_found = true; 
            ++rule_itr;
          }
          else
          {
            throw std::exception();
          } 
          break;

        default:
          assert(false);      
      }
    }
  }

  // Now make sure that all the required arguments were specified:
  for(ArgumentRuleset::RuleVector::const_iterator rule_itr = ruleset.get_rules().begin(); rule_itr != ruleset.get_rules().end(); ++rule_itr)
  {
    // If it's REQUIRED and NOT FOUND, throw an exception:
    if(    (rule_itr->second == ArgumentRuleset::ONE || rule_itr->second == ArgumentRuleset::ONE_OR_MORE)
        && found.get_arguments(rule_itr->first).empty())
      throw std::exception();
  }
  
  // Return the arguments found instance:
  return found;
}

//==================================
//
//  process_commandline(...)
// 
//==================================
ArgumentsFound process_commandline(
        int                argc,   
        char            ** argv,   
  const std::string     &  program_name,
  const ArgumentRuleset &  ruleset,
        std::ostream    &  out,
        bool               exit_on_fail)
{
  // Try processing the commandline:
  try
  {
    return ArgumentParser::parse_commandline(ArgumentParser::convert_args(argc, argv), ruleset);
  }
  
  // On failure:
  catch(const std::exception & e)
  {
    // Display appropriate usage:
    ArgumentParser::display_usage(program_name, ruleset, out);
    
    // If flag was set, exit the program:
    if(exit_on_fail)
      exit(0);
      
    // Otherwise, just return a dummy instance:
    else
      return ArgumentsFound();
  }
}

} // End namespace APP
} // End namespace SAGE


