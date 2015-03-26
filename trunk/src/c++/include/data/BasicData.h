#ifndef DATA_BASIC_DATA_H
#define DATA_BASIC_DATA_H

#include <vector>
#include "data/SAGEdata.h"
#include "data/ModelParser.h"
#include "containers/CompleteProcessMgr.h"

namespace SAGE {
namespace APP  {

//===============================================
//
//  BasicData
//
//===============================================
class BasicData : public APP::SAGE_Simple_Data
{
  friend class BasicApp;

public:

  /// @name Registering analysis types and their parsers
  //@{

    ///
    /// Registers a parameter block name as corresponding to a particular analysis type.
    ///
    /// \param ANALYSIS_TYPE The object that the parser will populate with content
    /// \param analysis_block_name The name of this analysis' subblock
    /// \param PARSER_TYPE An object with the following two functions defined:
    /// \code
    /// PARSER_TYPE(const PARSER_TYPE &);  // Copy constructor
    /// void operator() (ANALYSIS_TYPE & f, const LSFBase * block) const
    /// \endcode
    ///
    /// \returns \c true if successful, \c false otherwise
    template<class ANALYSIS_TYPE, class PARSER_TYPE> bool registerAnalysis(const std::string & analysis_block_name, const PARSER_TYPE & parser);

  //@}
                                      
  /// @name Analysis access
  //@{
  
    ///
    /// Returns a const begin iterator across all analyses of type ANALYSIS_TYPE.
    template<class ANALYSIS_TYPE> UntypedSet::ConstIterator<ANALYSIS_TYPE> getAnalysisBegin() const;

    ///
    /// Returns a const end iterator across all analyses of type ANALYSIS_TYPE.
    template<class ANALYSIS_TYPE> UntypedSet::ConstIterator<ANALYSIS_TYPE> getAnalysisEnd() const;

    ///
    /// Returns a non-const begin iterator across all analyses of type ANALYSIS_TYPE.
    template<class ANALYSIS_TYPE> UntypedSet::Iterator<ANALYSIS_TYPE> getAnalysisBegin();

    ///
    /// Returns a non-const end iterator across all analyses of type ANALYSIS_TYPE.
    template<class ANALYSIS_TYPE> UntypedSet::Iterator<ANALYSIS_TYPE> getAnalysisEnd();

  //@}

  /// \internal
  /// @name Required by inheritance (but not used)
  virtual void process_input(int argc, char** argv) {}

  /// \internal
  /// @name Required by inheritance (but not used)
  virtual bool read_analysis () { return true; }

private:

  class LSFBaseClassifier { public:
    typedef std::string return_type;
    std::string operator() (const LSFBase * input) const { return input->name(); } };

  ///
  /// Constructor.
  BasicData(const string & program_name, bool debug_on);

  void readParameterFile(const std::string & filename);

  void readPedigreeFile(const std::string & filename);

  void readLocusDescriptionFile(const std::string & filename);

  CompleteProcessMgr<UntypedSet, const LSFBase *, LSFBaseClassifier> my_analysis_parsers;
};

//=================
// CONSTRUCTOR
//=================

inline BasicData::BasicData(const string& program_name, bool debug_on) :
  SAGE_Simple_Data(program_name, debug_on)
{}

//=====================================
// getAnalysisBegin(), getAnalysisEnd()
//=====================================

template<class ANALYSIS_TYPE> inline UntypedSet::ConstIterator<ANALYSIS_TYPE> BasicData::getAnalysisBegin() const
{ return my_analysis_parsers.getContainer().begin<ANALYSIS_TYPE>(); }

template<class ANALYSIS_TYPE> inline UntypedSet::ConstIterator<ANALYSIS_TYPE> BasicData::getAnalysisEnd() const
{ return my_analysis_parsers.getContainer().end<ANALYSIS_TYPE>(); }

template<class ANALYSIS_TYPE> inline UntypedSet::Iterator<ANALYSIS_TYPE> BasicData::getAnalysisBegin()           
{ return my_analysis_parsers.getContainer().begin<ANALYSIS_TYPE>(); }

template<class ANALYSIS_TYPE> inline UntypedSet::Iterator<ANALYSIS_TYPE> BasicData::getAnalysisEnd()
{ return my_analysis_parsers.getContainer().end<ANALYSIS_TYPE>(); }

//=======================
// registerAnalysis() ...
//=======================

template<class ANALYSIS_TYPE, class PARSER_TYPE> 
inline bool 
BasicData::registerAnalysis(const std::string & analysis_block_name, const PARSER_TYPE & parser)
{
  return my_analysis_parsers.addProcessor<ANALYSIS_TYPE>(analysis_block_name, parser);
}

//===============================
//  readLocusDescriptionFile(...)
//===============================
inline void BasicData::readLocusDescriptionFile(const std::string & filename)
{
  read_locus_description_file(filename);
}

//=======================
// readParameterFile(...)
//=======================
inline void 
BasicData::readParameterFile(const std::string & filename)
{
  read_parameter_file(filename);
  
  for(LSFList::const_iterator i = my_params->List()->begin(); i != my_params->List()->end(); ++i)
    if(*i)
      my_analysis_parsers.processInput(*i);
}

//=======================
//  readPedigreeFile(...)  
//=======================
inline void BasicData::readPedigreeFile(const std::string & filename)
{
  read_family_data_file (filename, true, false, false, false, true);
  evaluate_functions    ();  
}

}} /// End namespace

#endif
