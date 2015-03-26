#ifndef DATA_BASIC_APP_H
#define DATA_BASIC_APP_H

#include <vector>
#include "app/SAGEapp.h"
#include "data/SAGEdata.h"
#include "data/BasicData.h"

namespace SAGE {
namespace APP  {

//===============================================
//
//  BasicApp
//
//===============================================
class BasicApp : public APP::SAGEapp
{
  public:

    enum StandardFileType
    {
      PARAMETER    = 0,
      PEDIGREE     = 1,
      MARKER_LOCUS = 2
    };
    
    enum NonstandardFileType
    {
      SPECIAL1 = MARKER_LOCUS + 1,
      SPECIAL2,
      SPECIAL3,
      SPECIAL4,
      SPECIAL5,
      SPECIAL6,
      SPECIAL7,
      SPECIAL8,
      SPECIAL9
    };
    
    enum RequirementLevel
    {
      REQUIRED,
      OPTIONAL
    };

    struct FileReaderInputType
    {
      FileReaderInputType(int _filetype, const std::string & _filename) : filetype(_filetype), filename(_filename) {}
      int         filetype;
      std::string filename;
    };

  /// @name Construction and setup
  //@{

    ///
    /// Constructor
    /// \param
    BasicApp(const app_index_type program_index, int argc, char ** argv);

    ///
    /// For standard input files (see StandardFileType), there are already predefined file readers.
    /// If such a file is required by your program, you should use this function to indicate as such.
    ///
    /// \param FILETYPE The type of input file
    /// \param level The requirement level (required or optional)
    ///
    /// \returns \c true if successful, \c false otherwise
    bool requireStandardInputFile(StandardFileType type, RequirementLevel level);

    ///
    /// If a nonstandard file is required as input, it's up to you (the programmer) to indicate as such,
    /// as well as write a file reading object to handle it. The FILE_READER_TYPE class must have a
    /// the following public member functions defined:
    ///
    /// \code
    /// FILE_READER_TYPE(const FILE_READER_TYPE & other) // Copy constructor
    /// void operator() (const BasicApp::FileReaderInputType & input) const // operator()
    /// \endcode
    ///
    /// \returns \c true if successful, \c false otherwise (see note above)
    template<class FILE_READER_TYPE>
    bool requireNonstandardInputFile(
            NonstandardFileType   type,
            RequirementLevel      level,
      const std::string         & name,
      const FILE_READER_TYPE    & file_reader);

    ///
    /// Having instructed the BasicApp as to which input files are required & optional, this function
    /// will then process a given commandline on the basis of those options.
    void processCommandLine();

  //@}
  
  /// @name Associated data access
  //@{

    ///
    /// Returns the data object associated with this application (through which
    /// the various output streams and runtime input are available).
    BasicData & getData();

    ///
    /// Returns the data object associated with this application (through which
    /// the various output streams and runtime input are available).
    const BasicData & getData() const;

  //@}
  
  /// \internal
  /// @name Required virtual functions because of inheritance
  //@{

    /// \internal
    /// Required because of inheritance; not used
    virtual int main();

  //@}

  private:

    //===================
    // STRUCTS & TYPEDEFS
    //===================

    struct InputFileInfo
    {
      InputFileInfo(int _type, RequirementLevel _level, const std::string & _name)
        : type(_type), level(_level), name(_name) {}

      int              type;
      RequirementLevel level;
      std::string      name;
    };

    typedef std::vector<InputFileInfo> InputFileInfoVector;


    //========
    // CLASSES
    //========

    class ParameterFileReader
    { 
    public:
      explicit ParameterFileReader(BasicData & data) : my_data(data) {}
      ParameterFileReader(const ParameterFileReader & other) : my_data(other.my_data) {}

      void operator() (const BasicApp::FileReaderInputType & input) const
      { my_data.readParameterFile(input.filename); }

    private:
      BasicData & my_data;
    };   
    
    class PedigreeFileReader
    { 
    public:
      explicit PedigreeFileReader(BasicData & data) : my_data(data) {}
      PedigreeFileReader(const PedigreeFileReader & other) : my_data(other.my_data) {}

      void operator() (const BasicApp::FileReaderInputType & input) const
      { my_data.readPedigreeFile(input.filename); }

    private:
      BasicData & my_data;
    };   

    class LocusDescriptionFileReader
    { 
    public:
      explicit LocusDescriptionFileReader(BasicData & data) : my_data(data) {}
      LocusDescriptionFileReader(const LocusDescriptionFileReader & other) : my_data(other.my_data) {}

      void operator() (const BasicApp::FileReaderInputType & input) const
      { my_data.readLocusDescriptionFile(input.filename); }

    private:
      BasicData & my_data;
    };   

    class FileTypeClassifier
    {
    public:
      typedef int return_type;
      int operator() (const FileReaderInputType & i) const { return i.filetype; }
    };

    //==========
    // FUNCTIONS
    //==========

    void printUsage();

    //=============
    // DATA MEMBERS
    //=============

    BasicData                                             my_data;
    InputFileInfoVector                                   my_input_files;
    ProcessorMgr<FileReaderInputType, FileTypeClassifier> my_file_readers;
};

//=================
// INLINE FUNCTIONS
//=================

inline const BasicData & BasicApp::getData() const { return my_data; }
inline       BasicData & BasicApp::getData()       { return my_data; }
    
inline int  BasicApp::main         ()                 { return 0; }

}} /// End namespace

#endif
