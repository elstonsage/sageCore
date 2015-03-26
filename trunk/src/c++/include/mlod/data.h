#ifndef MLOD_DATA_H
#define MLOD_DATA_H

#include <string>
#include <iostream>
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "data/SAGEdata.h"
#include "mlod/analysis_parameters.h"
#include "mlod/parser.h"

namespace SAGE
{
namespace MLOD
{

/// \brief MLOD-specific instantiation of SAGE_Data
///
/// A derived type of the APP::SAGE_Data, includes custom methods to read
/// in the data files and analyses.  Also contains new data structures for storing
/// the marker loci, trait loci, and complete marker information separately.
///
/// \section Note
///
/// MLOD has to read multiple Marker Locus files (one for markers,
/// one for traits) and has to store them separately.  It therefore
/// doesn't use the standard methods for reading these files, but
/// codes its own. which are used in the input function.
class Data : public APP::SAGE_Data
{
  public:
    typedef RPED::genome_description::region_type region;

    /// Basic constructor.
    ///
    Data(const string& program_name, bool debug = false);
    
    /// Basic Destructor
    ///
    ~Data();
    
    /// Input based upon command line arguments.  Replaces the APP::SAGE_Data
    /// function
    void  input(int argc, char** const argv);

    /// Read in the analyses to be performed
    ///
    virtual bool  read_analysis();
    
    /// Returns an inheritance model map containing \b just the marker loci
    ///
    /// \returns The markers in an inheritance_model_map
    const MLOCUS::inheritance_model_map& get_marker_loci() const;

    /// Returns an inheritance model map containing \b just the marker loci
    ///
    /// \returns The markers in an inheritance_model_map
    const MLOCUS::inheritance_model_map& get_trait_loci() const;

    /// Returns the vector of analyses that are to be run.
    ///
    const vector<AnalysisParameters>&  get_analyses() const;
    
  private:
  
    //lint -e{1704} Disabled default construction
    /// Default Constructor disabled for this class
    ///
    Data();
  
    /// Reads the marker loci and trait loci files and stores them separately
    /// and in the ref multipedigree.
    ///
    /// \param _argv The command line arguments which specify files to read.
    void read_locus_models(char* const _argv[]);
    
    /// Reads the file specified by fname and inserts the loci in it into the
    /// MLOCUS::inheritance_model_map specified by imap.
    ///
    /// Since MLOD has two different sets of loci (marker and trait), we want to 
    /// read them in separately to separate maps.  Once read in, the 
    /// copy_locus_models_to_multipedigree function transfers them all to a common place.
    ///
    /// \param imap  The inheritance map to read into
    /// \param fname The name of the file to read
    void read_locus_description_file_to_map(MLOCUS::inheritance_model_map& imap, 
                                            const string&                  fname);
    
    
    /// Once the marker loci are read in, they must be copied to the
    /// multipedigree.  Markers are copied first, followed by traits.
    void copy_locus_models_to_multipedigree();
  
    /// A map containing just the marker loci
    ///
    MLOCUS::inheritance_model_map my_marker_loci;
    
    /// A map containing just the trait loci
    ///
    MLOCUS::inheritance_model_map my_trait_loci;
    
    /// Vector containing analysis to be performed
    ///
    vector<AnalysisParameters>  my_analyses;
};

}
}

#include "mlod/data.ipp"


#endif

