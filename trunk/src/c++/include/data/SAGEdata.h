#ifndef SAGEDATA_H
#define SAGEDATA_H

#include <list>
#include <vector>
#include <algorithm>
#include "mlocus/imodel.h"
#include "mlocus/mfile.h"
#include "rped/rped.h"
#include "rped/rpfile.h"
#include "rped/genome_description.h"
#include "LSF/LSFsymbol.h"
#include "app/output_streams.h"
#include "boost/smart_ptr.hpp"
#include "data/ArgumentRuleset.h"

namespace SAGE {
namespace APP  {

/** 
 *  Allows for iteration over markers in their original order after 
 *  RefMPedinfo's inheritance_model_map has changed due to marker filtration 
 *  (see SAGEdata::delete_marker()). 
 */
class marker_order
{
  public:
    marker_order();
    marker_order(const MLOCUS::inheritance_model_map* imodel);
    
    class iterator
    {
      public:
        friend class marker_order;
      
        iterator  operator ++();
        bool      operator ==(const iterator& other) const;
        bool      operator !=(const iterator& other) const;
        bool      operator  >(const iterator& other) const;
        const MLOCUS::inheritance_model*  operator *() const;
        size_t  imodel_index() const;

      private:
        iterator(const marker_order* mo, bool at_end = false);      
      
        // Data members
        const marker_order*  my_marker_order;
        vector<string>::const_iterator  my_position; 
    };
    
    friend class iterator;
    
    iterator  begin() const;
    iterator  end() const;
    iterator  find(size_t imodel_index) const;    
    
    void  remove_marker(const string& marker);    
    
  private:
    void  dump() const;
  
    // Data members
    const MLOCUS::inheritance_model_map*  my_imodel;
    vector<string>  my_marker_names;
};


/** \brief Basic data needed for most SAGE programs.
 *
 * The SAGE_Data class is the basic data needed for most SAGE programs. Not
 * all elements are needed for all programs, but the extra overhead in
 * having them included here is made up for in the ease of maintainability.
 */
class SAGE_Data
{
public:

    //lint --e{1712} no default constructor ok

  /// @name Constructor/destructor
  //@{

    ///
    /// Constructor
    /// \param program_name The name of the program
    /// \param debug Boolean option indicating whether or not to include debug information
    SAGE_Data(const string& program_name, bool debug);

    ///
    /// Destructor
    virtual ~SAGE_Data();

  //@}

  /// @name Required virtual interface
  //@{

    /**
      * Reads in all the analyses for this program.
      *
      * The implementation of this function should take the form:
      * \code
      * parser_object my_parser;
      * 
      * foreach LSFObject in my_params:
      *   if(LSFObject is an analysis block):
      *     my_parser.parse(LSFObject)
      * \endcode
      * 
      * Here is a real-world example from assoc:
      * 
      * \code
      * bool
      * assoc_data::read_analysis()
      * {
      *   // 1. Create the local model and parser
      *             
      *         assoc_parser my_parser(pedigrees(), my_output.messages(), my_output.info(), errors());
      *         assoc_model  my_model;
      *               
      *   // 2. Clear out the vector of analyses
      *          
      *         my_analyses.clear();
      *         
      *   // 3. Parse the parameter file
      *   
      *         for(LSFList::const_iterator i = my_params->List()->begin(); i != my_params->List()->end(); ++i)
      *         {
      *           if(!*i) continue;
      *  
      *           if(toUpper((*i)->name()) == "ASSOC_ANALYSIS" || toUpper((*i)->name()) == "ASSOC")
      *           {
      *             my_parser.parse_test_parameter_section(*i);
      *         
      *             my_model = my_parser.get_model();
      *   
      *             if(my_model.calc_valid(&(my_output.messages())))
      *               my_analyses.push_back(my_model);
      *           }
      *         }
      * 
      *         return true;
      * }        
      * \endcode
      */
    virtual bool read_analysis() = 0;

  //@}

  /// @name Input/parsing routines
  //@{
    /// Reads in the parameter file and parses for a few key blocks (MARKER blocks
    /// and marker delimiter parameters)
    ///
    /// \param fname The name of the file to read in
    bool read_parameter_file  (const string& fname);

    ///
    /// Reads in the pedigree data file.
    /// \param fname The name of the file to read in
    /// \param dump_trait ?
    /// \param dump_marker ?
    /// \param skip_traits ?
    /// \param skip_markers ?
    /// \param dynamic_markers ?
    bool read_family_data_file(const string& fname, bool dump_trait      = true,
                                                    bool dump_marker     = true,
                                                    bool skip_traits     = false,
                                                    bool skip_markers    = false,
                                                    bool dynamic_markers = false);

    ///
    /// Reads in the locus description file.
    ///
    /// This function uses the RPED::PhenoReaderInfo object found in the
    /// pedigrees for its parsing.
    ///
    /// \param fname The name of the file to read in
    bool read_locus_description_file (const string& fname);

    ///
    /// Reads in the genome description file.
    /// \param fname The name of the file to read in
    bool read_genome_description_file(const string& fname);

    ///
    /// Reads in the genome description file, with a few more options.
    /// \param fname The name of the file to read in
    /// \param multipoint Boolean value indicating whether or not to do a multipoint thingy
    /// \param distance Genetic distance thingy
    bool read_genome_description_file(const string& fname, bool multipoint, double distance);

  //@}

  /// @name Evaluation routines
  //@{

    ///
    /// Evaluate function blocks.
    void evaluate_functions();

    ///
    /// Evaluate function blocks.
    /// \param mped The data source to be used
    void evaluate_functions(RPED::RefMultiPedigree& mped);

  //@}

  /// @name Information printing routines
  //@{

    ///
    /// Print genome information file.
    /// \param markers ???
    void print_genome_info_file(const MLOCUS::inheritance_model_map & markers) const;

  //@}

  /// @name Error checking functions
  //@{

    ///
    /// Check for structural errors.
    /// \param mped The data source whose structure will be checked for errors
    void check_family_data(const RPED::RefMultiPedigree& mped);

  //@}

  /// @name Filtering functions
  //@{

    ///
    /// Filter out trait entries with no valid data.
    /// \param t The name of the trait to be filtered
    void filter_trait(size_t t);

    ///
    /// Filter out marker entries with no valid data.
    /// \param The number of the marker to be filtered.
    void filter_marker(size_t m);

    ///
    /// Filter out all traits with no valid data.
    void filter_traits();

    enum  FilterType
    {
      ft_NO_DATA,
      ft_SINGLE_ALLELE,
      ft_X_LINKED,
      ft_NONE
    };

    ///
    /// Filter out all markers with no valid data.
    void filter_markers();
    
    ///
    /// Filter out all markers with only one allele.
    void filter_single_allelic_markers();    
    
    ///
    /// Filter out all x-linked markers.
    void filter_x_linked_markers();

  //@}

  /// @name Member data access
  //@{
  
    const ArgumentsFound&  parsed_arguments();

    ///
    /// Returns a non-const pointer to the LSF parameters associated with this object.
    LSFBase * parameters();

    ///
    /// Returns a const pointer to the LSF parameters associated with this object.
    const LSFBase * parameters() const;

    ///
    /// Returns a non-const reference to the RPED::RefMultiPedigree associated with this object.
    RPED::RefMultiPedigree & pedigrees();

    ///
    /// Returns a const reference to the RPED::RefMultiPedigree associated with this object.
    const RPED::RefMultiPedigree & pedigrees() const;

    ///
    /// Returns a non-const reference to the structure containing information about the first ten individuals listed 
    /// in the RPED::RefMultiPedigree.
    ///
    /// It should be noted that, although the variable is \b called first_ten_ind, it may in fact represent a variable
    /// number of individuals. This is because the user may request a different number of individuals to be reported
    /// in the information file.
    ///
    /// The structure takes the form of a vector of pairs, where the first element in the pair is the pedigree name,
    /// and the second element in the pair is the individual name.
    vector<pair<string, string> >& first_ten_ind();

    ///
    /// Returns a const reference to the structure containing information about the first ten individuals list
    /// in the RPED::RefMultiPedigree. See the non-const version for more information.
    const vector<pair<string, string> >& first_ten_ind() const;

    ///
    /// Returns a const reference to the interhitance_model_map that describes the markers that have been read in.
    const MLOCUS::inheritance_model_map & markers() const;

    ///
    /// Returns a non-const pointer to the interhitance_model_map that describes the markers that have been read in.
    RPED::genome_description * genome();

    ///
    /// Returns a const pointer to the interhitance_model_map that describes the markers that have been read in.
    const RPED::genome_description * genome() const;

    ///
    /// Returns a non-const pointer to the LSFBase that describes the genomic regions that have been read in.
    LSFBase * regions();

    ///
    /// Returns a const pointer to the LSFBase that describes the genomic regions that have been read in.
    const LSFBase * regions() const;

  //@}
 
  /// @name Output stream accessing functions
  //@{

    ///
    /// Returns a non-const reference to this object's Output_Streams member.
    Output_Streams & get_ostreams() const;

    ///
    /// Returns this object's outputstream for error messages.
    cerrorstream & errors() const;

    ///
    /// Returns this object's outputstream for informative output.
    ostream & info() const;

    ///
    /// Returns this object's outputstream for screen output.
    ostream & screen() const;

    ///
    /// Returns this object's outputstream for messages.
    ostream &  messages() const;

  //@}
    const marker_order&  original_marker_order() const;
    
    const ArgumentsFound&  parsed_arguments() const;
    
  protected:
  
    void  parse_cmdline(int argc, char** argv);  
  
    /// Finds and parses the marker information within a parameter file.
    /// This includes both the marker block and the ALLELE_DELIMITER/MARKER_DELIMITER
    /// variables.
    void parse_pheno_info();

    void dump_new_traits(const RPED::RefMultiPedigree& mped, size_t org_trait_count) const;

    void delete_trait(size_t t);
    void  delete_marker(size_t m, FilterType filter = ft_NONE);

    bool trait_data_exist(size_t t) const;
    bool marker_data_exist(size_t m) const;
    bool only_one_allele(size_t m) const;

    string                              my_program_name;
    ArgumentRuleset                     my_cmdline_rules;
    ArgumentsFound                      my_parsed_arguments;
    LSF_ptr<LSFBase>                    my_params;
    RPED::RefMultiPedigree              my_pedigrees;
    mutable Output_Streams              my_output;
    vector<pair<string, string> >       my_first_ten_ind;
    char                                my_allele_delimiter;

    LSF_ptr<RPED::genome_description>   my_genome;
    LSF_ptr<LSFBase>                    my_regions;
    marker_order                        my_marker_order;

  private:
  
};

typedef SAGE_Data SAGE_Simple_Data;

} // End namespace APP
} // End namespace SAGE

#include "data/SAGEdata.ipp"

#endif
