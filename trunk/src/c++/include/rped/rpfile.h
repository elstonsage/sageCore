#ifndef REF_PED_FILE_H
#define REF_PED_FILE_H

#include <set>
#include <vector>
#include <list>
#include <iostream>
#include "LSF/LSF.h"
#include "rped/rped.h"
#include "error/errorstream.h"
#include "error/bufferederrorstream.h"

namespace SAGE {
namespace RPED {

/** \brief Describes meta-information for reading in a pedigree data file
  *
  * \par Introduction
  *
  * When reading in a pedigree data file, there are all sorts of options that need to 
  * be set (sex code, individual code, delimited vs. fortran, etc.). The
  * RefPedigreeFile is a base class that stores pedigree data importation options
  * common to both fortran-style and character-delimited data files.
  *
  * Please note that this is a \b pure virtual class, and as such cannot be directly
  * instantiated.
  *
  * \par Data importation 
  *
  * The first step in data importation is to set the various importation configurartion options.
  * The documentation for these various options is available in the description of each function.
  *
  * The next step is to invoke the input() method (see below for more details).
  *
  * \par Output generation 
  *
  * \anchor RefPedigreeFile_Validity_functions
  * \par Validity functions 
  *
  * \par Statistics on imported data 
  *
  * \par Imported fields 
  *
  * \par Quick data extraction functions 
  *
  */
class RefPedigreeFile
{
public:

  /// Uniquely identifies the type of a field.
  enum field_t 
  { 
    skip,          /*!< Skipped over                       */
    study_id,      /*!< Indicates the study id number      */
    pedigree_id,   /*!< Indicates the pedigree id number   */
    individual_id, /*!< Indicates the individual id number */
    parent_id,     /*!< Indicates a parent id number       */
    sex_code,      /*!< Indicates an individual's sex      */
    trait,         /*!< Indicates a trait-type field       */
    marker,        /*!< Indicates a marker                 */
    allele,        /*!< Indicates an allele                */
    string_field,  /*!< Indicates a string-type field      */
    allele_cov,    /*!< Indicates a covariate from marker data field      */
    marker_cov     /*!< Indicates a covariate from marker data field      */
  };
  
  /** \brief Describes a field read in during data importation
    *
    * Describes a field read in during data importation.
    */
  struct field
  {
    /// Constructor
    field(field_t t = skip, const std::string & f = "", const std::string & n  ="",  size_t i = (size_t)-1);

    /// Type of the field (see RefPedigreeFile::field_t
    field_t type;
    
    /// Name of the field as indicated by the pedigree block
    std::string field_name;

    /// The name of the field as indicated in the pedigree data file
    std::string name;
    
    /// Index number of the field
    size_t index;

    /// Marks this field as invalid.
    void invalidate();
  };

  /// List of fields
  typedef std::list<field> field_list_type;
  
  /// Map of fields (mapping field name to field instance)
  typedef std::map<std::string, field> field_map_type;

  /// @name Constructor / destructor
  //@{
  
    ///
    /// Constructor
    /// \param err Errorstream to which messages will be sent
    RefPedigreeFile(cerrorstream &err = sage_cerr);

    ///
    /// Destructor
    virtual ~RefPedigreeFile();
    
  //@}

  /// @name Option-setting for data importation (mutators)
  //@{
  
    ///
    /// Sets the errorstream to which messages will be directed.
    /// \param err The errorstream
    void set_error_sink(cerrorstream &err);

    ///
    /// Sets the number of individuals whose pedigree id and name will be reported by the get_ind_list() function.
    ///
    /// Please note that the default value for this feature (at construction) is 10.
    /// \param v Number of individuals
    void set_verbose_output(size_t v = (size_t)-1);

    ///
    /// Sets whether or not to reject individuals whose lineage is not complete (ie: only one parent specified).
    ///
    /// Please note that the default value for this feature (at construction) is \c true.
    /// \param r Boolean value indicating whether or not to reject the individual
    void set_reject_partial_lineage(bool r = true);
    
    ///
    /// Sets whether or not a record is required for every indicated individual. If this is set to false, records
    /// will be \b created for individuals lacking records. For instance, if individual #10 lists parents #1 and #2,
    /// but there are no records for #1 and #2, then the RefPedigreeFile will create empty entries for those parents if
    /// set_require_record() is set to false.
    ///
    /// Please note that the default value for this feature (at construction) is \c false.
    /// \param a Boolean value indicating whether or not to require a record
    void set_require_record(bool a = true);
    
    ///
    /// Sets whether or not to skip reading in traits.
    void set_skip_traits(bool s = true);
    
    ///
    /// Sets whether or not to skip marker reading.
    void set_skip_markers(bool s = true);
    
    ///
    /// Sets whether or not to do 'dynamic markers' ???
    void set_dynamic_markers(bool d = true);
    
    ///
    /// Sets whether or not to forcibly skip traits ???
    void set_force_skip_traits(bool s=true);
    
    ///
    ///
    void set_force_skip_markers(bool s = true);
    
    ///
    ///
    void set_force_dynamic_markers(bool d = true);
    
    ///
    ///
    void set_sex_linked_exist(bool s = true);

    ///
    ///
    void set_build_incremental(bool b = true);
    
    ///
    ///
    void set_sex_code_trait(bool a = false);

    ///
    ///
    void set_sex_field_name(const std::string & s);

    ///
    ///
    void set_pedigree_id_trait(bool a = false);

    ///
    ///
    void set_pedigree_id_name(const std::string & s);
    
    ///
    /// Records the "FORMAT" option from the pedigree block (see section 2.3.2.1 of the SAGE 4.6 manual).
    /// \param f The value of the "FORMAT" option
    void set_format(const std::string & f);
    
    ///
    /// Sets the treat_ped_id option.
    void set_treat_as_sibs(bool t);
    
    ///
    /// Sets the no_sex_field option.
    void set_no_sex_field(bool s);

    ///
    /// Sets the no_sex_ok option.
    void set_no_sex_ok_option(bool s);
    
  //@}

  /// @name Option-getting for data importation (accessors)
  //@{

    ///
    /// Returns the errorstream to which error messages from this object are directed.
    cerrorstream error_sink() const;

    ///
    /// Returns the number of individuals whose pedigree id / name will be stored (see set_verbose_output() )
    size_t verbose_output() const;
    
    ///
    /// Returns whether or not individuals with partial lineage will be rejected (see set_reject_partial_lineage() )
    bool reject_partial_lineage() const;

    ///
    /// Returns whether or not records are required for all individuals (see set_require_record() )
    bool require_record() const;

    ///
    ///
    bool skip_traits() const;

    ///
    ///
    bool skip_markers() const;

    ///
    ///
    bool dynamic_markers() const;

    ///
    ///
    bool force_skip_traits() const;

    ///
    ///
    bool force_skip_markers() const;

    ///
    ///
    bool force_dynamic_markers() const;

    ///
    ///
    bool sex_linked_exist() const;

    ///
    ///
    bool build_incremental() const;

    ///
    ///
    bool sex_code_trait() const;

    ///
    ///
    const std::string & sex_field_name() const;

    ///
    ///
    bool pedigree_id_trait() const;

    ///
    ///
    const std::string & pedigree_id_name() const;

    ///
    /// Returns the value of the "FORMAT" option from the pedigree block (see set_format() )
    const std::string & format() const;
    
    ///
    /// Returns the treat_ped_id option.
    bool get_treat_as_sibs() const;
    
    ///
    /// Returns the no sex field option.
    bool get_no_sex_field() const;

    ///
    /// Returns the no_sex_ok option.
    bool get_no_sex_ok_option() const;
  //@}

  /// @name Data importation
  //@{
  
    ///
    /// Assuming that the configuration settings have been set, this function reads in all data from 
    /// the indicated inputfile, parses it, and places it in the indicated RefMultiPedigree.
    ///
    /// \param p The RefMultiPedigree into which data will be placed
    /// \param filename The name of the input file
    /// \param messages The outputstream to which messages will be sent
    /// \retval true Data was successfully imported
    /// \retval false Data was \b not successfully imported
    virtual bool input(RefMultiPedigree &p, const std::string &filename, ostream &messages = cout);

    ///
    /// \internal
    /// Imports pedigree \b structure.
    virtual bool input_pedigree(RefMultiPedigree &p, const std::string &filename, ostream &messages = cout, bool quiet = false) = 0;

    ///
    /// \internal
    /// Imports individual trait/marker data.
    virtual bool input_data(RefMultiPedigree &p, const std::string &filename, ostream &messages = cout, bool quiet = false) = 0;


  //@}
  
  /// @name Validity
  //@{

    ///
    /// Sets the internal validity state to \c true (valid)
    ///
    /// See \ref RefPedigreeFile_Validity_functions "detailed description" for details.
    void validate();
    
    ///
    /// Sets the internal validity state of \c false (invalid)
    ///
    /// See \ref RefPedigreeFile_Validity_functions "detailed description" for details.
    void invalidate();
    
    ///
    /// Returns the internal validity state.
    ///
    /// See \ref RefPedigreeFile_Validity_functions "detailed description" for details.
    bool valid() const;

  //@}

  /// @name Output generation
  //@{

    ///
    /// Prints information on the first N individuals in the RefMultiPedigree (where n is given by verbosity() ).
    /// \param p The RefMultiPedigree from which the data values will be read
    /// \param filename The name of the original file from which the data was read
    /// \param messages The outputstream to which the output will be sent
    /// \param dump_trait Boolean value indicating whether or not to include trait values in the output
    /// \param dump_marker Boolean value indicating whether or not to include marker values in the output
    virtual void print_mped(const RefMultiPedigree &p, 
                            const std::string &filename, 
                            ostream &messages = cout,
                            bool dump_trait = true, 
                            bool dump_marker = false);

    ///
    /// Prints the entire RefMultiPedigree's contents into the indicated output file, using the same
    /// formatting options that were used to import the data.
    /// \param p The RefMultiPedigree from which the data values will be read
    /// \param filename The name of the file in which the output will be created
    /// \param messages The outputstream to which informative messages will be sent
    /// \param quiet Boolean value indicating whether or not to include verbose output
    virtual bool output(
            RefMultiPedigree & p, 
      const std::string      & filename, 
            ostream          & messages = cout, 
            bool               quiet    = false) = 0;
          
                              
  //@}
  
  /// @name Statistics on imported data
  //@{

    ///
    /// Returns the number of fields imported.
    size_t field_count() const;
    
    ///
    /// Returns the number of trait-type fields imported.
    size_t trait_count() const;
    
    ///
    /// Returns the number of invalid trait-type fields encountered during import.
    size_t invalid_trait_count() const;
    
    ///
    /// Returns the number of string-type fields imported.
    size_t string_count() const;
    
    ///
    /// Returns the number of invalid string-type fields encountered durint import.
    size_t invalid_string_count() const;
    
    ///
    /// Returns the number of markers imported.
    size_t marker_count() const;
    
    ///
    /// Returns the number of invalid marker fields encountered during import.
    size_t invalid_marker_count() const;
    
    ///
    /// Returns the number of marker covariates imported.
    size_t marker_cov_count() const;
    
    ///
    /// Returns the number of invalid marker covariate fields encountered during import.
    size_t invalid_marker_cov_count() const;
    
    ///
    /// Returns the number of sex-type fields encountered during import.
    size_t sex_count() const;
    
    ///
    /// Returns the number of fields skipped during import.
    size_t skip_count() const;
    
    ///
    /// Returns the number of study-id-type fields encountered during import.
    size_t study_id_count() const;
    
    ///
    /// Returns the number of pedigree-id-type fields encountered during import.
    size_t pedigree_id_count() const;
    
    ///
    /// Returns the number of individual-id-type fields encountered during import.
    size_t individual_id_count()  const;
    
    ///
    /// Returns the number of parent-id-type fields encountered during import.
    size_t parent_id_count() const;

  //@}
  
  /// @name Imported fields
  //@{

    ///
    /// Returns a const-reference to the list of imported fields.
    const field_list_type & field_list() const;

    ///
    /// Returns a non-const reference to the list of imported fields.
    field_list_type & field_list();

  //@}
  
  /// @name Summary data extraction functions
  //@{
  
    ///
    /// Returns a vector of pedigree id / individual id pairs.
    const vector<pair<std::string, std::string> > & get_ind_list() const;

  //@}
  

protected:

  void reset_counts();

  /// @name Printing family structure
  //@{

    void print_family_structure_header(ostream &messages, const std::string &filename) const;
    void print_family_structure_footer(ostream &messages) const;
    void print_family_structure       (ostream &messages,
                                       const std::string& pn,
                                       const std::string& in,
                                       const std::string& s,
                                       const std::string& p1,
                                       const std::string& p2) const;
  //@}
  
  /// @name Printing trait info
  //@{

    void print_trait_header (ostream &messages, const RefMPedInfo &mped_info,  const std::string &filename) const;
    void print_trait_footer (ostream &messages) const;
    void print_trait        (ostream &messages,
                             const RefMPedInfo &mped_info,
                             const std::string& pn,
                             const std::string& in,
                             const vector<pair<size_t,std::string> >& trait_values,
                             const vector<pair<size_t,std::string> >& string_values) const;

  //@}
  
  /// @name Printing marker info
  //@{
  
    void print_marker_header(ostream &messages, const RefMPedInfo &mped_info, const std::string &filename) const;
    void print_marker_footer(ostream &messages) const;
    void print_marker       (ostream &messages,
                             const RefMPedInfo &mped_info,
                             const std::string& pn,
                             const std::string& in,
                             const vector<pair<size_t, std::string> >& marker_values) const;

  //@}

  void parse_dynamic_markers_missing(RefMarkerInfo& d_marker, const std::string& missing);

  bool validate_fields(bool data_only = false, bool quiet = false);

  bool build_pedigree(RefMultiPedigree &p);
  bool build_data(RefMultiPedigree &p);

  void add_member(RefMultiPedigree &p, const std::string &pn,
                  const std::string &id, const std::string &sex,
                  const std::string &parent1, const std::string &parent2,
                  size_t line, size_t count=(size_t)-1);

  // Process allele information modification option for a marker.
  //
  void parse_allele_frequency_adjustment(RefMarkerInfo& marker, const LSFBase *param);

  void normalize_allele_freq(RefMarkerInfo& m);
  void complement_allele_freq(RefMarkerInfo& m);
  void adjust_min_allele_freq(RefMarkerInfo& m, double f);
  void adjust_max_allele_freq(RefMarkerInfo& m, double f);

  void update_marker_delimiter_info(RefMarkerInfo& marker, char sep);
  void update_marker_missing_info  (RefMarkerInfo& marker, const std::string& missing);

  // Update sex-linked marker information for members according to sex.
  //
  void update_sex_linked_marker_info(RefMultiPedigree &mp);

  /// Test for missing sexes of individuals in pedigree structures.  If there
  /// are individuals with unknown sex after teh build is complete (inferred 
  /// sexes are ok) that are linked into a pedigree structure, this function
  /// produces a warning message that indicates this to the user.
  ///
  /// NOTE:  This is a temporary item, created specifically to deal with ticket
  /// #1601 dealing with the fact that many of our algorithms require sex,
  /// even though the analysis results are not affected by it.  This warning 
  /// should be removed post 5.3 when sex checking should be done by the 
  /// specific algorithms, which can take actions appropriate to thier needs, 
  /// rather than globally where we don't have any ability to determine what 
  /// the needs are.
  bool do_no_sex_structural_test(const RefMultiPedigree& mp);

  bool report_pedigree_build_errors   (const Pedigree &p) const;
  bool report_pedigree_build_warnings (const Pedigree &p) const;

  // Data members
  //
  field_list_type my_fields;

  size_t     my_skip_count;
  size_t     my_trait_count;
  size_t     my_invalid_trait_count;
  size_t     my_string_count;
  size_t     my_invalid_string_count;
  size_t     my_marker_count;
  size_t     my_invalid_marker_count;
  size_t     my_marker_cov_count;
  size_t     my_invalid_marker_cov_count;
  size_t     my_sex_count;
  size_t     my_study_id_count;
  size_t     my_pedigree_id_count;
  size_t     my_individual_id_count;
  size_t     my_parent_id_count;
 
  mutable cerrorstream errors;

  std::vector<std::pair<std::string, std::string> > my_ind_list;

private:

  bool         my_treat_as_sibs;
  size_t       my_verbose_output;
  bool         my_reject_partial_lineage;
  bool         my_require_record;
  bool         my_skip_traits;
  bool         my_skip_markers;
  bool         my_dynamic_markers;
  bool         my_force_skip_traits;
  bool         my_force_skip_markers;
  bool         my_force_dynamic_markers;
  bool         my_sex_linked_exist;
  bool         my_build_incremental;
  bool         my_sex_code_trait;
  bool         my_pedigree_id_trait;
  bool         my_no_sex_field;
  bool         my_no_sex_ok_option;
  bool         my_valid;
  std::string  my_sex_field_name;
  std::string  my_pedigree_id_name;
  std::string  my_format;
  
  // This set is only used for checking for duplicate records.
  std::multiset<std::pair<std::string, std::string> > my_inds;
};

/** \brief Stores information for reading in a fortran-formatted file
  *
  * Stores information for reading in a fortran-formatted file.
  */
class RefFortranPedigreeFile : virtual public RefPedigreeFile
{
public:

  /// @name Constructor/destructor
  //@{
  
    explicit RefFortranPedigreeFile(cerrorstream &err = sage_cerr);

    ~RefFortranPedigreeFile();

  //@}

  virtual bool input_pedigree(RefMultiPedigree &p, const std::string &filename, ostream &messages = cout, bool quiet = false);
  virtual bool input_data(RefMultiPedigree &p, const std::string &filename, ostream &messages = cout, bool quiet = false);
  virtual bool output(RefMultiPedigree &p, const std::string &filename, ostream &messages = cout, bool quiet = false);

  void add_skip_field();
  void add_study_id();
  void add_pedigree_id();
  void add_individual_id();
  void add_parent_id();
  void add_sex_field();
  void add_trait_field(const std::string &trait_name);
  void add_string_field(const std::string &name);
  void add_marker_field(const std::string &marker_name);
  void add_allele_field(const std::string &marker_name);
  
protected:
  bool nondefault_field_order;

private:
  bool build_fields(const RefMPedInfo &mped_info, bool quiet = false);
};

/** \brief Stores information for reading in a character-delimited file
  *
  * Stores information for reading in a character-delimited file.
  */
class RefDelimitedPedigreeFile : virtual public RefPedigreeFile
{
public:
  
  RefDelimitedPedigreeFile(cerrorstream &err = sage_cerr);

  ~RefDelimitedPedigreeFile();

  virtual bool input(RefMultiPedigree &p, const std::string &filename, ostream &messages = cout);

  virtual bool  input_pedigree(RefMultiPedigree &p, const std::string &filename, ostream &messages = cout, bool quiet = false);

  virtual bool  input_data(RefMultiPedigree &p, const std::string &filename, ostream &messages = cout, bool quiet = false);

  virtual bool output(RefMultiPedigree &p, const std::string &filename, ostream &messages = cout, bool quiet = false);

  bool format_in_file() const;

  const std::string &delimiters() const;

  const std::string &whitespace() const;

  bool skip_consecutive_delimiters() const;

  bool skip_leading_delimiters()     const;

  bool skip_trailing_delimiters()    const;

  void set_format_in_file(bool f);

  void set_delimiters(const std::string &d);

  void set_whitespace(const std::string &w);

  void set_skip_consecutive_delimiters(bool skip=true);

  void set_skip_leading_delimiters(bool skip=true);

  void set_skip_trailing_delimiters(bool skip=true);

  void add_study_id_field(const std::string &field_name);

  void add_pedigree_id_field(const std::string &field_name);

  void add_individual_id_field(const std::string &field_name);

  void add_parent_id_field(const std::string &field_name);

  void add_sex_field(const std::string &field_name);

  void add_trait_field(const std::string &field_name, const std::string &trait_name);

  void add_string_field(const std::string &field_name, const std::string &string_name);

  void add_allele_field(const std::string &field_name, const std::string &marker_name);

  void add_marker_field(const std::string &field_name, const std::string &marker_name);

  void add_allele_cov_field(const std::string &field_name, const std::string &mcov_name);

  void add_marker_cov_field(const std::string &field_name, const std::string &mcov_name);

  const field_map_type & field_map() const;

  field_map_type & field_map();

protected:

  bool validate_file  (const std::string &filename);
  bool validate_format(const std::string &filename);
  
  void setup_tokenizer(string_tokenizer& tokenizer);


private:
  bool build_fields(string_tokenizer &header, const RefMPedInfo &mped_info,
                    bool quiet = false);
                    
  void build_field(const RefMPedInfo& mped_info, field_map_type::const_iterator f,
                   bool quiet);
                   
  void build_trait_field  (const RefMPedInfo& mped_info, field& current_field);
  void build_marker_field (const RefMPedInfo& mped_info, field& current_field);
  void build_string_field (const RefMPedInfo& mped_info, field& current_field);

  bool check_header_vs_field_map(string_tokenizer &header);

  virtual string get_marker_covariate_value(string m, const string& v1, const string& v2) { return ""; }

  bool   my_format_in_file;
  std::string my_delimiters;
  std::string my_whitespace;
  
  bool my_skip_leading_delimiters;
  bool my_skip_trailing_delimiters;
  bool my_skip_consecutive_delimiters;

  field_map_type  my_field_map;
};


/** \brief Reads in LSF parameters general to all kinds of pedigree files (delimited & fortran)
  *
  * Reads in LSF parameters general to all kinds of pedigree files (delimited & fortran).
  */
class RefLSFPedigreeFile : virtual public RefPedigreeFile
{
public:
  ///
  /// Constructor.
  RefLSFPedigreeFile(cerrorstream &errors = sage_cerr);

  ///
  /// Processes an entire "PEDIGREE" block.
  /// \param mped_info The RefMPedInfo instance associated with the RefMultiPedigree in question
  /// \param params The LSFBase pointer corresponding to the pedigree block
  virtual bool process_parameters(RefMPedInfo &mped_info, const LSFBase *params);

  ///
  /// Processes a single parameter in the "PEDIGREE" block.
  /// \param mped_info The RefMPedInfo instance associated with the RefMultiPedigree in question
  /// \param param The LSFBase pointer corresponding to the parameter
  virtual bool process_parameter(RefMPedInfo &mped_info, const LSFBase *param);

  bool process_sex_code                 (RefMPedInfo &mped_info, const LSFBase* param);

private:

  bool process_format                   (RefMPedInfo &mped_info, const LSFBase* param);
  bool process_verbose                  (RefMPedInfo &mped_info, const LSFBase* param);
  bool process_require_record           (RefMPedInfo &mped_info, const LSFBase* param);
  bool process_skip_traits              (RefMPedInfo &mped_info, const LSFBase* param);
  bool process_skip_markers             (RefMPedInfo &mped_info, const LSFBase* param);
  bool process_dynamic_markers          (RefMPedInfo &mped_info, const LSFBase* param);
  bool process_reject_partial_lineage   (RefMPedInfo &mped_info, const LSFBase* param);
  bool process_individual_missing_value (RefMPedInfo &mped_info, const LSFBase* param);
  bool process_no_sex_ok                (RefMPedInfo &mped_info, const LSFBase* param);
};

/** \brief Reads in LSF parameters specific to fortran-formatted pedigree files
  *
  * Reads in LSF parameters specific to fortran-formatted pedigree files.
  */
class RefLSFFortranPedigreeFile : virtual public RefFortranPedigreeFile, virtual public RefLSFPedigreeFile
{
public:
  ///
  /// Constructor.
  RefLSFFortranPedigreeFile(cerrorstream &errors = sage_cerr);

  ///
  /// Processes an entire "PEDIGREE" block.
  /// \param mped_info The RefMPedInfo instance associated with the RefMultiPedigree in question
  /// \param params The LSFBase pointer corresponding to the pedigree block
  virtual bool process_parameters(RefMPedInfo &mped_info, const LSFBase *params);

  ///
  /// Processes a single parameter in the "PEDIGREE" block.
  /// \param mped_info The RefMPedInfo instance associated with the RefMultiPedigree in question
  /// \param param The LSFBase pointer corresponding to the parameter
  virtual bool process_parameter(RefMPedInfo &mped_info, const LSFBase *param);

private:

  virtual string get_marker_covariate_value(string m, const string& v1, const string& v2) { return ""; }
};

/** \brief Reads in LSF parameters specific to character delimeted pedigree files
  *
  * Reads in LSF parameters specific to character delimeted pedigree files.
  */
class RefLSFDelimitedPedigreeFile : virtual public RefDelimitedPedigreeFile, virtual public RefLSFPedigreeFile
{
public:
  ///
  /// Constructor.
  RefLSFDelimitedPedigreeFile(cerrorstream &errors = sage_cerr);
  
  virtual bool process_parameters(RefMPedInfo &mped_info, const LSFBase *params);
  ///
  /// Processes a single parameter in the "PEDIGREE" block.
  /// \param mped_info The RefMPedInfo instance associated with the RefMultiPedigree in question
  /// \param param The LSFBase pointer corresponding to the parameter
  virtual bool process_parameter(RefMPedInfo &mped_info, const LSFBase *param);

  virtual bool input(RefMultiPedigree &p, const std::string &filename, ostream &messages = cout);

private:

  struct MarkerListElement
  {
    std::string start_marker;
    std::string end_marker;
    
    LSF_ptr<LSFBase> marker_params;
  };

  std::list<MarkerListElement> my_marker_lists;

  std::list<MarkerListElement>::iterator find_marker_list(string start);

  bool build_marker_list_parameters (RefMPedInfo &mped_info);
  bool build_marker_list ( string_tokenizer& header,
                           string_tokenizer::iterator i, 
                           RefMPedInfo &mped_info,
                           std::list<MarkerListElement>::iterator mlist_index);

  std::list<MarkerListElement> my_covariate_lists;

  std::list<MarkerListElement>::iterator find_covariate_list(string start);

  bool build_covariate_list_parameters (RefMPedInfo &mped_info);
  bool build_covariate_list ( string_tokenizer& header,
                              string_tokenizer::iterator i, 
                              RefMPedInfo &mped_info,
                              std::list<MarkerListElement>::iterator mlist_index);

  bool process_format2                     (RefMPedInfo &mped_info, const LSFBase *param);
  bool process_whitespace                  (RefMPedInfo &mped_info, const LSFBase *param);
  bool process_delimiter_mode              (RefMPedInfo &mped_info, const LSFBase *param);
  bool process_delimiters                  (RefMPedInfo &mped_info, const LSFBase *param);
  bool process_skip_leading_delimiters     (RefMPedInfo &mped_info, const LSFBase *param);
  bool process_skip_trailing_delimiters    (RefMPedInfo &mped_info, const LSFBase *param);
  bool process_skip_consecutive_delimiters (RefMPedInfo &mped_info, const LSFBase *param);
  bool process_study_id                    (RefMPedInfo &mped_info, const LSFBase *param);
  bool process_pedigree_id                 (RefMPedInfo &mped_info, const LSFBase *param);
  bool process_individual_id               (RefMPedInfo &mped_info, const LSFBase *param);
  bool process_parent_id                   (RefMPedInfo &mped_info, const LSFBase *param);
  bool process_sex_field                   (RefMPedInfo &mped_info, const LSFBase *param);
  bool process_marker                      (RefMPedInfo &mped_info, const LSFBase *param);
  bool process_phenotype                   (RefMPedInfo &mped_info, const LSFBase *param);
  bool process_string                      (RefMPedInfo &mped_info, const LSFBase *param);
  bool process_marker_list                 (RefMPedInfo &mped_info, const LSFBase *param);
  bool process_covariate_list              (RefMPedInfo &mped_info, const LSFBase *param);

  bool test_skip_marker(const LSFBase* param) const;
  bool test_allow_dynamic(const LSFBase* param) const;
  size_t setup_dynamic_marker(RefMPedInfo &mped_info, const LSFBase* param, const string& marker_name);  
  void parse_marker_parameters(RefMPedInfo &mped_info,const LSFBase* param, size_t marker_id);

  // Added for marker_covariate, i.e., read in marker data as covariate.
  //
  struct marker_covariate_info
  {
    string  type;
    string  field_name;
    string  mcov_name;
    char    allele_delimiter;
    string  allele_missing;
    string  covariate_moi;
    string  covariate_allele;
    bool    allow_hemizygote;
    size_t  trait_index;
  };

  std::map<string, marker_covariate_info> my_marker_covariates;

  bool use_as_covariate(const LSFBase *param);
  void add_marker_covariate(RefMPedInfo &mped_info, const LSFBase *param,
                            const string& fname,    const string& mname);
  bool process_marker_covariates(RefMPedInfo &mped_info);
  virtual string get_marker_covariate_value(string m, const string& v1, const string& v2);
};

//============================
//  INLINE FUNCTIONS
//
//  RefPedigreeFile
//============================

inline cerrorstream RefPedigreeFile::error_sink() const        { return errors; }
inline void RefPedigreeFile::set_error_sink(cerrorstream &err) { errors = err;  }

inline size_t RefPedigreeFile::verbose_output         () const { return my_verbose_output;                              }
inline bool   RefPedigreeFile::reject_partial_lineage () const { return my_reject_partial_lineage;                      }
inline bool   RefPedigreeFile::require_record         () const { return my_require_record;                              }
inline bool   RefPedigreeFile::skip_traits            () const { return my_force_skip_traits     || my_skip_traits;     }
inline bool   RefPedigreeFile::skip_markers           () const { return my_force_skip_markers    || my_skip_markers;    }
inline bool   RefPedigreeFile::dynamic_markers        () const { return my_force_dynamic_markers || my_dynamic_markers; }
inline bool   RefPedigreeFile::force_skip_traits      () const { return my_force_skip_traits;                           }
inline bool   RefPedigreeFile::force_skip_markers     () const { return my_force_skip_markers;                          }
inline bool   RefPedigreeFile::force_dynamic_markers  () const { return my_force_dynamic_markers;                       }
inline bool   RefPedigreeFile::sex_linked_exist       () const { return my_sex_linked_exist;                            }
inline bool   RefPedigreeFile::build_incremental      () const { return my_build_incremental;                           }
inline bool   RefPedigreeFile::sex_code_trait         () const { return my_sex_code_trait;                              }
inline bool   RefPedigreeFile::pedigree_id_trait      () const { return my_pedigree_id_trait;                              }

inline const std::string & RefPedigreeFile::sex_field_name()   const { return my_sex_field_name;   }
inline const std::string & RefPedigreeFile::pedigree_id_name() const { return my_pedigree_id_name; }
inline const std::string & RefPedigreeFile::format()           const { return my_format;           }

inline void RefPedigreeFile::set_treat_as_sibs          (bool   t) { my_treat_as_sibs          = t; }
inline void RefPedigreeFile::set_no_sex_field           (bool   s) { my_no_sex_field           = s; }
inline void RefPedigreeFile::set_verbose_output         (size_t v) { my_verbose_output         = v; }
inline void RefPedigreeFile::set_reject_partial_lineage (bool   r) { my_reject_partial_lineage = r; }
inline void RefPedigreeFile::set_require_record         (bool   a) { my_require_record         = a; }
inline void RefPedigreeFile::set_skip_traits            (bool   s) { my_skip_traits            = s; }
inline void RefPedigreeFile::set_skip_markers           (bool   s) { my_skip_markers           = s; }
inline void RefPedigreeFile::set_dynamic_markers        (bool   d) { my_dynamic_markers        = d; }
inline void RefPedigreeFile::set_force_skip_traits      (bool   s) { my_force_skip_traits      = s; }
inline void RefPedigreeFile::set_force_skip_markers     (bool   s) { my_force_skip_markers     = s; }
inline void RefPedigreeFile::set_force_dynamic_markers  (bool   d) { my_force_dynamic_markers  = d; }
inline void RefPedigreeFile::set_sex_linked_exist       (bool   s) { my_sex_linked_exist       = s; }
inline void RefPedigreeFile::set_build_incremental      (bool   b) { my_build_incremental      = b; }
inline void RefPedigreeFile::set_sex_code_trait         (bool   a) { my_sex_code_trait         = a; }
inline void RefPedigreeFile::set_pedigree_id_trait      (bool   a) { my_pedigree_id_trait      = a; }
inline void RefPedigreeFile::set_no_sex_ok_option       (bool   a) { my_no_sex_ok_option       = a; }

inline void RefPedigreeFile::set_sex_field_name   (const std::string &f) { my_sex_field_name = f;   }
inline void RefPedigreeFile::set_pedigree_id_name (const std::string &f) { my_pedigree_id_name = f; }
inline void RefPedigreeFile::set_format           (const std::string &f) { my_format = f;           }

inline void RefPedigreeFile::validate()    { my_valid = true;  }
inline void RefPedigreeFile::invalidate()  { my_valid = false; }
inline bool RefPedigreeFile::valid() const { return my_valid;  }

inline RefPedigreeFile::field::field(field_t t, const std::string &f, const std::string &n, size_t i) : 
  type(t), field_name(f), name(n), index(i) 
{ }

inline void RefPedigreeFile::field::invalidate()
{
  index = (size_t)-1;
  type  = skip;
}

inline size_t RefPedigreeFile::field_count          () const { return my_fields.size();        }
inline size_t RefPedigreeFile::trait_count          () const { return my_trait_count;          }
inline size_t RefPedigreeFile::invalid_trait_count  () const { return my_invalid_trait_count;  }
inline size_t RefPedigreeFile::string_count         () const { return my_string_count;         }
inline size_t RefPedigreeFile::invalid_string_count () const { return my_invalid_string_count; }
inline size_t RefPedigreeFile::marker_count         () const { return my_marker_count;         }
inline size_t RefPedigreeFile::invalid_marker_count () const { return my_invalid_marker_count; }
inline size_t RefPedigreeFile::marker_cov_count         () const { return my_marker_cov_count;         }
inline size_t RefPedigreeFile::invalid_marker_cov_count () const { return my_invalid_marker_cov_count; }
inline size_t RefPedigreeFile::sex_count            () const { return my_sex_count;            }
inline size_t RefPedigreeFile::skip_count           () const { return my_skip_count;           }
inline size_t RefPedigreeFile::study_id_count       () const { return my_study_id_count;       }
inline size_t RefPedigreeFile::pedigree_id_count    () const { return my_pedigree_id_count;    }
inline size_t RefPedigreeFile::individual_id_count  () const { return my_individual_id_count;  }
inline size_t RefPedigreeFile::parent_id_count      () const { return my_parent_id_count;      }

inline bool RefPedigreeFile::get_treat_as_sibs    () const { return my_treat_as_sibs;    }
inline bool RefPedigreeFile::get_no_sex_field     () const { return my_no_sex_field;     }
inline bool RefPedigreeFile::get_no_sex_ok_option () const { return my_no_sex_ok_option; }

inline const RefPedigreeFile::field_list_type & RefPedigreeFile::field_list() const { return my_fields; }
inline       RefPedigreeFile::field_list_type & RefPedigreeFile::field_list()       { return my_fields; }

inline const vector<pair<std::string, std::string> > & RefPedigreeFile::get_ind_list() const { return my_ind_list; }

//===============================
//  INLINE FUNCTIONS
//
//  RefFortranPedigreeFile
//===============================

inline void RefFortranPedigreeFile::add_skip_field    () { field_list().push_back( field(skip));          }
inline void RefFortranPedigreeFile::add_study_id      () { field_list().push_back( field(study_id));      }
inline void RefFortranPedigreeFile::add_pedigree_id   () { field_list().push_back( field(pedigree_id));   }  
inline void RefFortranPedigreeFile::add_individual_id () { field_list().push_back( field(individual_id)); }  
inline void RefFortranPedigreeFile::add_parent_id     () { field_list().push_back( field(parent_id));     }
inline void RefFortranPedigreeFile::add_sex_field     () { field_list().push_back( field(sex_code));      }

inline void RefFortranPedigreeFile::add_trait_field  (const std::string & trait_name) { field_list().push_back(field(trait,        trait_name,  trait_name));  }
inline void RefFortranPedigreeFile::add_string_field (const std::string & name)       { field_list().push_back(field(string_field, name,        name));        }
inline void RefFortranPedigreeFile::add_marker_field (const std::string & marker_name){ field_list().push_back(field(marker,       marker_name, marker_name)); }
inline void RefFortranPedigreeFile::add_allele_field (const std::string & marker_name){ field_list().push_back(field(allele,       marker_name, marker_name)); }

//==========================
//  INLINE FUNCTIONS
//
//  RefDelimitedPedigreeFile
//==========================

inline bool RefDelimitedPedigreeFile::format_in_file() const     { return my_format_in_file; }

inline const std::string & RefDelimitedPedigreeFile::delimiters() const { return my_delimiters; }
inline const std::string & RefDelimitedPedigreeFile::whitespace() const { return my_whitespace; }

inline bool RefDelimitedPedigreeFile::skip_consecutive_delimiters () const { return my_skip_consecutive_delimiters; }
inline bool RefDelimitedPedigreeFile::skip_leading_delimiters     () const { return my_skip_leading_delimiters;     }
inline bool RefDelimitedPedigreeFile::skip_trailing_delimiters    () const { return my_skip_trailing_delimiters;    }

inline void RefDelimitedPedigreeFile::set_format_in_file              (bool f)    { my_format_in_file              = f;    }
inline void RefDelimitedPedigreeFile::set_skip_consecutive_delimiters (bool skip) { my_skip_consecutive_delimiters = skip; }
inline void RefDelimitedPedigreeFile::set_skip_leading_delimiters     (bool skip) { my_skip_leading_delimiters     = skip; }
inline void RefDelimitedPedigreeFile::set_skip_trailing_delimiters    (bool skip) { my_skip_trailing_delimiters    = skip; }

inline const RefDelimitedPedigreeFile::field_map_type & RefDelimitedPedigreeFile::field_map() const { return my_field_map; }
inline       RefDelimitedPedigreeFile::field_map_type & RefDelimitedPedigreeFile::field_map()       { return my_field_map; }

inline void RefDelimitedPedigreeFile::add_study_id_field      (const std::string &field_name)                                  { my_field_map[toUpper(field_name)] = field(study_id,      field_name);              }
inline void RefDelimitedPedigreeFile::add_pedigree_id_field   (const std::string &field_name)                                  { my_field_map[toUpper(field_name)] = field(pedigree_id,   field_name);              }
inline void RefDelimitedPedigreeFile::add_individual_id_field (const std::string &field_name)                                  { my_field_map[toUpper(field_name)] = field(individual_id, field_name);              }
inline void RefDelimitedPedigreeFile::add_parent_id_field     (const std::string &field_name)                                  { my_field_map[toUpper(field_name)] = field(parent_id,     field_name);              }
inline void RefDelimitedPedigreeFile::add_sex_field           (const std::string &field_name)                                  { my_field_map[toUpper(field_name)] = field(sex_code,      field_name);              }
inline void RefDelimitedPedigreeFile::add_trait_field         (const std::string &field_name, const std::string & trait_name)  { my_field_map[toUpper(field_name)] = field(trait,         field_name, trait_name);  }
inline void RefDelimitedPedigreeFile::add_string_field        (const std::string &field_name, const std::string & string_name) { my_field_map[toUpper(field_name)] = field(string_field,  field_name, string_name); }
inline void RefDelimitedPedigreeFile::add_allele_field        (const std::string &field_name, const std::string & marker_name) { my_field_map[toUpper(field_name)] = field(allele,        field_name, marker_name); }
inline void RefDelimitedPedigreeFile::add_marker_field        (const std::string &field_name, const std::string & marker_name) { my_field_map[toUpper(field_name)] = field(marker,        field_name, marker_name); }
inline void RefDelimitedPedigreeFile::add_allele_cov_field    (const std::string &field_name, const std::string & mcov_name)   { my_field_map[toUpper(field_name)] = field(allele_cov,    field_name, mcov_name);   }
inline void RefDelimitedPedigreeFile::add_marker_cov_field    (const std::string &field_name, const std::string & mcov_name)   { my_field_map[toUpper(field_name)] = field(marker_cov,    field_name, mcov_name);   }

} // End namespace RPED
} // End namespace SAGE

#endif
