#ifndef MEMBER_DATA_SAMPLE_H
#define MEMBER_DATA_SAMPLE_H

#include "output/Output.h"
#include "error/internal_error.h"
#include "fped/fped.h"
#include "mped/mp_utilities.h"
#include "sampling/Field.h"
#include "sampling/IndividualValidator.h"
#include "sampling/TraitSet.h"

namespace SAGE {
namespace SAMPLING {

class TraitValueCalculator;

/** Assists in the importation of data from a RefMultiPedigree for use in a particular analysis.
  *
  * \par Purpose
  *
  * The MemberDataSample object assists in the importation of data from a RefMultiPedigree for use in a particular analysis. 
  * Furthermore, the MemberDataSample object is \b generalized; that is, it is a library-level feature that generalizes a 
  * piece of functionality that currently exists in proprietary, program-specific form in all SAGE programs.
  *
  * When a SAGE program prepares data for use, it takes a RefMultiPedigree as input. Then, based on information from 
  * the program's parameter file, it pulls out and transforms the RefMultiPedigree data into a form more usable by the 
  * program. Consider, for instance, the following parameter file:
  *
  * pedigree
  * {
  *  ...
  *  trait=sqrtdbh,missing=-1
  *  ...
  * }
  *
  * assoc_analysis
  * {
  *  ...
  *   primary_trait=sqrtdbh
  * }
  *
  * When assoc is run, the multipedigree library first processes the pedigree block and pulls in the indicated 
  * information from the data file. "sqrtdbh", for example, is read from the data file and stored as a continuous trait.
  *
  * Then, once the RefMultiPedigree is passed on to assoc, assoc in turn consults its analysis parameters and pulls in 
  * and translates required data from the RefMultiPedigree. In this case, assoc creates an assoc_sampledata object, where the 
  * get_primary_trait() function returns values gleaned from the RefMultiPedigree's "sqrtdbh" data.
  *
  * Although the data importation process works fine as things currently stand, it is implemented on a proprietary, 
  * program-by-program basis. We need to generalize that procedure, so that data importation can be standardized, and 
  * suite-wide features applied more easily.
  *
  * \par Basic interface
  *
  * An MemberDataSample has several components to its interface:
  *
  * (1) Data importation - functions that tell the MemberDataSample which data to import from the RefMultiPedigree and 
  * how to organize it.
  *
  * (2) Data finalization - functions that tell the MemberDataSample how to process the imported data
  *
  * (3) Data creation - functions for creating new fields not directly from the source data
  *
  * (4) User-created data finalization - functions that tell the MemberDataSample how to process the user-created data
  *
  * (5) Data extraction - functions that return requested sample values
  *
  * \par Data importation
  *
  * Data importation is accomplished through the importField() function. For every trait/covariate in your
  * RefMultiPedigree that you want imported into the sample for your analysis, you should use the importField()
  * function. If, for instance, the analysis' parameter specs look like this:
  *
  * assoc
  * {
  *  ...
  *    primary_trait = sqrtdbh
  *   ...
  * }
  *
  * Then your invocation of importField() would look something like this:
  * 
  * FieldID primary_trait_id = my_analysis_sample.importField("sqrtdbh", "general", "primary_trait");
  *
  * \par Finalizing data
  *
  * Once you have imported all the fields you need for a particular analysis, you should invoke finalizeData().
  * This causes the MemberDataSample object to review all individual data, process uninformative individuals, and
  * adjust field data according to each Field's import flags.
  *
  * The finalization process allows for an individual \b validation algorithm to be used. Although validity has
  * an analysis-specific meaning, the basic concept holds for all analyses: a "valid" individual has trait values
  * in the proper range for conducting an analysis.
  *
  * When finalizeData() is invoked, it validates each individual. If you do not specify your own validator, finalizeData()
  * will invalidate any individual with a missing value for any single trait. "Invalidation" means that that
  * individual's trait values will \b all be replaced with missing values. It is important to remember that invalidation
  * does not simply flag an individual as uninformative--it wipes out all trait values for that individual.
  *
  * It might occur to you at this point to wonder why it is not sufficient to simply flag the individual as
  * uninformative and leave it at that. The reason is because the post-importation Field adjustment depends on a
  * completely verified data set. If a field is to be mean-adjusted, then the correct mean must be calculated with
  * which to make the adjustment. Said mean must only be derived from the values of valid individuals. Consider the
  * following dataset:
  *
  * IND         COV1          COV2
  * 1           1.0           2.5
  * 2           1.1           2.4
  * 3           <missing>     2.0
  * 4           1.2           2.1
  * 5           0.9           <missing>
  * 6           0.8           2.2
  * 7           1.0           2.1
  *
  * If the Field "COV1" is supposed to be mean adjusted, then the mean for this Field must only include fully valid
  * individuals. If COV1's mean were calculated using the above data as it is shown, then COV1's mean would equal 1.0.
  * After the default validation algorithm is applied, the dataset looks like this:
  *
  * IND         COV1          COV2
  * 1           1.0           2.5
  * 2           1.1           2.4
  * 3           <missing>     <missing>
  * 4           1.2           2.1
  * 5           <missing>     <missing>
  * 6           0.8           2.2
  * 7           1.0           2.1
  *
  * Note that individuals #3 and #5 have been stripped of all their trait values. This is because, for this analysis,
  * any missing value for an individual renders that individual fully uninformative. With the "corrected" sample, the
  * mean of COV1 now equals 1.02. This is the \b correct mean.
  *
  * \par Individual validity
  *
  * It should be noted that unless you use the RPEDNEW::Field::ALLOW_AVERAGING flag when you import your data
  * (see RPEDNEW::Field::ImportFlagsEnum), then any individual with \b any single missing value for \b any field
  * will be considered invalid. When finalizeData() is invoked (after all Field's have been imported), the
  * MemberDataSample object reviews the data for all individuals. Any individual with missing data in a field imported
  * \b without the allow-averaging flag will have all of its field values marked as uninformative. Even if the
  * individual has informative values for other fields, those values will be changed to uninformative. Thus, the
  * mean and variance adjustment applied to those parameters will not reflect the individual's removed values.
  *
  * \par Trait presence
  *
  * You may find it useful at some point in your calculations to see if an individual actually \b has
  * a value for some trait. The easy way to do this is with SAGE::isnan() :
  *
  * if(SAGE::isnan(sample.getAdjValue(x, "group", "foo")) // ... do something
  *
  * If you want to avoid the unpleasant overhead of using SAGE::isnan(), however, you can simply use the 
  * Field::isAdjValuePresent() function:
  *
  * if(sample.getField("group", "foo").isAdjValuePresent(x))) // ... do something
  *
  * \par Data creation
  *
  * It is possible to create sample-based fields whose values are calculated by your program. For instance, if the
  * original sample includes a trait indicating height, you may wish to create a field indicating whether a value
  * for height is present for an individual in the sample. That is, the desired trait equals \c true if height is
  * indicated, \c false otherwise. (Please note that the absence of a source-indicated value is indicated by QNAN).
  *
  * Please note that you can only create new fields \b after you have finalized the imported data.
  *
  * The first step, obviously, is to import the height attribute into your sample and finalize it:
  *
  * sample.importField("height", "general", "height");
  * sample.finalizeData();
  *
  * Now, you should create a field entry for the presence of a height attribute:
  *
  * sample.createField("general", "height_present");
  *
  * Next, you populate the field with values:
  *
  * for(size_t i = 0; i < sample->getTotalIndividualCount(); ++i)
  * {
  *   double height = sample->getAdjValue(i, "general", "height");
  *   sample->getField("general", "height").setOrigValue(i, !SAGE::isnan(height));
  * }
  *
  * Lastly, you finalize the user created fields:
  *
  * sample->finalizeUserCreatedField("general", "height");
  *
  * \par ACCESSING SAMPLING
  *
  * Accessing imported data is relatively straightforward. To access any individual's field data, you can
  * use the isValid() and getAdjValue() functions, where the individual is identified by its 
  * RefMultiPedigree::member_type.
  *
  * You can also access the Field's themselves with the getField function, where Field's are identified either
  * by FieldID, group name and GroupID, or group name and field name.
  *
  * Lastly, if you need to iterate across Field (for a composite calculation, for instance), you can you the
  * getFieldBegin() and getFieldEnd() functions.
  */
class MemberDataSample
{
public:

  friend class PartitionedMemberDataSample;

  // Indicates the status of this MemberDataSample.
  enum ImportStatusEnum
  {
    READY_FOR_IMPORT_AND_FINALIZE = 0,
      /* User may invoke the importField() function zero or more times, then invoke
       * the finalizeData() function \b once.
       */
       
    READY_FOR_CREATE_AND_FINALIZE = 1,
      /* User may invoke the createField() function zero or more times, then invoke
       * the finalizeUserCreatedData() function \b once.
       */
       
    COMPLETELY_FINALIZED = 2
      /* User may only use const functions, and cannot use importField() or createField().
       */
  };

    MemberDataSample(const FPED::Multipedigree& rmp, cerrorstream& err = sage_cerr);
    MemberDataSample(const MemberDataSample& other); 
    ~MemberDataSample();
  
  // Data importation
    const FPED::Multipedigree& getMultipedigree() const;

     /* Resets the MemberDataSample object (ie: deletes all field entries).
      *
      * If you are using the same instance of an MemberDataSample object for multiple analyses, then you should
      * make sure to clear & setup the MemberDataSample prior to each analysis. Thus, your analysis control loop
      * should take the form
      * 
      * for cur_analysis in analysis:
      *   ...
      *   my_analysis_sample.reset();    // This is supposed to describe clearData()??? -djb
      *
      *   for cur_trait in cur_analysis:
      *     my_analysis_sample.importField(cur_trait..., ...);
      * 
      * \retval 0 Data cleared successfully.
      * \retval 1 Data \b not cleared successfully.
      */
    bool clearData();
    void  reset();

    bool addGroup(const string& group_name);

    // Imports a trait/covariate from a data source (the RefMultiPedigree indicated at construction) into the 
    // MemberDataSample object.
    //
    // Please note that this function only operates if the object has not been finalized (see finalizeData()).
    //
    // \param source_trait_name The name of trait as it appears in the RefMultiPedigree
    // \param group_name The name of the field group to which the new field will belong
    // \param new_field_name The name of the new field
    // \param flags Control flags for customizing the import process. See Field::ImportFlagsEnum for more information.
    // \retval The unique FieldID of the newly created Field
    const FieldID& importField(const string& source_trait_name, const string& group_name,
                                const string& new_field_name, unsigned long flags = Field::NO_FLAGS);

    // Flags the indicated individual as invalid.
    //
    // Please note that this function only operates if the object has not been finalized (see finalizeData()).
    //
    // \param i The id of the member to be flagged as invalid.
    // \retval 0 Successful
    // \retval 1 Unsuccessful
    int invalidateIndividual(size_t i);

  // Data creation
    //
    // Creates an entry for a field with the indicated configuration options.
    // \param group_name The name of the group to which the new field will belong
    // \param field_name The name of the new field
    // \param flags Control flags for customizing the finalization process. See Field::ImportFlagsEnum for more information.
    // \retval The unique FieldID of the newly created Field
    const FieldID& createField(string group_name, string field_name, unsigned long flags = Field::NO_FLAGS);

    // Populates a user-defined Field with values taken from the given functor.
    // For more information, see TraitValueCalculator.
    //
    // \param group_name The name of the group to which the Field belongs
    // \param field_name The name of the field
    // \param calc A TraitValueCalculator designed to calculate the field's value for a given individual
    void populateUserCreatedField(const string& group_name, const string& field_name, const TraitValueCalculator & calc);

  // Data finalization
    
    ImportStatusEnum getImportStatus() const;

    // Post-processes the imported data, checking for individual validity, centering/standardizing 
    // data as requested, and filling in missing data (as requested).
    //
    // This version of finalizeData() invalidates any single individual if that individual has \b any missing
    // values for any fields. If you want to specify your own validation scheme, create an invalidator class
    // (see validation section in detailed description) and use the appropriate version of finalizeData().
    //
    // Please note that this function only operates if the object has not been finalized (see finalizeData()).
    //
    // After you have imported all the fields you need (through importField()), you need to invoke
    // finalizeData().
    // \retval 0 Successful
    // \retval 1 Unsuccessful
    int finalizeData();

    // Post-processes the imported data, checking for individual validity, centering/standardizing 
    // data as requested, and filling in missing data (as requested).
    //
    // Please note that this function only operates if the object has not been finalized (see finalizeData()).
    //
    // After you have imported all the fields you need (through importField()), you need to invoke
    // finalizeData().
    // \param invalidator Any class with the validation interface defined (see detailed notes on this MemberDataSample)
    // \retval 0 Successful
    // \retval 1 Unsuccessful
    int finalizeData(const IndividualValidator& invalidator);

    // Post-process any fields created via the createField() functions. This will \b NOT alter any aspects of
    // individual validity; this function simply applies the adjustment options requested on the field.
    int finalizeUserCreatedData();

  // Individual counts

    size_t getTotalIndividualCount() const;

    // Please note that this function only operates if the object has been finalized (see finalizeData()).
    //
    size_t getValidIndividualCount() const;

    // Please note that this function only operates if the object has been finalized (see finalizeData()).
    //
    size_t getInvalidIndividualCount() const;
    size_t getValidSingletonCount() const;
    size_t getInvalidSingletonCount() const;


  // Individual information

    const FPED::Member& getIndividual(size_t i) const;

    // Indicates whether or not the individual is considered valid.
    //
    // Please note that this function only operates if the object has been finalized (see finalizeData()).
    bool isValid(size_t i) const;

    // Returns the individual's value for the indicated field.
    // \param i The individual in question
    // \param field_id The FieldID of the Field in question
    double getOrigValue(size_t i, const FieldID& field_id) const;

    // Returns the individual's value for the indicated field.
    // \param i The individual in question
    // \param group_name The name of the group to which the Field belongs
    // \param group_id The GroupID of the Field in question
    double getOrigValue(size_t i, const string& group_name, const GroupID& group_id) const;

    // Returns the individual's value for the indicated field.
    // \param i The individual in question
    // \param group_name The name of the group to which the Field belongs
    // \param field_name The name of the Field in question
    double getOrigValue(size_t i, const string& group_name, const string& field_name) const;

    // Returns the individual's value for the indicated field.
    // \param i The individual in question
    // \param field_id The FieldID of the Field in question
    double getReplacedValue(size_t i, const FieldID & field_id) const;

    // Returns the individual's value for the indicated field.
    // \param i The individual in question
    // \param group_name The name of the group to which the Field belongs
    // \param group_id The GroupID of the Field in question
    double getReplacedValue(size_t i, const string& group_name, const GroupID& group_id) const;

    // Returns the individual's value for the indicated field.
    // \param i The individual in question
    // \param group_name The name of the group to which the Field belongs
    // \param field_name The name of the Field in question
    double getReplacedValue(size_t i, const string& group_name, const string& field_name) const;

    // Returns the individual's value for the indicated field.
    // \param i The individual in question
    // \param field_id The FieldID of the Field in question
    double getAdjValue(size_t i, const FieldID& field_id) const;

    // Returns the individual's value for the indicated field.
    // \param i The individual in question
    // \param group_name The name of the group to which the Field belongs
    // \param group_id The GroupID of the Field in question
    double getAdjValue(size_t i, const string& group_name, const GroupID& group_id) const;

    // Returns the individual's value for the indicated field.
    // \param i The individual in question
    // \param group_name The name of the group to which the Field belongs
    // \param field_name The name of the Field in question
    double getAdjValue(size_t i, const string & group_name, const string & field_name) const;

  // Field counts

    // Returns the total number of fields present in the MemberDataSample.
    size_t getFieldCount() const;

    // Returns the total number of fields present in the named group.
    // \param group_name The name of the group whose number of fields will be returned
    size_t getFieldCount(const string & group_name) const;

    // Returns the number of field groups present in the MemberDataSample.
    size_t getGroupCount() const;

  // Field extraction

     // Returns \c true if the named parameter exists, \c false if it does not.
     // \param group_name The name of the group in question
     // \param field_name The name of the field in question
     bool fieldExists(const string& group_name, const string& field_name) const;
     
     // Returns \c true if the named parameter exists, \c false if it does not.
     // \param group_name The name of the group in question
     // \param group_id The group_id of the Field in question
     bool fieldExists(const string& group_name, const GroupID& group_id) const;

     // Returns a non-const reference to the numbered field.
     // \param id The FieldID that identifies the field.
     Field& getField(const FieldID& id);

     // Returns a non-const reference to the named & numbered field.
     // \param group_name The name of the group to which the field belongs
     // \param group_id The GroupID that identifies the field within the named group
     Field& getField(const string& group_name, const GroupID& group_id);

     // Returns a non-const reference to the named field.
     // \param group_name The name of the group to which the field belongs
     // \param field_name The name of the field
     Field& getField(const string& group_name, const string& field_name);

     // Returns a const reference to the numbered field.
     // \param id The FieldID that identifies the field.
     const Field& getField(const FieldID& id) const;

     // Returns a const reference to the named field.
     // \param group_name The name of the group to which the field belongs
     // \param field_name The name of the field
     const Field& getField(const string& group_name, const string& field_name) const;

     // Returns a const reference to the named & numbered field.
     // \param group_name The name of the group to which the field belongs
     // \param group_id The GroupID that identifies the field within the named group
     const Field& getField(const string& group_name, const GroupID& group_id) const;


  // Group extraction
  
     // Returns \c true if the named group exists, \c false if it does not.
     // \param group_name The name of the group in question
     bool groupExists(const string& group_name) const;

     // Returns the named field group.
     //
     // \param group_name The name of the requested group
     const FieldGroupType& getGroup(const string& group_name) const;

     // Returns the named field group.
     //
     // \param group_name The name of the requested group
     FieldGroupType& getGroup(const string& group_name);


  // name Field traversal

     // Returns a non-const begin iterator for all Field's.
     FieldIterator getFieldBegin();

     // Returns a non-const end iterator for all Field's.
     FieldIterator getFieldEnd();

     // Returns a non-const begin iterator for the named group.
     // \param group_name The name of the group whose iterator will be returned.
     FieldIterator getFieldBegin(const string& group_name);

     // Returns a non-const end iterator for the named group.
     // \param group_name The name of the group whose iterator will be returned.
     FieldIterator getFieldEnd(const string& group_name);

     // Returns a const begin iterator for all Field's.
     FieldConstIterator getFieldBegin() const;

     // Returns a const end iterator for all Field's.
     FieldConstIterator getFieldEnd() const;

     // Returns a const begin iterator for the named group.
     // \param group_name The name of the group whose iterator will be returned.
     FieldConstIterator getFieldBegin(const string& group_name) const;

     // Returns a const end iterator for the named group.
     // \param group_name The name of the group whose iterator will be returned.
     FieldConstIterator getFieldEnd(const string& group_name) const;


  // Output formatting
  
    // Returns summary information about the number of individuals as an
    // OUTPUT::Table instance.
    const OUTPUT::Table& getSummaryTable() const { return my_summary_table; }
  
  // Debugging

    // Dumps a table of all trait values to cout
    void dumpTraitValues() const;
    
protected:
    cerrorstream& getErrorstream();

private:

    // Adds a new field with the indicated configuration options.  Used by
    // importField() and createField()
    //
    // \param source_trait_name The name of trait as it appears in the RefMultiPedigree
    // \param group_name The name of the group to which the new field will belong
    // \param field_name The name of the new field
    // \param flags Control flags for customizing the finalization process. See Field::ImportFlagsEnum for more information.
    // \param user_defined Is the new Field user-defined?
    // \retval The unique FieldID of the newly created Field
    const FieldID& addField(const string& source_trait_name, const string& group_name,
                            const string& new_field_name, bool user_defined, unsigned long flags = Field::NO_FLAGS);

    void validateIndividuals(const IndividualValidator& invalidator);
    void updateCounts();

    MemberDataSample& operator=(const MemberDataSample&); 

  // Data members:

    // - Affected by clearData()
    //
    FieldVectorType                            my_field_vector;    // The actual data    
    std::map<pair<string, string>, FieldID>    my_field_id_table;  // Lookup table match group/field name to FieldID    
    mutable std::map<string, FieldGroupType>   my_groups;          // Field groups    
    ImportStatusEnum                           my_import_status;
    
    // - Affected by updateCounts()
    //
    // - Note a valid individual has no missing data for the main phenotype or any of the
    //   covariates specified in any of the models for an analysis.
    //
    int  my_total_individual_count;
    int  my_valid_individual_count;    // Includes both singletons and subpedigree members.
    int  my_invalid_individual_count;  // Includes both singletons and subpedigree members.
    int  my_valid_singleton_count;        // Added 6-11-7. djb
    int  my_invalid_singleton_count;      // Added 6-11-7. djb    
    
    vector<bool>  my_valid;            // individual validity's
    OUTPUT::Table  my_summary_table;        
    
    FieldID          my_bad_field_id;               // for returning failed addField operations
    FieldGroupType   my_bad_group;                  // for returning failed getGroup operations    

    mutable cerrorstream  my_errors;           
    const   FPED::Multipedigree&  my_mp;
};





// \brief Functor for quickly calculating user-defined trait values
//
// \par Introduction
//
// If you want to create a user-defined trait, it's relatively simple.
// First, you create the trait itself (with MemberDataSample::createTrait).
// Then, you loop through each individual and populate it with a value.
// For instance:
// \code
// sample.createField("general", "height_present");
//
// for(size_t i = 0; i < sample->getTotalIndividualCount(); ++i)
// {
//   double height = sample->getAdjValue(i, "general", "height");
//   sample->getField("general", "height").setOrigValue(i, !SAGE::isnan(height));
// }
// \endcode
//
// This works fine, but wouldn't it be nice if you could simply write a 
// functor to carry out the calculation, and then pass have the sampling
// library do the work for you?
//
// You can!
//
// With the TraitValueCalculator, you can write a functor to easily calculate
// user-defined traits. For instance, let's say you wanted to create some
// user-defined trait "bmi", where "BMI = weight / (height*height)".
//
// To do this with the TraitValueCalculator, you simply derived your own
// calculator class from the TraitValueCalculator, define the operator(), and
// voila!
//
// \code
// class BMICalculator : public SAMPLING::TraitValueCalculator 
// {
//   virtual double operator() (size_t i, const SAMPLING::MemberDataSample & sample) const 
//   { 
//     double weight = sample.getAdjValue(i, "general", "weight"),
//            height = sample.getAdjValue(i, "general", "height");
//            
//     return weight / (height * height);
//   }    
// }    
//
// sample.createField("special", "BMI");
// sample.populateUserCreatedField("special", "valid", BMICalculator());
// \endcode
class TraitValueCalculator
{
public:
  // Required to make compiler happy.
  virtual ~TraitValueCalculator() { }
  virtual double operator() (size_t i, const MemberDataSample& sample) const = 0;
};


} // End namespace SAMPLING
} // End namespace SAGE

#include "sampling/MemberDataSample.ipp"

#endif
