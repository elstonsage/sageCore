#ifndef FIELD_H
#define FIELD_H
//=====================================================
//
//  File:	Field.h
//
//  Author:	Stephen Gross
//
// Copyright 2004 R.C. Elston
//=====================================================


#include <string>
#include <vector>
#include "numerics/sinfo.h"

namespace SAGE     {
namespace SAMPLING {

class Field;
class FieldIterator;
class FieldConstIterator;

/** \class FieldID
  * \brief Uniquely identifies a field in an MemberDataSample.
  */
class FieldID
{
public:
  friend class Field;
  friend class FieldIterator;
  friend class FieldConstIterator;
  friend class MemberDataSample;
  friend class PartitionedMemberDataSample;
     
  ///
  /// Constructor
  FieldID() {}

  ///
  /// Copy constructor
  FieldID(const FieldID & other);
        
  ///
  /// Assignment operator
  FieldID& operator= (const FieldID & other);
  
  ///
  /// Returns \c true if it is a valid id, \c false otherwise
  bool isValid() const { return id != (size_t) -1; } 
  
  ///
  /// Sets the id to be invalid
  void invalidate() { id = (size_t) -1; }
  
  ///
  /// Compares two ids and says if they're equal
  bool operator==(const FieldID& f) { return f.id == id; }
     
private:
  FieldID(size_t value) : id(value) {}
  size_t id;
};

/** \class GroupID
  * \brief Uniquely identifies a field within a field group in an MemberDataSample
  */
class GroupID
{
public:
  friend class Field;
  friend class FieldIterator;
  friend class FieldConstIterator;
  friend class MemberDataSample;
  
  ///
  /// Copy constructor
  GroupID(const GroupID & other);
  
  ///
  /// Assignment operator
  GroupID& operator= (const GroupID & other);
  
private:
  GroupID() {}
  GroupID(size_t value) : id(value) {}
  size_t id;
};

/** \class Field
 *  \brief Provides descriptive information about a set of continuous values.
 *
 *  \par INTRODUCTION
 *
 *  Field stores imported values from a RefMultiPedigree. Through the MemberDataSample,
 *  Field is capable of changing the imported values (through missing-value replacement,
 *  mean adjustment, and std. dev. adjustment).
 *
 *  \par DESCRIPTIVE STATISTICS
 *
 *  The following information regarding the imported data is available:
 *
 *  Mean
 *
 *  Variance
 *
 *  Standard deviation
 *
 */
class Field
{
public:

  friend class MemberDataSample;
  friend class PartitionedMemberDataSample;
  friend class std::vector<Field>;
    
  ///
  /// Bitwise flags for use in controlling the importation process.
  enum ImportFlagsEnum
  {
    NO_FLAGS    = 0,
          /*!< \hideinitializer
           * No flags enabled.
           */
    MEAN_ADJUST = 1,
          /*!< \hideinitializer
           * Flags all values for mean adjustment when finalizeData() is invoked.
           * ( \f$ \hat x_i = x_i - \bar X \f$ )
           */
    STDEV_ADJUST = 2,
          /*!< \hideinitializer
           * Flags all values for standard deviation adjustment when finalizeData() is invoked.
           * ( \f$ \hat x_i = { x_i \over \sigma{X} } \f$ )
           */
    ALLOW_AVERAGING = 4
          /*!< \hideinitializer
           * Allows missing values to be replaced with the field's mean.
           */
  }; 

  /// @name Structs
  //@{
  
    struct SummaryInfo
    {
      std::string name;
      double mean;
      double stdev;
      double min;
      double max;
    };
  
  //@}

  /// @name Copy constructor & operators
  //@{
  
    ///
    /// Copy constructor.
    Field(const Field & other);
  
    ///
    /// Assignment operator.
    Field& operator= (const Field & other);

  //@}

  /// @name Modifying field contents
  //@{

    ///
    /// Sets the i'th (original, non-finalized) value to val.
    /// \param i The index number of the element to be changed (the member's mpindex())
    /// \param val The value to which the element will be changed.
    /// \retval 0 If successful.
    /// \retval 1 If \b not successful.
    int setOrigValue(size_t i, double val);

    ///
    /// Non-const version of the getValues() function.
    std::vector<double> & getAllOrigValues();

  //@}

  /// @name Basic information
  //@{

    ///
    /// Returns a SummaryInfo instance describing this field.
    SummaryInfo getSummaryInfo() const
    {
      SummaryInfo i;
      
      i.name  = getFieldName  ();
      i.mean  = getMean       ();
      i.stdev = getStdev      ();
      i.min   = getSampleInfo ().min();
      i.max   = getSampleInfo ().max();
      
      return i;
    }

    ///
    /// Returns the name of the trait originally used to read in data for this field.
    const std::string & getOrigName() const;

    ///
    /// Returns the name of the field.
    const std::string & getFieldName() const;
    
    ///
    /// Returns the name of the group to which the field belongs.
    const std::string & getGroupName() const;

    ///
    /// Returns this Field's FieldID.
    const FieldID & getFieldID() const;

    ///
    /// Returns this Field's GroupID.
    const GroupID & getGroupID() const;

    ///
    /// Returns the flags that were used to control the post-import adjustment of this field.
    unsigned long getInputFlags() const;

    ///
    /// If this field was created via MemberDataSample::createField(), then this function will return true.
    /// Otherwise, it returns false.
    bool getUserDefined() const;

  //@}

  /// @name Querying specific values
  //{

    ///
    /// Returns the original value for the i'th member (as read from the RefMultiPedigree)
    /// \param i The member's mpindex
    double getOrigValue(size_t i) const;

    ///
    /// Returns the replaced value for the i'th member (as read from the RefMultiPedigree)
    /// This will be simply the "orig" value if no mean replacement was requested in the 
    /// import flags; otherwise, it will be the mean of the orig values.
    /// \param i The member's mpindex
    double getReplacedValue(size_t i) const;

    ///
    /// Returns the adjusted (post-finalization) value for the i'th member.
    /// \param i The member's mpindex
    double getAdjValue(size_t i) const;

    ///
    /// Returns the original values as a single vector (as read from the RefMultiPedigree)
    const std::vector<double> & getAllOrigValues() const;

    ///
    /// Returns the replaced values as a single vector.
    const std::vector<double> & getAllReplacedValues() const;

    ///
    /// Returns the adjusted (post-finalization) values as a single vector.
    const std::vector<double> & getAllAdjValues() const;

    ///
    /// Returns \c true if the individual's trait value is available; \c false otherwise
    /// (that is, the individual's trait value is QNAN).
    /// \param i The member's mpindex
    bool isAdjValuePresent(size_t i) const;

  //@}

  /// @name Descriptive statistics
  //@{

    ///
    /// Returns the sampleinfo object for this field.
    const SampleInfo & getSampleInfo() const;

    ///
    /// Returns the mean of the sample's original values.
    /// Please note: This function is only available \b after the Field has been finalized.
    double getMean() const;

    ///
    /// Returns the variance of the sample's original values.
    /// Please note: This function is only available \b after the Field has been finalized.
    double getVariance() const;

    ///
    /// Returns the standard deviation of the sample's original values.
    /// Please note: This function is only available \b after the Field has been finalized.
    double getStdev() const;

  //@}

  /// @name Setting
  //@{

    /// \internal
    /// Sets the user_defined aspect.
    ///
    /// Note: This has been made available in the event that you want to first invoke createField()
    /// before later invoking importField(). Normally, importField() is all you need, but if you want
    /// to mix and match imported and created fields in the same field group, you'll need access to this function.
    int setUserDefined(bool user_defined);

  //@}

protected:

    ///
    /// Constructor
    Field();

    ///
    /// Constructor
    /// \param size The number of values that will be stored in this Field.
    explicit Field(size_t size);
    
    /// Initializes the data structure for importing data.  Clears out the
    /// value vectors, the sample info and the validity flags.  Options such
    /// as the names and ids are not cleared.
    void initializeForImport(size_t size);

    ///
    /// Sets the name of the trait originally used to read in data for this field.
    int setOrigName(const std::string & orig_name);

    ///
    /// Returns the name of the field.
    int setFieldName(const std::string & field_name);
    
    ///
    /// Returns the name of the group to which the field belongs.
    int setGroupName(const std::string & group_name);

    ///
    /// Sets this Field's FieldID (0 = success, 1 = failure)
    int setFieldID(const FieldID & id);

    ///
    /// Sets this Field's GroupID (0 = success, 1 = failure)
    int setGroupID(const GroupID & id);

    ///
    /// Sets the sampleinfo object for this field.
    int setSampleInfo(const SampleInfo & sinfo);

    ///
    /// Sets the import flags for controlling data importation.
    int setInputFlags(unsigned long flags);
 
    ///
    /// Copies the original values into the adjusted values, adjusting as necessary.
    /// Adjustment options are specified via the input flags.
    int finalizeData();
    
    ///
    /// Replaces all missing values (QNANs and Infinities) with the mean of the sample.
    /// \retval 0 Replaced values successfully.
    /// \retval 1 Did \b not replace values successfully.
    int replaceMissing();

    // Mean adjusts all values.
    // \retval 0 Values adjusted successfully.
    // \retval 1 Values \b not adjusted successfully.
    int meanAdjustValues();

    // Std. dev. adjusts all values.
    // \retval 0 Values adjusted successfully.
    // \retval 1 Values \b not adjusted successfully.
    int stdevAdjustValues();

    ///
    /// Non-const version of the getAllAdjValues() function.
    std::vector<double> & getAllAdjValues();

    ///
    /// Sets the i'th (adjusted) value to val.
    /// \param i The index number of the element to be changed (the member's mpindex())
    /// \param val The value to which the element will be changed.
    /// \retval 0 If successful.
    /// \retval 1 If \b not successful.
    int setAdjValue(size_t i, double val);

private:

  void copy(const Field &);
  int  updateStats();

  std::string         my_orig_name;
  std::string         my_field_name;
  std::string         my_group_name;
  FieldID             my_field_id;
  GroupID             my_group_id;
  unsigned long       my_input_flags;
  std::vector<double> my_orig_values;
  std::vector<double> my_replaced_values;
  std::vector<double> my_adj_values;
  std::vector<bool>   my_adj_value_presents;
  SampleInfo          my_sample_info;
  bool                my_valid;
  bool                my_finalized;
  bool                my_user_defined;
};

typedef std::vector<Field>   FieldVectorType;
typedef std::vector<FieldID> FieldGroupType;

/** \class FieldIterator
 *  \brief Non-const iterator across RPEDNEW::Field's
 *
 */
class FieldIterator
{
  public:
        // Friends
        friend class MemberDataSample;
        friend class FieldConstIterator;
        
        // Typedefs
        typedef FieldIterator                     iterator;
        typedef FieldConstIterator                const_iterator;
        typedef Field                           & reference;
        typedef Field                           * pointer;
        typedef const Field                     * const_pointer;
        typedef const Field                     & const_reference;
        typedef Field                             value_type;
        typedef FieldGroupType::difference_type   difference_type;
        typedef FieldGroupType::size_type         size_type;
        typedef std::random_access_iterator_tag        iterator_category;

        // Constructor
        FieldIterator();
        
        // Copy constructor
        FieldIterator(const FieldIterator &);
        
        // Dereference
        reference operator*  ();
        reference operator[] (difference_type i);
        
        pointer operator-> ();
        
        // Increment/Decrement
        
        iterator operator++ ();
        iterator operator++ (int);
        iterator operator-- ();
        iterator operator-- (int);
        iterator operator+= (difference_type i);
        iterator operator-= (difference_type i);
        iterator operator+  (difference_type i) const;
        iterator operator-  (difference_type i) const;
        
        difference_type operator-(const iterator & x) const;
        
        // Comparison
        
        bool operator== (const iterator & x) const;
        bool operator!= (const iterator & x) const;
        bool operator<  (const iterator & x) const;
        
  protected:
        
        // Constructor
        FieldIterator(FieldVectorType *, FieldGroupType::iterator);
        
  private:
        
        // Reference to the master list of Fields
        FieldVectorType * masterlist;
        
        // Iterator to the current element of this group:   
        FieldGroupType::iterator iter;
        
        void increment();
        void decrement();
};

/** \class FieldConstIterator
 *  \brief Const iterator across Field's
 *
 */
class FieldConstIterator
{
  public:
	// Friends
	friend class MemberDataSample;
  
	// Typedefs
	typedef FieldIterator                     iterator;
	typedef FieldConstIterator                const_iterator;
	typedef Field                           & reference;
	typedef Field                           * pointer;
	typedef const Field                     * const_pointer;
	typedef const Field                     & const_reference;
	typedef Field                             value_type;
	typedef FieldGroupType::difference_type   difference_type;
	typedef FieldGroupType::size_type         size_type;
	typedef std::random_access_iterator_tag        iterator_category;

	// Constructor
	FieldConstIterator();

	// Copy constructor
	FieldConstIterator(const FieldConstIterator &);
	FieldConstIterator(const FieldIterator       &);

	// Dereference
	const_reference operator*  ()                  const;
	const_reference operator[] (difference_type i) const;

	const_pointer operator-> ();

	// Increment/Decrement

	const_iterator operator++ ();
	const_iterator operator++ (int);
	const_iterator operator-- ();
	const_iterator operator-- (int);
	const_iterator operator+= (difference_type i);
	const_iterator operator-= (difference_type i);
	const_iterator operator+  (difference_type i) const;
	const_iterator operator-  (difference_type i) const;

	difference_type operator-(const const_iterator & x) const;

	// Comparison

	bool operator== (const const_iterator & x) const;
	bool operator!= (const const_iterator & x) const;
	bool operator<  (const const_iterator & x) const;

  protected:

	// Constructor
	FieldConstIterator(const FieldVectorType *, FieldGroupType::const_iterator);

  private:

	void increment();
	void decrement();

	// Reference to the master list of fields
	const FieldVectorType * masterlist;

	// Iterator to the current element of this group:
	FieldGroupType::const_iterator iter;
};

} // End namespace SAMPLING
} // End namespace SAGE

#include "sampling/Field.ipp"
#include "sampling/FieldItr.ipp"

#endif
