#ifndef PEDIGREE_DATA_SET_H
#define PEDIGREE_DATA_SET_H

#include <segreg/model.h>
#include <fped/fped.h>
#include <rped/rped.h>


namespace SAGE {
namespace SEGREG {

/// \brief Describes the dataset to be analyzed by a particular analysis
///
/// The dataset's responsibility is to create and provide access to the
/// analysis-specific pedigree information to be analyzed by SEGREG.  It verifies
/// that the dataset generated is consitent with the model being specified.
///
/// When the data is inconsistent with the model's requirements, it produces an
/// empty data set and appropriate error messages.
class PedigreeDataSet
{
  public:

    typedef std::list<FPED::SubpedigreeConstPointer> SubpedigreeList;

    /// \brief Provides iteration over all subpedigrees in the multipedigree
    ///
    class SubpedigreeIterator 
      : public boost::iterator_adaptor<SubpedigreeIterator,
                                       SubpedigreeList::const_iterator,
                                       const FPED::Subpedigree&,
                                       boost::use_default>
    {
      public:
        SubpedigreeIterator();
        explicit SubpedigreeIterator(const SubpedigreeList::const_iterator& i);
          
      private:
        friend class boost::iterator_core_access;
        
        const FPED::Subpedigree& dereference() const;
    };

    typedef std::pair<SubpedigreeIterator,
                      SubpedigreeIterator>           SubpedigreeCursor;

    typedef std::list<FPED::MemberConstPointer>      MemberList;

    /// \brief Provides iteration over members in the multipedigree
    ///
    class MemberIterator 
      : public boost::iterator_adaptor<MemberIterator,
                                       MemberList::const_iterator,
                                       const FPED::Member&,
                                       boost::use_default>
    {
      public:
        MemberIterator();
        explicit MemberIterator(const MemberList::const_iterator& i);
          
      private:
        friend class boost::iterator_core_access;
      
        const FPED::Member& dereference() const;
    };

    typedef std::pair<MemberIterator,
                      MemberIterator>           MemberCursor;

    PedigreeDataSet();
    PedigreeDataSet(const RPED::MultiPedigree& mped,
                    const model&               model,
                    APP::Output_Streams&       output);

    PedigreeDataSet(const PedigreeDataSet& other);

    PedigreeDataSet& operator=(const PedigreeDataSet& other);

    ~PedigreeDataSet();

    boost::shared_ptr<const FPED::Multipedigree> get_raw_data() const;
    
    bool is_valid() const;

    bool is_empty() const;

    /// \name Subpedigree functions
    //@{
    SubpedigreeIterator get_subpedigree_begin() const;
    SubpedigreeIterator get_subpedigree_end() const;

    SubpedigreeCursor   get_subpedigrees() const;

    size_t get_subpedigree_count() const;
    //@}
    
    /// \name Member functions
    ///
    /// These operations are for all members in the dataset.
    //@{
    MemberIterator get_member_begin() const;
    MemberIterator get_member_end() const;

    MemberCursor get_members() const;

    size_t get_member_count() const;
    //@}
    
    /// \name Unconnected member functions
    ///
    /// There operations are for only the members which are unconnected.
    //@{
    MemberIterator get_unconnected_begin() const;
    MemberIterator get_unconnected_end() const;

    MemberCursor get_unconnecteds() const;

    size_t get_unconnected_count() const;
    //@}
    
    MemberList getMemberList() const;

  private:

    void create_dataset(const RPED::MultiPedigree& mped, const model& model);

    bool is_valid_with(const model& model);
    
    bool has_sex_errors(const model& model);
    
    bool has_transmission_error(const model& model);

    void clear_dataset();
    
    void build_subpedigree_list();
    void build_member_list();
    void build_unconnected_list();
    
    void report_dataset_errors       (const model&         model,
                                      APP::Output_Streams& output);
    
    void report_empty_dataset_error  (APP::Output_Streams& output);
    void report_dataset_sex_errors   (const model&         model,
                                      APP::Output_Streams& output);
    void report_transmission_error   (APP::Output_Streams& output);

    boost::shared_ptr<FPED::Multipedigree> my_data;

    SubpedigreeList my_subpedigrees;
    MemberList      my_members;
    MemberList      my_unconnecteds;
};

}
}

#include <segreg/PedigreeDataSet.ipp>

#endif
