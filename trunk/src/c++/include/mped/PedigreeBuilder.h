#ifndef PEDIGREE_BUILDER_H
#define PEDIGREE_BUILDER_H

#include "mped/spbaseiter.h"

namespace SAGE {
namespace MPED {

class PedigreeBuilder
{
    friend class pedigree_base;

    /** \internal
      * \brief This struct is used to represent buffered meta-relationship information and its corresponding state.
      */
    struct name_pair
    {
        /// Describes if the pair has been used in the building of the pedigree yet.
        enum UseStatus 
        {
            unused = 0, 
            possibly_used,
            used
        };
    
        name_pair();
        name_pair(const string& f, UseStatus s=unused);
        name_pair(const string& f, const string& l, UseStatus s=unused);
    
        string     name1;
        string     name2;
        UseStatus  state;
    };
    //- These container types are used to hold meta-relationship 
    //  information in the pedigree.
    //
    typedef std::list<name_pair>        marriage_list;
    typedef std::list<name_pair>        sibship_list;
    typedef std::list<member_id>        sib_chain;
    typedef std::map<string,name_pair>  lineage_map;

    bool add_lineage(const string& child, const string& parent);
    bool add_lineage(const string& child, const string& parent1, const string& parent2);

    void add_lineage(const string& child, const name_pair& parents);
    
    bool add_marriage(const string& spouse1, const string& spouse2);
    
    bool add_sibship(const string& sib1, const string& sib2);

    // Building functions
    
    void        build_pedigree(pedigree_base& ped);
    
    void        process_marriages (pedigree_base&);
    void        process_sibships  (pedigree_base&);
    void        process_lineages  (pedigree_base&);

    member_id   build_sib_chain(pedigree_base& ped, member_id kid1,
                                member_id par1, sib_chain& sibs);
                                
    void        build_subpedigrees(pedigree_base& ped);
    void        mark_family(family_id f, subped_id s);
    void        mark_all(pedigree_base& ped,subped_id src, subped_id dst);

    void        build_indices(pedigree_base& ped);

    void        infer_sexes(pedigree_base& ped);
    void        check_parental_sex_consistency(pedigree_base& ped);
    void        check_parental_missing_sexes(pedigree_base& ped);

    void        assign_arbitrary_sexes(pedigree_base& ped);
    bool        assign_arbitrary_sex(member_id ind, member_id mate,
                                     SexCode sex,
                                     std::vector<std::pair<bool, SexCode> >& visit_vector);

    void        test_marriage_loop_consistency(pedigree_base& ped);
    bool        test_marriage_loop_consistency(member_id ind, member_id mate,
                                               SexCode sex,
                                               std::vector<std::pair<bool, SexCode> >& visit_vector);
                                               
    void        set_family_mother_father(pedigree_base& ped);

    //- Cleanup
    
    void cleanup();
    
    void flush_build_info();
                                               
    //- Buffering area for meta-relationship and error information.
    //

    marriage_list   my_marriages;
    lineage_map     my_lineages;
    sibship_list    my_sibships;

    error_list      my_errors;
    error_list      my_warnings;
};
  
}
}

#endif

