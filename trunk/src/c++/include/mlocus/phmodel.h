#ifndef PHMODEL_H
#define PHMODEL_H

//============================================================================
//  File:       phmodel.h
//
//  Author:     Bob Steagall
//
//  History:    Version 1.00.00
//
//  Notes:    
//
//  Copyright (c) 1999 R.C. Elston
//  All Rights Reserved
//============================================================================
//
#include "mlocus/phenotype.h"

#ifndef SAGE_ASSERT
    #define SAGE_ASSERT(A,B)
#endif

namespace SAGE   {
namespace MLOCUS {

//----------------------------------------------------------------------------
//  Class:      phenotype_model_info
//
//  Purpose:    This class...
//----------------------------------------------------------------------------
//
class phenotype_model_info
{
  public:
    friend class phenotype_model;

    typedef std::vector<phenotype>          phenotype_vector;
    typedef std::map<string, uint>          phenotype_map;

    phenotype_model_info();
    phenotype_model_info(const string& name);

    string              name;
    phenotype_vector    phenotypes;
    phenotype_map       phenotype_names;
    string              missing_ptname;
    uint                missing_ptid;
    
    bool                my_has_generated_phenotypes; /// Set to \c true if there are 
                                                     /// phenotypes which we generated
                                                     /// from genotype names, \c false otherwise
    bool                my_has_external_phenotypes;  /// Set to \c true if there are
                                                     /// phenotypes which were created
                                                     /// from external sources (add_phenotype()),
                                                     /// \c false otherwise.

};


//----------------------------------------------------------------------------
//  Class:      phenotype_model
//
//  Purpose:    This class...
//----------------------------------------------------------------------------
//
class phenotype_model
  {
  public:

    friend class penetrance_model;
    friend class phenotype;

    typedef std::vector<phenotype>              phenotype_vector;
    typedef std::map<string, uint>              phenotype_map;
    typedef phenotype_vector::const_iterator    phenotype_iterator;
    typedef phenotype_vector::const_iterator    phenotype_const_iterator;

    phenotype_model();
    phenotype_model(const genotype_model& gm);
    phenotype_model(const genotype_model& gm, const string& name);
    phenotype_model(const phenotype_model& m);
    ~phenotype_model();

    phenotype_model&   operator =(const phenotype_model& m);

    //- Counts of various contained elements.
    //
    uint        phenotype_count() const;
    uint        alias_count() const;     // Total number of acceptable aliases

    //- Iteration over phenotypes.
    //
    phenotype_iterator      phenotype_begin() const;
    phenotype_iterator      phenotype_end() const;

    //- Element indexing by ID.
    //
    const phenotype&    get_phenotype(uint id) const;

    //- Direct element lookup by name.
    //
    phenotype           get_phenotype(const string& name) const;
    const phenotype&    get_missing_phenotype() const;

    //- Element ID lookup by name.
    //
    uint        get_phenotype_id(const string& name) const;
    uint        get_missing_phenotype_id() const;

    //- Other properties.
    //
    phenotype_model         clone() const;
    const string&           name() const;
    const string&           missing_phenotype_name() const;

    bool        has_generated_phenotypes() const;
    bool        has_external_phenotypes() const;

    //- Modifiers.
    //
    void        add_phenotype(const string& pname);
    void        alias_phenotype(const string& pname, const string& alias);
    void        alias_phenotype(uint pid, const string& alias);

    void        add_genotypes(const genotype_model&);
    void        add_allele(const string& name, const genotype_model&);

    void        clear();
    void        set_name(const string& name);
    void        set_missing_phenotype_name(const string& name);

  private:
    boost::shared_ptr<phenotype_model_info>    my_info;
    
    void        add_autosomal_phenotypes_by_allele(const string& name, const genotype_model&);
    void        add_female_x_linked_phenotypes_by_allele(const string& name, const genotype_model&);
    void        add_male_x_linked_phenotypes_by_allele(const string& name, const genotype_model&);
    void        add_female_y_linked_phenotypes_by_allele(const string& name, const genotype_model&);
    void        add_male_y_linked_phenotypes_by_allele(const string& name, const genotype_model&);

    //lint -e{1704}
    phenotype_model(phenotype_model_info* info);

    void        push_phenotype(const string& name, const genotype_model&);
    void        build_phenotypes_symmetric(const genotype_model& gm);

    void        add_missing_phenotype();
    void        add_missing_phenotype_aliases(const genotype_model&);

    void        uniquify();
};

} // End namespace MLOCUS
} // End namespace SAGE

#include "mlocus/phmodel.ipp"

#endif
