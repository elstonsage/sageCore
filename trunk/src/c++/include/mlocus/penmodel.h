#ifndef PENMODEL_H
#define PENMODEL_H

//============================================================================
//  File:       penmodel.h
//
//  Author:     Bob Steagall
//
//  History:    Version 1.00.00
//              X, Y-linkage added   - yjs  Mar. 2002
//
//  Notes:    
//
//  Copyright (c) 1999 R.C. Elston
//  All Rights Reserved
//============================================================================
//
#include "mlocus/phmodel.h"
#include "mlocus/penetrance_matrix.h"
#include "output/Output.h"
#include "error/internal_error.h"

#ifndef SAGE_ASSERT
    #define SAGE_ASSERT(A,B)
#endif

namespace SAGE   {
namespace MLOCUS {

//----------------------------------------------------------------------------
//  Class:      penetrance_info
//
//  Purpose:    This class...
//----------------------------------------------------------------------------
//
class penetrance_info
{
  public:
    penetrance_info();
    penetrance_info(double pval, const string& p, const string& a1, const string& a2, Ordering o);

    double      pvalue;
    string      pname;
    string      a1name;
    string      a2name;
    Ordering    order;
};


//----------------------------------------------------------------------------
//  Class:      penetrance_model_info
//
//  Purpose:    This class...
//----------------------------------------------------------------------------
//
class penetrance_model_info
{
  public:

    typedef std::vector<phenotype>                 phenotype_vector;
    typedef std::vector<bool>                      codominance_vector;
    typedef penetrance_matrix<int,int>             phased_penetrance_matrix;
    typedef penetrance_matrix<int,int>             unphased_penetrance_matrix;
    typedef std::vector<uint>                      count_vector;

    penetrance_model_info();
    penetrance_model_info(const penetrance_model_info&);
    penetrance_model_info(const string& name);
    penetrance_model_info(const genotype_model& gm);
    penetrance_model_info(const genotype_model& gm, const string& name);
    penetrance_model_info(const genotype_model& gm, const phenotype_model& pm);
    penetrance_model_info(const genotype_model& gm, const phenotype_model& pm,
                          const string& name);

    genotype_model      gmodel;
    phenotype_model     phmodel;
    codominance_vector  codominant_phenotypes;
    codominance_vector  strict_phenotypes;
    string              name;
    phased_penetrance_matrix   phased_penetrance;
    unphased_penetrance_matrix unphased_penetrance;
    uint                codominant;           // Number of non-strict phenotypes that
                                              // are *not* codominant.
    uint                strict_codominant;    // Number of phenotypes that are *not*
                                              // codominant

    uint                last_alleles;

  protected:

    void init();
};




//----------------------------------------------------------------------------
//  Class:      penetrance_model
//
//  Purpose:    This class...
//----------------------------------------------------------------------------
//
class penetrance_model
{
  public:
    class phased_penetrance_iterator;
    class unphased_penetrance_iterator;

    friend class phenotype;
    friend class phased_penetrance_iterator;
    friend class unphased_penetrance_iterator;

    typedef std::vector<phenotype>                      phenotype_vector;
    typedef phenotype_vector::const_iterator            phenotype_iterator;
    typedef phenotype_vector::const_iterator            phenotype_const_iterator;

    typedef penetrance_matrix<int,int>             phased_penetrance_matrix;
    typedef penetrance_matrix<int,int>             unphased_penetrance_matrix;

    class phased_penetrance_iterator : public std::forward_iterator_tag
    {
      public:
        friend class penetrance_model;
        typedef double                      reference;
        typedef penetrance_model            host_type;
        typedef const penetrance_model*     host_pointer;
        typedef phased_penetrance_iterator  this_type;

        phased_penetrance_iterator()     : my_host(0), my_row(0) {}

        bool    operator ==(const phased_penetrance_iterator& i) const  { return my_host == i.my_host  &&  my_element == i.my_element; }
        bool    operator !=(const phased_penetrance_iterator& i) const  { return my_host != i.my_host  ||  my_element != i.my_element; }

        reference   operator *() const                                  { return my_element->second;     }
        this_type&  operator ++()                                       { ++my_element;  return *this;  }
        this_type   operator ++(int)                                    { this_type tmp(*this); ++my_element; return tmp; }

        int         geno_id() const                                     { return my_element->first.col; }
        uint        phenotype_id() const                                { return my_row;                }

        phased_genotype         phased_geno() const                 { return /*lint -e{613} */ my_host->get_phased_genotype(geno_id()); }
        const phenotype&        pheno() const                       { return /*lint -e{613} */ my_host->get_phenotype(my_row); }

      private:
        host_pointer    my_host;
        uint            my_row;
        phased_penetrance_matrix::row_iterator   my_element;

        //lint -e{1704} <-- intended private constructor
        phased_penetrance_iterator(const host_type& h, uint row, bool end);
    };

    class unphased_penetrance_iterator : public std::forward_iterator_tag
    {
      public:
        friend class penetrance_model;
        typedef double                          reference;
        typedef penetrance_model                host_type;
        typedef const penetrance_model*         host_pointer;
        typedef unphased_penetrance_iterator    this_type;

        unphased_penetrance_iterator() : my_host(0), my_row(0) {}

        bool    operator ==(const unphased_penetrance_iterator& i) const{ return my_host == i.my_host  &&  my_element == i.my_element; }
        bool    operator !=(const unphased_penetrance_iterator& i) const{ return my_host != i.my_host  ||  my_element != i.my_element; }

        reference   operator *() const                                  { return my_element->second;     }
        this_type&  operator ++()                                       { ++my_element;  return *this;  }
        this_type   operator ++(int)                                    { this_type tmp(*this); ++my_element; return tmp; }

        uint        geno_id() const                                     { /*lint -e{732} */ return my_element->first.col; }
        uint        phenotype_id() const                                { return my_row;            }

        unphased_genotype           unphased_geno() const           { return /*lint -e{613} */ my_host->get_unphased_genotype(geno_id()); }
        const phenotype&            pheno() const                   { return /*lint -e{613} */ my_host->get_phenotype(my_row);           }

        /// Often, we wish to determine if the current genotype phenotype
        /// combination is only consistent with one of its phased counterparts.
        /// This function returns \c true if exactly one of the phased genotypes
        /// has a penetrance.  In the case of homozygosity, returns \c false.
        bool is_phased() const;
        
      private:
        host_pointer    my_host;
        uint            my_row;
        unphased_penetrance_matrix::row_iterator   my_element;

        //lint -e{1704} <-- intended private constructor
        unphased_penetrance_iterator(const host_type& h, uint row, bool end);
    };

    penetrance_model();
    penetrance_model(const genotype_model& gm, bool add_phenotypes = false,
                     bool add_penetrances = false);
    penetrance_model(const genotype_model& gm, const phenotype_model& pm,
                     bool add_penetrances = false);
    penetrance_model(const genotype_model& gm, const string& name,
                     bool add_phenotypes = false, bool add_penetrances = false);
    penetrance_model(const genotype_model& gm, const phenotype_model& pm,
                     const string& name, bool add_penetrances = false);
    penetrance_model(const penetrance_model& m);
    ~penetrance_model();

    penetrance_model&   operator =(const penetrance_model& m);

    //- Counts of various contained elements.
    //
    uint        allele_count() const;
    uint        phased_genotype_count() const;
    uint        unphased_genotype_count() const;
    uint        phenotype_count() const;

    //- Iteration over alleles.
    //
    allele_iterator     allele_begin() const;
    allele_iterator     allele_end() const;

    //- Iteration over genotypes.
    //
    phased_genotype_iterator    phased_genotype_begin() const;
    phased_genotype_iterator    phased_genotype_end() const;
    unphased_genotype_iterator  unphased_genotype_begin() const;
    unphased_genotype_iterator  unphased_genotype_end() const;

    //- Iteration over phenotypes.
    //
    phenotype_iterator      phenotype_begin() const;
    phenotype_iterator      phenotype_end() const;

    //- Iteration over the penetrance matrices.
    //
    phased_penetrance_iterator      phased_penetrance_begin(uint id) const;
    phased_penetrance_iterator      phased_penetrance_end(uint id) const;
    phased_penetrance_iterator      phased_penetrance_begin(const phenotype& p) const;
    phased_penetrance_iterator      phased_penetrance_end(const phenotype& p) const;

    unphased_penetrance_iterator    unphased_penetrance_begin(uint id) const;
    unphased_penetrance_iterator    unphased_penetrance_end(uint id) const;
    unphased_penetrance_iterator    unphased_penetrance_begin(const phenotype& p) const;
    unphased_penetrance_iterator    unphased_penetrance_end(const phenotype& p) const;

    uint phased_penetrance_count(uint id) const;
    uint phased_penetrance_count(const phenotype& p) const;
    uint unphased_penetrance_count(uint id) const;
    uint unphased_penetrance_count(const phenotype& p) const;

    //- Element indexing by ID.
    //
    allele                      get_allele(uint id) const;
    phased_genotype             get_phased_genotype(int id) const;
    unphased_genotype           get_unphased_genotype(uint id) const;
    const phenotype&            get_phenotype(uint id) const;

    //- Direct element lookup by name.
    //
    allele              get_allele(const string& name) const;
    phased_genotype     get_phased_genotype(const string& name) const;
    unphased_genotype   get_unphased_genotype(const string& name) const;
    phenotype           get_phenotype(const string& name) const;
    phenotype           get_missing_phenotype() const;

    //- Element ID lookup by name.
    //
    uint        get_phenotype_id(const string& name) const;
    uint        get_missing_phenotype_id() const;

    //- Other properties.
    //
    penetrance_model        clone() const;
    const genotype_model&   gmodel() const;
    genotype_model&         gmodel();
    phenotype_model&        phmodel();
    const phenotype_model&  phmodel() const;
    const string&           name() const;
    const string&           missing_allele_name() const;
    const string&           missing_phenotype_name() const;

    /** Informativity options measure the informativity of the penetrance model.
     *  There are two types of informativity: genotypic and penetrance.
     *
     *  - Genotypic informativity detects the presence of both penetrant and
     *    non-penetrant genotypes.  Ie, there is at least one genotype for
     *    which the penetrance is zero and one for which the penetrance is >
     *    0.
     *
     *  - Penetrance informativity is a looser form of informativeness.  It
     *    detects the presence of non-equal genotype penetrance.  Ie, there
     *    is at least one genotype pair at a single phenotype for which
     *    P(p | g1) != P(p | g2).
     */
    //@{
    bool        genotype_informative_phenotype(const phenotype&) const;
    bool        genotype_informative_phenotype(uint id) const;

    bool        penetrance_informative_phenotype(const phenotype&) const;
    bool        penetrance_informative_phenotype(uint id) const;

    bool        genotype_informative() const;
    bool        penetrance_informative() const;
    //@}
    
    bool        strict_phenotype(const phenotype&) const;
    bool        strict_phenotype(uint id) const;

    bool        codominant(bool strict = true) const;
    bool        codominant(uint ptid, bool strict = true) const;

    double      phased_penetrance(uint ptid, int gtid) const;
    double      phased_penetrance(const phenotype& p, const phased_genotype& g) const;
    double      unphased_penetrance(uint ptid, uint gtid) const;
    double      unphased_penetrance(const phenotype& p, const phased_genotype& g) const;

    // Added for X-linkage - yjs Mar. 2002
    //
    bool        is_x_linked() const;
    bool        is_y_linked() const;
    
    bool        is_autosomal() const;
    
    GenotypeModelType get_model_type() const;

    void        set_model_type(GenotypeModelType);

    //- Modifiers.
    //
    const string& set_name(const string&);

    void        add_allele(const string& name, double freq=0.0,
                           bool add_phenotypes = false, bool add_penetrances = false);
                           
    bool        add_genotype_dynamically(const string& genotype);

    uint        add_phenotype(const string& pname, bool strict = true);
    void        alias_phenotype(const string& pname, const string& alias);
    void        alias_phenotype(uint pid, const string& alias);

    void        set_phenotype_strict(uint pid, bool strict = true);

    void        copy_penetrance(uint psource, uint pdest, bool override = false);
    void        copy_penetrance(const penetrance_model& model, uint psource, uint pdest, bool override = false);

    void        copy_penetrance_sexed(uint psource, uint pdest, MPED::SexCode sex, bool override = false);
    void        copy_penetrance_sexed(const penetrance_model& model, uint psource, uint pdest, MPED::SexCode sex, bool override = false);

    /// Clears the entire row of penetrances
    void        clear_penetrance(uint pdest);

    void        add_penetrance(double pval, const string& pname, const string& gname);
    void        add_penetrance(double pval, const string& pname, const string& gname, Ordering order);
    void        add_penetrance(double pval, const string& pname, const string& gname, Ordering order, char sep);
    void        add_phased_penetrance(double pval, uint pid, int gid);
    void        add_unphased_penetrance(double pval, uint pid, uint gid);
    void        remove_phased_penetrance(uint pid, int gid, bool maintain_consistency = true);
    void        remove_unphased_penetrance(uint pid, uint gid, bool maintain_consistency = true);
    void        clear();
    void        set_missing_phenotype_name(const string& name);

    // Checks for consistency of the penetrance values for each phenotype,
    // genotype pair for which values exist
    void        make_consistent();

    // Remapping

    void        mark_for_remap(const string& name);
    void        remap();
    void        remap(penetrance_model& m) const;

  private:
    boost::shared_ptr<penetrance_model_info>    my_info;

    //lint -e{1704} <-- intended private constructor
    penetrance_model(penetrance_model_info* info);

    void        build_symmetric_penetrances();

    void        add_penetrance_with_missing_allele
                    (uint pid, allele a, Ordering ordering, double val);

    bool        check_for_codominance(uint);
    void        reset_codominance(uint, bool);

    void        create_missing_phenotype_penetrance(uint pid, double val = 1.0);
    
    void        resize_matrices();

    void        uniquify();

    static void remap(const penetrance_model& src, penetrance_model& dst);
};

} // End namespace MLOCUS
} // End namespace SAGE

#include "mlocus/penmodel.ipp"

#endif
