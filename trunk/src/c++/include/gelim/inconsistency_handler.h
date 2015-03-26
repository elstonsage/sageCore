#ifndef INCONSISTENCY_HANDLER_H
#define INCONSISTENCY_HANDLER_H

#include <map>
#include <set>
#include <list>
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "mped/mp.h"
#include "fped/fped.h"

namespace SAGE
{

class inconsistency_handler
{
  public:

    friend class genotype_eliminator;

    typedef inconsistency_handler                    handler;
    typedef FPED::FilteredMultipedigree              multipedigree_type;
    typedef FPED::Pedigree                           pedigree_type;
    typedef FPED::Family                             family_type;
    typedef FPED::Member                             individual_type;
    typedef pedigree_type::pedigree_const_pointer    ped_id;
    typedef pedigree_type::family_const_pointer      fam_id;
    typedef pedigree_type::member_const_pointer      ind_id;

    enum error_type { none, nuclear_family, mendelian, xy_linked };

    typedef std::vector<error_type>           error_vector;
    typedef std::map<size_t, error_type>      error_map;
    typedef std::list<fam_id>                 fam_list;
    typedef std::pair<ind_id, error_map>      ind_error_type;
    typedef std::list<ind_error_type>         family_error_type;

  protected:

    struct ped_incon
    {
      ped_incon() { }
      
      vector<int> inconsistent;   // Number of times marker i has come up inconsistent for ped.
      vector<int> informative;    // Number of times marker i has been informative for ped.
      vector<int> checked;        // Number of times marker i has been checked on this pedigree.
                                  // Likely to indicate number of pedigree sections, but not
                                  // certain
    };

    // Families are indexed by the first sibling in the family.
    typedef std::map<fam_id, family_error_type> incon_family_map;
    typedef std::map<ped_id, ped_incon>         incon_pedigree_map;

  public:

    typedef error_vector::size_type            size_type;
    typedef incon_family_map::const_iterator   incon_family_iterator;
    typedef incon_pedigree_map::const_iterator incon_pedigree_iterator;
    typedef fam_list::const_iterator           fam_iterator;

    inconsistency_handler();

    ~inconsistency_handler();

    void clear();

    size_t incon_family_count()                           const;

    const fam_list&          get_incon_family_list()                 const;
    const family_error_type& get_incon_family(const family_type& id) const;

    bool is_sex_linked_error_exist()                     const;

    // iteration - Over families

    incon_family_iterator   family_begin()   const;
    incon_family_iterator   family_end()     const;

    incon_pedigree_iterator pedigree_begin() const;
    incon_pedigree_iterator pedigree_end()   const;

  protected:

    bool add_error(const family_type& fam, size_type marker, error_type);
    bool add_error(const family_type& fam, ind_id inconsistent_child,
                   size_type marker, error_type);

    // Counts the number of times a pedigree is used with a given marker, and how
    // many are inconsistent.
    void mark_info(ped_id ped, size_type marker, bool inconsistent, bool informative);

    void build_family(const family_type&, family_error_type&);

    incon_family_map    incon_families;
    incon_pedigree_map  incon_pedigrees;  
    fam_list            incon_families_list;

    bool sex_linked_error_exist;
};

class inconsistency_printer
{
  public:

    typedef inconsistency_handler               handler;
    typedef handler::error_vector               error_vector;
    typedef handler::ind_error_type             ind_error_type;
    typedef handler::family_error_type          family_error_type;
    typedef MLOCUS::inheritance_model_map       imodel_map;

    inconsistency_printer(SAGE::cerrorstream& c = SAGE::sage_cerr) : err(c) { }

    void print_table(const handler&) const;

  protected:

    mutable SAGE::cerrorstream err;
};

class inconsistency_table_printer
{
  public:

    typedef inconsistency_handler               handler;
    typedef handler::error_vector               error_vector;
    typedef handler::ind_error_type             ind_error_type;
    typedef handler::family_error_type          family_error_type;
    typedef unsigned long                       size_type;
    typedef FPED::Multipedigree                 multipedigree;
    typedef multipedigree::subpedigree_const_pointer sped_id;
    typedef multipedigree::family_const_pointer      fam_id;
    typedef FPED::Subpedigree                   subped_type;
    typedef FPED::Member                        member_type;
    typedef MLOCUS::inheritance_model_map       imodel_map;

    inconsistency_table_printer(ostream& c = cerr, size_type col = 3)
        : columns(col), err(c) { }

    void print_pedigree_table(const handler&, bool errors_only = true) const;
    void print_marker_table(const handler&, const imodel_map&, bool errors_only = true) const;

    void print_inconsistency_table(const handler&, const imodel_map&,
                                   const vector<size_t>& s_id,
                                   bool consistent_out = false, bool sort = true) const;

    void print_inconsistency_table(const handler&, const imodel_map&,
                                   size_type start_col, size_type end_col,
                                   const vector<size_t>& s_id,
                                   bool consistent_out = false, bool sort = true) const;

  // Later
  //  void print_inconsistency_table(const handler&, vector<string> markers) const;

    // Column modification
    
    size_type set_columns(size_type s = 3) { return columns = s; }
    size_type get_columns() const { return columns; }

  protected:

    struct ped_incon
    {
      ped_incon() 
        : pedigree(), inconsistent(0), informative(0), total(0) { }
      
      ped_incon(const string& ped, unsigned int c, unsigned int i, unsigned int t)
        : pedigree(ped), inconsistent(c), informative(i), total(t) { }
    
      string        pedigree;
      unsigned int  inconsistent;
      unsigned int  informative;
      unsigned int  total;
    };

    struct ped_incon_compare
    {
      // Sort first by number of inconsistencies, then by number informative,
      // then by total number, then by pedigree name
      bool operator()(const ped_incon& l, const ped_incon& r) const
      { 
        if(l.inconsistent  != r.inconsistent)
          return (l.inconsistent  > r.inconsistent);
          
        if(l.informative != r.informative)
          return (l.informative > r.informative);

        if(l.total != r.total)
          return (l.total > r.total);

        return (l.pedigree < r.pedigree);
      }
    };

    struct marker_incon
    {
      marker_incon() : marker(0), inconsistent(0), informative(0), total(0) { }
      
      explicit marker_incon(unsigned int m)
        : marker(m), inconsistent(0), informative(0), total(0) { }

      marker_incon(unsigned int m, unsigned int c, unsigned int i, unsigned int t)
        : marker(m), inconsistent(c), informative(i), total(t) { }
    
      unsigned int  marker;
      unsigned int  inconsistent;
      unsigned int  informative;
      unsigned int  total;
    };

    struct marker_incon_compare
    {
      bool operator()(const marker_incon& l, const marker_incon& r) const
      { 
        if(l.inconsistent  != r.inconsistent)
          return (l.inconsistent  > r.inconsistent);
        
        if(l.informative != r.informative)
          return (l.informative > r.informative);

        return l.marker < r.marker;
      }
    };

    void make_marker_list(const handler&, const imodel_map&, bool) const;

    void print_header(const handler&, bool consistent_out = false) const;

    void print_table_section(const handler&, const imodel_map&,
                             size_type, size_type, const vector<size_t>&,
                             bool consistent_out = false) const;

    void print_table_header(const handler&, const imodel_map&,
                            size_type, size_type, const vector<size_t>&) const;

    bool print_inconsistency_member(stringstream& output,
                                    const member_type& mem,
                                    const vector<marker_incon>& incon_markers,
                                    const imodel_map& imap,
                                    size_type start_col, size_type end_col, size_t type,
                                    const vector<size_t>& sample_id) const;

    size_type                          columns;

    //lint -e(1725)
    ostream&                           err;

    mutable vector<marker_incon>       marker_results;
};

} // End of SAGE namespace

#include "gelim/inconsistency_handler.ipp"

#endif
