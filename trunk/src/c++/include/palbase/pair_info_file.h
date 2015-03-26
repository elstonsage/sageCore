#ifndef PAIR_INFO_FILE_H
#define PAIR_INFO_FILE_H

//****************************************************************************
//* File:      pair_info_file.h                                              *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0 Initial implementation                yjs Jan. 02 *
//*                                                                          *
//* Notes:     This file defines class for I/O of pair information file.     *
//*                                                                          *
//* Copyright (c) 2002 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "palbase/relative_pairs.h"

namespace SAGE    {
namespace PALBASE {

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
// ~ Class:   RefPairInfoFile                                                ~ 
// ~                                                                         ~ 
// ~ Purpose: This class implements streaming of pairs to and from an ASCII  ~
// ~          data file format. This file stores pairs with weight/covariate ~
// ~          values.                                                        ~
// ~                                                                         ~ 
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

class RefPairInfoFile
{
  public:

    RefPairInfoFile(cerrorstream &err = sage_cerr);

    virtual ~RefPairInfoFile();

    virtual bool  input(relative_pairs &p, const string &filename, ostream &messages = cout) = 0;

    cerrorstream error_sink() const;
    void set_error_sink(cerrorstream &err);

    size_t verbose_output() const;
    
    void set_verbose_output(size_t v = (size_t)-1);

    void validate();
    void invalidate();
    bool valid() const;

    // Provide feedback about what is read.  Some readers/writters may
    // extend this mechanism to specify data formats
    enum field_t { skip,  pedigree_id, pair_id, pair_covariate, pair_weight };
    
    struct field
    {
      field(field_t t=skip, const string &f="", const string &n="", 
            size_t i = (size_t)-1) : type(t), field_name(f), name(n), index(i) { }

      ~field() {}

      field_t type;
      string  field_name;
      string  name;
      size_t  index;
    };

    typedef std::list<field> field_list_type;
    typedef std::map<std::string, field> field_map_type;

    size_t field_count()             const;
    size_t skip_count()              const;
    size_t pedigree_id_count()       const;
    size_t pair_id_count()           const;
    size_t pair_covariate_count()    const;
    size_t pair_weight_count()       const;
    size_t invalid_covariate_count() const;
    size_t invalid_weight_count()    const;

    const field_list_type &field_list() const;

  private:

    size_t my_verbose_output;
    bool   my_valid;
      
  protected:

    void reset_counts();

    void print_pair_info_header(ostream &messages, const relative_pairs &rpairs, 
                                const string &filename) const;
    void print_pair_info_footer(ostream &messages) const;
    void print_pair_info(ostream &messages,
                         const relative_pairs &rpairs,
                         const string& pn,
                         const string& in1,
                         const string& in2,
                         const vector<pair<size_t,string> >& pair_weight_values,
                         const vector<pair<size_t,string> >& pair_covariate_values) const;

    bool validate_fields(bool quiet = false);

    void build_pair_info(relative_pairs &p);

    size_t pair_find(relative_pairs &p, const string &pn,
                     const string &id1, const string &id2,
                     size_t line);

    field_list_type &field_list();

    field_list_type my_fields;

    size_t     my_skip_count;
    size_t     my_pedigree_id_count;
    size_t     my_pair_id_count;
    size_t     my_pair_covariate_count;
    size_t     my_pair_weight_count;
    size_t     my_invalid_covariate_count;
    size_t     my_invalid_weight_count;

    cerrorstream errors;
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:   RefDelimitedPairInfoFile                                       ~
// ~                                                                         ~
// ~ Purpose: This class implements streaming of pairs to and from an ASCII  ~
// ~          data file format. This file stores pairs with weight/covariate ~
// ~          values.                                                        ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class RefDelimitedPairInfoFile : public RefPairInfoFile
{
  public:
    
    RefDelimitedPairInfoFile(cerrorstream &err = sage_cerr);

    virtual ~RefDelimitedPairInfoFile();

    virtual bool  input(relative_pairs &p, const string &filename, ostream &messages = cout);

    const string &format()     const;
    const string &delimiters() const;
    const string &whitespace() const;

    bool skip_consecutive_delimiters() const;
    bool skip_leading_delimiters()     const;
    bool skip_trailing_delimiters()    const;

    void set_format(const std::string &f);
    void set_delimiters(const std::string &d);
    void set_whitespace(const std::string &w);

    void set_skip_consecutive_delimiters(bool skip=true);
    void set_skip_leading_delimiters(bool skip=true);
    void set_skip_trailing_delimiters(bool skip=true);

    void add_pedigree_id_field();

    void add_pair_id_field();

    void add_pair_covariate_field(const std::string &field_name,
                                  const std::string &pcov_name);

    void add_pair_weight_field(const std::string &field_name,
                               const std::string &weight_name);

    const field_map_type  &field_map() const;
          field_map_type  &field_map();

  private:

    bool build_fields(string_tokenizer &header, relative_pairs &rpairs,
                      bool quiet = false);

    string my_format;
    string my_delimiters;
    string my_whitespace;
    
    bool my_skip_leading_delimiters;
    bool my_skip_trailing_delimiters;
    bool my_skip_consecutive_delimiters;

    field_map_type  my_field_map;
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:   RefLSFDelimitedPairInfoFile                                    ~
// ~                                                                         ~
// ~ Purpose: This class implements streaming of pairs to and from an ASCII  ~
// ~          data file format. This file stores pairs with weight/covariate ~
// ~          values.                                                        ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class RefLSFDelimitedPairInfoFile : public RefDelimitedPairInfoFile
{
  public:
    RefLSFDelimitedPairInfoFile(const LSFBase *params, 
                                relative_pairs &rpairs,
                                cerrorstream &errors = sage_cerr);
};

#include "pair_info_file.ipp"

void read_pair_info_file(const vector<LSF_ptr<LSFBase> > f_list, relative_pairs &p,
                         cerrorstream &errors, ostream& info_file);

} // end of namespace PALBASE
} // end of namespace SAGE

#endif
