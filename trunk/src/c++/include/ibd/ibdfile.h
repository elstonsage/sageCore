#ifndef IBDFILE_H
#define IBDFILE_H

//
//  IBD sharing file -- I/O of allele sharing IBD for General
//                      data structures
//
//  Author: Kevin Jacobs (jacobs@darwin.cwru.edu)
//
//  History:   0.1  kbj Initial implementation              Aug 25 98
//
//  Copyright (c) 1998  R.C. Elston
// -------------------------------------------------------------------------
// Introduction:
//   This class implements streaming of pairs to and from an ASCII data
//   file format.  This file stores only pairs and marker allele sharing
//   at a set of loci.
//
// Warning:  These classes are only a prototype.  The interface will change
//           and much stricter checking will be added.

#include "ibd/ibd.h"

namespace SAGE {

class RefIBDReadFile
{
  public:

    RefIBDReadFile(std::ostream &output_messages = std::cerr) 
               : messages(output_messages) 
    { 
      reset();
      my_require_pedigree   = true;
      my_skip_unused_pairs  = true;
      my_warn_invalid_pairs = true;
      my_max_errors         = 50; 
    }

    virtual ~RefIBDReadFile() { }

    virtual bool input(const std::string &fname, IBD *ibd) = 0;
    virtual bool input_ibd_state(const std::string &fname, IBD *ibd) = 0;

    size_t   max_errors()                   const { return my_max_errors;         }
    size_t   error_count()                  const { return my_error_count;        }
    bool     require_pedigree()             const { return my_require_pedigree;   }
    bool     skip_unused_pairs()            const { return my_skip_unused_pairs;  }
    bool     warn_invalid_pairs()           const { return my_warn_invalid_pairs; }

    void     set_max_errors(size_t m)             { my_max_errors  = m;           }
    void     set_unlimited_errors()               { my_max_errors  = (size_t)-1;  }
    void     reset_error_count()                  { my_error_count = 0;           } 

    void     set_require_pedigree(bool o)         { my_require_pedigree   = o;    }
    void     set_skip_unused_pairs(bool o)        { my_skip_unused_pairs  = o;    }
    void     set_warn_invalid_pairs(bool o)       { my_warn_invalid_pairs = o;    }

  protected:

    void reset(const std::string &fname)
    {
      if(fname.size())
        filename = fname;
      else
        filename = "unknown file";
      reset();
    }

    void reset()
    {
      my_line = 0;
      my_error_count = 0;
    }

    void read_error(const std::string &m)
    {
      ++my_error_count;
      
      if(my_error_count == my_max_errors)
      {
        messages << "IBD file input [" << filename << ":" << my_line << "] "
                 << " Too many errors (" << my_error_count << ").  Reporting suspended" 
                 << endl;
      }
      else if(my_error_count < my_max_errors)
        messages << "IBD file input [" << filename << ":" << my_line 
                                       << "] " << m << endl;
    }

    std::ostream &messages;
    std::string filename;
    size_t my_line;
    size_t my_max_errors;
    size_t my_error_count;

    bool my_skip_unused_pairs;
    bool my_require_pedigree;
    bool my_warn_invalid_pairs;
};

class RefIBDReadFileStdIO : public RefIBDReadFile
{
  public:

    RefIBDReadFileStdIO(std::ostream &output_messages = std::cerr) 
               : RefIBDReadFile(output_messages) { }
    virtual ~RefIBDReadFileStdIO() { }

    virtual bool input(const std::string &fname, IBD *ibd);
    virtual bool input_ibd_state(const std::string &fname, IBD *ibd);

  protected:

    typedef string_tokenizer::iterator tok_iterator;

    bool do_input(FILE *file, IBD *ibd);
    bool do_input_ibd_state(FILE *file, IBD *ibd);

    void input_next_line(FILE *file, std::string &line);

    bool read_header(FILE *file, IBD *ibd, std::string &line, ibd_option_type& ibd_option, vector<size_t>& markers);
    bool read_ibd_option(FILE *file, std::string &line, ibd_option_type& io);

    bool read_pair_ids(tok_iterator &begin, tok_iterator &end, 
                       std::string &ped_name,
                       std::string &ind1_name, std::string &ind2_name);

    bool read_marker_pair_ids(tok_iterator &begin, tok_iterator &end, 
                              std::string &marker_name,
                              std::string &ped_name,
                              vector< pair<string, string> > &ind_names);
};

//
// -------------------------------------------------------------------------
//

class RefIBDWriteFile
{
  public:

    RefIBDWriteFile(const std::string &fname,
                    std::ostream &output_messages = std::cerr);

    operator void*() const { return  file; }
    int  operator!() const { return !file; }
        
    virtual ~RefIBDWriteFile();

    virtual bool output_probability_header(const IBD *ibd);
    virtual bool output_state_header(const IBD *ibd);

    virtual bool output_ibd_probability(const IBD *ibd);
    virtual bool output_ibd_state(const IBD *ibd);

  protected:

    virtual bool is_valid_ibd(const IBD *ibd);
    virtual void output_header_info(const IBD *ibd);
    
    void reset(const std::string &fname)
    {
      if(fname.size())
        filename = fname;
      else
        filename = "unknown file";
      my_line = 0;
      my_marker_count = 0;
    }

    void write_error(const std::string &m)
    {
      messages << "IBD file output [" << filename << ":" << my_line 
                                     << "] " << m << endl;
    }
    
    std::ostream& messages;
    std::ofstream file;
    std::string   filename;
    size_t        my_line;
    size_t        my_marker_count;
};

}

#endif
