#ifndef MFILE_H
#define MFILE_H

//
//  INHERITANCE MODEL FILE OBJECTS 0.1 -- Input and output of Inheritance Models.
//
//  Author: Geoff Wedig (wedig@darwin.cwru.edu)
//
//  History:   0.1   gcw Initial Implementation      96 07 02
//             1.01  gcw Porting to new marker locus 00 12 07
//
//  Copyright (c) 2000  R.C. Elston
//

#include <string>
#include <iostream>

#include "globals/config.h"
#include "LSF/LSFfile.h"
#include "error/errorstream.h"
#include "mlocus/imodel.h"
#include "rped/rped.h"

namespace SAGE   {
namespace MLOCUS {

class InheritanceModelFile
{
public:
  InheritanceModelFile (const cerrorstream& err= sage_cerr); 

  virtual ~InheritanceModelFile();
  
  // Note at present that we only deal with unphased genotypes.
  virtual bool  input(inheritance_model_map& m, const string &filename);
  virtual bool  input(inheritance_model_map& m, const string &filename, char separator);
  virtual bool  input(inheritance_model_map& m, const string &filename, char separator, ostream &messages);
  virtual bool  input(inheritance_model_map& m, const string &filename, const RPED::PhenotypeReaderInfo&);
  virtual bool  input(inheritance_model_map& m, const string &filename, const RPED::PhenotypeReaderInfo&, ostream &messages);

  virtual bool output(const inheritance_model_map& m, ostream &o);
  virtual bool output(const inheritance_model_map& m, const string &filename);
  virtual bool output(const inheritance_model_map& m, const string &filename, ostream &messages);
  
  cerrorstream error_sink() const       { return errors; }
  void set_error_sink(const cerrorstream& err) { errors = err;  }

  size_t marker_verbose_output()   const { return my_marker_verbose_output; }
  size_t genotype_verbose_output() const { return my_genotype_verbose_output; }

  void   set_marker_verbose_output  (size_t v = (size_t)-1) { my_marker_verbose_output = v;   }
  void   set_genotype_verbose_output(size_t v = (size_t)-1) { my_genotype_verbose_output = v; }

  void validate()    { my_valid = true;  }
  void invalidate()  { my_valid = false; }
  bool valid() const { return my_valid;  }

private:

  uint my_marker_verbose_output;
  uint my_genotype_verbose_output;

  bool my_valid;

protected:

  bool read_distances;

  void get_alleles   (istream&, genotype_model&, const RPED::PhenotypeReaderInfo&);
  void get_phenotypes(istream&, penetrance_model&);

  bool test_codominant(istream&);
  
  void test_eof(const istream&);

  void output_markers(ostream& o, const inheritance_model_map& m);
  void output_marker (ostream& o, const penetrance_model& pm, uint verbosity) const;

  cerrorstream errors;
};

/// class for turning a LSF formed set of models into inheritance_models.
class LSFInheritanceModelParser
{
  public:

    LSFInheritanceModelParser(char sep, const cerrorstream& = sage_cerr);
  
    bool input_to(inheritance_model_map&, LSFBase*);

  protected:

    void parse_model(inheritance_model_map& m, LSFBase* model);

    bool get_alleles    (LSFBase* model, genotype_model& gm);
    bool get_phenotypes (LSFBase* model, penetrance_model& gm);

    char         separator;

    cerrorstream errors;
};

} // End namespace MLOCUS
} // End namespace SAGE

#endif
