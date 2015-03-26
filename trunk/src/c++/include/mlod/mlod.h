#ifndef __MLOD_H
#define __MLOD_H

#include <iostream>
#include <map>
#include "app/SAGEapp.h"
#include "error/errorstream.h"
#include "rped/genome_description.h"
#include "lvec/lvector.h"
#include "lvec/mpoint_like.h"
#include "mlod/data.h"
#include "mlod/PedigreeAnalysisSample.h"

namespace SAGE
{
namespace MLOD
{

class MLOD : public SAGE::APP::SAGEapp
{
public:

  MLOD(int argc=0, char **argv=NULL)
    : SAGE::APP::SAGEapp(SAGE::APP::APP_MLOD, true, argc, argv), 
      my_data("mlod",false)
      { }

  ~MLOD();

  virtual int main();
  
private:

  typedef Likelihood_Vector                     lvector;
  typedef lvector::size_type                    size_type;
  typedef mpoint_likelihood_data                ldata;
  typedef RPED::genome_description::region_type region_type;

  typedef std::vector<boost::shared_ptr<PedigreeAnalysisSample> > PASVectorType;

  void init();

  void run_analyses();

  Data my_data;
};

// ================
// Inline functions
// ================

inline MLOD::~MLOD()
{ }

}
}

#endif
