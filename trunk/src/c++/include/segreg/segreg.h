#ifndef SEGREG_H
#define SEGREG_H

#include <limits>
#include <iostream>
#include <string>
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "segreg/segreg_input.h"

namespace SAGE
{

namespace SEGREG
{

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     SEGREG                                                        ~
// ~                                                                         ~
// ~ Purpose:   Defines segreg application derived from SAGEapp.              ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class segreg : public APP::SAGEapp
{
  public:

    segreg(int argc=0, char **argv=NULL);
   
    // Run the application
    virtual int main();

  protected:

    bool evaluate_pedigree_validity(const segreg_data&) const;
  
};

}

} // end of namespace SAGE 

#endif
