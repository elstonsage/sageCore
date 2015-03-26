#ifndef PHENOTYPE_H
#define PHENOTYPE_H

//============================================================================
//  File:       phenotype.h
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
#include "mlocus/genotype.h"

namespace SAGE   {
namespace MLOCUS {

//----------------------------------------------------------------------------
//  Class:      phenotype_info
//
//  Purpose:    This class...
//----------------------------------------------------------------------------
//
class phenotype_info
{
  public:
    phenotype_info();
    phenotype_info(const string& name, uint id);

    bool    operator ==(const phenotype_info& p) const;
    bool    operator !=(const phenotype_info& p) const;

    string  name;
    uint    id;
};


//----------------------------------------------------------------------------
//  Class:      phenotype
//
//  Purpose:    This class...
//----------------------------------------------------------------------------
//
class phenotype
{
  public:
    friend class phenotype_model;

    phenotype();
    phenotype(const phenotype& p);
    phenotype(const string& name, uint index=NPOS);
    ~phenotype();

    const phenotype&    operator =(const phenotype& p);

    bool    operator ==(const phenotype& rhs) const;
    bool    operator !=(const phenotype& rhs) const;

    const string&   name() const;
    uint            id() const;

  private:
    boost::shared_ptr<phenotype_info>   my_info;

    //lint -e{1704}
    phenotype(phenotype_info* info);

    void    set_id(uint id);
    void    uniquify();
};

} // End namespace MLOCUS
} // End namespace SAGE

#include "mlocus/phenotype.ipp"

#endif
