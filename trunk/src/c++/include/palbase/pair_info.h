#ifndef PAIR_INFO_H
#define PAIR_INFO_H

//****************************************************************************
//* File:      pair_info.h                                                   *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Initial implementation                                        *
//*                                                                          *
//* Notes:     This header file defines the following classes.               *
//*              pair_marker_info                                            *
//*              pair_pheno_info                                             *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "palbase/definitions.h"

namespace SAGE    {
namespace PALBASE {

class base_info
{
  public:

    base_info();
    base_info(string s);
    
    const string& get_name()    const;

    void set_name(string name);

  protected:

    void init();

    string        my_name;
};
/*
class pair_marker_info : public base_info
{
  public:

    pair_marker_info();

    explicit pair_marker_info(string      s,
                              gmodel_type t,
                              double      d = std::numeric_limits<double>::quiet_NaN());

    gmodel_type   get_type()    const;
    double        get_distance()const;

    void set_type(gmodel_type t);
    void set_distance(double d);

  private:

    void init();

    gmodel_type   my_type;
    double        my_distance;
};
*/
class pair_pheno_info : public base_info
{
  public :

    enum info_use { unknown, mean, minimum };
    
    pair_pheno_info();
    
    explicit pair_pheno_info(string   name,
                             info_use u = unknown,
                             double   value = std::numeric_limits<double>::quiet_NaN());

    info_use      get_usage() const;
    double        get_value() const;
    
    double        get_numeric_missing_code() const;
    const string& get_string_missing_code()  const;
    
    void set_usage(info_use u);
    void set_value(double value);
    void set_numeric_missing_code(double d);
    void set_string_missing_code(string s);
  
  private:
  
    void init();
    
    info_use my_usage;
    double   my_value;
    double   my_numeric_missing_code;
    string   my_string_missing_code;
};

#include "pair_info.ipp"

} // end of namespace PALBASE
} // end of namespace SAGE

#endif
