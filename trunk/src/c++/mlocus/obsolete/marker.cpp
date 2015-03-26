//============================================================================
//  File:       marker.cpp
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
#ifdef _MSC_VER
    #include <app/SAGEconfig.h>
    #pragma hdrstop
#endif

#include <mlocus.new/marker.h>
 
namespace SAGE   {
namespace MLOCUS {

//============================================================================
//  IMPLEMENTATION: marker_info
//============================================================================
//
marker_info::marker_info()
  : name(), region(), location(-1.0), id(NPOS)
{}


marker_info::marker_info(const string& n, const string& r, double l, uint i)
  : name(n), region(r), location(l), id(i)
{}


bool
marker_info::operator ==(const marker_info& mi) const
{
    //lint -e{777}
    return name == mi.name  &&  location == mi.location  &&  
            region == mi.region  &&  id == mi.id;
}

bool
marker_info::operator !=(const marker_info& mi) const
{
    //lint -e{777}
    return name != mi.name  ||  location != mi.location  ||  
            region != mi.region  ||  id != mi.id;
}


//============================================================================
//  IMPLEMENTATION: marker
//============================================================================
//
marker::marker()
  : my_info(new marker_info())
{}


marker::marker(const marker& a)
  : my_info(a.my_info)
{}


marker::marker
(const string& n, const string& r, double l, uint i)
  : my_info(new marker_info(n, r, l, i))
{}


marker::marker(marker_info* info)
  : my_info(info)
{}


marker::~marker()
{}


marker&
marker::operator =(const marker& a)
{
    if (&a != this  &&  a.my_info != my_info)
    {
        my_info = a.my_info;
    }
    return *this;
}


//----------
//
marker
marker::clone() const
{
    marker  m1(new marker_info(*my_info));

    return m1;
}


//----------
//
void
marker::reset
(const string& n, const string& r, double l, uint i)
{
    uniquify();
    my_info->name     = n;
    my_info->region   = r;
    my_info->location = l;
    my_info->id       = i;
}


void
marker::set_name(const string& n)
{
    uniquify();
    my_info->name = n;
}


void
marker::set_region(const string& r)
{
    uniquify();
    my_info->region = r;
}


void
marker::set_location(double l)
{
    uniquify();
    my_info->location = l;
}


//lint --e{1762}
void
marker::set_id(uint i)
{
//  uniquify();
    my_info->id = i;
}

} // End namespace MLOCUS
} // End namespace SAGE

