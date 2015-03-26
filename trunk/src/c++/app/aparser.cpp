#include <string>
#include <functional>
#include <fstream>
#include <math.h>
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "error/errorbuf.h"
#include "error/internal_error.h"
#include "LSF/LSFsymbol.h"
#include "app/aparser.h"

using namespace std;

namespace SAGE {
namespace APP  {


/// Performs basic error testing (param ! NULL, attributes exist) and returns an error_t
inline LSFConvert::error_t basic_parse_test(const LSFBase* param, attr_id id)
{
  if(!param)
  {
    SAGE_internal_error();

    return LSFConvert::LSF_NULL;
  }

  AttrList* a_list = param->attrs();      

  if(!a_list || a_list->size() == 0)
  {
    return LSFConvert::NO_LIST;
  }

  AttrList::const_iterator iter;          

  iter = a_list->find(id);

  if(iter == a_list->end())
  {
    return LSFConvert::NO_ATTR;
  }

  return LSFConvert::GOOD;
}

LSFConvert::error_t LSFConvert::to_boolean
    (const LSFBase* param, attr_id id, bool& value)
{
  error_t err = basic_parse_test(param, id);

  if(!err)
  {
    AttrVal a = attr_value(param, id);

    if( a.has_value() )
    {
      string n = toUpper(a.String());

      if( n == "TRUE" || n == "YES" || n == "ON")
        value = true;
      else if( n == "FALSE" || n == "NO" || n == "OFF" )
        value = false;
      else                             // Value not valid
      {
        err = ATTR_INVALID;
      }
    }
    else                               // No value
    {
      err = ATTR_EMPTY;
    }
  }

  return err;
}

LSFConvert::error_t LSFConvert::to_boolean
    (const LSFBase* param, attr_id id, bool& value, cerrorstream& errors)
{
  error_t err = to_boolean(param, id, value);

  if(err != GOOD) // Most of the time, good will work.
  {

    if(err == ATTR_EMPTY)
    {
      errors << priority(warning) << parameter_title(param, id, true) 
             << " expects a value.  Parameter will be ignored."
             << endl;
    }
    else if(err == ATTR_INVALID)
    {
      errors << priority(error) << "Unknown value for "
             << parameter_title(param, id) << endl;
    }
    else if(err == LSF_NULL)
    {
      errors << priority(debug) << "Big badness.  Null LSF passed to conversion routines!" << endl;
    }
  }

  return err;
}

LSFConvert::error_t LSFConvert::to_string
    (const LSFBase* param, attr_id id, string& value)
{
  error_t err = basic_parse_test(param, id);

  if(!err)
  {
    AttrVal a = attr_value(param, id);

    if( a.has_value() )
    {
      if( a.String().size() )
        value = a.String();
      else
      {
        err = ATTR_INVALID;
      }
    }
    else
    {
      err = ATTR_EMPTY;
    }
  }

  return err;
}
LSFConvert::error_t LSFConvert::to_string
    (const LSFBase* param, attr_id id, string& value, cerrorstream& errors)
{
  error_t err = to_string(param, id, value);


  if(err != GOOD) // Most of the time, good will work.
  {
    if(err == ATTR_EMPTY)
    {
      errors << priority(warning) << parameter_title(param, id, true) 
             << " expects a value.  Parameter will be ignored."
             << endl;
    }
    else if(err == ATTR_INVALID)
    {
      errors << priority(error) << "Bad or missing value for "
             << parameter_title(param,id) << endl;
    }
    else if(err == LSF_NULL)
    {
      errors << priority(debug) << "Big badness.  Null LSF passed to conversion routines!" << endl;
    }
  }

  return err;
}

LSFConvert::error_t LSFConvert::to_real
    (const LSFBase* param, attr_id id, double& value)
{
  error_t err = basic_parse_test(param, id);

  if(!err)
  {
    AttrVal a = attr_value(param, id);

    if( a.has_value() )
    {
      if( finite(a.Real()))
        value = a.Real();
      else
      {
        err = ATTR_INVALID;
      }
    }
    else
    {
      err = ATTR_EMPTY;
    }
  }

  return err;
}

LSFConvert::error_t LSFConvert::to_real
    (const LSFBase* param, attr_id id, double& value, cerrorstream& errors)
{
  error_t err = to_real(param, id, value);

  if(err != GOOD) // Most of the time, good will work.
  {

    if(err == ATTR_EMPTY)
    {
      errors << priority(warning) << parameter_title(param, id, true) 
             << " expects a value.  Parameter will be ignored."
             << endl;
    }
    else if(err == ATTR_INVALID)
    {
      errors << priority(error) << "Bad or missing value for "
             << parameter_title(param,id) << endl;
    }
    else if(err == LSF_NULL)
    {
      errors << priority(debug) << "Big badness.  Null LSF passed to conversion routines!" << endl;
    }
  }

  return err;
}

LSFConvert::error_t LSFConvert::to_integer
    (const LSFBase* param, attr_id id, int& value)
{
  error_t err = basic_parse_test(param, id);

  if(!err)
  {
    AttrVal a = attr_value(param, id);

    if( a.has_value() )
    {
      if( finite(a.Real()) && 0 <= a.Int())
        value = a.Int();
      else
      {
        err = ATTR_INVALID;
      }
    }
    else
    {
      err = ATTR_EMPTY;
    }
  }

  return err;
}


LSFConvert::error_t LSFConvert::to_integer
    (const LSFBase* param, attr_id id, int& value, cerrorstream& errors)
{
  error_t err = to_integer(param, id, value);

  if(err != GOOD) // Most of the time, good will work.
  {

    if(err == ATTR_EMPTY)
    {
      errors << priority(warning) << parameter_title(param, id, true) 
             << " expects a value.  Parameter will be ignored."
             << endl;
    }
    else if(err == ATTR_INVALID)
    {
      errors << priority(error) << "Bad or missing value for "
             << parameter_title(param,id) << endl;
    }
    else if(err == LSF_NULL)
    {
      errors << priority(debug) << "Big badness.  Null LSF passed to conversion routines!" << endl;
    }
  }

  return err;
}

} // End namespace APP
} // End namespace SAGE

