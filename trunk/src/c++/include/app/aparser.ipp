namespace SAGE {
namespace APP  {

// These extract the default attributes only
inline
LSFConvert::error_t
    LSFConvert::to_boolean(const LSFBase* param, bool&   value)
{
  return to_boolean(param, 0, value);
}

inline
LSFConvert::error_t
    LSFConvert::to_integer(const LSFBase* param, int&    value)
{
  return to_integer(param, 0, value);
}

inline
LSFConvert::error_t
    LSFConvert::to_string(const LSFBase* param, string& value)
{
  return to_string(param, 0, value);
}

inline
LSFConvert::error_t
    LSFConvert::to_real(const LSFBase* param, double& value)
{
  return to_real(param, 0, value);
}


// These extract any attribute by name
inline
LSFConvert::error_t
    LSFConvert::to_boolean(const LSFBase* param, const string& attr, bool&   value)
{
  return to_boolean(param, AttrNameMgr.query(attr), value);
}

inline
LSFConvert::error_t
    LSFConvert::to_integer(const LSFBase* param, const string& attr, int&    value)
{
  return to_integer(param, AttrNameMgr.query(attr), value);
}

inline
LSFConvert::error_t
    LSFConvert::to_string(const LSFBase* param, const string& attr, string& value)
{
  return to_string(param, AttrNameMgr.query(attr), value);
}

inline
LSFConvert::error_t
    LSFConvert::to_real(const LSFBase* param, const string& attr, double& value)
{
  return to_real(param, AttrNameMgr.query(attr), value);
}
// These extract the default attributes only
inline
LSFConvert::error_t
  LSFConvert::to_boolean(const LSFBase* param, bool&   value,
    cerrorstream& errors)
{
  return to_boolean(param, 0, value, errors);
}

inline
LSFConvert::error_t
  LSFConvert::to_integer(const LSFBase* param, int&    value,
    cerrorstream& errors)
{
  return to_integer(param, 0, value, errors);
}

inline
LSFConvert::error_t
  LSFConvert::to_string(const LSFBase* param, string& value,
    cerrorstream& errors)
{
  return to_string(param, 0, value, errors);
}

inline
LSFConvert::error_t
  LSFConvert::to_real(const LSFBase* param, double& value,
    cerrorstream& errors)
{
  return to_real(param, 0, value, errors);
}


// These extract any attribute by name
inline
LSFConvert::error_t
  LSFConvert::to_boolean(const LSFBase* param, const string& attr, bool&   value,
    cerrorstream& errors)
{
  return to_boolean(param, AttrNameMgr.query(attr), value, errors);
}

inline
LSFConvert::error_t
  LSFConvert::to_integer(const LSFBase* param, const string& attr, int&    value,
    cerrorstream& errors)
{
  return to_integer(param, AttrNameMgr.query(attr), value, errors);
}

inline
LSFConvert::error_t
  LSFConvert::to_string(const LSFBase* param, const string& attr, string& value,
    cerrorstream& errors)
{
  return to_string(param, AttrNameMgr.query(attr), value, errors);
}

inline
LSFConvert::error_t
  LSFConvert::to_real(const LSFBase* param, const string& attr, double& value,
    cerrorstream& errors)
{
  return to_real(param, AttrNameMgr.query(attr), value, errors);
}

// These extract the default attributes only
inline
BasicParser::error_t
BasicParser::parse_boolean(const LSFBase* param, bool&   value, bool verbose)
{
  return parse_boolean(param, 0, value, verbose);
}

inline
BasicParser::error_t
BasicParser::parse_integer(const LSFBase* param, int&    value, bool verbose)
{
  return parse_integer(param, 0, value, verbose);
}

inline
BasicParser::error_t
BasicParser::parse_string(const LSFBase* param, string& value, bool verbose)
{
  return parse_string(param, 0, value, verbose);
}

inline
BasicParser::error_t
BasicParser::parse_real(const LSFBase* param, double& value, bool verbose)
{
  return parse_real(param, 0, value, verbose);
}


// These extract any attribute by name
inline
BasicParser::error_t
BasicParser::parse_boolean(const LSFBase* param, const string& attr, bool&   value, bool verbose)
{
  return parse_boolean(param, AttrNameMgr.query(attr), value, verbose);
}

inline
BasicParser::error_t
BasicParser::parse_integer(const LSFBase* param, const string& attr, int&    value, bool verbose)
{
  return parse_integer(param, AttrNameMgr.query(attr), value, verbose);
}

inline
BasicParser::error_t
BasicParser::parse_string(const LSFBase* param, const string& attr, string& value, bool verbose)
{
  return parse_string(param, AttrNameMgr.query(attr), value, verbose);
}

inline
BasicParser::error_t
BasicParser::parse_real(const LSFBase* param, const string& attr, double& value, bool verbose)
{
  return parse_real(param, AttrNameMgr.query(attr), value, verbose);
}


// These extract any attribute by name
inline
BasicParser::error_t
BasicParser::parse_boolean(const LSFBase* param, attr_id attr, bool&   value, bool verbose)
{
  if(verbose) return LSFConvert::to_boolean(param, attr, value, errors);
  else        return LSFConvert::to_boolean(param, attr, value);
}

inline
BasicParser::error_t
BasicParser::parse_integer(const LSFBase* param, attr_id attr, int&    value, bool verbose)
{
  if(verbose) return LSFConvert::to_integer(param, attr, value, errors);
  else        return LSFConvert::to_integer(param, attr, value);
}

inline
BasicParser::error_t
BasicParser::parse_string(const LSFBase* param, attr_id attr, string& value, bool verbose)
{
  if(verbose) return LSFConvert::to_string(param, attr, value, errors);
  else        return LSFConvert::to_string(param, attr, value);
}

inline
BasicParser::error_t
BasicParser::parse_real(const LSFBase* param, attr_id attr, double& value, bool verbose)
{
  if(verbose) return LSFConvert::to_real(param, attr, value, errors);
  else        return LSFConvert::to_real(param, attr, value);
}

/// Returns the parameter/attribute pair as a string for error messages
inline string parameter_title(const LSFBase* param, attr_id id, bool cap)
{
  string p = string("parameter '") + param->name() + "'";

  if(id == 0)
  {
    if(cap) p[0] = 'P';

    return p;
  }

  string a = string("attribute '") + AttrNameMgr.query(id) + "' of " + p;

  if(cap) a[0] = 'A';

  return a;  
}

// Not inline, but here because it's a template

template<class T>
size_t check_for_incorrect_attributes
    (LSFBase*      param,
     T             exp_sym_list,
     size_t        exp_count,
     cerrorstream& errors)
{
  if(!param) return 0;

  AttrList* a_list = param->attrs();

  if(!a_list || !a_list->size()) return 0;

  size_t count = 0;

  for(AttrList::const_iterator iter = a_list->begin(); iter != a_list->end(); ++iter)
  {
    if(iter->first == 0) continue;

    string  attr_name = AttrNameMgr.query(iter->first);

    bool found = false;

    for(size_t i = 0; !found && i < exp_count; ++i)
      if(attr_name == toUpper(exp_sym_list[i]))
        found = true;

    if(!found)
    {
      errors << priority(error) << parameter_title(param, iter->first, true)
             << " not understood.  Parameter will be ignored." << endl;

      ++count;
    }
  }

  return count;
}

} // End namespace APP
} // End namespace SAGE

