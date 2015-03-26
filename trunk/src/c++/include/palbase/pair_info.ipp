//////////////////////////////////////////////////////////////////////////
//                 Implementation of base_info (Inline)                 //
//////////////////////////////////////////////////////////////////////////

inline
base_info::base_info()
{
  init();
}

inline
base_info::base_info(string s)
{
  init();
  set_name(s);
}

inline const string&
base_info::get_name() const
{
  return my_name;
}

inline void
base_info::set_name(string name)
{
  my_name = name;
}

inline void
base_info::init()
{
  my_name = "";
}

//////////////////////////////////////////////////////////////////////////
//                 Implementation of pair_marker_info (Inline)          //
//////////////////////////////////////////////////////////////////////////
/*
inline
pair_marker_info::pair_marker_info()
{
  init();
}

inline
pair_marker_info::pair_marker_info(string s, gmodel_type t, double d)
{
  init();

  set_name(s);
  set_type(t);
  set_distance(d);
}

inline gmodel_type
pair_marker_info::get_type() const
{
  return my_type;
}

inline double
pair_marker_info::get_distance() const
{
  return my_distance;
}

inline void
pair_marker_info::set_type(gmodel_type t)
{
  my_type = t;
}

inline void
pair_marker_info::set_distance(double d)
{
  my_distance = d;
}

inline void
pair_marker_info::init()
{
  base_info::init();

  my_type = AUTOSOMAL;
  my_distance = std::numeric_limits<double>::quiet_NaN();
}
*/
//////////////////////////////////////////////////////////////////////////
//                 Implementation of pair_pheno_info (Inline)           //
//////////////////////////////////////////////////////////////////////////

inline
pair_pheno_info::pair_pheno_info()
{
  init();
}

inline
pair_pheno_info::pair_pheno_info(string s, info_use u, double v)
{
  init();

  set_name(s);
  set_usage(u);
  set_value(v);
}

inline pair_pheno_info::info_use
pair_pheno_info::get_usage() const
{
  return my_usage;
}

inline double
pair_pheno_info::get_value() const
{
  return my_value;
}

inline double
pair_pheno_info::get_numeric_missing_code() const
{
  return my_numeric_missing_code;
}

inline const string&
pair_pheno_info::get_string_missing_code() const
{
  return my_string_missing_code;
}

inline void
pair_pheno_info::set_usage(pair_pheno_info::info_use u)
{
  my_usage = u;
}

inline void
pair_pheno_info::set_value(double value)
{
  my_value = value;
}

inline void
pair_pheno_info::set_numeric_missing_code(double d)
{
  my_numeric_missing_code = d;
}

inline void
pair_pheno_info::set_string_missing_code(string s)
{
  my_string_missing_code = s;
}

inline void
pair_pheno_info::init()
{
  base_info::init();

  my_usage = unknown;
  my_value = std::numeric_limits<double>::quiet_NaN();

  my_numeric_missing_code = std::numeric_limits<double>::quiet_NaN();
  my_string_missing_code  = "";
}

