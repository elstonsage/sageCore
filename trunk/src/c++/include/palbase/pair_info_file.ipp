////////////////////////////////////////////////////////////////////////////
//             Implementation of pairinfo_file(Inline)                    //
////////////////////////////////////////////////////////////////////////////

//
//--------------------------------------------------------------------------
//

inline
cerrorstream
RefPairInfoFile::error_sink() const
{
  return errors;
}

inline
void
RefPairInfoFile::set_error_sink(cerrorstream &err)
{
  errors = err;
}

inline
size_t
RefPairInfoFile::verbose_output() const
{
  return my_verbose_output;
}
    
inline
void
RefPairInfoFile::set_verbose_output(size_t v)
{
  my_verbose_output = v;
}

inline
void
RefPairInfoFile::validate()
{
  my_valid = true;
}

inline
void
RefPairInfoFile::invalidate()
{
  my_valid = false;
}

inline
bool
RefPairInfoFile::valid() const
{
  return my_valid;
}

inline
size_t
RefPairInfoFile::field_count() const
{
  return my_fields.size();
}

inline
size_t
RefPairInfoFile::skip_count() const
{
  return my_skip_count;
}

inline
size_t
RefPairInfoFile::pedigree_id_count() const
{
  return my_pedigree_id_count;
}

inline
size_t
RefPairInfoFile::pair_id_count() const
{
  return my_pair_id_count;
}

inline
size_t
RefPairInfoFile::pair_covariate_count() const
{
  return my_pair_covariate_count;
}

inline
size_t
RefPairInfoFile::pair_weight_count() const
{
  return my_pair_weight_count;
}

inline
size_t
RefPairInfoFile::invalid_covariate_count() const
{
  return my_invalid_covariate_count;
}

inline
size_t
RefPairInfoFile::invalid_weight_count() const
{
  return my_invalid_weight_count;
}

inline
const RefPairInfoFile::field_list_type&
RefPairInfoFile::field_list() const
{
  return my_fields;
}

inline
RefPairInfoFile::field_list_type&
RefPairInfoFile::field_list()
{
  return my_fields;
}

//
//--------------------------------------------------------------------------
//

inline
const string&
RefDelimitedPairInfoFile::format() const
{
  return my_format;
}

inline
const string&
RefDelimitedPairInfoFile::delimiters() const
{
  return my_delimiters;
}

inline
const string&
RefDelimitedPairInfoFile::whitespace() const
{
  return my_whitespace;
}

inline
bool
RefDelimitedPairInfoFile::skip_consecutive_delimiters() const
{
  return my_skip_consecutive_delimiters;
}

inline
bool
RefDelimitedPairInfoFile::skip_leading_delimiters() const
{
  return my_skip_leading_delimiters;
}

inline
bool
RefDelimitedPairInfoFile::skip_trailing_delimiters() const
{
  return my_skip_trailing_delimiters;
}

inline
void
RefDelimitedPairInfoFile::set_format(const std::string &f)
{
  my_format = f;
}

inline
void
RefDelimitedPairInfoFile::set_delimiters(const std::string &d)
{
  my_delimiters = d;
}

inline
void
RefDelimitedPairInfoFile::set_whitespace(const std::string &w)
{
  my_whitespace = w;
}

inline
void
RefDelimitedPairInfoFile::set_skip_consecutive_delimiters(bool skip)
{
  my_skip_consecutive_delimiters = skip;
}

inline
void
RefDelimitedPairInfoFile::set_skip_leading_delimiters(bool skip)
{
  my_skip_leading_delimiters = skip;
}

inline
void
RefDelimitedPairInfoFile::set_skip_trailing_delimiters(bool skip)
{
  my_skip_trailing_delimiters = skip;
}

inline
void
RefDelimitedPairInfoFile::add_pedigree_id_field()
{
  field_list().push_back( field(pedigree_id) );
}  

inline
void
RefDelimitedPairInfoFile::add_pair_id_field()
{
  field_list().push_back( field(pair_id) );
}  

inline
void
RefDelimitedPairInfoFile::add_pair_covariate_field(const std::string &field_name,
                                                   const std::string &pcov_name)
{
  my_field_map[ toUpper(field_name) ] = field(pair_covariate, field_name, pcov_name);
}

inline
void
RefDelimitedPairInfoFile::add_pair_weight_field(const std::string &field_name,
                                                const std::string &weight_name)
{
  my_field_map[ toUpper(field_name) ] = field(pair_weight, field_name, weight_name);
}

inline
const RefDelimitedPairInfoFile::field_map_type&
RefDelimitedPairInfoFile::field_map() const
{
  return my_field_map;
}

inline
RefDelimitedPairInfoFile::field_map_type&
RefDelimitedPairInfoFile::field_map()
{
  return my_field_map;
}


