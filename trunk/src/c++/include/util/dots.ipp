// ==================
// Inlines for dots.h
// ==================

inline dot_formatter::dot_formatter()
  : my_trigger_count(0), my_current_value(0)
{ }

inline
size_t dot_formatter::set_trigger_count(size_t t)
{
  if(!done()) finish();

  if(my_current_value >= my_trigger_count)
    my_current_value = t;  

  return my_trigger_count = t;
}

inline
size_t dot_formatter::get_trigger_count() const
{
  return my_trigger_count;
}  
  
inline bool dot_formatter::done() const
{
  return my_trigger_count <= my_current_value;
}

inline
bool dot_formatter::trigger()
{
  if(my_current_value >= my_trigger_count)
  {
    if(my_trigger_count == 0) return false;

    my_current_value = 0;

    prefix();
  }
  
  ++my_current_value;
  
  resync();

  if(my_trigger_count <= my_current_value)
  {
    postfix();

    return true;
  }

  return false;
}

inline
bool dot_formatter::finish()
{
  if(my_current_value >= my_trigger_count) return true;

  my_current_value = my_trigger_count;
  
  resync();

  postfix();

  return true;
}

inline
bool dot_formatter::spammed()
{
  // If the current line is complete, being spammed is no problem
  if(my_current_value >= my_trigger_count) return true;

  prefix();

  reset();
  
  resync();

  return true;
}

inline
text_dot_formatter::text_dot_formatter(std::ostream& o)
   : my_stream(o), my_postfix("Done.")
{
  my_prefix_width = 0;
  my_postfix_width = 0;

  my_prefix_char = ' ';
  my_postfix_char = ' ';

  my_output_count = 20;
  
  my_output_char = '.';

  my_current_width = 0;
}

inline
size_t text_dot_formatter::set_output_character_count(size_t t)
{
  my_output_count = t;

  // We need to resynchronize the line to the new fill width.  If
  // the width is too large, the resync will fail.  In this case, we assume
  // the line is spammed, and reproduce the line.
  if(!done() && !resync())
    spammed();

  return my_output_count;
}

inline size_t text_dot_formatter::get_output_character_count() const
{ return my_output_count; }

inline char text_dot_formatter::set_output_char(char c)
{
  bool b = !done() && c != my_output_char;
  
  my_output_char = c;
  
  if(b)
  {
    spammed();
  }

  return my_output_char;
}

inline char text_dot_formatter::get_output_char() const
{
  return my_output_char;
}

inline
const text_dot_formatter::string_type& text_dot_formatter::set_prefix
    (const string_type& s)
{
  bool b = !done() && my_prefix != s;

  my_prefix = s;

  if(b) finish();

  return my_prefix;
}

inline
const text_dot_formatter::string_type& text_dot_formatter::get_prefix() const
{
  return my_prefix;
}

inline
size_t text_dot_formatter::set_prefix_width(size_t t, char c)
{
  my_prefix_width = t;
  my_prefix_char = c;

  return t;
}

inline size_t text_dot_formatter::get_prefix_width() const
{ return my_prefix_width; }

inline char text_dot_formatter::get_prefix_char() const
{ return my_prefix_char; }

inline
const text_dot_formatter::string_type& text_dot_formatter::set_postfix
    (const string_type& s)
{
  my_postfix = s;

  return my_postfix;
}

inline
const text_dot_formatter::string_type& text_dot_formatter::get_postfix() const
{
  return my_postfix;
}

inline
size_t text_dot_formatter::set_postfix_width(size_t t, char c)
{
  my_postfix_width = t;
  my_postfix_char = c;

  return t;
}

inline size_t text_dot_formatter::get_postfix_width() const
{ return my_postfix_width; }

inline char text_dot_formatter::get_postfix_char() const
{ return my_postfix_char; }

inline
bool text_dot_formatter::prefix()
{
  char fill = my_stream.fill();
  
  if(fill != my_prefix_char) my_stream.fill(my_prefix_char);

  if(my_prefix_width) my_stream << std::left << std::setw(my_prefix_width);
  
  my_stream << my_prefix << std::flush;

  if(fill != my_prefix_char) my_stream.fill(fill);

  return true;
}

inline
bool text_dot_formatter::postfix()
{
  char fill = my_stream.fill();
  
  if(fill != my_postfix_char) my_stream.fill(my_postfix_char);

  if(my_postfix_width) my_stream << std::setw(my_postfix_width);
  
  my_stream << my_postfix << std::endl;

  if(fill != my_postfix_char) my_stream.fill(fill);

  my_current_width = 0;

  return true;
}

inline
bool text_dot_formatter::resync()
{
  // We can't resync if this is the case.
  if(my_output_count < my_current_width) return false;

  double value = my_current_value * my_output_count / (double) my_trigger_count;

  while(my_current_width < value)
  {
    ++my_current_width;

    my_stream << my_output_char << std::flush;
  }

  return true;
}

inline
bool text_dot_formatter::reset()
{
  my_current_width = 0;

  return true;
}

inline
time_dot_formatter::time_dot_formatter(std::ostream& o)
   : my_stream(o), my_prefix(""), my_postfix("Done.")
{
  my_prefix_width = 32;
  my_postfix_width = 0;

  my_prefix_char = ' ';
  my_postfix_char = ' ';

  interval = 0;

  set_output_time_count();
  set_max_time_interval();
  set_min_time_interval();
  set_time_check_count();
}

inline
size_t time_dot_formatter::set_output_time_count(size_t t)
{
  interval = 0;

  return output_time_count = t;
}

inline size_t time_dot_formatter::get_output_time_count() const
{ return output_time_count; }

inline
size_t time_dot_formatter::set_max_time_interval(size_t t)
{
  interval = 0;

  return max_time_interval = t;
}

inline size_t time_dot_formatter::get_max_time_interval() const
{ return max_time_interval; }

inline
size_t time_dot_formatter::set_min_time_interval(size_t t)
{
  interval = 0;

  return min_time_interval = t;
}

inline size_t time_dot_formatter::get_min_time_interval() const
{ return min_time_interval; }

inline
size_t time_dot_formatter::set_time_check_count(size_t t)
{
  return time_check_count = t;
}

inline size_t time_dot_formatter::get_time_check_count() const
{ return time_check_count; }

inline
const time_dot_formatter::string_type& time_dot_formatter::set_prefix
    (const string_type& s)
{
  bool b = !done() && my_prefix != s;

  my_prefix = s;

  if(b)
  {
    finish();
    interval = 0;
  }

  return my_prefix;
}

inline
const time_dot_formatter::string_type& time_dot_formatter::get_prefix() const
{
  return my_prefix;
}

inline
size_t time_dot_formatter::set_prefix_width(size_t t, char c)
{
  my_prefix_width = t;
  my_prefix_char = c;

  return t;
}

inline size_t time_dot_formatter::get_prefix_width() const
{ return my_prefix_width; }

inline char time_dot_formatter::get_prefix_char() const
{ return my_prefix_char; }

inline
const time_dot_formatter::string_type& time_dot_formatter::set_postfix
    (const string_type& s)
{
  my_postfix = s;

  return my_postfix;
}

inline
const time_dot_formatter::string_type& time_dot_formatter::get_postfix() const
{
  return my_postfix;
}

inline
size_t time_dot_formatter::set_postfix_width(size_t t, char c)
{
  my_postfix_width = t;
  my_postfix_char = c;

  return t;
}

inline size_t time_dot_formatter::get_postfix_width() const
{ return my_postfix_width; }

inline char time_dot_formatter::get_postfix_char() const
{ return my_postfix_char; }

inline
bool time_dot_formatter::prefix()
{
  if(my_current_value) return false;

  char fill = my_stream.fill();
  
  if(fill != my_prefix_char) my_stream.fill(my_prefix_char);

  my_stream << "        Percent       Time         Estimated Time"  << std::endl
            << "        Complete      Left          of Completion"   << std::endl
            << "        ========   ==========   =====================" << std::endl;

  if(fill != my_prefix_char) my_stream.fill(fill);

  if(my_current_value == 0)
  {
    time(&start_time);
    return true;
  }

  return true;
}

inline
bool time_dot_formatter::postfix()
{
  char fill = my_stream.fill();
  
  if(fill != my_postfix_char) my_stream.fill(my_postfix_char);

  if(interval) my_stream << std::setw(my_prefix_width) << my_postfix_char;

  if(my_postfix_width) my_stream << std::setw(my_postfix_width);
  
  my_stream << my_postfix;

  time(&current_time);

  size_t elapsed_time = current_time - start_time;

  my_stream << " (Time used: " << convert_time(elapsed_time, true) << ")" << std::endl;

  if(fill != my_postfix_char) my_stream.fill(fill);

  return true;
}

inline
bool time_dot_formatter::resync()
{
  if(my_current_value % time_check_count != 0) return true;

  time(&current_time);

  size_t elapsed_time = current_time - start_time;

  if(elapsed_time < 10) return true;

  if(!interval)
  {
    interval = std::max((size_t) 1, my_trigger_count / time_check_count / output_time_count);

    output_time(elapsed_time, true);

    return true;
  }

  if((my_current_value / time_check_count) % interval == 0)
    output_time(elapsed_time, false);

  return true;
}

inline void time_dot_formatter::output_time(size_t elapsed, bool first)
{
  double amount_complete = (double) my_current_value / (double) my_trigger_count;

  size_t time_remaining = (size_t) ((double) elapsed * (1.0 - amount_complete)
                                                     / amount_complete);

  if(time_remaining == 0) return;

  my_stream << std::setw(8) << ' ';

  my_stream << std::setw(7) << setiosflags(ios::right)
            << doub2str(amount_complete * 100, 4, 1, ios::fixed) 
            << "%    ";

  time_t end_time = current_time + time_remaining;

  string end = ctime(&end_time);

  // Remove day of week and ending \n
  end = end.substr(4, end.size() - 5);

  tm now  = *localtime(&current_time);
  tm then = *localtime(&end_time);

  // Blank out year only if year is not relevant
  if(now.tm_year == then.tm_year)
    std::fill(end.begin() + end.size() - 5, end.end(), ' ');

  // Blank out date only if is not relevant
  if(now.tm_mday == then.tm_mday && now.tm_mon == then.tm_mon)
    std::fill(end.begin(), end.begin() + 7, ' ');

  my_stream << std::setw(8) << setiosflags(ios::right)
            << convert_time(time_remaining, true)
            << resetiosflags(ios::right)
            << "     " << end << std::endl;
}

inline
bool time_dot_formatter::reset()
{
  return true;
}
