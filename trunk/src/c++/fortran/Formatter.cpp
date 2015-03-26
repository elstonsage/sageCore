#include <fstream>
#include "LSF/Attr.h"
#include "LSF/LSF.h"
#include "LSF/LSFfactory.h"
#include "fortran/Formatter.h"

// set_format makes the format string for the parser equal to the input
//     string and resets the parser to the beginning.

void FortranFormatter::set_format (const string& s)
{
  fmt_str = s;
  fmt_length = s.size();
  reset();
}

// reset - resets the parser to the beginning of the format string

void FortranFormatter::reset()
{
  fmt_pos   = 0;
  state     = 0;
  flag      = 0;
  input_pos = 0;

  while (!st.empty()) st.pop();
  st.push(stack_elt(0,1,0));
}

// get_action returns the next action based upon the format string
//     state 0: beginning of format string
//     state 1: Formatting element
//     state 2: In between elements (bound)
//     state 3: Determine if next state is state 1 or 2

FortranFormatter::action FortranFormatter::get_action()
{
  action a;

  while (true)
  {
    ignore_ws();
    switch (state)
    {
      case 0: if (fmt_chr()=='(')
                fmt_pos++;

              if (fmt_chr() == '/' || fmt_chr() == ')')
                state = 2;
              else
                state = 1;

              break;

      case 1: state = 2;
              if (isdigit(fmt_chr()))
                get_fmt_digit();
              else
                if (next_format_elt(a))
                  return a;
              break;

      case 2: if (next_bound_elt(a))
                return a;
              break;

      case 3: if(strchr("),/", fmt_chr())==NULL)
                state = 1;
              else
                state = 2;
    }
  }
}

// FORTRAN_ERROR returns error message - This probably needs restructuring.

void FortranFormatter::FORTRAN_ERROR(action& a, const char* err)
{
  char err_str[2000];

  sprintf(err_str, "Error: %s in format string \"%s\" at position %i",
          err, fmt_str.c_str(), fmt_pos);

  a.action_string = err_str;
  a.data_type     = action::Error;
  a.Action        = action::msg;
}

// next_format_elt - Parses a block of text

bool FortranFormatter::next_format_elt(action& a)
{
  a.data_type = action::None;
  switch (fmt_chr())
  {
    case 'I': case 'i': a.data_type = action::Integer;
                        return get_value(a);

    case 'A': case 'a': a.data_type = action::String;
                        return get_value(a);

    case 'F': case 'f': a.data_type = action::Float;
                        return get_value(a);

    case 'T': case 't': return get_tab(a);

    case 'B': case 'b': return get_BZflag(a);

    case '\''         : return get_quote(a);

    default:
        FORTRAN_ERROR(a, "Unknown type");
        return 1;
  }
}

// next_bound_elt figures out what to do between elements

bool FortranFormatter:: next_bound_elt(action& a)
{
  state = 3;

  char c = fmt_chr();

  switch (c)
  {
    case ',' :
               if (!popstack())
                 fmt_pos++;
               return 0;

    case '/' : if (!popstack())
               {
                 fmt_pos++;
                 input_pos = 0;

                 a.Action = action::flush;
                 return 1;
               }
               return 0;

    case '\0':
    case ')' : if ((--st.top().rpt) == 0)
                 return get_end_fmt(a);
               fmt_pos = st.top().start;
               return 0;

    default  : FORTRAN_ERROR(a, "Misplaced character");
               return 1;
  }
}

// get_value reads off the values associated with an int, float or string,
//   and returns the action associated.

bool FortranFormatter::get_value(action& a)
{
  fmt_pos++;								    

  int value = get_number();

  if (value < 0)
  {
    FORTRAN_ERROR(a, "Missing value");
    return 1;
  }									   

  if (a.data_type == action::Float)					   
  {									   
    if (fmt_chr() != '.')
    {									   
      a.decimal = -1;
    }
    else
    {
      fmt_pos++;								   

      if ((a.decimal = get_number()) < 0)
      {									   
        FORTRAN_ERROR(a, "Second value in floating point missing");		   
        return 1;								   
      }
    }
  }									   
  else
  {
    a.decimal = 0;
  }

  a.Action = action::read;
  a.start  = input_pos;						   
  a.end    = (input_pos += value) - 1;
  a.flag  ^= flag;						   

  return 1;								   
}

// get_tab parses a tab element (ex 'T5')

bool FortranFormatter::get_tab(action& a)
{
  int dir      =  1;
  int base_pos = -1;

  fmt_pos++;

  switch (fmt_chr() )
  {
    case 'l': case 'L': dir = -1;

    case 'r': case 'R': base_pos = input_pos;
                        fmt_pos++;

    default           : break;
  }

  int tab_amt = get_number();

  if (tab_amt < 0)
  {
    FORTRAN_ERROR(a, "Number must follow tab");
    return 1;
  }

  input_pos = base_pos + (dir * tab_amt);
  return 0;
}

// get_BZflag sets the zeros flag.

bool FortranFormatter::get_BZflag(action& a)
{
  fmt_pos++;

  switch(fmt_chr())
  {
    case 'N': case 'n': flag &=  !action::BZflag;
                        break;

    case 'Z': case 'z': flag |= action::BZflag;
                        break;

    default           : FORTRAN_ERROR(a, "Incorrect option");
                        return 1;
  }

  fmt_pos++;
  return 0;
}

// get_end_fmt deals with ) and \0 terminators

bool FortranFormatter::get_end_fmt( action& a )
{  
  if (st.size() == 1)  
  {
    if(fmt_length - 1 > fmt_pos )
    {   
      FORTRAN_ERROR(a, "Early end of format");
      return 1;   
    }   

    a.Action    = action::msg; 
    a.data_type = action::EndFormat;  

    reset(); 

    return 1;
  }

  if (!(st.top().end)) 
    fmt_pos++; 

  st.pop();  

  return 0;  
}  

bool FortranFormatter::get_quote (action& a)
{
  unsigned int i = fmt_pos + 1;

  for (; i < fmt_length && fmt_str[i] != '\''; i++);

  if(i == fmt_length || fmt_chr() != '\'')
  {
    FORTRAN_ERROR(a, "Unterminated literal string");
    return 1;
  }

  a.Action        = action::output;
  a.action_string = fmt_str.substr(fmt_pos + 1, i - fmt_pos - 1);
  a.start         = input_pos;
  a.end           = (input_pos += a.action_string.length()) - 1;
  a.decimal       = 0;

  fmt_pos = i + 1;

  return 1;
}

// get_fmt_digit reads in elements which begin with numbers (includes
//     pushing to the stack of repeating elements)

void FortranFormatter::get_fmt_digit()
{  
  int repeat = get_number();

  switch (fmt_chr())
  {
    case 'X' :
    case 'x' : input_pos += repeat;
               fmt_pos++;
               return;

    case '(' : st.push( stack_elt( ++fmt_pos, repeat, 0 ) );    

               if (fmt_chr() != '/' && fmt_chr() != ')')
                 state = 1;

               return;

    default  : state = 1;
               st.push( stack_elt(   fmt_pos, repeat, 1 ) );
  }
}

// get_number returns the value of the number at fmt_pos and moves fmt_pos to
//   the character after the end of the number

int FortranFormatter::get_number()
{
  unsigned int t = fmt_pos;

  for(; t < fmt_length && isdigit(fmt_str[t]); t++);

  if (fmt_pos == t)
    return -1;

  int value = atoi(fmt_str.c_str() + fmt_pos);

  fmt_pos = t;

  return value;
}

// popstack checks the top of the stack for repeats and either (if done)
//     pops the last element, or (if needs to repeat) returns the parser
//     to the beginning of the repeated string.

bool FortranFormatter::popstack()
{
  if (st.top().end)
    if ((--st.top().rpt) != 0)
    {
      fmt_pos = st.top().start;
      return 1;
    }
    else
      st.pop();

  return 0;
}
