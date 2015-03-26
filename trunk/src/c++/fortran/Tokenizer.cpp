#include <fstream>
#include <math.h>
#include "LSF/Attr.h"
#include "LSF/LSF.h"
#include "LSF/LSFfactory.h"
#include "fortran/Tokenizer.h"
#include "fortran/Formatter.h"

Tokenizer::Tokenizer(BaseFormatter* _p, istream& _i)
    : i(&_i), p(_p), _eof(false), count(0)
{
  if (p) p->grab();
  if (!i || i->fail() || !p)
    setstate( LSFBase::badbit );
}

Tokenizer::Tokenizer(BaseFormatter* _p, char* _n)
    : p(_p), _eof(false), count(0)
{
  if (p) p->grab();
  if (_n && *_n) i = new ifstream(_n);
  if (!i || i->fail() || !p)
    setstate( LSFBase::badbit );
}

Tokenizer::~Tokenizer()
{
  if (p) p->release();
}

bool Tokenizer::eol() const
{
  if(p)
    return (buffer.size() <= (unsigned) p->get_pos());
  else
    return true;
}

bool Tokenizer::eof() const
{
  return (_eof || (eol() && i->eof()));
}

void Tokenizer::read_line()
{
  if (!i->fail())
  {
    if(i->eof())
    {
      buffer = "";
      _eof = true;
    }
    else
    {
      getline( *i, buffer );
      if(i->eof() && buffer == "")
        _eof = true;
    }
  }
  else
    setstate(badbit);

  ++count;
}

Tokenizer::token Tokenizer::get_token()
{
  token  tok(action::None, 0);
  action act;

  while (true)                
  {
    if (_eof)
    {
      tok.first = action::EndOfFile;
      return tok;
    }
    act = p->get_action();    
    switch (act.Action)
    {
      case action::move   : break;

      case action::flush  : read_line();
                            break;

      case action::read   : return read_token(act, tok);

      case action::msg    : tok.first  = act.data_type;
                            tok.second = act.action_string;
                            return tok;

      case action::output :

      default             : break;
    }
  }
}

Tokenizer::token Tokenizer::read_token(action& act, token& tok)
{
  string s;

  int i = act.end + 1;

  if ((unsigned) i > buffer.length()) i = buffer.length();	          

  if (i < act.start) act.start = i;				          

  convert_string( s, buffer.substr(act.start, i - act.start),	          
                     act.data_type, act.flag );			          

  tok.second = s;						          

  switch (act.data_type)					          
  {								          
    case action::Integer : tok.second = tok.second.Int();				          
                           break;							          

    case action::Float   : convert_float(act, tok);
                           break;							          

    case action::String  :					          

    default              : break;							          
  }								          
  tok.first = action::None;						          

  return tok;							          
}

void Tokenizer::convert_string(string&   dest, const string& source,
                               const d_type dtype, const int     flag) 
{
  if ( dtype == action::Integer || dtype == action::Float )
  {
    if (flag & action::BZflag)
    {
      dest = source;
      for (unsigned int i = 0; (i < dest.size()); i++)
        if (strchr(" \t", dest[i])) dest[i] = '0';
    }
    else
    {
      for (unsigned int i = 0; (i < source.size()); i++)
        if(!strchr(" \t", source[i])) dest += source[i];
    }
  }
  else
    dest = source;
}

// Convert float follows the Fortran conversion process.  If a decimal point
//   is not located in the input, the input is converted based on the 
//   number of decimal places there 'should' be.

void Tokenizer::convert_float(action& act, token& tok)
{
  if (act.decimal > 0 && 
      (unsigned long) tok.second.String().find_last_of(".") == 
      (unsigned long) -1)      
    tok.second = tok.second.Real() / pow(10.0, act.decimal);
  else
    tok.second = tok.second.Real();
}

OutputTokenizer::OutputTokenizer(BaseFormatter* _p, ostream& _o)
    : o(&_o), p(_p)
{
  if (p) p->grab();
  if (!o || o->fail() || !p)
    setstate( LSFBase::badbit );
}

OutputTokenizer::OutputTokenizer(BaseFormatter* _p, char* n) : p(_p)
{
  if (p) p->grab();
  if (n && *n) o = new ofstream(n);
  if (!o || o->fail() || !p)
    setstate( LSFBase::badbit );
}

OutputTokenizer::~OutputTokenizer() { if (p) p->release(); }

void OutputTokenizer::write_line()
{
  if (!o->fail())
  {
    *o << buffer << endl; 
    buffer = "";
  }
  else
    setstate(badbit);
}

OutputTokenizer::token OutputTokenizer::put_token(const token& tok)
{
  action act;
  act.flag = 0;
  token t(action::None, 0);

  while (true)
  {
    act = p->get_action();

    switch (act.Action)
    {
      case action::move   : break;

      case action::flush  : write_line();
                            break;

      case action::read   : if(tok.second.Type() != AttrVal::Attr_None)
                              write_token(act, tok);
                            return t;

      case action::msg    : t.first  = act.data_type;
                            t.second = act.action_string;
                            return t;

      case action::output : if (buffer.length() <= (unsigned) act.end)                         
                              buffer.resize(act.end+1, ' ');
                            copy(act.action_string.begin(),
                                 act.action_string.end(),
                                 buffer.begin() + act.start);
                            break;

      default             : break;
    }
  }
}

void OutputTokenizer::write_token(action& act, const token& tok)
{
  string s;

  unsigned int size = act.end - act.start + 1;


  if(act.data_type != action::String)
  {
    if(act.flag & action::BZflag)
      s = tok.second.String(size, act.decimal, ios::fixed, '0');
    else
      s = tok.second.String(size, act.decimal, ios::fixed);

    if((signed) s.find_last_not_of(" ") == -1) 
      s = "";

    if(s.size() > size)
    {
      s = s.substr(0, size);
      asterik(s);
    }
  }
  else
    s = tok.second.String().substr(0, size);

  if (s.size() && buffer.length() < (unsigned) act.start + s.size())
    buffer.resize(act.start+s.size(), ' ');

  copy(s.begin(), s.end(), buffer.begin() + act.start);
}
