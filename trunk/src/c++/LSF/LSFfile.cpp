//
//  LSF FILE OBJECTS 1.0 -- Low level logical structure building block objects
//
//  Author: Kevin Jacobs (jacobs@darwin.cwru.edu)
//
//  History:   1.0   kbj  Initial implementation        Apr 5 1996
//
//  Copyright (c) 1996  R.C. Elston
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <string>
#include <ctype.h>
#include "LSF/parse_ops.h"
#include "LSF/LSF.h"
#include "LSF/LSFfactory.h"
#include "LSF/LSFfile.h"


// LSF_input constructor
LSF_input::LSF_input (istream &i, ostream &o) 
{ 
  lsf_file = &i; 
  free_file = false;
  err = &o;
  _lines = 1;
  if(!i) 
    setstate( LSFBase::badbit );
}

// LSF_input constructor
LSF_input::LSF_input (const char *f, ostream &o)
{ 
  _lines = 1;
  free_file = false;
  if(f && *f) 
  {
    lsf_file = new ifstream(f); 
    free_file = true;
  }

  err = &o;
  if( !lsf_file->good() ) 
    setstate( LSFBase::badbit );
}

LSF_input::~LSF_input () 
{
  if(lsf_file && free_file)
    delete lsf_file;
}

void LSF_input::input_to( LSFBase *root, bool typed )
{
  if( !lsf_file || !*lsf_file )
  {
    setstate( LSFBase::badbit );
    return;
  }

  char c;
  LSFBase *current = NULL;

  _path.erase(_path.begin(), _path.end() );
  assert( !root || root->List(TRUE) != NULL );

  if(root) _path.push_back(root);

  while( *lsf_file && !lsf_file->eof() )
  {
    _lines += kill_ws(*lsf_file);

    if( lsf_file->eof() ) break;

    c = lsf_file->peek();
    if( strchr("{}<>;#", c) )
    {
      lsf_file->get();

      switch(c)
      {
        case '#' : comment(); break;
        case ';' : break;
        case '<' :
        case '{' : if( current && current->List(TRUE) != NULL && 
                         _path.back() != current  )
                     _path.push_back( current );
                   else
                     (*err) << "Warning (LSF_input) Line # " << _lines 
                              << ": Miss-matched braces." << endl << flush;
                   break;
        case '>' :
        case '}' : if( _path.size() == 1 )              // Don't pop the root
                     (*err) << "Warning (LSF_input) Line # " << _lines 
                          << ": Miss-matched braces." << endl << flush;
                   else 
                   {
                     current = _path.back();
                     _path.pop_back();
                   }
                   break;
      }
      continue;
    }

    current = getLSF(typed);

    if(current)
    {
      if( _path.size() != 0 ) _path.back()->List(TRUE)->add( current );
      else (*err) << "Warning (LSF_input) Line # " << _lines 
                  << ": Internal error." << endl << flush;
    }
  }

  if( _path.size() != 1 )
    (*err) << "Warning (LSF_input) Line # " << _lines 
           << ": Unexpected end of file."   << endl << flush;


}

void LSF_input::error(char *m)
{
  if(m)
    (*err) << "Warning (" << name() << ") Line # " << _lines << ": " << m 
           << endl << flush;
  getString( *lsf_file, "{}<>=;\r\n", "");
  _lines += kill_ws(*lsf_file);

  do
  {
    if(!lsf_file->eof() && lsf_file->peek() == ',')
    {
      lsf_file->get();
      getString( *lsf_file, "{}<>=;\r\n", "");
      _lines += kill_ws(*lsf_file);
    }
    else
      break;
  } while( !lsf_file->eof() );
}  

void LSF_input::comment()
{
  getString( *lsf_file, "\r\n", "");
  _lines += kill_ws(*lsf_file);
}  

LSFBase *LSF_input::getLSF(bool typed)
{
  char c;

  def.erase( def.begin(), def.end() );

  string name, type;
 
  kill_ws(*lsf_file, " \t","");

  name = getString( *lsf_file, "=,#{}<>;\n\r", " \t" );

  if( name.length() && !(isalnum(name[0]) || strchr("$_",name[0])) )
  {
    error("Invalid LSF name");
    (*err) << "  Name = '" << name << "'" << endl << flush;
    return NULL;
  }

  kill_ws(*lsf_file, " \t","");

  c = lsf_file->peek();
  if (lsf_file->eof() || !strchr( "\n\r;#,={}<>", c ) )
  {
    error("Unexpected character");
    (*err) << "  Name = '" << name << "' Char = '" << c << "'" << endl;
    return NULL;
  }

  if(typed)
  {
    if( strchr("=", c) )
    {
      lsf_file->get();
      kill_ws(*lsf_file, " \t", "");
    }

    type = getString( *lsf_file, "=,#;{}<>\n\r", " \t" );

    if( type.length() && !isalpha(type[0]) )
    {
      error("Invalid LSF type");
      (*err) << "  Type = '" << type << "'" << endl << flush;
      return NULL;
    }
  }

  kill_ws(*lsf_file, " \t","");

  if (!lsf_file->eof())
  {
    c = lsf_file->peek();
    if( !strchr( "\n\r#;$,={}<>", c ) )
    {
      error("Unexpected character in LSF argument list");
      (*err) << "Name = '" << name << "' Type = '" << type 
               << "' Char = '" << c << "'" << endl;
      return NULL;
    }

    if(c == ',')                  // Allow for <name>==<value> syntax
    {
      lsf_file->get();
      _lines += kill_ws(*lsf_file); 
    }
    else
      kill_ws(*lsf_file, " \t","");

    while( !lsf_file->eof() )
    {
      getPair(*lsf_file);
      kill_ws(*lsf_file, " \t","");

      if(!lsf_file->eof() && lsf_file->peek() == ',')
      {
        lsf_file->get();
        _lines += kill_ws(*lsf_file);
        if( lsf_file->eof() )
          error("Unexpected end of file");
      }
      else
        break;
    }
  }
  return Factory->build(name.c_str(), type.c_str(), &def);
}

void LSF_input::getPair( istream &i, const char *delim,
                                     const char *ws )
{
  string left, right;
  char c;

  kill_ws(i, ws,"");
  left = getString( i, delim, ws );
  kill_ws(i, ws, "");
  if(i.eof() || i.peek() != '=')
  {
    if( left.length() )
    {
      def.set(left,right);
    }
    return;
  }

  i.get();
  kill_ws(i, ws,"");
  c=i.peek();

  if (isalpha(c) || strchr("\"$_./!?@#%^&*()",c) )
  {
    right = getQString( i, delim, ws, "\r\n", &_lines );
    def.set(left, right);
    return;
  }

  right = getString( i, delim, ws );

  if( !right.length() )
  {
    def.set(left,right);
    return;
  }
  c=right[0];
  
 // Floats have . or e/E 
  if ( strchr(".",c) || 
     ((strchr("-+",c) || isdigit(c)) 
      && (long)right.find_first_of(".eE",0) > 0)) 
  {
    double d = str2doub(right);
    def.set(left, d);
  }
  else if( strchr("-+",c) || isdigit(c) )
  {
    int i = str2long(right);
    def.set(left, i);
  }
  else
  {
    def.set(left, right);
  }
}

// Load an LSF file and insert it into a new LSF base with name = desc
LSFBase *loadLSFfile(const string& file, const string& desc,
                     ostream& errors, bool typed)
{
  LSF_input load_lsf(file.c_str(), errors);
  if(!load_lsf)
    return NULL;
  assert( load_lsf.good() );
  LSFBase *dest = new LSFBase(desc.c_str());
  load_lsf.input_to(dest, typed);
  return dest;
}

void SimpleLSFDump::output(ostream &o, LSFBase *i)
{
  out = &o;
  level = 0;
  i->Accept(this);
}

void SimpleLSFDump::def(LSFBase *c) { output(c); }

void SimpleLSFDump::Process(LSFBase *c)
{
  output(c);
  if(c->List())
  {
    level+=2;
    c->AcceptAll(this);
    level-=2;
  }
}

void SimpleLSFDump::output(LSFBase *c)
{
  indent(level);

  string type;
  if(Factory)
    type = Factory->type_map( c->type() );
  
  (*out) << type  << setw(20 - level - Factory->type_map(c->type()).length() ) 
         << right << ": "      << c->name();

  if( c->List() )
    (*out) << setw(20 - c->name().length() ) << right
           << " " << "[LIST " << setw(3) << c->List()->size() << ']';

  (*out) << endl;

}
void SimpleLSFDump::indent(int n) { for(int i=0; i<n; i++) (*out) << ' '; }

LSFDump::LSFDump(char *name) 
{ 
  out = new ofstream(name); 
  free_stream = true; 
}

void LSFDump::dump(LSFBase *i)
{
  level = 0;
  i->Accept(this);
}

void LSFDump::def(LSFBase *c) { output(c); }

void LSFDump::Process(LSFBase *c)
{
  output(c);
  if( c->List() && c->List()->size() )
  {
    indent(level);
    (*out) << '{' << endl;

    level+=2;
    c->AcceptAll(this);
    level-=2;

    indent(level);  
    (*out) << '}' << endl;
  }
}

void LSFDump::output(LSFBase *c)
{
  string val;
  bool quoted;

  string type;

  if(Factory)
    type = Factory->type_map(c->type());

  indent(level);
  (*out) << c->name();

  
  if((type.length() && !c->name().size())
    || (c->name().size() && type != "base"))
    (*out) << '=' <<  type;

  int mark = level + c->name().length() + 1;
  int pos  = mark  + type.length();
  int len1, len2;   // len1 - length of attrib name + format chars,
                    // len2 - length of atribute value as string.   

  if(c->attrs())
  {
    LSFBase::AList::iterator i;
    for( i=c->attrs()->begin(); i!= c->attrs()->end(); i++)
    {
      len2 = 0;                                    // set initially to 0
      quoted = false;                              // default - no quotes

      (*out) << ',';                               // output separator
      ++pos;

      len1 = c->attrs()->Name( (*i).first ).length();  // set len1

      val = (*i).second.String();      

      if( val.length() )
      {
        len1++;
        len2 = val.length();                              // set len2

        if( ( (*i).second.Type()==AttrVal::Attr_String && 
           (int)val.find_first_of(" \"\\,.()={};\n\r\t") >= 0) ||
            mark+len1+len2 > 78 )
        {
          quoted = true;
          len2+=2;                                 // Add quotes to len1    
        }
      }

      if( !quoted && pos + len1 + len2 > 78 && pos != mark )
      {
        (*out) << endl;
        indent( (pos = mark) );
      }
      else if(quoted && pos + len1 + len2 > 78)
      {
        int m, n = 0;             // m - current pos, n - last separator

        for(m = 5; m < (len2-2) && pos + len1 + m + 2 < 78; ++m)   // Find separator
          if(strchr(" \t\n;.,*-+/|?", val[m]))
            n = m; 

        if ((len2 - n < len2/4 || n == 0) && pos != mark)  // If there isn't
        {                                                  // much left
          (*out) << endl;
          indent( (pos = mark) );
        }
      }  

      if(len1)
      {      
        (*out) << c->attrs()->Name ((*i).first);
        pos+= len1;
      }

      if(len2 || !len1)
      {
        (*out) << '=';
        ++pos;
      }

      if(len2)
      {
        if (quoted) pos = printQuoted(*out,val,pos,mark);
        else        
        {
          (*out) << val;
          pos += len2;
        }
      }
    }
  }
  (*out) << endl;
}

void LSFDump::indent(int n) { for(int i=0; i<n; i++) (*out) << ' '; }
