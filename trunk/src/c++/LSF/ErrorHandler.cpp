#include <iostream>
#include "LSF/ErrorHandler.h"
#include "LSF/LSF.h"

#ifdef NAMESPACE_STD
NAMESPACE_STD
#endif

namespace SAGE  // Needed to avoid namespace collision
{

ErrorHandler* GlobalError;

void ErrorHandler::out_error(ostream* o, error_t e, const string& s)
{
  size_t indent = 0, _ind = 0;

  if (!o || !o->good()) return;

  (*o) << endl << '%';

  if (pname != "")
  {
    (*o) << pname << '-';
    _ind += pname.length() + 2;
  }

  switch (e)
  {
    case Warning     : (*o) << "W: ";
                       _ind += 3;
                       break;
    case Information : (*o) << "I: ";
                       _ind += 3;
                       break;

    case Error       : (*o) << "E: ";
                       _ind += 3;
  }

  size_t pos = 0, l_begin;
  int l_blank=-1;

  while (pos < s.length())
  {
    for (size_t j = 0; j < indent; j++) (*o) << ' ';
    indent = _ind;

    l_begin = l_blank+1;
    if (indent + s.length() - pos < 78)
    {
      (*o) << s.substr(l_begin, s.length()) << endl;
      pos = s.length();
    }
    else
    {
      for (; pos - l_begin + indent <= 78 && pos < s.length() - 10; pos++)
        if (s[pos] == ' ') l_blank=pos;
      (*o) << s.substr(l_begin, l_blank - l_begin) << endl;
    }
  }

}

}
