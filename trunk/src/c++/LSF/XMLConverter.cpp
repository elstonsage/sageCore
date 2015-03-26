#include "LSF/XMLConverter.h"

namespace SAGE {
namespace LSF  {

//======================================
//
//  toXML(...)
//
//======================================
std::string
XMLConverter::toXML(const LSFBase * data)
{
  if(!data)
    return "";

  if(data->name() == "")
    return "";

  std::ostringstream txt;

  // Open bracket + name:

  txt << "<" << data->name();

  // Attributes:

  if(data->attrs())
    for(AttrList::const_iterator attr_itr = data->attrs()->begin(); attr_itr != data->attrs()->end(); ++attr_itr)
    {
      txt << " " << AttrNameMgr.query(attr_itr->first)             << "=\""
                 << convertAttributeVal(attr_itr->second.String()) << "\"";
    }

  txt << ">" << std::endl;

  // Child nodes:

  if(data->List())
    for(LSFList::const_iterator itr = data->List()->begin(); itr != data->List()->end(); ++itr)
      txt << toXML(*itr);
  
  // Close tag:

  txt << "</" << data->name() << ">" << std::endl;

  return txt.str();
}

//================================================
//
//  convertAttributeVal(...)
//
//================================================
std::string
XMLConverter::convertAttributeVal(const std::string & src_txt)
{
  std::ostringstream txt;

  for(size_t i = 0; i < src_txt.length(); ++i)
  {
    char c = src_txt[i];

         if(c == '<') txt << "&lt;";
    else if(c == '>') txt << "&gt;";
    else              txt << c;
  }

  return txt.str();
}

} // End namespace LSF
} // End namespace SAGE
