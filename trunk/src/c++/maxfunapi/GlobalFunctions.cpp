
#include "maxfun/maxfun.h"
#include "maxfunapi/Datatypes.h"

#include <sstream>

namespace SAGE   {
namespace MAXFUN {

//=======================================================================
//
//  str2paramname
// 
//=======================================================================
pair<string, string>
str2paramname(const string & combined_name)
{   
  pair<string, string> pair_names;
    
  string::size_type split_idx = combined_name.find(':');
   
  if(split_idx == string::npos)
  {       
    pair_names.first  = combined_name;
    pair_names.second = combined_name;
  }
  else
  {
    pair_names.first  = combined_name.substr(0,             split_idx);
    pair_names.second = combined_name.substr(split_idx + 1, combined_name.length() - split_idx - 1);
  }

  return pair_names;
} 

//=======================================================================
//
//  paramname2str   
//
//=======================================================================
string
paramname2str(const string & group_name, const string & param_name)
{
  ostringstream os;
  
  os << group_name << ":" << param_name;
    
  return os.str();
}  


}} // End namespace
