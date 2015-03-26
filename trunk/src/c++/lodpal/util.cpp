#include <string>
#include "LSF/parse_ops.h"
#include "lodpal/util.h"

namespace SAGE {

string pval(double p, size_t w, int prec)
{
  if(prec < 0)
    prec = w - 3 - 2;

  string pv = fp(p,w-3,prec);
  if   (!finite(p)) pv = pv.substr(0,w) + "   ";
  else if(p < 0.01) pv = pv.substr(0,w) + " **";
  else if(p < 0.05) pv = pv.substr(0,w) + " * ";
  else              pv = pv.substr(0,w) + "   ";
  return pv;
}

}
