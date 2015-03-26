#include "maxfunapi/CovarianceMatrix.h"

using namespace std;

namespace SAGE   {
namespace MAXFUN {

//===========================================================
//  CONSTRUCTOR
//===========================================================
CovarianceMatrix::CovarianceMatrix()
{
  my_available = false;
  my_values    . clear();
  my_names     . clear();
}

//===========================================================
//  COPY CONSTRUCTOR
//===========================================================
CovarianceMatrix::CovarianceMatrix(const CovarianceMatrix & other)
{
  copy(other);
}

//===========================================================
//  operator=()
//===========================================================
CovarianceMatrix& 
CovarianceMatrix::operator= (const CovarianceMatrix& other)
{
  if(&other != this)
  {
    copy(other);
  }

  return *this;
}

//===========================================================
//  copy()
//===========================================================
void 
CovarianceMatrix::copy(const CovarianceMatrix & other)
{
  my_available = other.my_available;
  my_values    = other.my_values;
  my_names     = other.my_names;
}

//===========================================================
// isAvailable()
//===========================================================
bool CovarianceMatrix::isAvailable() const { return my_available; }

//===========================================================
// getSize()
//===========================================================
int CovarianceMatrix::getSize() const { return my_names.size(); }

//===========================================================
// getNames()
//===========================================================
const vector<pair<string, string> > & CovarianceMatrix::getNames() const { return my_names; }

//===========================================================
// getCovariance()
//===========================================================
double CovarianceMatrix::getCovariance(int param_id1, int param_id2) const { return my_values[param_id1][param_id2]; }

//===========================================================
// setCovariance()
//===========================================================
int 
CovarianceMatrix::setCovariance(int param_id1, int param_id2, double val)
{
  my_values[param_id1][param_id2] = val;

  return 0;
}

//===========================================================
//  inputData()
//===========================================================
void 
CovarianceMatrix::inputData(const Maxfun_Data & data)
{
  // 0. Set up local variables:

	int N = data.ne();

  // 1. Populate name list:

	my_names.clear();

	for(size_t i = 0; i < (size_t)data.nt(); i++)
	  if(data.ist(i) != 4)
	    my_names.push_back(str2paramname(data.label(i)));

  // 2. Populate matrix:
  
        my_values.clear();

        if(0 < N && N <= data.nt())
        {
          my_values.resize(N);
  
          for(int i = 0; i < N; i++)
            for(int j = 0; j < N; j++)
              my_values[i].push_back(data.av(i, j));

  // 3. Set available:

          if(N)
            my_available = true;
        }
}

} 
} 


