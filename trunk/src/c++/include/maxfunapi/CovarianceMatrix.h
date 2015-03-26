#ifndef COVARIANCEMATRIX_H
#define COVARIANCEMATRIX_H

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include "util/get_mem.h"
#include "maxfun/maxfun.h"
#include "maxfunapi/Datatypes.h"

using namespace std;

namespace SAGE   {
namespace MAXFUN {

/** \class CovarianceMatrix
 *  \brief Stores information about the variance-covariance matrix.
 *
 */
class CovarianceMatrix
{
  public:
	friend class ParameterMgr;
	friend class Parameter;
	friend class Results;
	MEM_FRIEND(SAGE::MAXFUN::CovarianceMatrix);

	// Default constructor
	CovarianceMatrix();

	// Copy constructor
	CovarianceMatrix(const CovarianceMatrix &);

	// Assignment operator
	CovarianceMatrix & operator= (const CovarianceMatrix &);

	// copy operation
	void copy(const CovarianceMatrix &);
	
	///
	/// Returns whether or not the covariance matrix is available.
	/// \retval true The matrix \b is available.
	/// \retval false The matrix is \b not available.
	bool isAvailable() const;

	///
	/// Returns the size of the matrix N, where the matrix is of N x N size.
	int getSize() const;

	///
	/// Returns a vector of names of parameters with available covariances.
	/// Please note that that pair is made up of the group name and parameter name (in that order).
	/// Also, please note that the index of the the name corresponds to its index in the matrix of
	/// covariance values.
	const vector<pair<string, string> > & getNames() const;

	///
	/// Returns the covariance between param1 and param2.
	/// Please note that the parameters are identified by their VarIndex(), not by their Index().
	double getCovariance(int param_id1, int param_id2) const;

	///
	/// Sets the indicated covariance to a new value.
	/// \param param_id1 The VarIndex() of the first parameter.
	/// \param param_id2 The VarIndex() of the second parameter.
	/// \param val The value to which the variance should be set.
	int setCovariance(int param_id1, int param_id2, double val);



  protected:

	// input operation
	void inputData(const Maxfun_Data &);

  private:
        bool                          my_available;
	vector<pair<string, string> > my_names;
	vector<vector<double> >       my_values;
};

}} /// End namespace

MEM_COUNT_BEGIN(SAGE::MAXFUN::CovarianceMatrix)
{
          size_t v = 0;
          
          for(vector<pair<string, string> >::const_iterator i = t.my_names.begin(); i != t.my_names.end(); ++i)
            v += i->first.size() + i->second.size();
            
          for(vector<vector<double> >::const_iterator i = t.my_values.begin(); i != t.my_values.end(); ++i)
            v += get_mem(*i);
            
          v += sizeof(SAGE::MAXFUN::CovarianceMatrix);
          
          return v;
        }
MEM_COUNT_END


#endif

