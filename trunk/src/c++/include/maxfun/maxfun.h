#ifndef MAXFUN_H
#define MAXFUN_H

#include <cstring>
#include <stdlib.h>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "util/get_mem.h"
#include "numerics/isnan.h"
#include "globals/config.h"

using namespace std;

namespace SAGE {

/// This advance declaration of the SAGE::MAXFUN::Maximizer class is necessary in order to get the 
/// compiler to correctly understand the 'friend class SAGE::MAXFUN::Maximizer' line.
namespace MAXFUN { class Maximizer; }

using std::vector;
using std::string;

// C++ Binding for MAXFUN, a Fortran subroutine which performs slow, but
// extremely reliable maximization of complex functions (such as likelihood
// surfaces in multivariate maximum likelihood estimation)

class MaxFunction
{
  public:

    typedef vector<double>   parameter_vector;

    virtual ~MaxFunction();

    inline double operator()(parameter_vector& theta)
    {
      return evaluate(theta);
    }

    inline int depar(parameter_vector& theta)
    {
      return update_bounds(theta);
    }

    size_t nfe;

  protected:

    virtual double evaluate(parameter_vector& theta)      = 0;
    virtual int    update_bounds(parameter_vector& theta) = 0;

};

/** The maxfun externally available data.
 *  This is mainly of use for copying out final estimates after a maxfun run
 *  has completed.  Note that it only provides const access to the data
 *  structures.  This is to avoid inadvertantly placing Maxfun in an
 *  inconsistent state.
 */
class Maxfun_Data
{
  public:

    friend class Maxfun;

    MEM_FRIEND(SAGE::Maxfun_Data);

    // Constructors

    Maxfun_Data();
    Maxfun_Data(const Maxfun_Data&);

    Maxfun_Data& operator=(const Maxfun_Data&);

    // Destructor
    ~Maxfun_Data() { }
  
    // MAXFUN parameters
    int NP()                const   { return param_NP;  }
    int NPV()               const   { return param_NPV; }

    // MAXFUN data bindings.  See MAXFUN manual for variable descriptions

    // Input Variables to Maxfun

    // due to JA, for score test purposes

    // Convergence Criteria and constants
    double   epsc1()                   const  { return maxf1_.epsc1;     }
    double   epsc2()                   const  { return maxf1_.epsc2;     }
    double   epsc3()                   const  { return maxf1_.epsc3;     }
    double   epsd()                    const  { return maxf1_.epsd;      }
    double   epst()                    const  { return maxf1_.epst;      }
    double   yota()                    const  { return maxf1_.yota;      }

    // Maximization/Computation options       
    int      ihit()                    const  { return maxf1_.ihit;      }
    int      ixvc()                    const  { return maxf1_.ixvc;      }
    int      maxit()                   const  { return maxf1_.maxit;     }
    int      method()                  const  { return maxf1_.method;    }

    // Parameter options                      
    int      nt()                      const  { return maxf1_.nt;        }
    int      istin(int i)              const  { return maxf1_.istin[i];  }
    double   stpin(int i)              const  { return maxf1_.stpin[i];  }
    double   thin(int i)               const  { return maxf1_.thin[i];   }
    double   thl(int i)                const  { return maxf1_.thl[i];    }
    double   thu(int i)                const  { return maxf1_.thu[i];    }

    // Read the i'th label
    const string& label(int i)         const  { return maxflb_.label[i]; }

    // Interface to parameter options
    int        param_count()           const  { return maxf1_.nt;        }
    double     initial_estimate(int i) const  { return maxf1_.thin[i];   }
    double     lower_bound(int i)      const  { return maxf1_.thl[i];    }
    double     upper_bound(int i)      const  { return maxf1_.thu[i];    }

    // Automatic Output
    double     value()                 const  { return f;        }
    double     param(int i)            const  { return theta[i]; }
    int        last_error()            const  { return lfl;      }
    int        evaluations()           const  { return nfe;      }

    // Optional Output Information
    double     av(int i, int j)        const  { return maxf2_.av[j*NP()+i]; }
    double     h(int i, int j)         const  { return maxf2_.h[j*NP()+i];  }
    double     g(int i)                const  { return maxf2_.g[i];         }
    double     stde(int i)             const  { return maxf2_.stde[i];      }
    double     stp(int i)              const  { return maxf2_.stp[i];       }
    double     thpr(int i)             const  { return maxf2_.thpr[i];      }
    double     cth(int i)              const  { return maxf2_.cth[i];       }
    int        ist(int i)              const  { return maxf2_.ist[i];       }

    double     difmax()                const  { return maxf2_.difmax;       }
    double     erm()                   const  { return maxf2_.erm;          }
    double     fch()                   const  { return maxf2_.fch;          }
    double     fpr()                   const  { return maxf2_.fpr;          }
    int        idif()                  const  { return maxf2_.idif;         }
    int        igage()                 const  { return maxf2_.igage;        }
    int        igfl()                  const  { return maxf2_.igfl;         }
    int        impbnd()                const  { return maxf2_.impbnd;       }
    int        it()                    const  { return maxf2_.it;           }
    long       ivage()                 const  { return maxf2_.ivage;        }
    int        ivfl()                  const  { return maxf2_.ivfl;         }
    int        nb()                    const  { return maxf2_.nb;           }
    int        nd()                    const  { return maxf2_.nd;           }
    int        ne()                    const  { return maxf2_.ne;           }
    int        ni()                    const  { return maxf2_.ni;           }
    int        nsurf2()                const  { return maxf2_.nsurf2;       }
    int        nv()                    const  { return maxf2_.nv;           }

    int        ireturn()               const  { return maxf2_.ireturn;      }
    // Status labels                   
    string     maxfst(int i)           const  { return maxfst_[i];          }

    /// @name Debugging
    //@{

    ///
    /// Returns the memory usage, in bytes, ROUGHLY.
  
    //@}

  protected:

    // Data members

    // Local storage of function value evaluted at theta
    double f;

    // Error flag storage for results of maximization.  See MAXFUN manual.
    int lfl;

    // Number of times _fun has been evaluated.  Reset only by calling reset()
    int nfe;  

  public: // SGROSS

    // Array of parameter estimates
    vector<double>      theta;

    // due to JA
    vector<double>      score_vec;
    vector<double>      score_vec2;
    vector<double>      score_stepsize;
    vector<int>         need_score;
    vector<int>         score_deriv_conv;
    vector<std::string> flat_dir;

  protected: // SGROSS
    /// Clear out the results of prior maxfun runs from maxf2_
    void clear_results();  

  private:

    static char maxfst_[][80];

  public:

    struct maxf1_struct
    {
       MEM_FRIEND(SAGE::Maxfun_Data::maxf1_struct);

       vector<double> thin;
       vector<double> thl;
       vector<double> thu;
       vector<double> stpin;

       double epsd;
       double yota;
       double epst;
       double epsc1;
       double epsc2;
       double epsc3;

       vector<int> istin;

       int nt;
       int maxit;
       int method;
       int ixvc;
       int ihit;
    };

    struct maxf2_struct
    {
       MEM_FRIEND(SAGE::Maxfun_Data::maxf2_struct);

       vector<double> thpr;
       vector<double> cth;
       vector<double> stp;
       vector<double> g;
       vector<double> h;
       vector<double> v;
       vector<double> av;
       vector<double> stde;
       vector<double> ult;
       vector<double> diag;
       vector<double> pdir;

       double erm;
       double fpr;
       double fch;
       double difmax;
       double gtg;
       double ptg;
       double tstep;

       vector<int> ist;

       int ne;
       int nd;
       int ni;
       int nb;
       int nv;
       int impbnd;
       int it;
       int nsurf2;
       int igfl;
       int ivfl;
       int igage;
       int ivage;
       int idif;

       int ireturn;
    };

    struct maxflb_struct
    {
      MEM_FRIEND(SAGE::Maxfun_Data::maxflb_struct);

      vector<string> label;
    };

  private:

    int                              param_NP;
    int                              param_NPV;
    int                              param_NPV2; // due to JA
    maxf1_struct                     maxf1_;
    maxf2_struct                     maxf2_;
    maxflb_struct                    maxflb_;
};

class Maxfun
{
  friend class SAGE::MAXFUN::Maximizer;

  private:  

    // Initialization for MAXFUN runtime
    void init();

    // Copy construction is forbidden
    Maxfun(const Maxfun&);

  public:

    // Constructor:  A function must be supplied.  If d is not specified then
    // a dummy routine is used to validate parameters.
    Maxfun(MaxFunction& mfun)
    {
      fun = &mfun;
      init();
      reset();
    }

    // Destructor
    ~Maxfun() {}
  
    int run()
    {
      return maxfun_();
    }

    int evaluate(vector<double>& theta, double& f, int& nfe2, int& lex)
    {
      lex = 0;

      depar(theta, lex);

      if(lex)
        return lex;

      int n = fun->nfe;

      f = (*fun)(theta);

      nfe2 += fun->nfe - n;

      if(SAGE::isnan(f))
        lex = 1;

      return lex;
    }

    int evaluate(double& f, int& nfe, int& lex)
    {
      return evaluate(theta, f, nfe, lex);
    }

    int depar(vector<double>& theta, int& lex)
    {
      lex = fun->depar(theta);
      return lex;
    }

    void reset()
    {
      f           = std::numeric_limits<double>::quiet_NaN();
      lfl = 0;
      nfe = 0;
    }

    // introduced by JA for score test 
    void set_score_calc(int param_indic, int val) 
    {
      my_data.need_score[param_indic] = val; 
      return;
    }

    void reset_score_vec()
    {
      unsigned currentsize = my_data.score_vec.size();

      for (unsigned i = 0; i != currentsize; i++)
      {
        my_data.score_vec[i] = \
        std::numeric_limits<double>::quiet_NaN();
      }

      return ;
    }

    void reset_score_deriv_conv()
    {
      unsigned currentsize = my_data.score_deriv_conv.size();
      for (unsigned i = 0; i != currentsize; i++)
      {
        my_data.score_vec[i] = 9;
      }

      return ;
    }

    // due to JA top help djb with score test execution
    double obtain_score();

    // MAXFUN parameters
    int NP()                const   { return my_data.param_NP;  }
    int NPV()               const   { return my_data.param_NPV; }

    // MAXFUN data bindings.  See MAXFUN manual for variable descriptions
    // Required Input Section

    // Convergance Criteria and constants
    double&   epsc1()             { return my_data.maxf1_.epsc1;     }
    double&   epsc2()             { return my_data.maxf1_.epsc2;     }
    double&   epsc3()             { return my_data.maxf1_.epsc3;     }
    double&   epsd()              { return my_data.maxf1_.epsd;      }
    double&   epst()              { return my_data.maxf1_.epst;      }
    double&   yota()              { return my_data.maxf1_.yota;      }

    // Maximization/Computation options
    int&      ihit()              { return my_data.maxf1_.ihit;      }
    int&      ixvc()              { return my_data.maxf1_.ixvc;      }
    int&      maxit()             { return my_data.maxf1_.maxit;     }
    int&      method()            { return my_data.maxf1_.method;    }

    // Parameter options
    int&      nt()                { return my_data.maxf1_.nt;        }
    int&      istin(int i)        { return my_data.maxf1_.istin[i];  }
    double&   stpin(int i)        { return my_data.maxf1_.stpin[i];  }
    double&   thin(int i)         { return my_data.maxf1_.thin[i];   }
    double&   thl(int i)          { return my_data.maxf1_.thl[i];    }
    double&   thu(int i)          { return my_data.maxf1_.thu[i];    }

    // Optional Input

    // Read the i'th label
    const string& label(int i)         const  { return my_data.maxflb_.label[i]; }

    // Set the i'th label
    void set_label(int i, const string& l) 
    {
      my_data.maxflb_.label[i] = l;
    }

    // Interface to parameter options
    int        param_count()           const  { return my_data.maxf1_.nt;        }
    double     initial_estimate(int i) const  { return my_data.maxf1_.thin[i];   }
    double     lower_bound(int i)      const  { return my_data.maxf1_.thl[i];    }
    double     upper_bound(int i)      const  { return my_data.maxf1_.thu[i];    }

    // Automatic Output
    double     value()                 const  { return f;        }
    double     param(int i)            const  { return theta[i]; }
    int        last_error()            const  { return lfl;      }
    int        evaluations()           const  { return nfe;      }
    
    // Optional Output Information
    double     av(int i, int j)        const  { return my_data.maxf2_.av[j*NP()+i]; }
    double     h(int i, int j)         const  { return my_data.maxf2_.h[j*NP()+i];  }
    double     g(int i)                const  { return my_data.maxf2_.g[i];         } 
    double     score(int i)            const  { return my_data.score_vec[i];        } 

    // introduced for score_test 
    int        need_score(int i)       const  { return my_data.need_score[i];       }
    int        score_deriv_conv(int i) const  { return my_data.score_deriv_conv[i]; }

    double     stde(int i)             const  { return my_data.maxf2_.stde[i];      }
    double     stp(int i)              const  { return my_data.maxf2_.stp[i];       }
    double     thpr(int i)             const  { return my_data.maxf2_.thpr[i];      }
    double     cth(int i)              const  { return my_data.maxf2_.cth[i];       }
    int        ist(int i)              const  { return my_data.maxf2_.ist[i];       }

    double     difmax()                const  { return my_data.maxf2_.difmax;       }
    double     erm()                   const  { return my_data.maxf2_.erm;          }
    double     fch()                   const  { return my_data.maxf2_.fch;          }
    double     fpr()                   const  { return my_data.maxf2_.fpr;          }
    int        idif()                  const  { return my_data.maxf2_.idif;         }
    int        igage()                 const  { return my_data.maxf2_.igage;        }
    int        igfl()                  const  { return my_data.maxf2_.igfl;         }
    int        impbnd()                const  { return my_data.maxf2_.impbnd;       }
    int        it()                    const  { return my_data.maxf2_.it;           }
    long       ivage()                 const  { return my_data.maxf2_.ivage;        }
    int        ivfl()                  const  { return my_data.maxf2_.ivfl;         }
    int        nb()                    const  { return my_data.maxf2_.nb;           }
    int        nd()                    const  { return my_data.maxf2_.nd;           }
    int        ne()                    const  { return my_data.maxf2_.ne;           }
    int        ni()                    const  { return my_data.maxf2_.ni;           }
    int        nsurf2()                const  { return my_data.maxf2_.nsurf2;       }
    int        nv()                    const  { return my_data.maxf2_.nv;           }

    // Status labels
    string     maxfst(int i)           const  { return my_data.maxfst_[i];          }

    const Maxfun_Data& get_results()   const;

    static bool iteration_ended;

  protected:

    // Data members

    // Local storage of function value evaluted at theta
    double f;

    // Error flag storage for results of maximization.  See MAXFUN manual.
    int lfl;

    // Number of times _fun has been evaluated.  Reset only by calling reset()
    int nfe;  

  public: // SGROSS
    // Array of parameter estimates
    vector<double> theta;

  protected: // SGROSS    
    // Pointers to the function and parameter validation functions
    MaxFunction* fun;

  private:

    int maxfun_();

    /* Other function prototypes */
    int    bsrch_(int& ifl);
    int    nsrch_(vector<double>& dth, int& ifl);
    int    psrch_(vector<double>& dth, int& ifl);
    int    comb_(int ntotal, int nchoos, int& nsw, vector<int>& list);
    int    nrstep_(int& ifl);
    int    binit_(int& ihess);
    int    bupdt_(vector<double>& gpr, int& lex);
    int    direct_(int& lex);
    int    lsrch_(int& ifirst, int& ifl);
    int    prepd_(int& lex);
    int    fixbnd_(int& lex);
    int    bndchk_(double& th, double& thb, int& lex);
    int    bcnvch_(int& lex);
    int    vcnvch_(int& lex);
    int    itest_(int& i, double& eloc, int& lex);
    int    endit_();
    int    augv_(int& ih);
    int    deriv1_();
    int    deriv2_(int& ih, int& lex);
    int    fitder_(double& d1, double& d2, double& d3, double& dd, double& prmu, int& lex);
    double dfn_(double ath, double sf);
    double dfninv_(double ath, double delth);
    double efn_(double ath, double eloc);
    int    mmult_(const vector<double>& a, int mra, int ma, int namb,
                  const vector<double>& b, int mrb, int nb, vector<double>& c, int mrc);

    bool   score_deriv1();
    bool   score_deriv2();
    double calc_score();
    bool   score_mxnvrt_(vector<double>& a, int mra, int m, vector<double>& b, int mrb, int& lex);

    int    vcmx_deriv2_(int& ih, int& lex);
    int    augv_fitder_(double& d1, double& d2, double& d3, double& dd, double& prmu, int& lex);
    int    vcmx_fitder_(double& d1, double& d2, double& d3, double& d4, double& dd, double& prmu, int& lex);

    int    vcmx_(int& ih);
    int    mxnvrt_(vector<double>& a, int mra, int m, vector<double>& b, int mrb, int& lex);

    //int    mout_(vector<double>& a, int *mra, int *m, int *n, int *ipr);
    //int    lbd_();
    //int    lbv_(int& lin);
    //int    litd_();
    //int    litv_();
    //int    lf_();

    mutable Maxfun_Data              my_data;


    // static crud from comb_
    vector<int>                      maxlst;
    int                              kv;
};

} // End of namespace SAGE

MEM_COUNT_BEGIN(SAGE::Maxfun_Data::maxf1_struct)
{
  return get_mem(t.thin) + get_mem(t.thl) + get_mem(t.thu) + get_mem(t.stpin) + get_mem(t.istin) + sizeof(t);
}
MEM_COUNT_END

MEM_COUNT_BEGIN(SAGE::Maxfun_Data::maxf2_struct)
{
  return get_mem(t.thpr) + get_mem(t.cth) + get_mem(t.stp) + get_mem(t.g) + get_mem(t.h) + get_mem(t.v) + get_mem(t.av) + get_mem(t.stde) + get_mem(t.ult) + get_mem(t.diag) + get_mem(t.pdir) + get_mem(t.ist) + sizeof(t);
}
MEM_COUNT_END

MEM_COUNT_BEGIN(SAGE::Maxfun_Data::maxflb_struct)
{
  return get_mem(t.label) + sizeof(t);
}
MEM_COUNT_END


MEM_COUNT_BEGIN(SAGE::Maxfun_Data)
{
//  return get_mem(t.param_NP) + get_mem(t.maxf1_) + get_mem(t.maxf2_) + get_mem(t.maxflb_) + sizeof(t);
// appended below to take into account new vectors
  return get_mem(t.param_NP) + get_mem(t.maxf1_) + get_mem(t.maxf2_) + get_mem(t.maxflb_)  
         + get_mem(t.param_NP) + get_mem(t.param_NP) + get_mem(t.param_NP) + get_mem(t.param_NP) 
         + get_mem(t.param_NP) + get_mem(t.param_NPV2) + get_mem(t.param_NP) + sizeof(t);
}
MEM_COUNT_END



#endif

