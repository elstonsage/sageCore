#include <cstring>
#include "maxfun/maxfun.h"

using namespace std;
using namespace SAGE;

MaxFunction::~MaxFunction()
{}

Maxfun_Data::Maxfun_Data()
{
  param_NP = param_NPV = 100;
  param_NPV2 = param_NPV*param_NPV; // due to JA

  maxf1_.nt = 0;

  maxf1_.epsc1 = 0.001;
  maxf1_.epsc2 = 1e-15;
  maxf1_.epsc3 = 0.0;

  maxf1_.epsd  = 0.0;
  maxf1_.epst  = 0.0;
  maxf1_.yota  = 0.0;
  
  maxf1_.ihit  = 0;
  maxf1_.ixvc  = 2;
  maxf1_.maxit = 0;

  maxf1_.thin.resize (param_NP, 0.0);
  maxf1_.thl.resize  (param_NP,-std::numeric_limits<double>::infinity());
  maxf1_.thu.resize  (param_NP, std::numeric_limits<double>::infinity()); 
  maxf1_.stpin.resize(param_NP, 0.1); 
  maxf1_.istin.resize(param_NP, 1);

  maxf2_.thpr.resize (param_NP);
  maxf2_.cth.resize  (param_NP);
  maxf2_.stp.resize  (param_NP);
  maxf2_.g.resize    (param_NPV);
  maxf2_.h.resize    (param_NPV*param_NPV);
  maxf2_.v.resize    (param_NPV*param_NPV);
  maxf2_.av.resize   (param_NP*param_NP);
  maxf2_.stde.resize (param_NP);
  maxf2_.ult.resize  (param_NPV*param_NPV);
  maxf2_.diag.resize (param_NPV);
  maxf2_.pdir.resize (param_NPV);
  maxf2_.ist.resize  (param_NP);

  maxflb_.label.resize(param_NP);

  // Set all the key variables
  // new vectors due to JA
  theta.resize(NP(), std::numeric_limits<double>::quiet_NaN());
  score_vec.resize(param_NPV, std::numeric_limits<double>::quiet_NaN());
  score_vec2.resize(param_NPV2, std::numeric_limits<double>::quiet_NaN());
  score_stepsize.resize(param_NPV, std::numeric_limits<double>::quiet_NaN());
  need_score.resize(param_NPV,0);
  score_deriv_conv.resize(param_NPV,9);
  flat_dir.resize(param_NP,""); // due to JA for flagging flat directions
  
  f   = std::numeric_limits<double>::quiet_NaN();
  lfl = 0;
  nfe = 0;
}

Maxfun_Data::Maxfun_Data(const Maxfun_Data& md)
  : f        (md.f),
    lfl      (md.lfl),
    nfe      (md.nfe),
    theta    (md.theta),
    score_vec (md.score_vec),
    score_vec2 (md.score_vec2),
    need_score (md.need_score),
    score_deriv_conv (md.score_deriv_conv), 
    param_NP (md.param_NP),
    param_NPV(md.param_NPV),
    maxf1_   (md.maxf1_),
    maxf2_   (md.maxf2_),
    maxflb_  (md.maxflb_)
{}

Maxfun_Data& Maxfun_Data::operator=(const Maxfun_Data& md)
{
  f         = md.f;
  lfl       = md.lfl;
  nfe       = md.nfe;
  theta     = md.theta;
  score_vec = md.score_vec; 
  score_vec2 = md.score_vec2; 
  need_score = md.need_score;
  score_deriv_conv = md.score_deriv_conv;
  param_NP  = md.param_NP;
  param_NPV = md.param_NPV;
  maxf1_    = md.maxf1_;
  maxf2_    = md.maxf2_;
  maxflb_   = md.maxflb_;

  return *this;
}


// Allocate memory for maxfun variables on the heap
void Maxfun::init()
{
  theta.resize(NP());

  maxlst.resize      (my_data.param_NPV);
  kv = 0;
}

const Maxfun_Data& Maxfun::get_results() const
{
  my_data.f         = f;
  my_data.lfl       = lfl;
  my_data.nfe       = nfe;
  my_data.theta     = theta;

  return my_data;
}

void Maxfun_Data::clear_results()
{
  maxf2_.thpr.clear ();
  maxf2_.cth.clear  ();
  maxf2_.stp.clear  ();
  maxf2_.g.clear    ();
  maxf2_.h.clear    ();
  maxf2_.v.clear    ();
  maxf2_.av.clear   ();
  maxf2_.stde.clear ();
  maxf2_.ult.clear  ();
  maxf2_.diag.clear ();
  maxf2_.pdir.clear ();
  maxf2_.ist.clear  ();

  maxf2_.thpr.resize (param_NP,0);
  maxf2_.cth.resize  (param_NP,0);
  maxf2_.stp.resize  (param_NP,0);
  maxf2_.g.resize    (param_NPV,0);
  maxf2_.h.resize    (param_NPV*param_NPV,0);
  maxf2_.v.resize    (param_NPV*param_NPV,0);
  maxf2_.av.resize   (param_NP*param_NP,0);
  maxf2_.stde.resize (param_NP,0);
  maxf2_.ult.resize  (param_NPV*param_NPV,0);
  maxf2_.diag.resize (param_NPV,0);
  maxf2_.pdir.resize (param_NPV,0);
  maxf2_.ist.resize  (param_NP,0);
}
