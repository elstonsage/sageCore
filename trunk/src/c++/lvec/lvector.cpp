#include <iostream>
#include <cmath>
#include <iomanip>
#include "numerics/kahan.h"
#include "lvec/codom_ivgen.h"
#include "lvec/iv_generator.h"
#include "lvec/lvector.h"

namespace SAGE
{

class lv_acceptor : public iv_acceptor
{
public:

  lv_acceptor(Likelihood_Vector& lv) : l(lv) { }

  virtual void accept(equivalence_class e, equivalence_class d, double p)
  { 
// This is useful for testing, sometimes
#if 0
    for(int i = 1; i <= e; i*=2)
      cout << (int)(bool) (e & i);
    cout << " = " << p << endl; 
#endif
    l.increment_value(e, d, p);
  }

protected:

  Likelihood_Vector& l;
};

Likelihood_Vector::Likelihood_Vector(mmap* mm, const SAGE::MLOCUS::inheritance_model& pm)
    : first(0), f(false), bits(0), max_bits(0), my_scale(1.0)
{
  if(!mm) 
  {
    set_valid(false);
    return;
  }

  set_valid(true);

  resize_and_clear_bits(mm->nonfounder_meiosis_count());

  fixed = (1 << bits)-1;

  boost::shared_ptr<lv_acceptor> lva(new lv_acceptor(*this));

  if(pm.codominant(false))
  {
    codominant_iv_generator ivg(lva);

    if(!ivg.build(*mm, pm))
    {
      set_valid(false);
    }
  }
  else
  {
    iv_generator ivg(lva);

    if(!ivg.build(*mm, pm))
    {
      set_valid(false);
    }
  }
}

/*
Likelihood_Vector::Likelihood_Vector
    (const Likelihood_Vector& l, const Marker_Recombination& m)
  : first(l.first), f(l.f), fixed(l.fixed), bits(0), max_bits(0), my_scale(1.0)
{
  if(!l.is_valid())
  {
    set_valid(false);
    return;
  }

  set_valid(true);

  resize(l.bits);

  // If the vector is empty, we don't need to do anything more.
  if(!f)
  {
    clear_bits();

    return;
  }

  bits = max_bits = l.bits;

  if(fixed_bits())
  {
    for(equivalence_class i = start(); i < size(); bump_up(i))
      storage[i] = l.storage[i];
  }
  else
  {
    ::memcpy(&*storage.begin(), &*l.storage.begin(), size() * sizeof(likelihood));
  }

  *this *= m;
}
*/

Likelihood_Vector::Likelihood_Vector
    (const Likelihood_Vector& l, const mmap* mm, double theta)
  : first(l.first), fixed(l.fixed), f(l.f), my_scale(1.0)
{
  if(!l.is_valid() || !mm)
  {
    set_valid(false);
    return;
  }
  set_valid(true);

  resize(l.bits);

  // If the vector is empty, we don't need to do anything more.
  if(!f)
  {
    clear_bits();

    return;
  }

  if(fixed_bits())
  {
    for(equivalence_class i = start(); i < size(); bump_up(i))
      storage[i] = l.storage[i];
  }
  else
  {
    ::memcpy(&*storage.begin(), &*l.storage.begin(), size() * sizeof(likelihood));
  }
  (*this)(mm, theta);
}

Likelihood_Vector& Likelihood_Vector::operator()
    (mmap* mm, const SAGE::MLOCUS::inheritance_model& pm, bool use_pf)
{
  first = 0;
  f     = false;

  if(!mm)
  {
    set_valid(false);
    return *this;
  }

  set_valid(true);

  resize_and_clear_bits(mm->nonfounder_meiosis_count());

  fixed = (1 << bits) - 1;
  
  boost::shared_ptr<lv_acceptor> lva(new lv_acceptor(*this));
  
  bool codom = pm.codominant(false);
  
  if(codom)
  {
    codominant_iv_generator ivg(lva);

    if(!ivg.build(*mm, pm, use_pf))
    {
      set_valid(false);
    }
  }
  else
  {
    iv_generator ivg(lva);


    if(!ivg.build(*mm, pm, use_pf))
    {
      set_valid(false);
    }
  }

  return *this;
}

/*
Likelihood_Vector& Likelihood_Vector::operator*=(const Marker_Recombination& m)
{
  // If the vector is empty, we don't need to do anything more.
  if(!f) return *this;

  const mmap* mm = m.meiosis_map();

  // If we don't have a good mmap, we can't continue
  if(!is_valid() || !mm || !mm->is_valid())
  {
    setstate(failbit);
    return *this;
  }

  // Preprocess:

  bool non_zero = false;

  vector<double> th_p(mm->meiosis_count()); // Theta primes

  double th_tp = 1;                       // theta tilde'

  for(int i = 0; i < mm->founder_meiosis_count(); ++i)
  {
    likelihood d  = m.meiosis(i);

    if(d != 0.0) non_zero = true;

    // If d is 1, we have an error.
    if(d == 1.0)
    {
      setstate(failbit);
      return *this;
    }

    th_p[i] = d / (1 - d); 

    th_tp *= (1 - d);
  }

  for(int i = 0; i < mm->nonfounder_meiosis_count(); ++i)
  {
    likelihood d  = m.meiosis(i + mm->meiosis_bits);

    if(d != 0.0) non_zero = true;

    // If d is 1, we have an error.
    if(d == 1.0)
    {
      setstate(failbit);
      return *this;
    }

    th_p[i + mm->founder_meiosis_count()] = d / (1 - d);

    th_tp *= (1 - d);
  }

  // If all our theta's are zero, we're done.
  if(!non_zero)
    return *this;

  Process(mm, th_p.begin(), th_tp);

  return *this;
}
*/

Likelihood_Vector& Likelihood_Vector::operator() (const mmap* mm, double theta)
{
#if 0
  cout << "Likelihood_Vector::operator().." << theta << endl;
#endif

  // If the vector is empty, we don't need to do anything more.
  if(!f) return *this;

  // If we don't have a good mmap or theta is 1, we have a problem.
  if(!is_valid() || !mm)
  {
    set_valid(false);
    return *this;
  }

  // If theta is 0, the vector remains the same.
  if(theta == 0.0) return *this;

  // Preprocess:

  vector<double> th_p(mm->meiosis_count()); // Theta primes

  double th_tp = 1;                       // theta tilde'

  double thp = theta / (1 - theta), tht = 1 - theta;

#if 0
  cout << "1. th_tp = " << th_tp << ", thp = " << thp << ", tht = " << tht << endl;
  cout << "th_p :" << endl;
  for( size_t i = 0; i < th_p.size(); ++i )
  {
    cout << "th_p[" << i << "] = " << th_p[i] << endl;
  }
#endif

  for(size_t i = 0; i < mm->founder_meiosis_count(); ++i)
  {
    th_p[i] = thp; 

    th_tp *= tht;
  }

#if 0
  cout << "2. th_tp = " << th_tp << ", thp = " << thp << endl;
  cout << "th_p :" << endl;
  for( size_t i = 0; i < th_p.size(); ++i )
  {
    cout << "th_p[" << i << "] = " << th_p[i] << endl;
  }
#endif

  for(size_t i = 0; i < mm->nonfounder_meiosis_count(); ++i)
  {
    th_p[i + mm->founder_meiosis_count()] = thp;

    th_tp *= tht;
  }

#if 0
  cout << "3. th_tp = " << th_tp << ", thp = " << thp << endl;
  cout << "th_p :" << endl;
  for( size_t i = 0; i < th_p.size(); ++i )
  {
    cout << "th_p[" << i << "] = " << th_p[i] << endl;
  }
#endif

  Process(mm, th_p.begin(), th_tp);

#if 0
  cout << "4. th_tp = " << th_tp << ", thp = " << thp << endl;
  cout << "th_p :" << endl;
  for( size_t i = 0; i < th_p.size(); ++i )
  {
    cout << "th_p[" << i << "] = " << th_p[i] << endl;
  }
#endif

  return *this;
}

Likelihood_Vector& Likelihood_Vector::operator() (const fft& f, double theta)
{
  if(f.get_meiosis_map()->nonfounder_meiosis_count() != bit_count() ||
     theta < 0.0 || theta >= 1.0)
  {
    set_valid(false);

    return *this;
  }

  if(theta == 0.0) return *this;

  // Do algorithm as given in Kruglyak and Lander, 1998

  double temp;
  size_t s = size();  // Moved to make Intel Compiler happy
  
  for(equivalence_class i = 1; i < size(); i *= 2)
  {
    for(equivalence_class j = 0; j < size(); ++j)
    {
      equivalence_class jp = j ^ i;
      
      if(j < jp)
      {
        temp        = storage[j];
        storage[j] += storage[jp];
        storage[jp] = temp - storage[jp];
      }
    }
  }

  vector<double> t(f.get_meiosis_map()->meiosis_count()+1);
  
  t[0] = 1.0;
  
  for(size_t i = 1; i < f.get_meiosis_map()->meiosis_count()+1; ++i)
    t[i] = t[i-1] * (1 - 2 * theta);
  
  
  for(equivalence_class j = 0; j < s; ++j)
    storage[j] *= t[f[j]];
  
  for(equivalence_class i = 1; i < size(); i *= 2)
  {
    for(equivalence_class j = 0; j < size(); ++j)
    {
      equivalence_class jp = j ^ i;
      
      if(j < jp)
      {
        temp        = storage[j];
        storage[j] += storage[jp];
        storage[jp] = temp - storage[jp];
      }
    }
  }

  for(equivalence_class j = 0; j < size(); ++j)
    storage[j] /= size();

  first = 0;
  fixed = 0;
  
  return *this;
}


Likelihood_Vector& Likelihood_Vector::operator*= (const lvector& l)
{
  if(!is_valid() || !l.is_valid() || size() != l.size())
  {
    set_valid(false);
    return *this;
  }

  if(!f) return *this;

  if(!l.f)
  {
    *this = l;

    return *this;
  }

  likelihood scaling_factor = total() * l.total() / size() / l.size();

//cout << "my_scale = " << my_scale << ", l.log_scale = " << l.log_scale() << endl;

  my_scale *= log_double(l.log_scale().get_double() * scaling_factor);

//cout << "my_scale = " << my_scale << endl;

  scaling_factor = ((likelihood) 1.0) / scaling_factor;

  equivalence_class i = start(), j = l.start();

  equivalence_class new_first = 0;
  bool new_f = false;

  // kbj: added j<=l.size() check, since j+1<=l.size() should never happen
  while(i < size() && j <= l.size())
  {
    if(i < j)
    {
      while(i < j) { storage[i] = 0.0; bump_up(i); }
    }
    else if(i > j)
    {
      while(i > j) l.bump_up(j);
    }
    else
    {
      if(!new_f) { new_first = i; new_f = true; }

//      storage[i] *= l.storage[j];
      storage[i] *= l.storage[j] * scaling_factor;

      bump_up(i);
      l.bump_up(j);
    }
  }

  fixed = fixed | l.fixed;
  first = new_first;

  f = new_f;

  return *this;
}

Likelihood_Vector& Likelihood_Vector::operator+= (const lvector& l)
{
  if(!is_valid() || !l.is_valid() || size() != l.size())
  {
    set_valid(false);
    return *this;
  }

  // If l is 0.0, then we don't do anything.
  if(!l.f)
  {
    return *this;
  }

  double delta = (l.my_scale / my_scale).get_double();

  equivalence_class j = l.start();

  while(j < l.size())
  {
    increment_value(j, 0, l.storage[j] * delta);

    l.bump_up(j);
  }

  first = min(start(), l.start());

  return *this;
}

void Likelihood_Vector::normalize()
{
  if(!is_valid() || !f) return;

  likelihood l = total();

  if(l == 0) return;

  for(equivalence_class i = start(); i < size(); bump_up(i))
    storage[i] /= l;

  my_scale = 1.0;
}

Likelihood_Vector::likelihood Likelihood_Vector::total() const
{
  // GCW - Add a storage of total so that we only recalcualte it when
  // changed.

  if(!f || !is_valid()) return 0.0;

  KahanAdder<likelihood> l = 0.0;

  for(equivalence_class i = start(); i < size(); bump_up(i))
  {
    //cout << "eq " << i << " = " << storage[i] << ", ";
    l += storage[i];
  }
  //cout << "total = " << l << endl;

  return l;
}

double Likelihood_Vector::information() const
{
  if(!f || !is_valid() ) return 0.0;

  likelihood tot = total();

  if(tot == 0.0) return 0.0;

  likelihood l = 0.0;

  for(equivalence_class i = start(); i < size(); bump_up(i))
  {
    double d = storage[i]/tot;

    if(d == 0.0) continue;

    l += d * log(d) / log((double) 2.0);
  }

  l /= bit_count();

  l += 1;

  return l;
}

void Likelihood_Vector::Process
    (const mmap* mm, vector<double>::iterator th_p, double th_tp)
{
#if 0
  cout << "Likelihood_Vector::Process().."
       << th_tp << ", first = " << first << ", fixed = " << fixed << endl;
#endif

  if(!mm || !is_valid()) return;

  size_t s = size();

  if(fixed_bits())
  {
    for(equivalence_class i = start(); i < s; bump_up(i))
      storage[i] *= th_tp;

    Process_Fixed  (0, bit_count(), th_p + mm->founder_meiosis_count(), mm);
  }
  else
  {
    for(equivalence_class i = 0; i < s; ++i)
      storage[i] *= th_tp;

    Process_Natural(0, bit_count(), th_p + mm->founder_meiosis_count(), mm);
  }

  Process_Founders(mm, th_p);

  // Post Process:

  first = 0;
  fixed = 0;
}

void Likelihood_Vector::Process_Fixed
    (equivalence_class eq, size_t bit, vector<double>::iterator th, const mmap* mm)
{
#if 0
  cout << "Likelihood_Vector::Process_Fixed().."
       << eq << ", " << bit << ", " << *th
       << ", first = " << first << ", fixed = " << fixed << endl;
#endif

  if((int) (--bit) == -1) return;

  equivalence_class bit_p = 1 << bit;  // bit_p == bit'

  if(bit_p & fixed)
  {
    if(bit_p & first)
    {
      Process_Fixed(eq + bit_p, bit, th, mm);

      if( mm->is_x_linked() && mm->is_father_bit(bit_p) )
        ;//cout << "x_linked 1 " << bit << ", " << bit_p << endl;
      else
        for(equivalence_class j = eq; j < eq + bit_p; ++j)
          storage[j] = th[bit] * storage[j + bit_p];
    }
    else
    {
      Process_Fixed(eq, bit, th, mm);

      if( mm->is_x_linked() && mm->is_father_bit(bit_p) )
        ;//cout << "x_linked 2 " << bit << ", " << bit_p << endl;
      else
        for(equivalence_class j = eq; j < eq + bit_p; ++j)
        {
          storage[j + bit_p] = th[bit] * storage[j];
        }
    }
    return;
  }

  Process_Fixed(eq,         bit, th, mm);
  Process_Fixed(eq + bit_p, bit, th, mm);

  if( mm->is_x_linked() && mm->is_father_bit(bit_p) )
  {
    //cout << "x_linked 3 " << bit << ", " << bit_p << endl;
    return;
  }

  likelihood temp;
  for(equivalence_class j = eq; j < eq + bit_p; ++j)
  {
    temp = storage[j];
 
    storage[j]         += th[bit] * storage[j + bit_p];
    storage[j + bit_p] += th[bit] * temp;
  }

  return;
}

void Likelihood_Vector::Process_Natural
    (equivalence_class eq, size_t bit, vector<double>::iterator th, const mmap* mm)
{
#if 0
  cout << "Likelihood_Vector::Process_Natural().." << endl;
#endif

  if(((int) --bit) == -1) return;

  equivalence_class bit_p = 1 << bit;  // bit_p == bit'

  Process_Natural(eq,         bit, th, mm);
  Process_Natural(eq + bit_p, bit, th, mm);

  if( mm->is_x_linked() && mm->is_father_bit(bit_p) )
    return;

  likelihood temp;
  for(equivalence_class j = eq; j < eq + bit_p; ++j)
  {
    temp = storage[j];

    storage[j]         += th[bit] * storage[j + bit_p];
    storage[j + bit_p] += th[bit] * temp;
  }

  return;
}

void Likelihood_Vector::Process_Founders
    (const mmap* mm, vector<double>::iterator th)
{
#if 0
  cout << "Likelihood_Vector::Process_Founders().." << *th << endl;
#endif

  if(!mm) return;

  likelihood temp;

  for(size_t i = 0; i < mm->founder_meiosis_count(); ++i)
  {
    equivalence_class mask = mm->mask(i);

    if( mm->is_x_linked() && mm->is_father_bit(mask) )
      continue;

    equivalence_class j_p;  //j_p == j'

    if(mask)
      for(equivalence_class j = 0; j < size(); ++j)
      {
        j_p = j ^ mask;

        if(j_p < j)
        {
          temp = storage[j];

          storage[j]   += th[i] * storage[j_p];
          storage[j_p] += th[i] * temp;
        }
      }
    else
      for(equivalence_class j = 0; j < size(); ++j)
        storage[j] *= (1 + th[i]);
  }

  temp = 0.0;

  return;
}

}
