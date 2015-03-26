#ifndef DECIPHER_SHUFFLER_H
#define DECIPHER_SHUFFLER_H

//============================================================================
// File:      sub_pop_shuffler.h
//                                                                          
// Author:    Yeunjoo Song
//                                                                          
// History:   Initial implementation.                                Feb. 05
//                                                                          
// Notes:     This header file contains a class for generating and storing
//            shuffled sub_populations information for permutation test.
//                                                                          
// Copyright (c) 2004 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "decipher/em.h"

namespace SAGE
{

namespace DECIPHER
{

class sub_pop_shuffler
{
  public:

    typedef vector<vector<member> >  sub_pops;
  
    // Constructor/destructor.
    sub_pop_shuffler();

    sub_pop_shuffler(const vector<const member_em_phenotype_map*>& sub_pops);

    void set_sub_pops(const vector<const member_em_phenotype_map*>& sub_pops);

    bool do_shuffle(int seed = 0);

    const vector<member>& get_whole_pop()    const;
    const sub_pops&       get_new_sub_pops() const;

    void dump(ostream& out) const;    
    
  private:

    void pool_members(const set<member, member_order<member> >& mems,  vector<member>& new_mems);
    void get_new_members(vector<member>& new_sub, size_t start, size_t count);

    // Data members.
    //
    struct MTRandomizer
    {
      unsigned int operator()(unsigned int N)
      {
        /*
        double tmp = static_cast<double>(rand())
                   / static_cast<double>(RAND_MAX);

        unsigned int r = static_cast<unsigned int>(tmp*N);

        double mt_tmp1 = static_cast<double>(mt.uniform_integer())
                      / static_cast<double>(RAND_MAX);

        unsigned int mt_r1 = static_cast<unsigned int>(mt_tmp1*N);
        */
        double mt_tmp2 = mt.uniform_real();

        unsigned int mt_r2 = static_cast<unsigned int>(mt_tmp2*N);
        
        //cout << N << " " << tmp << " " << r;
        //cout << "	" << mt_tmp1 << " " << mt_r1;
        //cout << "	" << mt_tmp2 << " " << mt_r2 << endl;

        return mt_r2;
      }

      MersenneTwister mt;
    };

    bool            my_valid_data;

    MTRandomizer    my_randomizer;

    vector<member>  my_whole_pop;
    sub_pops        my_new_sub_pops;
};


#include "decipher/sub_pop_shuffler.ipp"

}
} 

#endif

