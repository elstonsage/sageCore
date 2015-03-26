//============================================================================
//  File:       parental_penetrance_internals.h
//
//  Author:     Kai He
//
//  History:    Version 1.00.00
//
//  Notes:      for interanl definitions that should only be used for testing
//
//  Copyright (c) 2003 R.C. Elston
//  All Rights Reserved        
//============================================================================

#include "geno_elim.new/parental_penetrance.h"

namespace SAGE
{

//-----------------------------------------------------------
//
// Class:       int_pen_model_cursor       
//
// Is: an internal cursor from normal penetrance
//
//-----------------------------------------------------------

class int_pen_model_cursor : public int_pen_cursor
{
  public:

    friend class parental_penetrance_cursor_storage;  

    typedef int_pen_cursor  this_type;

    typedef MLOCUS::penetrance_model::phased_penetrance_iterator phased_penetrance_iterator;
    
    typedef parental_penetrance_cursor_storage cursor_storage;

    int_pen_model_cursor& operator=(const int_pen_model_cursor& that)
    {
       if(this!=&that)
       {
          my_location = that.my_location;
          
          int_pen_cursor::operator=(that);
       }
       return *this;
    }

    int_pen_model_cursor(){}
    
    int_pen_model_cursor(phased_penetrance_iterator ppi)
    {
       my_location = ppi;
    }
    void setup(const parental_penetrance_model* ppm, uint j)
    {
       int_pen_cursor::setup(ppm,j);
       
       my_location = pen_model().phased_penetrance_begin(j + 1);
    } 

    bool  is_begin() const
    {
       if(my_location == pen_model().phased_penetrance_begin(my_id + 1) )

                return true;

       else     return false;
    }

    bool  is_end()   const
    {
       if(my_location == pen_model().phased_penetrance_end(my_id + 1) )

                return true;

       else     return false;

    }

    this_type* copy_me() const;
    void       release_me();

    void             operator++()       { ++my_location;}
    uint             geno_id() const    { return my_location.geno_id();}
    MLOCUS::phased_genotype  phased_geno() const{ return my_location.phased_geno();}

    double     penetrance() const;
    double     penetrance(uint mo_gid, uint fa_gid) const;

  //protected:
  //private:
  
    phased_penetrance_iterator my_location;

};
//-----------------------------------------------------------
//
// Class:       int_par_pen_cursor       
//
// Is: a second level internal cursor base class
//
//-----------------------------------------------------------

class int_par_pen_cursor : public int_pen_cursor
{
  public:

    friend class parental_penetrance__cursor_storage;    

    typedef int_pen_cursor                     this_type;

    typedef parental_penetrance_cursor_storage cursor_storage;
          
    int_par_pen_cursor(const int_par_pen_cursor&);
    int_par_pen_cursor() : int_pen_cursor() {}
  
    int_par_pen_cursor(uint mi)
    {
       my_member = mi;
    }

    int_par_pen_cursor& operator=(const int_par_pen_cursor&);


    void setup(const parental_penetrance_model* ppm, uint j)
    {
       int_pen_cursor::setup(ppm,j);
       my_member = j;
    }

    virtual bool is_begin   () const = 0;
    virtual bool is_end     () const = 0;
    virtual void operator++ ()       = 0;
    virtual uint geno_id    () const = 0;

    virtual MLOCUS::phased_genotype  phased_geno() const=0;

    virtual double penetrance() const=0;
    virtual double penetrance(uint mo_gid,uint fa_gid) const=0;
  
    virtual this_type*  copy_me() const = 0;
    virtual void        release_me()    = 0;

  //protected:
  //private:
  
    uint  my_member;

};
//-----------------------------------------------------------
//
// Class:       int_par_cursor       
//
// Is: an internal cursor
//
//-----------------------------------------------------------

class int_par_cursor : public int_par_pen_cursor
{
  public:

    friend class parental_penetrance_cursor_storage;    

    typedef int_pen_cursor                     this_type;
    typedef parental_penetrance_cursor_storage cursor_storage;
          
    int_par_cursor& operator=(const int_par_cursor& that)
    {    
       if(this!=&that)
       {
          my_geno = that.my_geno;

          int_pen_cursor::operator=(that);
       }
       return *this;
    }
    int_par_cursor():int_par_pen_cursor(){}
    
    int_par_cursor(MLOCUS::phased_genotype_iterator pgi)
    {
       my_geno = pgi;
    }
    virtual ~int_par_cursor() { }
  
    void setup(const parental_penetrance_model* ppm, uint j)
    {
       int_par_pen_cursor::setup(ppm,j);

       my_geno = geno_model().phased_genotype_begin();
    }

    bool  is_begin() const
    {
       if(my_geno == geno_model().phased_genotype_begin())
         return true; 
       else
         return false;
    }

    bool  is_end()   const
    { 
       if(my_geno == geno_model().phased_genotype_end())
         return true; 
       else
         return false;
    }

    void operator++()       { ++my_geno;   }

    uint geno_id() const    { return my_geno->get_id(); }

    MLOCUS::phased_genotype   phased_geno() const{ return *my_geno; }

    this_type* copy_me() const;
    void release_me();

    double     penetrance() const;
    double     penetrance(uint mo_gid, uint fa_gid) const;
  
  //protected:
  //private:
    
    MLOCUS::phased_genotype_iterator  my_geno;

};
//-----------------------------------------------------------
//
// Class:       int_par_geno_cursor       
//
// Is: an internal cursor
//
//-----------------------------------------------------------

class int_par_geno_cursor : public int_par_pen_cursor
{
  public:

    friend class parental_penetrance_cursor_storage;    

    typedef int_pen_cursor                      this_type;
    typedef parental_penetrance_cursor_storage  cursor_storage;
          
    int_par_geno_cursor& operator=(const int_par_geno_cursor& that)
    {
      if(this!=&that)
      {
        int_pen_cursor::operator=(that);

        count         = that.count;
        mo_gid        = that.mo_gid;
        fa_gid        = that.fa_gid;

        valid_geno[0] = that.valid_geno[0];
        valid_geno[1] = that.valid_geno[1];
        valid_geno[2] = that.valid_geno[2];
        valid_geno[3] = that.valid_geno[3];
        valid_geno[4] = that.valid_geno[4];
      }
      return *this;
    }

    int_par_geno_cursor()
      : int_par_pen_cursor(),
        count(0)
    {
         valid_geno[0] = 
         valid_geno[1] = 
         valid_geno[2] = 
         valid_geno[3] = 
         valid_geno[4] = MLOCUS::NPOS;
    }
    
    int_par_geno_cursor(uint vg[5], uint c);

    void setup(const parental_penetrance_model* ppm,
               MLOCUS::phased_genotype moth,
               MLOCUS::phased_genotype fath,
               const member_type& id)
    {
       int_pen_cursor::setup(ppm,id.subindex());

       mo_gid = moth.get_id();
       fa_gid = fath.get_id();

       count = 0;

       MLOCUS::child_genotype_set cgs(moth,fath,true);
       
       uint vg = 0;

       for(uint i = 0; i < cgs.size(); ++i)
       {
         valid_geno[vg] = cgs[i].get_id();
         ++vg;
       }

       for(; vg<5;++vg)
         valid_geno[vg] = MLOCUS::NPOS;
    }

    bool is_begin() const   
    { 
      if(valid_geno[count] == 0) return true;
      else                       return false;
    }
    bool is_end()   const   
    { 
      if(valid_geno[count]== MLOCUS::NPOS) return true;
      else                         return false;
    }

    void operator++() { ++count; }

    uint geno_id() const { return valid_geno[count]; }

    MLOCUS::phased_genotype  phased_geno() const
    { 
        return geno_model().get_phased_genotype(geno_id()); 
    }
    
    this_type* copy_me() const;
    void release_me();

    double    penetrance() const;
    double    penetrance(uint mo_gid, uint fa_gid) const;
    double    penetrance(const MLOCUS::phased_genotype mo_gid, const MLOCUS::phased_genotype fa_gid) const;
    double    penetrance(const this_type& mo_gid, const this_type& fa_gid) const;
  
  //protected:
  //private:
  
    uint  valid_geno[5];
    uint  count;
    uint  mo_gid;
    uint  fa_gid;

};    

//-----------------------------------------------------------
//
// Class:       parental_penetrance_cursor_storage       
//
// Is: a internal cursor pointer storage with static methods
//
//-----------------------------------------------------------

class parental_penetrance_cursor_storage
{
  public:

    typedef std::list<int_pen_model_cursor*>	list_1;//int_pen_model_cursor;
    typedef std::list<int_par_cursor*>		list_2;//int_par_cursor;
    typedef std::list<int_par_geno_cursor*>	list_3;//int_par_geno_cursor;
    

    static int_pen_model_cursor*  get_int_pen_model_cursor()
    {
       if(my_int_pen_model_cursors.size() > 0 )
       {
 	  int_pen_model_cursor* X = my_int_pen_model_cursors.top();
 	  my_int_pen_model_cursors.pop();
 	  return X;
       }
       else
       {
          return new int_pen_model_cursor();
       }
    }

    static int_par_cursor*  get_int_par_cursor()
    {
       if(my_int_par_cursors.size() > 0 )
       {
          int_par_cursor* X = my_int_par_cursors.top();
          my_int_par_cursors.pop();
          return X;
       }
       else
       {
          return new int_par_cursor();
       }
    }
    
    static int_par_geno_cursor*  get_int_par_geno_cursor()
    {
       if(my_int_par_geno_cursors.size() > 0 )
       {
          int_par_geno_cursor* X = my_int_par_geno_cursors.top();
          my_int_par_geno_cursors.pop();
          return X;
       }
       else
       {
          return new int_par_geno_cursor();
       }
    }

  //private:
  //protected:
    
    static void put_int_pen_model_cursor(int_pen_model_cursor* ipmc) 
    { 
       //list_1 lst1;
       //lst1.push_back(ipmc);
       my_int_pen_model_cursors.push(ipmc); 
    }

    static void put_int_par_cursor(int_par_cursor* ipc)       
    { 
       //list_2 lst2;
       //lst2.push_back(ipc);
       my_int_par_cursors.push(ipc);       
    }

    static void put_int_par_geno_cursor(int_par_geno_cursor* ipgc)  
    { 
       //list_3 lst3;
       //lst3.push_back(ipgc);
       my_int_par_geno_cursors.push(ipgc);  
    }

    parental_penetrance_cursor_storage(){}
            
    static std::stack <int_pen_model_cursor*,list_1 >  my_int_pen_model_cursors;
    static std::stack <int_par_cursor*, list_2 >	       my_int_par_cursors;
    static std::stack <int_par_geno_cursor*,list_3 >   my_int_par_geno_cursors;
};

}
