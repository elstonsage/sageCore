//============================================================================
// File:      relpair.ipp      
//                                                                          
// Author:    Dan Baechle & Kevin Jacobs                                    
//                                                                          
// History:   7/00   created.  - djb
//            4/4/01 modified so that relative_pair class includes connecting 
//                   members.  -djb
//                                                                          
// Notes:     Inline implementation for the following classes -
//              pair_generator 
//                relative pair
//                base_iterator
//                iterator
//                const_iterator
//              ind_filter_trait
//              pair_filter_trait
//              ind_filter                                 
//              pair_filter      
//              filtering_pair_generator
//                base_iterator
//                iterator
//                const_iterator
//                                                                          
// Copyright (c) 2000 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


// - Print a relative pair.
//
inline std::ostream&
operator <<(std::ostream& out, pair_generator::relative_pair pr) 
{
  out << "pair type:  " << std::setw(10) << left << pair_generator::pair_type_to_string(pr.type()) << "      ";
  out << "pedigree:  " << pr.member_one()->pedigree()->name() << "   ";
  out << "member id's:  " << std::setw(6) << std::right << pr.member_one()->name() << "  ";
  out << std::setw(6) << std::right << pr.member_two()->name() << std::endl;
  return out;   
}


//============================================================================
// IMPLEMENTATION:  pair_generator
//============================================================================
//
inline std::string
pair_generator::pair_type_to_string(pair_type type)
{
  switch(type)
  {
    case PARENTAL:
      return "PARENTAL";
    case SIBSIB:
      return "SIBSIB";
    case SISSIS:
      return "SISSIS";
    case BROBRO:
      return "BROBRO";
    case BROSIS:
      return "BROSIS";
    case GRANDP:
      return "GRANDP";
    case AVUNC:
      return "AVUNC";
    case HALFSIB:
      return "HALFSIB";
    case COUSIN:
      return "COUSIN";
    case EVERY:
      return "EVERY";      
    default:
      return "NULL_TYPE";
  }
}

inline pair_generator::pair_type
pair_generator::mask_to_pair_type(mask m)
{
  switch(m)
  {
    case PARENTAL_MASK:
      return PARENTAL;
    case SIBSIB_MASK:
      return SIBSIB;
    case SISSIS_MASK:
      return SISSIS;
    case BROBRO_MASK:
      return BROBRO;
    case BROSIS_MASK:
      return BROSIS;
    case GRANDP_MASK:
      return GRANDP;
    case AVUNC_MASK:
      return AVUNC;
    case HALFSIB_MASK:
      return HALFSIB;
    case COUSIN_MASK:
      return COUSIN;
    case EVERY_MASK:
      return EVERY;      
    default:
      return NULL_TYPE;
  }  
}

inline pair_generator::mask
pair_generator::pair_type_to_mask(pair_type p)
{
  switch(p)
  {
    case PARENTAL:
      return PARENTAL_MASK;
    case SIBSIB:
      return SIBSIB_MASK;
    case SISSIS:
      return SISSIS_MASK;
    case BROBRO:
      return BROBRO_MASK;
    case BROSIS:
      return BROSIS_MASK;
    case GRANDP:
      return GRANDP_MASK;
    case AVUNC:
      return AVUNC_MASK;
    case HALFSIB:
      return HALFSIB_MASK;
    case COUSIN:
      return COUSIN_MASK;
    case EVERY:
      return EVERY_MASK;      
    default:
      return NULL_TYPE_MASK;
  }  
}

inline
pair_generator::pair_generator(unsigned int types)
      : my_p(0), my_types(types)
{
  my_first_type = first_type();
}

inline
pair_generator::pair_generator(RPED::RefPedigree* p, unsigned int types) 
      : my_p(p), my_types(types)
{
  my_first_type = first_type();
}
               
inline pair_generator::iterator
pair_generator::begin()
{
  if(!my_p)
  {
    iterator temp;  // An end iterator.
    return temp;
  }
  else
  {
    iterator temp(this, my_p->family_begin());
    return temp;
  }
}

inline pair_generator::iterator
pair_generator::end()
{
  if(!my_p)
  {
    iterator temp;  // An end iterator.
    return temp;
  }
  else
  {
    pair_generator::iterator temp(this, my_p->family_end());
    return temp;
  }
}

inline pair_generator::const_iterator
pair_generator::begin() const
{
  if(!my_p)
  {
    const_iterator temp;  // An end iterator.
    return temp; 
  }
  else
  {
    const_iterator temp(this, my_p->family_begin());
    return temp;  
  }
}

inline pair_generator::const_iterator
pair_generator::end() const
{
  if(!my_p)
  {
    const_iterator temp;  // An end iterator.
    return temp;
  }
  else
  {
    const_iterator temp(this, my_p->family_end());
    return temp;
  }
}

inline unsigned int      
pair_generator::types() const
{
  return my_types;
}

inline bool
pair_generator::check_type(mask m) const
{
  return (my_types & m);
}

inline pair_generator::pair_type
pair_generator::first_type() const  
{ 
  if(my_types & EVERY_MASK)
  {
    return  EVERY;
  }
                                     
  for(int i = 1; i <= COUSIN; ++i)
  {
    if((my_types >> (i - 1)) & 1U)
    {
      return pair_type(i);
    } 
  }
  
  return pair_type(NULL_TYPE);   
}

inline pair_generator::pair_type
pair_generator::next_type(pair_generator::pair_type current_type) const
{
  for(int i = current_type + 1; i <= COUSIN; ++i)
  {
    if((my_types >> (i - 1)) & 1U)
    {
      return pair_type(i); 
    }
  }
  
  return first_type();
}

inline void 
pair_generator::set_pedigree(RPED::RefPedigree* pedigree)
{
  my_p = pedigree;
}

inline void
pair_generator::set_types(unsigned int types)
{
  my_types = types;
  my_first_type = first_type();
}

inline void
pair_generator::set_type(pair_type p)
{
  my_types = pair_type_to_mask(p);
  my_first_type = first_type();
}

inline void
pair_generator::add_type(mask m)
{
  my_types = my_types | m;
}

inline bool
pair_generator::operator ==(const pair_generator& other) const
{
  return (my_p == other.my_p && my_types == other.my_types);
}

inline bool
pair_generator::operator !=(const pair_generator& other) const
{
  return !(operator ==(other));
}


//============================================================================
// IMPLEMENTATION:  pair_generator::relative_pair      
//============================================================================
//
inline pair_generator::pair_type
pair_generator::relative_pair::type() const
{
  return my_type;
} 

inline pair_generator::relative_pair::member             
pair_generator::relative_pair::member_one() const
{
  return my_member_one;
}

inline pair_generator::relative_pair::member             
pair_generator::relative_pair::member_two() const
{
  return my_member_two;
}

inline pair_generator::relative_pair::member             
pair_generator::relative_pair::connector_one() const
{
  return my_connector_one;
}

inline pair_generator::relative_pair::member             
pair_generator::relative_pair::connector_two() const
{
  return my_connector_two;
}

inline bool
pair_generator::relative_pair::operator ==(const relative_pair& other) const
{
  return       my_member_one     == other.my_member_one   
            && my_member_two     == other.my_member_two
            && my_connector_one  == other.my_connector_one
            && my_connector_two  == other.my_connector_two
            && my_type           == other.my_type;
}

inline bool
pair_generator::relative_pair::operator !=(const relative_pair& other) const
{
  return !(operator ==(other));
}

inline 
pair_generator::relative_pair::relative_pair()
{}

inline 
pair_generator::relative_pair::relative_pair(pair_generator::relative_pair::member member_one, 
                                             pair_generator::relative_pair::member member_two,
                                             pair_generator::relative_pair::member connector_one,
                                             pair_generator::relative_pair::member connector_two,
                                             pair_generator::pair_type type) 
      : my_member_one(member_one), my_member_two(member_two),
        my_connector_one(connector_one), my_connector_two(connector_two), my_type(type) 
{}
  
//============================================================================
// IMPLEMENTATION:  pair_generator::base_iterator
//============================================================================
//
inline
pair_generator::base_iterator::base_iterator(const base_iterator& other)
{
  std::auto_ptr<pair_generator>  ap(new pair_generator(*(other.my_generator)));
  my_generator            = ap;
  my_pair                 = other.my_pair;
  my_family               = other.my_family;
  my_children             = other.my_children;
  my_fullsibs             = other.my_fullsibs;
  my_cousins              = other.my_cousins;
  my_halfsibs             = other.my_halfsibs;
  my_parents_mates        = other.my_parents_mates;
  my_avuncs_mates         = other.my_avuncs_mates;
  my_avuncs               = other.my_avuncs;
  
  my_subped               = other.my_subped;
  my_member_one           = other.my_member_one;
  my_member_two           = other.my_member_two;
  
  my_parent               = other.my_parent;
  my_grandp               = other.my_grandp;
  my_parent_number        = other.my_parent_number;
  my_grandp_number        = other.my_grandp_number;
  my_current_type         = other.my_current_type;
  pair_valid              = other.pair_valid;
  at_end                  = other.at_end;
  halfsibs_init           = other.halfsibs_init;
  cousins_init            = other.cousins_init;
}

inline pair_generator::base_iterator&
pair_generator::base_iterator::operator  =(const base_iterator& other)
{
  if(&other != this)
  {
    std::auto_ptr<pair_generator>  ap(new pair_generator(*(other.my_generator)));
    my_generator            = ap;
    my_pair                 = other.my_pair;
    my_family               = other.my_family;
    my_children             = other.my_children;
    my_fullsibs             = other.my_fullsibs;
    my_cousins              = other.my_cousins;
    my_halfsibs             = other.my_halfsibs;
    my_parents_mates        = other.my_parents_mates;
    my_avuncs_mates         = other.my_avuncs_mates;
    my_avuncs               = other.my_avuncs;
    
    my_subped               = other.my_subped;
    my_member_one           = other.my_member_one;
    my_member_two           = other.my_member_two;
    
    my_parent               = other.my_parent;
    my_grandp               = other.my_grandp;
    my_parent_number        = other.my_parent_number;
    my_grandp_number        = other.my_grandp_number;
    my_current_type         = other.my_current_type;
    pair_valid              = other.pair_valid;
    at_end                  = other.at_end;
    halfsibs_init           = other.halfsibs_init;
    cousins_init            = other.cousins_init;
  }
  
  return *this;
} 

inline
pair_generator::base_iterator::~base_iterator()
{}

inline bool
pair_generator::base_iterator::operator ==(const base_iterator& other) const
{
  if(at_end || other.at_end)
  {
    return at_end && other.at_end;
  }
  else
  {
    return my_pair == other.my_pair;
  }
}

inline bool
pair_generator::base_iterator::operator !=(const base_iterator& other) const
{
  return !operator ==(other);
}

inline
pair_generator::base_iterator::base_iterator()
      : at_end(true)
{}

inline
pair_generator::base_iterator::base_iterator(pair_generator* generator,
                                             RPED::RefPedigree::family_iterator family) 
      : my_family(family)
{
  std::auto_ptr<pair_generator> ap(new pair_generator(*generator));
  my_generator = ap;
  my_subped = my_generator->my_p->subpedigree_begin();
}

inline
pair_generator::base_iterator::base_iterator(const pair_generator* generator,
                                             RPED::RefPedigree::family_iterator family) 
      : my_family(family)
{
  std::auto_ptr<pair_generator> ap(new pair_generator(*generator));
  my_generator = ap;
  my_subped = my_generator->my_p->subpedigree_begin();  
}
    
    
//============================================================================
// IMPLEMENTATION:  pair_generator::iterator
//============================================================================
//  
inline
pair_generator::iterator::iterator(const iterator& other)
      : base_iterator(other)
{}

inline pair_generator::iterator&
pair_generator::iterator::operator  =(const iterator& other)
{
  if(&other != this)
  {
    std::auto_ptr<pair_generator>  ap(new pair_generator(*(other.my_generator)));
    my_generator            = ap;
    my_pair                 = other.my_pair;
    my_family               = other.my_family;
    my_children             = other.my_children;
    my_fullsibs             = other.my_fullsibs;
    my_cousins              = other.my_cousins;
    my_halfsibs             = other.my_halfsibs;
    my_parents_mates        = other.my_parents_mates;
    my_avuncs_mates         = other.my_avuncs_mates;
    my_avuncs               = other.my_avuncs;
    
    my_subped               = other.my_subped;
    my_member_one           = other.my_member_one;
    my_member_two           = other.my_member_two;
    
    my_parent               = other.my_parent;
    my_grandp               = other.my_grandp;
    my_parent_number        = other.my_parent_number;
    my_grandp_number        = other.my_grandp_number;
    my_current_type         = other.my_current_type;
    pair_valid              = other.pair_valid;
    at_end                  = other.at_end;
    halfsibs_init           = other.halfsibs_init;
    cousins_init            = other.cousins_init;
  }
  
  return *this;
}

inline
pair_generator::iterator::~iterator()
{}

inline pair_generator::relative_pair&                    
pair_generator::iterator::operator *() 
{
  return my_pair;
}

inline pair_generator::relative_pair*
pair_generator::iterator::operator ->() 
{
  return &my_pair;
}

inline
pair_generator::iterator::iterator()
{}
                                                   
//============================================================================
// IMPLEMENTATION:  pair_generator::const_iterator
//============================================================================
//  
inline
pair_generator::const_iterator::const_iterator(const const_iterator& other)
      : base_iterator(other)
{}

inline pair_generator::const_iterator&
pair_generator::const_iterator::operator  =(const const_iterator& other)
{
  if(&other != this)
  {
    std::auto_ptr<pair_generator>  ap(new pair_generator(*(other.my_generator)));
    my_generator            = ap;
    my_pair                 = other.my_pair;
    my_family               = other.my_family;
    my_children             = other.my_children;
    my_fullsibs             = other.my_fullsibs;
    my_cousins              = other.my_cousins;
    my_halfsibs             = other.my_halfsibs;
    my_parents_mates        = other.my_parents_mates;
    my_avuncs_mates         = other.my_avuncs_mates;
    my_avuncs               = other.my_avuncs;
    
    my_subped               = other.my_subped;
    my_member_one           = other.my_member_one;
    my_member_two           = other.my_member_two;    
    
    my_parent               = other.my_parent;
    my_grandp               = other.my_grandp;
    my_parent_number        = other.my_parent_number;
    my_grandp_number        = other.my_grandp_number;
    my_current_type         = other.my_current_type;
    pair_valid              = other.pair_valid;
    at_end                  = other.at_end;
    halfsibs_init           = other.halfsibs_init;
    cousins_init            = other.cousins_init;
  }
  
  return *this;
}

inline
pair_generator::const_iterator::~const_iterator()
{}

const inline pair_generator::relative_pair&                    
pair_generator::const_iterator::operator *() 
{
  return my_pair;
}

const inline pair_generator::relative_pair*
pair_generator::const_iterator::operator ->() 
{
  return &my_pair;
}

inline
pair_generator::const_iterator::const_iterator()
{}  

//============================================================================
// IMPLEMENTATION:  ind_filter_trait
//============================================================================
//
inline
ind_filter_trait::ind_filter_trait(size_t t, double min, double max)
      : my_trait(t), my_min(min), my_max(max) 
{
  // - Determines affection status for continuous traits.
  //   Set defaults such that individuals are effectively not filtered
  //   for effection status.
  my_affected_min    =  -std::numeric_limits<double>::infinity();
  my_affected_max    =   std::numeric_limits<double>::infinity();
  my_unaffected_min  =  -std::numeric_limits<double>::infinity();
  my_unaffected_max  =   std::numeric_limits<double>::infinity();
}
    
inline size_t           
ind_filter_trait::get_trait() const
{
  return my_trait;
}
    
inline double
ind_filter_trait::get_min() const
{
  return my_min;
}

inline double 
ind_filter_trait::get_max() const
{
  return my_max;
}

inline double
ind_filter_trait::get_affected_min() const
{
  return my_affected_min;
}

inline double 
ind_filter_trait::get_affected_max() const
{
  return my_affected_max;
}

inline double
ind_filter_trait::get_unaffected_min() const
{
  return my_unaffected_min;
}

inline double 
ind_filter_trait::get_unaffected_max() const
{
  return my_unaffected_max;
}

inline void 
ind_filter_trait::clear()
{
  my_trait           =  (size_t)(-1);
  my_min             = -std::numeric_limits<double>::infinity();
  my_max             =  std::numeric_limits<double>::infinity();
  my_affected_min    = -std::numeric_limits<double>::infinity();
  my_affected_max    =  std::numeric_limits<double>::infinity();
  my_unaffected_min  = -std::numeric_limits<double>::infinity();
  my_unaffected_max  =  std::numeric_limits<double>::infinity();
}

inline void 
ind_filter_trait::set_trait(size_t trait)
{
  my_trait = trait;
}

inline void 
ind_filter_trait::set_min(double min)
{
  my_min = min;
} 

inline void 
ind_filter_trait::set_max(double max)
{
  my_max = max;
}

inline void 
ind_filter_trait::set_affected_min(double affected_min)
{
  my_affected_min = affected_min;
}

inline void 
ind_filter_trait::set_affected_max(double affected_max)
{
  my_affected_max = affected_max;
}

inline void 
ind_filter_trait::set_unaffected_min(double unaffected_min)
{
  my_unaffected_min = unaffected_min;
}

inline void 
ind_filter_trait::set_unaffected_max(double unaffected_max)
{
  my_unaffected_max = unaffected_max;
}

inline void
ind_filter_trait::set_affected_range(double min, double max)
{
  my_affected_min = min;
  my_affected_max = max;
}

inline void
ind_filter_trait::set_unaffected_range(double min, double max)
{
  my_unaffected_min = min;
  my_unaffected_max = max;
}

inline void
ind_filter_trait::set_threshold(double threshold)
{
  my_unaffected_min = -std::numeric_limits<double>::infinity();
  my_unaffected_max =  threshold;
  my_affected_min   =  threshold;
  my_affected_max   =  std::numeric_limits<double>::infinity();
}

inline bool
ind_filter_trait::operator ==(const ind_filter_trait& other) const
{
  return my_trait            == other.my_trait          &&
         my_min              == other.my_min            &&
         my_max              == other.my_max            &&
         my_affected_min     == other.my_affected_min   &&
         my_affected_max     == other.my_affected_max   &&
         my_unaffected_min   == other.my_unaffected_min &&
         my_unaffected_max   == other.my_unaffected_max;
}

inline bool
ind_filter_trait::operator !=(const ind_filter_trait& other) const
{
  return !(*this == other);
}

//============================================================================
// IMPLEMENTATION:  pair_filter_trait
//============================================================================
//
inline pair_filter_trait::affection_status
pair_filter_trait::mask_to_status(mask m)
{
  switch(m)
  {
    case CONCORD_UNAFF_MASK:
      return CONCORD_UNAFF;
    case DISCORD_MASK:
      return DISCORD;
    case CONCORD_AFF_MASK:
      return CONCORD_AFF;
    case UNINFORM_MASK:
      return UNINFORM;
    default:
      return NULL_STATUS;
  }  
}

inline pair_filter_trait::mask
pair_filter_trait::status_to_mask(affection_status s)
{
  switch(s)
  {
    case CONCORD_UNAFF:
      return CONCORD_UNAFF_MASK;
    case DISCORD:
      return DISCORD_MASK;
    case CONCORD_AFF:
      return CONCORD_AFF_MASK;
    case UNINFORM:
      return UNINFORM_MASK;
    default:
      return NULL_STATUS_MASK;
  }  
}

inline
pair_filter_trait::pair_filter_trait(size_t t, unsigned int a, double min, double max)
      : ind_filter_trait(t, min, max), my_status(a) 
{
  valid_status  =     (my_status & UNINFORM_MASK) || (my_status & CONCORD_UNAFF_MASK)
                   || (my_status & DISCORD_MASK)  || (my_status & CONCORD_AFF_MASK);
}
    
inline unsigned int
pair_filter_trait::get_status() const
{
  return my_status;
}

inline void 
pair_filter_trait::clear()
{
  ind_filter_trait::clear();
  my_status     = ALL_INFORMATIVE_STATUSES;
  valid_status  = true;
}

inline void 
pair_filter_trait::set_statuses(unsigned int status)
{
  my_status     = status;
  valid_status  =    (my_status & UNINFORM_MASK) || (my_status & CONCORD_UNAFF_MASK)
                  || (my_status & DISCORD_MASK)  || (my_status & CONCORD_AFF_MASK);
}

inline void
pair_filter_trait::set_status(affection_status s)
{
  my_status = status_to_mask(s);
  valid_status  =    (my_status & UNINFORM_MASK) || (my_status & CONCORD_UNAFF_MASK)
                  || (my_status & DISCORD_MASK)  || (my_status & CONCORD_AFF_MASK);
}

inline bool
pair_filter_trait::operator ==(const pair_filter_trait& other) const
{
  return ind_filter_trait::operator ==(other) &&
         my_status     == other.my_status     &&
         valid_status  == other.valid_status;
}

inline bool
pair_filter_trait::operator !=(const pair_filter_trait& other) const
{
  return !(*this == other);
}

//============================================================================
// IMPLEMENTATION:  ind_filter
//============================================================================
//
inline
ind_filter::ind_filter() 
{}
    
inline
ind_filter::ind_filter(ind_filter_trait trait)
{
  my_traits.push_back(trait);
}
  
inline void 
ind_filter::add_trait(ind_filter_trait trait)
{
  my_traits.push_back(trait);
}

inline void 
ind_filter::set_trait(ind_filter_trait trait)
{
  my_traits.remove(trait);
  my_traits.push_back(trait);
}

inline void 
ind_filter::clear_traits()
{
  my_traits.clear();
}

//============================================================================
// IMPLEMENTATION:  pair_filter
//============================================================================
//
inline
pair_filter::pair_filter() 
{}
    
inline
pair_filter::pair_filter(pair_filter_trait trait)
{
  my_traits.push_back(trait);
}
  
inline void 
pair_filter::add_trait(pair_filter_trait trait)
{
  my_traits.push_back(trait);
}

inline void 
pair_filter::set_trait(pair_filter_trait trait)
{
  my_traits.remove(trait);
  my_traits.push_back(trait);
}

inline void 
pair_filter::clear_traits()
{
  my_traits.clear();
}

//============================================================================
// IMPLEMENTATION:  filtering_pair_generator_rep
//============================================================================
//
inline
filtering_pair_generator_rep::filtering_pair_generator_rep()
{}

inline
filtering_pair_generator_rep::filtering_pair_generator_rep(const pair_generator& generator, 
                                                           pair_filter filter) 
      : my_generator(generator), my_filter(filter)
{
  my_end_iterator = my_generator.end();
}

inline
filtering_pair_generator_rep::filtering_pair_generator_rep(const filtering_pair_generator_rep& other)
{
  my_generator      = other.my_generator;
  my_end_iterator   = other.my_end_iterator;
  my_filter         = other.my_filter;
  my_cache          = other.my_cache;
}

inline pair_generator*             
filtering_pair_generator_rep::get_generator() 
{
  return &my_generator;
}

inline pair_generator::iterator*             
filtering_pair_generator_rep::get_end() 
{
  return &my_end_iterator;
}
    
inline pair_filter*             
filtering_pair_generator_rep::get_filter() 
{
  return &my_filter;
}

inline void                                  
filtering_pair_generator_rep::set_generator(const pair_generator& generator)
{
  my_generator = generator;
  my_end_iterator = my_generator.end();
}

inline void                                 
filtering_pair_generator_rep::set_filter(const pair_filter& filter)
{
  my_filter = filter;
}

//============================================================================
// IMPLEMENTATION:  filtering_pair_generator_rep::base_iterator
//============================================================================
//
inline
filtering_pair_generator_rep::base_iterator::base_iterator(const base_iterator& other)
      : my_rep(other.my_rep), my_iterator(other.my_iterator), initializing(other.initializing)
{ }

inline filtering_pair_generator_rep::base_iterator&
filtering_pair_generator_rep::base_iterator::operator =(const base_iterator& other)
{
  if(&other != this)
  {
    my_rep = other.my_rep;       // Attach to new representation.
    my_iterator = other.my_iterator;
    initializing = other.initializing;
  }
  return *this;
}
  
inline
filtering_pair_generator_rep::base_iterator::~base_iterator()
{ }

inline bool            
filtering_pair_generator_rep::base_iterator::operator ==(const base_iterator& other) const
{
  return (my_iterator == other.my_iterator);
}

inline bool            
filtering_pair_generator_rep::base_iterator::operator !=(const base_iterator& other) const
{
  return !(operator ==(other));
}
      
inline
filtering_pair_generator_rep::base_iterator::base_iterator(rep_ptr_type rep,
                                                           pair_generator::iterator iter) 
      : my_rep(rep), my_iterator(iter), initializing(true)
{ }
      
        
//============================================================================
// IMPLEMENTATION:  filtering_pair_generator_rep::iterator
//============================================================================
//
inline
filtering_pair_generator_rep::iterator::iterator(const iterator& other)
      : base_iterator(other)
{}

inline filtering_pair_generator_rep::iterator&
filtering_pair_generator_rep::iterator::operator =(const iterator& other)
{
  if(&other != this)
  {
    my_rep = other.my_rep;       // Attach to new representation.
    my_iterator = other.my_iterator;
    initializing = other.initializing;
  }
  return *this;                
}

inline
filtering_pair_generator_rep::iterator::~iterator()
{}

inline pair_generator::relative_pair&      
filtering_pair_generator_rep::iterator::operator *() 
{
  return *my_iterator;
}

inline pair_generator::relative_pair*      
filtering_pair_generator_rep::iterator::operator ->()
{
  return &(*my_iterator);
}

inline void
filtering_pair_generator_rep::iterator::print_references()
{
  std::cout << "my_rep:  " << ((void*)my_rep.get()) << "     ";
  std::cout << "no. of references: " << my_rep.use_count();
}

inline        
filtering_pair_generator_rep::iterator::iterator(rep_ptr_type rep,
                                                 pair_generator::iterator iter) 
      : base_iterator(rep, iter)
{
  // Move to first pair which passes filters.
  if(my_iterator != *(my_rep->get_end()))
  {
    ++(*this);
  }
}

//============================================================================
// IMPLEMENTATION:  filtering_pair_generator_rep::const_iterator
//============================================================================
//
inline
filtering_pair_generator_rep::const_iterator::const_iterator(const const_iterator& other)
      : base_iterator(other)
{}

inline filtering_pair_generator_rep::const_iterator&
filtering_pair_generator_rep::const_iterator::operator =(const const_iterator& other)
{
  if(&other != this)
  {
    my_rep = other.my_rep;       // Attach to new representation.
    my_iterator = other.my_iterator;
    initializing = other.initializing;
  }
  return *this;                
}

inline
filtering_pair_generator_rep::const_iterator::~const_iterator()
{}

inline const pair_generator::relative_pair&      
filtering_pair_generator_rep::const_iterator::operator *() 
{
  return *my_iterator;
}

inline const pair_generator::relative_pair*      
filtering_pair_generator_rep::const_iterator::operator ->()
{
  return &(*my_iterator);
}

inline void
filtering_pair_generator_rep::const_iterator::print_references()
{
  std::cout << "my_rep:  " << ((void*)my_rep.get()) << "     ";
  std::cout << "no. of references: " << my_rep.use_count();
}

inline        
filtering_pair_generator_rep::const_iterator::const_iterator(rep_ptr_type rep,
                                                             pair_generator::iterator iter) 
      : base_iterator(rep, iter)
{
  // Move to first pair which passes filters.
  if(my_iterator != *(my_rep->get_end()))
  {
    ++(*this);
  }
}
        
//============================================================================
// IMPLEMENTATION:  filtering_pair_generator
//============================================================================
//
inline
filtering_pair_generator::filtering_pair_generator()
{
  my_rep = rep_ptr_type(new filtering_pair_generator_rep());
}

inline
filtering_pair_generator::filtering_pair_generator(const pair_generator& generator, 
                                                   pair_filter filter)
{
  my_rep = rep_ptr_type(new filtering_pair_generator_rep(generator, filter));
}

inline
filtering_pair_generator::filtering_pair_generator(const filtering_pair_generator& other)
      : my_rep(other.my_rep)
{ }

inline filtering_pair_generator&
filtering_pair_generator::operator =(const filtering_pair_generator& other)
{
  if(&other != this)
  {
    my_rep = other.my_rep;         // Attach to new representation.
  }
  return *this;                  
}

inline
filtering_pair_generator::~filtering_pair_generator()
{ }

inline filtering_pair_generator_rep::iterator           
filtering_pair_generator::begin()
{
  iterator temp(my_rep, my_rep->get_generator()->begin());
  return temp;
}

inline filtering_pair_generator_rep::iterator           
filtering_pair_generator::end()
{
  iterator temp(my_rep, my_rep->get_generator()->end());
  return temp;
}

inline filtering_pair_generator_rep::const_iterator           
filtering_pair_generator::begin() const
{
  const_iterator temp(my_rep, my_rep->get_generator()->begin());
  return temp;
}

inline filtering_pair_generator_rep::const_iterator           
filtering_pair_generator::end() const
{
  const_iterator temp(my_rep, my_rep->get_generator()->end());
  return temp;
}

inline pair_generator
filtering_pair_generator::get_pair_generator() const
{
  return *(my_rep->get_generator());
}
    
inline pair_filter            
filtering_pair_generator::get_filter() const
{
  return *(my_rep->get_filter());
}

inline void
filtering_pair_generator::set_pair_generator(pair_generator& generator)
{
  make_unique();
  my_rep->set_generator(generator);
}
    
inline void                    
filtering_pair_generator::set_filter(const pair_filter& filter)
{
  make_unique();
  my_rep->set_filter(filter);
}

inline void
filtering_pair_generator::print_references()
{
  std::cout << "my_rep:  " << ((void*)my_rep.get()) << "     ";
  std::cout << "no. of references: " << my_rep.use_count();
}

// - If generator representation is shared, make a copy and refer to it.
//  
inline filtering_pair_generator::rep_ptr_type      
filtering_pair_generator::make_unique()
{
  if(my_rep.use_count() > 1)
  {
    my_rep = rep_ptr_type(new filtering_pair_generator_rep(*my_rep));  // Attach to new representation.
  }
  return my_rep;
}
    
//============================================================================
// IMPLEMENTATION:  Helper functions not belonging to a single class.
//============================================================================
//
inline bool
operator ==(const pair_generator::iterator& first, const pair_generator::const_iterator& second) 
{
  return first == second;
}

inline bool
operator !=(const pair_generator::iterator& first, const pair_generator::const_iterator& second) 
{
  return first != second;
}

inline bool
operator ==(const pair_generator::const_iterator& first, const pair_generator::iterator& second) 
{
  return first == second;
}

inline bool
operator !=(const pair_generator::const_iterator& first, const pair_generator::iterator& second) 
{
  return first != second;
}

inline bool
operator ==(const filtering_pair_generator::iterator& first, const filtering_pair_generator::const_iterator& second) 
{
  return first == second;
}

inline bool
operator !=(const filtering_pair_generator::iterator& first, const filtering_pair_generator::const_iterator& second) 
{
  return first != second;
}

inline bool
operator ==(const filtering_pair_generator::const_iterator& first, const filtering_pair_generator::iterator& second) 
{
  return first == second;
}

inline bool
operator !=(const filtering_pair_generator::const_iterator& first, const filtering_pair_generator::iterator& second) 
{
  return first != second;
}



