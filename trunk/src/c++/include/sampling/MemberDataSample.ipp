namespace SAGE {
namespace SAMPLING {

//==============================================================
//  getErrorstream()
//==============================================================
inline cerrorstream& 
MemberDataSample::getErrorstream()
{
  return my_errors;
}

//==============================================================
//  getRefMultiPedigree()
//==============================================================
inline const FPED::Multipedigree& 
MemberDataSample::getMultipedigree() const
{
  return my_mp;
}

//=================================================================
//  fieldExists() #1
//=================================================================
inline bool 
MemberDataSample::fieldExists(const string& group_name, const string& field_name) const
{
  if(my_field_id_table.find(make_pair(group_name, field_name)) == my_field_id_table.end())
    return false;
  else
    return true;
}
     
//=================================================================
//  fieldExists() #2
//=================================================================
inline bool 
MemberDataSample::fieldExists(const string& group_name, const GroupID& group_id) const
{
  if(groupExists(group_name))
  {
    if(group_id.id < getGroup(group_name).size())
      return true;
    else
      return false;
  }
  else
  {
    return false;
  }
}

//=================================================================
//  groupExists()
//=================================================================
inline bool
MemberDataSample::groupExists(const string& group_name) const
{
  return my_groups.find(group_name) != my_groups.end();
}

//=================================================================
//  getGroup() CONST
//=================================================================
inline const FieldGroupType& 
MemberDataSample::getGroup(const string& group_name) const
{
  std::map<string, FieldGroupType>::const_iterator group_itr = my_groups.find(group_name);

  if(group_itr == my_groups.end())
  {
    my_errors << "Error: Group '" << group_name << "' requested but doesn't exist." << endl;

    return my_bad_group;
  }
  else
  {
    return group_itr->second;
  }
}

//=================================================================
//  getGroup() NON-CONST
//=================================================================
inline FieldGroupType&
MemberDataSample::getGroup(const string & group_name)
{
  // You may at this point be wondering why in the world I chose to implement
  // this function as I did. I could have simply copied over the code from the
  // const version of getGroup() (above). That would have worked perfectly fine.
  // 
  // But duplicating code is not a good idea. Whenever possible, you should avoid
  // any code duplication.
  // 
  // Ok, so that basically means that the non-const version of getGroup() will
  // have to invoke the const version of getGroup(). This isn't so much of a
  // problem, as you can simply const_cast() the return type into a non-const
  // reference. That's easy.
  // 
  // The tricky thing is making sure the correct version of getGroup() is
  // invoked! If the code simply const_cast'ed getGroup(), then the compiler
  // would invoke the non-const version, and you would end up with infinite
  // recursion. The solution, then, is to construct a const reference to the MemberDataSample object
  // itself, and then invoke getGroup() on the const MemberDataSample reference.
  // 
  // Cool, huh? --sag 5 Aug 2004

  const MemberDataSample& temp = *this;
  return const_cast<FieldGroupType&> (temp.getGroup(group_name));
}

//=========================================
//  Individual information
//=========================================

inline const FPED::Member& 
MemberDataSample::getIndividual(size_t i) const 
{ 
  return getMultipedigree().member_index(i); 
}

inline size_t  
MemberDataSample::getTotalIndividualCount() const 
{ 
  return my_total_individual_count;   
}

inline size_t  
MemberDataSample::getValidIndividualCount() const 
{ 
  return my_valid_individual_count;   
}

inline size_t  
MemberDataSample::getInvalidIndividualCount() const 
{ 
  return my_invalid_individual_count; 
}

inline size_t  
MemberDataSample::getValidSingletonCount() const 
{ 
  return my_valid_singleton_count;   
}

inline size_t  
MemberDataSample::getInvalidSingletonCount() const 
{ 
  return my_invalid_singleton_count; 
}

inline bool 
MemberDataSample::isValid(size_t i) const 
{ 
  return my_valid[i];                 
}

inline double 
MemberDataSample::getOrigValue(size_t i, const FieldID& field_id) const
{
  return getField(field_id).getOrigValue(i);
}

inline double
MemberDataSample::getOrigValue(size_t i, const string& group_name, const GroupID& group_id) const
{
  return getField(group_name, group_id).getOrigValue(i);
}

inline double
MemberDataSample::getOrigValue(size_t i, const string& group_name, const string& field_name) const
{
  return getField(group_name, field_name).getOrigValue(i);
}

inline double 
MemberDataSample::getReplacedValue(size_t i, const FieldID& field_id) const
{
  return getField(field_id).getReplacedValue(i);
}

inline double
MemberDataSample::getReplacedValue(size_t i, const string& group_name, const GroupID& group_id) const
{
  return getField(group_name, group_id).getReplacedValue(i);
}

inline double
MemberDataSample::getReplacedValue(size_t i, const string& group_name, const string& field_name) const
{
  return getField(group_name, field_name).getReplacedValue(i);
}

inline double 
MemberDataSample::getAdjValue(size_t i, const FieldID& field_id) const
{
  return getField(field_id).getAdjValue(i);
}

inline double
MemberDataSample::getAdjValue(size_t i, const string& group_name, const GroupID& group_id) const
{
  return getField(group_name, group_id).getAdjValue(i);
}

inline double
MemberDataSample::getAdjValue(size_t i, const string& group_name, const string& field_name) const
{
  return getField(group_name, field_name).getAdjValue(i);
}

inline size_t
MemberDataSample::getFieldCount() const
{
  return my_field_vector.size();
}

inline size_t
MemberDataSample::getFieldCount(const string& group_name) const
{
  return getGroup(group_name).size();
}

inline size_t
MemberDataSample::getGroupCount() const
{
  return my_groups.size();
}

//=====================================================================
//  Field extraction
//=====================================================================

inline Field& 
MemberDataSample::getField(const FieldID& id) 
{ 
  return my_field_vector[id.id]; 
}

inline const Field& 
MemberDataSample::getField(const FieldID &id) const 
{ 
  return my_field_vector[id.id]; 
}

inline Field& 
MemberDataSample::getField(const string& group_name, const GroupID& group_id) 
{ 
  if(!fieldExists(group_name, group_id))
    SAGE_internal_error();

  return my_field_vector[getGroup(group_name)[group_id.id].id]; 
}

inline const Field& 
MemberDataSample::getField(const string& group_name, const GroupID& group_id) const 
{ 
  if(!fieldExists(group_name, group_id))
    SAGE_internal_error();

  return my_field_vector[getGroup(group_name)[group_id.id].id]; 
}

inline Field& 
MemberDataSample::getField(const string& group_name, const string& field_name) 
{ 
  if(!fieldExists(group_name, field_name))
    SAGE_internal_error();

  return my_field_vector[my_field_id_table.find(make_pair(group_name, field_name))->second.id]; 
}

inline const Field& 
MemberDataSample::getField (const string& group_name, const string& field_name) const 
{ 
  if(!fieldExists(group_name, field_name))
    SAGE_internal_error();

  return my_field_vector[my_field_id_table.find(make_pair(group_name, field_name))->second.id]; 
}


//======================================================================
//  Field traversal
//======================================================================

inline FieldIterator 
MemberDataSample::getFieldBegin() 
{ 
  return  getFieldBegin("__ALL__"); 
}

inline FieldIterator 
MemberDataSample::getFieldEnd() 
{ 
  return  getFieldEnd("__ALL__"); 
}

inline FieldIterator 
MemberDataSample::getFieldBegin(const string& group_name) 
{ 
  return  FieldIterator(&my_field_vector, getGroup(group_name).begin()); 
}

inline FieldIterator 
MemberDataSample::getFieldEnd(const string& group_name) 
{ 
  return  FieldIterator(&my_field_vector, getGroup(group_name).end()); 
}

inline FieldConstIterator 
MemberDataSample::getFieldBegin() const 
{ 
  return  getFieldBegin ("__ALL__"); 
}

inline FieldConstIterator 
MemberDataSample::getFieldEnd() const 
{ 
  return  getFieldEnd("__ALL__"); 
}

inline FieldConstIterator 
MemberDataSample::getFieldBegin(const string & group_name) const 
{ 
  return  FieldConstIterator(&my_field_vector, getGroup(group_name).begin()); 
}

inline FieldConstIterator 
MemberDataSample::getFieldEnd(const string & group_name) const 
{ 
  return FieldConstIterator(&my_field_vector, getGroup(group_name).end()); 
}

} // End namespace SAMPLING
} // End namespace SAGE
