#include "sampling/MemberDataSample.h"

namespace SAGE {
namespace SAMPLING {

MemberDataSample::MemberDataSample(const FPED::Multipedigree& rmp, cerrorstream& err) 
    : my_errors(err), my_mp(rmp)
{
	my_bad_field_id.invalidate();
	my_valid.resize(my_mp.member_count());
	clearData();
	updateCounts();
}


MemberDataSample::MemberDataSample(const MemberDataSample& other) 
    : my_summary_table(other.my_summary_table), my_bad_field_id(other.my_bad_field_id),
      my_errors(other.my_errors), my_mp(other.my_mp) 
{
  my_bad_group                = other.my_bad_group;
  my_field_id_table           = other.my_field_id_table;
  my_field_vector             = other.my_field_vector;
  my_groups                   = other.my_groups;
  my_import_status            = other.my_import_status;

  my_total_individual_count   = other.my_total_individual_count;
  my_valid                    = other.my_valid;
  my_valid_individual_count   = other.my_valid_individual_count;
  my_invalid_individual_count = other.my_invalid_individual_count;
  my_valid_singleton_count    = other.my_valid_singleton_count;
  my_invalid_singleton_count  = other.my_invalid_singleton_count;  
}


MemberDataSample::~MemberDataSample()
{}


bool 
MemberDataSample::clearData()
{
  my_field_vector.clear();
  my_field_id_table.clear();
  my_groups.clear();  
  my_import_status = READY_FOR_IMPORT_AND_FINALIZE;

  addGroup("__ALL__");

  return true;
}


void
MemberDataSample::reset()
{
  clearData();
  my_valid.resize(my_mp.member_count(), false);
  updateCounts();
  my_summary_table = OUTPUT::Table();
}


const FieldID& 
MemberDataSample::importField(const string& source_trait_name, const string& group_name,
                              const string& new_field_name, unsigned long flags)
{
  // Check status:
  if(getImportStatus() != READY_FOR_IMPORT_AND_FINALIZE)
    return my_bad_field_id;

  if(my_mp.info().trait_exists(source_trait_name) == false)
  {
    return my_bad_field_id;
  }
  
  // Add a new field to our set
  
  const FieldID& fid = 
    fieldExists(group_name, new_field_name) ? 
      getField(group_name, new_field_name).getFieldID() :
      addField(source_trait_name, group_name, new_field_name, false, flags);
                
  // Check status to make sure it actually did add it.
  if(!fid.isValid())
    return my_bad_field_id;
  
  Field & field = getField(fid);
  
  // 2.3. Import data into field:

  size_t trait_id = my_mp.info().trait_find(source_trait_name);

  for(uint i = 0; i < my_mp.member_count(); ++i)
  {
    size_t pedigree_member_id = my_mp.member_index(i).index();

    field.setOrigValue(i, my_mp.member_index(i).pedigree()->info().trait(pedigree_member_id, trait_id));
  }
  
  return fid;
}


const FieldID& 
MemberDataSample::createField(string group_name, 
                              string field_name, unsigned long flags)
{
  return addField("__USER__", group_name, field_name, true, flags);
}


bool
MemberDataSample::addGroup(const string & group_name)
{
  if(groupExists(group_name) == true)
  {
    return false;
  }
  else
  {
    my_groups.insert(make_pair(group_name, FieldGroupType()));

    return true;
  }
}


int 
MemberDataSample::invalidateIndividual(size_t i)
{
  my_valid[i] = false;

  for(FieldIterator field_itr = getFieldBegin(); field_itr != getFieldEnd(); ++field_itr)
    field_itr->setOrigValue(i, std::numeric_limits<double>::quiet_NaN());

  return 0;
}


MemberDataSample::ImportStatusEnum
MemberDataSample::getImportStatus() const
{
  return my_import_status;
}


int 
MemberDataSample::finalizeData()
{
  return finalizeData(DefaultValidator());
}


int
MemberDataSample::finalizeData(const IndividualValidator& validator)
{
  // Check status:
  if(getImportStatus() != READY_FOR_IMPORT_AND_FINALIZE)
    return 1;

  // Validate all individuals:
  validateIndividuals(validator);

  // Iterate across fields and process their input flags:
  for(FieldIterator field_itr = getFieldBegin(); field_itr != getFieldEnd(); ++field_itr)
  {
    if(field_itr->getUserDefined() == false)
    {
      field_itr->finalizeData();
    }
  }

  // Update various counts:
  updateCounts();

  // Set new status:
  my_import_status = READY_FOR_CREATE_AND_FINALIZE;

  // Return success:
  return 0;
}


int
MemberDataSample::finalizeUserCreatedData()
{
  // Check status:
  if(getImportStatus() != READY_FOR_CREATE_AND_FINALIZE)
    return 1;

  for(FieldIterator field_itr = getFieldBegin(); field_itr != getFieldEnd(); ++field_itr)
    if(field_itr->getUserDefined() == true)
      field_itr->finalizeData();

  // Set status:
  my_import_status = COMPLETELY_FINALIZED;

  // Populate summary table:
  my_summary_table = (OUTPUT::Table("Sample description")
    << (OUTPUT::TableRow() << "Number of pedigrees in dataset"           << getMultipedigree().info().get_source_rped()->pedigree_count())
    << (OUTPUT::TableRow() << "Number of analyzable pedigrees"           << getMultipedigree().pedigree_count())
    <<  OUTPUT::Table::INSERT_BLANK_ROW()
    << (OUTPUT::TableRow() << "Number of individuals in dataset"         << getMultipedigree().info().get_source_rped()->member_count())
    << (OUTPUT::TableRow() << "Number of analyzable individuals"         << getMultipedigree().member_count())
    << (OUTPUT::TableRow() << "Number of analyzable invalid individuals" << getInvalidIndividualCount())
    << (OUTPUT::TableRow() << "Number of analyzable valid individuals"   << getValidIndividualCount()));

  return 0;
}


void
MemberDataSample::validateIndividuals(const IndividualValidator& validator)
{
  // 0. Set up the IndividualTraitData (with empty values)

        IndividualTraitData trait_data;

        for(FieldConstIterator field_itr = getFieldBegin(); field_itr != getFieldEnd(); ++field_itr)
        {
          GroupInfo null_group;
          null_group.name = field_itr->getGroupName();

          GroupInfoIterator current_group = trait_data.find(field_itr->getGroupName());

          if(current_group == trait_data.data.end())
          {
            trait_data.data.push_back(null_group);

            current_group = trait_data.find(field_itr->getGroupName());
          }

          TraitInfo null_trait;

          null_trait.name            = field_itr->getFieldName();
          null_trait.allow_averaging = field_itr->getInputFlags() & Field::ALLOW_AVERAGING;
          null_trait.value           = std::numeric_limits<double>::quiet_NaN();
          null_trait.user_defined    = field_itr->getUserDefined();

          TraitInfoIterator current_trait = current_group->find(field_itr->getFieldName());

          if(current_trait == current_group->data.end())
          {
            current_group->data.push_back(null_trait);

            current_trait = current_group->find(field_itr->getFieldName());
          }
        }

  // 1. Loop through all individuals:

        for(size_t i = 0; i < getTotalIndividualCount(); ++i)
        {
  // 1.1. For this individual, populate the TraitSet with values by looping through the TraitSet:

          for(GroupInfoIterator group_info_itr  = trait_data.begin ();
                                group_info_itr != trait_data.end   (); ++group_info_itr)
          {
            for(TraitInfoIterator trait_info_itr  = group_info_itr->begin ();
                                  trait_info_itr != group_info_itr->end   (); ++trait_info_itr)
            {
              trait_info_itr->value = getField(group_info_itr->name, trait_info_itr->name).getOrigValue(i);
            }
          }

  // 1.2. If this individual is invalid, erase all his trait values:

          if(validator.isValid(i, trait_data) == false)
            invalidateIndividual(i);
	        else
	          my_valid[i] = true;
        }
}

const FieldID & 
MemberDataSample::addField(
                  const string&     source_name,
                  const string&     group_name,
                  const string&     field_name,
                  bool              user_def,
                  unsigned long     flags)
{
  // 2.1. Add group if necessary:
	if(groupExists(group_name) == false)
	  addGroup(group_name);

  // 2.2. Create Field and set basic options:
	my_field_vector      . push_back(Field(my_mp.member_count()));

	Field& field = my_field_vector.back();

	field.setOrigName    (source_name);
        field.setUserDefined (user_def);
	field.setFieldName   (field_name);
	field.setGroupName   (group_name);
	field.setFieldID     (FieldID(getFieldCount()-1));
	field.setGroupID     (GroupID(getFieldCount(group_name)-1));
	field.setInputFlags  (flags);

  // 2.6. Add the field to the appropriate group lists:

	getGroup("__ALL__")  . push_back(field.getFieldID());
	getGroup(group_name) . push_back(field.getFieldID());

  // 2.7. Add the name lookup for this field:

	my_field_id_table.insert(make_pair(make_pair(group_name, field_name), field.getFieldID()));

  // 2.8. Return field id:

        return field.getFieldID();
}

//==============================================================
//
//  dumpTraitValues()
//
//==============================================================
void 
MemberDataSample::dumpTraitValues() const
{
  // Print individuals:

  OUTPUT::Table t("Trait values");
  
  t << OUTPUT::TableColumn("ID") 
    << OUTPUT::TableColumn("Ped")
    << OUTPUT::TableColumn("Name")
    << OUTPUT::TableColumn("Valid");

  for(FieldConstIterator field_itr = getFieldBegin(); field_itr != getFieldEnd(); ++field_itr)
  {
    t << OUTPUT::Table::BEGIN_COLUMN_GROUP(field_itr->getGroupName())
      << OUTPUT::TableColumn(field_itr->getFieldName() + " [Orig]")
      << OUTPUT::TableColumn(field_itr->getFieldName() + " [Repl.]")
      << OUTPUT::TableColumn(field_itr->getFieldName() + " [Adj]");
  }

  OUTPUT::TableRow mean_row,
                   mean_adj_row,
                   stdev_row,
                   stdev_adj_row,
                   aa_row;

  mean_row      << "" << "" << "" << "";
  mean_adj_row  << "" << "" << "" << "";
  stdev_row     << "" << "" << "" << "";
  stdev_adj_row << "" << "" << "" << "";
  aa_row        << "" << "" << "" << "";
  
  for(FieldConstIterator field_itr = getFieldBegin(); field_itr != getFieldEnd(); ++field_itr)
  {
    mean_row      << "mean      =" <<  field_itr->getMean       ()                                           << "";
    stdev_row     << "stdev     =" <<  field_itr->getStdev      ()                                           << "";
    mean_adj_row  << "mean_adj  =" << (field_itr->getInputFlags () & Field::MEAN_ADJUST     ? "YES" : "NO ") << "";
    stdev_adj_row << "stdev_adj =" << (field_itr->getInputFlags () & Field::STDEV_ADJUST    ? "YES" : "NO ") << "";
    aa_row        << "avgng     =" << (field_itr->getInputFlags () & Field::ALLOW_AVERAGING ? "YES" : "NO ") << "";
  }
    
  t << OUTPUT::Table::BEGIN_ROW_GROUP("Summary information")
    << mean_row 
    << mean_adj_row 
    << stdev_row 
    << stdev_adj_row 
    << aa_row 
    << OUTPUT::Table::BEGIN_ROW_GROUP("Individual values");

  for(size_t i = 0; i < getTotalIndividualCount(); ++i)
  {
    OUTPUT::TableRow r;
    
    r << i 
      << getMultipedigree().member_index(i).pedigree()->name()
      << getMultipedigree().member_index(i).name()
      << (isValid(i) ? "Valid" : "Invalid");
    
    for(FieldConstIterator field_itr = getFieldBegin(); field_itr != getFieldEnd(); ++field_itr)
    {
      r << field_itr->getOrigValue(i);

      if(field_itr->isAdjValuePresent(i))
        r << field_itr->getReplacedValue(i) << field_itr->getAdjValue(i);
      else
        r << "<missing>" << "<missing>";
    }
    
    t << r;
  }

  cout << t;
}


void 
MemberDataSample::updateCounts()
{
  // - 6-11-7 added singleton counts. djb
  //
  my_total_individual_count   = my_valid.size();
  my_valid_individual_count   = 0;
  my_invalid_individual_count = 0;
  my_valid_singleton_count    = 0;
  my_invalid_singleton_count  = 0;

  for(size_t i = 0; i < my_valid.size(); ++i)
  {
    bool  singleton = MPED::mp_utilities::nuclear_family_count(&(my_mp.member_index(i))) == 0;
    bool  valid = my_valid[i];
    if(valid)
    {
      ++my_valid_individual_count;
    
      if(singleton)
      {
        ++my_valid_singleton_count;
      }
    }
    else
    {
      ++my_invalid_individual_count;
      
      if(singleton)
      {
        ++my_invalid_singleton_count;
      }
    }
  }
}


void 
MemberDataSample::populateUserCreatedField(const string & group_name, const string & field_name, const TraitValueCalculator & calc)
{
  Field & field = getField(group_name, field_name);

  for(size_t i = 0; i < getTotalIndividualCount(); ++i)
  {
    field.setOrigValue(i, isValid(i) ? 
      calc(i, (const MemberDataSample)*this) : 
      std::numeric_limits<double>::quiet_NaN());
  }
}

} // End namespace SAMPLING
} // End namespace SAGE

