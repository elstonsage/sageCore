#include "segreg/PedigreeDataSet.h"
#include "mped/mp_utilities.h"

namespace SAGE {
namespace SEGREG {

/// Helper function which returns \c true if the member is in a valid
/// structure, \c false otherwise.  Valid structures are defined as either
/// being unconnected (no structure), or the subpedigree is loopless.
bool
  member_in_valid_structure
    (const RPED::Member& mem)
{
  return mem.is_unconnected() ||
         MPED::mp_utilities::no_loops(*mem.subpedigree());
}

/// Assigns arbitrary sexes to individuals without sexes in the FPED::Multipedigree.
/// Note that sexes are assigned in consistent fashion, such that parents are
/// always sexed differently.
void
  assign_arbitrary_sexes
    (FPED::Multipedigree& fped)
{
  // Assign arbitrary sexes, so analyses won't break.
  for(size_t i = 0; i < fped.member_count(); ++i)
  {
    if(fped.member_index(i).is_sex_unknown())
    {
      fped.member_index(i).set_sex(MPED::SEX_XMALE);
      fped.member_index(i).pedigree()->build();
    }
  }
}

/// Constructor for the PedigreeDataSet.
///
/// \param mped  The Raw multipedigree data
/// \param model The analysis model
PedigreeDataSet::PedigreeDataSet
  (const RPED::MultiPedigree& mped,
   const model&               model,
   APP::Output_Streams&       output)
{
  create_dataset(mped,model);

  if(this->is_valid_with(model))
  {
    build_subpedigree_list();
    build_member_list();
    build_unconnected_list();
    
    assign_arbitrary_sexes(*my_data);
  }
  else
  {
    report_dataset_errors(model, output);
    
    clear_dataset();
  }
}

/// Given multipedigree and model, creates appropriate dataset
///
/// \param mped  The Raw multipedigree data
/// \param model The analysis model
void
  PedigreeDataSet::create_dataset
    (const RPED::MultiPedigree& mped,
     const model&               model)
{
  my_data = boost::shared_ptr<FPED::Multipedigree>
                          (new FPED::Multipedigree(mped));

  FPED::MPFilterer::add_multipedigree_filtered_by_members
          (*my_data, mped, boost::bind(member_in_valid_structure, _1));

  my_data->construct();
}

/// Tests if the current dataset is valid with the model given.
/// The dataset might be invalid if:
///   - it's empty,
///   - it has bad sex information when sex is required by the analysis, or
///   - it doesn't have structure, required by the analysis.
bool
  PedigreeDataSet::is_valid_with
    (const model& model)
{
  if(this->is_empty())
    return false;

  if(this->has_sex_errors(model))
    return false;
  
  if(this->has_transmission_error(model))
    return false;
  
  return true;
}

/// Helpter function returns \c true if there are parents in the dataset
/// which have missing sex, \c false otherwise
bool has_parental_missing_sex(const FPED::Multipedigree& fped)
{
  for(size_t mem = 0; mem != fped.member_count(); ++mem)
  {
      if(fped.member_index(mem).mate_count() && 
         fped.member_index(mem).is_sex_unknown()) return true;
  }
  return false;
}

/// Returns \c true if the model requires sex to be available, but there
/// are members in the dataset which are unsexed.
bool
  PedigreeDataSet::has_sex_errors
    (const model& model)
{
  if(!model.transm_sub_model.has_sex_effect() &&
     !model.resid_sub_model.has_sex_effect())
    return false;
  
  return has_parental_missing_sex(*my_data);
}

/// Returns \c true if the model requires pedigree structures to make transmission
/// work.  This occurs whenever there are parameters in the transmission model which
/// can be maximized.
bool transmission_requires_structure(const model& model)
{
  transmission_sub_model::sm_option opt = model.transm_sub_model.option();
  
  return model.get_trans_missing()                    ||
         opt == transmission_sub_model::homog_general ||
         opt == transmission_sub_model::general       ||
         opt == transmission_sub_model::tau_ab_free;
}

bool pedigree_contains_structure(const FPED::Multipedigree& fped)
{
  for(FPED::PedigreeConstIterator ped = fped.pedigree_begin();
      ped != fped.pedigree_end();
      ++ped)
  {
    if(ped->family_count()) return true;
  }
  
  return false;
}

bool
  PedigreeDataSet::has_transmission_error
    (const model& model)
{
  return transmission_requires_structure(model) &&
         !pedigree_contains_structure(*my_data);
}

/// Clears the dataset
void
  PedigreeDataSet::clear_dataset()
{
  my_data = boost::shared_ptr<FPED::Multipedigree>();
  
  my_subpedigrees = SubpedigreeList();
  my_members      = MemberList();
  my_unconnecteds = MemberList();
}

/// Builds the list of subpedigrees in the dataset
///
void
  PedigreeDataSet::build_subpedigree_list()
{
  for(FPED::PedigreeConstIterator ped = my_data->pedigree_begin();
      ped != my_data->pedigree_end(); ++ped)
  {
    for(FPED::SubpedigreeConstIterator sped = ped->subpedigree_begin();
        sped != ped->subpedigree_end(); ++sped)
    {
      my_subpedigrees.push_back(&*sped);
    }
  }
}

/// Builds the list of members in the dataset
///
void
  PedigreeDataSet::build_member_list()
{
  for(FPED::PedigreeConstIterator ped = my_data->pedigree_begin();
      ped != my_data->pedigree_end(); ++ped)
  {
    for(FPED::MemberConstIterator mem = ped->member_begin();
        mem != ped->member_end(); ++mem)
    {
      my_members.push_back(&*mem);
    }
  }
}

/// Builds the list of unconnecteds in the dataset
///
void
  PedigreeDataSet::build_unconnected_list()
{
  for(FPED::PedigreeConstIterator ped = my_data->pedigree_begin();
      ped != my_data->pedigree_end(); ++ped)
  {
    for(FPED::MemberConstIterator unc = ped->unconnected_begin();
        unc != ped->unconnected_end(); ++unc)
    {
      my_unconnecteds.push_back(&*unc);
    }
  }
}

/// Reports errors in constructing the dataset to the errors stream in the output object
void
  PedigreeDataSet::report_dataset_errors
    (const model&         model,
     APP::Output_Streams& output)
{
  if(this->is_empty())
    report_empty_dataset_error(output);
  
  if(this->has_sex_errors(model))
    report_dataset_sex_errors(model, output);

  if(this->has_transmission_error(model))
    report_transmission_error(output);
}

/// Creates a list of strings that are the missing sex parents in the pedigree,
/// for output
void
  generate_missing_sex_parent_list
    (std::list<std::string>& missing_sex_parents,
     const FPED::Multipedigree& fped)
{
  for(size_t mem = 0; mem != fped.member_count(); ++mem)
  {
    const FPED::Member& member = fped.member_index(mem);
    
    if(member.mate_count() && member.is_sex_unknown())
      missing_sex_parents.push_back(member.pedigree()->name() + ":" + member.name());
  }
}  

/// Converts a list of strings into a single pretty-print string with commans and
/// "and" at the end.
std::string compose_printable_list(const std::list<std::string>& elts)
{
  std::string result     = "";
  std::string result_end = "";  

  std::list<std::string>::const_iterator endpoint = elts.end();
  
  if(elts.size() > 0) result_end = *(--endpoint);
  if(elts.size() > 1) result_end = *(--endpoint) + " and " + result_end;
  
  // Create a comma delimited list of the names still on the list
  for(std::list<std::string>::const_iterator i = elts.begin();
      i != endpoint; ++i)
  {
    result = result + *i + ", ";
  }
  
  return result + result_end;
}

/// Creates and returns a pretty-print string that lists all members with
/// missing sex that are also parents.
std::string get_missing_sex_parent_list(const FPED::Multipedigree& ped_data)
{
  std::list<std::string> missing_sex_parents;
  
  generate_missing_sex_parent_list(missing_sex_parents, ped_data);

  std::string result = compose_printable_list(missing_sex_parents);
  
  return result;
}

/// Report dataset errors with sex coding to the error stream from out
void
  PedigreeDataSet::report_dataset_sex_errors
    (const model&           model,
     APP::Output_Streams&   out)
{
  out.errors() << priority(error);
  
  bool transm_err = model.transm_sub_model.has_sex_effect();
  bool resid_err  = model.resid_sub_model.has_sex_effect();

  if( transm_err &&  resid_err) out.errors() << "The transmission and residual options require";
  if( transm_err && !resid_err) out.errors() << "The transmission option requires";
  if(!transm_err &&  resid_err) out.errors() << "The residual option requires";

  out.errors() << " that sex be available for all parents in this analysis. The "
               << "following individuals have missing sex and are parents: "
               << get_missing_sex_parent_list(*my_data) << ".  Analysis is "
               << "invalid and will be discontinued." << endl;
}

void
  PedigreeDataSet::report_transmission_error
    (APP::Output_Streams&   out)
{
  out.errors() << priority(error);
  
  out.errors() << "The transmission/segregation model given requires that there be "
               << "pedigree structures in the sample to be analyzed, but no valid "
               << "structures have been detected in this sample.  Analysis is "
               << "invalid and will be discontinued." << endl;
}

/// Reports an error when the dataset is empty to the error stream from out
void
  PedigreeDataSet::report_empty_dataset_error
    (APP::Output_Streams&   out)
{
  out.errors()
        << priority(error)
        << "No valid data after filtering.  Analysis will be aborted." << endl;
}

}
}

