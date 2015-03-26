#include "mlod/PedigreeAnalysisSample.h"

namespace SAGE {
namespace MLOD {

PedigreeAnalysisSample::PedigreeAnalysisSample
    (const RPED::MultiPedigree& source_rped,
     const AnalysisParameters&  analysis_params)
  : my_original_rped(source_rped),
    my_fped         (source_rped)
{
  // Initialize our data structure
  initialize_member_type_vector();

  // Do a first pass on member informativity
  classify_members_informativity(analysis_params);

  // Create the fped based upon informativity 
  create_fped(analysis_params);

  // Count the speds for statistic output sections.
  count_speds();
}

/// Initializes the vector of member types, with all members labeled as
/// MIT_UNDEFINED
void PedigreeAnalysisSample::initialize_member_type_vector()
{
  my_member_types.clear();
  my_member_types.resize(my_original_rped.member_count(), MIT_UNDEFINED);
}

/// For each member in the data set, determines (based upon the AnalysisParameters)
/// their informatity.  Members can be classified as MIT_SINGLETON,
/// MIT_DIRECTLY_INFORMATIVE, MIT_STRUCTURALLY_INFORMATIVE or
/// MIT_UNINFORMATIVE at this stage.
///
/// \param ap The analysis paramters which determine what loci are relevant.
void PedigreeAnalysisSample::classify_members_informativity(const AnalysisParameters&  ap)
{
  // Create some filter functions.  One will determine basic informativity, and one
  // for structural informativity
  typedef FPED::has_informative_loci<RPED::Member>            IsDataInf;
  typedef FPED::is_inf_within_sped_t<RPED::Member, IsDataInf> IsStructInf;

  IsDataInf is_data_inf(my_original_rped, false);
  
  // For the region, add each marker locus to the included set of checked markers 
  for(size_t i = 0; i < ap.get_region().locus_count(); ++i)
  {
    is_data_inf.set_check_status_for_locus(ap.get_region().locus(i).marker_index(), true);
  }
 
  // For each trait locus, add it to the included set of checked markers 
  for(AnalysisParameters::TraitModelList::const_iterator tr_iter  = ap.get_trait_list().begin();
      tr_iter != ap.get_trait_list().end(); ++tr_iter)
  {
    is_data_inf.set_check_status_for_locus(tr_iter->first, true);
  }

  // Create an is_inf object based upon locus_test.
  IsStructInf is_struct_inf(is_data_inf);
  
  for(size_t i = 0; i < my_member_types.size(); ++i)
  {
    const RPED::Member& mem = my_original_rped.member_index(i);
    
    // Look for singletons first.  If a singleton, it doesn't matter if it's
    // informative
    if(mem.is_unconnected())        my_member_types[i] = MIT_SINGLETON;
    else if(is_data_inf(mem))       my_member_types[i] = MIT_DIRECTLY_INFORMATIVE;
    else if(is_struct_inf(mem))     my_member_types[i] = MIT_STRUCTURALLY_INFORMATIVE;
    else                            my_member_types[i] = MIT_UNINFORMATIVE;
  }
}

/// \internal
///
/// \brief Functor for calculating if a member is informative within the PedigreeAnalysisSample's constraints
///
/// This predicate functor returns \c true if the member is directly or
/// structurally informative and \c false otherwise.  In other words, a member
/// is member PAS informative if it could contribute information to this analysis.
class is_member_PAS_informative : public unary_function<RPED::Member, bool>
{
  public:
  
    //lint --e{1509}  We know unary function has no destructor
    //lint --e{1712}  No default constructor for this
  
    /// Constructor
    ///
    /// \param memtypes The vector of member types following classify_members_informativity
    is_member_PAS_informative(const vector<PedigreeAnalysisSample::MemberInclusionType>& memtypes)
      : my_member_types(memtypes)
    {}

    /// This function should return true if m is to be included in the
    /// filtered subset, false otherwise
    ///
    /// \param m The member we're checking.
    bool operator()(const RPED::Member& m) const
    {
      PedigreeAnalysisSample::MemberInclusionType t = my_member_types[m.mpindex()];
      
      return t == PedigreeAnalysisSample::MIT_DIRECTLY_INFORMATIVE ||
             t == PedigreeAnalysisSample::MIT_STRUCTURALLY_INFORMATIVE;
    }
    
    private:
    
      /// A reference to the members' informativities.
      const vector<PedigreeAnalysisSample::MemberInclusionType>& my_member_types;
};

/// Given the AnalysisParameters, creates the FPED::Multipedigree which
/// represents the analyzable subset of the data which is analyzable
/// under the constraints contained in the AnalysisParameters.
///
/// \param ap The AnalysisParameters constraining the use of the data.
void PedigreeAnalysisSample::create_fped(const AnalysisParameters& ap)
{
  // 1. Create the temporary fped.  This contains all individuals which are 
  //    informative.
  
  FPED::Multipedigree temp_fped(my_original_rped);
  
  //lint -e{1032} -e{534} Spurious, since function *is* const; We also don't care about the results
  FPED::MPFilterer::add_multipedigree_filtered_by_members(temp_fped, my_original_rped,
           is_member_PAS_informative(my_member_types));

  temp_fped.construct();
  
  // 2. Reclassify individuals
  
  // Iterate through the temp_fped and determine which subpedigrees and members 
  // we can use based upon subpedigree size.  Change the my_member_types vector 
  // appropriately.
  
  // While performing this step, also keep track of the largest number of bits in
  // the valid data set.
  my_largest_bit_count = 0;
  
  for(FPED::PedigreeConstIterator ped = temp_fped.pedigree_begin();
      ped != temp_fped.pedigree_end(); ++ped)
  {
    // Iterate through subpedigrees and check their sizes:
    
    for(FPED::SubpedigreeConstIterator sped = ped->subpedigree_begin();
        sped != ped->subpedigree_end(); ++sped)
    {
      size_t bit_count = calculate_sped_bit_count(*sped);
      
      if(bit_count < 1)
      {
        // Subpedigree too small.  Set all members to MIT_SPED_TOO_SMALL

        for(size_t i = 0; i < sped->member_count(); ++i)
        {
          const RPED::Member& mem = *(sped->member_index(i).info().get_source_member());
          
          my_member_types[mem.mpindex()] = MIT_SPED_TOO_SMALL;
        }
      }
      else if(ap.get_max_ped_size() < bit_count)
      {
        // Subpedigree too large.  Set all members to MIT_SPED_TOO_LARGE

        for(size_t i = 0; i < sped->member_count(); ++i)
        {
          const RPED::Member& mem = *(sped->member_index(i).info().get_source_member());
          
          my_member_types[mem.mpindex()] = MIT_SPED_TOO_LARGE;
        }
      }
      else // Sped within bounds
      {
        // Check if sped is the largest we've found so far.
        
        if(my_largest_bit_count < bit_count) my_largest_bit_count = bit_count;
      }
    }

    // Iterate through unconnecteds and mark them MIT_SINGLETON_AFTER
    
    for(FPED::MemberConstIterator ucon = ped->unconnected_begin();
        ucon != ped->unconnected_end(); ++ucon)
    {
      my_member_types[ucon->info().get_source_member()->mpindex()] = MIT_SINGLETON_AFTER;
    }
  }

  // 3. Build the final FPED

  //lint -e{1032} -e{534} Spurious, since function *is* const; We also don't care about the results
  FPED::MPFilterer::add_multipedigree_filtered_by_members(my_fped, my_original_rped,
           is_member_PAS_informative(my_member_types));

  my_fped.construct();
}

void PedigreeAnalysisSample::count_speds()
{
  my_original_sped_count = 0;
  
  for(RPED::PedigreeConstIterator ped = my_original_rped.pedigree_begin();
      ped != my_original_rped.pedigree_end(); ++ped)
  {
    my_original_sped_count += ped->subpedigree_count();
  }

  my_filtered_sped_count = 0;
  
  for(FPED::PedigreeConstIterator ped = my_fped.pedigree_begin();
      ped != my_fped.pedigree_end(); ++ped)
  {
    my_filtered_sped_count += ped->subpedigree_count();
  }

  // Iterate through the speds and put them into the vector
  
  my_spedigrees.reserve(my_filtered_sped_count);

  for(FPED::PedigreeConstIterator ped = my_fped.pedigree_begin();
      ped != my_fped.pedigree_end(); ++ped)
  {
    for(FPED::SubpedigreeConstIterator sped = ped->subpedigree_begin();
        sped != ped->subpedigree_end(); ++sped)
    {
      my_spedigrees.push_back(&*sped);
    }
  }
}

/// Given supedigree, give the number of bits in the subpedigree
size_t PedigreeAnalysisSample::calculate_sped_bit_count(const FPED::Subpedigree& sp) const
{
  size_t founder_count    = 0;
  size_t nonfounder_count = 0;
  for(size_t i = 0; i < sp.member_count(); ++i)
  {
    if(sp.member_index(i).is_founder()) ++founder_count;
    else                                ++nonfounder_count;
  }
  
  return 2 * nonfounder_count - founder_count;
}

void PedigreeAnalysisSample::print_summary_statistics_table(ostream& o)
{
  OUTPUT::Table sum_tbl("Analysis Data Sample (summary)");
  
  sum_tbl << OUTPUT::TableColumn("");
  
  sum_tbl.beginColumnGroup("To Be Analyzed");
  
  sum_tbl << OUTPUT::TableColumn("");
  
  sum_tbl.beginColumnGroup("In Data Set");
  
  sum_tbl << OUTPUT::TableColumn("")
          << OUTPUT::TableColumn("")
          << OUTPUT::TableColumn("");
  
  sum_tbl <<
    (OUTPUT::TableRow() << "Number of Individuals:"
                        << my_fped.member_count()
                        << "("
                        << my_original_rped.member_count()
                        << ")"
    );
  
  sum_tbl <<
    (OUTPUT::TableRow() << "Number of Pedigrees:"
                        << my_fped.pedigree_count()
                        << "("
                        << my_original_rped.pedigree_count()
                        << ")"
    );

  sum_tbl <<
    (OUTPUT::TableRow() << "Number of Constituent Pedigrees:"
                        << my_filtered_sped_count
                        << "("
                        << my_original_sped_count
                        << ")"
    );

  o << sum_tbl;
}

void PedigreeAnalysisSample::print_member_table(ostream& o, bool detailed)
{
  // 1. Initialize the table and columns
  
  string title = "Individuals Removed from Analysis";
  if(detailed)
    title = "Individual Status for Analysis";
  
  OUTPUT::Table mem_table(title);

  mem_table   << OUTPUT::TableColumn("Ped")
              << OUTPUT::TableColumn("Ind");
            
  if(detailed)
    mem_table << OUTPUT::TableColumn("Status");
  
  mem_table   << OUTPUT::TableColumn("Reason");
  
  // 2. Print Rows

  // We only need mem_inf if we're printing detailed, but we want it out
  // of the loop.
  is_member_PAS_informative mem_inf(my_member_types);
  
  for(size_t i = 0; i < my_member_types.size(); ++i)
  {
    const RPED::Member& mem = my_original_rped.member_index(i);
    
    if(!detailed && mem_inf(mem)) continue;
    
    OUTPUT::TableRow current_row;
    
    current_row << mem.pedigree()->name() << mem.name();
    
    if(detailed)
    {
      if(mem_inf(mem)) current_row << "Kept     ";
      else             current_row << "Removed  ";
    }
    
    switch(my_member_types[i])
    {
      case MIT_SINGLETON                 : current_row << "Singleton";                         break;
      case MIT_DIRECTLY_INFORMATIVE      : current_row << "Informative (data)";                break;
      case MIT_STRUCTURALLY_INFORMATIVE  : current_row << "Informative (structure)";           break;
      case MIT_UNINFORMATIVE             : current_row << "Uninformative";                     break;
      case MIT_SINGLETON_AFTER           : current_row << "Singleton (after initial removal)"; break;
      case MIT_SPED_TOO_LARGE            : current_row << "Constituent Pedigree too large";    break;
      case MIT_SPED_TOO_SMALL            : current_row << "Constituent Pedigree too small";    break;
      case MIT_UNDEFINED                 : 
      case MIT_ERROR                     :
        SAGE_internal_error(); break;
    }

    mem_table << current_row;
  }

  o << mem_table << flush;
}

}
}
