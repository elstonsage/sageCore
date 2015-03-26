//============================================================================
// File:      tasks.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   10/18/2 - created.                                djb
//                                                                          
// Notes:     
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


// - This function called by all hypothesis test calculate() functions
//   with the appropriate parameters.  T is test type.  R is result type.
//
template<class T, class R> void
do_task_calculations(T& task)
{
  if(task.my_instructions.valid)
  {
    const FPED::FilteredMultipedigreeInfo&  mp_info = task.my_mped.info();
    size_t  trait_index = mp_info.marker_find(task.my_instructions.trait);
    
    for(size_t marker_index = 0; marker_index < mp_info.marker_count(); ++marker_index)
    {
      if(marker_index == trait_index)
      {
        continue;
      }
      
      R*  result = new R;
      result->trait = task.my_instructions.trait;
      result->marker = mp_info.marker_info(marker_index).name();
      
      task.calculate_alt(trait_index, marker_index, *result);
      task.calculate_null(trait_index, marker_index, *result);      
      
      task.my_results.push_back(result_ptr(result)); 
    }
    
    task.completed = true;
  }
}



//============================================================================
// IMPLEMENTATION:  task
//============================================================================
//
inline
task::task(cerrorstream& errors, const FPED::FilteredMultipedigree& mped, const instructions& instr)
      : my_errors(errors), completed(false), my_mped(mped), my_instructions(instr)
{}

inline
task::~task()
{}

inline void
task::write(ostream& summary, ostream& detail) const
{
  assert(completed);

  write_summary(summary);
  write_detail(detail);
}


//============================================================================
// IMPLEMENTATION:  non_ss_smiths_faraways_test
//============================================================================
//
inline
non_ss_smiths_faraways_test::non_ss_smiths_faraways_test(cerrorstream& errors, const FPED::FilteredMultipedigree& mped, 
                                                         const instructions& instr,
                                                         sf_type t)
      : task(errors, mped, instr), my_type(t)
{}

inline void  
non_ss_smiths_faraways_test::announce_start() const
{
  cout << "Performing " << (my_type == sf_LINKAGE ? "test for linkage under Smith's model ..." :
                                                    "Smith's homogeneity test ..."    )
                        << endl;
}

//============================================================================
// IMPLEMENTATION:  ss_smiths_faraways_test
//============================================================================
//
inline
ss_smiths_faraways_test::ss_smiths_faraways_test(cerrorstream& errors, const FPED::FilteredMultipedigree& mped, 
                                                 const instructions& instr,
                                                 sf_type t)
      : task(errors, mped, instr), my_type(t)
{}

inline void  
ss_smiths_faraways_test::announce_start() const
{
  cout << "Performing " << (my_type == sf_LINKAGE ? "test for linkage under Smith's model ..." :
                                                    "Smith's homogeneity test ..."    )
                        << endl;
}


