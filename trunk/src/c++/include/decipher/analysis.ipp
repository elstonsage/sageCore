//============================================================================
// File:      analysis.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   4/8/4 created        -djb
//                                                                          
// Notes:     Inline implementation of analysis.
//                                                                          
// Copyright (c) 2004 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


//============================================================================
// IMPLEMENTATION:  analysis
//============================================================================
//
inline
analysis::analysis(APP::Output_Streams& streams, ostream& detail_file, ostream& summary_file,
                   ostream& dump_file, 
                   const RefMultiPedigree& mped, const genome_description* gen, 
                   const instructions& instr)
      : my_streams(streams), my_errors(streams.errors()), my_messages(streams.messages()), 
        my_detail_file(detail_file), my_summary_file(summary_file), my_dump_file(dump_file), 
        my_mped(mped), my_genome_description(gen), my_instructions(instr), my_filtered_mped(mped),
        my_partitioner(instr, my_filtered_mped, my_errors)
{
  build_filtered_mped();
}

inline void
analysis::report_inconsistencies(const family_generator::iterator& iter)
{
  const vector<size_t>&  inconsistencies = iter.get_inconsistent_loci();

  vector<size_t>::const_iterator  locus_iter = inconsistencies.begin();
  vector<size_t>::const_iterator  locus_end_iter = inconsistencies.end();
  for(; locus_iter != locus_end_iter; ++locus_iter)
  {
    my_errors << priority(warning) << "In pedigree " << iter.get_pedigree()->name() 
              << ", constituent pedigree containing member "
              << iter.get_member()->name() << " contains a Mendelian inconsistency "
              << "at marker " << my_mped.info().marker_info(*locus_iter).name() << ".  "
              << "Genotypes treated as missing at this marker for all members of "
              << "this constituent pedigree. " << endl;
  } 
}

inline const string&
analysis::title() const
{
  return  my_instructions.title;
}




