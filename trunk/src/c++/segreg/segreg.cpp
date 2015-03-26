#include <string>
#include <functional>
#include <fstream>
#include <iomanip>
#include <cassert>
#include "LSF/LSFinit.h"
#include "LSF/LSFfile.h"
#include "app/SAGEapp.h"
#include "rped/rpfile.h"
#include "rped/rped.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "segreg/segreg.h"
#include "segreg/segreg_input.h"
#include "segreg/AnalysisController.h"
#include "segreg/iterative_models.h"

using namespace std;
namespace SAGE
{
namespace SEGREG
{


segreg::segreg(int argc, char** argv)
     : APP::SAGEapp(APP::APP_SEGREG, true, argc, argv)
{
  LSFInit();
}

int segreg::main()
{
  segreg_data data(name, debug());

  print_title(data.info());

  cerrorstream& errors = data.errors();

  data.input(argc, argv);

  if(!data.analyses().size())
  {
    errors << priority(fatal) << "No valid analyses specified.  Program cannot continue." 
           << endl;

    exit(EXIT_FAILURE);
  }
  else if(!evaluate_pedigree_validity(data))
  {
    errors << priority(fatal) << "No valid pedigrees specified. Program cannot continue."
           << endl;

    exit(EXIT_FAILURE);
  }
  else
  {
    const vector<model>& analyses = data.analyses();
    
   
       if (data.locus_indic_vec.size() > 0) { // due to JA
       if (analyses.size() != data.locus_indic_vec.size() ){ 
         errors << priority(fatal) << "Miscounting of models ! "
           << endl;
          exit(EXIT_FAILURE);
        }
       } // making sure fpmm models are read in synch by the parser
     
    unsigned ii = 0; // due to JA
      
    for(vector<model>::const_iterator i = analyses.begin(); i != analyses.end(); ++i)
    {
      
      AnalysisController controller(data.get_ostreams()); 
      controller.like_cutoff = data.like_cutoff;
        if (data.locus_indic_vec[ii] == 1) {
            controller.fpmm_and_no_poly_loci = true; }
        else {controller.fpmm_and_no_poly_loci = false;}
       controller.like_cutoff = data.like_cutoff;
       controller.do_analysis(data.pedigrees(), *i);
       ii++;
       Iterative_Models::intermax.clear(); // due to JA needed to avoid seg_fault
    }
  }

  print_inf_banner(cout);

  return EXIT_SUCCESS;
}

//----------------------------------------------------------------------------------
//        
// evaluate_pedigree_validity(...)
//
//----------------------------------------------------------------------------------
inline
bool segreg::evaluate_pedigree_validity(const segreg_data& data) const
{
   const RPED::RefMultiPedigree& RMP = data.pedigrees();

   int valid_subpedigree_count       = 0;
   int total_valid_subpedigree_count = 0;
   int loop_count                    = 0;
   int total_loop_count              = 0;
        
   data.screen() << "Evaluating pedigrees........." << flush;
          
   bool valid = true;

   // Determine if there are loops.  Pedigrees without loops cannot be processed
   for(RPED::PedigreeConstIterator current_pedigree  = RMP.pedigree_begin ();
       valid && current_pedigree != RMP.pedigree_end();
       ++current_pedigree)
   {
     RPED::RefMultiPedigree::subpedigree_const_iterator sub_itr  = current_pedigree->subpedigree_begin ();
     RPED::RefMultiPedigree::subpedigree_const_iterator sub_end  = current_pedigree->subpedigree_end   ();

     for( ; valid && sub_itr != sub_end; ++sub_itr)
     {
       if(MPED::mp_utilities::has_loops(*sub_itr))
       {
         valid = false;
       }
     }
   }

   // If the valid flag is still set, there are no pedigree problems.
   if(valid)
   {
     data.screen() << "...Done." << endl << endl;

     return RMP.pedigree_count() > 0; 
   }

   // Otherwise, we must report the problems
   data.screen() << endl << endl;
 
   data.info() << "SEGREG has detected the following pedigree problems: "
               << endl << endl;

   data.info() << "         Connected pedigree sections" << endl
               << endl
               << "Pedigree      Total      Valid    w. Loops"  << endl
               << "----------  ---------  ---------  ---------" << endl;

   for(RPED::RefMultiPedigree::pedigree_const_iterator
           current_pedigree  = RMP.pedigree_begin ();
           current_pedigree != RMP.pedigree_end   (); ++current_pedigree)
   {
     valid_subpedigree_count = 0;
     loop_count              = 0;
          
     for(RPED::RefMultiPedigree::subpedigree_const_iterator sub_itr  = current_pedigree->subpedigree_begin ();
                                    sub_itr != current_pedigree->subpedigree_end   (); sub_itr++)
     {
       if(MPED::mp_utilities::has_loops(*sub_itr))
       {
         ++loop_count;
       }

       if(MPED::mp_utilities::no_loops(*sub_itr))
       {
         ++valid_subpedigree_count;
       }
     }
         
     total_loop_count              += loop_count;
     total_valid_subpedigree_count += valid_subpedigree_count;

     // If there are loops, we must report the problems.
     if(loop_count)
     {
         data.info() << setw(10) << current_pedigree->name()              << "  "
                     << setiosflags(ios::right)
                     << setw(9)  << current_pedigree->subpedigree_count() << "  "
                     << setw(9)  << valid_subpedigree_count               << "  "
                     << setw(9)  << loop_count << endl;
     }
   }
       
   cerrorstream screen(data.screen());

   screen << "There are a total of " << total_valid_subpedigree_count
          << " valid connected pedigree(s) in this analysis.";

   screen << "  There are " << total_loop_count << " connected pedigrees "
          << "that will be ignored due to loops.";

   screen << "  See .inf file for detailed information regarding pedigrees being ignored."
          << endl;

   data.screen() << endl << endl;

   data.info() << endl 
               << "These connected pedigrees will be ignored for analyses."
               << endl << endl;

   if(total_valid_subpedigree_count) return true;
   else                              return false;
}         

}
}

int main(int argc, char* argv[])
{
  free(malloc(1));

  SAGE::SEGREG::segreg segreg_inst(argc, argv);

  segreg_inst.main();

  exit(EXIT_SUCCESS);
}

