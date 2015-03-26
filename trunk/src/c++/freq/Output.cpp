#include "freq/Output.h"
#include "app/SAGEapp.h"

namespace SAGE {
namespace FREQ {

//================================
//
//  generateOutput(...)
//
//================================
void
Output::generateOutput(const Results & results, const FPED::FilteredMultipedigree & fp)
{ 
  generateAnalysisOutput   (results, fp, true);
  generateAnalysisOutput   (results, fp, false);
  generateLocusDescription (results, fp);
}

//================================
//
//  generateAnalysisOutput(...)
//  
//================================  
void
Output::generateAnalysisOutput(const Results & results, const FPED::FilteredMultipedigree & fp, bool include_detailed)
{ 
  OUTPUT::Section s("Freq Analysis");
  
  // Set up configuration table:

  OUTPUT::Table config_table("Configuration");
  
  config_table << (OUTPUT::TableRow() << "Output file:" << results.config.ofilename);

  OUTPUT::TableRow founder_weight_row;
  
  founder_weight_row << "Founder weight:";
  
  if(results.config.use_founder_weight)
  {
    founder_weight_row << results.config.founder_weight;
  }
  else
  {
    founder_weight_row << "Disabled";
  }
    
  config_table << founder_weight_row;
  
  config_table << (OUTPUT::TableRow() << "Max. likelihood estimation:" << (results.config.skip_mle ? "Disabled" : "Enabled"));

  config_table << (OUTPUT::TableRow() << "Inbreeding coefficient:" << (results.config.estimate_inbreeding ? "Estimated" : "Not estimated"));

  s << config_table;
  
  // Note about non-convergences:
  
  s << OUTPUT::NamedString("Note", "A '*' by an MLE column indicates that the corresponding maximization may not have converged. Please review the detailed output file for more information.");
  
  // Per-marker analysis:

  for(MarkerAnalyses::const_iterator analysis_itr = results.analyses.begin(); analysis_itr != results.analyses.end(); ++analysis_itr)
  {
    // Fetch the marker info:
    const RPED::RefMarkerInfo & minfo = fp.info().marker_info(analysis_itr->marker_id);
    
    // Add the composite allele freqs table:
    OUTPUT::Table allele_freqs_table("Allele frequencies '" + minfo.name() + "'");
    
    // Add columns:
    allele_freqs_table << OUTPUT::TableColumn("Allele") 
                       << OUTPUT::TableColumn("Founders only") 
                       << OUTPUT::TableColumn("Entire dataset");
                       
    std::string mle_name1 = std::string(analysis_itr->maxfun_results.getConverged() ? "" : "*") +
                            std::string(results.config.estimate_inbreeding ? "MLE (no inbreeding)" : "MLE");
                       
    allele_freqs_table << OUTPUT::TableColumn(mle_name1);

    if(!results.config.skip_mle && results.config.estimate_inbreeding)
    {
      std::ostringstream column_name;
      
      column_name << (analysis_itr->maxfun_results_inbreeding.getConverged() ? "" : "*")
                  << "MLE (inbreeding = "
                  << analysis_itr->maxfun_results_inbreeding.getParameterMgr().getParameter("Inbreeding coeff.", "Inbreeding coeff.").getFinalEstimate()
                  << ")";

      allele_freqs_table << OUTPUT::TableColumn(column_name.str());
    }
      
    // Loop across all alleles:
    for(MLOCUS::allele_iterator al = minfo.allele_begin(); al != minfo.allele_end(); ++al)
    {
      OUTPUT::TableRow r;
      std::string      name = al->name();
        
      r << name;
        
      if(minfo.codominant())
      {
        r << analysis_itr->allele_freqs.founder_freqs[al->id()] << analysis_itr->allele_freqs.all_freqs[al->id()];
      }
      else // marker non-codominant
      {
        r << "Not computed" << "Not computed";
      }
        
      if(results.config.skip_mle)
      {
        r << "Not computed";
      }
      else // Mle computed
      {
        r << analysis_itr->maxfun_results.getParameterMgr().getParameter("Alleles", name).getFinalEstimate();
          
        if(results.config.estimate_inbreeding)
        {
          r << analysis_itr->maxfun_results_inbreeding.getParameterMgr().getParameter("Alleles", name).getFinalEstimate();
        }
      }

      allele_freqs_table << r;
    }
      
    s << allele_freqs_table;
                    
    if(include_detailed && !results.config.skip_mle)
    {
      s << MAXFUN::OutputFormatter::convertEstimates(analysis_itr->maxfun_results, true, false);
      
      if(results.config.estimate_inbreeding)
      {
        s << MAXFUN::OutputFormatter::convertEstimates(analysis_itr->maxfun_results_inbreeding, true, false)
          << MAXFUN::JointTest(analysis_itr->maxfun_results, analysis_itr->maxfun_results_inbreeding).summarizeAsTable("Test of Marker Specific Inbreeding Coefficient");
      }
    }
  }

  // Send to file:
  std::ofstream ofile;
  
  ofile.open(std::string(results.config.ofilename + (include_detailed ? ".det" : ".sum")).c_str());
  
  ofile << APP::SAGEapp::getReleaseString() << s << std::endl;
}

//================================
//
//  generateLocusDescription(...)
//  
//================================  
void
Output::generateLocusDescription(const Results & results, const FPED::FilteredMultipedigree & fp)
{
  // Set up output file:
  
  std::ofstream o;
  
  o.open(std::string(results.config.ofilename + ".loc").c_str());

  // Loop across markers:
  
  for(MarkerAnalyses::const_iterator analysis_itr = results.analyses.begin(); analysis_itr != results.analyses.end(); ++analysis_itr)
  {
    // Fetch the marker info:
    const RPED::RefMarkerInfo & minfo = fp.info().marker_info(analysis_itr->marker_id);
    
    // If the allele freqs weren't counted (not codominant) and maximum likelihood was skipped, then skip this guy:
    if(!minfo.codominant() && results.config.skip_mle)
      continue;
      
    // Print marker name:
    o << minfo.gmodel().name() << std::endl;
  
    // Print allele frequencies:
    if(results.config.skip_mle) // If mle was skipped, use the initial counts:
    {
      for(MLOCUS::allele_iterator al = minfo.allele_begin(); al != minfo.allele_end(); ++al)
      {
        o << "  " << std::setw(4) << left << al->name() << " =  " 
          << std::fixed << std::setprecision(6)
          << analysis_itr->allele_freqs.all_freqs[al->id()] << std::endl;
      }
    }
    else // mle enabled
    {
      for(MAXFUN::ParameterConstIterator param_itr  = analysis_itr->maxfun_results.getParameterMgr().getParamBegin ("Alleles"); 
                                         param_itr != analysis_itr->maxfun_results.getParameterMgr().getParamEnd   ("Alleles"); ++param_itr)
      {
        o << "  " << std::setw(4) << left << param_itr->getName() << " =  " << std::fixed << std::setprecision(6) << param_itr->getFinalEstimate() << std::endl;
      }
    }

    o << ";" << std::endl;
    
    // Print phenotype info:

    if(!minfo.codominant())
    {
      // Loop across phenotypes:
      for(MLOCUS::penetrance_model::phenotype_const_iterator ph_itr = minfo.phenotype_begin(); ph_itr != minfo.phenotype_end(); ++ph_itr)
      {
        o << ph_itr->name() << " = { ";
      
        bool at_least_one_genotype_printed = false;

        // Loop across valid penetrances:
        for(MLOCUS::penetrance_model::unphased_penetrance_iterator pen_itr  = minfo.unphased_penetrance_begin (*ph_itr);
                                                                   pen_itr != minfo.unphased_penetrance_end   (*ph_itr); ++pen_itr)
        {
          if(at_least_one_genotype_printed)
          {
            o << ",";
          }
          else
          {
            at_least_one_genotype_printed = true;
          }
          
          o << pen_itr.unphased_geno().name();

        } // End genotype loop
      
        o << " }" << std::endl;

      } // End phenotype loop

    } // End if-marker-not-codominant

    o << ";" << endl;
    
  } // End of marker loop
}

 }// End namespace FREQ
} // End namespace SAGE
