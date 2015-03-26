#include "tdtex/tdtex.h"

namespace SAGE  {
namespace TDTEX {

MersenneTwister mt;

/* ************************************************************************** **
** ***********              Function Specification              ************* **
** ************************************************************************** **
   Standard constructor
----------------------------------------------------------------------------- */
AppData::AppData(const string& program_name, bool debug_flag) :
   APP::SAGE_Simple_Data(program_name, debug_flag)
{
   my_cmdline_rules.add_rule(
      APP::ArgumentRuleset::Rule(APP::PARAMETER_FILE, APP::ArgumentRuleset::ONE)
   );
   my_cmdline_rules.add_rule(
      APP::ArgumentRuleset::Rule(APP::PEDIGREE_FILE, APP::ArgumentRuleset::ONE)
   );
   return;
}


/* ************************************************************************** **
** ***********              Function Specification              ************* **
** ************************************************************************** **
   Load input file
----------------------------------------------------------------------------- */
void AppData::process_input(int argc, char** argv) {
  parse_cmdline(argc, argv);

  // Load parameter file & check for validity.
  //
  read_parameter_file(my_parsed_arguments.get_arguments(APP::PARAMETER_FILE)[0]);

  // Create RPED::RefMultiPedigree & read pedigree file & check for validity.
  // Check for the errors in pedigrees.
  // Sort pedigrees.
  //
  read_family_data_file(
    my_parsed_arguments.get_arguments(APP::PEDIGREE_FILE)[0],
    true,    // dump_trait
    true,    // dump_marker
    false,   // skip_traits
    false,   // skip_markers
    true     // dynamic_markers
  );


  // - Create traits specified by function blocks.  -djb 5/2/01
  //
  evaluate_functions();

  // Read pedinfo_analysis block.
  //
  read_analysis();
}



//====================================================================
//
// batch_steps(...)
//
//====================================================================
size_t
batch_steps(double p, size_t min_steps = 400, size_t max_steps = 1000000)
{
  size_t steps = min_steps;

  if(!finite(p) && p > 0)
  {
    steps = max_steps;
  }
  else if(p > 0)
  {
    steps = std::max<size_t>(steps, (size_t)((1 - p) / p / 10));
  }

  return std::min<size_t>(max_steps, steps);
}

//====================================================================
//
// Constructor
//
//====================================================================
TdtexApp::TdtexApp(int argc, char** argv)
  : APP::SAGEapp(APP::APP_TDTEX, true, argc, argv)
{
  LSFInit();
  sage_cerr << prefix("%%TDTEX-%P: ");
}

//====================================================================
//
//  run_single(...)
//
//====================================================================
void
TdtexApp::run_single(const Configuration & config, const RPED::RefMultiPedigree& mp, MPED::SexCode sex,  cerrorstream& err, OUTPUT::Section & section)
{
  // Create Sampler:
  Sampler sampler(config, mp, sex, err);

  // Get a (alphabetically) ordered matrix.
  TransmissionTable::ReorderedMatrix rmatrix = sampler.get_table()->generate_reordered_matrix();

  bool        uses_sib_pairs = sampler.get_configuration().get_max_sib_pairs () != 0,
              uses_children  = sampler.get_configuration().get_max_children  () != 0;
  std::string title          = sampler.get_table()->units() + " transmissions to affected ";

  if(uses_children && uses_sib_pairs) title += "children and sibling pairs";
  else if(uses_children)              title += "children";
  else if(uses_sib_pairs)             title += "sibling pairs";
  else                                title += "nobody";

  OUTPUT::Table config_table;

  config_table << (OUTPUT::TableRow() << "Transmission type:"       <<  sampler.get_table()->units())
               << (OUTPUT::TableRow() << "Marker:"                  <<  sampler.get_gmodel().name())
               << (OUTPUT::TableRow() << "Trait:"                   <<  sampler.get_multi_pedigree().info().trait_info(sampler.get_configuration().get_trait()).name())
               << (OUTPUT::TableRow() << "Parental trait:"          << (sampler.get_configuration().get_parent_trait() == (size_t)-1 ?
                                                                        "Not specified" :
                                                                        sampler.get_multi_pedigree().info().trait_info(sampler.get_configuration().get_parent_trait()).name()))
               << (OUTPUT::TableRow() << "Parental sex(es) scored:" << (MPED::is_male(sampler.get_sex())   ? "Paternal" :
                                                                       (MPED::is_female(sampler.get_sex()) ? "Maternal" : "Paternal/maternal")));

  // Max children per family:
  OUTPUT::TableRow row = (OUTPUT::TableRow() << "Max. children/family:");

  if(sampler.get_configuration().get_max_children() == (size_t)-1)
    row << "Unlimited";
  else
    row << sampler.get_configuration().get_max_children();

  config_table << row;

  // Max sib pairs per family:
  row = (OUTPUT::TableRow() << "Max. sib pairs/family: ");

  if(sampler.get_configuration().get_max_sib_pairs() == (size_t)-1)
    row << "Unlimited";
  else
    row << sampler.get_configuration().get_max_sib_pairs();

  config_table << row;

  // Any errors in the sample?
  if(sampler.getErrorCount())
  {
    std::ostringstream s;

    s << sampler.getErrorCount() << " error(s) were found in your sample!  See the Information Output File for a detailed description of each problem.";

    config_table << OUTPUT::NamedString("Note", s.str());
  }

  section << config_table;

  // Figure out which indices have valid values:
  std::vector<bool> indices_informative(rmatrix.matrix.rows(), false);
  size_t            informative_index_cnt = 0,
                    sample_size           = 0,
                    diag_count            = 0;

  for(size_t i = 0; i < rmatrix.matrix.rows(); ++i)
  {
    bool index_informative = false;

    for(size_t j = 0; j < rmatrix.matrix.cols(); ++j)
      index_informative |= (rmatrix.matrix(i, j) || rmatrix.matrix(j, i));

    indices_informative[i] = index_informative;

    if(indices_informative[i])
      ++informative_index_cnt;
  }

  if(informative_index_cnt > 0)
  {
    // Create table:
    OUTPUT::Table transmission_table; // (sampler.get_table()->units() + " transmissions");

    // Display note about removed rows/columns, if applicable:
    if(informative_index_cnt != rmatrix.matrix.rows())
    {
      std::ostringstream s;
      s << rmatrix.matrix.rows() - informative_index_cnt << " empty rows/columns not shown";

      transmission_table << OUTPUT::NamedString("Note", s.str());
    }

    // Add columns:
    transmission_table << OUTPUT::TableColumn("T/NT");

    for(size_t i = 0; i < rmatrix.matrix.cols(); ++i)
      if(indices_informative[i])
        transmission_table << OUTPUT::TableColumn(rmatrix.headings[i]);

    // Add rows:
    for(size_t i = 0; i < rmatrix.matrix.rows(); ++i)
    {
      // Skip uninformative rows:
      if(indices_informative[i] == false)
        continue;

      // Increment diag_countP
      diag_count += rmatrix.matrix(i, i);

      // Set up row:
      OUTPUT::TableRow row = (OUTPUT::TableRow() << rmatrix.headings[i]);

      // Add all the cells:
      for(size_t j = 0; j < rmatrix.matrix.cols(); ++j)
      {
        if(indices_informative[j])
        {
          row << rmatrix.matrix(i, j);
          sample_size += rmatrix.matrix(i, j);
        }
      }

      // Add the row:
      transmission_table << row;
    }

    // Add the table:
    section << transmission_table;
  }

  // Sample information:
  section << (OUTPUT::Table()
          <<    OUTPUT::TableColumn("Group")
          <<    OUTPUT::TableColumn("Informative count")
          <<    OUTPUT::TableColumn("Total count")
          <<    OUTPUT::TableColumn("% Informative")

          <<   (OUTPUT::TableRow() << "Pedigrees" << sampler.getInfPedigreeCount() << sampler.getPedigreeCount()
          <<     ((sampler.getInfPedigreeCount() && sampler.getPedigreeCount()) ? (double)sampler.getInfPedigreeCount() / (double)sampler.getPedigreeCount() : 0.0))

          <<   (OUTPUT::TableRow() << "Families" << sampler.getInfFamilyCount() << sampler.getFamilyCount()
          <<     ((sampler.getInfFamilyCount() && sampler.getFamilyCount()) ? (double)sampler.getInfFamilyCount() / (double)sampler.getFamilyCount() : 0.0))

          <<   (OUTPUT::TableRow() << "Affected children" << sampler.getInfChildCount() << sampler.getChildCount()
          <<     ((sampler.getInfChildCount() && sampler.getChildCount()) ? (double)sampler.getInfChildCount() / (double)sampler.getChildCount() : 0.0))

          <<   (OUTPUT::TableRow() << "Affected sib pairs" << sampler.getInfPairCount() << sampler.getPairCount()
          <<     ((sampler.getInfPairCount() && sampler.getPairCount()) ? (double)sampler.getInfPairCount() / (double)sampler.getPairCount() : 0.0))

          <<   (OUTPUT::TableRow() << "Sample size" << (sample_size - diag_count) << sample_size
          <<     (((sample_size - diag_count) && sample_size) ? ((double)sample_size - (double)diag_count) / (double)sample_size : 0.0)));

  // Exact test statistics:
  chi_squared xx1 = mcnemar_statistic                      (sampler.get_table()->get_counts(), false),
              xx2 = mcnemar_statistic                      (sampler.get_table()->get_counts(), true),
              xx3 = pearson_marginal_homogeneity_statistic (sampler.get_table()->get_counts());

  // Set up table:
  OUTPUT::Table exact_tests_table;

  exact_tests_table
        <<  OUTPUT::TableColumn("Exact test statistic")
        <<  OUTPUT::TableColumn("P-value")
        <<  OUTPUT::TableColumn("Std. err.")
        << (OUTPUT::TableRow() << "Exact McNemar test"
        <<   (config.get_skip_permutation_test() ? "(skipped)" : pval(McNemarExact::mcnemar_exact(sampler.get_table()->get_counts()), 13)));

  // MonteCarlo McNemar test:
  OUTPUT::TableRow mcmn_row = (OUTPUT::TableRow() << "Monte Carlo McNemar test");

  if(config.get_skip_mc_test())
  {
    mcmn_row << "(skipped)";
  }
  else
  {
    simulated_pvalue sp = McNemarMonteCarlo::mcnemar_monte_carlo(sampler.get_table()->get_counts(), 1000, batch_steps(xx2.pvalue()));

    mcmn_row << pval(sp.pvalue(), 13) << fp(sp.standard_error(), 10, 8);
  }

  exact_tests_table << mcmn_row;

  // Monte Carlo Marginal Homogeneity test:

  OUTPUT::TableRow mcmh_row = (OUTPUT::TableRow() << "Monte Carlo Marginal Homogeneity");

  if(config.get_skip_mcmh_test())
  {
    mcmh_row << "(skipped)";
  }
  else
  {
    simulated_pvalue sp = MargHomoMonteCarlo::mh_monte_carlo(sampler.get_table()->get_counts(), 1000, batch_steps(xx3.pvalue()));

    mcmh_row << pval(sp.pvalue(), 13) << fp(sp.standard_error(), 10, 8);
  }

  exact_tests_table << mcmh_row;

  section << exact_tests_table;

  // Asymptotic tests:
  section << (OUTPUT::Table()
          <<  OUTPUT::TableColumn("Asymptotic test statistic")
          <<  OUTPUT::TableColumn("P-value")
          << (OUTPUT::TableRow() << "McNemar test"                      << pval(xx1.pvalue(), 13))
          << (OUTPUT::TableRow() << "Continuity corrected McNemar test" << pval(xx2.pvalue(), 13))
          << (OUTPUT::TableRow() << "Marginal homogeneity test"         << pval(xx3.pvalue(), 13)));
}

//====================================================================
//
//  runOneMarkerOneTrait(...)
//
//====================================================================
void
TdtexApp::runOneMarkerOneTrait(const Configuration & config, const RPED::RefMultiPedigree& mped,  cerrorstream& err, OUTPUT::Section & file_specific_section)
{
  // Figure out name for analysis section:
  std::ostringstream name;

  name << "Analysis #" << (file_specific_section.getVector().count<OUTPUT::Section>() + 1);

  // Create analysis section with this name:
  OUTPUT::Section analysis_section(name.str());

  // If sex differential is enabled and we are analyzing alleles:
  if(config.get_sex_differential() && (config.get_method() == Configuration::ALLELES))
  {
    run_single(config, mped, MPED::SEX_MISSING, err, analysis_section);
    run_single(config, mped, MPED::SEX_MALE,    err, analysis_section);
    run_single(config, mped, MPED::SEX_FEMALE,  err, analysis_section);
  }
  // Otherwise, just run the simple analysis:
  else
  {
    run_single(config, mped, MPED::SEX_MISSING, err, analysis_section);
  }

  // Add analysis section to file_specific_section:
  file_specific_section << analysis_section;
}

//====================================================================
//
//  runAllMarkersAllTraits(...)
//
//====================================================================
void
TdtexApp::runAllMarkersAllTraits(
        Configuration       & config,
  const RPED::MultiPedigree & mped,
        cerrorstream        & err,
        OUTPUT::Section     & file_specific_section)
{
  for(size_t t = 0; t < mped.info().trait_count(); ++t)
  {
    if(    mped.info().trait_info(t).name() == "SEX_CODE"
        || mped.info().trait_info(t).name() == "FAMILIAL_INDICATOR"
        || mped.info().trait_info(t).name() == "FOUNDER_INDICATOR"
        || mped.info().trait_info(t).name() == "PEDIGREE_SIZE" )
      continue;

    config.set_trait(t);

    for(size_t m = 0; m < mped.info().marker_count(); ++m)
    {
      config.set_marker(m);

      runOneMarkerOneTrait(config, mped, err, file_specific_section);
    }
  }
}

//====================================================================
//
// main (TdtexApp)
//
//====================================================================
int
TdtexApp::main()
{
  // If the number of arguments is invalid:
  if((argc != 3) && (argc != 5)) {
     exit(EXIT_SUCCESS);
  }

  // Create and set up data object:
  AppData tdtdata(name, debug());
  print_title(tdtdata.info());
  tdtdata.process_input(argc, argv);

  // Make sure there is at least one trait and at least one marker

  if(!tdtdata.pedigrees().info().trait_count())
  {
    tdtdata.errors() << priority(critical) << "No traits to analyze... aborting." << std::endl;
    exit(EXIT_FAILURE);
  }

  if(!tdtdata.pedigrees().info().marker_count())
  {
    tdtdata.errors() << priority(critical) << "No markers to analyze... aborting." << std::endl;
    exit(EXIT_FAILURE);
  }

  // Runtime output:
  std::cout << "Performing TDT analysis .................." << flush;

  // Each entry indicates an output filename (eg: "tdtex1.out") that has been opened/closed (made empty).
  typedef std::map<std::string, OUTPUT::Section> SectionMap;
  SectionMap sections;

  for(LSFList::iterator i = tdtdata.parameters()->List()->begin(); i != tdtdata.parameters()->List()->end(); ++i)
  {
    if(!*i)
      continue;

    // Is this a TDT block?
    std::string name = toUpper((*i)->name());

    if(name == "TDT_ANALYSIS" || name == "TDT" || name == "TDTEX")
    {
      // Set up a Configuration:
      Configuration config;

      // Try reading in the configuration:
      try { config = Parser::parse_parameters(tdtdata.pedigrees().info(), *i); }

      catch(const std::exception & e) { continue; }

      // Figure out which analyses to run, and run them:
      bool all_markers = config.get_marker () == (size_t)-1,
           all_traits  = config.get_trait  () == (size_t)-1;

      if(all_markers)
      {
        if(all_traits) // All markers, all traits
        {
          runAllMarkersAllTraits(config, tdtdata.pedigrees(), tdtdata.errors(), sections[config.get_ofilename()]);
        }
        else // All markers, one trait
        {
          for(size_t m = 0; m < tdtdata.pedigrees().info().marker_count(); ++m)
          {
            config.set_marker(m);

            runOneMarkerOneTrait(config, tdtdata.pedigrees(), tdtdata.errors(), sections[config.get_ofilename()]);
          }
        }
      }
      else // one marker...
      {
        if(all_traits) // One marker, all traits...
        {
          for(size_t t = 0; t < tdtdata.pedigrees().info().trait_count(); ++t)
          {
            if(    tdtdata.pedigrees().info().trait_info(t).name() == "SEX_CODE"
                || tdtdata.pedigrees().info().trait_info(t).name() == "FAMILIAL_INDICATOR"
                || tdtdata.pedigrees().info().trait_info(t).name() == "FOUNDER_INDICATOR"
                || tdtdata.pedigrees().info().trait_info(t).name() == "PEDIGREE_SIZE" )
              continue;

            config.set_trait(t);

            runOneMarkerOneTrait(config, tdtdata.pedigrees(), tdtdata.errors(), sections[config.get_ofilename()]);
          }
        }
        else // One marker, one trait...
        {
          runOneMarkerOneTrait(config, tdtdata.pedigrees(), tdtdata.errors(), sections[config.get_ofilename()]);
        }
      }

    } // End if-this-is-an-analysis-block

  } // Done with analysis block loop

  // If there were no analysis blocks encountered, then just run everything!
  if(sections.size() == 0)
  {
    Configuration config;

    runAllMarkersAllTraits(config, tdtdata.pedigrees(), tdtdata.errors(), sections[config.get_ofilename()]);
  }

  // Generate output files:

  for(SectionMap::iterator s = sections.begin(); s != sections.end(); ++s)
  {
    std::ofstream ofile(s->first.c_str());

    ofile << getReleaseString();

    ofile << s->second << std::flush;
  }

  // Runtime output:
  std::cout << "done." << std::endl;

  print_inf_banner(cout);

  // Success!
  return EXIT_SUCCESS;
}

} // End namespace TDTEX
} // End namespace SAGE

//====================================================================
//
//  main(...)
//
//====================================================================
int
main(int argc, char* argv[])
{
  return SAGE::TDTEX::TdtexApp(argc, argv).main();
}

