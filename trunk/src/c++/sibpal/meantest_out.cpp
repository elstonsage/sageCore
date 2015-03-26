#include "sibpal/meantest_out.h"

using namespace std;

namespace SAGE   {
namespace SIBPAL {

void
meantest_textfile::print_results(const meantest_parameters& reg)
{
  ostream& out = output_stream();
  const relative_pairs& pairs = sib_pairs();

  if(!out)
    return;

  if( wide_output() )
    out << "====================================================================================" << endl;
  else
    out << "======================================================================" << endl;
  out << "  Test of Mean Allele Sharing IBD";

  if( reg.get_use_full_sibs() && !(reg.get_use_half_sibs()) )
    out << " for Full Sib Pairs" << endl;
  else if( reg.get_use_half_sibs() && !(reg.get_use_full_sibs()) )
    out << " for Half Sib Pairs" << endl;

  if( wide_output() )
    out << "====================================================================================" << endl;
  else
    out << "======================================================================" << endl;
  out << endl
      << "  Estimates:"                                                      << endl
      << "    pi - Average proportion of alleles shared IBD."                << endl;
  if( reg.get_use_full_sibs() && !(reg.get_use_half_sibs()) )
    out << "    fi - Estimated proportion of sib pairs sharing i alleles IBD." << endl;
  out << endl;

  if( reg.subset_count() )
  {
    out << "             Subset: ";

    for(size_t i=0; i < reg.subset_count(); ++i)
    {
       if(i)
         out << ", ";
       out << reg.subsets(i).name(pairs);
    }
    out << endl << endl;
  }

  if( reg.get_use_full_sibs() && !(reg.get_use_half_sibs()) )
    out << "  w1: " << fp(reg.get_w(), 4) << endl << endl;

  if( wide_output() )
    out << "====================================================================================" << endl
        << "Marker            Pairs        Estimate      Std Error   T-value       P-value      " << endl;
  else
    out << "======================================================================" << endl
        << "Marker            Pairs        Estimate      Std Error   P-value      " << endl;

  for(size_t i=0; i < reg.marker_count(); ++i)
    print_results(reg, i);

  out << "======================================================================";
  if( wide_output() )
    out << "==============";
  out << endl << endl << endl;
}

void
meantest_textfile::print_results(const meantest_parameters& reg, size_t i)
{
  ostream& out = output_stream();
  const relative_pairs& pairs = sib_pairs();

  if(!out)
    return;

  out << "----------------------------------------------------------------------";
  if( wide_output() )
    out << "--------------";
  out << endl;

  const marker_parameter& param = reg[i];

  string marker = pairs.marker_name(param.marker);
  if(!marker.size())
    marker = "Marker #" + long2str(param.marker);

  int size = marker.size();
  out << marker << setw(18-size) << "" << setw(5) << param.pair_count;
  int  w=12;
  int pw=16;
  int pr=8;
  long df = param.pair_count - 1;

  double pi_mean = 0.25 + reg.get_w() * 0.5;
  if( reg.get_use_half_sibs() && !(reg.get_use_full_sibs()) )
    pi_mean = 0.25;  

  out << "  pi  ";
  print_result_line(param.estimate.pi(), param.estimate.standard_error(),
                    param.t_value(pi_mean), param.p_value(df, pi_mean),
                    w, pr, pw, reg.get_pvalues_scientific_notation(), false); 

  if( print_only_pi() )
    return;

  if( reg.get_use_half_sibs() && !(reg.get_use_full_sibs()) )
    return;

  for( size_t j = 0; j < 3; ++j )
  {
    double fi_mean = 0.5;

    if( reg.get_use_full_sibs() && !(reg.get_use_half_sibs()) )
    {
      switch( j )
      {
        case 0: fi_mean = 0.25; break;
        case 1: fi_mean = 0.5;  break;
        case 2: fi_mean = 0.25; break;
      }

      out << setw(23) << ""
          << "  f" << j << "  ";

      if( j == 1 )
      {
        print_result_line(param.estimate.f(j), param.estimate.standard_error(j),
                          numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN(),
                          w, pr, pw, reg.get_pvalues_scientific_notation(), false); 
      }
      else
      {
        print_result_line(param.estimate.f(j), param.estimate.standard_error(j),
                          param.t_value(j, fi_mean), param.p_value(j, df, fi_mean),
                          w, pr, pw, reg.get_pvalues_scientific_notation(), false); 
      }
    }
    else if( reg.get_use_half_sibs() && !(reg.get_use_full_sibs()) )
    {
      if( j > 1 )
        continue;

      out << setw(23) << ""
          << "  f" << j << "  ";

      //print_result_line(param.estimate.f(j), param.estimate.standard_error(j),
      //                  param.t_value(j, fi_mean), param.p_value(j, df, fi_mean),
      //                  w, pr, pw, reg.get_pvalues_scientific_notation()); 
      print_result_line(param.estimate.f(j), param.estimate.standard_error(j),
                        param.t_value(j, fi_mean), numeric_limits<double>::quiet_NaN(),
                        w, pr, pw, reg.get_pvalues_scientific_notation(), false); 
    }
  }
}

void
meantest_textfile::print_results(const meantest_param_vector& reg, const string& trait)
{
  ostream& out = output_stream();
  const relative_pairs& pairs = sib_pairs();

  if( wide_output() )
    out << "=====================================================================================================" << endl;
  else
    out << "=======================================================================================" << endl;
  out << "  Test of Mean Allele Sharing IBD";

  if( reg[0].get_use_full_sibs() && !(reg[0].get_use_half_sibs()) )
    out << " for Full Sib Pairs" << endl;
  else if( reg[0].get_use_half_sibs() && !(reg[0].get_use_full_sibs()) )
    out << " for Half Sib Pairs" << endl;

  out << "    Conditional on Binary Trait '" 
      << reg[0].trait().name(pairs) << "'."
      << endl;
  if( wide_output() )
    out << "=====================================================================================================" << endl;
  else
    out << "=======================================================================================" << endl;
  out << endl;

  if(!out || reg.size() != 3)
  {
    out << "Cannot print regression results.  Bad size == " << reg.size() << endl;
    return;
  }

  if( reg[0].get_use_full_sibs() && !(reg[0].get_use_half_sibs()) )
    out << "  Estimates:                                    #Aff:"                          << endl
        << "    pi - Mean proportion of alleles shared IBD.   0 - Concordantly unaffected." << endl
        << "    fi - Proportion of sib pairs sharing i        1 - Discordant."   << endl
        << "         alleles IBD.                             2 - Concordantly affected."   << endl
        << endl;
  else
    out << "  Estimates:                                    #Aff:"                          << endl
        << "    pi - Mean proportion of alleles shared IBD.   0 - Concordantly unaffected." << endl
        << "                                                  1 - Discordant."   << endl
        << "                                                  2 - Concordantly affected."   << endl
        << endl;

  out << "  Exp LINK: Value of estimate expected if there is linkage." << endl
      << endl;

  out << "  Trait: " << reg[0].trait().name(pairs) << endl;

  if( reg[0].subset_count() )
  {
    out << "  Subset: ";

    for(size_t i=0; i < reg[0].subset_count(); ++i)
    {
       if(i)
         out << ", ";
       out << reg[0].subsets(i).name(pairs);
    }
    out << endl;
  }
  out << endl;

  if( reg[0].get_use_full_sibs() && !(reg[0].get_use_half_sibs()) )
    out << "  w1: " << fp(reg[0].get_w(), 4) << endl << endl;

  if( wide_output() )
    out << "=====================================================================================================" << endl
        << "Marker                #Aff  Pairs    Estimate      Std Error   T-value       P-value         Exp LINK" << endl;
  else
    out << "=======================================================================================" << endl
        << "Marker                #Aff  Pairs    Estimate      Std Error   P-value         Exp LINK" << endl;

  for(size_t i=0; i < reg[0].marker_count(); ++i)
    print_results(reg, i, trait);

  if( wide_output() )
    out << "=====================================================================================================" << endl;
  else
    out << "=======================================================================================" << endl;
  out << endl << endl << endl;
}

void
meantest_textfile::print_results(const meantest_param_vector& reg, size_t i, const string& trait)
{
  ostream& out = output_stream();
  const relative_pairs& pairs = sib_pairs();

  if(!out)
    return;

  if( wide_output() )
    out << "-----------------------------------------------------------------------------------------------------" << endl;
  else
    out << "---------------------------------------------------------------------------------------" << endl;

  const marker_parameter& param0 = reg[0][i];
  const marker_parameter& param1 = reg[1][i];
  const marker_parameter& param2 = reg[2][i];

  string marker = pairs.marker_name(param0.marker);
  if(!marker.size())
    marker = "Marker #" + long2str(param0.marker);

  int size = marker.size();
  out << marker;
  out << setw(15-size) << "";

  int  w=12;
  int pw=16;
  int pr=8;

  long df0 = param0.pair_count - 1;
  long df1 = param1.pair_count - 1;
  long df2 = param2.pair_count - 1;

  double pi_mean = 0.25 + reg[0].get_w() * 0.5;
  if( reg[0].get_use_half_sibs() && !(reg[0].get_use_full_sibs()) )
    pi_mean = 0.25;  

  // Print pi 0 line.
  out << "  pi     " << param0.affecteds            << "   " 
      << setw(5) << param0.pair_count               << "  ";

  print_result_line(param0.estimate.pi(), param0.estimate.standard_error(),
                    param0.t_value(pi_mean), param0.p_value(df0, pi_mean),
                    w, pr, pw, reg[0].get_pvalues_scientific_notation(), true, pi_mean, false); 

  // Print pi 1 line.
  out << setw(24) << "" << param1.affecteds         << "   "   
      << setw(5)  << param1.pair_count              << "  ";

  print_result_line(param1.estimate.pi(), param1.estimate.standard_error(),
                    param1.t_value(pi_mean), param1.p_value(df1, pi_mean),
                    w, pr, pw, reg[1].get_pvalues_scientific_notation(), true, pi_mean, true); 

  // Print pi 2 line.
  out << setw(24) << "" << param2.affecteds         << "   "   
      << setw(5)  << param2.pair_count              << "  ";

  print_result_line(param2.estimate.pi(), param2.estimate.standard_error(),
                    param2.t_value(pi_mean), param2.p_value(df2, pi_mean),
                    w, pr, pw, reg[2].get_pvalues_scientific_notation(), true, pi_mean, false); 

  if( print_only_pi() )
    return;

  if( reg[0].get_use_half_sibs() && !(reg[0].get_use_full_sibs()) )
    return;

  for( size_t j = 0; j < 3; ++j )
  {
    double fi_mean = 0.5;

    if( reg[j].get_use_full_sibs() && !(reg[j].get_use_half_sibs()) )
    {
      switch( j )
      {
        case 0: fi_mean = 0.25; break;
        case 1: fi_mean = 0.5;  break;
        case 2: fi_mean = 0.25; break;
      }

      if( wide_output() )
        out << setw(18) << "" << "    ----  -----  ------------  ------------  ------------  -------------    -------" << endl;
      else
        out << setw(18) << "" << "    ----  -----  ------------  ------------  -------------    -------" << endl;

      // param0 - 0 aff
      out << setw(15) << "" << "  f" << j << "     " << param0.affecteds << "   "
          << setw(5)  << param0.pair_count            << "  ";

      if( j == 0 )
      {
        print_result_line(param0.estimate.f(j), param0.estimate.standard_error(j),
                          param0.t_value(j, fi_mean), param0.p_value(j, df0, fi_mean),
                          w, pr, pw, reg[0].get_pvalues_scientific_notation(), true, fi_mean, true); 
      }
      else if( j == 1 )
      {

        print_result_line(param0.estimate.f(j), param0.estimate.standard_error(j),
                          numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN(),
                          w, pr, pw, reg[0].get_pvalues_scientific_notation(), false); 
      }
      else
      {
        print_result_line(param0.estimate.f(j), param0.estimate.standard_error(j),
                          param0.t_value(j, fi_mean), param0.p_value(j, df0, fi_mean),
                          w, pr, pw, reg[0].get_pvalues_scientific_notation(), true, fi_mean, false); 
      }

      // param1
      out << setw(24) << "" << param1.affecteds       << "   "
          << setw(5)  << param1.pair_count            << "  ";

      if( j == 0 )
      {
        print_result_line(param1.estimate.f(j), param1.estimate.standard_error(j),
                          param1.t_value(j, fi_mean), param1.p_value(j, df1, fi_mean),
                          w, pr, pw, reg[0].get_pvalues_scientific_notation(), true, fi_mean, false); 
      }
      else if( j == 1 )
      {
        print_result_line(param1.estimate.f(j), param1.estimate.standard_error(j),
                          numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN(),
                          w, pr, pw, reg[1].get_pvalues_scientific_notation(), false); 
      }
      else
      {
        print_result_line(param1.estimate.f(j), param1.estimate.standard_error(j),
                          param1.t_value(j, fi_mean), param1.p_value(j, df1, fi_mean),
                          w, pr, pw, reg[0].get_pvalues_scientific_notation(), true, fi_mean, true);
      }

      // param2
      out << setw(24) << "" << param2.affecteds       << "   "
          << setw(5)  << param2.pair_count            << "  ";

      if( j == 0 )
      {
        print_result_line(param2.estimate.f(j), param2.estimate.standard_error(j),
                          param2.t_value(j, fi_mean), param2.p_value(j, df2, fi_mean),
                          w, pr, pw, reg[2].get_pvalues_scientific_notation(), true, fi_mean, true);
      }
      else if( j == 1 )
      {
        print_result_line(param2.estimate.f(j), param2.estimate.standard_error(j),
                          numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN(),
                          w, pr, pw, reg[2].get_pvalues_scientific_notation(), false); 
      }
      else
      {
        print_result_line(param2.estimate.f(j), param2.estimate.standard_error(j),
                          param2.t_value(j, fi_mean), param2.p_value(j, df2, fi_mean),
                          w, pr, pw, reg[2].get_pvalues_scientific_notation(), true, fi_mean, false);
      }
    }
    else if( reg[j].get_use_half_sibs() && !(reg[j].get_use_full_sibs()) )
    {
      if( j > 1 )
        continue;

      if( wide_output() )
        out << setw(18) << "" << "    ----  -----  ------------  ------------  ------------  -------------    -------" << endl;
      else
        out << setw(18) << "" << "    ----  -----  ------------  ------------  -------------    -------" << endl;

      // param0 - 0 aff
      out << setw(15) << "" << "  f" << j << "     " << param0.affecteds << "   "
          << setw(5)  << param0.pair_count            << "  ";


      if( j == 0 )
        print_result_line(param0.estimate.f(j), param0.estimate.standard_error(j),
                          param0.t_value(j, fi_mean), param0.p_value(j, df0, fi_mean),
                          w, pr, pw, reg[0].get_pvalues_scientific_notation(), true, fi_mean, true);
      else
        print_result_line(param0.estimate.f(j), param0.estimate.standard_error(j),
                          param0.t_value(j, fi_mean), param0.p_value(j, df0, fi_mean),
                          w, pr, pw, reg[0].get_pvalues_scientific_notation(), true, fi_mean, false);

      // param1
      out << setw(24) << "" << param1.affecteds       << "   "
          << setw(5)  << param1.pair_count            << "  ";

      if( j == 0 )
        print_result_line(param1.estimate.f(j), param1.estimate.standard_error(j),
                          param1.t_value(j, fi_mean), param1.p_value(j, df1, fi_mean),
                          w, pr, pw, reg[1].get_pvalues_scientific_notation(), true, fi_mean, false); 
      else
        print_result_line(param1.estimate.f(j), param1.estimate.standard_error(j),
                          param1.t_value(j, fi_mean), param1.p_value(j, df1, fi_mean),
                          w, pr, pw, reg[1].get_pvalues_scientific_notation(), true, fi_mean, true); 

      // param2
      out << setw(24) << "" << param2.affecteds       << "   "
          << setw(5)  << param2.pair_count            << "  ";

      if( j == 0 )
        print_result_line(param2.estimate.f(j), param2.estimate.standard_error(j),
                          param2.t_value(j, fi_mean), param2.p_value(j, df2, fi_mean),
                          w, pr, pw, reg[2].get_pvalues_scientific_notation(), true, fi_mean, true); 
      else
        print_result_line(param2.estimate.f(j), param2.estimate.standard_error(j),
                          param2.t_value(j, fi_mean), param2.p_value(j, df2, fi_mean),
                          w, pr, pw, reg[2].get_pvalues_scientific_notation(), true, fi_mean, false);
    }
  }
}

void
meantest_textfile::print_result_line(double est, double se, double t_value, double p_value,
                                      int w, int pr, int pw, bool scientific_notation,
                                      bool prt_exp, double exp_val, bool less)
{
  ostream& out = output_stream();

  out << fp(est, w, pr)   << "  "
      << fp(se, w, pr)    << "  ";

  if( wide_output() )
    out << fp(t_value, w) << "  ";

  if( scientific_notation )
    out << pval_scientific(p_value, 11, 4, 3);
  else
    out << pval(p_value, pw, 11, 3);

  if( prt_exp )
  {
    if( less )
      out << "  < ";
    else
      out << "  > ";

    out << fp(exp_val, 4);
  }

  out << endl;
}

//
//------------------------------------------------------------------
//

void
meantest_csvfile::print_results(const meantest_parameters& reg)
{
  ostream& out = output_stream();

  if(!out)
    return;

//  if(printed)
//   return;

  out << "Marker\tPairs\tpi_Estimate\tpi_Std.Error";
  if( wide_output() )
    out << "\tpi_T-value";
  out << "\tpi_P-value";

  if( !print_only_pi() )
  {
    if( reg.get_use_full_sibs() && !(reg.get_use_half_sibs()) )
    {
      out << "\tf0_Estimate\tf0_Std.Error";
      if( wide_output() )
        out << "\tf0_T-value";
      out << "\tf0_P-value";

      out << "\tf1_Estimate\tf1_Std.Error";
      if( wide_output() )
        out << "\tf1_T-value";
      out << "\tf1_P-value";

      out << "\tf2_Estimate\tf2_Std.Error";
      if( wide_output() )
        out << "\tf2_T-value";
      out << "\tf2_P-value";
    }
  }

  out << endl;

//  printed = true;

  for( size_t i = 0; i < reg.marker_count(); ++i )
    print_results(reg, i);
}

void
meantest_csvfile::print_results(const meantest_parameters& reg, size_t i)
{
  ostream& out = output_stream();
  const relative_pairs& pairs = sib_pairs();

  if(!out)
    return;

  const marker_parameter& param = reg[i];

  string marker = pairs.marker_name(param.marker);
  if(!marker.size())
    marker = "Marker #" + long2str(param.marker);

  long df  = param.pair_count - 1;

  double pi_mean = 0.25 + reg.get_w() * 0.5;
  if( reg.get_use_half_sibs() && !(reg.get_use_full_sibs()) )
    pi_mean = 0.25;

  out << marker
      << "\t" << param.pair_count
      << "\t" << fp(param.estimate.pi(),11,9)
      << "\t" << fp(param.estimate.standard_error(),12,10);
  if( wide_output() )
    out << "\t" << fp(param.t_value(pi_mean),10,8);

  if( reg.get_pvalues_scientific_notation() )
  {
    string s = doub2str(param.p_value(df, pi_mean), 11, 4, ios::scientific | ios::showpoint);
    out << "\t" << s;
  }
  else
    out << "\t" << fp(param.p_value(df, pi_mean),10,8);

  if( print_only_pi() )
  {
    out << endl;
    return;
  }

  if( reg.get_use_half_sibs() && !(reg.get_use_full_sibs()) )
  {
    out << endl;
    return;
  }

  for( int j = 0; j < 3; ++j )
  {
    double fi_mean = 0.5;

    if( reg.get_use_full_sibs() && !(reg.get_use_half_sibs()) )
    {
      switch( j )
      {
        case 0: fi_mean = 0.25; break;
        case 1: fi_mean = 0.5;  break;
        case 2: fi_mean = 0.25; break;
      }

      out << "\t"<< fp(param.estimate.f(j),11,9)
          << "\t"<< fp(param.estimate.standard_error(j),12,10);

      if( j == 1 )
      {
        if( wide_output() )
          out << "\t" << fp(numeric_limits<double>::quiet_NaN(),10,8);
        out << "\t" << fp(numeric_limits<double>::quiet_NaN(),10,8);
      }
      else
      {
        if( wide_output() )
          out << "\t" << fp(param.t_value(j, fi_mean),10,8);

        if( reg.get_pvalues_scientific_notation() )
        {
          string s = doub2str(param.p_value(j, df, fi_mean), 11, 4, ios::scientific | ios::showpoint);
          out << "\t" << s;
        }
        else
          out << "\t" << fp(param.p_value(j, df, fi_mean),10,8);
      }
    }
    else if( reg.get_use_half_sibs() && !(reg.get_use_full_sibs()) )
    {
      if( j > 1 )
        continue;

      out << "\t"<< fp(param.estimate.f(j),11,9)
          << "\t"<< fp(param.estimate.standard_error(j),12,10);

      if( wide_output() )
        out << "\t" << fp(param.t_value(j, fi_mean),10,8);
      //out << "\t" << fp(param.p_value(j, df, fi_mean),10,8);
      out << "\t" << fp(numeric_limits<double>::quiet_NaN(),10,8);
    }
  }

  out << endl;
}

void
meantest_csvfile::print_results(const meantest_param_vector& reg, const string& trait)
{
  ostream& out = output_stream();

  if(!out || reg.size() != 3)
  {
    out << "Cannot print regression results.  Bad size == " << reg.size() << endl;
    return;
  }

  out << "Marker\tTrait";
  out << "\tpi_0_Pairs\tpi_0_Estimate\tpi_0_Std.Error";
  if( wide_output() )
    out << "\tpi_0_T-value";
  out << "\tpi_0_P-value";

  out << "\tpi_1_Pairs\tpi_1_Estimate\tpi_1_Std.Error";
  if( wide_output() )
    out << "\tpi_1_T-value";
  out << "\tpi_1_P-value";

  out << "\tpi_2_Pairs\tpi_2_Estimate\tpi_2_Std.Error";
  if( wide_output() )
    out << "\tpi_2_T-value";
  out << "\tpi_2_P-value";

  if( !print_only_pi() )
  {
    if( reg[0].get_use_full_sibs() && !(reg[0].get_use_half_sibs()) )
    {
      out << "\tf0_0_Pairs\tf0_0_Estimate\tf0_0_Std.Error";
      if( wide_output() )
        out << "\tf0_0_T-value";
      out << "\tf0_0_P-value";  

      out << "\tf0_1_Pairs\tf0_1_Estimate\tf0_1_Std.Error";
      if( wide_output() )
        out << "\tf0_1_T-value";
      out << "\tf0_1_P-value";  

      out << "\tf0_2_Pairs\tf0_2_Estimate\tf0_2_Std.Error";
      if( wide_output() )
        out << "\tf0_2_T-value";
      out << "\tf0_2_P-value";  

      out << "\tf1_0_Pairs\tf1_0_Estimate\tf1_0_Std.Error";
      if( wide_output() )
        out << "\tf1_0_T-value";
      out << "\tf1_0_P-value";  

      out << "\tf1_1_Pairs\tf1_1_Estimate\tf1_1_Std.Error";
      if( wide_output() )
        out << "\tf1_1_T-value";
      out << "\tf1_1_P-value";  

      out << "\tf1_2_Pairs\tf1_2_Estimate\tf1_2_Std.Error";
      if( wide_output() )
        out << "\tf1_2_T-value";
      out << "\tf1_2_P-value";  

      out << "\tf2_0_Pairs\tf2_0_Estimate\tf2_0_Std.Error";
      if( wide_output() )
        out << "\tf2_0_T-value";
      out << "\tf2_0_P-value";  

      out << "\tf2_1_Pairs\tf2_1_Estimate\tf2_1_Std.Error";
      if( wide_output() )
        out << "\tf2_1_T-value";
      out << "\tf2_1_P-value";  

      out << "\tf2_2_Pairs\tf2_2_Estimate\tf2_2_Std.Error";
      if( wide_output() )
        out << "\tf2_2_T-value";
      out << "\tf2_2_P-value";  
    }
  }

  out << endl;

  for( size_t i = 0; i < reg[0].marker_count(); ++i )
    print_results(reg, i, trait);
}

void
meantest_csvfile::print_results(const meantest_param_vector& reg, size_t i, const string& trait)
{
  ostream& out = output_stream();
  const relative_pairs& pairs = sib_pairs();

  if(!out)
    return;

  const marker_parameter& param0 = reg[0][i];
  const marker_parameter& param1 = reg[1][i];
  const marker_parameter& param2 = reg[2][i];

  string marker = pairs.marker_name(param0.marker);
  if(!marker.size())
    marker = "Marker #" + long2str(param0.marker);

  long df0 = param0.pair_count - 1;
  long df1 = param1.pair_count - 1;
  long df2 = param2.pair_count - 1;

  double pi_mean = 0.25 + reg[0].get_w() * 0.5;
  if( reg[0].get_use_half_sibs() && !(reg[0].get_use_full_sibs()) )
    pi_mean = 0.25;

  out << marker
      << "\t"<< trait;

  out << "\t"<< param0.pair_count
      << "\t"<< fp(param0.estimate.pi(),13,11)
      << "\t"<< fp(param0.estimate.standard_error(),14,12);
  if( wide_output() )
    out << "\t" << fp(param0.t_value(pi_mean),12,10);
  if( reg[0].get_pvalues_scientific_notation() )
  {
    string s = doub2str(param0.p_value(df0, pi_mean), 11, 4, ios::scientific | ios::showpoint);
    out << "\t" << s;
  }
  else
    out << "\t" << fp(param0.p_value(df0, pi_mean),12,10);

  out << "\t"<< param1.pair_count
      << "\t"<< fp(param1.estimate.pi(),13,11)
      << "\t"<< fp(param1.estimate.standard_error(),14,12);
  if( wide_output() )
    out << "\t" << fp(param1.t_value(pi_mean),12,10);
  if( reg[0].get_pvalues_scientific_notation() )
  {
    string s = doub2str(param1.p_value(df1, pi_mean), 11, 4, ios::scientific | ios::showpoint);
    out << "\t" << s;
  }
  else
    out << "\t" << fp(param1.p_value(df1, pi_mean),12,10);

  out << "\t"<< param2.pair_count
      << "\t"<< fp(param2.estimate.pi(),13,11)
      << "\t"<< fp(param2.estimate.standard_error(),14,12);
  if( wide_output() )
    out << "\t" << fp(param2.t_value(pi_mean),12,10);
  if( reg[0].get_pvalues_scientific_notation() )
  {
    string s = doub2str(param2.p_value(df1, pi_mean), 11, 4, ios::scientific | ios::showpoint);
    out << "\t" << s;
  }
  else
    out << "\t" << fp(param2.p_value(df2, pi_mean),12,10);

  if( print_only_pi() )
  {
    out << endl;
    return;
  }

  if( reg[0].get_use_half_sibs() && !(reg[0].get_use_full_sibs()) )
  {
    out << endl;
    return;
  }

  for( int j = 0; j < 3; ++j )
  {
    double fi_mean = 0.5;
    
    if( reg[j].get_use_full_sibs() && !(reg[j].get_use_half_sibs()) )
    {
      switch( j )
      {
        case 0: fi_mean = 0.25; break;
        case 1: fi_mean = 0.5;  break;
        case 2: fi_mean = 0.25; break;
      }

      out << "\t"<< param0.pair_count
          << "\t"<< fp(param0.estimate.f(j),13,11)
          << "\t"<< fp(param0.estimate.standard_error(j),14,12);

      if( j == 1 )
      {
        if( wide_output() )
          out << "\t" << fp(numeric_limits<double>::quiet_NaN(),10,8);
        out << "\t" << fp(numeric_limits<double>::quiet_NaN(),10,8);
      }
      else
      {
        if( wide_output() )
          out << "\t" << fp(param0.t_value(j, fi_mean),12,10);
        if( reg[j].get_pvalues_scientific_notation() )
        {
          string s = doub2str(param0.p_value(j, df0, fi_mean), 11, 4, ios::scientific | ios::showpoint);
          out << "\t" << s;
        }
        else
          out << "\t" << fp(param0.p_value(j, df0, fi_mean),12,10);
      }

      out << "\t"<< param1.pair_count
          << "\t"<< fp(param1.estimate.f(j),13,11)
          << "\t"<< fp(param1.estimate.standard_error(j),14,12);

      if( j == 1 )
      {
        if( wide_output() )
          out << "\t" << fp(numeric_limits<double>::quiet_NaN(),10,8);
        out << "\t" << fp(numeric_limits<double>::quiet_NaN(),10,8);
      }
      else
      {
        if( wide_output() )
          out << "\t" << fp(param1.t_value(j, fi_mean),12,10);
        if( reg[j].get_pvalues_scientific_notation() )
        {
          string s = doub2str(param1.p_value(j, df1, fi_mean), 11, 4, ios::scientific | ios::showpoint);
          out << "\t" << s;
        }
        else
          out << "\t" << fp(param1.p_value(j, df1, fi_mean),12,10);
      }

      out << "\t"<< param2.pair_count
          << "\t"<< fp(param2.estimate.f(j),13,11)
          << "\t"<< fp(param2.estimate.standard_error(j),14,12);

      if( j == 1 )
      {
        if( wide_output() )
          out << "\t" << fp(numeric_limits<double>::quiet_NaN(),10,8);
        out << "\t" << fp(numeric_limits<double>::quiet_NaN(),10,8);
      }
      else
      {
        if( wide_output() )
          out << "\t" << fp(param2.t_value(j, fi_mean),12,10);
        if( reg[j].get_pvalues_scientific_notation() )
        {
          string s = doub2str(param2.p_value(j, df2, fi_mean), 11, 4, ios::scientific | ios::showpoint);
          out << "\t" << s;
        }
        else
          out << "\t" << fp(param2.p_value(j, df2, fi_mean),12,10);
      }
    }
    else if( reg[j].get_use_half_sibs() && !(reg[j].get_use_full_sibs()) )
    {
      if( j > 1 )
        continue;

      out << "\t"<< param0.pair_count
          << "\t"<< fp(param0.estimate.f(j),13,11)
          << "\t"<< fp(param0.estimate.standard_error(j),14,12);
      if( wide_output() )
        out << "\t" << fp(param0.t_value(j, fi_mean),12,10);
      if( reg[j].get_pvalues_scientific_notation() )
      {
        string s = doub2str(param0.p_value(j, df0, fi_mean), 11, 4, ios::scientific | ios::showpoint);
        out << "\t" << s;
      }
      else
        out << "\t" << fp(param0.p_value(j, df0, fi_mean),12,10);

      out << "\t"<< param1.pair_count
          << "\t"<< fp(param1.estimate.f(j),13,11)
          << "\t"<< fp(param1.estimate.standard_error(j),14,12);
      if( wide_output() )
        out << "\t" << fp(param1.t_value(j, fi_mean),12,10);
      if( reg[j].get_pvalues_scientific_notation() )
      {
        string s = doub2str(param1.p_value(j, df1, fi_mean), 11, 4, ios::scientific | ios::showpoint);
        out << "\t" << s;
      }
      else
        out << "\t" << fp(param1.p_value(j, df1, fi_mean),12,10);

      out << "\t"<< param2.pair_count
          << "\t"<< fp(param2.estimate.f(j),13,11)
          << "\t"<< fp(param2.estimate.standard_error(j),14,12);
      if( wide_output() )
        out << "\t" << fp(param2.t_value(j, fi_mean),12,10);
      if( reg[j].get_pvalues_scientific_notation() )
      {
        string s = doub2str(param2.p_value(j, df2, fi_mean), 11, 4, ios::scientific | ios::showpoint);
        out << "\t" << s;
      }
      else
        out << "\t" << fp(param2.p_value(j, df2, fi_mean),12,10);
    }
  }

  out << endl;
}

} // end of namespace SIBPAL
} // end of namespace SAGE
