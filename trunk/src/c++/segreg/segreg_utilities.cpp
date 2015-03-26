//===================================================================
//
//  File:	segreg_utilities.cpp
//
//  Author:	Stephen Gross
//
//  History:	sag Initial implementation.		Aug 03 2001
//
//  Copyright (c) 2001 R. C. Elston
//  All rights reserved
//===================================================================

#include "segreg/segreg_utilities.h"
#include "output/Output.h"

using namespace std;

namespace SAGE {
namespace SEGREG {

bool
segreg_utilities::almost_equals(double x, double y)
{
  if(x < y * 0.999999 || x > y * 1.000001) 
    return false;
  else 
    return true;
}

//===================================================================
//  test_continuous_MC(...)
//
//===================================================================
int
segreg_utilities::test_continuous_MC
   (const FPED::Multipedigree& ped_data, const model & mod, bool asc)
{
  cout << "Performing test of the member covariate calculator...";
  cout << std::endl;

  // Dump basic information

  mod.mean_sub_model.dump(cout);
  mod.var_sub_model.dump(cout);

  cout << endl;
  
  continuous_member_calculator mcc(ped_data,mod,asc);
  mcc.update();

  for(FPED::PedigreeConstIterator
          pedigree  = ped_data.pedigree_begin(); 
          pedigree != ped_data.pedigree_end(); ++pedigree)
  {
    for(FPED::MemberConstIterator
            member  = pedigree->member_begin();
            member != pedigree->member_end(); ++member)
    {
      for(size_t genotype = 0; genotype < 3; genotype++)
      {
        // Print out member info:

        cout << "Ped:" << pedigree->name();
        cout << " Member:" << member->name();
        cout << " Genotype:" << genotype;
        cout << std::endl;
 
        // Print out analysis trait info:

        cout << "Trait:\t";
        cout << OUTPUT::Double(mcc.get_composite_trait(*member)) << " ";

#ifdef DEBUG_HIGH

        cout << "= ";

        cout << mcc.get_covariate_trait_data(*member).get_analysis_trait();

        for(size_t analysis_trait_loop = 0; 
                   analysis_trait_loop < 
                   mod.comp_trait_sub_model.covariates().size();
                 ++analysis_trait_loop)
        {
          cout << " + ";
          cout << mod.comp_trait_sub_model.covariates()[analysis_trait_loop].coefficient;
          cout << " * ";  
          cout << mcc.get_covariate_trait_data(*member).
                          get_composite_cov((size_t)analysis_trait_loop);
        }

        cout << std::endl;
#endif

        // Print out expected mean info:

        cout << "Exp Trait:\t";
        cout << OUTPUT::Double(mcc.get_expected_mean(*member,(genotype_index)genotype)) << " ";

#ifdef DEBUG_HIGH

        cout << "= ";

        cout << mod.mean_sub_model.parameter( (genotype_index)genotype);

        for(size_t analysis_trait_loop = 0; 
                   analysis_trait_loop < mod.mean_cov_sub_model.covariates().size();
                 ++analysis_trait_loop)
        {
          cout << " + ";
          cout << mod.mean_cov_sub_model.covariates()[analysis_trait_loop].coefficient;
          cout << " * ";  
          cout << mcc.get_covariate_trait_data(*member).get_mean_cov(
                  (size_t)analysis_trait_loop);
        }

        cout << std::endl;
#endif

        // Print out expected variance info:

        cout << "Var:\t\t";
        cout << OUTPUT::Double(mcc.get_expected_variance(*member,(genotype_index)genotype)) << " ";

#ifdef DEBUG_HIGH

        cout << " = ";

        cout << mod.var_sub_model.parameter( (genotype_index)genotype);
        cout << " + ";

        for(size_t analysis_trait_loop = 0; 
                   analysis_trait_loop < mod.var_cov_sub_model.covariates().size();
                 ++analysis_trait_loop)
        {
          cout << mod.var_cov_sub_model.covariates()[analysis_trait_loop].coefficient;
          cout << " * ";  
          cout << mcc.get_covariate_trait_data(*member).get_variance_cov((size_t)analysis_trait_loop);
          if(analysis_trait_loop < mod.var_cov_sub_model.covariates().size()-1)
            cout << " + ";
        }

        cout << std::endl;
#endif

        // Print out expected standard deviation info:

        cout << "SD:\t\t";
        cout << OUTPUT::Double(mcc.get_expected_sd(*member,(genotype_index)genotype)) << " ";

#ifdef DEBUG_HIGH

        cout << " = sqrt(" << mcc.get_expected_variance(*member,(genotype_index)genotype);

        cout << ")";
        cout << std::endl;

#endif

        // Print out the correlation info:

        cout << "(Du(i)) = ";
        cout << mcc.get_standardization(*member,(genotype_index)genotype);

        cout << std::endl;

      } // End of genotype loop
    } // End of member loop
  } // End of subpedigree loop
  return 0;
}

void
  create_member_penetrance_table_columns
    (OUTPUT::Table&          table,
     const TypeDescription&  tdesc,
     const FPED::Member&     mem)
{
  if(mem.is_nonfounder())
    table << OUTPUT::TableColumn("M. State")
          << OUTPUT::TableColumn("F. State");

  for(TypeDescription::StateIterator state = tdesc.begin(); state != tdesc.end(); ++state)
  {
    table << OUTPUT::TableColumn(state->get_name());
  }
        
  for(FPED::MateConstIterator spouse_loop  = mem.mate_begin();
      spouse_loop != mem.mate_end();
      ++spouse_loop)
  {
    for(TypeDescription::StateIterator spstate = tdesc.begin(); spstate != tdesc.end(); ++spstate)
    {
      table << OUTPUT::Table::BEGIN_COLUMN_GROUP(std::string("Sp = ") + 
                                                 spouse_loop->name() + " " + spstate->get_name() );
      for(TypeDescription::StateIterator state = tdesc.begin(); state != tdesc.end(); ++state)
      {
        table << OUTPUT::TableColumn(state->get_name());
      }
    }
  }
}

void
  add_member_penetrance_group_to_table_row
    (OUTPUT::TableRow& row,
     const RegPenetranceCalculator& pcalc,
     const FPED::Member& mem,
     const PenetranceContext& context)
{
  RegPenetranceCalculator::penetrance_info pen_indiv;
  
  pen_indiv.member = &mem;
  
  pen_indiv.genotype = index_AA; row << pcalc.get_penetrance(pen_indiv, context);
  pen_indiv.genotype = index_AB; row << pcalc.get_penetrance(pen_indiv, context);
  pen_indiv.genotype = index_BB; row << pcalc.get_penetrance(pen_indiv, context);
}

void
  add_member_penetrance_table_row_elements
    (OUTPUT::TableRow&              row,
     const RegPenetranceCalculator& pcalc,
     const TypeDescription&         tdesc,
     const FPED::Member&            mem,
     const PenetranceContext&       context)
{
  add_member_penetrance_group_to_table_row(row, pcalc, mem, context);
  
  for(FPED::MateConstIterator spouse_loop  = mem.mate_begin();
      spouse_loop != mem.mate_end();
      ++spouse_loop)
  {
    PenetranceContext temp = context;
    
    temp.set_link_couple(mem, spouse_loop->mate());

    for(TypeDescription::StateIterator spstate = tdesc.begin(); spstate != tdesc.end(); ++spstate)
    {
      temp.set_link_spouse_state(*spstate);

      add_member_penetrance_group_to_table_row(row, pcalc, mem, temp);
    } 
  }
}

void
  add_member_penetrance_table_no_parent_row
    (OUTPUT::Table&                 table,
     const RegPenetranceCalculator& pcalc,
     const TypeDescription&         tdesc,
     const FPED::Member&            mem)
{
  OUTPUT::TableRow row;
  
  if(mem.is_nonfounder())
    row << "--" << "--";
  
  PenetranceContext context = pcalc.get_context();
  
  add_member_penetrance_table_row_elements(row, pcalc, tdesc, mem, context);

  table << row;
}

void
  add_member_penetrance_table_with_parents_rows
    (OUTPUT::Table&                 table,
     const RegPenetranceCalculator& pcalc,
     const TypeDescription&         tdesc,
     const FPED::Member&            mem)
{
  if(mem.is_founder()) return;
  
  PenetranceContext context = pcalc.get_context();

  context.set_nuclear_family(*mem.family()); 

  for(TypeDescription::StateIterator mstate = tdesc.begin(); mstate != tdesc.end(); ++mstate)
  {
    context.set_mother_state(*mstate);
    
    for(TypeDescription::StateIterator fstate = tdesc.begin(); fstate != tdesc.end(); ++fstate)
    {
      context.set_father_state(*fstate);

      OUTPUT::TableRow row;
  
      row << mstate->get_name() << fstate->get_name();
  
      add_member_penetrance_table_row_elements(row, pcalc, tdesc, mem, context);

      table << row;
    }
  }
}


OUTPUT::Table 
  generate_member_penetrance_table
    (const RegPenetranceCalculator& pcalc,
     const TypeDescription&         tdesc,
     const FPED::Member&            mem)
{
  OUTPUT::Table table(mem.pedigree()->name() + ":" + mem.name());

  create_member_penetrance_table_columns(table, tdesc, mem);
  
  add_member_penetrance_table_no_parent_row(table, pcalc, tdesc, mem);

  add_member_penetrance_table_with_parents_rows(table, pcalc, tdesc, mem);
  
  return table;
}

//===================================================================
//  test_PENETRANCE(...)
//
//===================================================================

int
segreg_utilities::test_PENETRANCE
   (const FPED::Multipedigree & ped_data, const model & mod, bool ascer)
{
  cout << "Performing test of the penetrance calculator...";
  cout << std::endl;

  RegPenetranceCalculator pcalc(ped_data,mod,ascer);

  pcalc.update();
  
  TypeDescription tdesc;
  
  tdesc.add_state("A/A");
  tdesc.add_state("A/B");
  tdesc.add_state("B/B");

  // Start pedigree loop

  for(FPED::PedigreeConstIterator
          pedigree  = ped_data.pedigree_begin(); 
          pedigree != ped_data.pedigree_end(); ++pedigree)
  {
    // Start subpedigree loop
    for(FPED::SubpedigreeConstIterator 
        spedigree  = pedigree->subpedigree_begin();
        spedigree != pedigree->subpedigree_end(); 
      ++spedigree)
    {
      // Make sure subpedigree is valid

      if(MPED::mp_utilities::has_loops(*spedigree)) 
      {
        cout << std::endl << 
          "Constituent pedigree has loops, skipping analysis...";
        continue;
      }

      // Start member loop

      for(size_t mem = 0; mem < spedigree->member_count();  ++mem)
      {
        const FPED::Member& member = spedigree->member_index(mem);
        
        cout << generate_member_penetrance_table(pcalc, tdesc, member) << endl;

      } // End of member loop
    } // End of subpedigree loop
  } // End of pedigree loop
  return 0;
}

//===================================================================
//  test_POLYGENIC_TRANSITION(...)
//
//===================================================================

int
segreg_utilities::test_POLYGENIC_TRANSITION(const FPED::Multipedigree & ped_data, const model & mod)
{
  cout << "Performing test of the polygenic transition calculator...";
  cout << std::endl;

  polygenic_transition_calculator poly_calc(mod);
  double tempvalue;

  for(size_t l=0; l<mod.fpmm_sub_model.loci(); ++l)
  {
    cout << endl << "=======================" << endl << "For l=" << l;
    for(size_t Vm=0; Vm<=((l+1)*2); ++Vm)
    {
      for(size_t Vf=0; Vf<=((l+1)*2); ++Vf)
      {
        tempvalue=0;
        cout << endl << "Summation for Vm=" << Vm << ",Vf=" << Vf << ",l=" << l;
        for(size_t Vi=0; Vi<=((l+1)*2); ++Vi)
          tempvalue += poly_calc.prob(Vi,Vm,Vf,l+1);
        cout << " = " << tempvalue;
      }
    }
  }
/*
  for(size_t l=0; l<mod.fpmm_sub_model.loci(); ++l)
  {
    cout << endl << "=======================" << endl << "For l=" << l << endl;
    bool end_of_line = true;
    for(size_t x=0; x<=((l+1)*2); ++x)
    {
      for(size_t y=0; y<=((l+1)*2); ++y)
      {
        for(size_t z=0; z<=((l+1)*2); ++z)
        {
          if(end_of_line) 
          {
            cout << endl;
            end_of_line = false;
          }
          else
          {
           cout << "   ";
            end_of_line = true;
          }

          cout << "(Vi=" << setw(2) << x << ",";
          cout << "Vm="  << setw(2) << y << ",";
          cout << "Vf="  << setw(2) << z << ") = ";
          cout << setw(13) << poly_calc.prob(x,y,z,l+1);
        }
      }
    }
  }
*/
  return 0;
}

//===================================================================
//  test_A(...)
//
//===================================================================
int
segreg_utilities::test_Regressive(
	const FPED::Multipedigree & ped_data, const model & mod)
{
  cout << std::endl;
  cout << "Performing test of the SL calculator and peeler...";

  LikelihoodElements like_elt(ped_data,mod,false);

  const TypeDescription& tdesc = like_elt.get_type_description();

  cout << std::endl;
  cout << "Performing test of the Regressive_SL calculator...";

  // Test the calculator

  for(FPED::PedigreeConstIterator
          pedigree  = ped_data.pedigree_begin();
          pedigree != ped_data.pedigree_end(); 
        ++pedigree)
  {
    for(FPED::SubpedigreeConstIterator 
        spedigree  = pedigree->subpedigree_begin();
        spedigree != pedigree->subpedigree_end();
      ++spedigree)
    {
      if(MPED::mp_utilities::has_loops(*spedigree))
      {
        cout << std::endl;
        cout << "Bad pedigree!";
        continue;
      }

      // Construct peeler, and set up bidirectional communications
      // between the peeler and the SL calculator

      regressive_peeler reg_peeler(*spedigree, like_elt);
                        
      like_elt.update();

      // Test peeler

      double antpossum = numeric_limits<double>::quiet_NaN();

      double tempvalue;

      cout << std::endl << "Testing regressive peeler...";

      cout << endl << "Summation test: ";

      FPED::MemberConstIterator ind = spedigree->member_begin();

      for( ; ind != spedigree->member_end(); ++ind)
      {
        tempvalue = 0;
        for(TypeDescription::StateIterator state = tdesc.begin();
            state != tdesc.end(); ++state)
          tempvalue += reg_peeler.anterior(*ind,*state).get_double() *
                       reg_peeler.posterior(*ind,*state).get_double();

        if(SAGE::isnan(antpossum)) 
        {
          antpossum = tempvalue;
        }
        else if(!almost_equals(tempvalue,antpossum))
//        else if(1)
//        if(1)
        {
          cout << "FAILED at indiv " << ind->name()
               << " " << tempvalue << " - " << antpossum << " = " << tempvalue - antpossum
               << std::endl 
               << "            ant          pos          ant            pos          ant          pos            Sum:"
               << endl
               << "Indiv       AA           AA           AB             AB           BB           BB";

          FPED::MemberConstIterator ind2 = spedigree->member_begin();

          for( ; ind2 != spedigree->member_end(); ++ind2)
          {
            tempvalue = 0;
            cout << endl << setw(2) << ind2->name() << "       ";
            for(TypeDescription::StateIterator state2 = tdesc.begin();
                state2 != tdesc.end(); ++state2)
            {
              cout << setw(11)
                   << reg_peeler.anterior(*ind2,*state2).get_double()
                   << " * "
                   << setw(11)
                   << reg_peeler.posterior(*ind2,*state2).get_double();

              tempvalue += reg_peeler.anterior(*ind2,*state2).get_double() *
                           reg_peeler.posterior(*ind2,*state2).get_double();

              if(state2->get_index() < 2) cout << " + "; else cout << " = ";
            }
            cout << setw(11) << tempvalue;
          }
          return 0;
        }
      }
      cout << " PASSED with value " << antpossum << endl;
    } // End of subpedigree loop
  } // End of pedigree loop
  return 0;
}

//===================================================================
//  test_FPMM(...)
//
//===================================================================
int
segreg_utilities::test_FPMM(
	const FPED::Multipedigree & ped_data, const model & mod)
{
  cout << std::endl;
  cout << "Performing test of the SL calculator and peeler..." << flush;

  polygenic_transition_calculator poly_trans_calc(mod);

  FPMM_SL  SL_calc(ped_data,mod,false);

  SL_calc.update();

  polygenic_penetrance_calculator::penetrance_info indiv;

  cout << std::endl;
  cout << "Performing test of the FPMM calculator..." << flush;

  // Test the calculator

  for(FPED::PedigreeConstIterator
          pedigree  = ped_data.pedigree_begin();
          pedigree != ped_data.pedigree_end(); 
        ++pedigree)
  {
    for(FPED::SubpedigreeConstIterator 
        spedigree  = pedigree->subpedigree_begin();
        spedigree != pedigree->subpedigree_end();
      ++spedigree)
    {
      if(MPED::mp_utilities::has_loops(*spedigree))
      {
        cout << std::endl;
        cout << "Bad pedigree!" << flush;
        continue;
      }

      // Construct peeler, and set up bidirectional communications
      // between the peeler and the SL calculator

      FPMM_peeler reg_peeler(*spedigree, mod);

      SL_calc.set_peeler(&reg_peeler);

      reg_peeler.set_SL(&SL_calc);

      // Test peeler

      double antpossum = numeric_limits<double>::quiet_NaN();

      double tempvalue;
      double tempvalue2;

      cout << std::endl << "Testing FPMM peeler..." << flush;

      cout << endl << "Summation test: " << flush;

      FPED::MemberConstIterator ind = spedigree->member_begin();

      for( ; ind != spedigree->member_end(); ++ind)
      {
        indiv.member = &*ind;

        cout << " i" << indiv.member->name() << " " << flush;
        tempvalue = 0;
        for(size_t genotype = 0; genotype < 3; ++genotype)
        {
          for(size_t polygenotype = 0; polygenotype < mod.fpmm_sub_model.max_pgt(); ++polygenotype)
          {
            tempvalue += reg_peeler.anterior(*indiv.member,
                         genetic_info((genotype_index)genotype,polygenotype)).get_double() *
                         reg_peeler.posterior(*indiv.member,
                         genetic_info((genotype_index)genotype,polygenotype)).get_double();
//            cout << "(" << tempvalue << ") " << flush;
          }
        }
        cout << "(" << tempvalue << ") " << flush;
        if(SAGE::isnan(antpossum)) 
        {
          antpossum = tempvalue;
        }
        else if(!almost_equals(tempvalue,antpossum))
        {
          cout << "FAILED at indiv " << indiv.member
               << " " << tempvalue << " - " << antpossum << " = " << tempvalue - antpossum
               << std::endl 
               << "Indiv       E(Ui=AA) {ant*pos}    E(Ui=AB) {ant*pos}    E(Ui=BB) {ant*pos}    Sum"
               << endl << flush;

          FPED::MemberConstIterator ind = spedigree->member_begin();

          for( ; ind != spedigree->member_end(); ++ind)
          {
            indiv.member = &*ind;

            tempvalue = 0;
            cout << endl << indiv.member->name() << "       " << flush;
            for(size_t genotype = 0; genotype < 3; ++genotype)
            {
              tempvalue2 = 0;
              for(size_t polygenotype = 0; polygenotype < mod.fpmm_sub_model.max_pgt(); ++polygenotype)
              {
                tempvalue += reg_peeler.anterior(*indiv.member,
                             genetic_info((genotype_index)genotype,polygenotype)).get_double() *
                             reg_peeler.posterior(*indiv.member,
                             genetic_info((genotype_index)genotype,polygenotype)).get_double();
                tempvalue2+= reg_peeler.anterior(*indiv.member,
                             genetic_info((genotype_index)genotype,polygenotype)).get_double() *
                             reg_peeler.posterior(*indiv.member,
                             genetic_info((genotype_index)genotype,polygenotype)).get_double();
              }
              cout << tempvalue2 << "   ";

              if(genotype < 2) cout << " + "; else cout << " = ";
              cout << flush;
            }
            cout << setw(11) << tempvalue << flush;
          }
          return 0;
        }
      }
      cout << "PASSED with value " << antpossum;
    } // End of subpedigree loop
  } // End of pedigree loop
  cout << " Test finished " << flush;
  return 0;
}

//===================================================================
//  test_SEGREG_CALCULATOR(...)
//
//===================================================================
int
segreg_utilities::test_SEGREG_CALCULATOR(const FPED::Multipedigree & ped_data, const model & mod)
{
  time_t start, end;  double dif;  time(&start);

  cout<<endl;
  cout<<"segreg_calculator testing start."<<endl;
  //segreg_calculator  segcal(ped_data, mod);

  //segcal.calculate().get_double();
  cout<<endl<<"segreg_calculator testing end."<<endl;

  cout<<"=========================="<<endl;

  subped_list	 pslist;

  cout<<"model_class="
      <<model_class_2_string(mod.get_model_class())<<endl; 
  cout<<"primary_trait_type="
      <<primary_type_2_string(mod.get_primary_trait_type())<<endl;   

  cout<<"primary_trait="     <<mod.get_primary_trait()<<endl;

  {
    const prevalence_sub_model& psm = mod.prev_sub_model;

    for(size_t j = 0; j < psm.get_estimate_count(); ++j)
    {        
      cout << "\nPrevalence - \n";
      for(size_t i = 0; i < psm.get_estimate_susc_covariate_count(j); ++i)
      {
        cout << "\n" << "Name(S): " << psm.get_estimate_susc_covariate_name(j, i) << "\n"
                 << "Coefficient: " << psm.get_estimate_susc_covariate_value(j, i) << "\n"
                 << std::endl;
      }
      for(size_t i = 0; i < psm.get_estimate_mean_covariate_count(j); ++i)
      {
        cout << "\n" << "Name(M): " << psm.get_estimate_mean_covariate_name(j, i) << "\n"
                 << "Coefficient: " << psm.get_estimate_mean_covariate_value(j, i) << "\n"
                 << std::endl;
      }
      for(size_t i = 0; i < psm.get_estimate_var_covariate_count(j); ++i)
      {
        cout << "\n" << "Name(V): " << psm.get_estimate_var_covariate_name(j, i) << "\n"
                 << "Coefficient: " << psm.get_estimate_var_covariate_value(j, i) << "\n"
                 << std::endl;
      }
      if(mod.get_primary_trait_type() == pt_ONSET)
        cout << (OUTPUT::Table() << (OUTPUT::TableRow() << "AGE" << psm.get_estimate_age(j)));
    }
    for(size_t j = 0; j < psm.get_constraint_count(); ++j)
    {
      cout << "P Constraint: " << j << endl;
      for(size_t i = 0; i < psm.get_constraint_susc_covariate_count(j); ++i)
      {
        cout << "\n" << "Name(S): " << psm.get_constraint_susc_covariate_name(j, i) << "\n"
                 << "Coefficient: " << psm.get_constraint_susc_covariate_value(j, i) << "\n"
                 << std::endl;
      }
      for(size_t i = 0; i < psm.get_constraint_mean_covariate_count(j); ++i)
      {
        cout << "\n" << "Name(M): " << psm.get_constraint_mean_covariate_name(j, i) << "\n"
                 << "Coefficient: " << psm.get_constraint_mean_covariate_value(j, i) << "\n"
                 << std::endl;
      }
      for(size_t i = 0; i < psm.get_constraint_var_covariate_count(j); ++i)
      {
        cout << "\n" << "Name(V): " << psm.get_constraint_var_covariate_name(j, i) << "\n"
                 << "Coefficient: " << psm.get_constraint_var_covariate_value(j, i) << "\n"
                 << std::endl;
      }
    
      if(mod.get_primary_trait_type() == pt_ONSET)
        cout << (OUTPUT::Table() << (OUTPUT::TableRow() << "AGE" << psm.get_constraint_age(j)));
    
      cout << "R:  " << psm.get_constraint_number_affected(j) << endl;
      cout << "N:  " << psm.get_constraint_sample_size(j) << endl;
    }
    
    cout << endl;
  }

  for(FPED::PedigreeConstIterator ped=ped_data.pedigree_begin(); ped!=ped_data.pedigree_end(); ++ped)
  {
    cout<<"Pedigree: "<<ped->index()<<endl;

    for(FPED::SubpedigreeConstIterator sped =ped->subpedigree_begin();
        			   sped!=ped->subpedigree_end(); ++sped)
    {
       cout<<"Subpedigree: "      <<sped->index()<<endl;
       cout<<"number of family = "<<sped->family_count()<<endl;
       cout<<"number of member = "<<sped->member_count()<<endl;

       for(FPED::FamilyConstIterator fi2=sped->family_begin(); fi2!=sped->family_end(); ++fi2)
       {
           if(fi2->offspring_count() > 9) 
              cout<<"family: "<<fi2->name()<<setw(5)<<" "
		  <<fi2->offspring_count()<<" child"<<endl;
       }
       
       if(MPED::mp_utilities::has_loops(*sped))
       {  
	  cout << "Bad pedigree!" <<endl;  continue; 
       }

       pslist.push_back(&*sped);
    }
  }

  cout<<endl<<endl<<"<<< test_MLM result begin >>>"<<endl;

  for(subped_list::iterator ps = pslist.begin(); ps!= pslist.end(); ps++)
  {
       FPED::SubpedigreeConstPointer sped = *ps;

       time_t start1, end1; double dif1; time(&start1);

       binary_member_calculator bmc(ped_data,mod,false);

       bmc.update();
       
       MlmLikelihoodElements lelt(ped_data,mod,false);

       mlm_peeler mpl(*sped,&bmc, ped_data,mod,lelt,false);


       cout<<endl<<"Summary of testing:"<<endl
                 <<"-------------------"<<endl;

       cout<<"ind loc Yi   "
           <<setfill(' ')<<setw(13)<<"antAA"
           <<setfill(' ')<<setw(13)<<"posAA"
           <<setfill(' ')<<setw(13)<<"antAB"
	   <<setfill(' ')<<setw(13)<<"posAB"
	   <<setfill(' ')<<setw(13)<<"antBB"
	   <<setfill(' ')<<setw(13)<<"posBB"
	   <<setfill(' ')<<setw(10)<<"sum"
	   <<setfill(' ')<<setw(5)<<"time"<<endl<<endl;

       for(size_t indiv = 0; indiv < sped->member_count(); ++indiv)
       {
           const FPED::Member& ind = sped->member_index(indiv);

           time_t start2, end2; double dif2; time(&start2);

           double tempvalue   = 0.0;
          
           double trait_value = 0.0;
           if( bmc.is_member_valid(ind))
               trait_value = bmc.get_aff_status(ind);
           else
               trait_value = -9;

           cout<<setw(3)<<ind.name()
               <<setw(3)<<indiv
               <<setw(3)<<trait_value;

           //if(indiv==0) {

           for(TypeDescription::StateIterator state = lelt.get_type_description().begin();
               state != lelt.get_type_description().end(); ++state)
           //for(size_t genotype = 0; genotype < 3; ++genotype)
           {
               cout << setw(2) <<" "
		    << setw(9)<<setprecision(4)
                    << mpl.anterior (ind,*state).get_double()
                    << " * " 
                    << setw(9)
                    << mpl.posterior(ind,*state).get_double();
           
               tempvalue += mpl.anterior (ind,*state).get_double() *
                            mpl.posterior(ind,*state).get_double();
         
               if(state->get_index() != index_BB)  cout << " + "; 
               else                                cout << " = ";

           }// for genotype
           //}// if indiv=0

           cout<<setw(9)<< tempvalue;

           time(&end2); dif2 = difftime(end2,start2);
           if(dif2 < 60.0)
                 cout<<setw(5)<<dif2   <<" sec"<<endl;
           else    
                 cout<<setw(5)<<setprecision(3)<<dif2/60<<" min"<<setprecision(4)<<endl;

       }// for indiv

       time(&end1);
       dif1 = difftime(end1,start1);
       if(dif1<60.0) cout<<endl<<"Pedigree "<<sped->pedigree()->index()<<" takes "<<dif1   <<" sec"<<endl;
       else	     cout<<endl<<"Pedigree "<<sped->pedigree()->index()<<" takes "<<dif1/60<<" min"<<endl;  

  }// pslist

  time(&end);
  dif = difftime(end,start);     
  if(dif<60.0)	cout<<endl<<"The data set takes "<<dif   <<" sec"<<endl;
  else          cout<<endl<<"The data set takes "<<dif/60<<" min"<<endl;
  
  return 0;
}

//===================================================================
//  test_MLM(...)
//===================================================================
int
segreg_utilities::test_MLM(const FPED::Multipedigree & ped_data, const model & mod, bool ascer)
{
  cout << "Performing test of the MLM Anterior/Posterior...";
  cout << std::endl << endl;

  // Dump basic MLM information

  mod.susc_sub_model.dump(cout);

  cout << endl;

  binary_member_calculator bmc(ped_data,mod,ascer);

  bmc.update();

  for(FPED::PedigreeConstIterator ped=ped_data.pedigree_begin(); ped!=ped_data.pedigree_end(); ++ped)
  {
    for(FPED::SubpedigreeConstIterator sped =ped->subpedigree_begin();
                                   sped!=ped->subpedigree_end(); ++sped)
    {
       if(MPED::mp_utilities::has_loops(*sped))
         continue;

       MlmLikelihoodElements lelt(ped_data,mod,false);

       mlm_peeler mpl(*sped,&bmc, ped_data,mod,lelt,false);

       regressive_peeler::penetrance_info indiv;

       cout<<"<<< test_MLM result begin >>>"<<endl;

       //-----------------------------------------------------------------

       for(FPED::MemberConstIterator mi=sped->member_begin(); mi!= sped->member_end(); ++mi)
       //for(indiv.member = 0; indiv.member < ps.individual_count(); ++indiv.member)
       {
           indiv.member = &*mi;

      	   log_double mlm_prob(0.0);
           log_double ant_prob(0.0);
           log_double pos_prob(0.0);

           cout<<"name="<<mi->name()<<"  ";
           cout<<"location(i)="<< indiv.member->index()<<endl;
           cout<<"child_count = "<<mi->offspring_count()<<endl;

           for(TypeDescription::StateIterator state = lelt.get_type_description().begin();
               state != lelt.get_type_description().end(); ++state)
      	   {
               ant_prob = mpl.anterior (*indiv.member, *state);
               pos_prob = mpl.posterior(*indiv.member, *state);
               mlm_prob += ant_prob * pos_prob;
               
               cout << (OUTPUT::Table()
                    << (OUTPUT::TableRow() << "genotype" << state->get_index())
                    << (OUTPUT::TableRow() << "ant_prob" << ant_prob.get_double())
                    << (OUTPUT::TableRow() << "pos_prob" << pos_prob.get_double())
                    << (OUTPUT::TableRow() << "-----------------------------------L(P," << state->get_index() << ")=" << mlm_prob.get_double()));
      	   }
           cout<<" indiv="<<mi->name();

           cout << (OUTPUT::Table() << (OUTPUT::TableRow() << "--------------------------------------------L(P)  =" << mlm_prob.get_double()));

       } //end of ind loop
       cout<<"<<< test_MLM result end >>>"<<endl;       
    } //end of subped
  } // end of ped

  return 0;
}

//===================================================================
//  test_binary_MC(...)
//
//===================================================================
int
segreg_utilities::test_binary_MC
    (const FPED::Multipedigree & ped_data, const model & mod, bool ascer)
{
  cout << "Performing test of the binary member calculator...";
  cout << std::endl << endl;

  // Dump basic binary information

  mod.susc_sub_model.dump(cout);

  cout << endl;

  binary_member_calculator mcc(ped_data,mod,ascer);
  mcc.update();

  for(FPED::PedigreeConstIterator pedigree  = ped_data.pedigree_begin(); 
                              pedigree != ped_data.pedigree_end(); ++pedigree)
  {
    for(FPED::MemberConstIterator member  = pedigree->member_begin();
                              member != pedigree->member_end(); ++member)
    {
      for(size_t genotype = 0; genotype < 3; genotype++)
      {
        // Print out member info:

        cout << "Ped: "       << pedigree->name() << "\tMember: "   << member->name()
             << "\tGenotype: " << genotype << endl
             << "--------------------------------------------------------" << endl;
 
        // Print out analysis trait info:

        cout << (OUTPUT::Table() 
             << (OUTPUT::TableRow()
             << "Valid"
             << mcc.is_member_valid   (*member)
             << "Affection"  
             << mcc.get_aff_status    (*member)
             << "Exp Susc"   
             << mcc.get_expected_susc (*member, (genotype_index) genotype)
             << "Penetrance" 
             << mcc.get_penetrance    (*member, (genotype_index) genotype)));

        if(mod.get_model_class() == model_FPMM)
        {
          // Print out expected polygenic mean info:

          for(size_t polygenotype = 0; polygenotype < mod.fpmm_sub_model.max_pgt(); ++polygenotype)
          {
            cout << "(Poly = " << polygenotype << ") ";
            cout << mcc.get_penetrance(*member,(genotype_index)genotype,polygenotype);
            cout << "\t";
          }

          cout << endl;
        }
      } // End of genotype loop
    } // End of member loop
  } // End of subpedigree loop
  return 0;
}

ostream& operator<<(ostream& o, const TypeDescription::State& s)
{
  o << s.get_name() << ':' << s.get_index();

  return o;
}

int
segreg_utilities::test_binary_fam_resid
    (const FPED::Multipedigree & ped_data, const model & mod, bool ascer)
{
  cout << "Performing test of the binary residual adjustment calculator...";
  cout << std::endl << endl;

  // Dump basic binary information

  mod.susc_sub_model.dump(cout);

  cout << endl;

  binary_member_calculator mcc(ped_data,mod,ascer);
  mcc.update();
  
  PED_CALC::ApproximateFamResidAdj<genotype_index, FPED::Multipedigree,
                                   boost::counting_iterator<int> >
      approx_adj
        (ResidualGetter(mod.resid_sub_model),
         boost::bind(&binary_member_calculator::get_aff_status,
                     boost::ref(mcc), _1),
         boost::bind(&binary_member_calculator::get_expected_susc,
                     boost::ref(mcc), _1, _2),
         boost::bind(&transmission_sub_model::prob, boost::ref(mod.transm_sub_model),
                     _4, _2, _3),
         boost::counting_iterator<int>(index_AA),
         boost::counting_iterator<int>(index_INVALID));

  for(FPED::PedigreeConstIterator p = ped_data.pedigree_begin(); p != ped_data.pedigree_end(); ++p)
  {
    std::cout << PED_CALC::generate_fra_test_output
                    (approx_adj, *p,
                     boost::counting_iterator<int>(index_AA),
                     boost::counting_iterator<int>(index_INVALID));
  }

  return 0;
}



//===================================================================
//  test_onset_MC(...)
//
//===================================================================
int
segreg_utilities::test_onset_MC
    (const FPED::Multipedigree & ped_data, const model & mod, bool ascer)
{
  cout << "Performing test of the onset member calculator...";
  cout << std::endl << endl;

  // Dump basic onset information

  mod.ons_sub_model.dump(cout);
  mod.mean_sub_model.dump(cout);
  mod.var_sub_model.dump(cout);
  mod.susc_sub_model.dump(cout);
  mod.fpmm_sub_model.dump(cout);

  cout << endl;

  for(size_t i = 0; i < mod.fpmm_sub_model.max_pgt(); ++i)
    cout << "Poly " << i << "\t: " << mod.fpmm_sub_model.mean(i) << endl;

  cout << endl << endl;  

  onset_member_calculator mcc(ped_data,mod,ascer);

  mcc.update();

  for(FPED::PedigreeConstIterator pedigree  = ped_data.pedigree_begin(); 
                              pedigree != ped_data.pedigree_end(); ++pedigree)
  {
    for(FPED::MemberConstIterator member  = pedigree->member_begin();
                              member != pedigree->member_end(); ++member)
    {
      for(size_t genotype = 0; genotype < 3; genotype++)
      {
        size_t mem       = member->index();
        size_t abs_index = mcc.pedigree_index_map.find(&*pedigree)->second + mem;

        // Print out member info:

        cout << "Ped: "       << pedigree->name() << "\tMember: "   << member->name()
             << "\tGenotype: " << genotype << endl
             << "--------------------------------------------------------" << endl;
 
        // Print out analysis trait info:

        cout << (OUTPUT::Table() 
             << (OUTPUT::TableRow() 
             << "Affection"
             << mcc.get_aff_status(*member)
             << "Onset"   
             << mcc.age_onsets[abs_index] 
             << " -> "
             << mcc.get_age_onset(*member)
             << "Exam"   
             << mcc.age_exams[abs_index] 
             << " -> "
             << mcc.get_age_exam(*member)));

        // Print out expected polygenic mean info:

        for(size_t polygenotype = 0; polygenotype < mod.fpmm_sub_model.max_pgt(); ++polygenotype)
        {
          cout << "(Mean v=" << polygenotype << ") ";
          cout << OUTPUT::Double(mcc.get_expected_age_onset(*member,(genotype_index)genotype,polygenotype));

          cout << "\t";
        }

        cout << endl;

        // Print out expected alpha info:

        cout << "(alpha_i) ";
        cout << OUTPUT::Double(mcc.get_alpha_i(*member,(genotype_index)genotype));

        cout << " = pi / sqrt(3 * (";

        cout << mod.var_sub_model.parameter( (genotype_index)genotype);

        for(size_t analysis_trait_loop = 0; 
                   analysis_trait_loop < mod.var_cov_sub_model.covariates().size();
                 ++analysis_trait_loop)
        {
          cout << " + ";

          cout << mod.var_cov_sub_model.covariates()[analysis_trait_loop].coefficient;
          cout << " * ";  
          cout << mcc.var_cov_data[analysis_trait_loop][abs_index];
        }

        cout << "))" << std::endl;

        // Print out expected polygenic susceptibility info:

        for(size_t polygenotype = 0; polygenotype < mod.fpmm_sub_model.max_pgt(); ++polygenotype)
        {
          cout << "(Susc v=" << polygenotype << ") ";
          cout << OUTPUT::Double(mcc.get_expected_susc(*member,(genotype_index)genotype,polygenotype))
               << '\t';
          
        }
        cout << endl;

      } // End of genotype loop
    } // End of member loop
  } // End of subpedigree loop
  return 0;
}

int
segreg_utilities::test_polygenic_penetrance
   (const FPED::Multipedigree & ped_data, const model & mod, bool ascer)
{
  cout << "Performing test of the polygenic penetrance calculator...";
  cout << std::endl;

  polygenic_penetrance_calculator pcalc(ped_data,mod,ascer);

  pcalc.update();

  polygenic_transition_calculator poly_trans_calc(mod);

  FPMM_SL translation_calc(ped_data,mod,true);

  // Start pedigree loop

  for(FPED::PedigreeConstIterator pedigree  = ped_data.pedigree_begin(); 
                              pedigree != ped_data.pedigree_end(); ++pedigree)
  {
    // Start subpedigree loop
    for(FPED::SubpedigreeConstIterator 
        spedigree  = pedigree->subpedigree_begin();
        spedigree != pedigree->subpedigree_end(); 
      ++spedigree)
    {
      // Make sure subpedigree is valid

      if(MPED::mp_utilities::has_loops(*spedigree))
      {
        cout << std::endl << 
          "Invalid subpedigree (looped), skipping analysis...";
        continue;
      }

      // Start member loop

      polygenic_penetrance_calculator::penetrance_info pen_indiv;

      for(FPED::MemberConstIterator member  = pedigree->member_begin();
                                member != pedigree->member_end(); ++member)
      {
        pen_indiv.member = &*member;

        cout << endl << "Member " << member->name() << endl;

        OUTPUT::Table t;

        // (polygenic with just the individual)
        for(size_t genotype = 0; genotype < 3; genotype++)
        {
          OUTPUT::TableRow r;
          
          r << "G=" << genotype << " : ";

          pen_indiv.genotype = (genotype_index)genotype;

          for(size_t polygenotype = 0; polygenotype < mod.fpmm_sub_model.max_pgt(); ++polygenotype)
          {
            pen_indiv.polygenotype = polygenotype;

            r << "V=" << polygenotype << " : " << pcalc.get_polygenic_penetrance(pen_indiv);

          } // End of polygenotype loop
          
          t << r;

        } // End of genotype loop
        
        cout << t;
        
      } // End of member loop
    } // End of subpedigree loop
  } // End of pedigree loop
  return 0;
}

void
segreg_utilities::test_prev_calculation(const model& mod)
{
  if(mod.prev_sub_model.get_constraint_count() == 0) return;

  // Dump out the current set of prevalence constraints

  cout << "Prevalence Constraints: " << endl << endl;

  const prevalence_sub_model& psm = mod.prev_sub_model;

  for(size_t j = 0; j < psm.get_constraint_count(); ++j)
  {
    cout << "P Constraint: " << j << endl;
    for(size_t i = 0; i < psm.get_constraint_susc_covariate_count(j); ++i)
    {
      cout << "\n" << "Name(S): " << psm.get_constraint_susc_covariate_name(j, i) << "\n"
               << "Coefficient: " << psm.get_constraint_susc_covariate_value(j, i) << "\n"
               << std::endl;
    }
    for(size_t i = 0; i < psm.get_constraint_mean_covariate_count(j); ++i)
    {
      cout << "\n" << "Name(M): " << psm.get_constraint_mean_covariate_name(j, i) << "\n"
               << "Coefficient: " << psm.get_constraint_mean_covariate_value(j, i) << "\n"
               << std::endl;
    }
    for(size_t i = 0; i < psm.get_constraint_var_covariate_count(j); ++i)
    {
      cout << "\n" << "Name(V): " << psm.get_constraint_var_covariate_name(j, i) << "\n"
               << "Coefficient: " << psm.get_constraint_var_covariate_value(j, i) << "\n"
               << std::endl;
    }
  
    if(mod.get_primary_trait_type() == pt_ONSET)
      cout << (OUTPUT::Table() << (OUTPUT::TableRow() << "AGE" << psm.get_constraint_age(j)));

    cout << "R:  " << psm.get_constraint_number_affected(j) << endl;
    cout << "N:  " << psm.get_constraint_sample_size(j) << endl;
  }
  
  cout << endl;

  cout << "Calculated Penalty: " << psm.get_prevalence_penalty() << endl;

}

void
segreg_utilities::test_mlm_pairs(const PedigreeDataSet& ped_data, const model & mod)
{
  cout << "Testing MLM Correlation Verifier" << endl;
  cout << "================================" << endl;

  cout << endl << mod.get_title() << endl << endl;

  binary_member_calculator bmc(*ped_data.get_raw_data(),mod,false);

  bmc.update();

  MlmCorrelationVerifier mpc(ped_data.get_subpedigrees(), mod.resid_sub_model, bmc);

  mpc.dump_valid_pairs();

  cout << mpc.current_correlation_estimates_are_valid() << endl;
}
void
segreg_utilities::test_mlm_correlations(const PedigreeDataSet& ped_data, const model & mod)
{
  cout << "Testing MLM correlations" << endl;
  cout << "========================" << endl;

  cout << endl << mod.get_title() << endl << endl;

  binary_member_calculator bmc(*ped_data.get_raw_data(),mod,false);

  bmc.update();

  MlmResidCorrelationCalculator mrc(bmc, ped_data.get_subpedigrees(), mod);

  mrc.dump_data(mod);
}

}
}
