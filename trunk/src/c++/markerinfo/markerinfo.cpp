//==================================================================
//  File:     markerinfo.cpp                         
//
//  Author:   Yeunjoo Song                      
//
//  History:  Initial implementation.  Mar. 2002
//
//  Copyright (c) 2002 R. C. Elston
//==================================================================

#include <string>
#include <cassert>
#include <functional>
#include <vector>
#include <iomanip>
#include "LSF/LSFinit.h"
#include "error/bufferederrorstream.h"
#include "gelim/geno_eliminate.h"
#include "gelim/pedigree_region.h"
#include "markerinfo/markerinfo.h"
#include "markerinfo/input.h"

using namespace SAGE;

MARKERINFO::MARKERINFO(int argc, char **argv)
          : APP::SAGEapp(APP::APP_MARKERINFO, true, argc, argv),
            my_consistent_out(false), my_pedigree_out(false),
            my_summary_file(NULL), my_pedigree_file(NULL), my_table_file(NULL)
{
  LSFInit();
}

MARKERINFO::~MARKERINFO()
{}

int
MARKERINFO::main()
{
  // Create input data.
  //
  markerinfo_data mdata(name, debug());

  print_title(mdata.info());

  // Get the error stream and make it available locally
  //
  cerrorstream& errors = mdata.errors();

  mdata.input(argc, argv);

  // Start markerinfo analysis.
  //
  cout << endl << "MARKERINFO analysis......." << endl << endl << flush;

  if( mdata.analysis().size() )
    parse_analyses(mdata.analysis()[0], mdata.pedigrees(), errors);

  boost::scoped_ptr<ofstream> s_file;
  boost::scoped_ptr<ofstream> p1_file;
  boost::scoped_ptr<ofstream> p2_file;
  boost::scoped_ptr<ofstream> t_file;

  string   sum_filename;
  string   ped_filename;
  string   par_filename;
  string   tab_filename;

  ofstream sum_file;
  ofstream ped_file;
  ofstream par_file;
  ofstream tab_file;

  if( my_outfile_name.size() )
  {
    sum_filename = my_outfile_name + ".out";
    sum_file.open( sum_filename.c_str() );
    if(sum_file)
      print_title(sum_file);
    my_summary_file = &sum_file;

    if( my_pedigree_out )
    {
      ped_filename = my_outfile_name + "_clean.ped";
      ped_file.open( ped_filename.c_str() );
      my_pedigree_file = &ped_file;

      par_filename = my_outfile_name + "_clean.par";
      par_file.open( par_filename.c_str() );
      my_parameter_file = &par_file;

      //tab_filename = my_outfile_name + ".tab";
      //tab_file.open( tab_filename.c_str() );
      //if(tab_file)
      //  print_title(tab_file);
      //my_table_file = &tab_file;
    }
  }
  else
  {
    sum_filename = "markerinfo.out";
    if( s_file.get() == NULL )
    {
      s_file.reset( new ofstream(sum_filename.c_str()) );
      if(*s_file)
        print_title(*s_file);
    }      
    my_summary_file = s_file.get();

    if( my_pedigree_out )
    {
      ped_filename = "markerinfo_clean.ped";
      if( p1_file.get() == NULL )
      {
        p1_file.reset( new ofstream(ped_filename.c_str()) );
      }      
      my_pedigree_file = p1_file.get();

      par_filename = "markerinfo_clean.par";
      if( p2_file.get() == NULL )
      {
        p2_file.reset( new ofstream(par_filename.c_str()) );
      }      
      my_parameter_file = p2_file.get();

      //tab_filename = "markerinfo.tab";
      //if( t_file.get() == NULL )
      //{
      //  t_file.reset( new ofstream(tab_filename.c_str()) );
      //  if(*t_file)
      //    print_title(*t_file);
      //}      
      //my_table_file = t_file.get();
    }
  }

  run_analyses(mdata);

  cout << endl << "Analysis complete!" << endl << endl;

  print_inf_banner(cout);

  return EXIT_SUCCESS;
}

void
MARKERINFO::parse_analyses(LSFBase* params, const RPED::RefMultiPedigree& mp, cerrorstream& errors)
{
  AttrVal v = attr_value(params, "out");
  if( !v.has_value() || !v.String().size() )
    v = attr_value(params, "output");
  if( v.has_value() && v.String().size() )
    my_outfile_name = v.String();

  if( !params->List() )
    return;

  LSFList::const_iterator i;
  AttrVal a;
  for( i = params->List()->begin(); i != params->List()->end(); ++i)
  {
    if( !*i || !((*i)->name().size()) ) continue;

    string name = toUpper( (*i)->name() );

    if( name == "SAMPLE_ID" )
    {
      a = attr_value(*i, "SAMPLE_ID", 0);
      if( a.has_value() )
      {
        size_t string_index = mp.info().string_find(a.String());

        if( string_index == (size_t) - 1 )
          errors << priority(information) << "Invalid sample_id name : "
                 << a.String() << "\n             Skipped... " << endl;
        else
          my_sample_ids.push_back(string_index);
      }
    }
    else if( name == "CONSISTENT_OUT" )
    {
      a = attr_value(*i, "CONSISTENT_OUT", 0);
      if( a.has_value() )
      {
        if(    toUpper(a.String()) == "TRUE"
            || toUpper(a.String()) == "YES" )
          my_consistent_out = true;
      }
    }
    else if( name == "PEDIGREE_SKIP" )
    {
      a = attr_value(*i, "PEDIGREE_SKIP", 0);
      if( a.has_value() )
      {
        if( mp.pedigree_find(a.String() ) )
        {
          size_t pedigree_index = mp.pedigree_find(a.String())->index();
          my_pedigree_skips.insert(pedigree_index);
        }
        else
          errors << priority(information) << "Invalid pedigree name : "
                 << a.String() << "\n             Skipped... " << endl;
      }
    }
    else if( name == "PEDIGREE_OUT" )
    {
      a = attr_value(*i, "PEDIGREE_OUT", 0);
      if( a.has_value() )
      {
        if(    toUpper(a.String()) == "TRUE"
            || toUpper(a.String()) == "YES" )
          my_pedigree_out = true;
      }
    }
    else
      errors << priority(information) << "Unknown option : "
             << name << "\n             Skipped... " << endl;
  }

  return;
}

void
MARKERINFO::run_analyses(const markerinfo_data& mdata)
{
  genotype_eliminator gelim;

  const RPED::RefMultiPedigree& ref_mp = mdata.pedigrees();

  FPED::Multipedigree fped(ref_mp);
  FPED::MPFilterer::add_multipedigree(fped, ref_mp);

  fped.construct();

  FPED::PedigreeConstIterator ped = fped.pedigree_begin();

  vector< string > ped_name_sorted;

  for( ; ped != fped.pedigree_end(); ++ped )
  {
    const FPED::Pedigree& pedigree = *ped;

    ped_name_sorted.push_back(pedigree.name());
  }

  sort(ped_name_sorted.begin(), ped_name_sorted.end());

  for( size_t ped_name = 0; ped_name < ped_name_sorted.size(); ++ped_name )
  {
    const FPED::Pedigree& pedigree = *fped.pedigree_find(ped_name_sorted[ped_name]);

    char old_fill = cout.fill('.');

    string s = "'" + pedigree.name() + "'";

    if( my_pedigree_skips.find(pedigree.index()) != my_pedigree_skips.end() )
    {
      cout << "Skipping pedigree " << s << "........." << endl << flush;
      continue;
    }

    cout << "Processing pedigree " << s << flush;

#if 0
    cout << "member in mped : " << endl;
    RPED::filtered_multipedigree::pedigree_type::member_const_iterator m_i = pedigree.member_begin();
    for( ; m_i != pedigree.member_end(); ++m_i )
      cout << "name = " << m_i->name() << ", index = " << m_i->index() << endl;
    cout << endl;
#endif

    FPED::SubpedigreeConstIterator subped     = pedigree.subpedigree_begin();
    FPED::SubpedigreeConstIterator subped_end = pedigree.subpedigree_end();
    for( ; subped != subped_end; ++subped)
    {
      //string c;
      //cin >> c;

      pedigree_region ped_r(*subped, mdata.get_region(), mdata.errors(), false, false);

      gelim.set_subpedigree(*subped);

      for( size_t i = 0; i < ped_r.inheritance_model_count(); ++i )
      {
        MLOCUS::inheritance_model& i_model = ped_r[i];

        if( !i_model.penetrance_informative() )
        {
          gelim.mark_uninformative(i_model, i);
          continue;
        }

        // Generate errors
        //  - only family errors first.
        FPED::FamilyConstIterator f       = subped->family_begin();
        FPED::FamilyConstIterator fam_end = subped->family_end();
        for( ; f != fam_end; ++f )
        {
            gelim.process_family(i_model, i, *f, inconsistency_handler::nuclear_family,
                                 genotype_eliminator::none, false);
        }

        //  - mendelian errors.
        gelim.process(i_model, i);
      }
    }

    cout << left << setw(22-s.size()) << "" << "done."<< endl;
    cout.fill(old_fill);
  }

  // Produce output
  ostream& analysis_out = *my_summary_file;

  analysis_out << endl
               << "=====================================================================" << endl
               << "  Markerinfo Analysis Output" << endl
               << "=====================================================================" << endl << endl;

  if(gelim.get_errors().family_begin() != gelim.get_errors().family_end())
  {
    inconsistency_table_printer printer(analysis_out);

    //Not yet!
    //printer.set_columns();

    analysis_out << "---------------------------------------------------------------------" << endl
                 << "    Part 1.1: Number of Inconsistencies per pedigree" << endl
                 << "---------------------------------------------------------------------" << endl << endl;

    printer.print_pedigree_table(gelim.get_errors());

    analysis_out << endl;

    analysis_out << "---------------------------------------------------------------------" << endl
                 << "    Part 1.2: Number of Inconsistencies per marker" << endl
                 << "---------------------------------------------------------------------" << endl << endl;

    printer.print_marker_table(gelim.get_errors(), fped.info().markers());

    analysis_out << endl << endl;

    analysis_out << "---------------------------------------------------------------------" << endl
                 << "    Part 2: Inconsistencies" << endl
                 << "---------------------------------------------------------------------" << endl << endl;

    // find the name for missing allele.
    const MLOCUS::inheritance_model_map& imap = fped.info().markers();
    string mv  = imap[0].missing_allele_name();
    if( mv == "*missing" )
      mv = "?";

    analysis_out << "    missing allele code = " << mv << endl << endl;

    printer.print_inconsistency_table(gelim.get_errors(), fped.info().markers(), my_sample_ids, my_consistent_out);

    if( my_pedigree_out )
    {
      print_resetted_pedigree(gelim.get_errors(), fped);
      //print_marker_pedigree_table(gelim.get_errors(), fped);
    }
  }
  else
  {
    inconsistency_table_printer printer(analysis_out);


    analysis_out << "---------------------------------------------------------------------" << endl
                 << "    Part 1.1: Number of Inconsistencies per pedigree" << endl
                 << "---------------------------------------------------------------------" << endl << endl;

    printer.print_pedigree_table(gelim.get_errors(), false);

    analysis_out << endl;

    analysis_out << "---------------------------------------------------------------------" << endl
                 << "    Part 1.2: Number of Inconsistencies per marker" << endl
                 << "---------------------------------------------------------------------" << endl << endl;

    printer.print_marker_table(gelim.get_errors(), fped.info().markers(), false);

    analysis_out << endl << endl;

    analysis_out << "  No inconsistencies detected!" << endl;
  }

  return;
}

void
MARKERINFO::print_resetted_pedigree(const inconsistency_handler& h, const FPED::Multipedigree& fped)
{
  const MLOCUS::inheritance_model_map& imap = fped.info().markers();

  // Construct subpedigree inconsistency map.
  //
  map< const FPED::SubpedigreeConstPointer, vector<bool> > sped_incon_map;

  inconsistency_handler::incon_family_iterator fi = h.family_begin();

  for( ; fi != h.family_end(); ++fi )
  {
    const FPED::FamilyConstPointer      fam_id  = fi->first;
    const FPED::SubpedigreeConstPointer sped_id = fam_id->subpedigree();

    if( !sped_incon_map[sped_id].size() )
      sped_incon_map[sped_id].resize(imap.size(), false);

    const inconsistency_handler::family_error_type& fam = fi->second;

    inconsistency_handler::family_error_type::const_iterator fam_member = fam.begin();

    for( ; fam_member != fam.end(); ++fam_member )
    {
      if( fam_member->second.size() )
      {
        inconsistency_handler::error_map::const_iterator e = fam_member->second.begin();

        for( ; e != fam_member->second.end(); ++e )
        {
          size_t m =  e->first;

          sped_incon_map[sped_id][m] = true;
        }
      }
    }
  }

#if 0
  cout << "Dump sped_incon_map:" << endl;
  map< FPED::SubpedigreeConstPointer, vector<bool> >::const_iterator sp = sped_incon_map.begin();
  for( ; sp != sped_incon_map.end(); ++sp )
  {
    cout << "sped " << sp->first->pedigree()->name() << "-s" << sp->first->name();
    for( size_t m = 0; m < sp->second.size(); ++m )
    {
      cout << ", " << sp->second[m];
    }
    cout << endl;
  }
#endif

  // Print out cleaned pedigree file & parameter file.
  //
  ostream& out = *my_pedigree_file;
  ostream& par = *my_parameter_file;

  // id fields
  //
  out << "pedid\tindid\tdad\tmom\tsex";

  par << "pedigree\n"
      << "{\n"
      << "  delimiters = \"\\t\"\n"
      << "  delimiter_mode = single\n"
      << "  individual_missing_value = \"" << fped.info().individual_missing_code() << "\"\n\n"
      << "  sex_code, male = \"" << fped.info().sex_code_male() << "\""
      << ", female = \"" << fped.info().sex_code_female() << "\""
      << ", missing = \"" << fped.info().sex_code_unknown() << "\"\n\n"
      << "  pedigree_id = pedid\n"
      << "  individual_id = indid\n"
      << "  parent_id = dad\n"
      << "  parent_id = mom\n"
      << "  sex_field = sex\n\n";

  // trait field
  //
  for( size_t t = 0; t < fped.info().trait_count(); ++t )
  {
    const RPED::RefTraitInfo& t_info = fped.info().trait_info(t);

    if(    t_info.name() == "SEX_CODE"
        || t_info.name() == "FAMILIAL_INDICATOR"
        || t_info.name() == "FOUNDER_INDICATOR"
        || t_info.name() == "PEDIGREE_SIZE" )
      continue;

    out << "\t" << t_info.name();

    par << "  trait = \"" << t_info.name() << "\"";
    if( t_info.type() == RPED::RefTraitInfo::binary_trait )
      par << ", binary, affected = \"" << t_info.string_affected_code() << "\""
          << ", unaffected = \"" << t_info.string_unaffected_code() << "\"";
    par << ", missing = \"" << t_info.string_missing_code() << "\"\n";
  }
  par << "\n";

  // string field
  //
  for( size_t s = 0; s < fped.info().string_count(); ++s )
  {
    out << "\t" << fped.info().string_info(s).name();

    par << "  string = " << fped.info().string_info(s).name() << "\n";
  }
  par << "\n";

  // marker field
  //

  //  check first to see autosomal, sex-linked, or mixed
  set<MLOCUS::GenotypeModelType> m_type;
  set<string> m_miss;

  for( size_t m = 0; m < imap.size(); ++m )
  {
    const MLOCUS::inheritance_model& model = imap[m];

    out << "\t" << model.name();

    m_type.insert(model.get_model_type());
    m_miss.insert(model.missing_allele_name());
  }
  out << "\n";

  if( m_type.size() > 1 || m_miss.size() > 1 )
  {
    for( size_t m = 0; m < imap.size(); ++m )
    {
      const MLOCUS::inheritance_model& model = imap[m];

      par << "  marker = \"" << model.name() << "\"";

      if( m_miss.size() > 1 )
        par << ", missing = \"" << model.missing_allele_name() << "\"";

      if( model.get_model_type() == MLOCUS::X_LINKED )
        par << ", x_linked";
      else if( model.get_model_type() == MLOCUS::Y_LINKED )
        par << ", y_linked";

      par << "\n";
    }

    par << "}\n\n"
        << "marker\n"
        << "{\n";
    if( m_miss.size() == 1 )
      par << "  allele_missing = \"" << imap[0].missing_allele_name() << "\"\n";
    par << "  allele_delimiter = \"/\"\n"
        << "}\n";
  }
  else
  {
    par << "  marker_list, start = \"" << imap[0].name() << "\"";
    par << ", end = \"" << imap[imap.size()-1].name() << "\"";
    if( imap[0].get_model_type() == MLOCUS::X_LINKED )
      par << ", x_linked";
    else if( imap[0].get_model_type() == MLOCUS::Y_LINKED )
      par << ", y_linked";

    par << "\n}\n\n"
        << "marker\n"
        << "{\n"
        << "  allele_missing = \"" << imap[0].missing_allele_name() << "\"\n"
        << "  allele_delimiter = \"/\"\n"
        << "}\n";
  }

  // pedigree data
  //
  FPED::PedigreeConstIterator ped = fped.pedigree_begin();

  for( ; ped != fped.pedigree_end(); ++ped )
  {
    const FPED::FilteredPedigreeInfo& ped_info = ped->info();

    FPED::MemberConstIterator mem = ped->member_begin();

    for( ; mem != ped->member_end(); ++mem )
    {
      const FPED::SubpedigreeConstPointer sped = mem->subpedigree();

      out << ped->name() << "\t" << mem->name();

      if( mem->get_father() != NULL )
        out << "\t" << mem->get_father()->name();
      else
        out << "\t" << fped.info().individual_missing_code();

      if( mem->get_mother() != NULL )
        out << "\t" << mem->get_mother()->name();
      else
        out << "\t" << fped.info().individual_missing_code();

      if( mem->is_male() )
        out << "\t" << fped.info().sex_code_male();
      else if( mem->is_female() )
        out << "\t" << fped.info().sex_code_female();
      else
        out << "\t" << fped.info().sex_code_unknown();

      for( size_t t = 0; t < fped.info().trait_count(); ++t )
      {
        const RPED::RefTraitInfo& t_info = fped.info().trait_info(t);

        if(    t_info.name() == "SEX_CODE"
            || t_info.name() == "FAMILIAL_INDICATOR"
            || t_info.name() == "FOUNDER_INDICATOR"
            || t_info.name() == "PEDIGREE_SIZE" )
          continue;

        if( ped_info.trait_missing(mem->index(), t) )
          out << "\t" << fped.info().trait_info(t).string_missing_code();
        else
        {
          if( t_info.type() == RPED::RefTraitInfo::binary_trait )
          {
            if( ped_info.trait(mem->index(), t) == 1.0 )
              out << "\t" << t_info.string_affected_code();
            else
              out << "\t" << t_info.string_unaffected_code();
          }
          else
            out << "\t" << ped_info.trait(mem->index(), t);
        }
      }

      for( size_t s = 0; s < fped.info().string_count(); ++s )
      {
        out << "\t" << ped_info.get_string(mem->index(), s);
      }

      bool print_asis = false;

      if( sped_incon_map.find(sped) == sped_incon_map.end() )
        print_asis = true;

#if 0
      cout << ped->name() << ",s";
      if( sped != NULL )
        cout << setw(4) << sped->name();
      else
        cout << "    ";
      cout << ",m" << mem->name();
      if( print_asis )
        cout << " no incon";
      else
        cout << "    incon";
#endif

      for( size_t m = 0; m < imap.size(); ++m )
      {
        const MLOCUS::inheritance_model& model = imap[m];

        size_t p = ped_info.phenotype(mem->index(), m);

        if( !print_asis )
        {
          bool men_error = sped_incon_map[sped][m];

          if( men_error )
            p = model.get_missing_phenotype_id();
        }

        string pheno = get_pheno_string(p, model);
        out << "\t" << pheno;

#if 0
        cout << "\t" << p << "," << pheno;
#endif
      }

      out << endl;
#if 0
      cout << endl;
#endif
    }
  }
}

void
MARKERINFO::print_marker_pedigree_table(const inconsistency_handler& h, const FPED::Multipedigree& fped)
{
}

string
MARKERINFO::get_pheno_string(size_t p, const MLOCUS::inheritance_model& model) const
{
  string pheno_name = model.get_phenotype(p).name();

  if( p == model.get_missing_phenotype_id() )
  {
    string a = model.missing_allele_name();
    if( a == "*missing" )
      a = "?";

    pheno_name = a + "/" + a;
  }
  else
  {
    string al1, al2;
    MLOCUS::Ordering order;

    string geno_name = model.unphased_penetrance_begin(p).unphased_geno().name();

    model.gmodel().parse_genotype_name(geno_name, al1, al2, order);

    if( !al1.size() || al1.substr(0,1) == "~" ) al1 = model.missing_allele_name();
    if( !al2.size() || al2.substr(0,1) == "~" ) al2 = model.missing_allele_name();

    if( al1 == "*missing" )
      al1 = "?";

    if( al2 == "*missing" )
      al2 = "?";

    pheno_name = al1 + "/" + al2;
  }

  return pheno_name;
}

int main(int argc, char **argv)
{
  free(malloc(1));

  MARKERINFO *markerinfo = new MARKERINFO(argc,argv);
  assert(markerinfo != NULL);

  int r = markerinfo->main();

  delete markerinfo;

  return r;
}
