#include "gelim/inconsistency_handler.h"

namespace SAGE
{

void inconsistency_printer::print_table
   (const handler& h) const
{
  handler::incon_family_iterator i = h.family_begin();

  for( ; i != h.family_end(); ++i)
  {
    const imodel_map& markers =
        i->second.begin()->first->multipedigree()->info().markers();

    const family_error_type& fam = i->second;

    family_error_type::const_iterator c = fam.begin();

    const ind_error_type& mother = *c;

    ++c;

    const ind_error_type& father = *c;

    ++c;

    // Note.  Each vector of errors is only large enough to store the
    // largest inheritance model which has an inconsistency.  However, since
    // parents are always labeled inconsistent, this is not a problem.
    //
    // XXX This is somewhat inefficient, as there may only be 1-2
    // inconsistencies with high indices.  This should be fixed.

    for(inconsistency_handler::error_map::const_iterator j = mother.second.begin(); j != mother.second.end(); ++j)
    {
      handler::error_type e = j->second;
      if(e != handler::none)
      {
        err << SAGE::priority(SAGE::error) 
            << "At Marker '" << markers.name(j->first) << "' parents '" 
            << mother.first->name() << "' and '"
            << father.first->name() << "' in pedigree '"
            << mother.first->pedigree()->name()
            << "' are inconsistent with the following child(ren) of their nuclear "
            << " family: ";

        bool comma = false;

        for(family_error_type::const_iterator child = c; child != fam.end(); ++child)
        {
          if(child->second.find(j->first) != child->second.end())
          {
            handler::error_type e = child->second.find(j->first)->second;
            if(e != handler::none)
            {
              if(comma) err << ", ";
              else      comma = true;
  
              err << child->first->name();
            }
          }
        }
      }
      err << "." << endl;
    }
  }
}

//
//--------------------------------------------------------------
//

void inconsistency_table_printer::print_pedigree_table
  (const handler& h, bool errors_only) const
{
  vector<ped_incon> ped_results;

  handler::incon_pedigree_iterator i = h.pedigree_begin();

  if(i == h.pedigree_end()) return;

  for( ; i != h.pedigree_end(); ++i)
  {
    size_t incon  = 0;
    size_t inform = 0;
    size_t total  = 0;
    
    for(size_t j = 0; j < i->second.checked.size(); ++j)
    {
      if(i->second.inconsistent[j]) ++incon;
      if(i->second.informative[j])  ++inform;
      if(i->second.checked[j])      ++total;
    }

    // Print marker count to be consistent with old version
    //
    //total = i->second.checked.size();

    ped_results.push_back(ped_incon(i->first->name(), incon, inform, total));
  }

  std::sort(ped_results.begin(), ped_results.end(), ped_incon_compare());

  err << "=====================================================================" << endl;

  err << "                                          Number of Markers" << endl
      << "Pedigree                          Incon.     with Data       Total " << endl
      << "----------                      ----------  -------------  ----------" << endl;

  for(size_t i = 0; i < ped_results.size(); ++i)
  {
    if( ped_results[i].inconsistent <= 0 && errors_only )
      continue;

    err << left  << setw(30) << ped_results[i].pedigree     << "  "
        << right << setw(8)  << ped_results[i].inconsistent << "    "
                 << setw(11) << ped_results[i].informative  << "    "
                 << setw(8)  << ped_results[i].total        << endl;
  }
  err << "=====================================================================" << endl;
}

void inconsistency_table_printer::print_marker_table
    (const handler& h, const imodel_map& imap, bool errors_only) const
{
  make_marker_list(h, imap, true);

  err << "=====================================================================" << endl;

  err << "                                         Number of Pedigrees" << endl
      << "Marker                            Incon.     with Data       Total " << endl
      << "----------                      ----------  -------------  ----------" << endl;

  for(size_t i = 0; i < marker_results.size(); ++i)
  {
    if( marker_results[i].inconsistent <= 0 && errors_only )
      continue;

    err << left  << setw(30) << imap.name(marker_results[i].marker)  << "  "
        << right << setw(8)  << marker_results[i].inconsistent       << "    "
                 << setw(11) << marker_results[i].informative        << "    "
                 << setw(8)  << marker_results[i].total              << endl;
  }
  err << "=====================================================================" << endl;
}

// This is the one MARKERINFO currently uses.
//
void inconsistency_table_printer::print_inconsistency_table
   (const handler& h, const imodel_map& imap,
    const vector<size_t>& sample_id, bool c_out, bool sort) const
{
  if(h.family_begin() == h.family_end()) return;

  print_header(h, c_out);

  make_marker_list(h, imap, sort);

  for( size_t c = 0; c < imap.size(); c += columns )
  {
    string marker_section_list;

    size_t end_col = c;

    for( ; end_col != c + columns && end_col < marker_results.size(); ++end_col )
    {
      if( marker_results[end_col].inconsistent <= 0 )
        continue;

      marker_section_list += imap.name(marker_results[end_col].marker);
      marker_section_list += "  ";
    }
    
    if( !marker_section_list.size() )
      continue;

    err << "=====================================================================";
    if( sample_id.size() )
    {
      for( size_t i = 0; i < sample_id.size(); ++i )
      {
        if( sample_id[i] != (size_t)-1 )
          err << "===========";
      }
    }
    err << endl;

    err << "Table " << (c/columns)+1 << " with Marker(s)  " << marker_section_list << endl;

    print_table_section(h, imap, c, end_col, sample_id, c_out);

    err << "=====================================================================";
    if( sample_id.size() )
    {
      for( size_t i = 0; i < sample_id.size(); ++i )
      {
        if( sample_id[i] != (size_t)-1 )
          err << "===========";
      }
    }
    err << endl << endl;
  }
}

void inconsistency_table_printer::print_inconsistency_table
   (const handler& h, const imodel_map& imap,
    size_type start_col, size_type end_col,
    const vector<size_t>& sample_id, bool c_out, bool sort) const
{
  if(h.family_begin() == h.family_end() || start_col <= end_col) return;

  make_marker_list(h, imap, sort);

  if(!marker_results[start_col].inconsistent) return;

  print_header(h, c_out);

  print_table_section(h, imap, start_col, end_col, sample_id, c_out);

  return;
}

bool
inconsistency_table_printer::print_inconsistency_member(stringstream& output,
                                                        const member_type& mem,
                                                        const vector<marker_incon>& incon_markers,
                                                        const imodel_map& imap,
                                                        size_type start_col, size_type end_col, size_t type,
                                                        const vector<size_t>& sample_id) const
{
  bool write_family = false;

  const FPED::FilteredPedigreeInfo& ped_info = mem.pedigree()->info();

  output << left << setw(10) << mem.pedigree()->name() << "  "
                 << setw(10) << mem.name() << " ";

  if( sample_id.size() )
  {
    for( size_t i = 0; i < sample_id.size(); ++i )
    {
      if( sample_id[i] != (size_t)-1 )
        output << setw(10) << ped_info.get_string(mem.index(), sample_id[i]) << " ";
    }
  }

  switch(type)
  {
    case 0 : output << "Mother "; break;
    case 1 : output << "Father "; break;
    default: output << "       "; break;
  }

  for(size_t j = start_col; j != end_col && j < incon_markers.size(); ++j)
  {
    const MLOCUS::inheritance_model& model = imap[incon_markers[j%columns].marker];
    
    uint p = model.get_missing_phenotype_id();

    if( !ped_info.phenotype_missing(mem.index(), incon_markers[j%columns].marker, model) )
      p = ped_info.phenotype(mem.index(), incon_markers[j%columns].marker);

    string pheno_name = model.get_phenotype(p).name();
   
    if(p == model.get_missing_phenotype_id())
    {
      string a = model.missing_allele_name();
      if( a == "*missing" )
        a = "?";

      pheno_name = a + model.gmodel().unphased_separator() + a;
    }
    else if(imap[incon_markers[j%columns].marker].codominant(p))
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

      pheno_name = al1 + model.gmodel().unphased_separator() + al2;
    }

    if( incon_markers[j%columns].inconsistent > 0 )
    {
      output << "* " << setw(7) << pheno_name << "     ";
      write_family = true;
    }
    else
      output << "  " << setw(12) << " ";
  }
  output << endl;

  return write_family;
}

//
//--------------------------------------------------------------
//

void inconsistency_table_printer::print_header
   (const handler& h, bool c_out) const
{
  err << "    * More than one nuclear family must be examined to detect" << endl
      << "      the inconsistency." << endl << endl;

  if( c_out)
    err << "   [] Indicates marker alleles that are not inconsistent." << endl << endl;

  if( h.is_sex_linked_error_exist() )
    err << "    # Sex-linked marker alleles inconsistent with the person's sex." << endl << endl;
}

void inconsistency_table_printer::print_table_section
   (const handler& h, const imodel_map& imap,
    size_type start_col, size_type end_col,
    const vector<size_t>& sample_id, bool c_out) const
{
  bool header = false;

  bool io_flag = err.flags() & ios::right;

  err << right;

  // To make sure to print out in order by pedigree name.
  //
  const handler::fam_list& incon_id_list = h.get_incon_family_list();

  for(handler::fam_iterator l = incon_id_list.begin(); l != incon_id_list.end(); ++l)
  {
    bool write_family = false;

    stringstream output;

    if(header) output << endl;

    // family_error_type
    // - list< pair< family_pointer, map<size_t, error_type> > >
    const family_error_type& fam = h.get_incon_family(**l);

    // Check to see if an inconsistency exist within family, then
    //  print all members of family when asked.

    vector<bool> print_worthy;

    if( c_out )
    {
      print_worthy.resize(columns, false);

      for(size_t j = start_col; j != end_col && j < marker_results.size(); ++j)
      {
        if( marker_results[j].inconsistent <= 0 )
          continue;

        family_error_type::const_iterator temp_fam_member = fam.begin();

        for( ; temp_fam_member != fam.end(); ++temp_fam_member )
        {
          if(    temp_fam_member->second.size()
              && temp_fam_member->second.find(marker_results[j].marker) != temp_fam_member->second.end() )
          {
            print_worthy[j%columns] = true;
            break;
          }
        }
      }
    }

    family_error_type::const_iterator fam_member = fam.begin();

    const FPED::FilteredPedigreeInfo& ped_info = fam_member->first->pedigree()->info();

    vector<marker_incon> mendelian_found;
    bool mendelian_exist = false;

    mendelian_found.resize(columns);
    for( size_t j = start_col; j != end_col && j%columns < mendelian_found.size(); ++j )
    {
      mendelian_found[j%columns].marker       = marker_results[j].marker;
      mendelian_found[j%columns].inconsistent = 0;
      mendelian_found[j%columns].informative  = 0;
      mendelian_found[j%columns].total        = 0;
    }

    for( size_t t = 0; fam_member != fam.end(); ++fam_member, ++t)
    {
      FPED::Multipedigree::member_const_pointer mem = fam_member->first;

      output << left << setw(10) << mem->pedigree()->name() << "  "
                     << setw(10) << mem->name() << " ";

      if( sample_id.size() )
      {
        for( size_t i = 0; i < sample_id.size(); ++i )
        {
          if( sample_id[i] != (size_t)-1 )
            output << setw(10) << ped_info.get_string(mem->index(), sample_id[i]) << " ";
        }
      }

      switch(t)
      {
        case 0 : output << "Mother "; break;
        case 1 : output << "Father "; break;
        default: output << "       "; break;
      }

      for(size_t j = start_col; j != end_col && j < marker_results.size(); ++j)
      {
        if( marker_results[j].inconsistent <= 0 )
          continue;

        const MLOCUS::inheritance_model& model = imap[marker_results[j].marker];

        uint p = model.get_missing_phenotype_id();

        if( !ped_info.phenotype_missing(mem->index(), marker_results[j].marker, model) )
          p = ped_info.phenotype(mem->index(), marker_results[j].marker);

        string pheno_name = model.get_phenotype(p).name();
       
        if(p == model.get_missing_phenotype_id())
        {
          string a = model.missing_allele_name();
          if( a == "*missing" )
            a = "?";

          pheno_name = a + model.gmodel().unphased_separator() + a;
        }
        else if(imap[marker_results[j].marker].codominant(p))
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

          pheno_name = al1 + model.gmodel().unphased_separator() + al2;
        }

        if( fam_member->second.size() )
        {
          if( fam_member->second.find(marker_results[j].marker) != fam_member->second.end() )
          {
            handler::error_type e = fam_member->second.find(marker_results[j].marker)->second;
            switch(e)
            {
              case handler::nuclear_family : output << "  " << setw(7) << pheno_name << "     "; write_family = true; break;
              case handler::mendelian      : output << "* " << setw(7) << pheno_name << "     "; write_family = true;
                                             ++mendelian_found[j%columns].inconsistent;
                                             mendelian_exist = true;

                                             break;

              case handler::xy_linked      : output << "# " << setw(7) << pheno_name << "     "; write_family = true; break;
              case handler::none           : SAGE_internal_error();
            }
          }
          else
            if( c_out && print_worthy[j%columns] )
              output << " [" << setw(7) << pheno_name << "]    ";
            else
              output << "  " << setw(12) << " ";
        }
        else
          if( c_out && print_worthy[j%columns] )
            output << " [" << setw(7) << pheno_name << "]    ";
          else
            output << "  " << setw(12) << " ";
      }
      output << endl;
    }

    // print the rest of pedigree members if mendelian inconsistency found.
    //
    if( mendelian_exist )
    {
      const subped_type& subped = *(fam.begin()->first->subpedigree());

      FPED::Multipedigree::family_const_iterator mfam = subped.family_begin();

      for( ; mfam != subped.family_end(); ++mfam )
      {
        if( (++(++fam.begin()))->first->family()->name() == mfam->name() )
          continue;

        if( mfam->parent1()->is_male() )
        {
          bool w_family = print_inconsistency_member(output, *(mfam->parent2()), mendelian_found,
                                                    imap, start_col, end_col, 0, sample_id);

               w_family = print_inconsistency_member(output, *(mfam->parent1()), mendelian_found,
                                                    imap, start_col, end_col, 1, sample_id);
        }
        else
        {
          bool w_family = print_inconsistency_member(output, *(mfam->parent1()), mendelian_found,
                                                    imap, start_col, end_col, 0, sample_id);

               w_family = print_inconsistency_member(output, *(mfam->parent2()), mendelian_found,
                                                    imap, start_col, end_col, 1, sample_id);
        }

        FPED::Multipedigree::offspring_const_iterator mem = mfam->offspring_begin();

        for( ; mem != mfam->offspring_end(); ++mem )
        {
          print_inconsistency_member(output, *mem, mendelian_found,
                                                    imap, start_col, end_col, 3, sample_id);
        }
      }
    }

    if(write_family)
    {
      if(!header)
      {
        print_table_header(h, imap, start_col, end_col, sample_id);
        header = true;
      }
      err << output.str();
    }
  }

  if(!io_flag)
    err << right;

  err << endl;

  return;
}

void inconsistency_table_printer::print_table_header
   (const handler& h, const imodel_map& imap,
    size_type start_col, size_type end_col, const vector<size_t>& sample_id) const
{
  err << "=====================================================================";
  if( sample_id.size() )
  {
    for( size_t i = 0; i < sample_id.size(); ++i )
    {
      if( sample_id[i] != (size_t)-1 )
        err << "===========";
    }
  }
  err << endl;

  err << "Pedigree    Individual ";
  if( sample_id.size() )
  {
    for( size_t i = 0; i < sample_id.size(); ++i )
    {
      err << left;

      if( sample_id[i] != (size_t)-1 )
        err << setw(11)
            << h.get_incon_family_list().front()->multipedigree()->info().string_info(sample_id[i]).name();
    }
  }
  err << "      ";

  for(size_t j = start_col; j != end_col && j < marker_results.size(); ++j)
  {
    if( marker_results[j].inconsistent <= 0 )
      continue;

    err << right << setw(11) << imap.name(marker_results[j].marker)
        << "   ";
  }
  err << endl;

  err << "----------  ---------- ";
  if( sample_id.size() )
  {
    for( size_t i = 0; i < sample_id.size(); ++i )
    {
      if( sample_id[i] != (size_t)-1 )
        err << "---------- ";
    }
  }
  err << "      ";

  for(size_t j = start_col; j != end_col && j < marker_results.size(); ++j)
  {
    if( marker_results[j].inconsistent <= 0 )
      continue;
    err << "-----------   ";
  }
  err << endl;

  return;
}

void inconsistency_table_printer::make_marker_list
    (const handler& h, const imodel_map& imap, bool sort) const
{
  marker_results.resize(0);

  size_t msize = 0;

  handler::incon_pedigree_iterator i = h.pedigree_begin();

  for( ; i != h.pedigree_end(); ++i)
    if(msize < i->second.checked.size())
      msize = i->second.checked.size();

  marker_results.resize(msize);
  
  for(uint j = 0; j < msize; ++j)
    marker_results[j].marker = j;

  // Count up inconsistencies and such.

  i = h.pedigree_begin();
  for( ; i != h.pedigree_end(); ++i)
  {
    for(size_t j = 0; j < i->second.checked.size(); ++j)
    {
      if(i->second.inconsistent[j]) ++marker_results[j].inconsistent;
      if(i->second.informative[j])  ++marker_results[j].informative;
      if(i->second.checked[j])      ++marker_results[j].total;
    }
  }

  // Sort the markers.

  if(sort)
  {
    std::sort(marker_results.begin(), marker_results.end(), marker_incon_compare());
  }
  else
  {
    // Compress all the bad markers to the beginning.

    size_t i = 0;

    for( ; i < imap.size() && marker_results[i].inconsistent; ++i);

    size_t j = i;

    for(++i; i < imap.size(); ++i)
    {
      if(marker_results[i].inconsistent)
      {
        std::swap(marker_results[j], marker_results[i]);

        ++j;
      }
    }
  }

  return;
}

}
