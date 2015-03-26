#include "gelim/geno_eliminate.h"

namespace SAGE
{

genotype_eliminator::genotype_eliminator()
    : err()
{ }

bool
genotype_eliminator::set_subpedigree(const subpedigree_type& sped)
{
  subpedigree = &sped;

  return true;
}

size_t
genotype_eliminator::process(imodel& model, size_t marker, removal_type remove, bool do_x_data_check)
{
#if 0
  cout << "process()..." << endl;
#endif
  // Set the initial type of inconsistency.

  //error_type e_type = handler::nuclear_family;
  error_type e_type = handler::mendelian;

  // If either set isn't empty, we operate on those people only.

#if 0
  cout << endl
       << "m = " << marker
       << ", set1 size = " << set1.size()
       << ", set2 size = " << set2.size()
       << endl;
#endif

  if(set1.empty())
  {
    if(!set2.empty())
    {
      set1.swap(set2);

      e_type = handler::mendelian;
    }
    else
      generate_list(model);
  }
  else
    e_type = handler::mendelian;

#if 0
  cout << endl
       << "m = " << marker
       << ", set1 size = " << set1.size()
       << ", set2 size = " << set2.size()
       << endl;
#endif


  bool inconsistent = false;

  // How the algorithm works:
  //
  // There are two sets of families.  Initially, set1 contains the set we
  // want to check, while set2 is empty.  Whenever a family is checked, for
  // each individual in that family whose genoset is altered by the
  // elimination, the nuclear families to which that individual belongs are
  // removed from set1 (if present) and added to set2 (if propagation is
  // active).  This delays the evaluation of those nuclear families until
  // set1 has emptied.  This has the effect of allowing more information to
  // propagate to those families from other places in the subpedigree, meaning
  // that those families that are delayed will have more information when
  // they are eventually processed.  This reduces the number of times
  // individual families are processed.
  //
  // When set1 is empty, if there are delayed families in set2, the sets
  // are swapped, and processing continues as before.  When both sets are
  // empty, we have obtained our final state and the algorithm can exit.

  while(!set1.empty())
  {
    const family_type& family = *set1.back();   // Get the last family from the list
    set1.pop_back();

    size_t error_code = 0;

    error_code = process_family(model, marker, family, e_type, remove, true);

    if( error_code == 3 )
    {
      inconsistent = true;
    }

    if( set1.empty() )
    {
      set1.swap(set2);  // Move the new set of families in.
      e_type = handler::mendelian;
    }
  }

  err.mark_info(subpedigree->pedigree(), marker, inconsistent, true);
#if 0
  cout << "end of process()..." << endl;
#endif
  if(inconsistent) return 2; // Inconsistent
  
  return 0;                  // Everything is fine.
}

void dump_parental_genotypes(const FPED::Family& fam,
                             const valid_parental_genotypes& parent_genotypes)
{
  // Test the parental genotypes
  cout << endl << "dump_parental_genotypes..." << endl
       << fam.parent1()->name() << " X "
       << fam.parent2()->name() << endl;

  cout << parent_genotypes.genotype_pair_count()  << " M: "
       << parent_genotypes.mother_valid_genotype_count()   << ", "
       << parent_genotypes.mother_invalid_genotype_count() << " F: "
       << parent_genotypes.father_valid_genotype_count()   << ", "
       << parent_genotypes.father_invalid_genotype_count() << endl;

  valid_parental_genotypes::genotype_pair_iterator i = 
          parent_genotypes.genotype_pair_begin();

  for( ; i != parent_genotypes.genotype_pair_end(); ++i)
  {
    cout << i->first << ' ' << i->second << endl;
  }

  valid_parental_genotypes::genotype_iterator j =
          parent_genotypes.mother_valid_genotype_begin();

  cout << "MV:  ";

  for( ; j != parent_genotypes.mother_valid_genotype_end(); ++j)
  {
    cout << (*j) << ' ';
  }

  cout << endl << "MI:  ";

  j = parent_genotypes.mother_invalid_genotype_begin();

  for( ; j != parent_genotypes.mother_invalid_genotype_end(); ++j)
  {
    cout << (*j) << ' ';
  }

  j = parent_genotypes.father_valid_genotype_begin();

  cout << endl << "FV:  ";

  for( ; j != parent_genotypes.father_valid_genotype_end(); ++j)
  {
    cout << (*j) << ' ';
  }

  cout << endl << "FI:  ";

  j = parent_genotypes.father_invalid_genotype_begin();

  for( ; j != parent_genotypes.father_invalid_genotype_end(); ++j)
  {
    cout << (*j) << ' ';
  }

  cout << endl;
}

void dump_all_parental_genotypes(const FPED::Family&    fam,
                                 const valid_parental_genotypes::parental_genotype_pair_vector& parent_genotypes)
{
  // Test the parental genotypes
  cout << endl << "dump_all_parental_genotypes..." << endl
       << fam.parent1()->name() << " X "
       << fam.parent2()->name() << endl;

  cout << parent_genotypes.size() << endl;

  valid_parental_genotypes::genotype_pair_iterator i = parent_genotypes.begin();

  for( ; i != parent_genotypes.end(); ++i)
  {
    cout << i->first << ' ' << i->second << endl;
  }
}

/// Processes a single nuclear family.

/// process_family() works in two passes.  In the first pass, the valid
/// parental genotypes are generated.  A pair of parental genotypes is
/// considered valid if for every child there is at least one genotype which
/// is valid for the child and consistent with the pair of parental
/// genotypes.  This pass finds the genotypes that are valid for the parents
/// given the children.
///
/// In the second pass, if the parents are invalid, the child genotypes are
/// generated as if all pairs were valid (to find invalid children that
/// caused parents to be invalid).  If the parents are valid, only genotypes
/// using valid parental genotype pairs are used.  This pass removes child
/// genotypes that are not consistent with the parental genotypes and their
/// siblings.
  
size_t
genotype_eliminator::process_family(imodel& model, size_t marker, const family_type& fam,
                                    error_type e, removal_type remove, bool prop)
{
  // Build our data structures.  If this doesn't work, returns an error. 

  size_t error = build_family_data(model, fam);
#if 0
  cout << endl
       << "start of process_family()..." << fam.name()
       << endl << "error from build_family_data() = " << error << endl;
#endif
  if(error) 
    return error;  // If state 1 or 2 (see header)

  // Begin genotype elimination

  // First pass - Generate valid genotypes of parents

  par_genotypes parent_genotypes;

  generate_valid_parental_genotypes(model, parent_genotypes);

#if 0
  dump_parental_genotypes(fam, parent_genotypes);
#endif

  // Second pass - generate valid states for children by either:
  //   A: All parent pairs (if parents already invalid)
  //   B: Only valid parent pairs (if parents still have valid genotypes)
  
  // Determine if there were valid parental genotype pairs.
  if( parent_genotypes.genotype_pair_count() )
  {
    // If pairs exist, we only check the valid pairs against the children.
    // We know there aren't any errors at this point, because the children
    // have all been checked, but we want to eliminate invalid states.

    generate_valid_child_genotypes(model,
                                   parent_genotypes.genotype_pair_begin(),
                                   parent_genotypes.genotype_pair_end());
  }
  else
  {
    // If there aren't any valid parental genotype pairs, we know that we're
    // inconsistent, but we still want to check the children to determine
    // which, if any, specific children might be the culprit.  We'd like to
    // narrow it down to the smallest number possible.

    par_pair_vector all_pair_vector;

    generate_all_parental_genotypes(model, all_pair_vector);

    generate_valid_child_genotypes(model, all_pair_vector.begin(), all_pair_vector.end()); 

#if 0
  dump_all_parental_genotypes(fam, all_pair_vector);
#endif
  }

  // Produce any errors that might have been found

  bool error_added = generate_errors(model, parent_genotypes, marker, e);

  if( error_added )
    err.mark_info(subpedigree->pedigree(), marker, true, true);

  // Remove invalid genotypes

  switch(remove)
  {
    case genotype: 
      remove_genotypes(model, parent_genotypes, prop);
      break;

    case all:
      remove_all(model, parent_genotypes, prop);
      break;

    case none:
    default:
      break;
  }

#if 0
  //model.print_info_sparse_matrix();
  cout << "end of process_family()..." << endl;
#endif

  // If there weren't any valid parental pairs, then we know the subpedigree is
  // inconsistent.
  if(parent_genotypes.genotype_pair_count() == 0)
  {
    return 3;
  }

  return 0;
}

// ----------------------------------------

//
//-----------------------------------------
// 

void genotype_eliminator::generate_list(imodel& model)
{
  // We want set1 = all nuclear families with genetic information.
  //         set2 = empty.

  if(!set1.empty()) set1 = fam_set();
  if(!set2.empty()) set2 = fam_set();     // Make sure the sets are empty

  subpedigree_type::family_const_iterator fam     = subpedigree->family_begin();
  subpedigree_type::family_const_iterator fam_end = subpedigree->family_end();

  for( ; fam != fam_end; ++fam)
  {
    // If there's genetic information, add it to set1
    if(informative_family(model, *fam))
      set1.push_back( &*fam );
  }
}

void genotype_eliminator::generate_valid_parental_genotypes(imodel& model, par_genotypes& pg)
{
  pg.generate_valid_parental_genotypes(*my_child_info[0].ptr->family(), model);
}

void genotype_eliminator::generate_all_parental_genotypes(imodel& model, par_pair_vector& dest)
{
  unphased_iterator gm_begin = model.unphased_penetrance_begin (mother_ptr->subindex()+1);
  unphased_iterator gm_end   = model.unphased_penetrance_end   (mother_ptr->subindex()+1);

  unphased_iterator gf_begin = model.unphased_penetrance_begin (father_ptr->subindex()+1);
  unphased_iterator gf_end   = model.unphased_penetrance_end   (father_ptr->subindex()+1);

  for(unphased_iterator gm = gm_begin; gm != gm_end; ++gm)
  {
    for(unphased_iterator gf = gf_begin; gf != gf_end; ++gf)
    {
      dest.push_back(make_pair(gm.geno_id(), gf.geno_id()));
    }
  }
}

void
genotype_eliminator::generate_valid_child_genotypes(imodel& model,
                                                    par_pair_iterator begin,
                                                    par_pair_iterator end)
{
  for(par_pair_iterator p_genotypes = begin; p_genotypes != end; ++p_genotypes)
  {
    MLOCUS::unphased_genotype gm = model.get_unphased_genotype(p_genotypes->first);
    MLOCUS::unphased_genotype gf = model.get_unphased_genotype(p_genotypes->second);

    MLOCUS::child_genotype_set cg(gm, gf, false);

    child_vector::iterator child_iter;   // Iterator over the children

    for(child_iter = my_child_info.begin(); child_iter != my_child_info.end(); ++child_iter)
    {
      phased_iterator    gc = child_iter->begin;
      bit_type::iterator bc = child_iter->data.begin();

      for( ; gc != child_iter->end; ++bc, ++gc)
        if(cg.contains(gc.phased_geno()))
        {
          *bc = true;
        }
    }
  }

#if 0
  child_vector::iterator child_iter;
  for(child_iter = my_child_info.begin(); child_iter != my_child_info.end(); ++child_iter)
  {
    cout << "child " << child_iter->ptr->name()
         << ", bit size = " << child_iter->data.size() << endl;

    for( size_t i = 0; i < child_iter->data.size(); ++i )
      cout << " " << child_iter->data[i];
    cout << endl;

    phased_iterator    gc = child_iter->begin;
    for( ; gc != child_iter->end; ++gc )
      cout << " " << gc.geno_id() << "," << gc.phenotype_id();
    cout << endl;
  }
#endif
}

bool
genotype_eliminator::generate_errors(imodel& model, const par_genotypes& pg,
                                     size_t marker, error_type e)
{
#if 0
  cout << "start of generate_errors()..." << endl;
#endif
  // Check for errors

  bool errors_detected = false;

  // Child Errors
  child_vector::iterator child;         // Iterator over the children

  for(child = my_child_info.begin(); child != my_child_info.end(); ++child)
  {
#if 0
  cout << child->ptr->name() << " "
       << model.phased_penetrance_count(child->ptr->subindex()+1) << ", "
       << child->data.empty()
       << endl;
#endif
    if(model.phased_penetrance_count(child->ptr->subindex()+1) && child->data.empty())
    {
      errors_detected = true;

      err.add_error(*child->ptr->family(), child->ptr, marker, e);
#if 0
  cout << endl
       << "#ped = " << subpedigree->pedigree()->name()
       << ", sib = " << my_child_info.begin()->ptr->index()
       << ", id = " << child->ptr->name()
       << ", m = " << marker << endl;
#endif
    }
  }

  // Parental Errors
  if(!errors_detected && !pg.father_valid_genotype_count())
  {
    errors_detected = true;

    err.add_error(*my_child_info.begin()->ptr->family(), marker, e);

#if 0
  cout << endl
       << "$ped = " << subpedigree->pedigree()->name()
       << ", sib = " << my_child_info.begin()->ptr->name()
       << ", m = " << marker << endl;
#endif
  }

#if 0
  cout << "end of generate_errors()..." << endl;
#endif

  return errors_detected;
}

void
genotype_eliminator::remove_genotypes(imodel& model, const par_genotypes& pg, bool prop)
{
  if(pg.mother_invalid_genotype_count())
  {
    par_iterator mother_begin = pg.mother_invalid_genotype_begin();
    par_iterator mother_end   = pg.mother_invalid_genotype_end();

    remove_genotypes(model, *mother_ptr, mother_begin, mother_end);

    if(prop) parent_set_move(*mother_ptr, *father_ptr);
  }
  
  if(pg.father_invalid_genotype_count())
  {
    par_iterator father_begin = pg.father_invalid_genotype_begin();
    par_iterator father_end   = pg.father_invalid_genotype_end();

    remove_genotypes(model, *father_ptr, father_begin, father_end);

    if(prop) parent_set_move(*mother_ptr, *father_ptr);
  }
  
  child_vector::iterator child;   // Iterator over the children

  for(child = my_child_info.begin(); child < my_child_info.end(); ++child)
  {
    if(child->data.size())
    {
      bool child_genotypes_removed = remove_genotypes(model, *child);

      if(child_genotypes_removed && prop)
      {
        child_set_move(*child->ptr, *child->ptr);
      }
    }
  }
}

void
genotype_eliminator::remove_all(imodel& model, const par_genotypes& pg, bool prop)
{
  if(!pg.mother_valid_genotype_count())
  {
    par_iterator mother_begin = pg.mother_invalid_genotype_begin();
    par_iterator mother_end   = pg.mother_invalid_genotype_end();

    remove_genotypes(model, *mother_ptr, mother_begin, mother_end);

    if(prop) parent_set_move(*mother_ptr, *father_ptr);
  }
  
  if(!pg.father_valid_genotype_count())
  {
    par_iterator father_begin = pg.father_invalid_genotype_begin();
    par_iterator father_end   = pg.father_invalid_genotype_end();

    remove_genotypes(model, *father_ptr, father_begin, father_end);

    if(prop) parent_set_move(*mother_ptr, *father_ptr);
  }

  child_vector::iterator child;   // Iterator over the children

  for(child = my_child_info.begin(); child < my_child_info.end(); ++child)
  {
    if(child->data.empty() && child->data.size())
    { 
      bool child_genotypes_removed = remove_genotypes(model, *child);

      if(child_genotypes_removed && prop)
      {
        child_set_move(*child->ptr, *child->ptr);
      }
    }
  }
}

void
genotype_eliminator::parent_set_move(const individual_type& i, const individual_type& mate)
{
  if(i.family()) // if there are parents
  {
    move_family(*i.family());
  }

  child_set_move(i, mate);
}

void
genotype_eliminator::child_set_move(const individual_type& index, const individual_type& mate)
{
  individual_type::mate_const_iterator mt     = index.mate_begin();
  individual_type::mate_const_iterator mt_end = index.mate_end();

  for( ; mt != mt_end; ++mt)
    if(&mt->mate() != &mate)
    {
      move_family(mt->family());
    }
}

void
genotype_eliminator::move_family(const family_type& fam)
{
  fam_set::iterator i = set1.begin();

  for( ; i != set1.end(); ++i)
    if(*i == &fam) break;

  if(i != set1.end()) set1.erase(i);

  set2.push_back(&fam);
}

/// A nuclear family is informative if any individual within is individually informative.

bool
genotype_eliminator::informative_family(imodel& model, const family_type& i) const
{
#if 0
  cout << "start of informative_family()..." << i.name() << endl;
#endif
  size_t al_count = model.allele_count();

  size_t phased_size   = al_count * al_count;
  size_t unphased_size = al_count * (al_count - 1) / 2;

  size_t parent1 = i.parent1()->subindex();

  size_t ind_genotypes_p1 = model.unphased_penetrance_count(parent1+1);

#if 0
  cout << i.parent1()->name() << " p1 ind_genotypes = " << ind_genotypes_p1 << " "
       << "unphased_size = " << unphased_size << endl;
#endif

  if( ind_genotypes_p1 != unphased_size ) return true;

  size_t parent2 = i.parent2()->subindex();

  size_t ind_genotypes_p2 = model.unphased_penetrance_count(parent2+1);

#if 0
  cout << i.parent2()->name() << " p2 ind_genotypes = " << ind_genotypes_p2 << " "
       << "unphased_size = " << unphased_size << endl;
#endif

  if( ind_genotypes_p2 != unphased_size ) return true;

  // Check children
  family_type::offspring_const_iterator off     = i.offspring_begin();
  family_type::offspring_const_iterator off_end = i.offspring_end();

  for( ; off != off_end; ++off)
  {
    size_t ind_genotypes_c = model.phased_penetrance_count(off->subindex() + 1);

#if 0
  cout << off->name() << " off ind_genotypes = " << ind_genotypes_c << " "
       << "phased_size = " << phased_size << endl;
#endif

    if( ind_genotypes_p1 == 1 && ind_genotypes_p2 == 1 && ind_genotypes_c == phased_size )
      return true;

    if( ind_genotypes_c != phased_size ) return true;
  }

#if 0
  cout << "end of informative_family()..." << endl;
#endif
  return false;
}

//
//-----------------------------------------------
//

void
genotype_eliminator::remove_genotypes(imodel& model, const individual_type& par,
                                      par_iterator begin, par_iterator end)
{
  for(par_iterator i = begin; i != end; ++i)
    model.remove_unphased_penetrance(par.subindex()+1, *i);
}

bool
genotype_eliminator::remove_genotypes(imodel& model, const child_info_type& child)
{
  phased_iterator i = model.phased_penetrance_begin(child.ptr->subindex()+1);

  vector<int> genotypes(child.data.size());

  for(uint c = 0; i != model.phased_penetrance_end(child.ptr->subindex()+1); ++i, ++c)
  {
    genotypes[c] = i.geno_id();
  }

  bool elim = false;

  bit_field::const_iterator b = child.data.begin();

  for(uint c = 0; b != child.data.end(); ++b, ++c)
  {
    if(*b == 0)
    {
      model.remove_phased_penetrance(child.ptr->subindex()+1, genotypes[c]);

      elim = true;
    }
  }

  return elim;
}

} // End namespace SAGE

#undef __geno

