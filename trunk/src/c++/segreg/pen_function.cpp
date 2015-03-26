#include "segreg/pen_function.h"

namespace SAGE
{
namespace SEGREG
{

void pg::dump_post_genos(ostream& o, const post_geno_map& pens)
{
  string ped_name = "";

  for(post_geno_map::const_iterator p  = pens.begin(); p != pens.end(); ++p)
  {
    member_pointer         m = p->first;
    const member_pen_list& l = p->second;

    // If the pedigree name is different from ped_name, we've begun a new
    // pedigree.
    if(m->pedigree()->name() != ped_name)
    {
      if(p != pens.begin())
        o << "  }" << endl;

      ped_name = m->pedigree()->name();

      o << "  pedigree=\"" << ped_name << "\"" << endl;
      o << "  {" << endl;
    }

    o << "    individual=\"" << m->name() << "\", type=";

    if(MPED::mp_utilities::is_founder(m))
      o << "founder" << endl;
    else 
      o << "nonfounder" << endl;

    o << "    {" << endl;

    for(member_pen_list::const_iterator i = l.begin(); i != l.end(); ++i)
    {
      o << "      genotype=\"" << i->geno << "\", penetrance=" << i->pen 
        << ", posterior_genotype=" << i->post_geno << endl;
    }

    o << "    }" << endl;
  }

  o << "  }" << endl; // End of last pedigree block
}

void pf::dump_pen_func(ostream& o, const pen_func_map& pens)
{
  string ped_name = "";

  for(pen_func_map::const_iterator p  = pens.begin(); p != pens.end(); ++p)
  {
    member_pointer         m = p->first;
    const member_pen_list& l = p->second;

    // If the pedigree name is different from ped_name, we've begun a new
    // pedigree.
    if(m->pedigree()->name() != ped_name)
    {
      if(p != pens.begin())
        o << "  }" << endl;

      ped_name = m->pedigree()->name();

      o << "  pedigree=\"" << ped_name << "\"" << endl;
      o << "  {" << endl;
    }

    o << "    individual=\"" << m->name() << "\", type=";

    if(MPED::mp_utilities::is_founder(m))
      o << "founder" << endl;
    else 
      o << "nonfounder" << endl;

    o << "    {" << endl;

    for(member_pen_list::const_iterator i = l.begin(); i != l.end(); ++i)
    {
      o << "      genotype=\"" << i->geno;

      if(!MPED::mp_utilities::is_founder(m))
      {
        o << "\", mother_genotype=\"" << i->mo_geno
          << "\", father_genotype=\"" << i->fa_geno;
      }

      o << "\", penetrance=" << i->pen
        << endl;
    }

    o << "    }" << endl;
  }

  o << "  }" << endl; // End of last pedigree block
}

}}
