#define protected public

#include "mlocus/imodel.h"

namespace SAGE   {
namespace MLOCUS {

void inheritance_model_test(ostream& o, const penetrance_model& m, const char* str)
{
    uint        diff;
    uint        cnt;
    bool        good;

    o << "************************************" << endl;
    o << "Testing penetrance model " << str << endl;

    genotype_model_test(o, m.gmodel());

    //lint -e{56} -e{732}
    diff  = m.phenotype_end() - m.phenotype_begin();
    cnt = m.phenotype_count();
    if (cnt == diff)
        o << "..........phenotype count OK" << endl;
    else
        o << "..........phenotype count FAIL" << endl;

    penetrance_model::phenotype_iterator  pf = m.phenotype_begin();       
    penetrance_model::phenotype_iterator  pl = m.phenotype_end();     

    good = true;

    for (uint i = 0;  pf != pl;  ++pf, ++i)
    {
        if (*pf != m.get_phenotype(pf->id()))
        {
            o << "..........phenotype index FAIL at: " << i << endl;
            good = false;
        }

        if (pf->id() != m.get_phenotype_id(pf->name()))
        {
            o << ".........phenotype lookup FAIL at: " << i << endl;
            o << "  pf->id: " << pf->id() << "  lookup->id:" << m.get_phenotype_id(pf->name()) << endl;
            good = false;
        }

        if (m.get_phenotype(pf->name()).id() != m.get_phenotype_id(pf->name()))
        {
            o << ".........phenotype lookup FAIL at: " << i << endl;
            good = false;
        }
    }
    if (good)
    {
        o << ".........phenotype lookup OK" << endl;
    }
    o << endl;
}


void inheritance_model_print(ostream& o, const penetrance_model& m, const char* str)
{
    //lint --e{522}

    o << "************************************" << endl;
    o << "Printing penetrance model " << str << endl;
    o << "  symmetric: " << boolalpha << m.phmodel().has_generated_phenotypes() << endl;
    o << " codominant: " << boolalpha << m.codominant(true) << endl;
    
    if(m.codominant(true) != m.codominant(false))
      o << "lcodominant: " << boolalpha << m.codominant(false) << endl;

    genotype_model_print(o, m.gmodel());

    o << " phenotypes: ";
    penetrance_model::phenotype_iterator  pf = m.phenotype_begin();
    penetrance_model::phenotype_iterator  pl = m.phenotype_end();
    for (;  pf != pl;  ++pf)
    {
        o << " " << pf->name();
    }
    o << endl << endl;

    uint i;

    //lint -e{534}
    o.setf(ios::left);
 
    o << "unphased penetrance matrix:" << endl << "           ";
    unphased_genotype_iterator ugf = m.unphased_genotype_begin();
    unphased_genotype_iterator ugl = m.unphased_genotype_end();
    for (;  ugf != ugl;  ++ugf)
    {
        o << " " << (*ugf).name();
    }
    o << endl;
    for (i = 0;  i <= m.phenotype_count();  ++i)
    {
        o << setw(12) << left << m.get_phenotype(i).name().c_str();
        ugf = m.unphased_genotype_begin();
        for (;  ugf != ugl;  ++ugf)
        {
            o << setprecision(1);
            o << setw(3) << m.unphased_penetrance(i, ugf->get_id()) << " ";
        }
        o << endl;
    }
    o << endl;

    o << "phased penetrance matrix:" << endl << "           ";
    phased_genotype_iterator pgf = m.phased_genotype_begin();
    phased_genotype_iterator pgl = m.phased_genotype_end();
    for (;  pgf != pgl;  ++pgf)
    {
        o << " " << (*pgf).name();
    }
    o << endl;
    for (i = 0;  i <= m.phenotype_count();  ++i)
    {
        o << setw(12) << m.get_phenotype(i).name().c_str(); 
        for (pgf = m.phased_genotype_begin();  pgf != pgl;  ++pgf)
        {
            o << setprecision(1);
            o << setw(3) << m.phased_penetrance(i, pgf->get_id()) << " ";
        }
        o << endl;
    }
    o << endl;

    penetrance_model::unphased_penetrance_iterator    udf, udl;
    o << "unphased penetrance matrix iteration:" << endl;
    for (i = 0;  i <= m.phenotype_count();  ++i)
    {
        o << setw(12) << left << m.get_phenotype(i).name().c_str() << right;
        udf = m.unphased_penetrance_begin(i);
        udl = m.unphased_penetrance_end(i);
        for (;  udf != udl;  ++udf)
        {
            o << "  " << *udf << " [" << udf.unphased_geno().name() << "]";
        }
        o << endl;
    }
    o << endl;

    penetrance_model::phased_penetrance_iterator      pdf, pdl;
    o << "phased penetrance matrix iteration:" << endl;
    for (i = 0;  i <= m.phenotype_count();  ++i)
    {
        o << setw(12) << left << m.get_phenotype(i).name().c_str() << right;
        pdf = m.phased_penetrance_begin(i);
        pdl = m.phased_penetrance_end(i);
        for (;  pdf != pdl;  ++pdf)
        {
            o << "  " << *pdf << " [" << pdf.phased_geno().name() << "]";
        }
        o << endl;
    }
    o << endl;
}

} // End namespace MLOCUS
} // End namespace SAGE

