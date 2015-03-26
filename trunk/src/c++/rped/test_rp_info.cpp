#include "rped/rped.h"

namespace SAGE {
namespace RPED {

RefMarkerInfo create_marker(MLOCUS::GenotypeModelType mt)
{
  MLOCUS::genotype_model g("M1");
  
  g.set_model_type(mt);
  
  g.add_allele("A",0.5);
  g.add_allele("B",0.5);
  
  g.set_missing_allele_name("~M");
  
  RefMarkerInfo marker(g, true, true);
  
  return marker;
}
  
void test_ped_info_set_phenotypes()
{
    RefMPedInfo mpinfo;
    
    mpinfo.add_marker("M1");
    
    RefPedInfo pedinfo;
  
    MPED::pedigree_base   P("P6");

    P.add_member("m01");
    P.add_member("m02");
    P.add_member("m03");
    P.add_member("m04");
    P.add_member("m05");
    P.add_member("m06");
    P.add_lineage("m03", "m02", "m01");
    P.add_lineage("m04", "m01", "m02");
    P.add_lineage("m05", "m01");
    P.add_lineage("m05", "m02");
    P.add_lineage("m06", "m01");
    P.add_lineage("m06", "m02");

    P.build();
  
    pedinfo.build(P);
    
    pedinfo.resize_markers(1, mpinfo);

    // Basic setting
    assert(pedinfo.set_phenotype(0, 0, 0) == true);
    assert(pedinfo.set_phenotype(0, 1, 0) == false);
    assert(pedinfo.set_phenotype(5, 0, 0) == true);
    assert(pedinfo.set_phenotype(6, 0, 0) == false);
    
    // Setting using first marker
    RefMarkerInfo marker = create_marker(MLOCUS::AUTOSOMAL);
    
    assert(pedinfo.set_phenotype(0, 0, "A", "A", marker) == 0);
    assert(pedinfo.set_phenotype(0, 1, "A", "A", marker) == 3);
    assert(pedinfo.set_phenotype(5, 0, "A", "A", marker) == 0);
    assert(pedinfo.set_phenotype(6, 0, "A", "A", marker) == 3);

    assert(pedinfo.set_phenotype(5, 0, "A/B",   "",  marker) == 0);
    assert(pedinfo.set_phenotype(5, 0, "A/B",   "A", marker) == 2);
    assert(pedinfo.set_phenotype(5, 0, "A/C",   "",  marker) == 2);
    assert(pedinfo.set_phenotype(5, 0, "",      "",  marker) == 1);

    assert(pedinfo.set_phenotype(0, 0, "A",  "",   marker) == 2);
    assert(pedinfo.set_phenotype(0, 0, "A",  "~M", marker) == 2);
    assert(pedinfo.set_phenotype(0, 0, "~M", "~M", marker) == 1);

    assert(pedinfo.set_phenotype(5, 0, "*missing", "",  marker) == 1);
    
    // Try dynamics
    marker.gmodel().set_dynamic_alleles(true);

    assert(pedinfo.set_phenotype(0, 0, "A", "D", marker) == 0); // Adds D
    assert(pedinfo.set_phenotype(0, 0, "E", "D", marker) == 0); // Adds E
    assert(pedinfo.set_phenotype(0, 0, "F", "G", marker) == 0); // Adds F and G
    
    marker.gmodel().set_dynamic_alleles(false);

    // X-linked
    marker = create_marker(MLOCUS::X_LINKED);
    
    // Female
    assert(pedinfo.set_phenotype(0, MPED::SEX_FEMALE, 0, "A", "A", marker) == 0);
    assert(pedinfo.set_phenotype(0, MPED::SEX_FEMALE, 1, "A", "A", marker) == 3);
    assert(pedinfo.set_phenotype(5, MPED::SEX_FEMALE, 0, "A", "A", marker) == 0);
    assert(pedinfo.set_phenotype(6, MPED::SEX_FEMALE, 0, "A", "A", marker) == 3);

    assert(pedinfo.set_phenotype(5, MPED::SEX_FEMALE, 0, "A/B", "",  marker) == 0);
    assert(pedinfo.set_phenotype(5, MPED::SEX_FEMALE, 0, "A/B", "A", marker) == 2);
    assert(pedinfo.set_phenotype(5, MPED::SEX_FEMALE, 0, "A/C", "",  marker) == 2);
    assert(pedinfo.set_phenotype(5, MPED::SEX_FEMALE, 0, "",    "",  marker) == 1);

    assert(pedinfo.set_phenotype(0, MPED::SEX_FEMALE, 0, "A", "", marker) == 2);

    assert(pedinfo.set_phenotype(0, MPED::SEX_FEMALE, 0, "A",  "",   marker) == 2);
    assert(pedinfo.set_phenotype(0, MPED::SEX_FEMALE, 0, "A",  "~M", marker) == 5);
    assert(pedinfo.set_phenotype(0, MPED::SEX_FEMALE, 0, "~M", "~M", marker) == 1);

    // Male
    assert(pedinfo.set_phenotype(0, MPED::SEX_MALE, 0, "A", "A", marker) == 0);
    assert(pedinfo.set_phenotype(0, MPED::SEX_MALE, 1, "A", "A", marker) == 3);
    assert(pedinfo.set_phenotype(5, MPED::SEX_MALE, 0, "A", "A", marker) == 0);
    assert(pedinfo.set_phenotype(6, MPED::SEX_MALE, 0, "A", "A", marker) == 3);

    assert(pedinfo.set_phenotype(5, MPED::SEX_MALE, 0, "A/B", "",  marker) == 5);
    assert(pedinfo.set_phenotype(5, MPED::SEX_MALE, 0, "A/B", "A", marker) == 2);
    assert(pedinfo.set_phenotype(5, MPED::SEX_MALE, 0, "A/C", "",  marker) == 2);
    assert(pedinfo.set_phenotype(5, MPED::SEX_MALE, 0, "",    "",  marker) == 1);
    
    assert(pedinfo.set_phenotype(0, MPED::SEX_MALE, 0, "A", "", marker) == 2);

    assert(pedinfo.set_phenotype(0, MPED::SEX_MALE, 0, "A",  "",   marker) == 2);
    assert(pedinfo.set_phenotype(0, MPED::SEX_MALE, 0, "A",  "~M", marker) == 0);
    assert(pedinfo.set_phenotype(0, MPED::SEX_MALE, 0, "~M", "~M", marker) == 1);

    // Male vs. Female
    assert(pedinfo.set_phenotype(0, MPED::SEX_FEMALE, 0, "A", "A", marker) == 0);
    assert(pedinfo.set_phenotype(1, MPED::SEX_MALE, 0, "A", "A", marker) == 0);
    
    assert(pedinfo.phenotype(0, 0) != pedinfo.phenotype(1,0));
    
    // Y-linked
    marker = create_marker(MLOCUS::Y_LINKED);
    
    assert(pedinfo.set_phenotype(0, MPED::SEX_FEMALE, 0, "~X", "~X", marker) == 0);
    assert(pedinfo.set_phenotype(0, MPED::SEX_FEMALE, 0, "A", "A",   marker) == 5);
    assert(pedinfo.set_phenotype(0, MPED::SEX_FEMALE, 1, "A", "A",   marker) == 3);
    assert(pedinfo.set_phenotype(5, MPED::SEX_FEMALE, 0, "A", "A",   marker) == 5);
    assert(pedinfo.set_phenotype(6, MPED::SEX_FEMALE, 0, "A", "A",   marker) == 3);

    assert(pedinfo.set_phenotype(5, MPED::SEX_FEMALE, 0, "A/B", "",  marker) == 2);
    assert(pedinfo.set_phenotype(5, MPED::SEX_FEMALE, 0, "A/B", "A", marker) == 2);
    assert(pedinfo.set_phenotype(5, MPED::SEX_FEMALE, 0, "A/C", "",  marker) == 2);
    assert(pedinfo.set_phenotype(5, MPED::SEX_FEMALE, 0, "",    "",  marker) == 1);

    assert(pedinfo.set_phenotype(0, MPED::SEX_FEMALE, 0, "A", "", marker) == 2);

    assert(pedinfo.set_phenotype(0, MPED::SEX_FEMALE, 0, "A",  "",   marker) == 2);
    assert(pedinfo.set_phenotype(0, MPED::SEX_FEMALE, 0, "A",  "~M", marker) == 5);
    assert(pedinfo.set_phenotype(0, MPED::SEX_FEMALE, 0, "~M", "~M", marker) == 1);

    assert(pedinfo.set_phenotype(0, MPED::SEX_MALE, 0, "~X", "~X", marker) == 5);
    assert(pedinfo.set_phenotype(0, MPED::SEX_MALE, 0, "A", "A",   marker) == 0);
    assert(pedinfo.set_phenotype(0, MPED::SEX_MALE, 1, "A", "A",   marker) == 3);
    assert(pedinfo.set_phenotype(5, MPED::SEX_MALE, 0, "A", "A",   marker) == 0);
    assert(pedinfo.set_phenotype(6, MPED::SEX_MALE, 0, "A", "A",   marker) == 3);

    assert(pedinfo.set_phenotype(5, MPED::SEX_MALE, 0, "A/B", "",  marker) == 2);
    assert(pedinfo.set_phenotype(5, MPED::SEX_MALE, 0, "A/B", "A", marker) == 2);
    assert(pedinfo.set_phenotype(5, MPED::SEX_MALE, 0, "A/C", "",  marker) == 2);
    assert(pedinfo.set_phenotype(5, MPED::SEX_MALE, 0, "",    "",  marker) == 1);
    
    assert(pedinfo.set_phenotype(0, MPED::SEX_MALE, 0, "A", "", marker) == 2);

    assert(pedinfo.set_phenotype(0, MPED::SEX_MALE, 0, "A",  "",   marker) == 2);
    assert(pedinfo.set_phenotype(0, MPED::SEX_MALE, 0, "A",  "~M", marker) == 0);
    assert(pedinfo.set_phenotype(0, MPED::SEX_MALE, 0, "~M", "~M", marker) == 1);

}

}
}

int main()
{
  SAGE::RPED::test_ped_info_set_phenotypes();
  
  
  return 0;
}
