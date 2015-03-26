#include "freq/Sample.h"

namespace SAGE {
namespace FREQ {

//======================================================
//
//  Sample CONSTRUCTOR
//
//======================================================
Sample::Sample(const FPED::FilteredMultipedigree & mp, size_t marker)
  :
  my_mp        (mp),
  my_marker_id (marker)
{
  // Populate the mate const itr vectors:
  
  my_mate_begin_itrs .resize(my_mp.member_count());
  my_mate_end_itrs   .resize(my_mp.member_count());
  
  for(size_t i = 0; i < my_mp.member_count(); ++i)
  {
    my_mate_begin_itrs [i] = my_mp.member_index(i).mate_begin ();
    my_mate_end_itrs   [i] = my_mp.member_index(i).mate_end   ();
  }

  // Resize the master vector:
  my_genotype_infos.resize(my_mp.member_count());
  
  // Clear the valid subpedigree list and unconnecteds list:
  my_valid_speds  . clear();
  my_unconnecteds . clear();

  // Create the pedigree imodel generator and set it up:
  pedigree_imodel_generator generator;
  
  generator.set_prior_remap          (true);
  generator.set_post_remap           (true);
  generator.set_genotype_elimination (true);

  // Populate sample:
  
  // Loop across pedigrees:
  for(FPED::PedigreeConstIterator ped = my_mp.pedigree_begin(); ped != my_mp.pedigree_end(); ++ped)
  {
    // Loop across unconnecteds:
    for(FPED::MemberConstIterator unc = ped->unconnected_begin(); unc != ped->unconnected_end(); ++unc)
    {
      // If this individual has a missing phenotype, skip him!
      if(!phenotypePresent(*unc))
        continue;
        
      // Set up validity:
      bool valid = false;

      // Loop across genotypes:
      for(MLOCUS::unphased_genotype_iterator genotype = getMarkerInfo().unphased_genotype_begin();
          genotype != getMarkerInfo().unphased_genotype_end(); ++genotype)
      {
        double penetrance = getMarkerInfo().unphased_penetrance(getPhenotype(*unc), genotype->get_id());
          
        // If the penetrance is non-zero, add this genotype to the individual's genotype vector:
        if(penetrance)
        {
          GenotypeInfo g_info;
            
          g_info.genotype   = *genotype;
          g_info.penetrance = penetrance;
            
          my_genotype_infos[unc->mpindex()].push_back(g_info);
          
          // Validate the sample:
          valid = true;
        };
          
      } // End of genotype loop
      
      // Add this dude to the list:
      if(valid)
      {
        my_unconnecteds.push_back(&*unc);
      }

    } // End of loop across unconnecteds
    
    for(FPED::SubpedigreeConstIterator sped = ped->subpedigree_begin(); sped != ped->subpedigree_end(); ++sped)
    {
      // Create an imodel for this subpedigree:
      MLOCUS::inheritance_model imodel = generator(*sped, my_marker_id);
      
      // If there are any inconsistencies, throw out this subpedigree:
      if(generator.inconsistent())
      {
        continue;
      }
      else // Otherwise, add this subpedigree to the list of valids:
      {
        SpedInfo sped_info;

        sped_info.sped             = &*sped;
        sped_info.remapped_gmodel  =   imodel.gmodel();
        sped_info.allele_remapping .   setup(getMarkerInfo().gmodel(), imodel.gmodel());
        
        my_valid_speds.push_back(sped_info);
      }
        
      // Loop across all individuals:
      for(FPED::MemberConstIterator ind = sped->member_begin(); ind != sped->member_end(); ++ind)
      {
        // Loop across genotypes:
        for(MLOCUS::unphased_genotype_iterator genotype = imodel.unphased_genotype_begin();
            genotype != imodel.unphased_genotype_end(); ++genotype)
        {
          double penetrance = imodel.unphased_penetrance(ind->subindex() + 1, genotype->get_id());
          
          // If the penetrance is non-zero, add this genotype to the individual's genotype vector:
          if(penetrance)
          {
            GenotypeInfo g_info;
            
            g_info.genotype   = *genotype;
            g_info.penetrance = penetrance;
            
            my_genotype_infos[ind->mpindex()].push_back(g_info);
          };
          
        } // End of genotype loop
        
      } // End of member loop

    } // End subpedigree loop

  } // End pedigree loop
}

//======================================================
//
//  Sample COPY CONSTRUCTOR
//
//======================================================
Sample::Sample(const Sample & other) :
  my_mp               (other.my_mp),
  my_marker_id        (other.my_marker_id),
  my_genotype_infos   (other.my_genotype_infos),
  my_valid_speds      (other.my_valid_speds),
  my_unconnecteds     (other.my_unconnecteds),
  my_mate_begin_itrs  (other.my_mate_begin_itrs),
  my_mate_end_itrs    (other.my_mate_end_itrs)
{ }

//======================================================
//
//  dump()
//
//======================================================
void Sample::dump() const
{
  std::cout << "SAMPLE DUMP"               << std::endl
            << "==========="               << std::endl
            <<                                std::endl
            << "Marker = " << my_marker_id << std::endl
            <<                                std::endl;

  // Dump singletons:

  if(getUnconnecteds().size())
  {
    OUTPUT::Table st("Singletons");
  
    st << OUTPUT::TableColumn("mpindex") << OUTPUT::TableColumn("ind name") << OUTPUT::TableColumn("phenotype");
      
    for(MLOCUS::unphased_genotype_iterator genotype = getMarkerInfo().gmodel().unphased_genotype_begin();
        genotype != getMarkerInfo().gmodel().unphased_genotype_end(); ++genotype)
    {
      st << OUTPUT::TableColumn(genotype->allele1().name() + "/" + genotype->allele2().name());
    }

    for(MemberVector::const_iterator ind = getUnconnecteds().begin(); ind != getUnconnecteds().end(); ++ind)
    {
      OUTPUT::TableRow r;

      r << (*ind)->mpindex() << (*ind)->name() << getPhenotype(**ind);
    
      for(MLOCUS::unphased_genotype_iterator genotype = getMarkerInfo().unphased_genotype_begin();
          genotype != getMarkerInfo().unphased_genotype_end(); ++genotype)
      {
        bool g_present = false;
      
        for(GenotypeInfoVector::const_iterator g_itr = getGenotypeInfoVector(**ind).begin(); g_itr != getGenotypeInfoVector(**ind).end(); ++g_itr)
          if((*genotype) == g_itr->genotype)
            g_present = true;
      
        r << (g_present ? "X" : "");

      } // End individual genotype loop
      
      st << r;

    } // End individual loop

    std::cout << st;
  }
  
  if(getValidSpedInfos().size())
  {
    // Dump subpedigrees:

    for(SpedInfoVector::const_iterator sped_itr = getValidSpedInfos().begin(); sped_itr != getValidSpedInfos().end(); ++sped_itr)
    {
      OUTPUT::Table t(sped_itr->sped->member_begin()->pedigree()->name());
      
      t << OUTPUT::TableColumn("mpindex") << OUTPUT::TableColumn("ind name") << OUTPUT::TableColumn("phenotype");
      
      for(MLOCUS::unphased_genotype_iterator genotype = sped_itr->remapped_gmodel.unphased_genotype_begin();
          genotype != sped_itr->remapped_gmodel.unphased_genotype_end(); ++genotype)
      {
        t << OUTPUT::TableColumn(genotype->allele1().name() + "/" + genotype->allele2().name());
      }

      for(FPED::MemberConstIterator ind = sped_itr->sped->member_begin(); ind != sped_itr->sped->member_end(); ++ind)
      {
        OUTPUT::TableRow r;

        r << ind->mpindex() << ind->name();
      
        if(phenotypePresent(*ind))
        {
          r << getPhenotype(*ind);
        }
        else
        {
          r << "<Missing>";
        }
    
        for(MLOCUS::unphased_genotype_iterator genotype = sped_itr->remapped_gmodel.unphased_genotype_begin();
            genotype != sped_itr->remapped_gmodel.unphased_genotype_end(); ++genotype)
        {
          bool g_present = false;
      
          for(GenotypeInfoVector::const_iterator g_itr = getGenotypeInfoVector(*ind).begin(); g_itr != getGenotypeInfoVector(*ind).end(); ++g_itr)
            if(*genotype == g_itr->genotype)
              g_present = true;
      
          r << (g_present ? "X" : "");

        } // End individual genotype loop
      
        t << r;

      } // End individual loop

      std::cout << t;

      sped_itr->allele_remapping.dump();

    } // End subpedigree loop
    
  } // End if-has-valid-spedinfos
}

} // End namespace FREQ
} // End namespace SAGE
