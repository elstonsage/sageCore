#ifndef AO_MODEL_H
#define AO_MODEL_H
//======================================================================
//
//  File:	Model.h
//
//  Author:	Stephen Gross
//
//  Copyright 2002 R. C. Elston
//======================================================================


#include <vector>
#include "LSF/parse_ops.h"
#include "mped/mp.h"
#include "mped/sp.h"
#include "rped/rped.h"
#include "sampling/sampling.h"
#include "maxfun/maxfun.h"
#include "maxfun/transf_sub_model.h"
#include "maxfunapi/maxfunapi.h"
#include "ageon/Datatypes.h"
#include "ageon/Validator.h"
#include "ageon/PoolingCfg.h"

namespace SAGE {
namespace AO   {

//======================================================================
//
// class Model
//
//======================================================================
class Model
{
  public:

    ///
    /// \internal
    /// @name Index lookup variables
    /// These variables store the index numbers for parameter lookup (in Maxfunapi)
    /// and Field lookup (in MemberDataSample). They are populated with values when
    /// the respective components are setup.
    //@{
    
    int lambda1_mxid;
    int lambda2_mxid;
    
    SAMPLING::FieldID AO_id;    // age of onset id
    SAMPLING::FieldID AE_id;    // age at exam id
    SAMPLING::FieldID aff_id;   // affectedness id
    SAMPLING::FieldID class_id; // classification id
    
    //@}

  public:

    //==================================================================
    // Constructors & operators:
    //==================================================================

    Model            ();
    Model            (const Model &);
    Model& operator= (const Model &);
    void      copy   (const Model &);

    //==================================================================
    // Public utility functions:
    //==================================================================

    void verify                        ();
    void reset                         ();
    void update_after_parse            ();
    int  update                        (vector<double>     & params, 
                                        int                  t);
    void setupSample                   (SAMPLING::PartitionedMemberDataSample & sample);
    size_t getClassType                (const SAMPLING::PartitionedMemberDataSample & sample,
                                        const FPED::Member & ind);
    void calculate_initial_est         (const SAMPLING::PartitionedMemberDataSample & sample,
                                        int                  t);
    void add_trait                     (trait_type           t, 
                                        string               name,
                                        bool                 new_initial_est,
                                        double               initial_est,
                                        bool                 new_fixed,
                                        bool                 fixed);
    void translate_maxfun_output (const Maxfun_Data        & data,
                                        int                  t);
    void dump_model(ostream& o) const;

    //==================================================================
    // Public accessors:  
    //==================================================================

    const string                     & get_title                 () const;
    const string                     & get_affectedness_trait    () const;
    const string                     & get_age_of_onset_trait    () const;
    const string                     & get_age_of_exam_trait     () const;
    const string                     & get_class_trait           () const;

          double                       get_epsilon               () const;
          double                       get_min_denominator       () const;
          size_t                       get_num_of_classes        () const;
          string                       get_ofilename             () const;
          bool                         is_valid                  () const;
          bool                         get_truncate              () const;
          bool                         get_adjustment            () const;
          bool                         get_allow_averaging       () const;
          bool                         get_pool_class            () const;
	  bool                         get_debug                 () const;

    const map<size_t, size_t>&         get_class_type_map        () const;

    const PoolingCfg & get_pooling_cfg() const { return my_pooling_cfg; }

    const MAXFUN::DebugCfg           & getDebugCfg               () const;

    const MAXFUN::ParameterMgr       & GetParameterMgr           () const;

    //==================================================================
    // Public mutators:
    //==================================================================

          string                     & title                  ();
          string                     & affectedness_trait     ();
          string                     & age_of_onset_trait     ();
          string                     & age_of_exam_trait      ();
          string                     & class_trait            ();

          double                     & epsilon                ();
          double                     & min_denominator        () const;
          size_t                     & num_of_classes         ();
          string                     & ofilename              ();
          bool                       & validity               ();
          bool                       & truncate               ();
          bool                       & adjustment             ();
          bool                       & allow_averaging        ();
          bool                       & pool_class             ();
	  bool                       & debug                  ();
	  
	  PoolingCfg & get_pooling_cfg() { return my_pooling_cfg; }

	  MAXFUN::ParameterMgr       & GetParameterMgr        ();

          MAXFUN::DebugCfg           & getDebugCfg            ();
  private:
    //==================================================================
    // Data members:
    //==================================================================

	    // Basic information about trait names:

            string                     my_title;
            string                     my_affectedness_trait;
            string                     my_age_of_onset_trait;
            string                     my_age_of_exam_trait;
            string                     my_class_trait;

            // Parameter data:

		MAXFUN::ParameterMgr my_parameter_mgr;

		MAXFUN::DebugCfg my_debug_cfg;

            // Assorted variables:

            transformation_sub_model   my_transf_sub_model;
            PoolingCfg                 my_pooling_cfg;
	    string                     my_ofilename;
    mutable double                     my_min_denominator;
            double                     my_epsilon;
            size_t                     my_num_of_classes;
            bool                       my_valid;
            bool                       my_truncate;
            bool                       my_use_adjustment;
            bool                       my_allow_averaging;
            bool                       my_pool_class;

	    bool                       my_debug;

    map<size_t, size_t>                my_class_type_map;
};


// Accessors:

//=================================================================================================================//
//                                                                                                                 //
//                                                  TRAIT_NAMES                                                    //
//                                                                                                                 //
inline const std::string                & Model::get_title              () const { return my_title;              }
inline const std::string                & Model::get_affectedness_trait () const { return my_affectedness_trait; }
inline const std::string                & Model::get_age_of_onset_trait () const { return my_age_of_onset_trait; }
inline const std::string                & Model::get_age_of_exam_trait  () const { return my_age_of_exam_trait;  }
inline const std::string                & Model::get_class_trait        () const { return my_class_trait;        }
//                                                                                                                 //
//=================================================================================================================//
//                                                                                                                 //
//                                                  MISCELLANY                                                     //
//                                                                                                                 //
inline const MAXFUN::ParameterMgr         & Model::GetParameterMgr          () const { return my_parameter_mgr; }
inline const MAXFUN::DebugCfg             & Model::getDebugCfg              () const { return my_debug_cfg;     }

//
inline double Model::get_epsilon               () const { return my_epsilon;               }
inline double Model::get_min_denominator       () const { return my_min_denominator;       }
inline size_t Model::get_num_of_classes        () const { return my_num_of_classes;        }
inline string Model::get_ofilename             () const { return my_ofilename;             }
inline bool   Model::is_valid                  () const { return my_valid;                 }
inline bool   Model::get_truncate              () const { return my_truncate;              }
inline bool   Model::get_adjustment            () const { return my_use_adjustment;        }
inline bool   Model::get_allow_averaging       () const { return my_allow_averaging;       }
inline bool   Model::get_pool_class            () const { return my_pool_class;            }
inline bool   Model::get_debug                 () const { return my_debug;                 }
inline const map<size_t, size_t>& Model::get_class_type_map() const { return my_class_type_map; }

//                                                                                                                 //
//=================================================================================================================//

// Mutators:

//=================================================================================================================//
//                                                                                                                 //
//                                                  TRAIT_NAMES                                                    //
//                                                                                                                 //
inline       std::string                & Model::title                  ()       { return my_title;              }
inline       std::string                & Model::affectedness_trait     ()       { return my_affectedness_trait; }
inline       std::string                & Model::age_of_onset_trait     ()       { return my_age_of_onset_trait; }
inline       std::string                & Model::age_of_exam_trait      ()       { return my_age_of_exam_trait;  }
inline       std::string                & Model::class_trait            ()       { return my_class_trait;        }
//                                                                                                                 //
//=================================================================================================================//
//                                                                                                                 //
//                                                  MISCELLANY                                                     //
//                                                                                                                 //
inline       double                     & Model::epsilon                ()       { return my_epsilon;            }
inline       double                     & Model::min_denominator        () const { return my_min_denominator;    }
inline       size_t                     & Model::num_of_classes         ()       { return my_num_of_classes;     }
inline       string                     & Model::ofilename              ()       { return my_ofilename;          }
inline       bool                       & Model::validity               ()       { return my_valid;              }
inline       bool                       & Model::truncate               ()       { return my_truncate;           }
inline       bool                       & Model::adjustment             ()       { return my_use_adjustment;     }
inline       bool                       & Model::allow_averaging        ()       { return my_allow_averaging;    }
inline       bool                       & Model::pool_class             ()       { return my_pool_class;         }
inline       MAXFUN::ParameterMgr       & Model::GetParameterMgr        ()       { return my_parameter_mgr;      }
inline       MAXFUN::DebugCfg           & Model::getDebugCfg            ()       { return my_debug_cfg;          }
inline       bool                       & Model::debug                  ()       { return my_debug;              }
//                                                                                                                 //
//=================================================================================================================//

}} // End namespace

#endif
