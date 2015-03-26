#ifndef LODPAL_UTIL_H
#define LODPAL_UTIL_H

#include "maxfun/maxfun.h"
#include "maxfunapi/maxfunapi.h"
#include "numerics/matrix.h"
#include "numerics/sinfo.h"
#include "numerics/print_util.h"
#include "palbase/pal_ibd.h"
#include "palbase/pair_filter.h"
#include "palbase/pair_info_file.h"

namespace SAGE   {
namespace LODPAL {

typedef MAXFUN::ParameterMgr            maxfun_param_mgr;
typedef MAXFUN::SequenceCfg             maxfun_seq_cfg;
typedef MAXFUN::RunCfg                  maxfun_run_cfg;
typedef MAXFUN::Maximizer               maxfun_maximizer;
typedef MAXFUN::DebugCfg                maxfun_debug;
typedef MAXFUN::Results                 maxfun_results;

typedef PALBASE::rel_pair_data          rel_pair_data;
typedef PALBASE::relative_pairs         RelativePairs;
typedef PALBASE::rel_pair               rel_pair;
typedef PALBASE::pair_filter            pair_filter;

} // end of namespace LODPAL
} // end of namespace SAGE

#endif
