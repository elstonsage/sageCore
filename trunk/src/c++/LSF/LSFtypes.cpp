#include <string>
#include "LSF/LSF.h"
#include "LSF/LSFtypes.h"

LSF_mapping lsf_mappings[] =
{ { LSF_NONE,      "None"      },  /* LSF objects: formal names */
  { LSF_BASE,      "base"      },
  { LSF_COMPONENT, "item"      },
  { LSF_COMPOSITE, "list"      },
  { LSF_REF,       "ref"       },
  { LSF_MAP,       "map"       },
  { LSF_FACTORY,   "factory"   },
  { LSF_STRING,    "string"    },  /* LSF storage types */
  { LSF_INT,       "int"       },
  { LSF_REAL,      "real"      },
  { LSF_GUARD,     "guard"     },  /* LSF interpreted types */
  { LSF_SWITCH,    "switch"    },
  { LSF_EXPR,      "expr"      },
  { LSF_EXPR,      "expression"},
  { LSF_ITER,      "iterator"  },
  { (unsigned long) -1,     0           }  };

