#ifndef __LSF_TYPES_H
#define __LSF_TYPES_H

#define LSF_MASK    0xFFFF0000
#define LOCAL_MASK  0x0000FFFF
#define LSF_SHIFT       16
#define LOCAL_SHIFT      0

#define LSF_NONE         0
#define LSF_BASE         1
#define LSF_COMPONENT    2
#define LSF_COMPOSITE    3
#define LSF_REF          4
#define LSF_MAP          5
#define LSF_FACTORY      6
#define LSF_INT          7
#define LSF_REAL         8
#define LSF_STRING       9
#define LSF_ITER        10
#define LSF_GUARD       11
#define LSF_EXPR        12
#define LSF_SWITCH      13

#define LSF_DISPLAY   1000

struct Type_mapping      
{ 
  unsigned long value;
  char *name;
};

struct LSF_mapping      
{ 
  unsigned long value;
  char *name;
};

struct Local_mapping      
{ 
  unsigned long value;
  char *name;
};

extern LSF_mapping lsf_mappings[];

#endif
