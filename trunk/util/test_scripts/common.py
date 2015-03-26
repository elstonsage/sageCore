import os

###################################################
# OS Options
###################################################

def get_os() :
  osname = os.uname()[0]
  
  if osname == "Linux":
    if os.uname()[4] == "x86_64" :
      return "Linux_x86_64"
    elif os.uname()[4] == "ia64" :
      return "Linux_ia64"
    else:
      return "Linux_ix86"
  elif osname == "Darwin":
    if os.uname()[4] == "i386" :
      return "Darwin_ix86"
    else:
      return "Darwin_ppc"

  return osname

