###################################################
#
# The purpose of this file is to define commonly
# shared data and functions for the build_scripts
# system.
#
###################################################

###################################################
# Exit status flags
###################################################

EXIT_SUCCESS = 0     # Everything went as expected
EXIT_DIFFS   = 1     # Everything compiled, but there were diffs
EXIT_FAILURE = 2     # Compile failed.

###################################################
# Imports
###################################################

import os
import time
import getopt
import sys

###################################################
# OS Options
###################################################

os_trans = { "Linux_ix86"       : 'gcc',
             "Linux_x86_64"     : 'gcc64',
             "Linux_ia64"       : 'gcc',
             "SunOS"            : 'gcc',
             "CYGWIN_NT-5.1"    : 'gcc',
             "Darwin_ppc"       : 'gcc',
             "Darwin_ix86"      : 'gcc' }

# for i in range(0, len(os.uname())):
#    print "--> os.uname: %s" % repr(os.uname()[i])

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

def get_compiler():
  return os_trans[get_os()]

def get_machine():
  return os.uname()[1]

osystem   = get_os()
compiler  = get_compiler()
QUIET     = 0

###################################################
# Time options
###################################################

# Determine time and date for filenames

date = time.strftime("%Y.%m.%d", time.localtime(time.time()))

# Check the current directory.  This script can only be run from a 'sage'
# directory or a 'sage/build_scripts' directory
def check_directory() :

  dir_split = os.path.split(sys.path[-1])
  par_dir   = os.path.split(dir_split[0])

#  print(dir_split[1])
#  print(par_dir[1])

  if(dir_split[1] != "build_scripts" or par_dir[1] != "util") :
    print "  This script must be run in the util/build_scripts directory."
    sys.exit(EXIT_FAILURE)

  os.chdir("../..")


def remove_text_in_file(fname, text):
  os.system('sed -e"s/%s//g" <%s >%s' % (text, fname, fname + '~~~'))
  os.system("mv %s %s" % (fname + '~~~', fname))

def remove_absolute_path(fname):

  dir = os.getcwd()

  dir = dir[:dir.find('src/c++')]

  dir = dir.replace('/', '\\/', 999)

  remove_text_in_file(fname, dir)

def remove_recompiling_messages_and_flags(fname):

  f = open(fname, 'r')
  
  l = f.readlines()
  
  f.close()

  g = open(fname + '~~~', 'w')

  for i in l :
    if i[0:4] != 'make'       and \
       i[0:9] != 'Compiling'  and \
       i[0:5] != '	Reco' and \
       i[0:5] != '	cd "' and \
       i[0:2] != 'CC'         and \
       i[0:3] != "CXX"        and \
       i[0:6] != "CFLAGS"     and \
       i[0:6] != "LDLIBS"     and \
       i[0:7] != "LDFLAGS"    and \
       i[0:8] != "CXXFLAGS"   and \
       i[0:8] != "  LDLIBS"   and \
       i[0:9] != "  LDFLAGS"  and \
       i[0:10]!= "  CXXFLAGS" :
      g.write(i)
  
  g.close()
  
  os.system("mv %s %s" % (fname + '~~~', fname))

def remove_test_timing(fname):

  f = open(fname, 'r')
  
  l = f.readlines()
  
  f.close()

  g = open(fname + '~~~', 'w')

  for i in l :
    if i[0:4] != 'Ran ':
      g.write(i)
  
  g.close()
  
  os.system("mv %s %s" % (fname + '~~~', fname))

import smtplib

def mail_file(filename, subject, tolist):
  # Only one email address works.
  tolist = "djb22@case.edu, yxs30@case.edu"
  if tolist == "":
    return

  #fromaddr = "sagebuilder-%s@darwin.case.edu" % (get_os())
  # Itentionally leaving as broken since this is needed
  #  for daily build test at night.
  fromaddr = "sagebuilder-%s@darwin" % (get_os())

  f = open(filename)

  msg = f.read()

  f.close()

  msg = "From: %s\r\nTo: %s\r\n" % (fromaddr, tolist) +  \
        "Subject: %s\r\n" % (subject) + \
        msg

  print
  print "  Mailing report to:",tolist

  server = smtplib.SMTP('smtp.cwru.edu')
  server.sendmail(fromaddr, tolist, msg)
  server.quit()

def usage():

  print
  print "Usage: %s <options> [run]" % sys.argv[0]
  print
  print "Options:"
  print "  -h, --help          = this message"
  print "  -q, --quiet         = run in quiet mode"
  print
  print "  -n, --notify (file) = use file of emails as list of people to email reports"
  print
  print "  [run]               = The command line run option for make"

def get_options():
  quiet = 0
  nlist = ""

  try:
    opts, args = getopt.getopt(sys.argv[1:], "hqn:", ["help", 'quiet', 'notify='])
  except getopt.GetoptError:
    # print help information and exit:
    usage()
    sys.exit(EXIT_FAILURE)

  for o, a in opts:
    if o in ("-h", "--help"):
      usage()
      sys.exit()
    if o in ("-q", "--quiet"):
      quiet = 1
    if o in ("-n", "--notify"):
      f=open(a,'r')

      nlist = f.readline()

      print "setting notification list to:", nlist

      f.close()

  return args, quiet, nlist
