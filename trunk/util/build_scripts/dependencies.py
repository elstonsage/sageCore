#!/usr/local/bin/python

import os
import string

def generate_file_list(cdir):
  return "../include/" + cdir + "/*.h "   + \
         "../include/" + cdir + "/*.ipp " + \
         "./"                        + "*.cpp"

         
#determine where we are

(main_dir, current_dir) = os.path.split(os.path.abspath('./'))

# generate the set of files we'll be looking at.  Made as a function since
# it will likely be changed/enhanced later, but a simple command line is ok for now

file_list = generate_file_list(current_dir) 

# Generate the list of includes.  Note that it doesn't do any parsing, so it only
# finds includes at the beginning of a line, and doesn't care if that include is
# in a comment or preprocessed out.  This may need modified, but it's good enough
# for now.

os.system('grep "^#[ \\t]*include " ' + file_list + ' | ' + \
          'sed "s/[.\/0-9a-zA-Z_]*.[hcip]*://g"       | ' + \
          'sed "s/#[ \\t]*include[ \\t]*//g"            ' + \
          ' > include_files~~~')

# Now we have to read the file to pull out the appropriate text bits.

f = open('include_files~~~')

text = f.readlines()

tokens = {}

# For each line in the file, determine if it's a module.  IF so, store the module
# for later work.

for i in text:
  # Determine if there's a directory attached.  If so, we need to figure it out
  
  if string.find(i,'/') != -1:
    i=i[1:string.find(i,'/')]
  
    tokens[i] = i

f = open('include_files~~~', 'w')
for i in tokens :
  f.write(i + '\n')
