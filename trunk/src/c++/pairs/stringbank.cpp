//****************************************************************************
//* File:      relnamebank.h                                                 *
//*                                                                          *
//* Author:    Kevin Jacobs & Yeunjoo Song                                   *
//*                                                                          *
//* History:   Version 1.0                                                   *
//*                                                                          *
//* Notes:     Maps conventional names onto relationship objects.            *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include <string>
#include "pairs/stringbank.h"

namespace SAGE
{

const char* ordinal_names[] =
        { "zero",
          "one", 
          "two",
          "three",
          "four",
          "five",
          "six",
          "seven",
          "eight",
          "nine",
          "ten",
          "eleven",
          "twelve",
          "thirteen",
          "fourteen",
          "fifteen",
          "sixteen",
          "seventeen",
          "eighteen",
          "nineteen",
          "twenty"
        };

const int max_ordinal_names = 20;
     
const char* ordinal_th_names[] =
        { "",        // zero-th?
          "first",
          "second",
          "third",
          "fourth",
          "fifth",
          "sixth",
          "seventh",
          "eighth",
          "nineth",
          "tenth",
          "eleventh",
          "twelfth",
          "thirteenth",
          "fourteenth",
          "fifteenth",
          "sixteenth",
          "seventeenth",
          "eitheenth",
          "nineteenth",
          "twentieth"
       };

const int max_ordinal_th_names = 20;

const char* ordinal_tens_names[] =
       { "",    // tens
         "twenty",
         "thirty",
         "fourty",
         "fifty",
         "sixty",
         "seventy",
         "eighty",
         "ninety"
       };

const int max_ordinal_tens_names = 8;

const char* ordinal_power_names[] =
       { "ones",
         "tens",
         "hundred",
         "thousand",
         "million",
         "billion",
         "trillion",
         "quadrilian",
         "quintillian",
         "sextillion",
         "septillion",
         "octillion",
         "nonillian",
         "decillian"
       };

const int max_ordinal_power_names = 13;

const char* ordinal_st_names[] = 
       { "",   // zero-th
         "once",
         "twice",
         "thrice",
         "times"
       };

const int max_ordinal_st_names = 4;

const char* relative_words[] =
       {
         "male",		//[0]
         "female",		//[1]
         "unknown_gender",	//[2]
         "paternal",		//[3]
         "maternal",		//[4]
         "unknown_phase",	//[5]
         "self",		//[6]
         "father",		//[7]
         "mother",		//[8]
         "parent",		//[9]
         "son",			//[10]
         "daughter",		//[11]
         "offspring",		//[12]
         "brother",		//[13]
         "sister",		//[14]
         "mztwin",		//[15]
         "dztwin",		//[16]
         "half",		//[17]
         "grand",		//[18]
         "great",		//[19]
         "uncle",		//[20]
         "aunt",		//[21]
         "avuncular",		//[22]
         "nephew",		//[23]
         "niece",		//[24]
         "nephewship",		//[25]
         "cousin",		//[26]
         "removed",		//[27]
         "none",		//[28]
         "invalid",		//[29]
         "sibling",		//[30]
         "marital",		//[31]
         "through",		//[32]
         "'s",    		//[33]
         "removed"    		//[34]
       };

const int max_relative_words = 34;

} // end of namespace SAGE
