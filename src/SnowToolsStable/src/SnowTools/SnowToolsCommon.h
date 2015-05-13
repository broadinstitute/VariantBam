#ifndef SNOWTOOLS_COMMON_H__
#define SNOWTOOLS_COMMON_H__

/*! \mainpage Snowtools 0.1
 *
 * \section intro_sec Introduction
 *
 * This is the introduction.
 *
 * \section install_sec Installation
 *
 * \subsection step1 Step 1: Opening the box
 *
 * etc...
 */

#include <string>
#include <vector>
#include <cstdint>

namespace SnowTools {
  
  static const std::vector<std::string> CHR_NAME {"1", "2", "3", "4", "5", "6", "7", "8", "9",
      "10", "11", "12", "13", "14", "15", "16", "17", 
      "18", "19", "20", "21", "22", "X", "Y", "M"};
  static const std::vector<std::string> CHR_NAME_NUM {"1", "2", "3", "4", "5", "6", "7", "8", "9",
      "10", "11", "12", "13", "14", "15", "16", "17", 
      "18", "19", "20", "21", "22", "23", "24"};
  
  static const int CHR_LEN [25] = {249250621, 243199373, 198022430, 191154276, //1-4
				   180915260, 171115067, //5-6
				   159138663, 146364022, 141213431, 135534747, 135006516, 133851895, //7-12
				   115169878, 107349540, 102531392, 90354753,  81195210,  78077248, //13-18
				   59128983,  63025520,  48129895,  51304566,  155270560, 59373566, //19-24
				   16571}; //25
  
  static const uint32_t CHR_CLEN [25] = {0, 249250621, 492449994,  690472424, 881626700, 1062541960, 1233657027,
					 1392795690,1539159712,1680373143,1815907890,1950914406,2084766301,
					 2199936179, 2307285719, 2409817111, 2500171864, 2581367074, 2659444322,                                                                                                                   
					 2718573305, 2781598825, 2829728720, 2881033286, 3036303846, 3095677412};
  
  static const uint32_t genome_size_XY = 3095677411;
  
  static std::string const REFHG19 = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";
  
  static const int NONCENT_CHR [44] = {1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,
				       11,11,12,12,13,14,15,16,16,17,17,18,18,19,19,20,20,21,21,22,23,23,24,24};
  
  
}

#endif
