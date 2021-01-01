/** @file */
/** 
 * This module contains all the other functions which are common for the library and do not fit elsewhere 
 **/

#ifndef OTHERFUNCTIONS_H
#define OTHERFUNCTIONS_H

#include "ANNIncludes.h"



/// Converts a given datatype to a std::string;
template<typename var>
std::string varToString (var input){
    std::stringstream ss;
    ss << input;
    std::string output;
    output= ss.str();
    return output;
};

/// Converts argument number from std::string format to integer 
int stringToInt(std::string input);

/// Converts argument number from string format to floating point single precision
float stringToFloat(std::string input);

/// Converts argument number from string format to floating pointdouble precision 
double stringToDouble(std::string input);

/// Returns 1 if the file exists else 0	
bool fileExists(const std::string& filename);

/// Get current working directory
std::string getCurrentDir();


#endif


