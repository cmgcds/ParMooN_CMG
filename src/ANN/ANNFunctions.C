#include <ANNFunctions.h>


/// Converts string to int
int stringToInt(std::string input)
{
    std::stringstream convert(input);
    int output;
    convert >> output;
    return output;
};

/// Converts string to float
float stringToFloat(std::string input)
{
    std::stringstream convert(input);
    float output;
    convert >> output;
    return output;
};

/// Converts string to double
double stringToDouble(std::string input)
{
    std::stringstream convert(input);
    double output;
    convert >> output;
    return output;
};


/// Checks whether the given file exists
bool fileExists(const std::string& filename){
    std::ifstream ifile(filename.c_str());
    return (bool)ifile;
};

/// Get current working directory
std::string getCurrentDir(){
   char buff[FILENAME_MAX]; //create string buffer to hold path
   getcwd( buff, FILENAME_MAX );
   std::string current_working_dir(buff);
   return current_working_dir;
};

