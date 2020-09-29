#include <Database.h>
#include <MooNMD_Io.h>

std::ofstream OutFile;

void OpenFiles()
{
  OutFile.open(TDatabase::ParamDB->OUTFILE);
  OutFile.setf(std::ios::scientific);
}

void CloseFiles()
{
  OutFile.close();
}


