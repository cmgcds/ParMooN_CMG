/** ************************************************************************ 
* @brief     source file for GPU Routines of Particle Deposition within the problem
             
* @brief  Description :

            Contains all the functions that are needed to calculate the particle position within the cells.

            Functions and formulas are same as implemented in the sequential Implementation of ParMooN.
            They are just simplified ( in Datastructure wise ) for implementing in GPU

* @author  Thivin Anandh D
* @date      20-Jun-2023
* @History   Jun 20  - Creation of Cuda Routines
 ************************************************************************  */
#include <cuda.h>
#include "helper_cuda.h"

