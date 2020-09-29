/** rearrange bytes in int array */
void SwapIntArray(int *intarray, int length)
{
  char *array, c;
  int i;

  for(i=0;i<length;i++)
  {
    array = (char *)(intarray+i);
    c = array[0];
    array[0] = array[3];
    array[3] = c;

    c = array[1];
    array[1] = array[2];
    array[2] = c;
  }
} // SwapIntArray

/** rearrange bytes in float array */
void SwapFloatArray(float *floatarray, int length)
{
  char *array, c;
  int i;

  for(i=0;i<length;i++)
  {
    array = (char *)(floatarray+i);
    c = array[0];
    array[0] = array[3];
    array[3] = c;

    c = array[1];
    array[1] = array[2];
    array[2] = c;
  }
} // SwapFloatArray

/** rearrange bytes in double array */
void SwapDoubleArray(double *doublearray, int length)
{
  char *array, c;
  int i;

  for(i=0;i<length;i++)
  {
    array = (char *)(doublearray+i);
    c = array[0];
    array[0] = array[7];
    array[7] = c;

    c = array[1];
    array[1] = array[6];
    array[6] = c;

    c = array[2];
    array[2] = array[5];
    array[5] = c;

    c = array[3];
    array[3] = array[4];
    array[4] = c;
  }
} // SwapDoubleArray
