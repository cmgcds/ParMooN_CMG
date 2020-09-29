// rotating rigid body   
// for rotating system
 

void ExampleFile()
{
  OutPut("Example: RotatingBody.h" << endl) ;
}

 
// ========================================================================
// boundary conditions
// ========================================================================
void GridBoundCondition(int i, double t, BoundCond &cond)
{
  cond = NEUMANN;
}

 

void ModifyCoords(double Y, double &y)
{
 double disp, t=TDatabase::TimeDB->CURRENTTIME;
  
  
  disp = 0.5*sin(Pi*t);
  y  = Y + disp;
    
 
  //    cout << "x = " << x << "  y = " <<  y << endl;
}