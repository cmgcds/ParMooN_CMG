#include "Controller.h"
#include "Database.h"
#include "stdlib.h"

Controller::Controller()
    : reject(false), olderr0(1.0e-3),olderr1(1.0e-3),
      oldtau0(TDatabase::TimeDB->TIMESTEPLENGTH),
      oldtau1(TDatabase::TimeDB->TIMESTEPLENGTH)
{
  K_P=TDatabase::TimeDB->TIMESTEPLENGTH_PARA_KK_P/(TDatabase::TimeDB->T0+1);
  K_I=TDatabase::TimeDB->TIMESTEPLENGTH_PARA_KK_I/(TDatabase::TimeDB->T0+1);
  K_D=TDatabase::TimeDB->TIMESTEPLENGTH_PARA_KK_D/(TDatabase::TimeDB->T0+1);

  control_safty=TDatabase::TimeDB->CONTROL_SAFTY;
  control_maxscale=TDatabase::TimeDB->CONTROL_MAXSCALE;
  control_minscale=TDatabase::TimeDB->CONTROL_MINSCALE;
}
Controller::~Controller()
{

}

void Controller::StepLengthControl(int &m, bool &acc, double err2, double &err1, double &err0, 
                                   double &tau, double &tauold, int &step_rej)
{
  double scale=0.;
  double taunew;
  int controller=TDatabase::TimeDB->CONTROL;
  double tol=TDatabase::TimeDB->TIMESTEPLENGTH_TOL;
  acc=false;
  
  switch(controller){
    case 0: // I Controller
      scale=tau*control_safty*pow(tol/err2,K_I);
      break;
    case 1: // PI controller
      scale=tau*control_safty*pow(tol/err2,K_I)*pow(err1/err2,K_P);
      break;
    case 2: // PID controller
      scale=tau*control_safty*pow(err1/err2,K_P)*pow(tol/err2,K_I)*pow((err1*err1)/(err2*err0),K_D);
      break;
    case 3: // PC controller
      scale=(tau*tau)/tauold*control_safty*pow(tol/err2,K_I)
              *pow(err1/err2,K_P);
      break;
  }
  
  taunew=std::min(tau*control_maxscale,
                  std::max(tau*control_minscale,scale));
  
  TDatabase::TimeDB->TIMESTEPLENGTH=taunew;
  TDatabase::TimeDB->CURRENTTIMESTEPLENGTH=taunew;


  if(err2<=tol){
      acc=true;
      err0=err1;
      err1=err2;
      tauold=taunew;
      tau=TDatabase::TimeDB->TIMESTEPLENGTH;
  }
  else{
      m--;
      acc=false;
      tau=TDatabase::TimeDB->TIMESTEPLENGTH;
    step_rej++;
  }
}

void Controller::StepLengthControlTest(int &m, bool &acc, double err_at_tn_p1, 
                                   double &err_at_tn, double &err_at_tn_m1, 
                                   double &current_time_step, double &old_time_step,
                                   int &step_rej)
{
  double tau_star=0.;
  double new_time_step;
  int controller=TDatabase::TimeDB->CONTROL;
  double tol=TDatabase::TimeDB->TIMESTEPLENGTH_TOL;
  acc=false;
  OutPut("KKP: " << K_P << "  KK_I:  " << K_I << "  :K_D:   " << K_D<<endl);
  switch(controller){
    case 0: // I Controller
      tau_star= old_time_step*control_safty
                * pow(tol/err_at_tn_p1,K_I);
      break;
    case 1: // PI controller
      tau_star=  old_time_step*control_safty
               * pow( (tol/err_at_tn_p1) , K_I)
               * pow( (err_at_tn/err_at_tn_p1),K_P);
      break;
    case 2: // PID controller
      tau_star= old_time_step * control_safty
                * pow( (err_at_tn/err_at_tn_p1), K_P)
                * pow( (tol/err_at_tn_p1), K_I)
                * pow( (err_at_tn*err_at_tn)/(err_at_tn_p1*err_at_tn_m1),K_D);
      break;
    case 3: // PC controller
      tau_star=(old_time_step*old_time_step)/old_time_step
               * control_safty*pow(tol/err_at_tn_p1,K_I)
               * pow( (err_at_tn/err_at_tn_p1),K_P);
      break;
  }
  
  new_time_step=std::min(old_time_step*control_maxscale,
                  std::max(old_time_step*control_minscale,tau_star));
  
  TDatabase::TimeDB->TIMESTEPLENGTH=new_time_step;
  TDatabase::TimeDB->CURRENTTIMESTEPLENGTH=new_time_step;


  if(err_at_tn_p1<=tol){
      acc=true;
      err_at_tn_m1=err_at_tn;
      err_at_tn=err_at_tn_p1;
      old_time_step=new_time_step;
      current_time_step=TDatabase::TimeDB->TIMESTEPLENGTH;
  }
  else{
      m--;
      acc=false;
      current_time_step=TDatabase::TimeDB->TIMESTEPLENGTH;
    step_rej++;
  }
}