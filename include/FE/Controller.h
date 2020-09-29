#ifndef CONTROLLER_H
#define CONTROLLER_H

/** Step size controller **/

class Controller
{
  protected:
    /** */
    double olderr0;
    double olderr1;

    /** */
    double oldtau0;
    double oldtau1;

    /** */
    bool reject;

    /** */
    double K_P;

    /** */
    double K_I;

    /** */
    double K_D;

    /** */
    double control_safty;

    /** */
    double control_maxscale;

    /** */
    double control_minscale;
public:
    Controller();
    ~Controller();

    /** @param:
     *  input: err = error norm ..
     *         tau = timesteplength
     *  on return: true if err<= 1.
     *             false otherwise
     *          tau will be update
    **/
    bool success(int &m, bool &acc, double err, double &tau, int &step_rej);

    /** @param:
     *  input: err = error norm ..
     *         tau = timesteplength
     *  on return: true if err<= 1.
     *             false otherwise
     *          tau will be update
    **/
    void StepLengthControl(int &m, bool &acc, double err, double tau, int &step_rej);
    
    void StepLengthControl(int &m, bool &acc, double err, double &errold, double &tau, double &tauold, int &step_rej);
    
    void StepLengthControl(int &m, bool &acc, double err2, double &err1, double &err0, 
                                   double &tau, double &tauold, int &step_rej);
    
    
    void PID_Controller(double &errtn, double &errtn_p1, double &errtn_m1, 
                        int &rej_step, int &m, double &tau, double &oldtau, bool accepted);
    
    void StepLengthControlTest(int &m, bool &acc, double err_at_tn_p1, 
                                   double &err_at_tn, double &err_at_tn_m1, 
                                   double &current_time_step, double &old_time_step,
                                   int &step_rej);
};

#endif // CONTROLLER_H
