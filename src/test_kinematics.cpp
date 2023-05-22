#define _USE_MATH_DEFINES 
#include "forward_kinematics.h"
#include"data.h"
#include "utils.h"

int main()
{
	// double para[] = {1.33997, 0.98804, -3.811, -1.77003, -2.93483, -0.321918, -1.28835};
	// double para[] = {1.82072 , 0.587827 , 1.0189 , 4.33098 , -3.20266 , 1.06863 , 5.27103};
	double para[] = {1.23848 , 0.400183 , -0.838996 , 2.94077 , 1.44408 , -2.0944 , -1.96748};
	
    double* eta_end = new double;
    double* xi_end = new double[3];
    double* Pe = new double[3];
    double* eta_b = new double;
    double* xi_b = new double[3];
    double* quaternion_end_desired = new double[4];
    // double* quaternion_base_initial = new double[4];
    double* eta_base_initial = new double;
    double* xi_base_initial = new double[3];
    double* eta_end_desired = new double;
    double* xi_end_desired = new double[3];
    double* delta_eta_base;
    double* delta_eta_end;
    double* delta_xi_base;
    double* delta_xi_end;
    double* delta_Pe_end;
    double* rpy_end = new double[3];
    double* rpy_base = new double[3];
    double* delta_euler_end = new double[3];
    double* delta_euler_base = new double[3];
    
    delta_eta_base = new double;
    delta_eta_end = new double;
    delta_xi_base = new double[3];
    delta_xi_end = new double[3];
    delta_Pe_end = new double[3];
    
    double K = 50;
    double cost_func;

    double delta_xi_end_mod_temp = 0;   
    double delta_Pe_end_mod_temp = 0;
    double delta_xi_base_mod_temp = 0;
    
    forward_kin(para, eta_end, xi_end, Pe, eta_b, xi_b);

    quaternion2euler(*eta_end, xi_end[0], xi_end[1], xi_end[2], rpy_end);
    
    cout << "quaternion_end:" << endl;
    cout << *eta_end << "  " << xi_end[0] << "  " << xi_end[1] << "  " << xi_end[2] << endl;
    
    cout << "rpy_end:" << endl;
    cout << rpy_end[0] * 180 / M_PI << "  " << rpy_end[1] * 180 / M_PI << "  " << rpy_end[2] * 180 / M_PI << endl;
    cout << "Pe:" << endl;
    cout << Pe[0] << "  " << Pe[1] << "  " << Pe[2] << endl;
    
    

    quaternion2euler(*eta_b, xi_b[0], xi_b[1], xi_b[2], rpy_base);
    
    zyx2quaternion( RPY_END_DESIRED[2],  RPY_END_DESIRED[1],  RPY_END_DESIRED[0],  quaternion_end_desired);
    
    *eta_end_desired = quaternion_end_desired[0];

    for(int i=0;i<3;i++)
    {
        xi_end_desired[i] = quaternion_end_desired[i+1];
    }
    
    cout << "quaternion_end_desired:" << endl;
    cout << *eta_end_desired << "  " << xi_end_desired[0] << "  " << xi_end_desired[1] << "  " << xi_end_desired[2] << endl;
    
    
    delta_var( Pe_DESIRED, eta_end_desired, xi_end_desired, Pe, *eta_end, xi_end, *eta_b, xi_b, delta_eta_end, delta_xi_end, delta_Pe_end, delta_eta_base, delta_xi_base);
    
    
    for(int i = 0; i < 3; i++)
    {
        // delta_xi_end_mod_temp += delta_xi_end[i] * delta_xi_end[i];
        delta_Pe_end_mod_temp += delta_Pe_end[i] * delta_Pe_end[i];
        // delta_xi_base_mod_temp += delta_xi_base[i] * delta_xi_base[i];
    }

    // delta_xi_end_mod_temp = sqrt(delta_xi_end_mod_temp);
    delta_Pe_end_mod_temp = sqrt(delta_Pe_end_mod_temp);
    // delta_xi_base_mod_temp = sqrt(delta_xi_base_mod_temp);
    

    double delta_rpy_end_norm = 0;
    for(int i=0; i<3; i++)
    {
        delta_euler_end[i] = fabs(rpy_end[i] - RPY_END_DESIRED[i]);
        delta_rpy_end_norm += delta_euler_end[i] * delta_euler_end[i];
    }
    delta_rpy_end_norm = sqrt(delta_rpy_end_norm);

    //double delta_euler_base_norm = 0;
    //for(int i = 0; i < 3; i++)
    //{
    //    delta_euler_base[i] = fabs(rpy_base[i] - RPY_BASE_INITIAL[i]);
    //    delta_euler_base_norm += delta_euler_base[i] * delta_euler_base[i];
    //}
    //delta_euler_base_norm = sqrt(delta_euler_base_norm);

    // ****************************************  RPY error of end-effector,  Pe error  ********************************************************************************
    cost_func = K * delta_rpy_end_norm + K * delta_Pe_end_mod_temp;

    // cost_func = Ka * delta_euler_end_norm + Kp * delta_Pe_end_mod_temp + Kb * delta_euler_base_norm;
    

    

    delete[] delta_eta_base;
    delete[] delta_eta_end;
    delete[] delta_xi_base;
    delete[] delta_xi_end;
    delete[] delta_Pe_end;
    delete[] xi_end_desired;
    delete[] eta_end_desired;
    delete[] quaternion_end_desired;
    delete[] delta_euler_base;
    delete[] eta_base_initial;
    delete[] xi_base_initial;
    // delete[] quaternion_base_initial;
    delete[] delta_euler_end;
    delete[] xi_b;
    delete[] eta_b;
    delete[] Pe;
    delete[] xi_end;
    delete[] eta_end;
    delete[] rpy_end;
    delete[] rpy_base;
    
    
	
}
