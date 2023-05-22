
#define _USE_MATH_DEFINES 
#define fscanf_s fscanf

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include "utils.h"
#include "openGJK.h"


using namespace utils;
using namespace openGJK;

namespace kinetics{

    int N = 7;
    double Ez[] = {0, 0, 1};
    double eye[] = {1, 0, 0, 0, 1, 0, 0, 0, 1}; 

    double rpy_joints[] = {
          0,            0, M_PI/2, 
          -M_PI/2, M_PI/2,      0,
          M_PI/2,       0,   M_PI,
          -M_PI/2,      0,      0,
          M_PI/2,       0,   M_PI,
          M_PI/2,       0,      0,
          -M_PI/2,      0,      0};

    double joint_angle_velocity_min_limit = -0.0873;
    double joint_angle_velocity_max_limit = 0.0873;  // angle velocity 5 deg/s
    double joint_angle_acceleration_min_limit = -0.00873;
    double joint_angle_acceleration_max_limit = 0.00873;

    double joint_angle_min_limit_degree[] = {-360, 0, -360, -360, -360, -120, -360};

    double joint_angle_min_limit_rad[] = {joint_angle_min_limit_degree[0] * M_PI / 180,
                                          joint_angle_min_limit_degree[1] * M_PI / 180,
                                          joint_angle_min_limit_degree[2] * M_PI / 180,
                                          joint_angle_min_limit_degree[3] * M_PI / 180,
                                          joint_angle_min_limit_degree[4] * M_PI / 180,
                                          joint_angle_min_limit_degree[5] * M_PI / 180,
                                          joint_angle_min_limit_degree[6] * M_PI / 180};


    double joint_angle_max_limit_degree[] = {360, 180, 360, 360, 360, 120, 360};

    double joint_angle_max_limit_rad[] = {joint_angle_max_limit_degree[0] * M_PI / 180,
                                          joint_angle_max_limit_degree[1] * M_PI / 180,
                                          joint_angle_max_limit_degree[2] * M_PI / 180,
                                          joint_angle_max_limit_degree[3] * M_PI / 180,
                                          joint_angle_max_limit_degree[4] * M_PI / 180,
                                          joint_angle_max_limit_degree[5] * M_PI / 180,
                                          joint_angle_max_limit_degree[6] * M_PI / 180};

    double m_b = 1474.1;

    double b_b[] = {0, 0, 1.0864};

    double m[] = {39.618, 13.427, 22.757, 13.427, 42.614, 13.427, 8.269};

    double a[] = {
    0, 0.013631, 0.13304,
    0, 0.018673, 0.08,
    0, 0.011017, 0.8422,
    0, 0.018673, 0.08,
    0, 0.113430, 1.1109,
    0, 0.018673, 0.08,
    0, 0, 0.10565};

    double b[] = {
    0, 0.136369, 0.056960,
    0, 0.1-0.018673, 0,
    0, 0.1-0.011017, 1.08-0.8422,
    0, 0.1-0.018673, 0,
    0, 0.08-0.11343, 1.26-1.1109,
    0, 0.1-0.018673, 0,
    0, 0, 0.16435};

    double I_b_body[] = {
    17388.34, 0, 0,
    0, 1340.43, 0,
    0, 0, 17981.26};

    double I_links_body[] = {
    0.30554, 0, 0,
    0, 0.25403, 0.030757,
    0, 0.030757, 0.15563,

    0.043493, 0, 0,
    0, 0.031896, 0,
    0, 0, 0.029035,

    2.69, 0, 0,
    0, 2.68, 0,
    0, 0, 0.06,

    0.043493, 0, 0,
    0, 0.031896, 0,
    0, 0, 0.029035,

    1.75, 0, 0,
    0, 1.47, 0.29,
    0, 0.29, 0.33,

    0.043493, 0, 0,
    0, 0.031896, 0,
    0, 0, 0.029035,

    0.04, 0, 0,
    0, 0.04, 0,
    0, 0, 0.01};


    double q_INITIAL[] = {
    0, 
    0, 
    0,
    0, 
    0, 
    0, 
    0};

    double RPY_END_INITIAL[] = {0, M_PI / 2, M_PI / 2};    

    // double delta_tau = 0.001;
    // double delta_tau = 0.01;
    double delta_tau = 0.001;

    double RPY_BASE_INITIAL[] = {0, 0, 0};

    int nvrtx = 8;

    double center2vertices[15][8][3] = 
    {
          {
                {-1.0000,  1.0000 ,  0.9864},                            
                {1.0000 ,  1.0000 ,  0.9864},
                {1.0000 ,  -1.0000,  0.9864},
                {-1.0000,  -1.0000,  0.9864},
                {-1.0000,  1.0000 , -1.0136},
                {1.0000 ,  1.0000 , -1.0136},
                {1.0000 ,  -1.0000, -1.0136},
                {-1.0000,  -1.0000, -1.0136}

          },

          {
                {1.5000, 1.0000,   0.0364},
                {1.5000 , 1.0000,  -0.0636},
                {-1.5000, 1.0000,   0.0364},
                {-1.5000, 1.0000,   -0.0636},
                {1.5000 , 7.2000,   0.0364},
                {1.5000 , 7.2000,   -0.0636},
                {-1.5000, 7.2000,   0.0364},
                {-1.5000, 7.2000,   -0.0636}

          },

          {
                {1.5000 , -7.2000 , 0.0364},
                {1.5000 , -7.2000 , -0.0636},
                {-1.5000,  -7.2000,  0.0364},
                {-1.5000,  -7.2000,  -0.0636},
                {1.5000 , -1.0000 , 0.0364},
                {1.5000 , -1.0000 , -0.0636},
                {-1.5000,  -1.0000,  0.0364},
                {-1.5000,  -1.0000,  -0.0636}

          },

          {
                {0.075000 , -0.088631 ,  -0.133040},
                {0.075000 , 0.061369  ,  -0.133040},
                {-0.075000,  -0.088631,  -0.133040},
                {-0.075000,  0.061369 ,  -0.133040},
                {0.075000 , -0.088631 ,   0.116960},
                {0.075000 , 0.061369  ,   0.116960},
                {-0.075000,  -0.088631,   0.116960},
                {-0.075000,  0.061369 ,   0.116960}

          },

          {
                {0.050000 ,  0.061369 , 0.006960},
                {0.050000 ,  0.061369 , 0.106960},
                {-0.050000,  0.061369 , 0.106960},
                {-0.050000,  0.061369 , 0.006960},
                {-0.050000,  0.136369 , 0.006960},
                {-0.050000,  0.136369 , 0.106960},
                {0.050000 ,  0.136369 , 0.106960},
                {0.050000 ,  0.136369 , 0.006960}

          },

          {
                {0.050000 , -0.068673 ,  -0.080000},
                {0.050000 , 0.031327  ,  -0.080000},
                {-0.050000,  0.031327 ,  -0.080000},
                {-0.050000,  -0.068673,  -0.080000},
                {0.050000 , -0.068673 ,   0.080000},
                {0.050000 , 0.031327  ,   0.080000},
                {-0.050000,  0.031327 ,   0.080000},
                {-0.050000,  -0.068673,   0.080000}

          },

          {
                {0.050000 , 0.031327 ,  -0.050000},
                {0.050000 , 0.081327 , -0.050000},
                {-0.050000,  0.081327,  -0.050000},
                {-0.050000,  0.031327,  -0.050000},
                {0.050000 , 0.031327 , 0.050000},
                {0.050000 , 0.081327 , 0.050000},
                {-0.050000,  0.081327,  0.050000},
                {-0.050000,  0.031327,  0.050000}

          },

          {
                {0.030000 , -0.143430 ,  -1.110900},
                {0.030000 , -0.083430 , -1.110900},
                {-0.030000,  -0.143430,  -1.110900},
                {-0.030000,  -0.083430,  -1.110900},
                {0.030000 , -0.143430 , -0.110900},
                {0.030000 , -0.083430 , -0.110900},
                {-0.030000,  -0.143430,  -0.110900},
                {-0.030000,  -0.083430,  -0.110900}

          },

          {
                {0.050000 , -0.163430 ,  -0.110900},
                {0.050000 , -0.063430 , -0.110900},
                {-0.050000,  -0.163430,  -0.110900},
                {-0.050000,  -0.063430,  -0.110900},
                {0.050000 , -0.163430 , 0.049100},
                {0.050000 , -0.063430 , 0.049100},
                {-0.050000,  -0.163430,  0.049100},
                {-0.050000,  -0.063430,  0.049100}

          },

          {
                {0.050000 , -0.063430 ,  -0.08090},
                {-0.050000,  -0.063430,  -0.08090},
                {0.050000 , -0.063430 ,  0.019100},
                {-0.050000,  -0.063430,  0.019100},
                {0.050000 , 0.146570  , -0.080900},
                {-0.050000,  0.146570 , -0.080900},
                {0.050000 , 0.146570  ,  0.019100},
                {-0.050000,  0.146570 ,  0.019100}

          },

          {
                {0.050000 , 0.016570 ,  0.019100},
                {0.050000 , 0.116570 , 0.019100},
                {-0.050000,  0.016570,  0.019100},
                {-0.050000,  0.116570,  0.019100},
                {0.050000 , 0.016570 , 0.229100},
                {0.050000 , 0.116570 , 0.229100},
                {-0.050000,  0.016570,  0.229100},
                {-0.050000,  0.116570,  0.229100}

          },

          {
                {0.050000 , 0.016570  ,  0.099100},
                {-0.050000,  0.016570 ,  0.099100},
                {-0.050000,  0.016570 ,  0.199100},
                {0.050000 , 0.016570  ,  0.199100},
                {0.050000 , -0.033430 ,  0.099100},
                {-0.050000,  -0.033430,  0.099100},
                {-0.050000,  -0.033430,  0.199100},
                {0.050000 , -0.033430 , 0.199100}

          },

          {
                {0.050000 , -0.068673 ,  -0.080000},
                {-0.050000,  -0.068673,  -0.080000},
                {-0.050000,  0.031327 ,  -0.080000},
                {0.050000 , 0.031327  , -0.080000},
                {0.050000 , -0.068673 ,  0.080000},
                {-0.050000,  -0.068673,  0.080000},
                {-0.050000,  0.031327 ,  0.080000},
                {0.050000 , 0.031327  ,  0.080000}

          },

          {
                {0.050000 , 0.031327 ,  -0.050000},
                {-0.050000,  0.031327,  -0.050000},
                {-0.050000,  0.031327,  0.050000},
                {0.050000 , 0.031327 , 0.050000},
                {0.050000 , 0.081327 , -0.050000},
                {-0.050000,  0.081327,  -0.050000},
                {-0.050000,  0.081327,  0.050000},
                {0.050000 , 0.081327 , 0.050000}

          },

          {
                {-0.040000,  -0.040000, -0.105650},
                {-0.040000,  0.040000 , -0.105650},
                {0.040000 , 0.040000  , -0.105650},
                {0.040000 , -0.040000 , -0.105650},
                {-0.040000,  -0.040000,  0.164350},
                {-0.040000,  0.040000 ,  0.164350},
                {0.040000 , 0.040000  ,  0.164350},
                {0.040000 , -0.040000 ,  0.164350}

          }
    };

    void A_b(double alpha_base, double beta_base, double gamma_base, double* A_base)
    {
        rpy2dc(alpha_base, beta_base, gamma_base, A_base);

    }


    int calc_binomial(int n, int k)
    {
        int dp[20][20];

    	if (k == 0 || k == n)
            return 1;

        else
    	    return dp[n][k] = calc_binomial(n - 1, k - 1) + calc_binomial(n - 1, k);

    }


    void calc_p(double* r, double* A_links_transform, double* p)
    {   
        double* a_tempt = new double[3];
        double* r_tempt = new double[3];
        double* A_links_transform_tempt_i = new double[9];
        double* A_links_transform_tempt_i_multi_a_tempt = new double[3];

        SetZeroMatrix_( 1,  3,  a_tempt);
        SetZeroMatrix_( 1,  3,  r_tempt);
        SetZeroMatrix_( 1,  9,  A_links_transform_tempt_i);
    
    	for(int i = 0; i < N; i++)
    	{		
    	    SetZeroMatrix_( 1, 3, A_links_transform_tempt_i_multi_a_tempt);
    
    		MatrixExtract_( 1, 3*N, 1, 1, i*3+1, i*3+3, a, a_tempt);
    
    		ScaleMatrix_(1, 3, (-1),  a_tempt, a_tempt);
    
    		MatrixExtract_( 1, 3*N, 1, 1, i*3+1, i*3+3, r, r_tempt);
    
    		MatrixExtract_( 1, 3*3*N, 1, 1, i*9+1, i*9+9, A_links_transform, A_links_transform_tempt_i);
    
    		/*
    		for(int j = 0; j < 3; j++)
    		{
    			a_tempt[j] = (-1) * a[i * 3 + j];
    			r_tempt[j] = r[i * 3 + j];
    
    			for(int k = 0; k < 3; k++)
    				A_links_transform_tempt_i[j * 3 + k] = A_links_transform[i * 9 + j * 3 + k];
    
    		}
    		*/

    		MatrixMulti_(3, 3, 1, A_links_transform_tempt_i, a_tempt, A_links_transform_tempt_i_multi_a_tempt);
    
    		for(int k = 0; k < 3; k++)
    			p[i * 3 + k] = r_tempt[k] + A_links_transform_tempt_i_multi_a_tempt[k];
    
    	}

        delete[] a_tempt;
        delete[] r_tempt;
        delete[] A_links_transform_tempt_i;
        delete[] A_links_transform_tempt_i_multi_a_tempt;
    
    }


    void calc_para_range(double* result_min, double* result_max)
    {

        double M1_temp;

        double M2_temp;

        double M1;

        double M2;

        double* result_1 = new double[N];

        double* result_1_temp = new double[N];

        double* result_2 = new double[N];

        double* result_2_temp = new double[N];

        M1_temp = pow((1-0.0001), 5) + 5 * pow((1-0.0001), 4) * 0.0001 + 10 * pow((1-0.0001), 3) * pow(0.0001, 2);

        M2_temp = 10 * pow((1-0.0001), 2) * pow(0.0001, 3) + 5 * (1-0.0001) * pow(0.0001, 4) + pow(0.0001, 5);

        for(int i = 0; i < N; i++)
        {
            result_1[i] = (joint_angle_min_limit_rad[i] - q_INITIAL[i] * M1_temp) / M2_temp;

            result_2[i] = (joint_angle_max_limit_rad[i] - q_INITIAL[i] * M1_temp) / M2_temp;

        }

        for(double tau = 0.001; tau < 1.0; tau += 0.0001)
        {
            M1 = pow((1 - tau), 5) + 5 * pow((1 - tau), 4) * tau + 10 * pow((1 - tau), 3) * pow(tau, 2);

            M2 = 10 * pow((1 - tau), 2) * pow(tau, 3) + 5 * (1 - tau) * pow(tau, 4) + pow(tau, 5);

            for(int i = 0; i < N; i++)
            {
                result_1_temp[i] = (joint_angle_min_limit_rad[i] - q_INITIAL[i] * M1) / M2;

                result_2_temp[i] = (joint_angle_max_limit_rad[i] - q_INITIAL[i] * M1) / M2;

                if(result_1_temp[i] >= result_1[i])
                {
                    result_1[i] = result_1_temp[i];
                }

                if(result_2_temp[i] <= result_2[i])
                {
                    result_2[i] = result_2_temp[i];
                }

            } 

        }   

        for(int i = 0; i < N; i++)
        {
            result_min[i] = result_1[i];

            result_max[i] = result_2[i];

        }

        delete[] result_1;
        delete[] result_1_temp;
        delete[] result_2;
        delete[] result_2_temp;

    }


    double calc_q_ddot(double tau)
    {
        return 5 * (calc_binomial(4, 2) * (-2) * (1-tau) * pow(tau, 2) + calc_binomial(4, 2) * pow((1-tau), 2) * 2 * tau);
    }


    double calc_q_dot(double tau)
    {
        // 根据 Beseir 曲线 计算得到
        return (5 * calc_binomial(4, 2) * pow((1 - tau), 2) * pow(tau, 2));
    }


    void calc_r(double* r_b, double* A_b, double* A_links_transform, double* r)
    {   
        double* a_tempt = new double[3];
        double* b_tempt = new double[3];
        double* A_links_transform_tempt = new double[9];
        double* A_links_transform_tempt_pre = new double[9];
        double* A_links_transform_multi_a_tempt = new double[3];
    	double* A_links_transform_tempt_pre_multi_b_tempt = new double[3];
        double* A_b_multi_b_b = new double[3];

        SetZeroMatrix_( 1, 3, a_tempt);
        SetZeroMatrix_( 1, 3, b_tempt);

    	for(int i = 0; i < N; i++)
    	{
            if(i == 0)
            {
                for(int j = 0; j < 3; j++)
                    a_tempt[j] = a[i * 3 + j];

                for(int j = 0; j< 3; j++)
    				for(int k = 0; k < 3; k++)
                        A_links_transform_tempt[j * 3 + k] = A_links_transform[i * 9 + j * 3 + k];

                MatrixMulti_(3, 3, 1, A_links_transform_tempt, a_tempt, A_links_transform_multi_a_tempt);

                MatrixMulti_(3, 3, 1, A_b, b_b, A_b_multi_b_b);

                for(int j = 0; j < 3; j++)
                    r[i * 3 + j] = r_b[j] + A_b_multi_b_b[j] + A_links_transform_multi_a_tempt[j];

            }
            else
            {
                for(int j = 0; j < 3; j++)
                {
                    a_tempt[j] = a[i * 3 + j];
                    b_tempt[j] = b[(i-1) * 3 + j];
                }

                for(int j = 0; j< 3; j++)
    				for(int k = 0; k < 3; k++)
    				{
    					A_links_transform_tempt[j * 3 + k] = A_links_transform[i * 9 + j * 3 + k];
    					A_links_transform_tempt_pre[j * 3 + k] = A_links_transform[(i-1) * 9 + j * 3 + k];
    				}

                MatrixMulti_(3, 3, 1, A_links_transform_tempt, a_tempt, A_links_transform_multi_a_tempt);

                MatrixMulti_(3, 3, 1, A_links_transform_tempt_pre, b_tempt, A_links_transform_tempt_pre_multi_b_tempt);

                for(int j = 0; j < 3; j++)
    			    r[i * 3 + j] = r[(i-1) * 3 + j] + A_links_transform_tempt_pre_multi_b_tempt[j] + A_links_transform_multi_a_tempt[j];

            }		    
    
    	}

        delete[] a_tempt;
        delete[] b_tempt;
        delete[] A_links_transform_tempt;
    	delete[] A_links_transform_tempt_pre;
    	delete[] A_links_transform_multi_a_tempt;
    	delete[] A_links_transform_tempt_pre_multi_b_tempt;
    	delete[] A_b_multi_b_b;	
    
    }


    void delta_var(double* Pe_desired, double* eta_end_desired, double* xi_end_desired, 
                    double* Pe, double eta_end, double* xi_end, double eta_b, double* xi_b,
                    double* delta_eta_end, double* delta_xi_end, double* delta_Pe_end, 
                    double* delta_eta_base, double* delta_xi_base)
    {
        double* quaternion_base_initial = new double[4];
        double eta_b_initial;
        double* xi_b_initial = new double[3];
        double* delta_xi_base_tempt = new double[3];
        double* cross_xi_b = new double[3*3];
        double* delta_eta_base_tempt = new double;
        double* cross_xi_end = new double[3*3];
        double* delta_xi_end_tempt = new double[3];
        double* delta_eta_end_tempt = new double;

        zyx2quaternion(RPY_BASE_INITIAL[2], RPY_BASE_INITIAL[1], RPY_BASE_INITIAL[0], quaternion_base_initial);
        eta_b_initial = quaternion_base_initial[0];
        for (int i = 0; i < 3; i++)
        {
            xi_b_initial[i] = quaternion_base_initial[i + 1];
        }


        cross(xi_b, cross_xi_b);
        MatrixMulti_(3, 3, 1, cross_xi_b, xi_b_initial, delta_xi_base_tempt);
        for(int i = 0; i < 3; i++)
        {
            delta_xi_base_tempt[i] = (-1) * delta_xi_base_tempt[i] - eta_b_initial * xi_b[i];     
            delta_xi_base_tempt[i] = delta_xi_base_tempt[i] + eta_b * xi_b_initial[i];
        }

        for (int i = 0; i < 3; i++)
        {
            delta_xi_base[i] = delta_xi_base_tempt[i];
        }

        MatrixMulti_(1, 3, 1, xi_b, xi_b_initial, delta_eta_base_tempt);
        (*delta_eta_base_tempt) = (*delta_eta_base_tempt) + eta_b * eta_b_initial;
        *delta_eta_base = (*delta_eta_base_tempt);

        cross(xi_end, cross_xi_end);
        MatrixMulti_(3, 3, 1, cross_xi_end, xi_end_desired, delta_xi_end_tempt);
        for(int i = 0; i < 3; i++)
        {
            delta_xi_end_tempt[i] = (-1) * delta_xi_end_tempt[i] - (*eta_end_desired) * xi_end[i];        
            delta_xi_end_tempt[i] = delta_xi_end_tempt[i] + eta_end * xi_end_desired[i];
        }

        for (int i = 0; i < 3; i++)
        {
            delta_xi_end[i] = delta_xi_end_tempt[i];
        }
        MatrixMulti_(1, 3, 1, xi_end, xi_end_desired, delta_eta_end_tempt);
        (*delta_eta_end_tempt) = (*delta_eta_end_tempt) + eta_end * (*eta_end_desired);    
        (*delta_eta_end) = (*delta_eta_end_tempt);

        for (int i = 0; i < 3; i++)
        {
            delta_Pe_end[i] = Pe_desired[i] - Pe[i];
        }


        delete[] quaternion_base_initial;
        delete[] xi_b_initial;
        delete[] delta_xi_base_tempt;
        delete[] cross_xi_b;
        delete[] delta_eta_base_tempt;
        delete[] cross_xi_end;
        delete[] delta_xi_end_tempt;
        delete[] delta_eta_end_tempt;

    }


    void forward_kin(double* para, double* eta_end, double* xi_end, double* Pe, 
    double* eta_b, double* xi_b, double* p_e_initial, double* locus, 
    double* delta_xi_b_distrb_max, double* manipl, double* T_min, double* collision)
    {
        double* q = new double[N];
        double* q_dot = new double[N];
        double* q_ddot = new double[N];
        double* quaternion_base_initial = new double[4];
        double* A_b = new double[3*3];
        double* xi_b_tempt = new double[3];
        double* xi_end_tempt = new double[3];
        double* A_links_transform  = new double[3*3*N];
        double* A_b_multi_b_b = new double[3];
        double* r_b_tempt_1 = new double[3];
        double* r_b_tempt_2 = new double[3];
        double* r_b_tempt_2_temp = new double[3];
        double* a_tempt_i = new double[3];
        double* r_b_tempt_3 = new double[3];
        double* r_b = new double[3];
        double* a_tempt_k = new double[3];
        double* b_tempt_k = new double[3];
        double* l_tempt_k = new double[3];    
        double* A_links_transform_tempt_i = new double[3*3];
        double* A_links_transform_tempt_i_multi_a_tempt_i = new double[3];
        double* A_links_transform_tempt_k = new double[3*3];
        double* A_links_transform_tempt_k_multi_l_tempt_k = new double[3];
        double* r_b_tempt_3_tempt = new double[3];
        double* r = new double[3 * N];
        double* r_e = new double[3];
        double* Pe_initial = new double[3];
        double* A_links_transform_N_minus_1 = new double[3*3];
        double* r_e_tempt = new double[3];
        double* b_tempt_N_minus_1 = new double[3];
        double* p = new double[N*3];
        double* A_links_transform_tempt_i_multi_Ez = new double[3];
        double* A_links_transform_tempt_i_multi_Ez_cross = new double[9];
        double* Jm_v = new double[3*N];
        double* Jm_v_tempt = new double[N*3];
        double* Jm_v_tempt_i = new double[3];
        double* Jm_w = new double[3*N];
        double* Jm_w_tempt = new double[N*3];
        double* p_tempt = new double[3];
        double* r_e_minus_p_i = new double[3];
        double* Jm = new double[6*N];
        double* Jm_tempt = new double[N*6];
        double* J_bm_w = new double[3*N];
        double* J_bm_v = new double[3*N];
        double* J_bm = new double[6*N];	
        double* J_bE = new double[6*6];
        double* J_g = new double[6*N];
        double* J_g_transpose = new double[N*6];
        double* J_g_multi_J_g_transpose = new double[6*6];
        double* J_bE_multi_J_bm = new double[6*N];
        double* J_g_v = new double[3*N];
        double* J_g_w = new double[3*N];    
        double* quaternion_end  = new double[4];
        double* eta_b_dot_tempt_1 = new double[N];
        double* eta_b_dot_tempt_2 = new double;  
        double* xi_b_dot = new double[3];
        double* cross_xi_b = new double[3*3];  
        double* xi_b_dot_tempt_1 = new double[3*3];
        double* xi_b_dot_tempt_2 = new double[3*N];
        double* xi_b_dot_tempt_3 = new double[3]; 
        double* eta_end_dot_tempt_1 = new double[N];
        double* eta_end_dot_tempt_2 = new double;
        double* xi_end_dot = new double[3];
        double* cross_xi_end = new double[3*3]; 
        double* xi_end_dot_tempt_1 = new double[3*3];
        double* xi_end_dot_tempt_2 = new double[3*N];
        double* xi_end_dot_tempt_3 = new double[3];
        double* v_e = new double[3];
        double* v_b = new double[3];
        double* A_b_expand = new double[2*3 * 2*3];
    	double* A_b_expand_transpose = new double[2*3 * 2*3];
    	double* A_b_expand_transpose_multi_J_g = new double[6 * N];
    	double* A_b_expand_transpose_multi_J_g_transpose = new double[N * 6];
    	double* manipl_temp = new double[6 * 6];
        double* xi_b_initial = new double[3];
        double* delta_xi_base = new double[3];
        double* delta_eta_base_tempt = new double;
        double* delta_eta_base = new double;
        double* delta_xi_base_tempt = new double[3];

        double eta_end_dot = 0;
        double eta_b_dot = 0;
        double total_mass_for_links = 0;
    	double total_mass = 0;	
    	double q_dot_temp = 0;
        double q_ddot_temp = 0;
        double eta_b_initial = 0;
        double delta_xi_base_norm_max = 0;
        double delta_xi_base_norm = 0;
        double q_dot_max = 0;
        double q_ddot_max = 0;
        double q_dot_max_vmax = 0;
        double q_ddot_max_alphamax = 0;
        double T_min_temp = 0;


        double center2vertices_temp_temp[nvrtx * 3];

        double** objects_temp = new double*[15];
        for(int i = 0; i < 15; i++)
        {
            objects_temp[i] = new double[nvrtx * 3];
        }

        double** A_links_transform_i = new double*[15];
        for(int i = 0; i < 15; i++)
        {
            A_links_transform_i[i] = new double[3*3];
        }
        // set A_links_transform_i to 0
        for ( int i = 0; i < 15; i++ )
        {
            for ( int j = 0; j < 3*3; j++ )
            {
                A_links_transform_i[i][j] = 0;
            }
        }

        double*** objects = new double**[15];
        for (int i = 0; i < 15; i++) 
        {
            objects[i] = new double*[nvrtx];

            for (int j = 0; j < nvrtx; j++) 
            {
                objects[i][j] = new double[3];
            }

        }

        gkSimplex s;
        s.nvrtx = 0;

        gkPolytope* bd_objects = new gkPolytope[15];

        double collision_test_base_link567 = 0;
        double collision_test_left_link567 = 0;
        double collision_test_right_link567 = 0;
        double collision_test_link_1_link567 = 0;
        double collision_test_link_2_link567 = 0;
        double collision_test = 0;

        double* dis_base2links = new double[8];

        double* dis_left2links = new double[7];

        double* dis_right2links = new double[7];

        double* dis_link_1_2_links = new double[14];

        double* dis_link_2_2_links = new double[14];

        double distance_limit = 0.001;


        double*** center2vertices_buff = new double**[15];
        for(int i = 0; i < 15; i++)
        {
            center2vertices_buff[i] = new double*[nvrtx];

            for(int j = 0; j < nvrtx; j++)
            {
                center2vertices_buff[i][j] = new double[3];
            }
        }

        for(int i = 0; i < 15; i++)
        {
            for(int j = 0; j < nvrtx; j++)
            {
                for(int k = 0; k < 3; k++)
                {
                    center2vertices_buff[i][j][k] = center2vertices[i][j][k];
                }
            }
        }

        double** center2vertices_temp = new double*[15];
        for(int i = 0; i < 15; i++)
        {
            center2vertices_temp[i] = new double[nvrtx * 3];
        }

        for (int i = 0; i < 15; i++)
        {
            for (int j = 0; j < nvrtx; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    center2vertices_temp[i][j * 3 + k] = center2vertices_buff[i][j][k];
                }
            }
        }

        for (int i = 0; i < 15; i++)
        {
            for (int j = 0; j < nvrtx * 3; j++)
            {
                center2vertices_temp_temp[j] = center2vertices_temp[i][j];
            } 

            // center2vertices_temp_temp is transformed into 3row * 8 column
            MatrixTranspose_(nvrtx, 3, center2vertices_temp_temp, center2vertices_temp_temp);

            for (int j = 0; j < nvrtx * 3; j++)
            {
                center2vertices_temp[i][j] = center2vertices_temp_temp[j];
            } 
        }


    //****************************************************************************************************
        for(int i = 0; i < N; i++)
            q[i] = q_INITIAL[i];

        zyx2quaternion(RPY_BASE_INITIAL[2], RPY_BASE_INITIAL[1], 
                        RPY_BASE_INITIAL[0], quaternion_base_initial);

        quaternion2dc(quaternion_base_initial[0], quaternion_base_initial[1], 
                        quaternion_base_initial[2], quaternion_base_initial[3], A_b);

        *eta_b = quaternion_base_initial[0];

        eta_b_initial = quaternion_base_initial[0];

        for(int i = 0; i < 3; i++)
        {
            xi_b[i] = quaternion_base_initial[i+1];

            xi_b_initial[i] = quaternion_base_initial[i+1];

        }

    	*locus = 0;


    // ****************************************************************************************************
        links_transform(A_b, q, A_links_transform);

        for(int i = 0; i < 3; i++)
        {
            r_b_tempt_2[i] = 0;
        }

        for(int i = 0; i < 3; i++)
        {
            r_b_tempt_3[i] = 0;
        }

        // 自由漂浮，初始动量为 0 ， 以 系统质心为 惯性系原点。
        for (int i = 0; i < N; i++)
        {
            total_mass_for_links += m[i];
        }

    	total_mass = m_b + total_mass_for_links;

    	MatrixMulti_(3, 3, 1, A_b, b_b, A_b_multi_b_b);	

    	ScaleMatrix_( 3, 1, total_mass_for_links, A_b_multi_b_b, r_b_tempt_1);
    
    	for(int i = 0; i < N; i++)
    	{
    		for(int j = 0; j < 3; j++)
    		{
    			a_tempt_i[j] = a[i*3+j];

    			r_b_tempt_3_tempt[j] = 0;	

    			for(int k = 0; k < 3; k++)
                {
                    A_links_transform_tempt_i[j * 3 + k] = A_links_transform[i * 9 + j * 3 + k];
                }
    		}

    		MatrixMulti_( 3, 3, 1, A_links_transform_tempt_i, a_tempt_i, 
                            A_links_transform_tempt_i_multi_a_tempt_i);	

    		ScaleMatrix_( 1, 3, m[i], 
                            A_links_transform_tempt_i_multi_a_tempt_i, r_b_tempt_2_temp);		

    		MatrixAdd_( 1, 3, r_b_tempt_2, r_b_tempt_2_temp, r_b_tempt_2);

    		for(int k = 0; k < i; k++)
    		{
    			for(int jj = 0; jj < 3; jj++)
    			{
    				a_tempt_k[jj] = a[k * 3 + jj];

    				b_tempt_k[jj] = b[k * 3 + jj];

    				l_tempt_k[jj] = a_tempt_k[jj] + b_tempt_k[jj];	

    				for(int kk = 0; kk < 3; kk++)
                    {
                        A_links_transform_tempt_k[jj * 3 + kk] = A_links_transform[k * 9 + jj * 3 + kk];
                    }

    			}			

    			MatrixMulti_(3, 3, 1, A_links_transform_tempt_k, 
                                l_tempt_k, A_links_transform_tempt_k_multi_l_tempt_k);	

    			MatrixAdd_( 1, 3, r_b_tempt_3_tempt, 
                            A_links_transform_tempt_k_multi_l_tempt_k, r_b_tempt_3_tempt);
    
    		}
    
    		ScaleMatrix_( 1, 3, m[i], r_b_tempt_3_tempt, r_b_tempt_3_tempt);

    		MatrixAdd_( 1, 3, r_b_tempt_3, r_b_tempt_3_tempt, r_b_tempt_3);
    
    	}

        //  计算 基座 位置 
        MatrixAdd_( 1, 3, r_b_tempt_1, r_b_tempt_2, r_b_tempt_2);

    	MatrixAdd_( 1, 3, r_b_tempt_2, r_b_tempt_3, r_b_tempt_3);

    	ScaleMatrix_( 1, 3, (-1) / total_mass, r_b_tempt_3, r_b);

    	// 计算各个连杆位置	
    	calc_r(r_b, A_b, A_links_transform, r);
    
        //  计算末端位置
    	for(int i=0;i<3;i++)
    	{
    		b_tempt_N_minus_1[i] = b[(N-1)*3+i];	

    		for(int j=0;j<3;j++)
            {
                A_links_transform_N_minus_1[i*3+j] = A_links_transform[(N-1)*9+i*3+j];
            }
    
    	}

    	MatrixMulti_(3, 3, 1, A_links_transform_N_minus_1, b_tempt_N_minus_1, r_e_tempt);

    	for(int i=0;i<3;i++)
    		r_e[i] = r[(N-1)*3+i] + r_e_tempt[i];    
    
        // ***************************** 初始末端位姿 ******************************************************   	
        for(int i=0;i<3;i++)
            Pe_initial[i] = r_e[i];  

        zyx2quaternion(RPY_END_INITIAL[2], RPY_END_INITIAL[1], RPY_END_INITIAL[0], quaternion_end);  

        *eta_end = quaternion_end[0];

        for(int i=0;i<3;i++)
            xi_end[i] = quaternion_end[i+1];

        for(int i=0;i<3;i++)
        {
            Pe[i] = Pe_initial[i];

            p_e_initial[i] = Pe_initial[i];
        }

        // *****************************************************************************************      
        // 计算 关节位置
    	calc_p(r, A_links_transform, p);    

        // ******  计算 初始 Jm_v, Jm_w, Jm, J_bm_w, J_bm_v, J_bm, J_g_v, J_g_w, J_g  ***************
        for(int i = 0; i < N; i++)
    	{

    		MatrixExtract_( 1, 3*N, 1, 1, i*3+1, i*3+3, p, p_tempt);

    		ScaleMatrix_( 1, 3, (-1), p_tempt, p_tempt);	

    		MatrixAdd_( 1, 3, r_e, p_tempt, r_e_minus_p_i);	

    		MatrixExtract_( 1, 3*3*N, 1, 1, i*9+1, i*9+9, A_links_transform, A_links_transform_tempt_i);	

    		MatrixMulti_(3, 3, 1, A_links_transform_tempt_i, Ez, A_links_transform_tempt_i_multi_Ez);	

    		cross(A_links_transform_tempt_i_multi_Ez, A_links_transform_tempt_i_multi_Ez_cross);	

    		MatrixMulti_(3, 3, 1, A_links_transform_tempt_i_multi_Ez_cross, r_e_minus_p_i, Jm_v_tempt_i);		
    
    		for(int k = 0; k < 3; k++)
    		{
    			Jm_v_tempt[i * 3 + k] = Jm_v_tempt_i[k];

    			Jm_w_tempt[i * 3 + k] = A_links_transform_tempt_i_multi_Ez[k];
    		}
    
    		for(int jj = 0; jj < 6; jj++)
    		{
    			if(jj<3)
    				Jm_tempt[i * 6 + jj] = Jm_v_tempt[i * 3 + jj];
    			else
    				Jm_tempt[i * 6 + jj] = Jm_w_tempt[i * 3 + (jj - 3)];
    		}
    	}

    	MatrixTranspose_(N, 3, Jm_v_tempt, Jm_v);

    	MatrixTranspose_(N, 3, Jm_w_tempt, Jm_w);

    	MatrixTranspose_(N, 6, Jm_tempt, Jm);

        calc_J_bm( r, r_b, A_b, A_links_transform, p, J_bm_w, J_bm_v);

    	for(int j = 0; j < 6; j++)
        {
    		for(int k = 0; k < N; k++)
    			{
    				if(j < 3)
    					J_bm[j * N + k] = J_bm_v[j * N + k];
    				else
    					J_bm[j * N + k] = J_bm_w[(j-3) * N + k];
    			}

        }

    	J_Base2EE(r_e, r_b, J_bE);

    	MatrixMulti_(6, 6, N, J_bE, J_bm, J_bE_multi_J_bm);

    	MatrixAdd_( 6, N, J_bE_multi_J_bm, Jm, J_g);		
    
    	for(int i=0;i<3;i++)
        {
    		for(int j=0;j<N;j++)
    		{
    			J_g_v[i*N+j] = J_g[i*N+j];

    			J_g_w[i*N+j] = J_g[(i+3)*N+j];
    		}
        }

        // ***************  进入迭代计算, 步长 delta_tau  **************************************

        for(double tau = 0; tau < 1; tau += delta_tau)
        {
            for(int i = 0; i < N; i++)
            {
                q_dot_temp = calc_q_dot(tau);

                q_ddot_temp = calc_q_ddot(tau);

                q_dot[i] = q_dot_temp * (para[i] - q_INITIAL[i]);

                q_ddot[i] = q_ddot_temp * (para[i] - q_INITIAL[i]);

                q[i] = q[i] + q_dot[i] * delta_tau;  

                if(fabs(q_dot[i]) > q_dot_max)
                    q_dot_max = fabs(q_dot[i]);  // 

                if(fabs(q_ddot[i]) > q_ddot_max)
                    q_ddot_max = fabs(q_ddot[i]);  // 
            }

            if(q_dot_max / joint_angle_velocity_max_limit > q_dot_max_vmax)
                q_dot_max_vmax = q_dot_max / joint_angle_velocity_max_limit;

            if(q_ddot_max / joint_angle_acceleration_max_limit > q_ddot_max_alphamax)
                q_ddot_max_alphamax = q_ddot_max / joint_angle_acceleration_max_limit;

            if(sqrt(q_ddot_max_alphamax > T_min_temp))
                T_min_temp = sqrt(q_ddot_max_alphamax);

            if(q_dot_max_vmax > T_min_temp)
                T_min_temp = q_dot_max_vmax;

            // ************************** eta_b_dot = (- xi_b.T @ J_bm_w @ q_dot) / 2    **************************************************************************

            for(int i = 0; i < 3; i++)
            {
                xi_b_tempt[i] = (-1) * xi_b[i];
            }   

            MatrixMulti_( 1, 3, N, xi_b_tempt, J_bm_w, eta_b_dot_tempt_1);

            MatrixMulti_( 1, N, 1, eta_b_dot_tempt_1, q_dot, eta_b_dot_tempt_2);

            eta_b_dot = (*eta_b_dot_tempt_2) / 2;

            //********************************************************************************************************************************************************

    		// ****************************************** xi_b_dot = ((eta_b * np.eye(3) - cross(xi_b)) @ J_bm_w @ q_dot) / 2  ***************************************

            for(int i = 0; i < 9; i++)
            {
                xi_b_dot_tempt_1[i] = (*eta_b) * eye[i];
            }

            cross(xi_b, cross_xi_b);

            for(int i = 0; i < 9; i++)
            {
                xi_b_dot_tempt_1[i] = xi_b_dot_tempt_1[i] - cross_xi_b[i];
            }

            MatrixMulti_( 3, 3, N, xi_b_dot_tempt_1, J_bm_w, xi_b_dot_tempt_2);

            MatrixMulti_( 3, N, 1, xi_b_dot_tempt_2, q_dot, xi_b_dot_tempt_3);  

            for(int i = 0; i < 3; i++)
            {
                xi_b_dot[i] = xi_b_dot_tempt_3[i] / 2;
            }

    		//****************************************************************************************************************************************************        

    		// ************************************** eta_end_dot = (- xi_end.T @ J_g_w @ q_dot) / 2  **************************************************************

            for(int i = 0; i < 3; i++)
            {
                xi_end_tempt[i] = (-1)*xi_end[i];
            }

            MatrixMulti_( 1, 3, N, xi_end_tempt, J_g_w, eta_end_dot_tempt_1);

            MatrixMulti_( 1, N, 1, eta_end_dot_tempt_1, q_dot, eta_end_dot_tempt_2);

            eta_end_dot = (*eta_end_dot_tempt_2) / 2;

    		//******************************************************************************************************        

    		// **** xi_end_dot = ((eta_end * np.eye(3) - cross(xi_end)) @ J_g_w @ q_dot) / 2  **********************

            for(int i = 0; i < 9; i++)
            {
                xi_end_dot_tempt_1[i] = (*eta_end) * eye[i];
            }

            cross(xi_end, cross_xi_end);

            for(int i = 0; i < 9; i++)
            {
                xi_end_dot_tempt_1[i] = xi_end_dot_tempt_1[i] - cross_xi_end[i];
            }

            MatrixMulti_( 3, 3, N, xi_end_dot_tempt_1, J_g_w, xi_end_dot_tempt_2);

            MatrixMulti_( 3, N, 1, xi_end_dot_tempt_2, q_dot, xi_end_dot_tempt_3);

            for(int i=0;i<3;i++)
            {
                xi_end_dot[i] = xi_end_dot_tempt_3[i] / 2;
            }

    		//*******************************************************************************************

    		// ********  next eta_b, xi_b, eta_end, xi_end, Pe, r_b, A_b, A_links_transform, r, r_e, p, Jm, J_bm, J_g ***************

            (*eta_b) = (*eta_b) + eta_b_dot * delta_tau;

            for(int i = 0; i < 3; i++)
            {
                xi_b[i] = xi_b[i] + xi_b_dot[i] * delta_tau;
            }

            //***************  Normalization *********************************************************
            (*eta_b) = (*eta_b) / sqrt((*eta_b) * (*eta_b) + xi_b[0] * xi_b[0] + xi_b[1] * xi_b[1] + xi_b[2] * xi_b[2]);

            for(int i = 0; i < 3; i++)
            {
                xi_b[i] = xi_b[i] / sqrt((*eta_b) * (*eta_b) + xi_b[0] * xi_b[0] + xi_b[1] * xi_b[1] + xi_b[2] * xi_b[2]);
            }

    		// *********************  计算基座扰动 RPY and delta_xi_b 最大值  ************************************

            cross(xi_b, cross_xi_b);

            MatrixMulti_(3, 3, 1, cross_xi_b, xi_b_initial, delta_xi_base_tempt);

            for(int i = 0; i < 3; i++)
            {
                delta_xi_base_tempt[i] = (-1) * delta_xi_base_tempt[i] - eta_b_initial * xi_b[i];    

                delta_xi_base_tempt[i] = delta_xi_base_tempt[i] + (*eta_b) * xi_b_initial[i];
            }

            for(int i=0;i<3;i++)
                delta_xi_base[i] = delta_xi_base_tempt[i];

            MatrixMulti_(1, 3, 1, xi_b, xi_b_initial, delta_eta_base_tempt);

            (*delta_eta_base_tempt) = (*delta_eta_base_tempt) + (*eta_b) * eta_b_initial;

            (*delta_eta_base) = (*delta_eta_base_tempt);

            delta_xi_base_norm = sqrt(delta_xi_base[0] * delta_xi_base[0] + delta_xi_base[1] * delta_xi_base[1] + delta_xi_base[2] * delta_xi_base[2]);

            if(delta_xi_base_norm > delta_xi_base_norm_max)
                delta_xi_base_norm_max = delta_xi_base_norm;

    		// *****************************************************************************************

            (*eta_end) = (*eta_end) + eta_end_dot * delta_tau;

            for(int i=0;i<3;i++)
            {
                xi_end[i] = xi_end[i] + xi_end_dot[i] * delta_tau;
            }

            //****************************  Normalization ***************************************************
            (*eta_end) = (*eta_end) / sqrt((*eta_end) * (*eta_end) + xi_end[0] * xi_end[0] + xi_end[1] * xi_end[1] + xi_end[2] * xi_end[2]);

            for(int i = 0; i < 3; i++)
            {
                xi_end[i] = xi_end[i] / sqrt((*eta_end) * (*eta_end) + xi_end[0] * xi_end[0] + xi_end[1] * xi_end[1] + xi_end[2] * xi_end[2]);
            }

            MatrixMulti_( 3, N, 1, J_g_v, q_dot, v_e);

    		// *****************  计算末端走过的路程  *******************************************************
    		(*locus) += sqrt(v_e[0] * v_e[0] + v_e[1] * v_e[1] + v_e[2] * v_e[2]) * delta_tau;
    		//*****************************************************************************************

            for(int i=0;i<3;i++)
                Pe[i] = Pe[i] + v_e[i] * delta_tau;

            MatrixMulti_( 3, N, 1, J_bm_v, q_dot, v_b);

            for(int i=0;i<3;i++)
                r_b[i] = r_b[i] + v_b[i] * delta_tau;

            quaternion2dc((*eta_b), xi_b[0], xi_b[1], xi_b[2], A_b);     

            links_transform(A_b, q, A_links_transform);     

            calc_r(r_b, A_b, A_links_transform, r);

            for(int i=0;i<3;i++)
    	    {
    		    b_tempt_N_minus_1[i] = b[(N-1)*3+i];	

    		    for(int j=0;j<3;j++)
    			    A_links_transform_N_minus_1[i*3+j] = A_links_transform[(N-1)*9+i*3+j];

    	    }

    	    MatrixMulti_(3, 3, 1, A_links_transform_N_minus_1, b_tempt_N_minus_1, r_e_tempt);

            //  计算末端位置
    	    for(int i=0;i<3;i++)
    		    r_e[i] = r[(N-1)*3+i] + r_e_tempt[i]; 

            // 计算 关节位置
    	    calc_p(r, A_links_transform, p);

            //  Jm, Jbm, Jg
            for(int i=0;i<N;i++)
    	    {
            
            //******************************************************************************************
    		    MatrixExtract_( 1, 3*N, 1, 1, i*3+1, i*3+3, p, p_tempt);	

    		    ScaleMatrix_( 1, 3, (-1), p_tempt, p_tempt);	

    		    MatrixAdd_( 1, 3, r_e, p_tempt, r_e_minus_p_i);		

    		    MatrixExtract_( 1, 3*3*N, 1, 1, i*9+1, i*9+9, A_links_transform, A_links_transform_tempt_i);
            //*******************************************************************************************
    
    		    MatrixMulti_(3, 3, 1, A_links_transform_tempt_i, Ez, A_links_transform_tempt_i_multi_Ez);	

    		    cross(A_links_transform_tempt_i_multi_Ez, A_links_transform_tempt_i_multi_Ez_cross);	

    		    MatrixMulti_(3, 3, 1, A_links_transform_tempt_i_multi_Ez_cross, r_e_minus_p_i, Jm_v_tempt_i);
    
    		    for(int k = 0; k < 3; k++)
    		    {
    			    Jm_v_tempt[i * 3 + k] = Jm_v_tempt_i[k];

    			    Jm_w_tempt[i * 3 + k] = A_links_transform_tempt_i_multi_Ez[k];
    		    }
    
    		    for(int jj = 0; jj < 6; jj++)
    		    {
    			    if(jj<3)
    				    Jm_tempt[i * 6 + jj] = Jm_v_tempt[i * 3 + jj];
    			    else
    				    Jm_tempt[i * 6 + jj] = Jm_w_tempt[i * 3 + (jj - 3)];
    		    }
    
    	    }
    
    	    MatrixTranspose_(N, 3, Jm_v_tempt, Jm_v);

    	    MatrixTranspose_(N, 3, Jm_w_tempt, Jm_w);

    	    MatrixTranspose_(N, 6, Jm_tempt, Jm);

    	    calc_J_bm(r, r_b, A_b, A_links_transform, p, J_bm_w, J_bm_v);

    	    for(int j = 0; j < 6; j++)
            {
    		    for(int k = 0; k < N; k++)
    			    {
    				    if(j < 3)
                        {
                            J_bm[j * N + k] = J_bm_v[j * N + k];
                        }
    				    else
                        {
                            J_bm[j * N + k] = J_bm_w[(j-3) * N + k];
                        }
    
    			    }
            }

    	    J_Base2EE(r_e, r_b, J_bE);

    	    MatrixMulti_(6, 6, N, J_bE, J_bm, J_bE_multi_J_bm);

    	    for(int i = 0; i < 6; i++)
            {
                for(int j = 0; j < N; j++)
                {
                    J_g[i * N + j] = J_bE_multi_J_bm[i * N + j] + Jm[i * N + j];
                }
            }
    
    
    
    	    for(int i=0;i<3;i++)
            {
    		    for(int j=0;j<N;j++)
    		    {
    			    J_g_v[i*N+j] = J_g[i*N+j];

    			    J_g_w[i*N+j] = J_g[(i+3)*N+j];
    		    }
            } 


            // ********************* calculate distances **********************************
            for (int i = 0; i < 15; i++)
            {
                if ( i < 3 )
                {
                    MatrixMulti_(3, 3, nvrtx, A_b, center2vertices_temp[i], objects_temp[i]);

                    MatrixTranspose_(3, nvrtx, objects_temp[i], objects_temp[i]);

                }

                else if ( i > 2 && i < 5 )
                {
                    MatrixExtract_( 1, 3*3*N, 1, 1, 0*9+1, 0*9+9, A_links_transform, A_links_transform_i[i]);

                    MatrixMulti_(3, 3, nvrtx, A_links_transform_i[i], center2vertices_temp[i], objects_temp[i]);

                    MatrixTranspose_(3, nvrtx, objects_temp[i], objects_temp[i]);

                }

                else if ( i > 4 && i < 7 )
                {
                    MatrixExtract_( 1, 3*3*N, 1, 1, 1*9+1, 1*9+9, A_links_transform, A_links_transform_i[i]);

                    MatrixMulti_(3, 3, nvrtx, A_links_transform_i[i], center2vertices_temp[i], objects_temp[i]);

                    MatrixTranspose_(3, nvrtx, objects_temp[i], objects_temp[i]);

                }

                else if ( i > 6 && i < 12 )
                {
                    MatrixExtract_( 1, 3*3*N, 1, 1, 4*9+1, 4*9+9, A_links_transform, A_links_transform_i[i]);

                    MatrixMulti_(3, 3, nvrtx, A_links_transform_i[i], center2vertices_temp[i], objects_temp[i]);

                    MatrixTranspose_(3, nvrtx, objects_temp[i], objects_temp[i]);

                }

                else if ( i > 11 && i < 14 )
                {
                    MatrixExtract_( 1, 3*3*N, 1, 1, 5*9+1, 5*9+9, A_links_transform, A_links_transform_i[i]);

                    MatrixMulti_(3, 3, nvrtx, A_links_transform_i[i], center2vertices_temp[i], objects_temp[i]);

                    MatrixTranspose_(3, nvrtx, objects_temp[i], objects_temp[i]);

                }

                else
                {
                    MatrixExtract_( 1, 3*3*N, 1, 1, 6*9+1, 6*9+9, A_links_transform, A_links_transform_i[i]);

                    MatrixMulti_(3, 3, nvrtx, A_links_transform_i[i], center2vertices_temp[i], objects_temp[i]);

                    MatrixTranspose_(3, nvrtx, objects_temp[i], objects_temp[i]);

                }

            }

            for ( int i = 0; i < 15; i++)
            {
                for ( int j = 0; j < nvrtx; j++)
                {
                    for ( int k = 0; k < 3; k++)
                    {
                        objects[i][j][k] = objects_temp[i][j * 3 + k];
                    }
                }
            }

            for ( int i = 0; i < 15; i++)
            {
                for ( int j = 0; j < nvrtx; j++)
                {
                    for ( int k = 0; k < 3; k++)
                    {
                        if ( i < 3)
                        {
                            objects[i][j][k] = objects[i][j][k] + r_b[k];
                        }

                        else if ( i == 3 || i == 4)
                        {
                            objects[i][j][k] = objects[i][j][k] + r[0 * 3 + k];
                        }

                        else if ( i == 5 || i == 6)
                        {
                            objects[i][j][k] = objects[i][j][k] + r[1 * 3 + k];
                        }

                        else if ( i == 7 || i == 8 || i == 9 || i == 10 || i == 11)
                        {
                            objects[i][j][k] = objects[i][j][k] + r[4 * 3 + k];
                        }

                        else if ( i == 12 || i == 13 )
                        {
                            objects[i][j][k] = objects[i][j][k] + r[5 * 3 + k];
                        }

                        else
                        {
                            objects[i][j][k] = objects[i][j][k] + r[6 * 3 + k];
                        }

                    }
                }
            }


            for ( int i = 0; i < 15; i++)
            {
                bd_objects[i].coord = objects[i];

                bd_objects[i].numpoints = nvrtx;

            }

            for ( int i = 0; i < 15; i++)
            {

                if ( i > 6 )
                {
                    // static int j = 0;  // bug

                    dis_base2links[ i - 7 ] = compute_minimum_distance( bd_objects[0],  bd_objects[i], &s);

                    // j++;

                }

                if ( i > 7 )
                {
                    // static int j = 0;

                    dis_left2links[ i - 8 ] = compute_minimum_distance( bd_objects[1],  bd_objects[i], &s);

                    dis_right2links[ i - 8 ] = compute_minimum_distance( bd_objects[2],  bd_objects[i], &s);

                    dis_link_1_2_links[ i - 8 ] = compute_minimum_distance( bd_objects[3],  bd_objects[i], &s);

                    dis_link_2_2_links[ i - 8 ] = compute_minimum_distance( bd_objects[5],  bd_objects[i], &s);

                    dis_link_1_2_links[ i - 1 ] = compute_minimum_distance( bd_objects[4],  bd_objects[i], &s);

                    dis_link_2_2_links[ i - 1 ] = compute_minimum_distance( bd_objects[6],  bd_objects[i], &s);

                    // j++;

                }

            }



            for(int ii = 0; ii < 8; ii++)
            {
                if(dis_base2links[ii] < distance_limit)
                {
                    collision_test_base_link567 += 1; 
                }
            }

            for(int ii = 0; ii < 7; ii++)
            {
                if(dis_left2links[ii] < distance_limit)
                {
                    collision_test_left_link567 += 1; 
                }

                if(dis_right2links[ii] < distance_limit)
                {
                    collision_test_right_link567 += 1; 
                }
            }


            for(int ii = 0; ii < 14; ii++)
            {
                if(dis_link_1_2_links[ii] < distance_limit)
                {
                    collision_test_link_1_link567 += 1; 
                }
            }

            for(int ii = 0; ii < 14; ii++)
            {
                if(dis_link_2_2_links[ii] < distance_limit)
                {
                    collision_test_link_2_link567 += 1; 
                }
            }

            collision_test = collision_test + collision_test_base_link567 + collision_test_left_link567 + collision_test_right_link567 + collision_test_link_1_link567 + collision_test_link_2_link567;

        }

        //*******************************************************************************************

        (*delta_xi_b_distrb_max) = delta_xi_base_norm_max;

        (*T_min) = T_min_temp;

        (*collision) = collision_test;    //  total collision times during whole motion process

    	// *************************    计算最终时刻的可操作度   ****************************************
    
        double manipl_temp_2 = 0;

    	// MatrixDiagExpand( A_b, 3, 3, A_b_expand);
    	// MatrixTranspose_(6, 6, A_b_expand, A_b_expand_transpose);
    	// MatrixMulti_(6, 6, N, A_b_expand_transpose, J_g, A_b_expand_transpose_multi_J_g);
    	// MatrixTranspose_(6, N, A_b_expand_transpose_multi_J_g, A_b_expand_transpose_multi_J_g_transpose);
    	// MatrixMulti_(6, N, 6, A_b_expand_transpose_multi_J_g, A_b_expand_transpose_multi_J_g_transpose, manipl_temp);
        // manipl_temp_2 = calc_determinantOfMatrix(manipl_temp, 6*6, 6, 6, 6);

        // for(int i = 0; i < 6; i++)
        // {
        //     for(int j = 0; j < N; j++)
        //     {
        //         cout << J_g[i * N + j] << " ";
        //     }
        //     cout << endl;
        // }
    

        MatrixMulti_(6, N, 6, J_g, J_g_transpose, J_g_multi_J_g_transpose);

        // // for(int i = 0; i < 6; i++)
        // // {
        // //     for(int j = 0; j < N; j++)
        // //     {
        // //         cout << J_g_multi_J_g_transpose[i * N + j] << "  " << "  " << "  ";
        // //     }
        // //     cout << endl;
        // //     cout << endl;
        // // }

    	manipl_temp_2 = calc_determinantOfMatrix(J_g_multi_J_g_transpose, 6 * 6, 6, 6, 6);

        // if manipl_temp_2 is near to 0, the manipulability is near to 0, set it to a small value
        if(manipl_temp_2 < 0.000001)
        {
            manipl_temp_2 = 0.000001;
        }

    	*manipl = sqrt(manipl_temp_2);




    // ***************************************************  delete memory allocation  ***********************************************
        for ( int i = 0; i < 15; i++ )
            delete[] objects_temp[i];
        delete[] objects_temp;

        for ( int i = 0; i < 15; i++ )
            delete[] A_links_transform_i[i];
        delete[] A_links_transform_i;

        // delete objects, note that the objects is a double*** type variable
        for ( int i = 0; i < 15; i++ )
        {
            for ( int j = 0; j < nvrtx; j++ )
            {
                delete[] objects[i][j];
            }
            delete[] objects[i];
        }
        delete[] objects;


        // delete center2vertices_temp, note that the center2vertices_temp is a double** type variable
        for ( int i = 0; i < 15; i++ )
            delete[] center2vertices_temp[i];
        delete[] center2vertices_temp;

        delete[] dis_left2links;
        delete[] dis_base2links;
        delete[] dis_right2links;
        delete[] dis_link_1_2_links;
        delete[] dis_link_2_2_links;

        // delete bd_objects, note that the bd_objects is a gkPolytope* type variable
        delete[] bd_objects;

        // delete center2vertices_buff, note that the center2vertices_buff is a double*** type variable
        for ( int i = 0; i < 15; i++ )
        {
            for ( int j = 0; j < nvrtx; j++ )
            {
                delete[] center2vertices_buff[i][j];
            }
            delete[] center2vertices_buff[i];
        }
        delete[] center2vertices_buff;

    	delete[] q ;
        delete[] q_dot ;
        delete[] q_ddot ;
        delete[] quaternion_base_initial ;
        delete[] A_b ;
        delete[] xi_b_tempt ;
        delete[] xi_end_tempt ;
        delete[] A_links_transform ;
        delete[] A_b_multi_b_b ;
        delete[] r_b_tempt_1 ;
        delete[] r_b_tempt_2 ;
        delete[] r_b_tempt_2_temp ;
        delete[] a_tempt_i ;
        delete[] r_b_tempt_3 ;
        delete[] r_b ;
        delete[] a_tempt_k ;
        delete[] b_tempt_k ;
        delete[] l_tempt_k ;    
        delete[] A_links_transform_tempt_i ;
        delete[] A_links_transform_tempt_i_multi_a_tempt_i;
        delete[] A_links_transform_tempt_k ;
        delete[] A_links_transform_tempt_k_multi_l_tempt_k ;
        delete[] r_b_tempt_3_tempt ;
        delete[] r ;
        delete[] r_e ;
        delete[] Pe_initial ;
        delete[] A_links_transform_N_minus_1;
        delete[] r_e_tempt;
        delete[] b_tempt_N_minus_1 ;
        delete[] p ;
        delete[] A_links_transform_tempt_i_multi_Ez ;
        delete[] A_links_transform_tempt_i_multi_Ez_cross;
        delete[] Jm_v ;
        delete[] Jm_v_tempt ;
        delete[] Jm_v_tempt_i ;
        delete[] Jm_w ;
        delete[] Jm_w_tempt ;
        delete[] p_tempt ;
        delete[] r_e_minus_p_i ;
        delete[] Jm ;
        delete[] Jm_tempt;
        delete[] J_bm_w ;
        delete[] J_bm_v ;
        delete[] J_bm ;	
        delete[] J_bE ;
        delete[] J_g ;
        delete[] J_g_transpose;
        delete[] J_g_multi_J_g_transpose;
        delete[] J_bE_multi_J_bm ;
        delete[] J_g_v ;
        delete[] J_g_w ;    
        delete[] quaternion_end ;
        delete[] eta_b_dot_tempt_1 ;
        delete eta_b_dot_tempt_2 ;  
        delete[] xi_b_dot ;
        delete[] cross_xi_b ;  
        delete[] xi_b_dot_tempt_1 ;
        delete[] xi_b_dot_tempt_2 ;
        delete[] xi_b_dot_tempt_3 ; 
        delete[] eta_end_dot_tempt_1 ;
        delete eta_end_dot_tempt_2 ;
        delete[] xi_end_dot ;
        delete[] cross_xi_end; 
        delete[] xi_end_dot_tempt_1 ;
        delete[] xi_end_dot_tempt_2 ;
        delete[] xi_end_dot_tempt_3 ;
        delete[] v_e ;
        delete[] v_b ;
        delete[] A_b_expand ;
    	delete[] A_b_expand_transpose ;
    	delete[] A_b_expand_transpose_multi_J_g ;
    	delete[] A_b_expand_transpose_multi_J_g_transpose ;
    	delete[] manipl_temp ;
        delete[] xi_b_initial ;
        delete[] delta_xi_base ;
        delete delta_eta_base_tempt;
        delete delta_eta_base ;
        delete[] delta_xi_base_tempt;
    }


    void J_Base2EE(double* r_e, double* r_b, double* J_bE)
    {
        double* tempt;
        tempt = new double[9];

        double* r_b2r_e;
        r_b2r_e = new double[3];
    
    
    	for(int i = 0; i < 6; i++)
        {
    		for(int j = 0; j < 6; j++)
    			{
    				if(i == j)
    					J_bE[i*6+j] = 1;
    
    				else
    					J_bE[i*6+j] = 0;
    
    			}
        }
    
    	// double r_b2r_e[3] = {0};
    
    	for(int i = 0; i < 3; i++)
    		r_b2r_e[i] = r_e[i] - r_b[i];
    
    	cross(r_b2r_e, tempt);
    
    	for(int i = 0; i < 9; i++)
    		tempt[i] = tempt[i] * (-1);    // 得到 -cross(r_e - r_b)
    
    	for(int i = 0; i < 3; i++)
    		for(int j = 0; j < 3; j++)
    			J_bE[i * 6 + (3 + j)] = tempt[i * 3 + j];


        delete[] tempt;
        delete[] r_b2r_e;
    
    }


    void calc_J_bm(double* r, double* r_b, double* A_b, double* A_links_transform, double* p, double* J_bm_w, double* J_bm_v)
    {
    	/*
        :param I_b_body: 在基座本体系中表示的基座惯量
        :param I_links_body: 在各个连杆本体系中表示的连杆惯量
    	*/

        double* I_b = new double[9];
        double* I_b_tempt_2 = new double[9];
        double* I_b_tempt_1 = new double[9];
        double* r_g = new double[3];

    	double total_mass = m_b;

    	double* a_tempt_1 = new double[3];
    	double* a_tempt_2 = new double[3];
    	double* A_links_transform_tempt = new double[9];
    	double* A_links_transform_tempt_2 = new double[9];    
    
    	double* I_links_body_tempt = new double[9];
    	double* I_links_tempt_1 = new double[9];
    	double* I_links_tempt_2 = new double[9];
    	double* I_links = new double[N * 3 * 3];
    	double* I_links_i_tempt = new double[9];
    
        double* Hw = new double[9];
        double* r_i_b = new double[3];
        double* r_i_b_cross_1 = new double[9];
        double* r_i_b_cross_1_T = new double[9];
        double* r_i_b_cross_1_tempt = new double[9];
        double* r_i_b_cross_2 = new double[9];
        double* joint_revolute_axis_1 = new double[3];
        double* joint_revolute_axis_2 = new double[3];
        double* JR = new double[N * N * 3];
        double* JR_i_tempt = new double[3 * N];

        for(int i = 0; i < N; i++)
            for(int j = 0; j < 3; j++)
                for(int k = 0; k < N; k++)
                    JR[i * 3*N + j * N + k] = 0;

        double* JT = new double[N * N * 3];
        double* JT_i_tempt = new double[3 * N];

        for(int i = 0; i < N; i++)
            for(int j = 0; j < 3; j++)
                for(int k = 0; k < N; k++)
                    JT[i * 3*N + j * N + k] = 0;
    
        double* J_Tw = new double[3*N];
        double* J_Tw_tempt = new double[3*N];

        for(int i = 0; i < 3; i++)
            for(int j = 0; j < N ; j++)
                J_Tw[i * N + j] = 0;
    
        double* J_bm_w_tempt = new double[3*N];
        double* joint_revolute_axis_cross_r_i_minus_p_i = new double[3];
        double* joint_revolute_axis_cross_r_i_minus_p_j = new double[3];
        double* joint_revolute_axis_2_cross = new double[3*3];
        double* r_i_minus_p_j = new double[3];
        double* r_g_minus_r_b = new double[3];
        double* r_g_minus_r_b_cross = new double[9];
        double* r_g_minus_r_b_cross_multi_r_g_minus_r_b_cross = new double[9];
        double* r_g_minus_r_b_cross_multi_J_Tw = new double[3*N];
    
        double* H_wq = new double[3*N];
        double* H_wq_tempt_1 = new double[3*N];
        double* H_wq_tempt_2 = new double[3*N];
        double* H_wq_tempt_3 = new double[3*N];

        for(int i = 0; i < 3; i++)
            for(int j = 0; j < N ; j++)
                H_wq[i * N + j] = 0;

    
        double* H_s = new double[3*3];
    	double* H_q = new double[3*N];
        double* r_g_minus_r_b_cross_multi_H_s_inv = new double[3*3];
        double* r_g_minus_r_b_cross_multi_H_s_inv_multi_H_q = new double[3*N];
        double* H_s_inv = new double[3*3];
        double* H_s_inv_tempt = new double[3*3];
    
    	// 计算惯性系中表示的基座惯量

    	MatrixMultiTranspose(I_b_body, A_b, I_b_tempt_1);
    
    	MatrixMulti_(3, 3, 3, A_b, I_b_tempt_1, I_b_tempt_2);
    
    	MatrixCopy_( 3, 3, I_b_tempt_2, I_b );
    
    	// 计算系统质心位置
    
    	ScaleMatrix_( 1, 3, m_b, r_b, r_g);
    
    
    	for(int i=0;i<N;i++)
    	{
    		total_mass += m[i];  //  total_mass 初始化为 m_b
    
    		for(int j=0;j<3;j++)
    			r_g[j] = r_g[j] + m[i] * r[i * 3 + j];
    	}
    
    	ScaleMatrix_( 1, 3, 1/total_mass, r_g, r_g);
    
        MatrixSub_( 1, 3, r_g, r_b, r_g_minus_r_b);
    
    	cross(r_g_minus_r_b, r_g_minus_r_b_cross);
    
    
    	/*   // OK
    	//-----------------------------------------------------------------------------------------------------------------------------------
    	cout << "r_g_minus_r_b_cross:" << "  " 
    	               << r_g_minus_r_b_cross[0] << "  "
    	               << r_g_minus_r_b_cross[1] << "  "
    	               << r_g_minus_r_b_cross[2] << endl;
    	               cout << r_g_minus_r_b_cross[3] << "  "
    	               << r_g_minus_r_b_cross[4] << "  "
    	               << r_g_minus_r_b_cross[5] << endl;
    	               cout << r_g_minus_r_b_cross[6] << "  "
    	               << r_g_minus_r_b_cross[7] << "  "
    	               << r_g_minus_r_b_cross[8] << endl;
    
    	//------------------------------------------------------------------------------------------------------------------------------------
    	*/
    
    

    
    
    	//  初始化 Hw
    	for(int j=0;j<3;j++)
    		for(int k=0;k<3;k++)
    			Hw[j * 3 + k] = I_b[j * 3 + k];
    
    
    	for(int i = 0; i < N; i++)
    	{
    		for(int j=0;j<3;j++)
    			r_i_b[j] = r[i * 3 + j] - r_b[j];
    
    		cross(r_i_b, r_i_b_cross_1);  // r_i_b cross
    
    		for(int j=0;j<3;j++)
            {
    			for(int k=0;k<3;k++)
    			{
    				r_i_b_cross_1_tempt[j * 3 + k] = r_i_b_cross_1[j * 3 + k];
    
    				A_links_transform_tempt[j * 3 + k] = A_links_transform[i * 9 + j * 3 + k];  // A_links_transform_tempt  第 i 个连杆的旋转矩阵
    
    				I_links_body_tempt[j * 3 + k] = I_links_body[i * 9 + j * 3 + k];
    			}
            }
    
    		MatrixTranspose_(3, 3, r_i_b_cross_1, r_i_b_cross_1_T);    // MatrixTranspose(r_i_b_cross_1);  // (r_i -r_b)叉 转置
    		MatrixMulti_(3, 3, 3, r_i_b_cross_1_T, r_i_b_cross_1_tempt, r_i_b_cross_2);  //  r_i_b_cross_2 = (r_i -r_b)叉.T * (r_i -r_b)叉
    
    		MatrixMulti_(3, 3, 3, A_links_transform_tempt, I_links_body_tempt, I_links_tempt_1);
    
    		MatrixMultiTranspose(I_links_tempt_1, A_links_transform_tempt, I_links_tempt_2);  // I_links_tempt_2   第 i 个连杆 在 惯性系 下的 惯量
    
    
    		for(int j=0;j<3;j++)
    		{	
    			for(int k=0;k<3;k++)
    		    {
    				I_links[i * 9 + j * 3 + k] = I_links_tempt_2[j * 3 + k];
    
    				I_links_i_tempt[j * 3 + k] = I_links[i * 9 + j * 3 + k];  //  缓存 第 i 个连杆的 惯量
    
    	            /*
                        Hw = I_b
                        for i in range(N):
                            Hw += I_links[i, :, :] + m[i] * cross(r[i, :, :] - r_b).T @ cross(r[i, :, :] - r_b)
                    */
               
    				Hw[j * 3 + k] += I_links_i_tempt[j * 3 + k] + m[i] * r_i_b_cross_2[j * 3 + k];
    
    			}
    		}
    
    		MatrixMulti_(3, 3, 1, A_links_transform_tempt, Ez, joint_revolute_axis_1);  //  joint_revolute_axis_1  第 i 个关节的旋转轴
    
    		for(int j = 0; j < (i+1); j++)
    		{
    			for(int jj = 0; jj < 3; jj++)
    			{
    				r_i_minus_p_j[jj] = r[i * 3 + jj] - p[j * 3 + jj];    // 求 r_i - p_j （r_i_minus_p_j） 
    
    				for(int kk = 0; kk < 3; kk++)
    					A_links_transform_tempt_2[jj * 3 + kk] = A_links_transform[j * 9 + jj * 3 + kk];  //  缓存 第 j 个连杆的旋转矩阵
    
    			}
    
    			MatrixMulti_(3, 3, 1, A_links_transform_tempt_2, Ez, joint_revolute_axis_2);  // joint_revolute_axis_2  第 j 个关节的旋转轴
    
                cross(joint_revolute_axis_2, joint_revolute_axis_2_cross);				

    			MatrixMulti_(3, 3, 1, joint_revolute_axis_2_cross, r_i_minus_p_j, joint_revolute_axis_cross_r_i_minus_p_j);
    
    			for(int k = 0; k < 3; k++)
    			{	
    				if(j == i)
    				{
    					JR[i * (3*N) + k * N + j] = joint_revolute_axis_1[k];
    					JT[i * (3*N) + k * N + j] = joint_revolute_axis_cross_r_i_minus_p_j[k];
    				}
    				else
    				{
    					JR[i * (3*N) + k * N + j] = JR[(i-1) * (3*N) + k * N + j];					
    					JT[i * (3*N) + k * N + j] = joint_revolute_axis_cross_r_i_minus_p_j[k];			
    				}
    			}

    		}
    
    		for(int j = 0; j < 3; j++)
            {
    			for(int k = 0; k < N; k++)
    			{
    				JR_i_tempt[j * N + k] = JR[i * (3 * N) + j * N + k];  //  缓存 JR_i
    				JT_i_tempt[j * N + k] = JT[i * (3 * N) + j * N + k];  //  缓存 JT_i
    
    			}
            }
    
    		MatrixMulti_(3, 3, N, I_links_i_tempt, JR_i_tempt, H_wq_tempt_1);
    
    		MatrixMulti_(3, 3, N, r_i_b_cross_1_tempt, JT_i_tempt, H_wq_tempt_2);
    
    		ScaleMatrix_( 3, N, m[i], H_wq_tempt_2, H_wq_tempt_2);		
    
    		MatrixAdd_( 3, N, H_wq_tempt_1, H_wq_tempt_2, H_wq_tempt_3);	
    
    		MatrixAdd_( 3, N, H_wq_tempt_3, H_wq, H_wq);
    
    		ScaleMatrix_( 3, N, m[i], JT_i_tempt, JT_i_tempt);
    
    		MatrixAdd_( 3, N, JT_i_tempt, J_Tw, J_Tw);
    
    
    	}
    
    	/*    // OK
    	//-----------------------------------------------------------------------------------------------------------------------------------
    	cout << "Hw:" << "  " 
    	              << Hw[0] << "  "
    	              << Hw[1] << "  "
    	              << Hw[2] << endl;
    	         cout << Hw[3] << "  "
    	              << Hw[4] << "  "
    	              << Hw[5] << endl;
    	         cout << Hw[6] << "  "
    	              << Hw[7] << "  "
    	              << Hw[8] << endl;
    
    	//------------------------------------------------------------------------------------------------------------------------------------
    	*/


    
    
    
    
    	MatrixMulti_(3, 3, 3, r_g_minus_r_b_cross, r_g_minus_r_b_cross, r_g_minus_r_b_cross_multi_r_g_minus_r_b_cross);
    
    	ScaleMatrix_( 3, 3, total_mass, r_g_minus_r_b_cross_multi_r_g_minus_r_b_cross, r_g_minus_r_b_cross_multi_r_g_minus_r_b_cross);
    
    	MatrixAdd_( 3, 3, r_g_minus_r_b_cross_multi_r_g_minus_r_b_cross, Hw, H_s);		
    
    
    	/*
    	//-----------------------------------------------------------------------------------------------------------------------------------
    	cout << "H_s:" << "  " 
    	              << H_s[0] << "  "
    	              << H_s[1] << "  "
    	              << H_s[2] << endl;
    	         cout << H_s[3] << "  "
    	              << H_s[4] << "  "
    	              << H_s[5] << endl;
    	         cout << H_s[6] << "  "
    	              << H_s[7] << "  "
    	              << H_s[8] << endl;
    
    	//------------------------------------------------------------------------------------------------------------------------------------
    	*/
    
    
    
    
    	// Hq = Hwq - cross(r_g - r_b) @ JTw
    	MatrixMulti_(3, 3, N, r_g_minus_r_b_cross, J_Tw, r_g_minus_r_b_cross_multi_J_Tw);
    
    	/*
    	for(int kk = 0; kk < 3; kk++)
    		for(int jj = 0; jj < N; jj++)
    			H_q[kk * N + jj] = H_wq[kk * N + jj] - r_g_minus_r_b_cross_multi_J_Tw[kk * N + jj];
    
    	*/		
    
    	MatrixSub_( 3, N, H_wq, r_g_minus_r_b_cross_multi_J_Tw, H_q);


        /*    // OK
    	//-----------------------------------------------------------------------------------------------------------------------------------
    	cout << "H_q:" << "  " 
    	              << H_q[0] << "  "
    	              << H_q[1] << "  "
    	              << H_q[2] << "  "
    	              << H_q[3] << "  "
    	              << H_q[4] << "  "
    	              << H_q[5] << "  "
    	              << H_q[6] << endl;
    	         cout << H_q[7] << "  "
    	              << H_q[8] << "  "
    	              << H_q[9] << "  "
    	              << H_q[10] << "  "
    	              << H_q[11] << "  "
    	              << H_q[12] << "  "
    	              << H_q[13] << endl;
    	         cout << H_q[14] << "  "
    	              << H_q[15] << "  "
    	              << H_q[16] << "  "
    	              << H_q[17] << "  "
    	              << H_q[18] << "  "
    	              << H_q[19] << "  "
    	              << H_q[20] << endl;
    
    	//------------------------------------------------------------------------------------------------------------------------------------
    	*/
    
    
    
    	//J_bm_w = (-1) * linalg.inv(Hs) @ Hq
        //J_bm_v = (-1) * (JTw / (m_b + sum(m)) + cross(r_g - r_b) @ linalg.inv(Hs) @ Hq)

    
    	LUP_solve_inverse(H_s, 3, H_s_inv);	
    
    	MatrixCopy_( 3, 3, H_s_inv, H_s_inv_tempt );
    
    	/*
    	for(int i=0;i<3;i++)
    		for(int j=0;j<3;j++)
    			H_s_inv_tempt[i * 3 + j] = H_s_inv[i * 3 + j];
    	*/
    
    	MatrixMulti_(3, 3, N, H_s_inv_tempt, H_q, J_bm_w_tempt);
    
    	ScaleMatrix_( 3, N, (-1), J_bm_w_tempt, J_bm_w);
    
    	ScaleMatrix_( 3, N, 1/total_mass, J_Tw, J_Tw_tempt);
    
    	MatrixMulti_(3, 3, 3, r_g_minus_r_b_cross, H_s_inv_tempt, r_g_minus_r_b_cross_multi_H_s_inv);
    	MatrixMulti_(3, 3, N, r_g_minus_r_b_cross_multi_H_s_inv, H_q, r_g_minus_r_b_cross_multi_H_s_inv_multi_H_q);
    
    	MatrixAdd_( 3, N, J_Tw_tempt, r_g_minus_r_b_cross_multi_H_s_inv_multi_H_q, J_bm_v);
    
    	ScaleMatrix_( 3, N, (-1), J_bm_v, J_bm_v);
    
    
    	/*	
    	//-------------------OK--------------------------------------------------------------------------------------------------
    	for(int i=0;i<3;i++)
    		for(int j=0;j<N;j++)
    			cout << J_bm_v[i * N + j] << endl; 
    	//-----------------------------------------------------------------------------------------------------------------------
    
    	*/
    
    

        delete[] I_b;
        delete[] I_b_tempt_2;
        delete[] I_b_tempt_1;
        delete[] r_g;
    	delete[] a_tempt_1;
    	delete[] a_tempt_2;
    	delete[] A_links_transform_tempt;
    	delete[] A_links_transform_tempt_2;
    	delete[] I_links_body_tempt;
    	delete[] I_links_tempt_1;
    	delete[] I_links_tempt_2;
    	delete[] I_links;
    	delete[] I_links_i_tempt;
        delete[] Hw;	
        delete[] r_i_b;
        delete[] r_i_b_cross_1;
        delete[] r_i_b_cross_1_T;
        delete[] r_i_b_cross_1_tempt;
        delete[] r_i_b_cross_2;
        delete[] joint_revolute_axis_1;
        delete[] joint_revolute_axis_2;
        delete[] JR;
        delete[] JR_i_tempt;
        delete[] JT;
        delete[] JT_i_tempt;
        delete[] J_Tw;
        delete[] J_Tw_tempt;
        delete[] J_bm_w_tempt;
        delete[] joint_revolute_axis_cross_r_i_minus_p_i;
        delete[] joint_revolute_axis_cross_r_i_minus_p_j;
        delete[] joint_revolute_axis_2_cross;
        delete[] r_i_minus_p_j;
        delete[] r_g_minus_r_b;
        delete[] r_g_minus_r_b_cross;
        delete[] r_g_minus_r_b_cross_multi_r_g_minus_r_b_cross;
        delete[] r_g_minus_r_b_cross_multi_J_Tw;
        delete[] H_wq;
        delete[] H_wq_tempt_1;
        delete[] H_wq_tempt_2;
        delete[] H_wq_tempt_3;
        delete[] H_s;
    	delete[] H_q;	
        delete[] r_g_minus_r_b_cross_multi_H_s_inv;
        delete[] r_g_minus_r_b_cross_multi_H_s_inv_multi_H_q;
        delete[] H_s_inv;
        delete[] H_s_inv_tempt;
    
    }


    void links_transform(double* A_base, double* q, double* A_links_transform)
    {
    	double* tempt_rpy2dc_joints_initial = new double[3*3];
    	double* tempt_rpy2dc_joints_var = new double[3*3];
    	double* tempt_A_links_transform_pre = new double[3*3];
    	double* tempt_A_links_transform_curr = new double[3*3];
    	double* tempt = new double[3*3];
    
    	for(int i = 0; i < N; i++)
    	{
    		if (i == 0)
    		{
    			rpy2dc(rpy_joints[i * 3], rpy_joints[i * 3 + 1], rpy_joints[i * 3 + 2], tempt_rpy2dc_joints_initial);

    			MatrixMulti_(3, 3, 3, A_base, tempt_rpy2dc_joints_initial, tempt);
    
    			rpy2dc(0, 0, q[i], tempt_rpy2dc_joints_var);
    
    			MatrixMulti_(3, 3, 3, tempt, tempt_rpy2dc_joints_var, tempt_A_links_transform_curr);
    
    			for(int j = 0; j< 3; j++)
                {
    				for(int k = 0; k < 3; k++)
    				{
    					A_links_transform[i * 9 + j * 3 + k] = tempt_A_links_transform_curr[j * 3 + k];

    				}
                }

    		}
    
    		else
    		{
    			rpy2dc(rpy_joints[i * 3], rpy_joints[i * 3 + 1], rpy_joints[i * 3 + 2], tempt_rpy2dc_joints_initial);
    
    			for(int j = 0; j < 3; j++)
    				for(int k = 0; k < 3; k++)
    					tempt_A_links_transform_pre[j * 3 + k] = A_links_transform[(i-1) * 9 + j * 3 + k];
    
    			MatrixMulti_(3, 3, 3, tempt_A_links_transform_pre, tempt_rpy2dc_joints_initial, tempt);
    
    			rpy2dc(0, 0, q[i], tempt_rpy2dc_joints_var);
    
    			MatrixMulti_(3, 3, 3, tempt, tempt_rpy2dc_joints_var, tempt_A_links_transform_curr);
    
    			for(int j = 0; j< 3; j++)
                {
    				for(int k = 0; k < 3; k++)
    				{
    					A_links_transform[i * 9 + j * 3 + k] = tempt_A_links_transform_curr[j * 3 + k];

    				}
                }

    		}
    
    	}
    
    	delete[] tempt_rpy2dc_joints_initial;
    	delete[] tempt_rpy2dc_joints_var;
    	delete[] tempt_A_links_transform_pre;
    	delete[] tempt_A_links_transform_curr;	
    	delete[] tempt;

    }


    int readinput(const char *inputfile, double ***pts, int *out) 
    {
        int npoints = 0;

        int idx = 0;

        FILE *fp;

        /* Open file. */
        /*
            C 库函数 :
            FILE *fopen(const char *filename, const char *mode) 
            使用给定的模式 mode 打开 filename 所指向的文件。
            下面是 fopen() 函数的声明:
            FILE *fopen(const char *filename, const char *mode)

            filename -- 字符串，表示要打开的文件名称。
            mode -- 字符串，表示文件的访问模式，可以是以下表格中的值：
            "r"	打开一个用于读取的文件。该文件必须存在。
            "w"	创建一个用于写入的空文件。如果文件名称与已存在的文件相同，则会删除已有文件的内容，文件被视为一个新的空文件。
            "a"	追加到一个文件。写操作向文件末尾追加数据。如果文件不存在，则创建文件。
            "r+"	打开一个用于更新的文件，可读取也可写入。该文件必须存在。
            "w+"	创建一个用于读写的空文件。
            "a+"	打开一个用于读取和追加的文件。

            该函数返回一个 FILE 指针。否则返回 NULL，且设置全局变量 errno 来标识错误。

        */

        fp = fopen(inputfile, "r");

        if ((fp) == NULL) 
        {
            fprintf(stdout, "ERROR: input file %s not found!\n", inputfile);
            fprintf(stdout, "  -> The file must be in the folder from which this " "program is launched\n\n");
            return 1;
        }

        /* Read number of input vertices. */
        if (fscanf_s(fp, "%d", &npoints) != 1)
            return 1;

        /* Allocate memory. */
        double **arr = (double **)malloc(npoints * sizeof(double *));
        for (int i = 0; i < npoints; i++)
            arr[i] = (double *)malloc(3 * sizeof(double));

        /* Read and store vertices' coordinates. */
        for (idx = 0; idx < npoints; idx++) 
        {
          if (fscanf_s(fp, "%lf %lf %lf\n", &arr[idx][0], &arr[idx][1], &arr[idx][2]) != 3)
              return 1;
        }

        fclose(fp);

        *pts = arr;
        *out = idx;

        return (0);
    }






}