#ifndef CALC_UTILS
#define CALC_UTILS


#include <stdio.h>
#include <cmath>
#include <iostream>
#include <ctime>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <random>
#include <chrono>


namespace utils{

    // generate random double number
    void random_double_generation(double* random, int num, double min, double max){
        std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
        std::uniform_real_distribution<double> distribution(min, max);
        for(int i = 0; i < num; ++i){
            random[i] = distribution(generator);
        }

        // srand((unsigned)time(NULL));
        // for (int i = 0; i < num; i++)
		//     random[i] = (rand()/double(RAND_MAX) * (max - min) + min);

    }

    void cxyz(double theta, int x, int y, int z, double* direction_cosines)
    {
        if (x == 1)
	    {
	    	direction_cosines[0] = 1; direction_cosines[1] = 0;          direction_cosines[2] = 0;
    
	    	direction_cosines[3] = 0; direction_cosines[4] = cos(theta); direction_cosines[5] = -sin(theta);
    
	    	direction_cosines[6] = 0; direction_cosines[7] = sin(theta); direction_cosines[8] = cos(theta);
    
	    }
	    else if (y == 1)
	    {
	    	direction_cosines[0] = cos(theta);  direction_cosines[1] = 0; direction_cosines[2] = sin(theta);
    
	    	direction_cosines[3] = 0;           direction_cosines[4] = 1; direction_cosines[5] = 0;
    
	    	direction_cosines[6] = -sin(theta); direction_cosines[7] = 0; direction_cosines[8] = cos(theta);
    
	    }
	    else if (z == 1)
	    {
	    	direction_cosines[0] = cos(theta); direction_cosines[1] = -sin(theta); direction_cosines[2] = 0;
    
	    	direction_cosines[3] = sin(theta); direction_cosines[4] = cos(theta);  direction_cosines[5] = 0;
    
	    	direction_cosines[6] = 0;          direction_cosines[7] = 0;           direction_cosines[8] = 1;
    
	    }
    }


    //两个 3x3 矩阵相乘
    void MatrixMulti(double* A, double* B, double* AB)
    {
	    for(int i = 0; i < 3; i++)
	    	for(int j = 0; j < 3; j++)
	    		AB[i * 3 + j] = 0;
    
	    for(int i = 0; i < 3; i++)
	    {
	    	for(int j = 0; j < 3; j++)
	    	{
	    		if (i == 0)
	    			AB[i + j] = A[i] * B[j] + A[i + 1] * B[j + 3] + A[i + 2] * B[j + 6];
    
	    		else if (i == 1)
	    			AB[i + 2 + j] = A[i + 2] * B[j] + A[i + 2 + 1] * B[j + 3] + A[i + 2 + 2] * B[j + 6];
    
	    		else
	    			AB[i + 2 + 2 + j] = A[i + 2 + 2] * B[j] + A[i + 2 + 2 + 1] * B[j + 3] + A[i + 2 + 2 + 2] * B[j + 6];
    
	    	}
    
	    }
    }


    void MatrixMulti_(int row_1, int col_1, int col_2, double* A, double* B, double* AB)
    {
        /*
         两个 任意维数 矩阵相乘
         row_1 第一个矩阵的行数
         col_1 第一个矩阵的列数
         col_1 第二个矩阵的行数和第一个矩阵的列数相等
         col_2 第二个矩阵的列数
        */

        double* A_tempt;
        A_tempt = new double[row_1 * col_1];

        double* B_tempt;
        B_tempt = new double[col_1 * col_2];

        double* AB_tempt;
        AB_tempt = new double[row_1 * col_2];


        for(int i = 0; i < row_1; i++)
        	for(int j = 0; j < col_1; j++)
        		A_tempt[i * col_1 + j] = A[i * col_1 + j];

        for(int i = 0; i < col_1; i++)
        	for(int j = 0; j < col_2; j++)
        		B_tempt[i * col_2 + j] = B[i * col_2 + j];


        for(int i = 0; i < row_1; i++)
        	for(int j = 0; j < col_2; j++)
        		AB_tempt[i * col_2 + j] = 0;


        for(int i = 0; i < row_1; i++)
        	for(int j = 0; j < col_2; j++)
        		for(int k = 0; k < col_1; k++)
        			AB_tempt[i * col_2 + j] += A_tempt[i * col_1 + k] * B_tempt[k * col_2 + j];

        for(int i = 0; i < row_1; i++)
        	for(int j = 0; j < col_2; j++)
        		AB[i * col_2 + j] = AB_tempt[i * col_2 + j];

        delete[] A_tempt;
        delete[] B_tempt;
        delete[] AB_tempt;
    }

    void zyx2quaternion(double psi, double theta, double phi, double* quaternion)
    {
        /*
	     欧拉角转四元素，其中 欧拉角 用 动系相继旋转 Z-Y-X法表示
         psi: 表示 Z 轴转动量
         theta: 表示 Y 轴转动量
         phi: 表示 X 轴转动量
         return: 四元素，q0标量，q1,2,3向量
	    */
    
	    quaternion[0] = cos(phi/2) * cos(theta/2) * cos(psi/2) + sin(phi/2) * sin(theta/2) * sin(psi/2);
	    quaternion[1] = sin(phi/2) * cos(theta/2) * cos(psi/2) - cos(phi/2) * sin(theta/2) * sin(psi/2);
	    quaternion[2] = cos(phi/2) * sin(theta/2) * cos(psi/2) + sin(phi/2) * cos(theta/2) * sin(psi/2);
	    quaternion[3] = cos(phi/2) * cos(theta/2) * sin(psi/2) - sin(phi/2) * sin(theta/2) * cos(psi/2);
    }

    void quaternion2dc(double q_0, double q_1, double q_2, double q_3, double* dc)
    {
        dc[0] = pow(q_0, 2) + pow(q_1, 2) - pow(q_2, 2) - pow(q_3, 2);
	    dc[1] = 2 * (q_1 * q_2 - q_0 * q_3);
	    dc[2] = 2 * (q_1 * q_3 + q_0 * q_2);
	    dc[3] = 2 * (q_1 * q_2 + q_0 * q_3);
	    dc[4] = pow(q_0, 2) - pow(q_1, 2) + pow(q_2, 2) - pow(q_3, 2);
	    dc[5] = 2 * (q_2 * q_3 - q_0 * q_1);
	    dc[6] = 2 * (q_1 * q_3 - q_0 * q_2);
	    dc[7] = 2 * (q_2 * q_3 + q_0 * q_1);
	    dc[8] = pow(q_0, 2) - pow(q_1, 2) - pow(q_2, 2) + pow(q_3, 2);
    }

    void rpy2dc(double roll, double pitch, double yaw, double* direction_cosine)
    {
        double* tempt_yaw = new double[9];
	    double* tempt_pitch = new double[9];
	    double* tempt_roll = new double[9];
	    double* tempt_yawMulpitch = new double[9];
	    cxyz(yaw, 0, 0, 1, tempt_yaw);
	    cxyz(pitch, 0, 1, 0, tempt_pitch);
	    cxyz(roll, 1, 0, 0, tempt_roll);
	    MatrixMulti_(3, 3, 3, tempt_yaw, tempt_pitch, tempt_yawMulpitch);
	    MatrixMulti_(3, 3, 3, tempt_yawMulpitch, tempt_roll, direction_cosine);
	    delete[] tempt_yaw;
	    delete[] tempt_pitch;
	    delete[] tempt_roll;
	    delete[] tempt_yawMulpitch;
    }

    void dc2euler(double* R, double* euler)
    {
        double theta_y_1;
	    double theta_y_2;
	    double theta_x_1;
	    double theta_x_2;
	    double theta_z_1;
	    double theta_z_2;
    
	    if(R[2 * 3 + 0] != 1 && R[2 * 3 + 0] != (-1))
	    {
	    	theta_y_1 = -asin(R[2 * 3 + 0]);
	    	theta_y_2 = M_PI - theta_y_1;
    
	    	theta_x_1 = atan2(R[2 * 3 + 1] / cos(theta_y_1), R[2 * 3 + 2] / cos(theta_y_1));
	    	theta_x_2 = atan2(R[2 * 3 + 1] / cos(theta_y_2), R[2 * 3 + 2] / cos(theta_y_2));
    
	    	theta_z_1 = atan2(R[1 * 3 + 0] / cos(theta_y_1), R[0 * 3 + 0] / cos(theta_y_1));
            theta_z_2 = atan2(R[1 * 3 + 0] / cos(theta_y_2), R[0 * 3 + 0] / cos(theta_y_2));
    
	    }
    
	    else if(R[2 * 3 + 0] == (-1))
	    {
	    	theta_z_1 = 0;
            theta_z_2 = 0;
	    	theta_y_1 = M_PI / 2;
	    	theta_y_2 = M_PI / 2;
    
	    	theta_x_1 = theta_z_1 + atan2(R[0 * 3 + 1], R[0 * 3 + 2]);
	    	theta_x_2 = theta_z_1 + atan2(R[0 * 3 + 1], R[0 * 3 + 2]);
	    	theta_x_1 = theta_z_1 + atan2(R[0 * 3 + 1], R[0 * 3 + 2]);
            theta_x_2 = theta_z_1 + atan2(R[0 * 3 + 1], R[0 * 3 + 2]);
    
	    }
    
	    else
	    {
	    	theta_z_1 = 0;
            theta_z_2 = 0;
            theta_y_1 = -M_PI / 2;
            theta_y_2 = -M_PI / 2;
            theta_x_1 = -theta_z_1 + atan2(-R[0 * 3 + 1], -R[0 * 3 + 2]);
            theta_x_2 = -theta_z_1 + atan2(-R[0 * 3 + 1], -R[0 * 3 + 2]);
	    }
    
        euler[0] = theta_z_1; 
        euler[1] = theta_y_1;
        euler[2] = theta_x_1;
        euler[3] = theta_z_2;
        euler[4] = theta_y_2;
        euler[5] = theta_x_2;
	
    }

    void euler_ZYX2dc(double alpha, double beta, double gamma, double* A)
    {
        // A 是一个 一维数组（3 * 3），用这个一维数组来表示二维数组，元素用 i * N + j 来查找
        // Z-Y-X 欧拉角[alpha, beta, gamma]定义为：frame_j分别绕其 Z 轴、旋转后的 Y 轴、再旋转后的 X 轴旋转alpha, beta, gamma角后, 与frame_i重合
	    A[0 * 3 + 0] = cos(alpha) * cos(beta);
	    A[0 * 3 + 1] = cos(alpha) * sin(beta) * sin(gamma) - sin(alpha) * cos(gamma);
        A[0 * 3 + 2] = cos(alpha) * sin(beta) * cos(gamma) + sin(alpha) * sin(gamma);

        A[1 * 3 + 0] = sin(alpha) * cos(beta);
        A[1 * 3 + 1] = sin(alpha) * sin(beta) * sin(gamma) + cos(alpha) * cos(gamma);
        A[1 * 3 + 2] = sin(alpha) * sin(beta) * cos(gamma) - cos(alpha) * sin(gamma);

        A[2 * 3 + 0] = -sin(beta);
        A[2 * 3 + 1] = cos(beta) * sin(gamma);
        A[2 * 3 + 2] = cos(beta) * cos(gamma);
	
    }

    void quaternion2euler(double w, double x, double y, double z, double* euler)
    {
        double Epsilon = 0.0009765625;
        double Threshold = 0.5 - Epsilon;

        double TEST = w * y - x * z;

        double gamma;
        double beta;
        double alpha;

        if(TEST < -Threshold || TEST > Threshold)
        {
            int sign;
            if(TEST < 0)
                sign = (-1);
            else
                sign = (1);

            gamma = -2 * sign * (double)atan2(x, w); // yaw

            beta = sign * (M_PI / 2.0); // pitch

            alpha = 0; // roll
        }
        else
        {
            alpha = atan2(2 * (y * z + w * x), w * w - x * x - y * y + z * z);
            beta = asin(-2 * (x * z - w * y));
            gamma = atan2(2 * (x * y + w * z), w * w + x * x - y * y - z * z);
        }

        euler[0] = alpha;
        euler[1] = beta;
        euler[2] = gamma; // 绕固定轴转
    
    }

    void MatrixMultiTranspose(double* A, double* B, double* AB)
    {
        // 矩阵 A 乘以 B.T， 两个矩阵都是 3 * 3 
	
	    for(int i = 0; i < 3; i++)
	    	for(int j = 0; j < 3; j++)
	    		AB[i * 3 + j] = 0;
    
	    for(int i = 0; i < 3; i++)
	    {
	    	for(int j = 0; j < 3; j++)
	    	{
	    		AB[i * 3 + j] = A[i * 3] * B[j * 3] + A[i * 3 + 1] * B[j * 3 + 1] + A[i * 3 + 2] * B[j * 3 + 2];
	    	}
    
	    }
				
    }

    void MatrixTranspose_(int row, int col, double* M, double* M_T)
    {
        // M 是一个 任意 的 m * n 的矩阵, 输出 M 矩阵 的 转置
	
        int m = row;
	    int n = col;

        double* M_tempt;
        M_tempt = new double[m * n];

        double* M_T_tempt;
        M_T_tempt = new double[n * m];

	    for(int i = 0; i < m; i++)
        	for(int j = 0; j < n; j++)
        		M_tempt[i * n + j] = M[i * n + j];
    
	    for(int p = 0; p < n; p++)
	    	for(int q = 0; q < m; q++)
	    		M_T_tempt[p * m + q] = M_tempt[q * n + p];

        for(int i = 0; i < n; i++)
        	for(int j = 0; j < m; j++)
        		M_T[i * m + j] = M_T_tempt[i * m + j];

        delete[] M_tempt;
        delete[] M_T_tempt;
    }

    void MatrixTranspose(double* M)
    {
        // double M_tempt[9] = {0};
        double* M_tempt;
        M_tempt = new double[3*3];
    
	    for(int i=0;i<3;i++)
	    	for(int j=0;j<3;j++)
	    		M_tempt[i * 3 + j] = M[i * 3 + j];
    
	    M[0] = M_tempt[0];
	    M[1] = M_tempt[3];
	    M[2] = M_tempt[6];
    
	    M[3] = M_tempt[1];
	    M[4] = M_tempt[4];
	    M[5] = M_tempt[7];
    
	    M[6] = M_tempt[2];
	    M[7] = M_tempt[5];
	    M[8] = M_tempt[8];

        delete[] M_tempt;
    }

    void MatrixAdd_(int row, int col, double* A, double* B, double* AB)
    {
        double* A_tempt = new double[row * col];
        double* B_tempt = new double[row * col];
        double* AB_tempt = new double[row * col];

        for(int i = 0; i < row; i++)
        	for(int j = 0; j < col; j++)
        		A_tempt[i * col + j] = A[i * col + j];

        for(int i = 0; i < row; i++)
        	for(int j = 0; j < col; j++)
        		B_tempt[i * col + j] = B[i * col + j];

        for(int i = 0; i < row; i++)
        	for(int j = 0; j < col; j++)
        		AB_tempt[i * col + j] = 0;

	    for(int i = 0; i < row; i++)
	    	for(int j = 0; j < col; j++)
	    		AB_tempt[i * col + j] = A_tempt[i * col + j] + B_tempt[i * col + j];

        for(int i = 0; i < row; i++)
        	for(int j = 0; j < col; j++)
        		AB[i * col + j] = AB_tempt[i * col + j];

        delete[] A_tempt;
        delete[] B_tempt;
        delete[] AB_tempt;
    
    }

    void MatrixSub_(int row, int col, double* A, double* B, double* AB)
    {
        double* A_tempt = new double[row * col];
        double* B_tempt = new double[row * col];
        double* AB_tempt = new double[row * col];

        for(int i = 0; i < row; i++)
        	for(int j = 0; j < col; j++)
        		A_tempt[i * col + j] = A[i * col + j];

        for(int i = 0; i < row; i++)
        	for(int j = 0; j < col; j++)
        		B_tempt[i * col + j] = B[i * col + j];

        for(int i = 0; i < row; i++)
        	for(int j = 0; j < col; j++)
        		AB_tempt[i * col + j] = 0;

	    for(int i = 0; i < row; i++)
	    	for(int j = 0; j < col; j++)
	    		AB_tempt[i * col + j] = A_tempt[i * col + j] - B_tempt[i * col + j];

        for(int i = 0; i < row; i++)
        	for(int j = 0; j < col; j++)
        		AB[i * col + j] = AB_tempt[i * col + j];

        delete[] A_tempt;
        delete[] B_tempt;
        delete[] AB_tempt;
    
    }

    void ScaleMatrix_(int row, int col, double scale, double* A, double* scale_A)
    {
        double* A_tempt = new double[row * col];
        double* AB_tempt = new double[row * col];

        for(int i = 0; i < row; i++)
        	for(int j = 0; j < col; j++)
        		A_tempt[i * col + j] = A[i * col + j];

        for(int i = 0; i < row; i++)
        	for(int j = 0; j < col; j++)
        		AB_tempt[i * col + j] = 0;

	    for(int i = 0; i < row; i++)
	    	for(int j = 0; j < col; j++)
	    		AB_tempt[i * col + j] = scale * A_tempt[i * col + j];

        for(int i = 0; i < row; i++)
        	for(int j = 0; j < col; j++)
        		scale_A[i * col + j] = AB_tempt[i * col + j];

        delete[] A_tempt;
        delete[] AB_tempt;
    
    }

    void SetZeroMatrix_(int row, int col, double* A)
    {
        for(int i = 0; i < row; i++)
    	    for(int j = 0; j < col; j++)
    		    A[i * col + j] = 0;

    }

    void MatrixExtract_(int m, int n, int i_s, int i_f, int j_s, int j_f, double *A, double *ans)
    {
        int row, col;
        row = i_f - (i_s-1);
        col = j_f - (j_s-1);
        for(int i = (i_s-1); i < i_f; i++)
        {
            for(int j = (j_s-1); j < j_f; j++ )
            {
                ans[col*(i-(i_s-1))+(j-(j_s-1))] = A[n*i+j];
            }
        }
    }

    void MatrixCopy_(int m, int n, double *a, double *u)
    {
        int max;
        max = m * n ;
        for ( int i = 0 ; i < max ; i++ )
	        u[ i ] = a[ i ];
    }

    void MatrixDiagExpand(double* A, int row, int col, double* A_expand)
    {
        for(int i = 0; i < 2*row; i++)
	    {
	    	for(int j = 0; j < 2*col; j++)
	    	{
	    		A_expand[2*col*i + j] = 0;
	    	}
	    }
    
	    for(int i = 0; i < 2*row; i++)
	    {
	    	if(i < row)
	    	{
	    		for(int j = 0; j < col; j++)
	    		{
	    			A_expand[2*col*i + j] = A[col*i + j];
	    		}
	    	}
	    	else
	    	{
	    		for(int j = col; j < 2*col; j++)
	    		{
	    			A_expand[2*col*i + j] = A[col*(i-row) + (j-col)];
	    		}
	    	}
    
	    }
	
    }

    void cross(double* V, double* M)
    {
        M[0 * 3 + 0] = 0;     M[0 * 3 + 1] = -V[2];  M[0 * 3 + 2] = V[1];
	    M[1 * 3 + 0] = V[2];  M[1 * 3 + 1] = 0;      M[1 * 3 + 2] = -V[0];
	    M[2 * 3 + 0] = -V[1]; M[2 * 3 + 1] = V[0];   M[2 * 3 + 2] = 0;
	
    }


    void LUP_Descomposition(double* A, double* L, double* U, int* P, int n)
    {
        int row = 0;

        for(int i = 0; i < n; i++)
            P[i] = i;
    
        for(int i = 0; i < n-1; i++)
        {
            double p = 0.0;

            for(int j = i; j < n; j++)
            {
                if(abs(A[j * n + i]) > p)
                {
                    p = abs(A[ j * n + i]);
                    row = j;
                }
            }
            if(0 == p)
            {
                // cout<< "矩阵奇异，无法计算逆" <<endl;
                return ;
            }

            //交换P[i]和P[row]
            int tmp = P[i];
            P[i] = P[row];
            P[row] = tmp;

            double tmp2 = 0.0;
            for(int j = 0; j < n; j++)
            {
                //交换A[i][j]和 A[row][j]
                tmp2 = A[i * n + j];
                A[i * n + j] = A[row * n + j];
                A[row * n + j] = tmp2;
            }

            //以下同LU分解
            double u = A[i * n + i], l = 0.0;
            for(int j = i+1; j < n; j++)
            {
                l = A[j * n + i] / u;
                A[j * n + i] = l;
                for(int k = i + 1; k < n; k++)
                    A[j * n + k] = A[j *n + k] - A[i * n + k] * l;

            }

        }

        //构造L和U
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j <= i; j++)
            {
                if(i != j)
                {
                    L[i * n + j] = A[i * n + j];
                }
                else
                {
                    L[i * n + j] = 1;
                }
            }
            for(int k = i; k < n; k++)
            {
                U[i * n + k] = A[i * n + k];
            }
        }

    }

    //LUP求解方程
    void LUP_Solve(double* L, double* U, int* P, double* b, int n, double* inv_A_each)
    {
        // double* x = new double[n];
        double* y = new double[n];

        //正向替换
        for(int i = 0; i < n; i++)
        {
            y[i] = b[P[i]];

            for(int j = 0; j < i; j++)
            {
                y[i] = y[i] - L[i * n + j] * y[j];
            }
        }

        //反向替换
        for(int i = n-1; i >= 0; i--)
        {
            inv_A_each[i] = y[i];

            for(int j = n - 1; j > i; j--)
            {
                inv_A_each[i] = inv_A_each[i] - U[i * n + j] * inv_A_each[j];

            }

            inv_A_each[i] /= U[i * n + i];

        }

        delete[] y;

    }

    /*****************矩阵原地转置BEGIN********************/

    /* 后继 */
    int getNext(int i, int m, int n)
    {
      return (i%n)*m + i/n;
    }

    /* 前驱 */
    int getPre(int i, int m, int n)
    {
      return (i%m)*n + i/m;
    }

    /* 处理以下标i为起点的环 */
    void movedata(double *mtx, int i, int m, int n)
    {
      double temp = mtx[i]; // 暂存
      int cur = i;    // 当前下标
      int pre = getPre(cur, m, n);
      while(pre != i)
      {
        mtx[cur] = mtx[pre];
        cur = pre;
        pre = getPre(cur, m, n);
      }
      mtx[cur] = temp;
    }

    /* 转置，即循环处理所有环 */
    void transpose(double *mtx, int m, int n)
    {
      for(int i = 0; i < m * n; ++i)
      {
        int next = getNext(i, m, n);
        while (next > i)                 // 若存在后继小于i说明重复,就不进行下去了（只有不重复时进入while循环）
        {
            next = getNext(next, m, n);
        }
        if (next == i)
        {
            movedata(mtx, i, m, n);     // 处理当前环
        }    

      }
    }
    /*****************矩阵原地转置END********************/

    //LUP求逆(将每列b求出的各列x进行组装)
    void LUP_solve_inverse(double* A, int n, double* inv_A)
    {
        // 创建矩阵A的副本，注意不能直接用A计算，因为LUP分解算法已将其改变
        // double *A_mirror = new double[n*n]();
        // double A_mirror[n*n];
        double* A_mirror = new double[n * n];

        // double *inv_A=new double[n*n]();  // 最终的逆矩阵（还需要转置）
        // double inv_A_tempt[n*n];
        double* inv_A_tempt = new double[n * n];

        // double *inv_A_each=new double[n]();  // 矩阵逆的各列
        // double inv_A_each[n];
        double* inv_A_each = new double[n];

        // double *B    =new double[N*N]();
        // double *b    =new double[n]();  // b阵为B阵的列矩阵分量
        // double b[n];
        double* b = new double[n];  // b阵为B阵的列矩阵分量


        for(int i=0;i<n;i++)
        {
            double* L = new double[n*n];
            double* U = new double[n*n];
            int* P = new int[n];

            //构造单位阵的每一列
            for(int i=0;i<n;i++)
            {
                b[i]=0;
            }
            b[i] = 1;

            //每次都需要重新将A复制一份
            for(int i=0;i<n*n;i++)
            {
                A_mirror[i]=A[i];
            }

            LUP_Descomposition(A_mirror, L, U, P, n);
            LUP_Solve (L, U, P, b, n, inv_A_each);
            memcpy(inv_A_tempt + i * n, inv_A_each, n * sizeof(double));  //  将各列拼接起来

            delete[] L;
       	 	delete[] U;
        	delete[] P;

        }


        transpose(inv_A_tempt, n, n);  // 由于现在根据每列b算出的x按行存储，因此需转置
        for(int i=0;i<n;i++)
        	for(int j=0;j<n;j++)
        		inv_A[i * n + j] = inv_A_tempt[i * n + j];

        delete[] A_mirror;
        delete[] inv_A_tempt;
        delete[] inv_A_each;
        delete[] b;

    }


    double calc_determinantOfMatrix(double* a, int num, int row, int col, int n)
    {
        /*
        n 表示 要求行列式的矩阵
        a 是一个一维数组；
        num 是数组 a 的元素个数；
        row 是 a 要转换成的二维数组的 行数；
        col 是 a 要转换成的二维数组的 列数；
        该程序中，因为是求 矩阵的行列式 ， 所以 row = col = n, num = n**2
        */
        double num1 =1;
	    double num2 = 1;
	    double det = 1;
	    int index = 1;
	    double total = 1; // Initialize result
        double result_tempt = 0;
        double result = 0;

        double** mat_tempt = new double* [n];

        // temporary array for storing row
        double* temp = new double[n];    

        for(int k = 0; k < row; k++){
            mat_tempt[k] = new double[col];
            for(int j = 0; j < (num / row); j++)
                mat_tempt[k][j] = a[k * col + j];
        }

	    // loop for traversing the diagonal elements
	    for (int i = 0; i < n; i++){
	    	index = i; // initialize the index

	    	// finding the index which has non zero value
	    	while (index < n && mat_tempt[index][i] == 0){
	    		index++;
	    	}
	    	if (index == n) {// if there is non zero element
	    	
	    		// the determinant of matrix as zero
	    		continue;
	    	}
	    	if (index != i){
	    		// loop for swapping the diagonal element row and
	    		// index row
	    		for (int j = 0; j < n; j++)
	    		{
	    			std::swap(mat_tempt[index][j], mat_tempt[i][j]);
	    		}
	    		// determinant sign changes when we shift rows
	    		// go through determinant properties
	    		det = det * pow(-1, index - i);
	    	}

	    	// storing the values of diagonal row elements
	    	for (int j = 0; j < n; j++){
	    		temp[j] = mat_tempt[i][j];
	    	}
	    	// traversing every row below the diagonal element
	    	for (int j = i + 1; j < n; j++){
	    		num1 = temp[i]; // value of diagonal element
	    		num2 = mat_tempt[j][i]; // value of next row element

	    		// traversing every column of row
	    		// and multiplying to every row
	    		for (int k = 0; k < n; k++)
	    		{
	    			// multiplying to make the diagonal
	    			// element and next row element equal
	    			mat_tempt[j][k] = (num1 * mat_tempt[j][k]) - (num2 * temp[k]);

	    		}

	    		total = total * num1; // Det(kA)=kDet(A);

	    	}
	    }

	    // multiplying the diagonal elements to get determinant
	    for (int i = 0; i < n; i++)
	    	det = det * mat_tempt[i][i];

        result_tempt = det / total;
        result = result_tempt;

        for (int i = 0; i < n; i++)
            delete[] mat_tempt[i];

        delete[] mat_tempt;

        delete[] temp;
        return result;
    }

    void calc_Euler_dev2angvl(double* Euler_angle, double* Euler_dev, double* angvl)
    {
        /*
        N_Phi = np.array([
            [0, -sin(alpha_z), cos(alpha_z) * cos(beta_y)],
            [0,  cos(alpha_z), sin(alpha_z) * cos(beta_y)],
            [1,  0,            -sin(beta_y)]
        ])
        Euler_dev = np.array([
            [alpha_dot],
            [beta_dot],
            [gamma_dot]
        ])
        omega = N_Phi @ Euler_dev
        */

        double* N_phi = new double[3*3];
        double alpha_z = Euler_angle[0];
        double beta_y  = Euler_angle[1];
        double gamma_x = Euler_angle[2];

        N_phi[0] = 0; N_phi[1] = -sin(alpha_z); N_phi[2] = cos(alpha_z) * cos(beta_y);
        N_phi[3] = 0; N_phi[4] = cos(alpha_z);  N_phi[5] = sin(alpha_z) * cos(beta_y);
        N_phi[6] = 1; N_phi[7] = 0;             N_phi[8] = -sin(beta_y);
        MatrixMulti_( 3,  3,  1,  N_phi,  Euler_dev,  angvl);

        delete[] N_phi;

    }


    double rand_range(const double& min_val, const double& max_val){
        std::default_random_engine rng(std::chrono::system_clock::now().time_since_epoch().count());
        std::uniform_real_distribution<double> distr(min_val, max_val);
        return distr(rng);
    }

    std::vector<double> rand_range_vector(const std::vector<double>& min_val, const std::vector<double>& max_val){
        std::vector<double> result;
        for (std::vector<double>::size_type i = 0; i < min_val.size(); ++i) {
            result.push_back(rand_range(min_val[i], max_val[i]));
        }
        return result;
    }

    double rand_normal(const double& mean_val, const double& std_val){
        std::default_random_engine rng(std::chrono::system_clock::now().time_since_epoch().count());
        std::normal_distribution<double> distr(mean_val, std_val);
        return distr(rng);
    }

    double rand_levy(const double& beta){
        // double beta = 1.5;
        double alpha_u = pow(tgamma(1 + beta) * sin(M_PI * beta / 2) / (tgamma((1 + beta) / 2) * beta * pow(2, (beta - 1) / 2)), 1 / beta);
        double alpha_v = 1;
        double u = rand_normal(0, alpha_u);
        double v = rand_normal(0, alpha_v);
        double z = u / pow(fabs(v), 1 / beta);
        return z;

    }

    double sign_func(const double& x)
    {
        if (x > 0) {
            return 1;
        } else if (x < 0) {
            return -1;
        } else {
            return 0;
        }
    }

    std::vector<double> vector_multiply(const double& scalar, const std::vector<double>& vector_1)
    {
        std::vector<double> result;
        for (std::vector<double>::size_type i = 0; i < vector_1.size(); ++i) {
            result.push_back(scalar * vector_1[i]);
        }
        return result;
    }

    std::vector<double> vector_add(const std::vector<double>& vector_1, const std::vector<double>& vector_2){
        std::vector<double> result;
        for (std::vector<double>::size_type i = 0; i < vector_1.size(); ++i) {
            result.push_back(vector_1[i] + vector_2[i]);
        }
        return result;
    }

    std::vector<double> vector_subtract(const std::vector<double>& vector_1, const std::vector<double>& vector_2)
    {
        std::vector<double> result;
        for (std::vector<double>::size_type i = 0; i < vector_1.size(); ++i) {
            result.push_back(vector_1[i] - vector_2[i]);
        }
        return result;
    }

    std::vector<double> vector_dot_vector(const std::vector<double>& vector1, const std::vector<double>& vector2)
    {
        std::vector<double> result;
        for (std::vector<double>::size_type i = 0; i < vector1.size(); ++i) {
            result.push_back(vector1[i] * vector2[i]);
        }
        return result;
    }

    std::vector<double> vector_min(const std::vector<double>& vector1, const std::vector<double>& vector2)
    {
        std::vector<double> result;
        for (std::vector<double>::size_type i = 0; i < vector1.size(); ++i) {
            result.push_back(std::min(vector1[i], vector2[i]));
        }
        return result;
    }

    std::vector<double> vector_max(const std::vector<double>& vector1, const std::vector<double>& vector2)
    {
        std::vector<double> result;
        for (std::vector<double>::size_type i = 0; i < vector1.size(); ++i) {
            result.push_back(std::max(vector1[i], vector2[i]));
        }
        return result;
    }

    std::vector<double> vector_fabs(const std::vector<double>& vector1)
    {
        std::vector<double> result;
        for (std::vector<double>::size_type i = 0; i < vector1.size(); ++i) {
            result.push_back(std::fabs(vector1[i]));
        }
        return result;
    }

    std::vector<double> vector_pow(const std::vector<double>& vector1, const double& power)
    {
        std::vector<double> result;
        for (std::vector<double>::size_type i = 0; i < vector1.size(); ++i) {
            result.push_back(std::pow(vector1[i], power));
        }
        return result;
    }

}



#endif 

