#ifndef PSO_WOA_H
#define PSO_WOA_H

#include "kinetics.h"

using namespace utils;
using namespace kinetics;
using namespace openGJK;


double RPY_END_DESIRED[] = {-30 * M_PI / 180, -50 * M_PI /180, 45 * M_PI /180 };
double Pe_DESIRED[] = {1.5, -1.4, 2.2};
int NUM_CALCULATIONS = 40;
int NUM_PARTICLES = 96;
int NUM_DIMENSIONS = 7;
int MAX_ITER = 2000;
// set the range of the parameters
std::vector<double> maxSpeed = {1.5,   0.628,  1.5,  1.5,  1.0,  0.8,  1.5};
std::vector<double> minSpeed = {-1.5, -0.628, -1.5, -1.5, -1.5, -0.8, -1.5};
std::vector<double> minPosition = {-6.28319,  0,       -6.28319, -6.28319, -6.28319, -2.0944, -6.28319};
std::vector<double> maxPosition = { 6.28319,  3.14159,  6.28319,  6.28319,  6.28319,  2.0944,  6.28319};
// set _rand_0 to 0.0, set _rand_1 to 1.0
std::vector<double> _vector_0(NUM_DIMENSIONS, 0);
std::vector<double> _vector_1(NUM_DIMENSIONS, 1);
// Define the PSO parameters
double inertGuideCoe = 0.9;  // 0.7;    // inertia weight
double localGuideCoe = 1.47;  //  1.2;   // cognitive weight
double globalGuideCoe = 1.47;  //  1.2;   // social weight


// define a class to store the position and fitness of Particle_pso and Particle_woa 
class Particle
{
public:
    std::vector<double> position;
    double fitness;

    Particle() {}

};

// define a function to compare the fitness of two particles 
bool particleComparator(const Particle& p1, const Particle& p2) 
{
    return p1.fitness < p2.fitness;
}



double calc_fitness_woa(const std::vector<double>& position) 
{
    // define the parameters
    double* para = new double[N];
	double* eta_end = new double;
	double* xi_end = new double[3];
	double* Pe = new double[3];
	double* eta_b = new double;
	double* xi_b = new double[3];
	double* quaternion_end_desired = new double[4];
	double* eta_end_desired = new double;
	double* xi_end_desired = new double[3];
	double* delta_eta_base = new double;
	double* delta_eta_end = new double;
	double* delta_xi_base = new double[3];
	double* delta_xi_end = new double[3];
	double* delta_Pe_end = new double[3];
    double* p_e_initial = new double[3];
    double* locus = new double;
    double* delta_xi_b_distrb_max = new double;
    double* manipl = new double;
    double* T_min = new double;
	double* collision_times = new double;

	for (int i = 0; i < N; i++)
		para[i] = position[i];

    double K_a = 1 / 0.0002;  
	double K_p = 1 / 0.002;  
    double K_s = 0;  
    double K_b = 1 / 0.0008;
    double K_M = 1;
    double K_t = 1 / 100;  
	double cost_func = 0;

	double delta_xi_end_mod_temp = 0;
	double delta_Pe_end_mod_temp = 0;
	double delta_xi_base_mod_temp = 0;

	forward_kin( para,  eta_end,  xi_end,  Pe,  eta_b,  xi_b,  p_e_initial,  locus,  delta_xi_b_distrb_max,  manipl, T_min, collision_times);

// ********************************************************  straight_line_locus  ************************************************************************************
    double delta_p_e = 0;
    double straight_line_locus = 0;
    for(int i = 0; i < 3; i++)
    {
        delta_p_e = p_e_initial[i] - Pe[i];
        straight_line_locus += delta_p_e * delta_p_e;
    }
    straight_line_locus = sqrt(straight_line_locus);
// **********************************************************************************************************************************************************************

	zyx2quaternion(RPY_END_DESIRED[2], RPY_END_DESIRED[1], RPY_END_DESIRED[0], quaternion_end_desired);

	*eta_end_desired = quaternion_end_desired[0];

	for (int i = 0; i < 3; i++)
	{
		xi_end_desired[i] = quaternion_end_desired[i + 1];
	}

	delta_var(Pe_DESIRED, eta_end_desired, xi_end_desired, Pe, *eta_end, xi_end, *eta_b, xi_b, delta_eta_end, delta_xi_end, delta_Pe_end, delta_eta_base, delta_xi_base);

	for (int i = 0; i < 3; i++)
	{
		delta_xi_end_mod_temp += delta_xi_end[i] * delta_xi_end[i];
		delta_Pe_end_mod_temp += delta_Pe_end[i] * delta_Pe_end[i];
		delta_xi_base_mod_temp += delta_xi_base[i] * delta_xi_base[i];
	}

	delta_xi_end_mod_temp = sqrt(delta_xi_end_mod_temp); 

    delta_Pe_end_mod_temp = sqrt(delta_Pe_end_mod_temp); 

    delta_xi_base_mod_temp = sqrt(delta_xi_base_mod_temp);


	// ****************************************  RPY error of end-effector,  Pe error  ********************************************************************************
    
	// cout << "manipl:" << "  " << (*manipl) << endl;

    cost_func = K_a * delta_xi_end_mod_temp 
                + K_p * delta_Pe_end_mod_temp 
                + K_b * (*delta_xi_b_distrb_max) 
                + K_s * fabs((*locus) - straight_line_locus)  
                + K_M * (1 / (*manipl))
                + K_t * (*T_min)
				+ (*collision_times);


    delete[] para;
    delete eta_end;
    delete[] xi_end ;
    delete[] Pe ;
    delete eta_b ;
	delete[] xi_b ;
	delete[] quaternion_end_desired;
	delete eta_end_desired;
	delete[] xi_end_desired ;
	delete delta_eta_base ;
	delete delta_eta_end ;
	delete[] delta_xi_base ;
	delete[] delta_xi_end ;
	delete[] delta_Pe_end ; 
	delete[] p_e_initial ;
	delete locus ;
	delete delta_xi_b_distrb_max;
	delete manipl;
	delete T_min;
	delete collision_times;


	return cost_func;

}



// define the fitness function
double calc_fitness_pso(const std::vector<double>& position) 
{
    // define the parameters
    double* para = new double[N];

	double* eta_end = new double;
	double* xi_end = new double[3];
	double* eta_end_desired = new double;
	double* xi_end_desired = new double[3];
	double* quaternion_end_desired = new double[4];
	double* delta_eta_end = new double;
	double* delta_xi_end = new double[3];
	double* delta_Pe_end = new double[3];
	double* p_e_initial = new double[3];
	double* Pe = new double[3];

	double* eta_b = new double;
	double* xi_b = new double[3];
	double* delta_eta_base = new double;
	double* delta_xi_base = new double[3];
	
    double* locus = new double;
    double* delta_xi_b_distrb_max = new double;
    double* manipl = new double;
    double* T_min = new double;
	double* collision_times = new double;

	for (int i = 0; i < N; i++)
		para[i] = position[i];

    double K_a = 1 / 0.0002;  
	double K_p = 1 / 0.002;  
    double K_s = 0;  
    double K_b = 1 / 0.0008;
    double K_M = 1;
    double K_t = 1 / 100;  
	double cost_func = 0;

	double delta_xi_end_mod_temp = 0;
	double delta_Pe_end_mod_temp = 0;
	double delta_xi_base_mod_temp = 0;

	// forward_kin_3( para,  eta_end,  xi_end,  Pe,  eta_b,  xi_b,  p_e_initial,  locus,  delta_xi_b_distrb_max,  manipl, T_min, collision_times);
	forward_kin( para,  eta_end,  xi_end,  Pe,  eta_b,  xi_b,  p_e_initial,  locus,  delta_xi_b_distrb_max,  manipl, T_min, collision_times);

    // *************************  straight_line_locus  **************************************************
    double delta_p_e = 0;

    double straight_line_locus = 0;

    for(int i = 0; i < 3; i++)
    {
        delta_p_e = p_e_initial[i] - Pe[i];

        straight_line_locus += delta_p_e * delta_p_e;

    }

    straight_line_locus = sqrt(straight_line_locus);

    // **************************************************************************************************

	zyx2quaternion(RPY_END_DESIRED[2], RPY_END_DESIRED[1], RPY_END_DESIRED[0], quaternion_end_desired);

	*eta_end_desired = quaternion_end_desired[0];

	for (int i = 0; i < 3; i++)
	{
		xi_end_desired[i] = quaternion_end_desired[i + 1];
	}

	delta_var(Pe_DESIRED, eta_end_desired, xi_end_desired, Pe, *eta_end, xi_end, 
				*eta_b, xi_b, delta_eta_end, delta_xi_end, delta_Pe_end, delta_eta_base, delta_xi_base);

	for (int i = 0; i < 3; i++)
	{
		delta_xi_end_mod_temp += delta_xi_end[i] * delta_xi_end[i];

		delta_Pe_end_mod_temp += delta_Pe_end[i] * delta_Pe_end[i];

		delta_xi_base_mod_temp += delta_xi_base[i] * delta_xi_base[i];
		
	}

	delta_xi_end_mod_temp = sqrt(delta_xi_end_mod_temp); 

    delta_Pe_end_mod_temp = sqrt(delta_Pe_end_mod_temp); 

    delta_xi_base_mod_temp = sqrt(delta_xi_base_mod_temp);


	// ****************************************  RPY error of end-effector,  Pe error  ********************************************************************************
    
	// cout << "manipl:" << "  " << (*manipl) << endl;

    cost_func = K_a * delta_xi_end_mod_temp + K_p * delta_Pe_end_mod_temp + K_b * (*delta_xi_b_distrb_max) + 
				K_s * fabs((*locus) - straight_line_locus)  + 
				K_M * (1 / (*manipl)) + 
				K_t * (*T_min) + (*collision_times);


    delete[] para;
    delete eta_end;
    delete[] xi_end ;
    delete[] Pe ;
    delete eta_b ;
	delete[] xi_b ;
	delete[] quaternion_end_desired;
	delete eta_end_desired;
	delete[] xi_end_desired ;
	delete delta_eta_base ;
	delete delta_eta_end ;
	delete[] delta_xi_base ;
	delete[] delta_xi_end ;
	delete[] delta_Pe_end ; 
	delete[] p_e_initial ;
	delete locus ;
	delete delta_xi_b_distrb_max;
	delete manipl;
	delete T_min;
	delete collision_times;


	return cost_func;

}




// Define a class to represent a particle
class Particle_pso {
public:
    std::vector<double> position;
    std::vector<double> velocity;
    std::vector<double> best_position;
    std::vector<double> r1;
    std::vector<double> r2;
    double fitness;
    double best_fitness;

    Particle_pso() 
    {
        position.resize(NUM_DIMENSIONS);
        velocity.resize(NUM_DIMENSIONS);
        best_position.resize(NUM_DIMENSIONS);
        r1.resize(NUM_DIMENSIONS);
        r2.resize(NUM_DIMENSIONS);

        velocity = rand_range_vector(minSpeed, maxSpeed);
        position = rand_range_vector(minPosition, maxPosition);

        fitness = calc_fitness_pso(position);
        best_fitness = fitness;
        best_position = position;
    }

   
    void update_fitness() 
    {
        fitness = calc_fitness_pso(position);
        if (fitness < best_fitness) 
        {
            best_fitness = fitness;
            best_position = position;
        }
    }
    
    void update_velocity(std::vector<double>& global_best_position) 
    {
        r1 = rand_range_vector(_vector_0, _vector_1);
        r2 = rand_range_vector(_vector_0, _vector_1);

        velocity = vector_add(
                        vector_add(
                            vector_multiply(inertGuideCoe, velocity), 
                            vector_multiply(localGuideCoe, vector_dot_vector(r1, vector_subtract(best_position, position)))
                        ), 
                        vector_multiply(globalGuideCoe, vector_dot_vector(r2, vector_subtract(global_best_position, position)))
                    );
    
        velocity = vector_min(vector_max(velocity, minSpeed), maxSpeed);
    }
    
    void update_position() 
    {
        position = vector_add(position, velocity);
        position = vector_min(vector_max(position, minPosition), maxPosition);
    }
};


// Define a class to represent a particle
class Particle_woa {
public:
    std::vector<double> position;
    std::vector<double> r1;
    std::vector<double> r2;
    std::vector<double> r3;
    double fitness;
    double beta;
    double b;

    Particle_woa() {
        position.resize(NUM_DIMENSIONS);
        r1.resize(NUM_DIMENSIONS);
        r2.resize(NUM_DIMENSIONS);
        r3.resize(NUM_DIMENSIONS);
        position = rand_range_vector(minPosition, maxPosition);
        fitness = calc_fitness_woa(position);
        beta = 1.5;
        b = 1;
    }

    // Define a function to update the fitness
    void update_fitness() 
    {
        fitness = calc_fitness_woa(position);
    }

    void update_position( double a,  double a2,  Particle_woa& rand_particle_woa,  Particle_woa& best_particle_woa) 
    {
        r1 = rand_range_vector(_vector_0, _vector_1);
        r2 = rand_range_vector(_vector_0, _vector_1);
        r3 = rand_range_vector(_vector_0, _vector_1);
        std::vector<double> A = vector_multiply(a, vector_subtract(vector_multiply(2, r1), _vector_1));
        std::vector<double> C = vector_multiply(2, r2);
        std::vector<double> l = vector_add(vector_multiply(a2 - 1, r3), _vector_1);
        double p = rand_range(0, 1);

        for (int i = 0; i < NUM_DIMENSIONS; i++) 
        {
            if (p < 0.5) 
            {
                if (std::fabs(A[i]) > 1) 
                {
                    position[i] = rand_particle_woa.position[i] - A[i] * (std::fabs(C[i] * rand_particle_woa.position[i] - position[i])) + 
                                    rand_range(0, 1) * sign_func(rand_range(0, 1) - 0.5) * rand_levy(beta);
                } 
                else 
                {
                    position[i] = best_particle_woa.position[i] - A[i] * (std::fabs(C[i] * best_particle_woa.position[i] - position[i])) + 
                                    rand_range(0, 1) * sign_func(rand_range(0, 1) - 0.5) * rand_levy(beta);
                }
            }
            else 
            {
                position[i] = std::fabs(best_particle_woa.position[i] - position[i]) * exp(b * l[i]) * cos(2 * M_PI * l[i]) + 
                                best_particle_woa.position[i] + rand_range(0, 1) * sign_func(rand_range(0, 1) - 0.5) * rand_levy(beta);
            }
            // clamp the position vector to the range [minPosition, maxPosition]
            // if the position is greater than maxPosition, set it to (position[j] + maxPosition[j]) / 2
            // if the position is less than minPosition, set it to (position[j] + minPosition[j]) / 2
            if (position[i] > maxPosition[i]) 
            {
                position[i] = (position[i] + maxPosition[i]) / 2;
            }
            if (position[i] < minPosition[i]) 
            {
                position[i] = (position[i] + minPosition[i]) / 2;
            }
        }
    }
};






#endif // PSO_WOA_H
