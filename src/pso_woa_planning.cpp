#include <thread>
#include <algorithm>
#include <numeric>
#include <limits>  
#include "pso_woa.h"


#define WMAX 0.9
#define WMIN 0.1


void evaluate_fitness(Particle_pso& particle_pso, Particle_woa& particle_woa);

void update_global_best(std::vector<Particle_pso>& particles_pso, std::vector<double>& pso_global_best_position, double& pso_global_best_fitness,
                        std::vector<Particle_woa>& particles_woa, std::vector<double>& woa_global_best_position, double& woa_global_best_fitness, 
                        Particle_woa& best_particle_woa);

void run_pso_woa(std::vector<Particle_pso>& particles_pso, std::vector<Particle_woa>& particles_woa, std::vector<Particle>& particles, int num_threads);

void reset_GuideCoe();


int main() 
{
    // set the start time of the program
    auto start = std::chrono::high_resolution_clock::now();

    // int num_threads = std::thread::hardware_concurrency();  
    // Get the number of threads supported by the hardware
    int num_threads = std::thread::hardware_concurrency();

    // Print the parameters
    std::cout << "NUM_CALCULATIONS:  " << NUM_CALCULATIONS << std::endl;
    std::cout << "NUM_PARTICLES:  " << NUM_PARTICLES << std::endl;
    std::cout << "NUM_DIMENSIONS:  " << NUM_DIMENSIONS << std::endl;
    std::cout << "MAX_ITER:  " << MAX_ITER << std::endl;
    std::cout << "num_threads:  " << num_threads << std::endl;

    for(int k = 0; k < NUM_CALCULATIONS; ++k)
    {
        reset_GuideCoe();
        std::cout << "The current calculation loop is: " << k << std::endl;  

        // Create the particles 
        std::vector<Particle_pso> particles_pso(NUM_PARTICLES);  
        std::vector<Particle_woa> particles_woa(NUM_PARTICLES);
        std::vector<Particle> particles(NUM_PARTICLES * 2);

        // Run the PSO algorithm using 12 threads
        run_pso_woa(particles_pso, particles_woa, particles, num_threads);

        particles_pso.clear();
        particles_woa.clear();
        particles.clear();

    }

    // set the end time of the program
    auto end = std::chrono::high_resolution_clock::now();

    // print the total time taken by the program
    std::cout << "Total time taken by the program: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;

    return 0;

}


// Define a function to run the PSO algorithm
void run_pso_woa(std::vector<Particle_pso>& particles_pso, std::vector<Particle_woa>& particles_woa, std::vector<Particle>& particles, int num_threads) 
{
    // define parameters for woa
    double a;
    double a2;
    int rand_particle_woa; 

    Particle_woa best_particle_woa;  // Create a particle to store the best particle for woa

    int num_particles_per_thread = NUM_PARTICLES / num_threads;  // Calculate the number of particles per thread

    std::vector<double> pso_global_best_position(NUM_DIMENSIONS); // Create a vector to store the global best position for pso and woa
    std::vector<double> woa_global_best_position(NUM_DIMENSIONS);  

    double pso_global_best_fitness = std::numeric_limits<double>::max();  // initialize the global best fitness to the maximum value
    double woa_global_best_fitness = std::numeric_limits<double>::max();
    
    std::vector<std::thread> threads;  // Create the threads vector

    for (int i = 0; i < num_threads; ++i)  // Create the threads and evaluate the fitness of each particle
    {
        threads.emplace_back([&particles_pso, &particles_woa, i, num_particles_per_thread]() {
            for ( int j = i * num_particles_per_thread; j < (i + 1) * num_particles_per_thread; j++ ) 
            {
                evaluate_fitness(particles_pso[j], particles_woa[j]);
            }
        });
    }

    for (auto& thread : threads)  // Wait for the threads to finish
    {
        thread.join();
    }

    for(int i = 0; i < NUM_PARTICLES; ++i)  // push the position and fitness of particles_pso and particles_woa into the particles
    {
        particles[i].position = particles_pso[i].position;
        particles[i].fitness = particles_pso[i].fitness;
        particles[i + NUM_PARTICLES].position = particles_woa[i].position;
        particles[i + NUM_PARTICLES].fitness = particles_woa[i].fitness;
    }

    std::sort(particles.begin(), particles.end(), particleComparator);  // Sort the particles in ascending order of fitness

    for(int i = 0; i < NUM_PARTICLES; ++i)
    {
        particles_pso[i].position = particles[i].position;
        particles_pso[i].fitness = particles[i].fitness;
        particles_woa[i].position = particles[i + NUM_PARTICLES].position;
        particles_woa[i].fitness = particles[i + NUM_PARTICLES].fitness;
    }

    // Update the global best position
    update_global_best(particles_pso, pso_global_best_position, pso_global_best_fitness, particles_woa, woa_global_best_position, woa_global_best_fitness, best_particle_woa);

    // Run the PSO_WOA algorithm
    for (int iter = 0; iter < MAX_ITER; ++iter) 
    {
        // update the inertial coefficient, please do not update the global parameters in the threads, otherwise, it will cause a mess on the global parameters.
        // Global objects, or pointers to objects passed into a thread’s initial “main” function, can be accessed by multiple threads. 
        // However, to have a fully functioning, bug-free program, you need to ensure threads share resources safely. 
        // If one thread reads a resource, like an array or chunk of memory, while another one is writing to it, that’s going to result in corrupted data. 
        inertGuideCoe = WMAX - (WMAX - WMIN) * (2* double(iter) / double(MAX_ITER) - pow(double(iter) / double(MAX_ITER), 2));  

        a = 2 * (1 - (double)iter / (double)MAX_ITER);  
        
		a2 = (-1) + (double)iter * ((-1) / (double)MAX_ITER);  

        rand_particle_woa = rand() % ((NUM_PARTICLES - 1) - 0 + 1) + 0;  

        threads.clear();  

        for (int i = 0; i < num_threads; ++i)
        {  
            threads.emplace_back([&particles_pso, &pso_global_best_position, &particles_woa, &woa_global_best_position, 
                                    a, a2, &best_particle_woa, rand_particle_woa, i, num_particles_per_thread]() {  

                for (int j = i * num_particles_per_thread; j < (i + 1) * num_particles_per_thread; ++j) 
                {
                    particles_pso[j].update_velocity(pso_global_best_position);
                    particles_pso[j].update_position();

                    auto get_other_particle = [j, rand_particle_woa]() -> int 
                    {
                        // Make a copy of rand_particle_woa
                        int temp_rand_particle_woa = rand_particle_woa; 
                        while (temp_rand_particle_woa == j) 
                        {
                            temp_rand_particle_woa = rand() % ((NUM_PARTICLES - 1) - 0 + 1) + 0;
                        }
                        return temp_rand_particle_woa;
                    };
                    particles_woa[j].update_position(a, a2, particles_woa[get_other_particle()], best_particle_woa);
                }
            });
        }

        // wait for the threads to finish
        for (auto& thread : threads) 
        {
            thread.join();  
        }

        threads.clear();

        // Evaluate the fitness of each particle
        for (int i = 0; i < num_threads; ++i) 
        {
            threads.emplace_back([&particles_pso, &particles_woa, i, num_particles_per_thread]() {
                for ( int j = i * num_particles_per_thread; j < (i + 1) * num_particles_per_thread; j++ ) 
                {
                    evaluate_fitness(particles_pso[j], particles_woa[j]);
                }
            });
        }

        // Wait for the threads to finish
        for (auto& thread : threads) 
        {
            thread.join();
        }

        for(int i = 0; i < NUM_PARTICLES; ++i)
        {
            particles[i].position = particles_pso[i].position;
            particles[i].fitness = particles_pso[i].fitness;
            particles[i + NUM_PARTICLES].position = particles_woa[i].position;
            particles[i + NUM_PARTICLES].fitness = particles_woa[i].fitness;
        }

        std::sort(particles.begin(), particles.end(), particleComparator);

        for(int i = 0; i < NUM_PARTICLES; ++i)
        {
            particles_pso[i].position = particles[i].position;
            particles_pso[i].fitness = particles[i].fitness;
            particles_woa[i].position = particles[i + NUM_PARTICLES].position;
            particles_woa[i].fitness = particles[i + NUM_PARTICLES].fitness;
        }

        // Update the global best position
        update_global_best(particles_pso, pso_global_best_position, pso_global_best_fitness, 
                            particles_woa, woa_global_best_position, woa_global_best_fitness, best_particle_woa);

        // Print the results every 50 iterations 
        if(iter % 50 == 0)
        {
            std::cout << "inertGuideCoe is :" << inertGuideCoe << std::endl;
            std::cout << "a is :" << a << std::endl;
            std::cout << "a2 is :" << a2 << std::endl;
			std::cout << "iteration is :" << iter << std::endl;
			std::cout << "particles[0].fitness is :" << particles[0].fitness << std::endl;

			for (auto x : particles[0].position) 
            {
                std::cout << x << " ";
            }	
        			
            std::cout << std::endl;

		}
		
        // Check if the global best fitness is less than 1
		if (particles[0].fitness < 1)
		{
			std::cout << "Find solution, iteration is:" << iter << std::endl;
            std::cout << "global_best_fitness is :" << particles[0].fitness << std::endl;

			for (auto x : particles[0].position) {
                std::cout << x << " ";
            }	
        			
            std::cout << std::endl;
			break;
		}

    }

    // Clear the threads vector 
    threads.clear();

    // Print the results
    std::cout << "Global best fitness: " << particles[0].fitness << std::endl;
    std::cout << "Global best position: ";
    for (auto x : particles[0].position) {
        std::cout << x << " ";
    }

    std::cout << std::endl;

}


// Define a function to evaluate the fitness of a particle
void evaluate_fitness(Particle_pso& particle_pso, Particle_woa& particle_woa) 
{
    particle_pso.update_fitness();
    particle_woa.update_fitness();
}


// Define a function to update the global best position
void update_global_best(std::vector<Particle_pso>& particles_pso, std::vector<double>& pso_global_best_position, double& pso_global_best_fitness,
                        std::vector<Particle_woa>& particles_woa, std::vector<double>& woa_global_best_position, double& woa_global_best_fitness, 
                        Particle_woa& best_particle_woa) 
{
    for (auto& particle_pso : particles_pso) 
    {
        if (particle_pso.best_fitness < pso_global_best_fitness) 
        {  
            // minimize the fitness
            pso_global_best_fitness = particle_pso.best_fitness;
            pso_global_best_position = particle_pso.best_position;
        }
    }

    for (auto& particle_woa : particles_woa) 
    {
        if (particle_woa.fitness < woa_global_best_fitness) 
        {  
            // minimize the fitness
            woa_global_best_fitness = particle_woa.fitness;
            woa_global_best_position = particle_woa.position;
            best_particle_woa = particle_woa;
        }
    }
}


void reset_GuideCoe()
{
    inertGuideCoe = 0.9;  
    localGuideCoe = 1.47;  
    globalGuideCoe = 1.47; 
}

