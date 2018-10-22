/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	//calculate needed number of particles; Hypothesis : the accuracy should be between 0  to 10 cm
	// so I will divide the GPS std deviation in x to sample every 10 cm multiplied by 2

	// convert from meters to cm
	int x_cm =  x * 100;

	// i.e no. of samples = 2(x_cm / 10)
	this->num_particles = 2* (x_cm / 10);

	// default random generator
	default_random_engine gen;

	// create normal Gussian distribution for GPS x pos
	normal_distribution<double> dist_x(x, std[0]);

	// Create normal distributions for y  .
    normal_distribution<double> dist_y(y, std[1]);

    // Create normal distributions for theta
    normal_distribution<double> dist_theta(theta, std[2]);

    // create and initialize particles

    for(int i = 0 ; i < this->num_particles ; i++)
    {
    	Particle particle;

    	particle.x = dist_x(gen);
    	particle.y = dist_y(gen);
    	particle.theta = dist_theta(gen);
    	particle.weight = 1.0;
    	particle.id = i;
    	this->particles.push_back(particle);

    }

    this->is_initialized = true;


}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	//Gussian distribution for the movement step
	// default random generator
	default_random_engine gen;

	// create normal Gussian distribution for GPS x pos
	normal_distribution<double> dist_x(0, std_pos[0]);

	// Create normal distributions for y  .
    normal_distribution<double> dist_y(0, std_pos[1]);

    // Create normal distributions for theta
    normal_distribution<double> dist_theta(0, std_pos[2]);

	// loop on all particles and apply the movement
	for(int i = 0 ; i < this->num_particles ; i++)
	{
		// check if the yaw_rate !=0
		if(fabs(yaw_rate) > 0.00001 )
		{
			// use the  motion model equation with non-zero yaw_rate

			this->particles[i].x = this->particles[i].x +( (velocity / yaw_rate) *
					(sin(this->particles[i].theta + (yaw_rate * delta_t )) - sin(this->particles[i].theta))) + dist_x(gen);

			this->particles[i].y = this->particles[i].y + ((velocity / yaw_rate) *
					(cos(this->particles[i].theta) - cos(this->particles[i].theta + (yaw_rate * delta_t)))) + dist_y(gen);

			this->particles[i].theta = this->particles[i].theta + (yaw_rate *delta_t) + dist_theta(gen);

		}
		else
		{

			// motion model equation if the yaw rate = 0

			this->particles[i].x = this->particles[i].x +
					( (velocity * delta_t) *  cos(this->particles[i].theta) ) + dist_x(gen);

			this->particles[i].y = this->particles[i].y +
					( (velocity * delta_t) *  sin(this->particles[i].theta) ) + dist_y(gen);

			this->particles[i].theta = this->particles[i].theta + dist_theta(gen);

		}


	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	//loop on all observations to associate it with one particular map landmark then update this
	// observation with the ID of this landmark
	for (int i = 0 ; i < observations.size() ; i++)
	{
		double temp_min_dist =  std::numeric_limits<double>::max();

		// loop on all landmarks to get the closest one

		for(int j = 0 ; j < predicted.size() ; j++)
		{
			//calculate the distance between two points
			double temp_dist = dist(predicted[j].x , predicted[j].y , observations[i].x , observations[i].y );

			// if the temp minimum distance value is bigger than the returned distance
			// then update the value of temp minimum distance with this value and update the observation ID.
			if(temp_min_dist > temp_dist)
			{
				temp_min_dist = temp_dist;

				observations[i].id = predicted[j].id ;

			}

		}

	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_nomal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html



	// temp vector for holding one particle transformed observations till it pushed back
	std::vector<LandmarkObs> tobservations;

	//temp vector for holding the map landmarks that is in the particle sensor range
	std::vector<LandmarkObs> inrange_map_landmarks;

	LandmarkObs tobservation;

	std::vector<int> association_ids;
	std::vector<double> sense_x;
	std::vector<double> sense_y;

	// loop on all particles to transform the observations from the vehicle coordinates to particle
	//Map coordinates , create vector with the interested map landmarks which are in sensor region
	// then calculate the particle weight and assign its associated landmark IDs and their x and y.

	for (int i = 0 ; i < this->num_particles ; i++)
	{

		// loop on observations and transform them from vehicle coordinates to particle map coordinates
		for(int j = 0 ; j < observations.size() ; j++)
		{

			tobservation.x = particles[i].x +
					(cos(particles[i].theta) * observations[j].x) -
					(sin(particles[i].theta) * observations[j].y);

			tobservation.y = particles[i].y +
					(sin(particles[i].theta) * observations[j].x) +
					(cos(particles[i].theta) * observations[j].y);


			//push the transformed observation to the observations vector
			tobservations.push_back(tobservation);

		}

		// create a vector with the landmarks in the sensor range according to the particle pos.

		for(int j = 0 ; j < map_landmarks.landmark_list.size() ; j++)
		{
			tobservation.x = map_landmarks.landmark_list[j].x_f;
			tobservation.y = map_landmarks.landmark_list[j].y_f;
			tobservation.id = map_landmarks.landmark_list[j].id_i;

			//check if this landmark in particle sensor range and if yes push this landmark into the
			// in range vector.

			if(dist(tobservation.x , tobservation.y , particles[i].x , particles[i].y) <= sensor_range )
			{
				inrange_map_landmarks.push_back(tobservation);

			}


		}

		// Associate the tobservations with the inrange_landmarks
		dataAssociation(inrange_map_landmarks,tobservations);

		// reset weight for particle before calc.
		 particles[i].weight = 1.0;

		// loop on all observations for this particle and update the weight of it.
		for(int j = 0 ; j < tobservations.size() ; j++)
		{

			double predicted_x, predicted_y, observed_x, observed_y, std_x, std_y;
			int id;

			observed_x = tobservations[j].x;
			observed_y = tobservations[j].y;
			id = tobservations[j].id;

			std_x = std_landmark[0];
			std_y = std_landmark[1];

			predicted_x = (double)map_landmarks.landmark_list[id-1].x_f;
			predicted_y = (double)map_landmarks.landmark_list[id-1].y_f;


			//Update the particle weights by multiplying the gaussian multivariate distributions
			particles[i].weight *= 1.0/(2* M_PI * std_x * std_y) *
					exp(-((pow(observed_x - predicted_x, 2) / (2 * pow(std_x,2))) +
							(pow(observed_y - predicted_y, 2) / (2 * pow(std_y,2)))));


			association_ids.push_back(tobservations[j].id);
			sense_x.push_back(tobservations[j].x);
			sense_y.push_back(tobservations[j].y);

		}

		//set the particle association ID and the x_sens and y_sense vectors

		SetAssociations(particles[i],association_ids,sense_x,sense_y);



		// clear the tobservations and inrange_map_landmarks vectors for the next particle
		//calculations.

		tobservations.clear();
		inrange_map_landmarks.clear();
		association_ids.clear();
		sense_x.clear();
		sense_y.clear();

	}


}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	default_random_engine gen;

	//This vector contains the rsampled particles
	std::vector<Particle> resampled_particles;

	//choose random index
	uniform_int_distribution<int> int_dist(0,this->num_particles-1);

	int index = int_dist(gen);

	double beta = 0.0;

	std::vector<double> weights;

	//collect all particles weights in the weights vector
	for (int i=0; i<this->num_particles; i++)
	{
		weights.push_back(particles[i].weight);
	}

	// extract the maximum weight
	double max_weight = *max_element(weights.begin(), weights.end());

	uniform_real_distribution<double> double_dist(0,max_weight);


	//run the resampling wheel
	for (int i = 0; i< this->num_particles; i++)
	{

		beta = beta + 2.0 * double_dist(gen);

		while (weights[index] < beta)
		{
			beta = beta - weights[index];
			index = (index+1) % this->num_particles;
		}

		resampled_particles.push_back(particles[index]);
	}

	particles = resampled_particles;


}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;

    return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
