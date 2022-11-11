/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles

  //pass x,y,theta noise
  double std_x = std[0];
  double std_y = std[1];
  double std_theta = std[2];

  //creat gaussian
  using std::normal_distribution;
  

  normal_distribution<double> dist_x(x, std_x);
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> dist_theta(theta, std_theta);

  //generate particles using gaussian distribution function
  std::default_random_engine gen;  

  for(int i=0;i<num_particles;i++){
    particles[i].id = i;
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);
    particles[i].weight = 1.;
  }
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  for(int i=0;i<num_particles;i++){
    particles[i].x = particles[i].x + velocity/yaw_rate*(sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta));
    particles[i].y = particles[i].y + velocity/yaw_rate*(cos(particles[i].theta)-cos(particles[i].theta)+yaw_rate*delta_t);
    particles[i].theta = particles[i].theta + yaw_rate*delta_t;
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

  double dist_= INFINITY;
  
  for(int i=0;i<observations.size();i++){
    for(int j=0;j<predicted.size();j++){
      if(dist(predicted[j].x,predicted[j].y,observations[i].x,observations[i].y)<dist_){
        dist_ = dist(predicted[j].x,predicted[j].y,observations[i].x,observations[i].y);
        observations[i].id = predicted[j].id;        
      }     
    }
  }  
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  
  vector<LandmarkObs> predicted;
  vector<LandmarkObs> t_obs;

  for (int i=0; i<particles.size();i++){

    //find landmarks with particle sense range
    for(int id=0;id<map_landmarks.landmark_list.size();id++){
      if(dist(map_landmarks.landmark_list[id].x_f,map_landmarks.landmark_list[id].y_f,particles[i].x,particles[i].y)<sensor_range){

        predicted[id].id = map_landmarks.landmark_list[id].id_i;
        predicted[id].x = map_landmarks.landmark_list[id].x_f;
        predicted[id].y = map_landmarks.landmark_list[id].y_f;
        
      }
    }

    //transform observations to map coordinate
    for (int j=0;j< observations.size();j++){
      t_obs[j].x=particles[i].x + cos(particles[i].theta)*observations[j].x - sin(particles[i].theta)*observations[j].y;
      t_obs[j].y=particles[i].y + sin(particles[i].theta)*observations[j].x + cos(particles[i].theta)*observations[j].y;
      t_obs[j].id=observations[j].id;      
    }

    

    //data associate
    dataAssociation(predicted,t_obs);
    
    
    //calculate particle weight
  
    
    
    
  for(int i=0;i<particles.size();i++){
      int id = particles[i].associations;
      double x_obs = observations[ass_id].x;
      double y_obs = observations[ass_id].y;
      double x_tobs = particles[i].sense_x[ass_id];
      double y_tobs = particles[i].sense_y[ass_id];
      double std_x = std_landmark[0];
      double std_y = std_landmark[1];

      particles[i].weight += (1/2*M_PI*std_x*std_y)*exp(-(pow(x_obs-x_tobs,2)/2*std_x*std_x + pow(y_obs-y_tobs,2)/2*std_y*std_y));
    }  
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}