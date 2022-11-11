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

  double std_x = std[0];
  double std_y = std[1];
  double std_theta = std[2];

  std::normal_distribution<double> dist_x(x, std_x);
  std::normal_distribution<double> dist_y(y, std_y);
  std::normal_distribution<double> dist_theta(theta, std_theta);
  std::default_random_engine gen;  

  for(unsigned int i=0;i<num_particles;i++){
    Particle particle;

    particle.id = i;
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1.0;

    particles.push_back(particle);
    weights.push_back(1.0);
  }

  is_initialized=true;
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

  std::default_random_engine gen;    

  //add noise to prediction
    
  std::normal_distribution<double> dist_x(0, std_pos[0]);
  std::normal_distribution<double> dist_y(0, std_pos[1]);
  std::normal_distribution<double> dist_theta(0, std_pos[2]);

  for(unsigned int i=0;i<num_particles;i++){   

    double x_pred, y_pred, theta_pred; 

    //assign current particle state
    double x_0=particles[i].x;
    double y_0=particles[i].y;
    double theta_0=particles[i].theta;   
    

    if(fabs(yaw_rate)>0.0001){
      x_pred = x_0 + (velocity/yaw_rate)*(std::sin(theta_0+yaw_rate*delta_t) - std::sin(theta_0));
      y_pred = y_0 + (velocity/yaw_rate)*(std::cos(theta_0) - std::cos(theta_0+yaw_rate*delta_t));
      theta_pred = theta_0 + yaw_rate*delta_t;
    }

    else{
      x_pred = x_0 + velocity*delta_t*cos(theta_0);
      y_pred = y_0 + velocity*delta_t*sin(theta_0);
      theta_pred = theta_0;
    }    
    
    particles[i].x = x_pred + dist_x(gen);
    particles[i].y = y_pred + dist_y(gen);
    particles[i].theta = theta_pred + dist_theta(gen); 
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

  
  //find closest dist_ to each observation
  for(unsigned int i=0;i<observations.size();i++){
    LandmarkObs obs=observations[i];
    int new_id = obs.id;

    double current_euclidian_distance = std::numeric_limits<double>::max();

    for(unsigned int j=0;j<predicted.size();j++){

      LandmarkObs pred = predicted[j];
            
      if(dist(obs.x,obs.y,pred.x,pred.y)<current_euclidian_distance){
        new_id = pred.id;
        current_euclidian_distance=dist(obs.x,obs.y,pred.x,pred.y);
      }     
    }
    observations[i].id = new_id;
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

  for (unsigned int i=0; i<num_particles;i++){

    double xp = particles[i].x;
    double yp = particles[i].y;
    double theta = particles[i].theta; 

    vector<LandmarkObs> tf_observations;

    for (unsigned int j=0;j< observations.size();j++){      
      LandmarkObs tf_observation;
      double xc = observations[j].x;
      double yc = observations[j].y;

      tf_observation.x = xp + (std::cos(theta)*xc) - (std::sin(theta)*yc);
      tf_observation.y = yp + (std::sin(theta)*xc) + (std::cos(theta)*yc);
      tf_observation.id = j;

      tf_observations.push_back(tf_observation);
    }    

    //find landmarks within particle's sensor range then assign to predicted    
    vector<LandmarkObs> landmarks_in_range; 

    for(unsigned int k=0;k<map_landmarks.landmark_list.size();k++){   
      double distance = dist(xp,yp,map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f);

      if(distance<=sensor_range){
        LandmarkObs landmark_in_range;
        landmark_in_range.id = map_landmarks.landmark_list[k].id_i;
        landmark_in_range.x = map_landmarks.landmark_list[k].x_f;
        landmark_in_range.y = map_landmarks.landmark_list[k].y_f;

        landmarks_in_range.push_back(landmark_in_range);        
      }
    }    

    //data associate
    dataAssociation(landmarks_in_range,tf_observations);      
    
    //calculate particle weight   

    double sig_x = std_landmark[0];
    double sig_y = std_landmark[1];
    double var_x = sig_x*sig_x;
    double var_y = sig_y*sig_y;
    double gaussian_norm = 1.0/(2.0*M_PI*sig_x*sig_y);

    particles[i].weight=1.0;
    double p_=1.0;   

    for (int j=0;j<tf_observations.size();j++){

      double tf_x = tf_observations[j].x;
      double tf_y = tf_observations[j].y;
      double tf_id = tf_observations[j].id;

      for (int k=0;k<landmarks_in_range.size();k++){
        double ldmk_x =landmarks_in_range[k].x;
        double ldmk_y =landmarks_in_range[k].y;
        double ldmk_id = landmarks_in_range[k].id;

        if (tf_id==ldmk_id){
          double mu_x = ldmk_x;
          double mu_y = ldmk_y;

          double exp_x = ((tf_x-mu_x)*(tf_x-mu_x))/(2*var_x);
          double exp_y = ((tf_y-mu_y)*(tf_y-mu_y))/(2*var_y);

          p_ = p_*gaussian_norm*exp(-(exp_x+exp_y));
          break;
        }
      }    
    }    
    particles[i].weight=p_;    
    weights[i] = p_;
  }

  //normalize weights
  double weight_norm = std::accumulate(weights.begin(),weights.end(),0.0f); 
  
  for(int p=0;p<num_particles;p++){
    particles[p].weight = particles[p].weight/weight_norm;
    weights[p] = particles[p].weight;
  }  
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

  vector<Particle> particles_resample;
  
  std:: default_random_engine gen;
  std::discrete_distribution<int> sample(weights.begin(),weights.end());

  for (int i=0;i<num_particles;i++){
    int k = sample(gen);
    particles_resample.push_back(particles[k]);
  }

  particles = particles_resample;
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