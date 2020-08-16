/*
 *  volume.hpp
 *  Contains declaration of N Dimensional Volume algorithm
 *
 *  Author :
 *  Sagar Vishwakarma (svishwa2@binghamton.edu)
 *  State University of New York, Binghamton
 */

#ifndef _VOLUME_HPP_
#define _VOLUME_HPP_


#include <cstdio>
#include <cmath>
#include <numeric>
#include <random>
#include <vector>
#include <iostream>
#include <assert.h>
#include <sys/time.h>
#include <sys/resource.h>


using namespace std;

vector<float> get_random_point(int radius, int dimension);
vector<vector<float>> generate_sample(vector<float> center, int radius, int num_sample, int* point_dist);
vector<int> get_dimensional_distance_array(int radius, int num_dimension, int num_sample);
int create_dimensional_distance(int radius, int num_dimension, int num_sample);


#endif /* _VOLUME_HPP_ */
