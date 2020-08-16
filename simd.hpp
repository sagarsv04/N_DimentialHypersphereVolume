/*
 *  simd.hpp
 *  Contains declaration of Line Segments in 4 Dimensional
 *
 *  Author :
 *  Sagar Vishwakarma (svishwa2@binghamton.edu)
 *  State University of New York, Binghamton
 */

#ifndef _SIMD_HPP_
#define _SIMD_HPP_


#include <cstdio>
#include <cmath>
#include <numeric>
#include <random>
#include <vector>
#include <iostream>

using namespace std;

vector<float> get_random_point(int range, int dimension);
vector<vector<float>> generate_line_points(int range, int dimension);
void print_distance_value(float* dist);

void calculate_dist_seq(vector<vector<float>>& a, vector<vector<float>>& b, float* dist);
void calculate_dist_simd(vector<vector<float>>& a, vector<vector<float>>& b, float* dist);

void compute_line_segment(int range, int dimension);


#endif /* _SIMD_HPP_ */
