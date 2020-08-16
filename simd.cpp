/*
 *  simd.cpp
 *  Contains implementation of Line Segment in 4 Dimensional
 *
 *  Author :
 *  Sagar Vishwakarma (svishwa2@binghamton.edu)
 *  State University of New York, Binghamton
 */

#include "simd.hpp"
#include <immintrin.h>
#include <xmmintrin.h>
#include <functional>
#include <algorithm>
#include <iterator>
#include <chrono>

#define DEBUG 0

const int NUM_SAMPLE = 16*1000000;
// const int NUM_SAMPLE = 16*10;

random_device rand_dev;
mt19937 generator(rand_dev());


// Computes the distance between two vectors
template <typename T>
float	vectors_distance(const vector<T>& a, const vector<T>& b) {
	vector<float> auxiliary;
	transform(a.begin(), a.end(), b.begin(), back_inserter(auxiliary), [](T element1, T element2) {
		return pow((element1-element2),2);
	});
	// auxiliary.shrink_to_fit();
	return  sqrt(accumulate(auxiliary.begin(), auxiliary.end(), 0.0));
} // end template vectors_distance


vector<float> get_random_point(int range, int dimension) {

	vector<float> r_num(dimension, 0.0);

	uniform_real_distribution<float>  distr(0.0, (float)range);

	if (DEBUG==2) {
		printf("Random Point in %dD: [ ",dimension);
		for (int idx = 0; idx < dimension; idx++) {
			r_num[idx] = distr(generator);
			printf("%f, ", r_num[idx]);
		}
		printf(" ]\n");
	}
	else {
		for (int idx = 0; idx < dimension; idx++) {
			r_num[idx] = distr(generator);
		}
	}
	return r_num;
}

vector<vector<float>> generate_line_points(int range, int dimension) {

  vector<vector<float>> line_points(NUM_SAMPLE, vector<float>(dimension, 0.0));

  for (int num_idx = 0; num_idx < NUM_SAMPLE; num_idx++) {
    line_points[num_idx] = get_random_point(range, dimension);
  }

  return line_points;
}


void print_distance_value(float* dist) {

	for (int i = 0; i < NUM_SAMPLE; i++) {
		printf("Point : %d, Distance : %f\n", i, dist[i]);
	}
}


void calculate_dist_seq(vector<vector<float>>& a, vector<vector<float>>& b, float* dist) {

	auto t1 = chrono::high_resolution_clock::now();
	for (int i = 0; i < NUM_SAMPLE; i++) {
    dist[i] = vectors_distance(a[i], b[i]);
  }
	auto t2 = chrono::high_resolution_clock::now();
  if (DEBUG) {
    print_distance_value(dist);
  }
	auto duration = chrono::duration_cast<chrono::microseconds>(t2-t1).count();
  cout << "Seq Execution Time :: " << duration << "\u03BC" << " sec" << '\n';
	cout << "Sequential: " << (NUM_SAMPLE/duration) << " Mops/" << "\u03BC" << "s" << endl;
}


void calculate_dist_simd(vector<vector<float>>& a, vector<vector<float>>& b, float* dist) {

	static float a_1[NUM_SAMPLE];
	static float a_2[NUM_SAMPLE];
	static float a_3[NUM_SAMPLE];
	static float a_4[NUM_SAMPLE];

	static float b_1[NUM_SAMPLE];
	static float b_2[NUM_SAMPLE];
	static float b_3[NUM_SAMPLE];
	static float b_4[NUM_SAMPLE];

	for (int i = 0; i < NUM_SAMPLE; i++) {
		a_1[i] = a[i][0];
		a_2[i] = a[i][1];
		a_3[i] = a[i][2];
		a_4[i] = a[i][3];

		b_1[i] = b[i][0];
		b_2[i] = b[i][1];
		b_3[i] = b[i][2];
		b_4[i] = b[i][3];
	}

	alignas(32) static float dist_val[NUM_SAMPLE];

	auto t1 = chrono::high_resolution_clock::now();

	for (int i = 0; i < NUM_SAMPLE/4; i++) {
		__m128 ax = _mm_load_ps(a_1 + 4*i);
		__m128 ay = _mm_load_ps(a_2 + 4*i);
		__m128 az = _mm_load_ps(a_3 + 4*i);
		__m128 aw = _mm_load_ps(a_4 + 4*i);

		__m128 bx = _mm_load_ps(b_1 + 4*i);
		__m128 by = _mm_load_ps(b_2 + 4*i);
		__m128 bz = _mm_load_ps(b_3 + 4*i);
		__m128 bw = _mm_load_ps(b_4 + 4*i);

		__m128 ax_sub_bx = _mm_sub_ps(ax, bx);
		__m128 ay_sub_by = _mm_sub_ps(ay, by);
		__m128 az_sub_bz = _mm_sub_ps(az, bz);
		__m128 aw_sub_bw = _mm_sub_ps(aw, bw);

		__m128 value = _mm_sqrt_ps(_mm_mul_ps(ax_sub_bx, ax_sub_bx)
															+ _mm_mul_ps(ay_sub_by, ay_sub_by)
															+ _mm_mul_ps(az_sub_bz, az_sub_bz)
															+ _mm_mul_ps(aw_sub_bw, aw_sub_bw));
		_mm_store_ps(dist_val + 4*i, value);
	}

  for (int i = 0; i < NUM_SAMPLE; i++) {
		dist[i]	= dist_val[i];
  }
	auto t2 = chrono::high_resolution_clock::now();
  if (DEBUG) {
    print_distance_value(dist);
  }
  auto duration = chrono::duration_cast<chrono::microseconds>(t2-t1).count();
  cout << "SIMD Execution Time :: " << duration << "\u03BC" << " sec" << '\n';
	cout << "Sequential: " << (NUM_SAMPLE/duration) << " Mops/" << "\u03BC" << "s" << endl;
}



void compute_line_segment(int range, int dimension) {

  vector<vector<float>> a;
  vector<vector<float>> b;
	float *dist = new float[NUM_SAMPLE];
  // float dist[NUM_SAMPLE];

  a = generate_line_points(range, dimension);
  b = generate_line_points(range, dimension);
  // a b together form line segments
  calculate_dist_seq(a, b, dist);
  calculate_dist_simd(a, b, dist);

	delete [] dist;
}


int main(int argc, char const* argv[]) {

	if (argc > 1) {
		cerr << "Help : Usage "<< argv[0] << endl;
		cerr << "NUM_SAMPLE : Variable To Change Number of Samples" << endl;
		exit(1);
	}
	else {
  	compute_line_segment(1, 4);
  }

  return 0;
}
