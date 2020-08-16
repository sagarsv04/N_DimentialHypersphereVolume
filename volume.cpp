/*
 *  volume.cpp
 *  Contains implementation of N Dimensional Volume algorithm
 *
 *  Author :
 *  Sagar Vishwakarma (svishwa2@binghamton.edu)
 *  State University of New York, Binghamton
 */

#include "volume.hpp"
#include <functional>
#include <algorithm>
#include <iterator>
#include <chrono>
#include <omp.h>

#define DEBUG 1

#define PREVIOUS_DIM 0

#define NUM_THREADS 8


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


vector<float> get_random_point(int radius, int dimension) {

	vector<float> r_num(dimension, 0.0);

	uniform_real_distribution<float>  distr(0.0, (float)radius);

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


vector<vector<float>> generate_sample(vector<float> center, int radius, int num_sample, int* point_dist) {

	int nthreads;
	int per_thr;
	int dimension = center.size();
	vector<vector<float>> sample_points(num_sample, vector<float>(dimension, 0.0));
	vector<int> points_dist_count(101, 0);
	float points_distance = 0.0;

	omp_set_num_threads(NUM_THREADS);
	#pragma omp parallel
	{
		vector<float> thr_r_num(dimension, 0.0);
		vector<float> thr_center(dimension, 0.0);
		vector<int> thr_dist_count(101, 0);
		float thr_points_dist = 0.0;
		float thr_point_dist;

		int start_idx;
		int end_idx;
		int thr_idx;
		int num_thr;

		num_thr = omp_get_num_threads();
		thr_idx = omp_get_thread_num();

		if (thr_idx == 0) {
			nthreads = num_thr;
		}
		per_thr = num_sample/num_thr;

		start_idx = thr_idx*per_thr;

		if (thr_idx==num_thr-1) {
			end_idx = num_sample;
		}
		else {
			end_idx = (thr_idx+1)*per_thr;
		}
		if (DEBUG==2) {
			printf("Thread %d, Start Index %d, End Index %d\n",thr_idx, start_idx, end_idx);
		}
		for (int num_idx = start_idx; num_idx < end_idx; num_idx++) {
			thr_r_num = get_random_point(radius, dimension);
			thr_point_dist = vectors_distance(thr_center, thr_r_num);
			thr_point_dist = fabs((float)radius - thr_point_dist);
			thr_points_dist += thr_point_dist;
			thr_point_dist = roundf(thr_point_dist*100)/100;
			if (DEBUG==2) {
				printf("%d, Thread Dist %f, Point Dist %f\n",thr_idx, thr_points_dist, thr_point_dist);
			}
			thr_point_dist = thr_point_dist*100.0;
			if (thr_point_dist > 100.0) {
				thr_point_dist = 100.0;
			}
			thr_dist_count[(int)thr_point_dist] += 1;
			sample_points[num_idx] = thr_r_num;
		}
		#pragma omp critical
		{
			for (int i = 0; i < 101; i++) {
				points_dist_count[i] += thr_dist_count[i];
			}
			points_distance += thr_points_dist;
		}
	}

	*point_dist = (int)points_distance;

	if (DEBUG==2) {
		printf("%4d Threads Requested, %4d Threads Executed\n",NUM_THREADS, nthreads);
	}

	if (DEBUG) {
		printf("Histogram of Points Distance From Surface in %dD\n",dimension);
		for (int i = 0; i < 101; i++) {
			printf("%2f :: %d\n", i*0.01, points_dist_count[i]);
		}
	}

	return sample_points;
}


vector<int> get_dimensional_distance_array(int radius, int num_dimension, int num_sample) {

	vector<vector<float>> sample_points;
	vector<int> dist_array(num_dimension, 0.0);
	int start_dim = num_dimension;
	int points_dist = 0;

	if (PREVIOUS_DIM) {
		start_dim = 1;
	}

	for (int dim = start_dim; dim < num_dimension + 1; dim++) {
		vector<float> center(dim, 0.0);
		// generate samples
		sample_points = generate_sample(center, radius, num_sample, &points_dist);
		dist_array[dim-1] = points_dist;
		if (DEBUG) {
			printf("Distance of Points in %2dD is: %d\n", dim, dist_array[dim-1]);
		}
	}

	return dist_array;
}


int create_dimensional_distance(int radius, int num_dimension, int num_sample) {

	vector<int> dist_array;

	dist_array = get_dimensional_distance_array(radius, num_dimension, num_sample);

	return 0;
}


int main(int argc, char const* argv[]) {

	if (argc != 4) {
		cerr << "Help : Usage "<< argv[0] << " radius num_dimension num_sample" << endl;
		exit(1);
	}
	else {
		int radius;
		int ret_val;
		int num_dimension;
		int num_sample;
		// assign argv values to variables
		radius	= atoi(argv[1]);
		num_dimension	= atoi(argv[2]);
		num_sample	= atoi(argv[3]);

		if (radius<=0) {
			cerr << "Error : Radius should be greater than 0" << endl;
			exit(1);
		}
		if (num_dimension<=0) {
			cerr << "Error : Number of dimension should be greater than 0" << endl;
			exit(1);
		}
		if (num_sample<=0) {
			cerr << "Error : Number of sample should be greater than 0" << endl;
			exit(1);
		}
		else {
			// code goes here
			struct rusage res_usage;
			printf("Radius :: %d\n", radius);
			printf("Number of dimension :: %d\n", num_dimension);
			printf("Number of sample :: %d\n", num_sample);

			auto t1 = chrono::high_resolution_clock::now();
			ret_val = create_dimensional_distance(radius, num_dimension, num_sample);
			auto t2 = chrono::high_resolution_clock::now();
			auto duration = chrono::duration_cast<chrono::microseconds>(t2-t1).count();
			cout << "Total Execution Time :: " << duration << "\u03BC" << " sec" << '\n';
			ret_val = getrusage(RUSAGE_SELF, &res_usage); assert(ret_val == 0);
			auto cv = [](const timeval &time_val) {
				return double(time_val.tv_sec) + double(time_val.tv_usec)/1000000;
			};
			cerr <<"Resource Usage: \n";
			cerr << "    User CPU Time: " << cv(res_usage.ru_utime) << '\n';
			cerr << "    Sys CPU Time: " << cv(res_usage.ru_stime) << '\n';
			cerr << "    Max Resident: " << res_usage.ru_maxrss << '\n';
			cerr << "    Page Faults: " << res_usage.ru_majflt << '\n';
		}
	}
	return 0;
}
