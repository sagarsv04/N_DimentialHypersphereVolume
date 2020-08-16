# Volume of N Dimential Hypersphere
Link: http://www.cs.binghamton.edu/~kchiu/cs547/prog/3/


## Author :

Sagar Vishwakarma (svishwa2@binghamton.edu)

State University of New York, Binghamton


## File :

1)	Makefile                     - Compile the program
2)	volume.cpp                   - Contains implementation of N Dimensional Volume algorithm
3)	volume.hpp                   - Contains declaration of N Dimensional Volume algorithm
4)	volume.py                    - Contains python implementation of algorithm


## Run :

- Open a terminal in project directory      : make (to build project)
- To run C++ volume algorithm               : ./volume radius num_dimension num_sample
- To run python volume algorithm            : python volume.py radius num_dimension num_sample



## Report :

I ran the code for 50 dimension with 1000000 sample points:

File : points_histogram.png shows how most of the points lies close to the surface of hypersphere.

File : distance_from_surface.png shows how total distance of all the points increase as we move in higher dimensions.

File : 2D_representation.png shows 2 dimensional representation of points in 50 dimensional hypersphere.

File : 3D_representation.png shows 2 dimensional representation of points in 50 dimensional hypersphere.


## Note :

Part A :
- Both sequential and parallel implementation (using OpenMP) works by changing NUM_THREADS to 1 to number of threads you want.
- This parameter is inside the volume.cpp file. Also use DEBUG to see the prints.
- I have also created a python file which computes the volume of hypersphere, and plots different graph.

Part B :
- Both sequential and simd implementation executing one after the other.
- Use DEBUG to see the prints. This parameter is inside the simd.cpp file.
- Use NUM_SAMPLE to change number of samples. This parameter is inside the simd.cpp file.
