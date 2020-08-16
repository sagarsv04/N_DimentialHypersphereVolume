#!/usr/bin/python
import os
import sys
import argparse
from argparse import RawTextHelpFormatter
from tqdm import tqdm
import numpy as np
from scipy.special import gammainc
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

np.random.seed(42)


def two_dim_plot_sample(sample_points, center, radius):

	fig1 = plt.figure(1)
	ax1 = fig1.gca()
	ax1.scatter(sample_points[:,0], sample_points[:,1], s=0.5)
	ax1.add_artist(plt.Circle(center, radius, fill=False, color='0.5'))
	ax1.set_xlim(-1.5,1.5)
	ax1.set_ylim(-1.5,1.5)
	ax1.set_aspect("equal")
	plt.title("A 2D Representation of {0}D points".format(center.size))
	plt.tight_layout()
	plt.show()
	return 0


def three_dim_plot_sample(sample_points, center, radius):

	phi, theta = np.mgrid[0.0:np.pi:100j, 0.0:2.0*np.pi:100j]
	x = radius*np.sin(phi)*np.cos(theta)
	y = radius*np.sin(phi)*np.sin(theta)
	z = radius*np.cos(phi)
	#Set colours and render
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.plot_surface(x, y, z,  rstride=1, cstride=1, color='c', alpha=0.3, linewidth=0)
	ax.scatter(sample_points[:,0], sample_points[:,1], sample_points[:,2], color="k", s=0.5)
	ax.set_xlim([-1,1])
	ax.set_ylim([-1,1])
	ax.set_zlim([-1,1])
	if sys.platform == "linux":
		ax.set_aspect("equal") # not working in matplotlib > 3.0.x
	plt.title("A 3D Representation of {0}D points".format(center.size))
	plt.tight_layout()
	plt.show()
	return 0


def plot_sample_histogram(sample_points, center, radius):

	dist_bucket = np.arange(0, radius+0.01, 0.01*radius)
	dist_count = np.array([0]*dist_bucket.shape[0])
	for idx in range(sample_points.shape[0]):
		# idx = 0
		dist = abs(radius-np.linalg.norm(center-sample_points[idx]))
		dist = round(dist, 2)
		if np.where(dist_bucket==dist)[0].shape[0]>0:
			dist_count[np.where(dist_bucket==dist)[0][0]] += 1

	fig = plt.figure()
	ax = plt.axes()
	# ax.bar(dist_bucket, dist_count, width = 0.5, color='#0504aa',alpha=0.7)
	ax.plot(dist_bucket, dist_count, marker="H", markersize=5);
	ax.set(xlabel="Distance from Surface", ylabel="Number Of Points", title="Histogram of Points Distance in {0} Dimensions".format(center.shape[0]));
	# plt.xticks(dist_bucket)
	plt.grid()
	plt.show()
	return 0


def plot_dimensional_distance_plot(dist_array, dimention):

	dimention_array = np.array(list(range(1, dimention+1)))
	fig = plt.figure()
	ax = plt.axes()
	# ax.bar(dimention_array, dist_array, width = 0.5, color='#0504aa',alpha=0.7)
	ax.plot(dimention_array, dist_array, marker="H", markersize=5);
	ax.set(xlabel="Num of Dimensions", ylabel="Distance from Surface", title="Distance between points and surface of Hypersphere");
	plt.grid()
	plt.show()
	return 0


def generate_sample(center, radius, num_sample):

	num_dimension = center.size
	random_num = np.random.normal(size=(num_sample, num_dimension))
	# random_num[0]
	squared_random_num = random_num**2
	# squared_random_num[0]
	squared_sum = np.sum(squared_random_num, axis=1)
	# squared_sum[0]
	fractional_radius = radius*gammainc(num_dimension/2, squared_sum/2)**(1/num_dimension)/np.sqrt(squared_sum)
	ndim_fractional_ratius = np.tile(fractional_radius.reshape(num_sample, 1), (1, num_dimension))
	sample_points = center + np.multiply(random_num, ndim_fractional_ratius)

	return sample_points


def create_hypersphere(sample_points, center, radius):
	plot_sample_histogram(sample_points, center, radius)
	two_dim_plot_sample(sample_points, center, radius)
	three_dim_plot_sample(sample_points, center, radius)
	return 0


def get_dimensional_distance_array(radius, dimention, num_sample):

	dist_array = np.array([0.0]*dimention)
	for dim in tqdm(range(1, dimention+1)):
		# dim = 1
		center = np.array([0]*dim)
		sample_points = generate_sample(center, radius, num_sample)
		dist = abs(radius-np.linalg.norm(center-sample_points))
		dist_array[dim-1] = round(dist, 4)
		if dim==dimention:
			create_hypersphere(sample_points, center, radius)

	return dist_array


def create_dimensional_distance(radius, dimention, num_sample):
	# radius, dimention, num_sample = 1, 50, 1000000
	dist_array = get_dimensional_distance_array(radius, dimention, num_sample)
	plot_dimensional_distance_plot(dist_array, dimention)

	return 0


def main():
	parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
	parser.add_argument("radius", type=int, help="Radius")
	parser.add_argument("num_dimension", type=int, help="Number of dimension")
	parser.add_argument("num_sample", type=int, help="Number of sample points")

	args = parser.parse_args()
	create_dimensional_distance(args.radius, args.num_dimension, args.num_sample)

	return 0



if __name__ == '__main__':
	main()
