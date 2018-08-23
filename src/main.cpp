//
// Created by nilsiism on 20/08/18.
//
#include "iostream"
#include <mav_tube_trajectory_generation/polynomial_optimization_linear.h>
#include <time.h>
#include <mav_tube_trajectory_generation/polynomial_optimization_linear_constrained.h>

int main () {
  mav_trajectory_generation::Vertex::Vector vertices;
  const int dimension = 3;
  const int derivative_to_optimize = mav_trajectory_generation::derivative_order::SNAP;
  mav_trajectory_generation::Vertex start(dimension), middle1(dimension), middle2(dimension), middle3(dimension), end(dimension);

  start.makeStartOrEnd(Eigen::Vector3d(0,0,0), derivative_to_optimize);
  vertices.push_back(start);

  middle1.addConstraint(mav_trajectory_generation::derivative_order::POSITION, Eigen::Vector3d(1,1,2));
  vertices.push_back(middle1);

  middle2.addConstraint(mav_trajectory_generation::derivative_order::POSITION, Eigen::Vector3d(2,2,1));
  vertices.push_back(middle2);

  middle3.addConstraint(mav_trajectory_generation::derivative_order::POSITION, Eigen::Vector3d(3,3,3));
  vertices.push_back(middle3);

  end.makeStartOrEnd(Eigen::Vector3d(3,4,2), derivative_to_optimize);
  vertices.push_back(end);

  std::vector<double> segment_times;
  const double v_max = 2.0;
  const double a_max = 2.0;
  segment_times = estimateSegmentTimes(vertices, v_max, a_max);

  std::vector<std::pair<double, double>> segment_radii;
  std::pair<double, double> r1(0.5, 0.5);
  std::pair<double, double> r2(0.5, 0.5);
  std::pair<double, double> r3(0.5, 0.5);
  std::pair<double, double> r4(0.5, 0.5);
  segment_radii.push_back(r1);
  segment_radii.push_back(r2);
  segment_radii.push_back(r3);
  segment_radii.push_back(r4);

  clock_t t;
  t = clock();

  const int N = 10;
  mav_trajectory_generation::PolynomialOptimizationConstrained<N> opt_constr(dimension);
  opt_constr.setupFromVertices(vertices, segment_times, segment_radii, derivative_to_optimize);
  opt_constr.solveLinear();
  opt_constr.solveQCQP();

  t = clock() - t;
  std::cout << "sec" << ((float)t)/CLOCKS_PER_SEC << std::endl;
  std::cout << "test" << std::endl;


  return 0;
};
