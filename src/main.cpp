//
// Created by nilsiism on 20/08/18.
//
#include "iostream"
#include <mav_tube_trajectory_generation/polynomial_optimization_linear.h>
#include <time.h>

#include <mav_tube_trajectory_generation/polynomial_optimization_qcqp.h>
#include <mav_tube_trajectory_generation/polynomial_optimization_nonlinear.h>

#include <se/octree.hpp>
#include <se/node.hpp>
#include <se/volume_traits.hpp>

int main () {

  std::string read_filename = "/home/nilsiism/workspace/prob_trajectory_planning/data/writen_data/fr_078_tidyup_supereight.bt";
  se::Octree<OFusion>* octree_ = new se::Octree<OFusion>;
  octree_->load(read_filename);

  mav_trajectory_generation::Vertex::Vector vertices;
  const int dimension = 3;
  const int derivative_to_optimize = mav_trajectory_generation::derivative_order::SNAP;
  mav_trajectory_generation::Vertex start(dimension), middle1(dimension), middle2(dimension), middle3(dimension), middle4(dimension), end(dimension);

  start.makeStartOrEnd(Eigen::Vector3d(9.8, 11.8, 12.45), derivative_to_optimize);
  vertices.push_back(start);

  middle1.addConstraint(mav_trajectory_generation::derivative_order::POSITION, Eigen::Vector3d(10.4474, 14.05, 12.7796));
  vertices.push_back(middle1);

  middle2.addConstraint(mav_trajectory_generation::derivative_order::POSITION, Eigen::Vector3d(13.0126, 14.7041, 12.3976));
  vertices.push_back(middle2);

  middle3.addConstraint(mav_trajectory_generation::derivative_order::POSITION, Eigen::Vector3d(13.6771, 15.358, 12.8908));
  vertices.push_back(middle3);

  middle4.addConstraint(mav_trajectory_generation::derivative_order::POSITION, Eigen::Vector3d(15.2574, 16.6208, 12.2477));
  vertices.push_back(middle4);

  end.makeStartOrEnd(Eigen::Vector3d(16.6, 16.6, 12.45), derivative_to_optimize);
  vertices.push_back(end);

  std::vector<double> segment_times;
  const double v_max = 2.0;
  const double a_max = 2.0;
  segment_times = estimateSegmentTimes(vertices, v_max, a_max);

  std::vector<std::pair<double, double>> segment_radii;
  std::pair<double, double> r1(0.8, 0.3);
  std::pair<double, double> r2(0.3, 0.3);
  std::pair<double, double> r3(0.3, 0.3);
  std::pair<double, double> r4(0.3, 0.3);
  std::pair<double, double> r5(0.3, 0.3);
  segment_radii.push_back(r1);
  segment_radii.push_back(r2);
  segment_radii.push_back(r3);
  segment_radii.push_back(r4);
  segment_radii.push_back(r5);

  clock_t t;
  t = clock();

  const int N = 10;

  mav_trajectory_generation::NonlinearOptimizationParameters parameters;
  parameters.solve_with_position_constraint = true;
  parameters.objective = mav_trajectory_generation::NonlinearOptimizationParameters::OptimizationObjective::kOptimizeFreeConstraintsAndCollision;
  parameters.print_debug_info = false;
  parameters.max_iterations = 200;
  parameters.max_time = -1.0;
  parameters.f_rel = 0.000001;
  parameters.x_rel = 0.01;
  parameters.use_soft_constraints = false;
  parameters.soft_constraint_weight = 100.0;
  parameters.time_penalty = 500.0;
  parameters.initial_stepsize_rel = 0.1;
  parameters.inequality_constraint_tolerance = 0.1;
  parameters.weights.w_d = 0.1;
  parameters.weights.w_c = 50;
  parameters.weights.w_t = 50;
  parameters.weights.w_sc = 1.0;
  parameters.use_numeric_grad = false;
  parameters.use_continous_distance = false;
  parameters.increment_time = 0.00001;
  parameters.epsilon = 0.3;
  parameters.coll_pot_multiplier = 20.0;
  parameters.is_collision_safe = true;
  parameters.is_simple_numgrad_time = true;
  parameters.is_simple_numgrad_constraints = true;
  parameters.coll_check_time_increment = 0.1;
  parameters.is_coll_raise_first_iter = true;
  parameters.add_coll_raise = 0.0000001;
  parameters.algorithm =
          static_cast<nlopt::algorithm>(11);
  parameters.random_seed = 12345678;
  parameters.map_resolution = 0.05;
  parameters.min_bound = Eigen::Vector3d(6.4, 5.8, 10.45) ;
  parameters.max_bound = Eigen::Vector3d(19.15, 19.75, 15.1);
  parameters.side = 0.5;


  mav_trajectory_generation::PolynomialOptimizationNonLinear<N> opt_nonlin(3, parameters);
  mav_trajectory_generation::PolynomialOptimizationConstrained<N> opt_constr(dimension);

  opt_nonlin.setOctree(octree_);
  opt_nonlin.setupFromVertices(vertices, segment_times, segment_radii, derivative_to_optimize);

  opt_nonlin.optimize();

  opt_constr.setupFromVertices(vertices, segment_times, segment_radii, derivative_to_optimize);
  opt_constr.solveLinear();
  opt_constr.solveQCQP();

  mav_trajectory_generation::Trajectory traj;
  opt_constr.getTrajectory(&traj);

  std::vector<Eigen::VectorXd> free_constraints;
  opt_constr.getFreeConstraints(&free_constraints);

  for (Eigen::VectorXd df : free_constraints)
    std::cout << df << "\n" << std::endl;

  t = clock() - t;
  std::cout << "Time to solve optimization [sec]: " << ((float)t)/CLOCKS_PER_SEC << std::endl;


  return 0;
};
