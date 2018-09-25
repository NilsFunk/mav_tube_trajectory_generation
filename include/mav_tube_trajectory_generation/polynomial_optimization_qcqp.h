/*
 * Probabilistic Trajectory Planning Quadraticlly Constraind Quadratic Program.
 *
 * Copyright (C) 2018 Imperial College London.
 * Copyright (C) 2018 ETH Zürich.
 *
 * @todo LICENSE
 * 
 * @original file polynomial_optimization_linear.h
 * Copyright (c) 2016, Markus Achtelik, ASL, ETH Zurich, Switzerland
 * Copyright (c) 2016, Michael Burri, ASL, ETH Zurich, Switzerland
 * Copyright (c) 2016, Helen Oleynikova, ASL, ETH Zurich, Switzerland
 * Copyright (c) 2016, Rik Bähnemann, ASL, ETH Zurich, Switzerland
 * Copyright (c) 2016, Marija Popovic, ASL, ETH Zurich, Switzerland
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * @file file polynomial_optimization_qcqp.h
 * @author Nils Funk
 * @date August, 2018
 */

#ifndef MAV_TUBE_TRAJECTORY_GENERATION_POLYNOMIAL_OPTIMIZATION_QCQP_H
#define MAV_TUBE_TRAJECTORY_GENERATION_POLYNOMIAL_OPTIMIZATION_QCQP_H

#include <mav_tube_trajectory_generation/polynomial_optimization_linear.h>
#include <iostream>
#include <stdio.h>
#include <h/fusion.h>
#include <Eigen/Cholesky>

extern "C"{
  #include "h/mosek.h"
}

using namespace mosek::fusion;
using namespace monty;

static void MSKAPI printstr(void *handle,
                            const char str[])
{
  printf("%s", str);
} /* printstr */

namespace mav_trajectory_generation {

  template <int _N = 10>
  class PolynomialOptimizationConstrained : public PolynomialOptimization<_N> {
  public:
    PolynomialOptimizationConstrained(size_t dimension);

    bool setupFromVertices(const Vertex::Vector& vertices, const std::vector<double>& times,
            const std::vector<std::pair<double, double>>& radii, int derivative_to_optimize);

    bool solveLinear();
    int solveQCQP();
    int solveSOCP();

    void getSegmentRadii(std::vector<std::pair<double, double>>* segment_radii) const {
      CHECK(segment_radii != nullptr);
      *segment_radii = segment_radii_;
    }

  private:
    void setupConstraintReorderingMatrixkDim();
    void constructRkDim(Eigen::SparseMatrix<double>* R_kDim) const;
    void setupInverseControlPointMappingMatrix(typename PolynomialOptimization<_N>::SquareMatrix* B, double T);

    void setupControlPointConstraints();
    void compute_sphere_constraints(std::vector<Eigen::MatrixXd> &control_pt_extraction_matrices, int segment_idx);
    void compute_tube_constraints(std::vector<Eigen::MatrixXd> &control_pt_extraction_matrices, int segment_idx);
    void compute_tube_end_constraints(std::vector<Eigen::MatrixXd> &control_pt_extraction_matrices, int segment_idx);

    int factorial(int n);
    int binomialCoeff(int n, int k);

    Eigen::VectorXd free_constraints_compact_kDim_;
    Eigen::VectorXd fixed_constraints_compact_kDim_;
    std::vector<std::pair<double, double>> segment_radii_;
    size_t n_all_constraints_kDim_;
    size_t n_fixed_constraints_kDim_;
    size_t n_free_constraints_kDim_;
    Eigen::SparseMatrix<double>  constraint_reordering_kDim_;
    typename PolynomialOptimization<_N>::SquareMatrixVector inverse_control_pt_mapping_matrices_;
    std::vector<Eigen::MatrixXd> constr_quad_;
    std::vector<Eigen::MatrixXd> constr_lin_;
    std::vector<double> constr_const_;

    using PolynomialOptimization<_N>::n_free_constraints_;
    using PolynomialOptimization<_N>::n_fixed_constraints_;
    using PolynomialOptimization<_N>::free_constraints_compact_;
    using PolynomialOptimization<_N>::fixed_constraints_compact_;
    using PolynomialOptimization<_N>::n_all_constraints_;
    using PolynomialOptimization<_N>::n_segments_;
    using PolynomialOptimization<_N>::n_vertices_;
    using PolynomialOptimization<_N>::N;
    using PolynomialOptimization<_N>::dimension_;
    using PolynomialOptimization<_N>::vertices_;
    using PolynomialOptimization<_N>::segments_;
    using PolynomialOptimization<_N>::constraint_reordering_;
    using PolynomialOptimization<_N>::derivative_to_optimize_;
    using PolynomialOptimization<_N>::segment_times_;
    using PolynomialOptimization<_N>::inverse_mapping_matrices_;
    using PolynomialOptimization<_N>::cost_matrices_;
    using PolynomialOptimization<_N>::kHighestDerivativeToOptimize;
  };

} //namespace mav_tube_trajectory_generation

#include <mav_tube_trajectory_generation/impl/polynomial_optimization_qcqp_impl.h>



#endif //MAV_TUBE_TRAJECTORY_GENERATION_POLYNOMIAL_OPTIMIZATION_QCQP_H
