//
// Created by nilsiism on 21/08/18.
//

#ifndef MAV_TUBE_TRAJECTORY_GENERATION_POLYNOMIAL_OPTIMIZATION_QCQP_H
#define MAV_TUBE_TRAJECTORY_GENERATION_POLYNOMIAL_OPTIMIZATION_QCQP_H

#include <mav_tube_trajectory_generation/polynomial_optimization_linear.h>
#include <iostream>
#include <stdio.h>

extern "C"{
  #include "h/mosek.h"
}

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


template <int _N>
PolynomialOptimizationConstrained<_N>::PolynomialOptimizationConstrained(size_t dimension)
    : PolynomialOptimization<_N>(dimension),
      n_all_constraints_kDim_(0),
      n_fixed_constraints_kDim_(0),
      n_free_constraints_kDim_(0) {
};

template <int _N>
void PolynomialOptimizationConstrained<_N>::setupConstraintReorderingMatrixkDim() {
  typedef Eigen::Triplet<double> Triplet;
  std::vector<Triplet> reordering_list_kDim;
  std::vector<Triplet> reordering_list;

  n_fixed_constraints_ = N;
  n_fixed_constraints_kDim_ = N * dimension_;

  n_free_constraints_ = (n_segments_-1) * N / 2;
  n_free_constraints_kDim_ = (n_segments_-1) * N / 2 * dimension_;

  n_all_constraints_ = n_fixed_constraints_ + 2 * n_free_constraints_;
  n_all_constraints_kDim_ = n_fixed_constraints_kDim_ + 2 * n_free_constraints_kDim_;

  reordering_list_kDim.reserve(n_all_constraints_kDim_);
  constraint_reordering_kDim_ = Eigen::SparseMatrix<double>(
          n_all_constraints_kDim_, n_fixed_constraints_kDim_ + n_free_constraints_kDim_);

  reordering_list.reserve(n_all_constraints_);
  constraint_reordering_ = Eigen::SparseMatrix<double>(
          n_all_constraints_, n_fixed_constraints_ + n_free_constraints_);

  fixed_constraints_compact_kDim_.resize(n_fixed_constraints_kDim_, Eigen::NoChange);
  for (Eigen::VectorXd& df : fixed_constraints_compact_)
    df.resize(n_fixed_constraints_, Eigen::NoChange);

  free_constraints_compact_kDim_.resize(n_free_constraints_kDim_, Eigen::NoChange);

  for (size_t ending = 0; ending < 2; ending++) {
    size_t vertex_idx = ending * n_segments_;
    const Vertex& vertex = vertices_[vertex_idx];
    for (size_t constraint_idx = 0; constraint_idx < N / 2;
         ++constraint_idx) {
      Constraint constraint;
      constraint.vertex_idx = vertex_idx;
      constraint.constraint_idx = constraint_idx;
      vertex.getConstraint(constraint_idx, &(constraint.value));
      for (size_t d = 0; d < dimension_; d++) {
        const Eigen::VectorXd constraint_value = constraint.value;
        size_t fixed_constraint_kDim_idx = ending * N / 2 + d * N + constraint_idx;
        fixed_constraints_compact_kDim_[fixed_constraint_kDim_idx] = constraint_value[d];

        size_t fixed_constraint_idx = ending * N / 2 + constraint_idx;
        Eigen::VectorXd& df = fixed_constraints_compact_[d];
        df[fixed_constraint_idx] = constraint_value[d];;
      }
    }
  }

  for (size_t d = 0; d<dimension_; d++){
    for (size_t k = 0; k < 2; k++) {
      size_t row = d * (2 * N / 2 + (n_segments_-1) * N) + k * (N / 2 + (n_segments_-1) * N);
      size_t col = d * 2 * N / 2 + k * N / 2;
      for (size_t i = 0; i< N / 2; i++) {
        reordering_list_kDim.emplace_back(Triplet(row, col, 1.0));
        row++;
        col++;
      }
    }
  }

  for (size_t d = 0; d < dimension_; d++) {
    for (size_t k = 0; k < n_segments_ - 1; k++) {
      size_t row = d * (2 * N / 2 + (n_segments_-1) * N) + k * N + N / 2;
      size_t col = dimension_ * 2 * N / 2 + d * (n_segments_-1) * N / 2 + k * N / 2;
      for (size_t i = 0; i < N / 2; i++) {
        reordering_list_kDim.emplace_back(Triplet(row, col, 1.0));
        reordering_list_kDim.emplace_back(Triplet(row + N / 2, col, 1.0));
        row++;
        col++;
      }
    }
  }

  constraint_reordering_kDim_.setFromTriplets(reordering_list_kDim.begin(),
                                         reordering_list_kDim.end());
  for (size_t k = 0; k < 2; k++) {
    size_t row = k * (N / 2 + (n_segments_-1) * N);
    size_t col = k * N / 2;
    for (size_t i = 0; i< N / 2; i++) {
      reordering_list.emplace_back(Triplet(row, col, 1.0));
      row++;
      col++;
    }
  }

  for (size_t k = 0; k < n_segments_ - 1; k++) {
    size_t row = k * N + N / 2;
    size_t col = 2 * N / 2 + k * N / 2;
    for (size_t i = 0; i < N / 2; i++) {
      reordering_list.emplace_back(Triplet(row, col, 1.0));
      reordering_list.emplace_back(Triplet(row + N / 2, col, 1.0));
      row++;
      col++;
    }
  }

  constraint_reordering_.setFromTriplets(reordering_list.begin(),
                                         reordering_list.end());
}


template <int _N>
bool PolynomialOptimizationConstrained<_N>::setupFromVertices(
        const Vertex::Vector& vertices, const std::vector<double>& times, const std::vector<std::pair<double, double>>& radii,
        int derivative_to_optimize) {
  CHECK(derivative_to_optimize >= 0 &&
        derivative_to_optimize <= kHighestDerivativeToOptimize)
  << "You tried to optimize the " << derivative_to_optimize
  << "th derivative of position on a " << N
  << "th order polynomial. This is not possible, you either need a higher "
     "order polynomial or a smaller derivative to optimize.";

  derivative_to_optimize_ = derivative_to_optimize;
  vertices_ = vertices;
  segment_times_ = times;
  segment_radii_ = radii;

  n_vertices_ = vertices.size();
  n_segments_ = n_vertices_ - 1;

  segments_.resize(n_segments_, Segment(N, dimension_));

  CHECK(n_vertices_ == times.size() + 1)
  << "Size of times must be one less than positions.";

  CHECK(n_vertices_ == radii.size() + 1)
  << "Size of radii must be one less than positions.";

  inverse_mapping_matrices_.resize(n_segments_);
  cost_matrices_.resize(n_segments_);
  inverse_control_pt_mapping_matrices_.resize(n_segments_);

  for (size_t i = 0; i < n_segments_; ++i) {
    const double segment_time = segment_times_[i];
    typename PolynomialOptimization<_N>::SquareMatrix inverse_control_pt_mapping_matrix;
    setupInverseControlPointMappingMatrix(&inverse_control_pt_mapping_matrix, segment_time);
    inverse_control_pt_mapping_matrices_[i] = inverse_control_pt_mapping_matrix;
  }
  // Iterate through all vertices and remove invalid constraints (order too
  // high).

  for (size_t vertex_idx = 0; vertex_idx < n_vertices_; ++vertex_idx) {
    Vertex& vertex = vertices_[vertex_idx];

    // Check if we have valid constraints.
    bool vertex_valid = true;
    Vertex vertex_tmp(dimension_);
    for (Vertex::Constraints::const_iterator it = vertex.cBegin();
         it != vertex.cEnd(); ++it) {
      if (it->first > kHighestDerivativeToOptimize) {
        vertex_valid = false;
        LOG(WARNING) << "Invalid constraint on vertex " << vertex_idx
                     << ": maximum possible derivative is "
                     << kHighestDerivativeToOptimize << ", but was set to "
                     << it->first << ". Ignoring constraint";
      } else {
        vertex_tmp.addConstraint(it->first, it->second);
      }
    }
    if (!vertex_valid) {
      vertex = vertex_tmp;
    }
  }
  this->updateSegmentTimes(times);
  setupConstraintReorderingMatrixkDim();
  return true;
}

template <int _N>
void PolynomialOptimizationConstrained<_N>::constructRkDim(
        Eigen::SparseMatrix<double>* R_kDim) const {

  CHECK_NOTNULL(R_kDim);
  typedef Eigen::Triplet<double> Triplet;
  std::vector<Triplet> cost_unconstrained_triplets;
  cost_unconstrained_triplets.reserve(dimension_* N * N * n_segments_);

  for (size_t i = 0; i < n_segments_; ++i) {
    const typename PolynomialOptimization<_N>::SquareMatrix& Ai = inverse_mapping_matrices_[i];
    const typename PolynomialOptimization<_N>::SquareMatrix& Q = cost_matrices_[i];
    const typename PolynomialOptimization<_N>::SquareMatrix H = Ai.transpose() * Q * Ai;

    for (size_t d = 0; d < dimension_; d++){
      const int start_pos = i * N + d * N * n_segments_;
      for (int row = 0; row < N; ++row) {
        for (int col = 0; col < N; ++col) {
          cost_unconstrained_triplets.emplace_back(
                  Triplet(start_pos + row, start_pos + col, H(row, col)));
        }
      }
    }
  }
  Eigen::SparseMatrix<double> cost_unconstrained(N * n_segments_ * dimension_,
          N * n_segments_ * dimension_);
  cost_unconstrained.setFromTriplets(cost_unconstrained_triplets.begin(),
                                     cost_unconstrained_triplets.end());

  // [1]: R = C^T * H * C. C: constraint_reodering_ ; H: cost_unconstrained,
  // assembled from the block-H above.
  *R_kDim = constraint_reordering_kDim_.transpose() * cost_unconstrained *
       constraint_reordering_kDim_;
}

template <int _N>
bool PolynomialOptimizationConstrained<_N>::solveLinear() {
  CHECK(derivative_to_optimize_ >= 0 &&
        derivative_to_optimize_ <= kHighestDerivativeToOptimize);
  // Catch the fully constrained case:
  if (n_free_constraints_kDim_ == 0) {
    LOG(WARNING)
            << "No free constraints set in the vertices. Polynomial can "
               "not be optimized. Outputting fully constrained polynomial.";
    //updateSegmentsFromCompactConstraints();
    return true;
  }

  // Compute cost matrix for the unconstrained optimization problem.
  // Block-wise H = A^{-T}QA^{-1} according to [1]
  Eigen::SparseMatrix<double> R_kDim;
  constructRkDim(&R_kDim);

  setupControlPointConstraints();

  // Extract block matrices and prepare solver.
  Eigen::SparseMatrix<double> R_kDim_pf = R_kDim.block(
          n_fixed_constraints_kDim_, 0, n_free_constraints_kDim_, n_fixed_constraints_kDim_);
  Eigen::SparseMatrix<double> R_kDim_pp =
          R_kDim.block(n_fixed_constraints_kDim_, n_fixed_constraints_kDim_, n_free_constraints_kDim_,
                  n_free_constraints_kDim_);
  Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
          solver;
  solver.compute(R_kDim_pp);

  // Compute dp_opt for every dimension.
  Eigen::VectorXd df =
          -R_kDim_pf * fixed_constraints_compact_kDim_;  // Rpf = Rfp^T
  free_constraints_compact_kDim_ = solver.solve(df);  // dp = -Rpp^-1 * Rpf * df
  free_constraints_compact_kDim_ = df;

  free_constraints_compact_[0] = free_constraints_compact_kDim_.topRows(n_free_constraints_);
  free_constraints_compact_[1] = free_constraints_compact_kDim_.middleRows(n_free_constraints_,n_free_constraints_);
  free_constraints_compact_[2] = free_constraints_compact_kDim_.bottomRows(n_free_constraints_);

  this->updateSegmentsFromCompactConstraints();
  return true;
}

template <int _N>
void PolynomialOptimizationConstrained<_N>::setupInverseControlPointMappingMatrix(typename PolynomialOptimization<_N>::SquareMatrix* B_inv, double T) {
  Eigen::SparseMatrix<double> control_point_mapping_matrix(N * n_segments_ * dimension_,
          N * n_segments_ * dimension_);
  // [x 0 0 0 0 0]
  // [x x 0 0 0 0]
  // [x x x 0 0 0]
  // [0 0 0 0 0 x]
  // [0 0 0 0 x x]
  // [0 0 0 x x x]
  // ==>
  // [ B_ul   B_ur=0]
  // [ B_ll=0 B_lr  ]
  // We make use of the sparsity of the matrix, so the inverse is:
  // [ inv(B_ul) B_ur=0   ]
  // [ B_ll=0    inv(B_lr)]

  Eigen::MatrixXd control_point_mapping_coefficients(N / 2, N / 2);

  control_point_mapping_coefficients.setZero();
  control_point_mapping_coefficients(0,0) = 1;

  int n = N - 1;
  for (int l = 1; l < N / 2; l++) {
    for (int j = 0; j < N / 2; j++) {
      if (j <= l) {
        int deriv = l;
        control_point_mapping_coefficients(l,j) = factorial(n) / factorial(n - deriv) * pow(-1,(deriv + j)) / pow(T,deriv) * binomialCoeff(deriv, j);
      }
    }
  }

  Eigen::Matrix<double, N / 2, N / 2> B_ul_inv = control_point_mapping_coefficients.inverse();

  for (int k = 0; k < N / 2; ++k) {
    for (int i = 0; i < N / 2; ++i) {
      if (B_ul_inv(k,i) > -0.00001 && B_ul_inv(k,i) < 0.00001) {
        B_ul_inv(k,i) = 0;
      }
    }
  }

  Eigen::VectorXd alternating_sign(N / 2,1);
  for (int i = 0; i < N / 2; i++) {
    alternating_sign(i) = pow(-1, i);
  }

  Eigen::Matrix<double, N / 2, N / 2> B_lr_inv = B_ul_inv.colwise().reverse() * alternating_sign.asDiagonal();

  B_inv->setZero();
  B_inv->topLeftCorner(N / 2, N / 2) = B_ul_inv;
  B_inv->bottomRightCorner(N / 2, N / 2) = B_lr_inv;
}

template <int _N>
void PolynomialOptimizationConstrained<_N>::setupControlPointConstraints() {
  constr_quad_.clear();
  constr_lin_.clear();
  constr_const_.clear();
  Eigen::MatrixXd F = Eigen::MatrixXd::Zero(dimension_, n_all_constraints_kDim_);

  Eigen::MatrixXd B_inv = Eigen::MatrixXd::Zero(dimension_ * n_segments_ * N, dimension_ * n_segments_ * N);

  for (int i = 0; i < n_segments_; ++i) {
    for (int d = 0; d < dimension_; ++d) {
      int start_pos = i * N + d * N * n_segments_;
      B_inv.block<N, N>(start_pos, start_pos) = inverse_control_pt_mapping_matrices_[i];
    }
  }

  for (int i = 0; i < n_segments_; i++) {
    std::vector<Eigen::MatrixXd> control_pt_extraction_matrices;
    //std::cout << "============================" << std::endl;
    for (int j = 0; j < N; ++j) {
      F.setZero();
      for (int k = 0; k < dimension_; ++k) {
        F(k, (k * N * n_segments_ + i * N + j)) = 1;
      }
      control_pt_extraction_matrices.push_back(F * B_inv * constraint_reordering_kDim_);
    }

    if (i < n_segments_ - 1) {
      compute_sphere_constraints(control_pt_extraction_matrices, i);
    }

    compute_tube_constraints(control_pt_extraction_matrices, i);
    compute_tube_end_constraints(control_pt_extraction_matrices, i);
  }
}

template <int _N>
void PolynomialOptimizationConstrained<_N>::compute_sphere_constraints(std::vector<Eigen::MatrixXd> &control_pt_extraction_matrices, int segment_idx) {
  Eigen::MatrixXd control_pt_extraction_matrix = control_pt_extraction_matrices.back();
  Eigen::VectorXd pos_vec;
  vertices_[segment_idx+1].getConstraint(0, &pos_vec);
  constr_const_.push_back(pos_vec.transpose() * pos_vec - pow(segment_radii_[segment_idx].second,2));
  constr_lin_.push_back(-2*pos_vec.transpose() * control_pt_extraction_matrix.rightCols(n_free_constraints_kDim_));
  constr_quad_.push_back(2*control_pt_extraction_matrix.rightCols(n_free_constraints_kDim_).transpose()*control_pt_extraction_matrix.rightCols(n_free_constraints_kDim_));
};



template <int _N>
void PolynomialOptimizationConstrained<_N>::compute_tube_constraints(std::vector<Eigen::MatrixXd> &control_pt_extraction_matrices, int segment_idx) {
  Eigen::VectorXd segment_start;
  Eigen::VectorXd segment_end;

  vertices_[segment_idx].getConstraint(0, &segment_start);
  vertices_[segment_idx + 1].getConstraint(0, &segment_end);

  Eigen::VectorXd segment_vec = segment_end - segment_start;
  segment_vec.normalize();

  double nx = segment_vec(0);
  double ny = segment_vec(1);
  double nz = segment_vec(2);

  double px = segment_start(0);
  double py = segment_start(1);
  double pz = segment_start(2);

  Eigen::Matrix3d A;
  Eigen::Vector3d b;
  A  << 1 - pow(nx, 2), -nx * ny, -nx * nz,
        -nx * ny, 1 - pow(ny, 2), -ny * nz,
        -nx * nz, -ny * nz, 1 - pow(nz, 2);

  for (int i = 0; i < dimension_; ++i) {
    for (int j = 0; j < dimension_; ++j) {
      if (A(i,j) > -0.000001 && A(i,j) < 0.000001) {
        A(i,j) = 0;
      }
    }
  }

  b <<  (pow(nx, 2) - 1) * px + nx * ny * py + nx * nz * pz,
        nx * ny * px + (pow(ny, 2) - 1) * py + ny * nz * pz,
        nx * nz * px + ny * nz * py + (pow(nz, 2) - 1) * pz;

  for (int j = 0; j < dimension_; ++j) {
    if (b(j) > -0.000001 && b(j) < 0.000001) {
      b(j) = 0;
    }
  }

  Eigen::Matrix3d LL = A.transpose()*A;
  Eigen::Vector3d L  = 2*b.transpose()*A;

  double mu = b.transpose()*b-pow(segment_radii_[segment_idx].first,2);

  for (int i = 1; i < N - 1; ++i) {
    Eigen::MatrixXd EE = (control_pt_extraction_matrices[i]).transpose() * LL * control_pt_extraction_matrices[i];
    Eigen::MatrixXd E  = L.transpose() * control_pt_extraction_matrices[i];

    constr_const_.push_back((fixed_constraints_compact_kDim_.transpose() * EE.topLeftCorner(n_fixed_constraints_kDim_, n_fixed_constraints_kDim_) *
                               fixed_constraints_compact_kDim_ + E.leftCols(n_fixed_constraints_kDim_) * fixed_constraints_compact_kDim_)(0) + mu);

    constr_lin_.push_back(2*(fixed_constraints_compact_kDim_.transpose() * EE.topRightCorner(n_fixed_constraints_kDim_, n_free_constraints_kDim_)) +
                           E.rightCols(n_free_constraints_kDim_));

    constr_quad_.push_back(2 * EE.bottomRightCorner(n_free_constraints_kDim_, n_free_constraints_kDim_));
  }
}

template <int _N>
void PolynomialOptimizationConstrained<_N>::compute_tube_end_constraints(std::vector<Eigen::MatrixXd> &control_pt_extraction_matrices, int segment_idx) {

  int n_mid_control_pts = N-2;

  Eigen::MatrixXd blk_0 = Eigen::MatrixXd::Zero(n_free_constraints_kDim_, n_free_constraints_kDim_);

  for (int k = 0; k < n_mid_control_pts; ++k) {
    Eigen::VectorXd segment_start;
    Eigen::VectorXd segment_end;

    vertices_[segment_idx].getConstraint(0, &segment_start);
    vertices_[segment_idx + 1].getConstraint(0, &segment_end);

    Eigen::Vector3d segment_vec = segment_end - segment_start;
    segment_vec.normalize();

    Eigen::Vector3d norm_vec_start = -segment_vec;
    Eigen::Vector3d norm_vec_end   = segment_vec;

    Eigen::Vector3d p_start;
    if (segment_idx == 0)
      p_start = segment_start + norm_vec_start * segment_radii_[0].first;
    else
      p_start = segment_start + norm_vec_start * segment_radii_[segment_idx-1].second;

    Eigen::Vector3d p_end = segment_end + norm_vec_end * segment_radii_[segment_idx].second;

    std::vector<Eigen::MatrixXd> side_lin;
    std::vector<Eigen::MatrixXd> side_con;

    side_lin.push_back(norm_vec_start.transpose() * control_pt_extraction_matrices[k+1]);
    side_con.push_back(-norm_vec_start.transpose() * p_start);

    side_lin.push_back(norm_vec_end.transpose() * control_pt_extraction_matrices[k+1]);
    side_con.push_back(-norm_vec_end.transpose() * p_end);

    for (int i = 0; i < 2; ++i) {
      constr_const_.push_back((side_con[i] + side_lin[i].leftCols(n_fixed_constraints_kDim_) * fixed_constraints_compact_kDim_)(0));
      constr_lin_.push_back(side_lin[i].rightCols(n_free_constraints_kDim_));
      constr_quad_.push_back(blk_0);
    }
  }
}

template <int _N>
int PolynomialOptimizationConstrained<_N>::solveQCQP() {
  CHECK(derivative_to_optimize_ >= 0 &&
        derivative_to_optimize_ <= kHighestDerivativeToOptimize);
  // Catch the fully constrained case:
  if (n_free_constraints_kDim_ == 0) {
    LOG(WARNING)
            << "No free constraints set in the vertices. Polynomial can "
               "not be optimized. Outputting fully constrained polynomial.";
    this->updateSegmentsFromCompactConstraints();
    return true;
  }

  // Compute cost matrix for the unconstrained optimization problem.
  // Block-wise H = A^{-T}QA^{-1} according to [1]
  Eigen::SparseMatrix<double> R_kDim;
  constructRkDim(&R_kDim);

  setupControlPointConstraints();

  // Extract block matrices and prepare solver.
  Eigen::SparseMatrix<double> R_kDim_fp = R_kDim.topRightCorner(
          n_fixed_constraints_kDim_, n_free_constraints_kDim_);
  Eigen::SparseMatrix<double> R_kDim_pp =
          R_kDim.bottomRightCorner(n_free_constraints_kDim_, n_free_constraints_kDim_);

  Eigen::MatrixXd f = 2 * fixed_constraints_compact_kDim_.transpose() * R_kDim_fp;

  int NUMCON = constr_const_.size();
  int NUMVAR = n_free_constraints_kDim_;
  int NUMQNZ = (R_kDim_pp.diagonal().nonZeros() + R_kDim_pp.nonZeros())/2;

  MSKrescodee  r;

  MSKboundkeye bkc  = MSK_BK_UP;

  double       blc  = -MSK_INFINITY;

  MSKboundkeye  bkx = MSK_BK_FR;

  double        blx = -MSK_INFINITY;

  double        bux = +MSK_INFINITY;

  Eigen::MatrixXd constr_lin_matrix(NUMCON, n_free_constraints_kDim_);

  for (int constr_idx = 0; constr_idx < NUMCON; ++constr_idx) {
    constr_lin_matrix.row(constr_idx) = constr_lin_[constr_idx];
  }

  MSKint32t aptrb[n_free_constraints_kDim_];
  MSKint32t aptre[n_free_constraints_kDim_];
  std::vector<int> asub_temp;
  std::vector<double> aval_temp;

  for (int variable_idx = 0; variable_idx < n_free_constraints_kDim_; ++variable_idx) {
    int n_non_zero = 0;
    Eigen::VectorXd constr_lin_matrix_col = constr_lin_matrix.col(variable_idx);

    //std::cout << "constr_lin_matrix_col " << constr_lin_matrix_col << std::endl;

    if (variable_idx == 0) {
      aptrb[variable_idx] = 0;
    } else {
      aptrb[variable_idx] = aptre[variable_idx-1];
    }

    for (int col_idx = 0; col_idx < NUMCON; ++col_idx) {
      if (constr_lin_matrix_col(col_idx) != 0) {
        n_non_zero++;
        asub_temp.push_back(col_idx);
        aval_temp.push_back(constr_lin_matrix_col(col_idx));
      }
    }

    aptre[variable_idx] = aptrb[variable_idx] + n_non_zero;

  }

  MSKint32t   asub[asub_temp.size()];

  for (int asub_idx = 0; asub_idx < asub_temp.size(); ++asub_idx) {
    asub[asub_idx] = asub_temp[asub_idx];
  }

  double      aval[aval_temp.size()];

  for (int aval_idx = 0; aval_idx < aval_temp.size(); ++aval_idx) {
    aval[aval_idx] = aval_temp[aval_idx];
  }

  MSKint32t   qsubi[NUMQNZ],
              qsubj[NUMQNZ];
  double      qval[NUMQNZ];

  MSKint32t   j, i;
  double      xx[NUMVAR];
  MSKenv_t    env;
  MSKtask_t   task;

  /* Create the mosek environment. */
  r = MSK_makeenv(&env, NULL);

  if ( r == MSK_RES_OK )
  {
    /* Create the optimization task. */
    r = MSK_maketask(env, NUMCON, NUMVAR, &task);

    if ( r == MSK_RES_OK )
    {
      r = MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr);

      /* Append 'NUMCON' empty constraints.
       The constraints will initially have no bounds. */
      if ( r == MSK_RES_OK )
        r = MSK_appendcons(task, NUMCON);

      /* Append 'NUMVAR' variables.
       The variables will initially be fixed at zero (x=0). */
      if ( r == MSK_RES_OK )
        r = MSK_appendvars(task, NUMVAR);

      /* Optionally add a constant term to the objective. */
      if ( r == MSK_RES_OK )
        r = MSK_putcfix(task, 0.0);
      for (j = 0; j < NUMVAR && r == MSK_RES_OK; ++j)
      {
        /* Set the linear term c_j in the objective.*/
        if (r == MSK_RES_OK)
          r = MSK_putcj(task, j, f(j));

        /* Set the bounds on variable j.
        blx[j] <= x_j <= bux[j] */
        if (r == MSK_RES_OK)
          r = MSK_putvarbound(task,
                              j,           /* Index of variable.*/
                              bkx,      /* Bound key.*/
                              blx,      /* Numerical value of lower bound.*/
                              bux);     /* Numerical value of upper bound.*/

        /* Input column j of A */
        auto test_a = aptre[j] - aptrb[j];
        auto test_b = *(asub + aptrb[j]);
        auto test_c = *(aval + aptrb[j]);

        if (r == MSK_RES_OK)
          r = MSK_putacol(task,
                          j,                 /* Variable (column) index.*/
                          aptre[j] - aptrb[j], /* Number of non-zeros in column j.*/
                          asub + aptrb[j],   /* Pointer to row indexes of column j.*/
                          aval + aptrb[j]);  /* Pointer to Values of column j.*/

      }

      /* Set the bounds on constraints.
         for i=1, ...,NUMCON : blc[i] <= constraint i <= buc[i] */
      for (i = 0; i < NUMCON && r == MSK_RES_OK; ++i) {
        r = MSK_putconbound(task,
                            i,           /* Index of constraint.*/
                            bkc,      /* Bound key.*/
                            blc,      /* Numerical value of lower bound.*/
                            -constr_const_[i]);     /* Numerical value of upper bound.*/
      }

      if ( r == MSK_RES_OK )
      {
        /*
         * The lower triangular part of the Q^o
         * matrix in the objective is specified.
         */
        int q_idx = 0;
        for (int row = 0; row < (n_free_constraints_) ; ++row) {
          for (int col = 0; col < (n_free_constraints_); ++col) {
            if (col <= row) {
              double cost_value = R_kDim_pp.coeff(row,col);
              if (cost_value != 0) {
                for (int d = 0; d < dimension_; ++d) {

                  qsubi[d * NUMQNZ / 3 + q_idx] = row + d * n_free_constraints_;
                  qsubj[d * NUMQNZ / 3 + q_idx] = col + d * n_free_constraints_;
                  qval[d * NUMQNZ / 3 + q_idx] = 2*cost_value;

                }
                q_idx++;
              }
            }
          }
        }

        r = MSK_putqobj(task, NUMQNZ, qsubi, qsubj, qval);
      }

      if ( r == MSK_RES_OK )
      {
        for (int k = 0; k < NUMCON; ++k) {

          Eigen::MatrixXd constr_quad_temp = constr_quad_[k];

          int q_idx = 0;
          for (int row = 0; row < n_free_constraints_kDim_; ++row) {
            for (int col = 0; col < n_free_constraints_kDim_; ++col) {
              if (col <= row) {
                double cost_value = constr_quad_temp(row,col);
                if (cost_value != 0) {
                  qsubi[q_idx] = row;
                  qsubj[q_idx] = col;
                  qval[q_idx] = cost_value;
                  q_idx++;
                }
              }
            }
          }

          int n_constr = q_idx;

          r = MSK_putqconk(task,
                           k,
                           n_constr,
                           qsubi,
                           qsubj,
                           qval
                          );
        }
      }

      if ( r == MSK_RES_OK )
        r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);

      if ( r == MSK_RES_OK )
      {
        MSKrescodee trmcode;

        /* Run optimizer */
        r = MSK_optimizetrm(task, &trmcode);

        /* Print a summary containing information
           about the solution for debugging purposes*/
        MSK_solutionsummary (task, MSK_STREAM_LOG);

        if ( r == MSK_RES_OK )
        {
          MSKsolstae solsta;
          int j;
          MSK_getsolsta (task, MSK_SOL_ITR, &solsta);
          switch (solsta)
          {
            case MSK_SOL_STA_OPTIMAL:
            case MSK_SOL_STA_NEAR_OPTIMAL:
              MSK_getxx(task,
                        MSK_SOL_ITR,    /* Request the interior solution. */
                        xx);
              std::cout << "d_p = [";
              for (j = 0; j < NUMVAR; ++j)
                if (j < NUMVAR-1) {
                  std::cout << xx[j] << "; ";
                } else {
                  std::cout << xx[j];
                }
              std::cout << "];" << "\n";
              break;

            case MSK_SOL_STA_DUAL_INFEAS_CER:
            case MSK_SOL_STA_PRIM_INFEAS_CER:
            case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
            case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
              printf("Primal or dual infeasibility certificate found.\n");
              break;

            case MSK_SOL_STA_UNKNOWN:
              printf("The status of the solution could not be determined.\n");
              break;

            default:
              printf("Other solution status.");
              break;
          }
        }
        else
        {
          printf("Error while optimizing.\n");
        }
      }

      if (r != MSK_RES_OK)
      {
        /* In case of an error print error code and description. */
        char symname[MSK_MAX_STR_LEN];
        char desc[MSK_MAX_STR_LEN];

        printf("An error occurred while optimizing.\n");
        MSK_getcodedesc (r,
                         symname,
                         desc);
        printf("Error %s - '%s'\n", symname, desc);
      }
    }

    MSK_deletetask(&task);
  }
  MSK_deleteenv(&env);

  for (int i = 0; i < NUMVAR; i++) {
    free_constraints_compact_kDim_(i) = xx[i];
  }

  free_constraints_compact_[0] = free_constraints_compact_kDim_.topRows(n_free_constraints_);
  free_constraints_compact_[1] = free_constraints_compact_kDim_.middleRows(n_free_constraints_,n_free_constraints_);
  free_constraints_compact_[2] = free_constraints_compact_kDim_.bottomRows(n_free_constraints_);

  this->updateSegmentsFromCompactConstraints();

  return ( r );
}


template <int _N>
int PolynomialOptimizationConstrained<_N>::factorial(int n) {
  if(n > 1)
    return n * factorial(n - 1);
  else
    return 1;
}

template <int _N>
int PolynomialOptimizationConstrained<_N>::binomialCoeff(int n, int k) {
  int res = 1;
  // Since C(n, k) = C(n, n-k)
  if ( k > n - k )
    k = n - k;

  // Calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
  for (int i = 0; i < k; ++i)
  {
    res *= (n - i);
    res /= (i + 1);
  }

  return res;
}

//#include <impl/polynomial_optimization_linear_constrained_impl.h>

} //namespace mav_trajectory_generation

#endif //MAV_TUBE_TRAJECTORY_GENERATION_POLYNOMIAL_OPTIMIZATION_QCQP_H
