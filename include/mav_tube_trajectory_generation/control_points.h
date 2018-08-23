//
// Created by nilsiism on 22/08/18.
//

#ifndef MAV_TUBE_TRAJECTORY_GENERATION_CONTROL_POINTS_H
#define MAV_TUBE_TRAJECTORY_GENERATION_CONTROL_POINTS_H

class ControlPoints {

  static void controlPointMappingCoeffs(int N, int derivative,
                                 Eigen::VectorXd* coeffs) {
    CHECK_LT(derivative, N);
    CHECK_GE(derivative, 0);

    coeffs->resize(N, 1);
    coeffs->setZero();
    // first coefficient doesn't get multiplied
    (*coeffs)[derivative] = base_coefficients_(derivative, derivative);

    if (std::abs(t) < std::numeric_limits<double>::epsilon()) return;

    // now multiply increasing power of t towards the right
    for (int j = derivative + 1; j < N; j++) {
      (*coeffs)[j] = base_coefficients_(derivative, j) * t_power;
      t_power = t_power * t;
    }
  }

  int factorial(int n) {
    if(n > 1)
      return n * factorial(n - 1);
    else
      return 1;
  }

  int binomialCoeff(int n, int k) {
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



};

#endif //MAV_TUBE_TRAJECTORY_GENERATION_CONTROL_POINTS_H
