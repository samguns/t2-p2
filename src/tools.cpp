#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse.setZero();

  int vector_size = estimations.size();
  if (vector_size == 0 ||
      (estimations.size() != ground_truth.size())) {
    return rmse;
  }

  for (int i = 0; i < vector_size; i++) {
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  rmse = rmse / vector_size;
  rmse = rmse.array().sqrt();

  return rmse;
}