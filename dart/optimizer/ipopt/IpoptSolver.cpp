/*
 * Copyright (c) 2014, Georgia Tech Research Corporation
 * All rights reserved.
 *
 * Author(s): Jeongseok Lee <jslee02@gmail.com>
 *
 * Georgia Tech Graphics Lab and Humanoid Robotics Lab
 *
 * Directed by Prof. C. Karen Liu and Prof. Mike Stilman
 * <karenliu@cc.gatech.edu> <mstilman@cc.gatech.edu>
 *
 * This file is provided under the following "BSD-style" License:
 *   Redistribution and use in source and binary forms, with or
 *   without modification, are permitted provided that the following
 *   conditions are met:
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 *   CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 *   INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *   MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 *   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 *   USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *   AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *   LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *   POSSIBILITY OF SUCH DAMAGE.
 */

#include "dart/optimizer/ipopt/IpoptSolver.h"

#include "dart/common/Console.h"
#include "dart/optimizer/Function.h"
#include "dart/optimizer/Problem.h"

namespace dart {
namespace optimizer {

//==============================================================================
IpoptSolver::IpoptSolver(Problem* _problem)
  : Solver(_problem)
{
  assert(_problem);

  // Create a new instance of nlp (use a SmartPtr, not raw)
  mNlp = new DartTNLP(_problem);

  // Create a new instance of IpoptApplication (use a SmartPtr, not raw). We are
  // using the factory, since this allows us to compile this with an Ipopt
  // Windows DLL.
  mIpoptApp = IpoptApplicationFactory();

  // Change some options
  // Note: The following choices are only examples, they might not be
  //       suitable for your optimization problem.
  mIpoptApp->Options()->SetNumericValue("tol", 1e-9);
  mIpoptApp->Options()->SetStringValue("mu_strategy", "adaptive");
  mIpoptApp->Options()->SetStringValue("output_file", "ipopt.out");
  mIpoptApp->Options()->SetStringValue("hessian_approximation", "limited-memory");

  // Intialize the IpoptApplication and process the options
  Ipopt::ApplicationReturnStatus status = mIpoptApp->Initialize();
  if (status != Ipopt::Solve_Succeeded)
    dterr << "Error during ipopt initialization.\n";
}

//==============================================================================
IpoptSolver::~IpoptSolver()
{
}

//==============================================================================
bool IpoptSolver::solve()
{
  // Ask Ipopt to solve the problem
  Ipopt::ApplicationReturnStatus status = mIpoptApp->OptimizeTNLP(mNlp);

  if (status == Ipopt::Solve_Succeeded)
    return true;
  else
    return false;
}

//==============================================================================
DartTNLP::DartTNLP(Problem* _problem)
  : Ipopt::TNLP(),
    mProblem(_problem)
{
  assert(_problem && "Null pointer is not allowed.");
}

//==============================================================================
DartTNLP::~DartTNLP()
{

}

//==============================================================================
bool DartTNLP::get_nlp_info(Ipopt::Index& n,
                            Ipopt::Index& m,
                            Ipopt::Index& nnz_jac_g,
                            Ipopt::Index& nnz_h_lag,
                            Ipopt::TNLP::IndexStyleEnum& index_style)
{
  // The problem described in HS071_NLP.hpp has 4 variables, x[0] through x[3]
  n = mProblem->getDimension();

  // one equality constraint and one inequality constraint
  m = mProblem->getNumEqConstraints() + mProblem->getNumIneqConstraints();

  // in this example the Jacobian is dense and contains 8 nonzeros
  nnz_jac_g = n * m;

  // the Hessian is also dense and has 16 total nonzeros, but we
  // only need the lower left corner (since it is symmetric)
  nnz_h_lag = n * n * m;

  // use the C style indexing (0-based)
  index_style = Ipopt::TNLP::C_STYLE;

  return true;
}

//==============================================================================
bool DartTNLP::get_bounds_info(Ipopt::Index n,
                               Ipopt::Number* x_l,
                               Ipopt::Number* x_u,
                               Ipopt::Index m,
                               Ipopt::Number* g_l,
                               Ipopt::Number* g_u)
{
  // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
  // If desired, we could assert to make sure they are what we think they are.
  assert(n == mProblem->getDimension());
  assert(m == mProblem->getNumEqConstraints()
              + mProblem->getNumIneqConstraints());

  // lower and upper bounds
  for (Ipopt::Index i = 0; i < n; i++)
  {
    x_l[i] = mProblem->getLowerBounds()[i];
    x_u[i] = mProblem->getUpperBounds()[i];
  }

  // Add inequality constraint functions
  size_t idx = 0;
  for (size_t i = 0; i < mProblem->getNumEqConstraints(); ++i)
  {
    g_l[idx] = g_u[idx] = 0.0;
    idx++;
  }

  for (size_t i = 0; i < mProblem->getNumIneqConstraints(); ++i)
  {
    // Ipopt interprets any number greater than nlp_upper_bound_inf as
    // infinity. The default value of nlp_upper_bound_inf and
    // nlp_lower_bound_inf is 1e+19 and can be changed through ipopt options.
    g_l[idx] = -2e+19;
    g_u[idx] = 0.0;
    idx++;
  }

  return true;
}

//==============================================================================
bool DartTNLP::get_starting_point(Ipopt::Index n,
                                  bool init_x,
                                  Ipopt::Number* x,
                                  bool init_z,
                                  Ipopt::Number* /*z_L*/,
                                  Ipopt::Number* /*z_U*/,
                                  Ipopt::Index /*m*/,
                                  bool init_lambda,
                                  Ipopt::Number* /*lambda*/)
{
  // If init_x is true, this method must provide an initial value for x.
  if (init_x)
  {
    for (size_t i = 0; i < n; ++i)
      x[i] = mProblem->getInitialGuess()[i];
  }

  // If init_z is true, this method must provide an initial value for the bound
  // multipliers z^L and z^U
  if (init_z)
  {
    // TODO(JS): Not implemented yet.
    dterr << "Initializing lower/upper bounds for z is not supported yet. "
          << "Ignored here.\n";
  }

  // If init_lambda is true, this method must provide an initial value for the
  // constraint multipliers, lambda.
  if (init_lambda)
  {
    // TODO(JS): Not implemented yet.
    dterr << "Initializing lambda is not supported yet. "
          << "Ignored here.\n";
  }

  return true;
}

//==============================================================================
bool DartTNLP::eval_f(Ipopt::Index _n,
                      const Ipopt::Number* _x,
                      bool _new_x,
                      Ipopt::Number& _obj_value)
{
  if (_new_x)
  {
    Eigen::Map<const Eigen::VectorXd> x(_x, _n);
    mObjValue = mProblem->getObjective()->eval(x);
  }

  _obj_value = mObjValue;

  return true;
}

//==============================================================================
bool DartTNLP::eval_grad_f(Ipopt::Index _n,
                           const Ipopt::Number* _x,
                           bool _new_x,
                           Ipopt::Number* _grad_f)
{
  if (_new_x)
  {
    Eigen::Map<const Eigen::VectorXd> x(_x, _n);
    Eigen::Map<Eigen::VectorXd> grad(_grad_f, _n);
    mProblem->getObjective()->evalGradient(x, grad);
  }

  return true;
}

//==============================================================================
bool DartTNLP::eval_g(Ipopt::Index _n,
                      const Ipopt::Number* _x,
                      bool _new_x,
                      Ipopt::Index _m,
                      Ipopt::Number* _g)
{
  assert(_m == mProblem->getNumEqConstraints()
               + mProblem->getNumIneqConstraints());

  // TODO(JS):
  if (_new_x)
  {
  }

  Eigen::Map<const Eigen::VectorXd> x(_x, _n);
  size_t idx = 0;

  // Evaluate function values for equality constraints
  for (size_t i = 0; i < mProblem->getNumEqConstraints(); ++i)
  {
    _g[idx] = mProblem->getEqConstraint(i)->eval(x);
    idx++;
  }

  // Evaluate function values for inequality constraints
  for (size_t i = 0; i < mProblem->getNumIneqConstraints(); ++i)
  {
    _g[idx] = mProblem->getIneqConstraint(i)->eval(x);
    idx++;
  }

  return true;
}

//==============================================================================
bool DartTNLP::eval_jac_g(Ipopt::Index _n,
                          const Ipopt::Number* _x,
                          bool _new_x,
                          Ipopt::Index _m,
                          Ipopt::Index _nele_jac,
                          Ipopt::Index* _iRow,
                          Ipopt::Index* _jCol,
                          Ipopt::Number* _values)
{
  // If the iRow and jCol arguments are not NULL, then IPOPT wants you to fill
  // in the sparsity structure of the Jacobian (the row and column indices
  // only). At this time, the x argument and the values argument will be NULL.

  if (_values == NULL)
  {
    // return the structure of the Jacobian

    // Assume the gradient is dense
    size_t idx = 0;
    for (size_t i = 0; i < _m; ++i)
    {
      for (size_t j = 0; j < _n; ++j)
      {
        _iRow[idx] = i;
        _jCol[idx] = j;
        idx++;
      }
    }
  }
  else
  {
    // return the values of the Jacobian of the constraints
    size_t idx = 0;
    Eigen::Map<const Eigen::VectorXd> x(_x, _n);
    Eigen::Map<Eigen::VectorXd> grad(NULL, 0);

    // Evaluate function values for equality constraints
    for (size_t i = 0; i < mProblem->getNumEqConstraints(); ++i)
    {
      new (&grad)Eigen::Map<Eigen::VectorXd>(_values + idx, _n);
      mProblem->getEqConstraint(i)->evalGradient(x, grad);
      idx += _n;
    }

    // Evaluate function values for inequality constraints
    for (size_t i = 0; i < mProblem->getNumIneqConstraints(); ++i)
    {
      new (&grad)Eigen::Map<Eigen::VectorXd>(_values + idx, _n);
      mProblem->getIneqConstraint(i)->evalGradient(x, grad);
      idx += _n;
    }
  }

  return true;
}

//==============================================================================
bool DartTNLP::eval_h(Ipopt::Index _n,
                      const Ipopt::Number* _x,
                      bool _new_x,
                      Ipopt::Number _obj_factor,
                      Ipopt::Index _m,
                      const Ipopt::Number* _lambda,
                      bool _new_lambda,
                      Ipopt::Index _nele_hess,
                      Ipopt::Index* _iRow,
                      Ipopt::Index* _jCol,
                      Ipopt::Number* _values)
{
  // TODO(JS): Not implemented yet.
  return TNLP::eval_h(_n, _x, _new_x, _obj_factor, _m, _lambda, _new_lambda,
                      _nele_hess, _iRow, _jCol, _values);
}

//==============================================================================
void DartTNLP::finalize_solution(Ipopt::SolverReturn _status,
                                 Ipopt::Index _n,
                                 const Ipopt::Number* _x,
                                 const Ipopt::Number* _z_L,
                                 const Ipopt::Number* _z_U,
                                 Ipopt::Index _m,
                                 const Ipopt::Number* _g,
                                 const Ipopt::Number* _lambda,
                                 Ipopt::Number _obj_value,
                                 const Ipopt::IpoptData* _ip_data,
                                 Ipopt::IpoptCalculatedQuantities* _ip_cq)
{
  // Store optimal and optimum values
  mProblem->setOptimumValue(_obj_value);
  Eigen::VectorXd x = Eigen::VectorXd::Zero(_n);
  for (size_t i = 0; i < _n; ++i)
    x[i] = _x[i];
  mProblem->setOptimalSolution(x);
}

}  // namespace optimizer
}  // namespace dart
