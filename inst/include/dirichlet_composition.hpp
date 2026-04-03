#ifndef DIRICHLET_COMPOSITION_HPP
#define DIRICHLET_COMPOSITION_HPP

#include "ATL/lib/Optimization.hpp"
#include "functors/dirichlet_default.hpp"
#include "functors/dirichlet_linear.hpp"
#include "functors/dirichlet_saturated.hpp"
#include <random>

template <typename T>
class DirichletDefault : public atl::ObjectiveFunction<T> {
public:
  size_t A = 0;
  std::vector<T> x;
  std::vector<T> x_observed;
  std::vector<T> x_predicted;
  std::vector<T> P;
  std::vector<T> P_observed;
  std::vector<T> P_predicted;
  T Nff;

  DirichletCompositionType type = DirichletCompositionType::DEFAULT;

  DirichletDefault() {
  }

  virtual void Initialize() {}

  virtual void Evaluate(atl::Variable<T> &f) {
    f = static_cast<T>(0);
    f = log_dirichlet_multinom(x, P_predicted);
  }

private:
};

template <typename T> class DirichletLinear : public atl::ObjectiveFunction<T> {
public:
  size_t A = 0;
  std::vector<T> x;
  std::vector<T> x_observed;
  std::vector<T> x_predicted;
  std::vector<T> P;
  std::vector<T> P_observed;
  std::vector<T> P_predicted;
  T Nff;

  DirichletCompositionType type = DirichletCompositionType::DEFAULT;

  DirichletLinear() {  }

  virtual void Initialize() {}

  virtual void Evaluate(atl::Variable<T> &f) {
    f = static_cast<T>(0);
    f =log_dirichlet_multinom_linear(x, P_predicted);
  }

private:
};

template <typename T>
class DirichletSaturated : public atl::ObjectiveFunction<T> {
public:
  size_t A = 0;
  std::vector<T> x;
  std::vector<T> x_observed;
  std::vector<T> x_predicted;
  std::vector<T> P;
  std::vector<T> P_observed;
  std::vector<T> P_predicted;
  T Nff;

  DirichletCompositionType type = DirichletCompositionType::DEFAULT;

  DirichletSaturated() {

  }

  virtual void Initialize() {}

  virtual void Evaluate(atl::Variable<T> &f) {
    f = static_cast<T>(0);
    f = log_dirichlet_multinom_saturated(x, P_predicted, Nff);
  }

private:
};

#endif // DIRICHLET_COMPOSITION_HPP
