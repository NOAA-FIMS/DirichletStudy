#ifndef DIRICHLET_COMPOSITION_HPP
#define DIRICHLET_COMPOSITION_HPP

#include "ATL/lib/Optimization.hpp"
#include <random>

template <typename T>
class DirichletComposition : public atl::FunctionMinimizer<T>
{
public:
    enum class DirichletCompositionType
    {
        DEFAULT = 0,
        LINEAR,
        SATURATED
    };

    DirichletCompositionType type = DirichletCompositionType::DEFAULT;

    DirichletComposition()
    {
    }

    virtual void Initialize()
    {
    }

    virtual void Evaluate(atl::Variable<T> &x)
    {
    }

private:
};

#endif // DIRICHLET_COMPOSITION_HPP