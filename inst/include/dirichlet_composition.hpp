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

    std::shared_ptr<FunctorBase<T> > functor;

    DirichletCompositionType type = DirichletCompositionType::DEFAULT;

    DirichletComposition()
    {
    }

    void setFunctor(std::shared_ptr<FunctorBase<T> > functor)
    {
        this->functor = functor;
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