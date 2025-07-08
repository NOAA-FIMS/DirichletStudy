#ifndef DIRICHLET_LINEAR_HPP
#define DIRICHLET_LINEAR_HPP

#include "functor_base.hpp"

template <typename Type>
class DirichletLinear : public FunctorBase<Type>
{
public:
    virtual void Initialize() override
    {
    }

    virtual void Evaluate() override
    {
    }
};

#endif // DIRICHLET_LINEAR_HPP
