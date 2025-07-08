#if !defined(DIRICHLET_DEFAULT_HPP)
#define DIRICHLET_DEFAULT_HPP

template <typename Type>
class DirichletDefault : public FunctorBase<Type>
{
public:
    virtual void Initialize() override
    {
    }

    virtual void Evaluate() override
    {
    }
};

#endif // DIRICHLET_DEFAULT_HPP