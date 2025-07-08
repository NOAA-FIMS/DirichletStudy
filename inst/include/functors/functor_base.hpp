#ifndef FUNCTOR_BASE_HPP
#define FUNCTOR_BASE_HPP

template <typename T>
class FunctorBase
{
public:
    virtual void Initialize() = 0;
    virtual void Evaluate() = 0;
    virtual ~FunctorBase() = default;

    void operator()()
    {
        this->Evaluate();
    }
};

#endif // FUNCTOR_BASE_HPP