#ifndef RCPP_INTERFACE_HPP
#define RCPP_INTERFACE_HPP

#include <Rcpp.h>
#include "dirichlet_fa.hpp"


class DirichletStudyInterfaceBase
{
public:
    virtual void setCompositionData(const Rcpp::NumericMatrix &data) = 0;
    virtual void setSimplexData(const Rcpp::NumericMatrix &simplex_data) = 0;
    virtual bool runAnalysis() = 0;
    virtual Rcpp::List getResults() = 0;

    virtual ~DirichletStudyInterfaceBase() {}
};


class DirichletDefaultInterface : public DirichletStudyInterfaceBase
{
public:
    DirichletDefaultInterface()
    {
    }           
    virtual ~DirichletDefaultInterface()
    {
    }   
    void setCompositionData(const Rcpp::NumericMatrix &data) override
    {
        this->data = data;
    }
    void setSimplexData(const Rcpp::NumericMatrix &simplex_data) override
    {
        this->simplex_data = simplex_data;
    }
    bool runAnalysis() override
    {
        // Placeholder for running the analysis
        // This would typically call the Dirichlet_Default class methods.
        return true; // Indicating success
    }
    Rcpp::List getResults() override
    {
        // Placeholder for returning results
        // This would typically return the results of the analysis.
        return Rcpp::List::create(Rcpp::Named("status") = "success");
    }
private:
    Rcpp::NumericMatrix data;
    Rcpp::NumericMatrix simplex_data;
};

class DirichletLinearInterface : public DirichletStudyInterfaceBase
{
public:
    DirichletLinearInterface()
    {
    }       
    virtual ~DirichletLinearInterface()
    {
    }   

    void setCompositionData(const Rcpp::NumericMatrix &data) override
    {
        this->data = data;
    }
    void setSimplexData(const Rcpp::NumericMatrix &simplex_data) override
    {
        this->simplex_data = simplex_data;
    }
    bool runAnalysis() override
    {
        // Placeholder for running the analysis
        // This would typically call the Dirichlet_Linear class methods.
        return true; // Indicating success
    }
    Rcpp::List getResults() override
    {
        // Placeholder for returning results        
        // This would typically return the results of the analysis.
        return Rcpp::List::create(Rcpp::Named("status") = "success");
    }
private:
    Rcpp::NumericMatrix data;
    Rcpp::NumericMatrix simplex_data;
};

class DirichletFischInterface : public DirichletStudyInterfaceBase
{
public:
    DirichletFischInterface()
    {
    }       
    virtual ~DirichletFischInterface()
    {
    }       
    void setCompositionData(const Rcpp::NumericMatrix &data) override
    {
        this->data = data;
    }
    void setSimplexData(const Rcpp::NumericMatrix &simplex_data) override
    {
        this->simplex_data = simplex_data;
    }
    bool runAnalysis() override
    {
        // Placeholder for running the analysis
        // This would typically call the Dirichlet_Fisch class methods.
        return true; // Indicating success
    }
    Rcpp::List getResults() override
    {
        // Placeholder for returning results        
        // This would typically return the results of the analysis.
        return Rcpp::List::create(Rcpp::Named("status") = "success");
    }
private:
    Rcpp::NumericMatrix data;
    Rcpp::NumericMatrix simplex_data;
};



class DirichletStudyInterface
{
public:
    DirichletStudyInterface()
    {
    }

    virtual ~DirichletStudyInterface()
    {
    }

    void setCompositionData(const Rcpp::NumericMatrix &data)
    {
        this->data = data;
    }

    void setSimplexData(const Rcpp::NumericMatrix &simplex_data)
    {
        this->simplex_data = simplex_data;
    }

    bool runAnalysis()
    {

        Dirichlet_Default<double> dirichlet_default;
        // set input values and by pass the Initialize method
        // dirichlet_default.Initialize();
        // dirichlet_default.Evaluate();

        Dirichlet_Fisch<double> dirichlet_fisch;
        // set input values and by pass the Initialize method
        // set theta values
        // dirichlet_fisch.Initialize();
        // dirichlet_fisch.Evaluate();

        Dirichlet_Linear<double> dirichlet_linear;
        // set input values and by pass the Initialize method
        // set theta values
        //  dirichlet_linear.Initialize();
        // dirichlet_linear.Evaluate();

        // Assuming the analysis is successful, we can return true.
        return true; // Indicating success
    }

    Rcpp::List getResults()
    {
        // Placeholder for returning results
        // This would typically return the results of the analysis.
        return Rcpp::List::create(Rcpp::Named("status") = "success");
    }

private:


    void prepare_inputs(DirichletStudyInterfaceBase &study, Rcpp::NumericMatrix &simplex_data)
    {
        // This function would prepare the inputs for the Dirichlet study
    }

    Rcpp::NumericMatrix data;
    Rcpp::NumericMatrix simplex_data;
};

#endif // RCPP_INTERFACE_HPP