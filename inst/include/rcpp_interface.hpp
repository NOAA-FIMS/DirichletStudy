#ifndef RCPP_INTERFACE_HPP
#define RCPP_INTERFACE_HPP

#include <Rcpp.h>

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


    void prepare_inputs(DirichletStudyBase<double> &study, Rcpp::NumericMatrix &simplex_data)
    {
        // This function would prepare the inputs for the Dirichlet study
    }

    Rcpp::NumericMatrix data;
    Rcpp::NumericMatrix simplex_data;
};

#endif // RCPP_INTERFACE_HPP