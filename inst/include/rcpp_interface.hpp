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
        dirichlet_default.Initialize();

        Diririchlet_Fisch<double> dirichlet_fisch;
        dirichlet_fisch.Initialize();

        Dirichlet_Linear<double> dirichlet_linear;
        dirichlet_linear.Initialize();

        // Placeholder for analysis logic
        // This would typically involve calling the DirichletComposition methods
        // and processing the data.
        return true; // Indicating success
    }

    Rcpp::List getResults()
    {
        // Placeholder for returning results
        // This would typically return the results of the analysis.
        return Rcpp::List::create(Rcpp::Named("status") = "success");
    }

private:
    Rcpp::NumericMatrix data;
    Rcpp::NumericMatrix simplex_data;
};



#endif // RCPP_INTERFACE_HPP