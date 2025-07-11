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
RCPP_EXPOSED_CLASS(DirichletStudyInterface)

RCPP_MODULE(dirichletstudy) {
    using namespace Rcpp;

    class_<DirichletStudyInterface("DirichletStudyInterface")
        .constructor()
        .method("setCompositionData", &DirichletStudyInterface::setCompositionData)
        .method("setSimplexData", &DirichletStudyInterface::setSimplexData)
        .method("runAnalysis", &DirichletStudyInterface::runAnalysis)
        .method("getResults", &DirichletStudyInterface::getResults);
}



#endif // RCPP_INTERFACE_HPP