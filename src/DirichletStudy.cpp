
#include <Rcpp.h>
#include "../inst/include/rcpp_interface.hpp"

using namespace Rcpp;


RCPP_EXPOSED_CLASS(DirichletStudyInterface)
RCPP_MODULE(ds) {
    Rcpp::class_<DirichletStudyInterface>("DirichletStudyInterface")
        .constructor()
        .method("setCompositionData", &DirichletStudyInterface::setCompositionData)
        .method("setSimplexData", &DirichletStudyInterface::setSimplexData)
        .method("runAnalysis", &DirichletStudyInterface::runAnalysis)
        .method("getResults", &DirichletStudyInterface::getResults);
}
