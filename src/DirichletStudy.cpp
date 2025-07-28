
#include <Rcpp.h>
#include "../inst/include/rcpp_interface.hpp"

using namespace Rcpp;

RCPP_EXPOSED_CLASS(DirichletStudyInterface)
RCPP_MODULE(ds)
{
    Rcpp::class_<DirichletStudyInterface>("DirichletStudyInterface")
        .constructor()
        .method("setCompositionData", &DirichletStudyInterface::setCompositionData)
        .method("setSimplexData", &DirichletStudyInterface::setSimplexData)
        .method("runAnalysis", &DirichletStudyInterface::runAnalysis)
        .method("getResults", &DirichletStudyInterface::getResults)
       // .method("addStudy", &DirichletStudyInterface::addStudy)
        .method("clearStudies", &DirichletStudyInterface::clearStudies);

    Rcpp::class_<DirichletLinearInterface>("DirichletLinearInterface")
        .constructor()
        .method("setCompositionData", &DirichletLinearInterface::setCompositionData)
        .method("setSimplexData", &DirichletLinearInterface::setSimplexData)
        .method("runAnalysis", &DirichletLinearInterface::runAnalysis)
        .method("getResults", &DirichletLinearInterface::getResults);

    Rcpp::class_<DirichletDefaultInterface>("DirichletDefaultInterface")
        .constructor()
        .method("setCompositionData", &DirichletDefaultInterface::setCompositionData)
        .method("setSimplexData", &DirichletDefaultInterface::setSimplexData)
        .method("runAnalysis", &DirichletDefaultInterface::runAnalysis)
        .method("getResults", &DirichletDefaultInterface::getResults);

    Rcpp::class_<DirichletFischInterface>("DirichletFischInterface")
        .constructor()
        .method("setCompositionData", &DirichletFischInterface::setCompositionData)
        .method("setSimplexData", &DirichletFischInterface::setSimplexData)
        .method("runAnalysis", &DirichletFischInterface::runAnalysis)
        .method("getResults", &DirichletFischInterface::getResults);
}
