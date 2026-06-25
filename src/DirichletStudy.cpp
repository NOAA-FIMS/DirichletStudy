
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
        .method("addStudy", &DirichletStudyInterface::addStudy)
        .method("clearStudies", &DirichletStudyInterface::clearStudies);


    Rcpp::class_<DirichletDefaultInterface>("DirichletDefaultInterface")
        .constructor()
        .method("setCompositionData", &DirichletDefaultInterface::setCompositionData)
        .method("setSimplexData", &DirichletDefaultInterface::setSimplexData)
        .method("setCounts", &DirichletDefaultInterface::setCounts)
        .method("runAnalysis", &DirichletDefaultInterface::runAnalysis)
        .method("getResults", &DirichletDefaultInterface::getResults)
        .method("getId", &DirichletDefaultInterface::getId)
        .field("write_output", &DirichletDefaultInterface::write_output);


    Rcpp::class_<DirichletLinearInterface>("DirichletLinearInterface")
        .constructor()
        .method("setCompositionData", &DirichletLinearInterface::setCompositionData)
        .method("setSimplexData", &DirichletLinearInterface::setSimplexData)
        .method("setCounts", &DirichletLinearInterface::setCounts)
        .method("runAnalysis", &DirichletLinearInterface::runAnalysis)
        .method("getResults", &DirichletLinearInterface::getResults)
        .method("getId", &DirichletLinearInterface::getId)
        .field("write_output", &DirichletLinearInterface::write_output)
        .field("theta", &DirichletLinearInterface::theta);

    Rcpp::class_<DirichletSaturatedInterface>("DirichletSaturatedInterface")
        .constructor()
        .method("setCompositionData", &DirichletSaturatedInterface::setCompositionData)
        .method("setSimplexData", &DirichletSaturatedInterface::setSimplexData)
        .method("setCounts", &DirichletSaturatedInterface::setCounts)
        .method("runAnalysis", &DirichletSaturatedInterface::runAnalysis)
        .method("getResults", &DirichletSaturatedInterface::getResults)
        .method("getId", &DirichletSaturatedInterface::getId)
        .field("write_output", &DirichletSaturatedInterface::write_output)
        .field("beta", &DirichletSaturatedInterface::beta);
}
