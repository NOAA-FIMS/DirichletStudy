#ifndef RCPP_INTERFACE_HPP
#define RCPP_INTERFACE_HPP

#include <Rcpp.h>
#include "dirichlet_fa.hpp"

class DirichletStudyInterfaceBase
{
public:
    static uint32_t next_id;
    uint32_t id;
    static std::map<uint32_t, DirichletStudyInterfaceBase *> instances;

    virtual void setCompositionData(const Rcpp::NumericMatrix &data) = 0;
    virtual void setSimplexData(const Rcpp::NumericMatrix &simplex_data) = 0;
    virtual bool runAnalysis() = 0;
    virtual Rcpp::List getResults() = 0;
    uint32_t getId() const
    {
        return id;
    }

    virtual ~DirichletStudyInterfaceBase() {}
};

std::map<uint32_t, DirichletStudyInterfaceBase *> DirichletStudyInterfaceBase::instances;

class DirichletDefaultInterface : public DirichletStudyInterfaceBase
{
public:
    Dirichlet_Default<double> dirichlet_default;

    DirichletDefaultInterface():DirichletStudyInterfaceBase
    {
        this->id = next_id++;
        instances[this->id] = this;
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

    uint32_t getId() const override
    {
        return this->id;
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
    Dirichlet_Linear<double> dirichlet_linear;

    DirichletLinearInterface():DirichletStudyInterfaceBase
    {
        this->id = next_id++;
        instances[this->id] = this;
    }
    virtual ~DirichletLinearInterface()
    {
    }

    uint32_t getId() const override
    {
        return this->id;
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
    Dirichlet_Fisch<double> dirichlet_fisch;

    DirichletFischInterface():DirichletStudyInterfaceBase
    {
        this->id = next_id++;
        instances[this->id] = this;
    }

    virtual ~DirichletFischInterface()
    {
    }

    uint32_t getId() const override
    {
        return this->id;
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
    std::vector<uint32_t> studies;

    DirichletStudyInterface()
    {
    }

    DirichletStudyInterface(const DirichletStudyInterface &other)
    {
        this->studies = other.studies;
    }

    virtual ~DirichletStudyInterface()
    {
    }

    void addStudy(uint32_t study)
    {
        studies.push_back(study);
    }

    void clearStudies()
    {
        studies.clear();
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