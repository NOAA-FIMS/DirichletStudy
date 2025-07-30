#ifndef RCPP_INTERFACE_HPP
#define RCPP_INTERFACE_HPP

#include <Rcpp.h>
#include "dirichlet_fa.hpp"

class DirichletStudyComponentBase
{
public:
    std::vector<std::vector<double>> fa_input_values;
    static uint32_t next_id;
    uint32_t id;
    static std::map<uint32_t, DirichletStudyComponentBase *> instances;

    virtual void setCompositionData(const Rcpp::NumericMatrix &data) = 0;
    virtual void setSimplexData(const Rcpp::NumericMatrix &simplex_data) = 0;
    virtual bool runAnalysis() = 0;
    virtual Rcpp::List getResults() = 0;
    virtual void makeInputValues()
    {
        //     std::cout << "Making input values..." << std::endl;
        //     this->fa_input_values.clear();
        //    for(size_t i = 0; i < this->setSimplexData.rows(); i++)
        //     {
        //         std::vector<double> input_row;
        //         for(size_t j = 0; j < this->setSimplexData.cols(); j++)
        //         {
        //             input_row.push_back(this->setSimplexData(i, j));
        //         }
        //         this->fa_input_values.push_back(input_row);
        //     }
        //     std::cout << "input_values size: " << this->fa_input_values.size() << "\n";
    }

    void test()
    {
        std::cout << "test" << std::endl;
    }

    uint32_t getId() const
    {
        return id;
    }

    virtual ~DirichletStudyComponentBase() {}
};

std::map<uint32_t, DirichletStudyComponentBase *> DirichletStudyComponentBase::instances;
uint32_t DirichletStudyComponentBase::next_id = 1;

class DirichletDefaultInterface : public DirichletStudyComponentBase
{
public:
    Dirichlet_Default<double> dirichlet_default;

    DirichletDefaultInterface() : DirichletStudyComponentBase()
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

    virtual void setSimplexData(const Rcpp::NumericMatrix &simplex_data)
    {
        this->simplex_data = simplex_data;
    }

    uint32_t getId() const
    {
        return this->id;
    }

    bool runAnalysis() override
    {
        this->makeInputValues();
        // dirichlet_default.input_values = this->fa_input_values;
        dirichlet_default.Initialize();
        dirichlet_default.Evaluate();
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

class DirichletLinearInterface : public DirichletStudyComponentBase
{
public:
    Dirichlet_Linear<double> dirichlet_linear;
    double theta = 1.0;

    DirichletLinearInterface() : DirichletStudyComponentBase()
    {
        this->id = next_id++;
        instances[this->id] = this;
    }
    virtual ~DirichletLinearInterface()
    {
    }

    uint32_t getId() const
    {
        return this->id;
    }

    void setCompositionData(const Rcpp::NumericMatrix &data)
    {
        this->data = data;
    }
    virtual void setSimplexData(const Rcpp::NumericMatrix &simplex_data)
    {
        this->simplex_data = simplex_data;
    }
    bool runAnalysis() override
    {
        this->makeInputValues();
        // dirichlet_linear.input_values = this->fa_input_values;
        dirichlet_linear.theta = this->theta;
        dirichlet_linear.Initialize();
        dirichlet_linear.Evaluate();
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

class DirichletFischInterface : public DirichletStudyComponentBase
{
public:
    Dirichlet_Fisch<double> dirichlet_fisch;
    double theta = 1.0;

    DirichletFischInterface() : DirichletStudyComponentBase()
    {
        this->id = next_id++;
        instances[this->id] = this;
    }

    virtual ~DirichletFischInterface()
    {
    }

    uint32_t getId() const
    {
        return this->id;
    }

    virtual void setCompositionData(const Rcpp::NumericMatrix &data)
    {
        this->data = data;
    }
    virtual void setSimplexData(const Rcpp::NumericMatrix &simplex_data)
    {
        this->simplex_data = simplex_data;
    }
    bool runAnalysis() override
    {
        this->makeInputValues();
        // dirichlet_fisch.input_values = this->fa_input_values;
        dirichlet_fisch.theta = this->theta;
        dirichlet_fisch.Initialize();
        dirichlet_fisch.Evaluate();
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
    std::set<uint32_t> studies;

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
        studies.insert(study);
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
        int num_successful = 0;
        // This function would run the analysis for all studies
        for (const auto &study_id : studies)
        {
            auto it = DirichletStudyComponentBase::instances.find(study_id);
            if (it != DirichletStudyComponentBase::instances.end())
            {
                num_successful++;
                DirichletStudyComponentBase *study = it->second;
                // study->setCompositionData(this->data);
                // study->setSimplexData(this->simplex_data);
                if (!study->runAnalysis())
                {
                    Rcpp::Rcerr << "Error running analysis for study ID: " << study_id << std::endl;
                    return false; // Indicating failure
                }
            }
            else
            {
                Rcpp::Rcerr << "Study ID not found: " << study_id << std::endl;
                return false; // Indicating failure
            }
        }

        std::cout << "Number of successful studies: " << num_successful << std::endl;
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
    void prepare_inputs(DirichletStudyComponentBase &study, Rcpp::NumericMatrix &simplex_data)
    {
        // This function would prepare the inputs for the Dirichlet study
    }

    Rcpp::NumericMatrix data;
    Rcpp::NumericMatrix simplex_data;
};

#endif // RCPP_INTERFACE_HPP