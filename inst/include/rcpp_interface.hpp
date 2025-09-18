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

    DirichletStudyComponentBase() {}
    DirichletStudyComponentBase(const DirichletStudyComponentBase &other)
        : id(other.id), data(other.data), simplex_data(other.simplex_data) {}

    virtual void setCompositionData(const Rcpp::NumericMatrix &data) = 0;
    virtual void setSimplexData(const Rcpp::NumericMatrix &simplex_data) = 0;
    virtual bool runAnalysis() = 0;
    virtual Rcpp::List getResults() = 0;
    virtual void makeInputValues()
    {
        std::cout << "Making input values..." << std::endl;
        this->fa_input_values.clear();
        for (size_t i = 0; i < this->simplex_data.rows(); i++)
        {
            std::cout << "Row " << i << ": ";
            std::vector<double> input_row;
            for (size_t j = 0; j < this->simplex_data.cols(); j++)
            {
                std::cout << "Col " << j << ": " << this->simplex_data(i, j) << std::endl;
                input_row.push_back(this->simplex_data(i, j));
            }
            this->fa_input_values.push_back(input_row);
        }
        std::cout << "input_values size: " << this->fa_input_values.size() << "\n";
    }

    void test()
    {
        std::cout << "test" << std::endl;
    }

    uint32_t getId() const
    {
        return id;
    }

    virtual Rcpp::List getResults() const
    {
        return Rcpp::List::create(Rcpp::Named("id") = id);
    }

    Rcpp::NumericMatrix ToRcppMatrix(const std::vector<std::vector<double>> &data) const
    {
        if (data.empty())
        {
            return Rcpp::NumericMatrix();
        }
        size_t nrows = data.size();
        size_t ncols = data[0].size();
        Rcpp::NumericMatrix mat(nrows, ncols);
        for (size_t i = 0; i < nrows; ++i)
        {
            for (size_t j = 0; j < ncols; ++j)
            {
                mat(i, j) = data[i][j];
            }
        }
        return mat;
    }

    Rcpp::NumericVector ToRcppVector(const std::vector<double> &data) const
    {
        if(data.empty())
        {
            std::cout << "Data is empty, returning empty NumericVector." << std::endl;
            return Rcpp::NumericVector();
        }
        Rcpp::NumericVector vec(data.size());
        for (size_t i = 0; i < data.size(); ++i)
        {
            vec[i] = data[i];
        }
        return vec;
    }

    virtual ~DirichletStudyComponentBase() {}

    Rcpp::NumericMatrix data;
    Rcpp::NumericMatrix simplex_data;
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
    virtual Rcpp::List getResults() override
    {
        std::vector<double> derivative_stochasticity;
        typename std::map<uint32_t, double>::iterator it;
        for (it = dirichlet_linear->stochasticity_of_derivatives.begin(); it != dirichlet_linear->stochasticity_of_derivatives.end(); ++it)
        {
            double stochasticity = it->second;
            derivative_stochasticity.push_back(it->second);

            
        }

        std::cout<<"dirichlet_linear->values size: "<<dirichlet_linear->values.size()<<"\n";

        // Placeholder for returning results
        // This would typically return the results of the analysis.
        return Rcpp::List::create(
            Rcpp::Named("LowerBoundsCovariance") = this->ToRcppMatrix(dirichlet_linear->lower_bound_covariance),
            Rcpp::Named("UpperBoundsCovariance") = this->ToRcppMatrix(dirichlet_linear->upper_bound_covariance),
            Rcpp::Named("CentralBoundsCovariance") = this->ToRcppMatrix(dirichlet_linear->central_bound_covariance),
            Rcpp::Named("LowerBoundsCorrelation") = this->ToRcppMatrix(dirichlet_linear->lower_bound_correlation),
            Rcpp::Named("UpperBoundsCorrelation") = this->ToRcppMatrix(dirichlet_linear->upper_bound_correlation),
            Rcpp::Named("CentralBoundsCorrelation") = this->ToRcppMatrix(dirichlet_linear->central_bound_correlation),
            Rcpp::Named("MeanParameterValues") = this->ToRcppVector(dirichlet_linear->mean_parameter_values),
            Rcpp::Named("DerivativesMatrix") = this->ToRcppMatrix(dirichlet_linear->derivatives_matrix),
            Rcpp::Named("Name") = dirichlet_linear->name,
            Rcpp::Named("MinValue") = dirichlet_linear->min_value,
            Rcpp::Named("MaxValue") = dirichlet_linear->max_value,
            Rcpp::Named("is_continuous") = dirichlet_linear->is_continuous,
            Rcpp::Named("input_values") = this->ToRcppMatrix(dirichlet_linear->input_values),
            Rcpp::Named("values") = this->ToRcppVector(dirichlet_linear->values),
            Rcpp::Named("stochasticity") = dirichlet_linear->stochasticity,
            Rcpp::Named("derivative_stochasticity") = this->ToRcppVector(derivative_stochasticity)
        );
    }

private:
};

class DirichletLinearInterface : public DirichletStudyComponentBase
{
public:
    std::shared_ptr<Dirichlet_Linear<double>> dirichlet_linear;
    double theta = 1.0;

    DirichletLinearInterface() : DirichletStudyComponentBase()
    {
        this->id = next_id++;
        this->dirichlet_linear = std::make_shared<Dirichlet_Linear<double>>();
        instances[this->id] = this;
    }
    DirichletLinearInterface(const DirichletLinearInterface &other)
        : DirichletStudyComponentBase(other),
          dirichlet_linear(other.dirichlet_linear)
    {
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
        dirichlet_linear->input_values = this->fa_input_values;
        dirichlet_linear->theta = this->theta;
        dirichlet_linear->build_parameter_sets = false;
        dirichlet_linear->Initialize();
        dirichlet_linear->Analyze();
        dirichlet_linear->Finalize();
        // Placeholder for running the analysis
        // This would typically call the Dirichlet_Linear class methods.
        return true; // Indicating success
    }
    virtual Rcpp::List getResults() override
    {
        std::vector<double> derivative_stochasticity;
        typename std::map<uint32_t, double>::iterator it;
        for (it = dirichlet_linear->stochasticity_of_derivatives.begin(); it != dirichlet_linear->stochasticity_of_derivatives.end(); ++it)
        {
            double stochasticity = it->second;
            derivative_stochasticity.push_back(it->second);

            
        }

        std::cout<<"dirichlet_linear->values size: "<<dirichlet_linear->values.size()<<"\n";

        // Placeholder for returning results
        // This would typically return the results of the analysis.
        return Rcpp::List::create(
            Rcpp::Named("LowerBoundsCovariance") = this->ToRcppMatrix(dirichlet_linear->lower_bound_covariance),
            Rcpp::Named("UpperBoundsCovariance") = this->ToRcppMatrix(dirichlet_linear->upper_bound_covariance),
            Rcpp::Named("CentralBoundsCovariance") = this->ToRcppMatrix(dirichlet_linear->central_bound_covariance),
            Rcpp::Named("LowerBoundsCorrelation") = this->ToRcppMatrix(dirichlet_linear->lower_bound_correlation),
            Rcpp::Named("UpperBoundsCorrelation") = this->ToRcppMatrix(dirichlet_linear->upper_bound_correlation),
            Rcpp::Named("CentralBoundsCorrelation") = this->ToRcppMatrix(dirichlet_linear->central_bound_correlation),
            Rcpp::Named("MeanParameterValues") = this->ToRcppVector(dirichlet_linear->mean_parameter_values),
            Rcpp::Named("DerivativesMatrix") = this->ToRcppMatrix(dirichlet_linear->derivatives_matrix),
            Rcpp::Named("Name") = dirichlet_linear->name,
            Rcpp::Named("MinValue") = dirichlet_linear->min_value,
            Rcpp::Named("MaxValue") = dirichlet_linear->max_value,
            Rcpp::Named("is_continuous") = dirichlet_linear->is_continuous,
            Rcpp::Named("input_values") = this->ToRcppMatrix(dirichlet_linear->input_values),
            Rcpp::Named("values") = this->ToRcppVector(dirichlet_linear->values),
            Rcpp::Named("stochasticity") = dirichlet_linear->stochasticity,
            Rcpp::Named("derivative_stochasticity") = this->ToRcppVector(derivative_stochasticity)
        );
    }

private:
};

class DirichletFischInterface : public DirichletStudyComponentBase
{
public:
    std::shared_ptr<Dirichlet_Fisch<double>> dirichlet_fisch;
    double theta = 1.0;

    DirichletFischInterface() : DirichletStudyComponentBase()
    {
        this->id = next_id++;
        this->dirichlet_fisch = std::make_shared<Dirichlet_Fisch<double>>();
        instances[this->id] = this;
    }
    DirichletFischInterface(const DirichletFischInterface &other)
        : DirichletStudyComponentBase(other),
          dirichlet_fisch(other.dirichlet_fisch)
    {
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
        dirichlet_fisch->input_values = this->fa_input_values;
        dirichlet_fisch->theta = this->theta;
        dirichlet_fisch->build_parameter_sets = false;
        dirichlet_fisch->Initialize();
        dirichlet_fisch->Analyze();
        dirichlet_fisch->Finalize();
        // Placeholder for running the analysis
        // This would typically call the Dirichlet_Fisch class methods.
        return true; // Indicating success
    }
     virtual Rcpp::List getResults() override
    {
        std::vector<double> derivative_stochasticity;
        typename std::map<uint32_t, double>::iterator it;
        for (it = dirichlet_linear->stochasticity_of_derivatives.begin(); it != dirichlet_linear->stochasticity_of_derivatives.end(); ++it)
        {
            double stochasticity = it->second;
            derivative_stochasticity.push_back(it->second);

            
        }

        std::cout<<"dirichlet_linear->values size: "<<dirichlet_linear->values.size()<<"\n";

        // Placeholder for returning results
        // This would typically return the results of the analysis.
        return Rcpp::List::create(
            Rcpp::Named("LowerBoundsCovariance") = this->ToRcppMatrix(dirichlet_linear->lower_bound_covariance),
            Rcpp::Named("UpperBoundsCovariance") = this->ToRcppMatrix(dirichlet_linear->upper_bound_covariance),
            Rcpp::Named("CentralBoundsCovariance") = this->ToRcppMatrix(dirichlet_linear->central_bound_covariance),
            Rcpp::Named("LowerBoundsCorrelation") = this->ToRcppMatrix(dirichlet_linear->lower_bound_correlation),
            Rcpp::Named("UpperBoundsCorrelation") = this->ToRcppMatrix(dirichlet_linear->upper_bound_correlation),
            Rcpp::Named("CentralBoundsCorrelation") = this->ToRcppMatrix(dirichlet_linear->central_bound_correlation),
            Rcpp::Named("MeanParameterValues") = this->ToRcppVector(dirichlet_linear->mean_parameter_values),
            Rcpp::Named("DerivativesMatrix") = this->ToRcppMatrix(dirichlet_linear->derivatives_matrix),
            Rcpp::Named("Name") = dirichlet_linear->name,
            Rcpp::Named("MinValue") = dirichlet_linear->min_value,
            Rcpp::Named("MaxValue") = dirichlet_linear->max_value,
            Rcpp::Named("is_continuous") = dirichlet_linear->is_continuous,
            Rcpp::Named("input_values") = this->ToRcppMatrix(dirichlet_linear->input_values),
            Rcpp::Named("values") = this->ToRcppVector(dirichlet_linear->values),
            Rcpp::Named("stochasticity") = dirichlet_linear->stochasticity,
            Rcpp::Named("derivative_stochasticity") = this->ToRcppVector(derivative_stochasticity)
        );
    }

private:
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