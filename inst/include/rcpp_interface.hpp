#ifndef RCPP_INTERFACE_HPP
#define RCPP_INTERFACE_HPP

#include <Rcpp.h>
#include "dirichlet_fa.hpp"
#include <cmath>
#include <sstream>

class DirichletStudyComponentBase
{
public:
    std::vector<std::vector<double>> fa_input_values;
    static uint32_t next_id;
    uint32_t id;
    static std::map<uint32_t, DirichletStudyComponentBase *> instances;
    bool verbose = false;

    DirichletStudyComponentBase() {}
    DirichletStudyComponentBase(const DirichletStudyComponentBase &other)
        : id(other.id), data(other.data), simplex_data(other.simplex_data) {}

    virtual void setCompositionData(const Rcpp::NumericMatrix &data) = 0;
    virtual void setSimplexData(const Rcpp::NumericMatrix &simplex_data) = 0;
    virtual bool runAnalysis() = 0;
    virtual Rcpp::List getResults() = 0;
    virtual void setCounts(const Rcpp::IntegerVector &counts) = 0;
    virtual void makeInputValues()
    {
        if (this->verbose)
        {
            Rcpp::Rcout << "Making input values..." << std::endl;
        }
        this->fa_input_values.clear();
        for (size_t i = 0; i < this->simplex_data.rows(); i++)
        {
            if (this->verbose)
            {
                Rcpp::Rcout << "Row " << i << ": ";
            }
            std::vector<double> input_row;
            for (size_t j = 0; j < this->simplex_data.cols(); j++)
            {
                if (this->verbose)
                {
                    Rcpp::Rcout << "Col " << j << ": " << this->simplex_data(i, j) << std::endl;
                }
                input_row.push_back(this->simplex_data(i, j));
            }
            this->fa_input_values.push_back(input_row);
        }
        if (this->verbose)
        {
            Rcpp::Rcout << "input_values size: " << this->fa_input_values.size() << "\n";
        }
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
        if (data.empty())
        {
            return Rcpp::NumericVector();
        }
        Rcpp::NumericVector vec(data.size());
        for (size_t i = 0; i < data.size(); ++i)
        {
            vec[i] = data[i];
        }
        return vec;
    }

    Rcpp::IntegerVector ToRcppIntegerVector(const std::vector<int> &data) const
    {
        Rcpp::IntegerVector vec(data.size());
        for (size_t i = 0; i < data.size(); ++i)
        {
            vec[i] = data[i];
        }
        return vec;
    }

    std::vector<int> CountsFromRcpp(const Rcpp::IntegerVector &counts) const
    {
        if (counts.size() == 0)
        {
            Rcpp::stop("counts must contain at least one category.");
        }
        std::vector<int> out(counts.size());
        int total = 0;
        for (R_xlen_t i = 0; i < counts.size(); ++i)
        {
            if (counts[i] < 0)
            {
                Rcpp::stop("counts must be non-negative.");
            }
            out[static_cast<size_t>(i)] = counts[i];
            total += counts[i];
        }
        if (total <= 0)
        {
            Rcpp::stop("counts must sum to a positive value.");
        }
        return out;
    }

    void ValidateSimplexData(size_t expected_columns) const
    {
        if (this->simplex_data.rows() == 0 || this->simplex_data.cols() == 0)
        {
            Rcpp::stop("simplex_data must be a non-empty numeric matrix.");
        }
        if (static_cast<size_t>(this->simplex_data.cols()) != expected_columns)
        {
            Rcpp::stop("simplex_data must have the same number of columns as counts.");
        }
        for (int i = 0; i < this->simplex_data.rows(); ++i)
        {
            double row_sum = 0.0;
            for (int j = 0; j < this->simplex_data.cols(); ++j)
            {
                double value = this->simplex_data(i, j);
                if (!std::isfinite(value))
                {
                    Rcpp::stop("simplex_data must contain only finite values.");
                }
                if (value <= 0.0)
                {
                    Rcpp::stop("simplex_data must be an open simplex: every entry must be greater than zero.");
                }
                row_sum += value;
            }
            if (std::fabs(row_sum - 1.0) > 1e-8)
            {
                Rcpp::stop("each simplex_data row must sum to 1 within tolerance 1e-8.");
            }
        }
    }

    double TotalCount(const std::vector<int> &counts) const
    {
        double total = 0.0;
        for (size_t i = 0; i < counts.size(); ++i)
        {
            total += static_cast<double>(counts[i]);
        }
        return total;
    }

    std::vector<double> ObservedProportions(const std::vector<int> &counts) const
    {
        double total = this->TotalCount(counts);
        std::vector<double> observed(counts.size(), 0.0);
        if (total <= 0.0)
        {
            return observed;
        }
        for (size_t i = 0; i < counts.size(); ++i)
        {
            observed[i] = static_cast<double>(counts[i]) / total;
        }
        return observed;
    }

    double EffectiveSampleSize(const std::string &variant,
                               const std::vector<int> &counts,
                               double dispersion) const
    {
        double total = this->TotalCount(counts);
        if (variant == "linear")
        {
            return (1.0 + dispersion * total) / (1.0 + dispersion);
        }
        if (variant == "saturated")
        {
            return total * (1.0 + dispersion) / (total + dispersion);
        }
        return std::numeric_limits<double>::quiet_NaN();
    }

    Rcpp::List PaperMetrics(const FunctionalAnalysis<double> &analysis,
                            const std::vector<int> &counts) const
    {
        size_t n = analysis.parameter_sets.size();
        std::vector<double> delta_loglik(n, std::numeric_limits<double>::quiet_NaN());
        std::vector<double> relative_likelihood(n, std::numeric_limits<double>::quiet_NaN());
        std::vector<double> gradient_norm(n, std::numeric_limits<double>::quiet_NaN());
        std::vector<double> max_abs_gradient(n, std::numeric_limits<double>::quiet_NaN());
        std::vector<double> boundary_distance(n, std::numeric_limits<double>::quiet_NaN());
        std::vector<double> l1_distance_from_observed(n, std::numeric_limits<double>::quiet_NaN());
        std::vector<double> l2_distance_from_observed(n, std::numeric_limits<double>::quiet_NaN());
        size_t n_derivatives = 0;
        for (size_t i = 0; i < analysis.derivatives_matrix.size(); ++i)
        {
            n_derivatives = std::max(n_derivatives, analysis.derivatives_matrix[i].size());
        }
        std::vector<double> derivative_total_variation(n_derivatives, 0.0);
        std::vector<double> mean_abs_derivative_step(n_derivatives, std::numeric_limits<double>::quiet_NaN());
        std::vector<int> derivative_step_counts(n_derivatives, 0);

        std::vector<double> observed = this->ObservedProportions(counts);
        double max_loglik = std::numeric_limits<double>::lowest();
        bool has_finite_loglik = false;
        for (size_t i = 0; i < analysis.values.size(); ++i)
        {
            if (std::isfinite(analysis.values[i]))
            {
                max_loglik = std::max(max_loglik, analysis.values[i]);
                has_finite_loglik = true;
            }
        }

        const double boundary_threshold = 0.05;
        double boundary_gradient_sum = 0.0;
        double interior_gradient_sum = 0.0;
        int boundary_count = 0;
        int interior_count = 0;

        for (size_t i = 0; i < n; ++i)
        {
            if (i < analysis.values.size() && std::isfinite(analysis.values[i]) && has_finite_loglik)
            {
                delta_loglik[i] = analysis.values[i] - max_loglik;
                relative_likelihood[i] = std::exp(delta_loglik[i]);
            }

            if (i < analysis.derivatives_matrix.size())
            {
                double squared_gradient_sum = 0.0;
                double largest_abs_gradient = 0.0;
                for (size_t j = 0; j < analysis.derivatives_matrix[i].size(); ++j)
                {
                    double derivative = analysis.derivatives_matrix[i][j];
                    if (std::isfinite(derivative))
                    {
                        squared_gradient_sum += derivative * derivative;
                        largest_abs_gradient = std::max(largest_abs_gradient, std::fabs(derivative));
                    }
                }
                gradient_norm[i] = std::sqrt(squared_gradient_sum);
                max_abs_gradient[i] = largest_abs_gradient;
            }

            if (i < analysis.parameter_sets.size())
            {
                double min_p = std::numeric_limits<double>::max();
                double l1 = 0.0;
                double l2_squared = 0.0;
                for (size_t j = 0; j < analysis.parameter_sets[i].size(); ++j)
                {
                    double p = analysis.parameter_sets[i][j];
                    min_p = std::min(min_p, p);
                    if (j < observed.size())
                    {
                        double diff = p - observed[j];
                        l1 += std::fabs(diff);
                        l2_squared += diff * diff;
                    }
                }
                boundary_distance[i] = min_p;
                l1_distance_from_observed[i] = l1;
                l2_distance_from_observed[i] = std::sqrt(l2_squared);

                if (std::isfinite(gradient_norm[i]))
                {
                    if (boundary_distance[i] < boundary_threshold)
                    {
                        boundary_gradient_sum += gradient_norm[i];
                        boundary_count++;
                    }
                    else
                    {
                        interior_gradient_sum += gradient_norm[i];
                        interior_count++;
                    }
                }
            }
        }

        for (size_t i = 1; i < analysis.derivatives_matrix.size(); ++i)
        {
            for (size_t j = 0; j < analysis.derivatives_matrix[i].size() && j < analysis.derivatives_matrix[i - 1].size(); ++j)
            {
                double current = analysis.derivatives_matrix[i][j];
                double previous = analysis.derivatives_matrix[i - 1][j];
                if (std::isfinite(current) && std::isfinite(previous))
                {
                    derivative_total_variation[j] += std::fabs(current - previous);
                    derivative_step_counts[j]++;
                }
            }
        }

        for (size_t j = 0; j < derivative_total_variation.size(); ++j)
        {
            if (derivative_step_counts[j] > 0)
            {
                mean_abs_derivative_step[j] = derivative_total_variation[j] /
                                              static_cast<double>(derivative_step_counts[j]);
            }
        }

        double boundary_gradient_ratio = std::numeric_limits<double>::quiet_NaN();
        if (boundary_count > 0 && interior_count > 0)
        {
            boundary_gradient_ratio = (boundary_gradient_sum / static_cast<double>(boundary_count)) /
                                      (interior_gradient_sum / static_cast<double>(interior_count));
        }

        return Rcpp::List::create(
            Rcpp::Named("observed_proportions") = this->ToRcppVector(observed),
            Rcpp::Named("delta_loglik") = this->ToRcppVector(delta_loglik),
            Rcpp::Named("relative_likelihood") = this->ToRcppVector(relative_likelihood),
            Rcpp::Named("gradient_norm") = this->ToRcppVector(gradient_norm),
            Rcpp::Named("max_abs_gradient") = this->ToRcppVector(max_abs_gradient),
            Rcpp::Named("boundary_distance") = this->ToRcppVector(boundary_distance),
            Rcpp::Named("boundary_threshold") = boundary_threshold,
            Rcpp::Named("boundary_gradient_ratio") = boundary_gradient_ratio,
            Rcpp::Named("l1_distance_from_observed") = this->ToRcppVector(l1_distance_from_observed),
            Rcpp::Named("l2_distance_from_observed") = this->ToRcppVector(l2_distance_from_observed),
            Rcpp::Named("derivative_total_variation") = this->ToRcppVector(derivative_total_variation),
            Rcpp::Named("mean_abs_derivative_step") = this->ToRcppVector(mean_abs_derivative_step));
    }

    Rcpp::List AnalysisResults(const FunctionalAnalysis<double> &analysis,
                               const std::vector<int> &counts,
                               const std::string &variant,
                               double dispersion) const
    {
        std::vector<double> derivative_volatility;
        typename std::map<uint32_t, double>::const_iterator it;
        for (it = analysis.stochasticity_of_derivatives.begin(); it != analysis.stochasticity_of_derivatives.end(); ++it)
        {
            derivative_volatility.push_back(it->second);
        }
        Rcpp::List paper_metrics = this->PaperMetrics(analysis, counts);

        return Rcpp::List::create(
            Rcpp::Named("variant") = variant,
            Rcpp::Named("name") = analysis.name,
            Rcpp::Named("description") = analysis.description,
            Rcpp::Named("counts") = this->ToRcppIntegerVector(counts),
            Rcpp::Named("dispersion") = dispersion,
            Rcpp::Named("effective_sample_size") = this->EffectiveSampleSize(variant, counts, dispersion),
            Rcpp::Named("runtime_seconds") = analysis.runtime,
            Rcpp::Named("n_parameter_sets") = static_cast<int>(analysis.parameter_sets.size()),
            Rcpp::Named("n_parameters") = static_cast<int>(analysis.parameters.size()),
            Rcpp::Named("is_continuous") = analysis.is_continuous,
            Rcpp::Named("min_value") = analysis.min_value,
            Rcpp::Named("max_value") = analysis.max_value,
            Rcpp::Named("parameter_set_min") = this->ToRcppVector(analysis.parameter_set_min),
            Rcpp::Named("parameter_set_max") = this->ToRcppVector(analysis.parameter_set_max),
            Rcpp::Named("mean_parameter_values") = this->ToRcppVector(analysis.mean_parameter_values),
            Rcpp::Named("input_values") = this->ToRcppMatrix(analysis.input_values),
            Rcpp::Named("parameter_sets") = this->ToRcppMatrix(analysis.parameter_sets),
            Rcpp::Named("values") = this->ToRcppVector(analysis.values),
            Rcpp::Named("derivatives_matrix") = this->ToRcppMatrix(analysis.derivatives_matrix),
            Rcpp::Named("paper_metrics") = paper_metrics,
            Rcpp::Named("derivative_volatility") = this->ToRcppVector(derivative_volatility),
            Rcpp::Named("lower_bound_covariance") = this->ToRcppMatrix(analysis.lower_bound_covariance),
            Rcpp::Named("central_bound_covariance") = this->ToRcppMatrix(analysis.central_bound_covariance),
            Rcpp::Named("upper_bound_covariance") = this->ToRcppMatrix(analysis.upper_bound_covariance),
            Rcpp::Named("lower_bound_correlation") = this->ToRcppMatrix(analysis.lower_bound_correlation),
            Rcpp::Named("central_bound_correlation") = this->ToRcppMatrix(analysis.central_bound_correlation),
            Rcpp::Named("upper_bound_correlation") = this->ToRcppMatrix(analysis.upper_bound_correlation),
            Rcpp::Named("lower_bound_derivative_covariance") = this->ToRcppMatrix(analysis.lower_bound_derivative_covariance),
            Rcpp::Named("central_bound_derivative_covariance") = this->ToRcppMatrix(analysis.central_bound_derivative_covariance),
            Rcpp::Named("upper_bound_derivative_covariance") = this->ToRcppMatrix(analysis.upper_bound_derivative_covariance),
            Rcpp::Named("lower_bound_derivative_correlation") = this->ToRcppMatrix(analysis.lower_bound_derivative_correlation),
            Rcpp::Named("central_bound_derivative_correlation") = this->ToRcppMatrix(analysis.central_bound_derivative_correlation),
            Rcpp::Named("upper_bound_derivative_correlation") = this->ToRcppMatrix(analysis.upper_bound_derivative_correlation));
    }

    virtual ~DirichletStudyComponentBase() {}

    void clear()
    {
    }

    Rcpp::NumericMatrix data;
    Rcpp::NumericMatrix simplex_data;
};

std::map<uint32_t, DirichletStudyComponentBase *> DirichletStudyComponentBase::instances;
uint32_t DirichletStudyComponentBase::next_id = 1;

class DirichletDefaultInterface : public DirichletStudyComponentBase
{
public:
    std::shared_ptr<Dirichlet_Default<double>> dirichlet_default;
    bool write_output = false;

    DirichletDefaultInterface() : DirichletStudyComponentBase()
    {
        this->id = next_id++;
        instances[this->id] = this;
        this->dirichlet_default = std::make_shared<Dirichlet_Default<double>>();
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

    void setCounts(const Rcpp::IntegerVector &counts) override
    {
        this->dirichlet_default->x = this->CountsFromRcpp(counts);
    }

    uint32_t getId() const
    {
        return this->id;
    }

    bool runAnalysis() override
    {
        this->makeInputValues();
        // dirichlet_default.input_values = this->fa_input_values;
        dirichlet_default->Initialize();
        dirichlet_default->Evaluate();
        // Placeholder for running the analysis
        // This would typically call the Dirichlet_Default class methods.
        return true; // Indicating success
    }
    virtual Rcpp::List getResults() override
    {
        std::vector<double> derivative_stochasticity;
        typename std::map<uint32_t, double>::iterator it;
        for (it = dirichlet_default->stochasticity_of_derivatives.begin(); it != dirichlet_default->stochasticity_of_derivatives.end(); ++it)
        {
            double stochasticity = it->second;
            derivative_stochasticity.push_back(it->second);
        }

        // Placeholder for returning results
        // This would typically return the results of the analysis.
        return Rcpp::List::create(
            Rcpp::Named("LowerBoundsCovariance") = this->ToRcppMatrix(dirichlet_default->lower_bound_covariance),
            Rcpp::Named("UpperBoundsCovariance") = this->ToRcppMatrix(dirichlet_default->upper_bound_covariance),
            Rcpp::Named("CentralBoundsCovariance") = this->ToRcppMatrix(dirichlet_default->central_bound_covariance),
            Rcpp::Named("LowerBoundsCorrelation") = this->ToRcppMatrix(dirichlet_default->lower_bound_correlation),
            Rcpp::Named("UpperBoundsCorrelation") = this->ToRcppMatrix(dirichlet_default->upper_bound_correlation),
            Rcpp::Named("CentralBoundsCorrelation") = this->ToRcppMatrix(dirichlet_default->central_bound_correlation),
            
            Rcpp::Named("LowerBoundsDeriviativeCovariance") = this->ToRcppMatrix(dirichlet_default->lower_bound_derivative_covariance),
            Rcpp::Named("UpperBoundsDeriviativeCovariance") = this->ToRcppMatrix(dirichlet_default->upper_bound_derivative_covariance),
            Rcpp::Named("CentralBoundsDeriviativeCovariance") = this->ToRcppMatrix(dirichlet_default->central_bound_derivative_covariance),
            Rcpp::Named("LowerBoundsDeriviativeCorrelation") = this->ToRcppMatrix(dirichlet_default->lower_bound_derivative_correlation),
            Rcpp::Named("UpperBoundsDeriviativeCorrelation") = this->ToRcppMatrix(dirichlet_default->upper_bound_derivative_correlation),
            Rcpp::Named("CentralBoundsDeriviativeCorrelation") = this->ToRcppMatrix(dirichlet_default->central_bound_derivative_correlation),

            Rcpp::Named("MeanParameterValues") = this->ToRcppVector(dirichlet_default->mean_parameter_values),
            Rcpp::Named("DerivativesMatrix") = this->ToRcppMatrix(dirichlet_default->derivatives_matrix),
            Rcpp::Named("Name") = dirichlet_default->name,
            Rcpp::Named("MinValue") = dirichlet_default->min_value,
            Rcpp::Named("MaxValue") = dirichlet_default->max_value,
            Rcpp::Named("is_continuous") = dirichlet_default->is_continuous,
            Rcpp::Named("input_values") = this->ToRcppMatrix(dirichlet_default->input_values),
            Rcpp::Named("values") = this->ToRcppVector(dirichlet_default->values),
            Rcpp::Named("stochasticity") = dirichlet_default->stochasticity,
            Rcpp::Named("derivative_stochasticity") = this->ToRcppVector(derivative_stochasticity));
    }
    void clear()
    {
        this->dirichlet_default->ClearData();
    }

private:
};

class DirichletLinearInterface : public DirichletStudyComponentBase
{
public:
    std::shared_ptr<Dirichlet_Linear<double>> dirichlet_linear;
    double theta = 1.0;
    bool write_output = false;

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

    void setCounts(const Rcpp::IntegerVector &counts) override
    {
        this->dirichlet_linear->x = this->CountsFromRcpp(counts);
    }

    bool runAnalysis() override
    {
        if (!std::isfinite(this->theta) || this->theta <= 0.0)
        {
            Rcpp::stop("theta must be finite and greater than zero.");
        }
        this->ValidateSimplexData(this->dirichlet_linear->x.size());
        this->makeInputValues();
        dirichlet_linear->ClearData();
        dirichlet_linear->input_values = this->fa_input_values;
        dirichlet_linear->theta = this->theta;
        dirichlet_linear->build_parameter_sets = false;
        dirichlet_linear->Initialize();
        dirichlet_linear->Analyze();
        if (this->write_output)
        {
            dirichlet_linear->Finalize();
        }
        // Placeholder for running the analysis
        // This would typically call the Dirichlet_Linear class methods.
        return true; // Indicating success
    }
    virtual Rcpp::List getResults() override
    {
        return this->AnalysisResults(*dirichlet_linear, dirichlet_linear->x, "linear", this->theta);
    }
    void clear()
    {
        this->dirichlet_linear->ClearData();
    }

private:
};

class DirichletSaturatedInterface : public DirichletStudyComponentBase
{
public:
    std::shared_ptr<Dirichlet_Saturated<double>> dirichlet_saturated;
    double beta = 1.0;
    bool write_output = false;

    DirichletSaturatedInterface() : DirichletStudyComponentBase()
    {
        this->id = next_id++;
        this->dirichlet_saturated = std::make_shared<Dirichlet_Saturated<double>>();
        instances[this->id] = this;
    }
    DirichletSaturatedInterface(const DirichletSaturatedInterface &other)
        : DirichletStudyComponentBase(other),
          dirichlet_saturated(other.dirichlet_saturated)
    {
    }

    virtual ~DirichletSaturatedInterface()
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

    void setCounts(const Rcpp::IntegerVector &counts) override
    {
        this->dirichlet_saturated->x = this->CountsFromRcpp(counts);
    }

    bool runAnalysis() override
    {
        if (!std::isfinite(this->beta) || this->beta <= 0.0)
        {
            Rcpp::stop("beta must be finite and greater than zero.");
        }
        this->ValidateSimplexData(this->dirichlet_saturated->x.size());
        this->makeInputValues();
        dirichlet_saturated->ClearData();
        dirichlet_saturated->input_values = this->fa_input_values;
        dirichlet_saturated->beta = this->beta;
        dirichlet_saturated->build_parameter_sets = false;
        dirichlet_saturated->Initialize();
        dirichlet_saturated->Analyze();
        if (this->write_output)
        {
            dirichlet_saturated->Finalize();
        }
        // Placeholder for running the analysis
        // This would typically call the Dirichlet_Fisch class methods.
        return true; // Indicating success
    }
    virtual Rcpp::List getResults() override
    {
        return this->AnalysisResults(*dirichlet_saturated, dirichlet_saturated->x, "saturated", this->beta);
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
        int n_studies = static_cast<int>(this->studies.size());
        Rcpp::List results(n_studies);
        Rcpp::CharacterVector names(n_studies);
        size_t index = 0;
        for (const auto &study_id : studies)
        {
            auto it = DirichletStudyComponentBase::instances.find(study_id);
            if (it == DirichletStudyComponentBase::instances.end())
            {
                Rcpp::stop("Study ID not found while collecting results.");
            }

            std::stringstream name;
            name << "study_" << study_id;
            results[index] = it->second->getResults();
            names[index] = name.str();
            index++;
        }
        results.attr("names") = names;
        return results;
    }

private:
    void prepare_inputs(DirichletStudyComponentBase &study, Rcpp::NumericMatrix &simplex_data)
    {
        // This function would prepare the inputs for the Dirichlet study
    }

    Rcpp::NumericMatrix data;
    Rcpp::NumericMatrix simplex_data;
};

void Clear()
{
    typename std::map<uint32_t, DirichletStudyComponentBase *>::iterator it;
    for (it = DirichletStudyComponentBase::instances.begin(); it != DirichletStudyComponentBase::instances.end(); ++it)
    {
        (*it).second->clear();
    }
}

#endif // RCPP_INTERFACE_HPP
