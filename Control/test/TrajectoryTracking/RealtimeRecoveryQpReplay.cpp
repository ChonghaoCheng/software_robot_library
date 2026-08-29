#include <Math/QPSolver.h>

#include <Eigen/Dense>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

Eigen::MatrixXd read_csv(const std::filesystem::path &path)
{
    std::ifstream stream(path);
    if(!stream) throw std::runtime_error("Unable to open " + path.string());
    std::vector<std::vector<double>> rows;
    std::string line;
    while(std::getline(stream, line))
    {
        if(line.empty()) continue;
        std::vector<double> row;
        std::stringstream parser(line);
        std::string field;
        while(std::getline(parser, field, ',')) row.push_back(std::stod(field));
        if(!rows.empty() && row.size() != rows.front().size())
            throw std::runtime_error("Ragged CSV: " + path.string());
        rows.push_back(std::move(row));
    }
    if(rows.empty()) return Eigen::MatrixXd(0, 0);
    Eigen::MatrixXd result(rows.size(), rows.front().size());
    for(Eigen::Index r = 0; r < result.rows(); ++r)
        for(Eigen::Index c = 0; c < result.cols(); ++c)
            result(r,c) = rows[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)];
    return result;
}

void write_vector(const std::filesystem::path &path, const Eigen::VectorXd &value)
{
    std::ofstream stream(path);
    stream << std::setprecision(17);
    for(Eigen::Index i = 0; i < value.size(); ++i) stream << value(i) << '\n';
}

struct SpecializedResult
{
    Eigen::VectorXd solution;
    unsigned int steps{0};
    bool converged{false};
    double finalStepSize{0.0};
    unsigned int activeSetChanges{0};
    unsigned int maximumActiveSetSize{0};
    unsigned int uniqueActiveConstraints{0};
    unsigned int cacheHits{0};
    unsigned int cacheMisses{0};
    unsigned int repeatedActivations{0};
};

SpecializedResult solve_box_aware(
    const Eigen::MatrixXd &H, const Eigen::VectorXd &f,
    const Eigen::MatrixXd &B, const Eigen::VectorXd &z,
    const Eigen::VectorXd &seed, const double tolerance, const unsigned int maxSteps)
{
    const int n = static_cast<int>(H.rows());
    if(B.rows() < 2 * n) throw std::invalid_argument("Expected 2*n box rows");
    const int m = static_cast<int>(B.rows());
    const int denseCount = m - 2 * n;
    for(int i = 0; i < n; ++i)
    {
        if(std::abs(B(i,i) - 1.0) > 1e-15 || std::abs(B(n+i,i) + 1.0) > 1e-15
           || (B.row(i).squaredNorm() - 1.0) > 1e-15
           || (B.row(n+i).squaredNorm() - 1.0) > 1e-15)
            throw std::invalid_argument("Frozen QP box-row layout is not +/- identity");
    }
    const Eigen::LDLT<Eigen::MatrixXd> decomposition = H.ldlt();
    if(decomposition.info() != Eigen::Success) throw std::runtime_error("H LDLT failed");
    const Eigen::VectorXd invHf = decomposition.solve(f);
    Eigen::VectorXd x = seed;
    std::vector<int> active, inactive;
    std::vector<bool> everActive(static_cast<std::size_t>(m), false);
    std::vector<unsigned int> activationCount(static_cast<std::size_t>(m), 0);
    for(int row = 0; row < m; ++row)
    {
        const bool onBoundary = z(row) - B.row(row).dot(x) <= 0.0;
        (onBoundary ? active : inactive).push_back(row);
        if(onBoundary) { everActive[static_cast<std::size_t>(row)] = true; ++activationCount[row]; }
    }
    if(static_cast<int>(active.size()) > n) throw std::runtime_error("Too many active rows");

    // H changes per QP. Cache only lazy columns within this solve; +/- box rows share a column.
    std::vector<std::optional<Eigen::VectorXd>> variableColumns(static_cast<std::size_t>(n));
    std::vector<std::optional<Eigen::VectorXd>> denseColumns(static_cast<std::size_t>(denseCount));
    SpecializedResult result;
    result.maximumActiveSetSize = static_cast<unsigned int>(active.size());
    auto cached_column = [&](const int row) -> Eigen::VectorXd
    {
        if(row < 2 * n)
        {
            const int variable = row % n;
            const double sign = row < n ? 1.0 : -1.0;
            auto &entry = variableColumns[static_cast<std::size_t>(variable)];
            if(entry) { ++result.cacheHits; return sign * *entry; }
            Eigen::VectorXd unit = Eigen::VectorXd::Zero(n); unit(variable) = 1.0;
            entry = decomposition.solve(unit); ++result.cacheMisses;
            return sign * *entry;
        }
        auto &entry = denseColumns[static_cast<std::size_t>(row - 2*n)];
        if(entry) { ++result.cacheHits; return *entry; }
        entry = decomposition.solve(B.row(row).transpose()); ++result.cacheMisses;
        return *entry;
    };

    double stepNorm = std::numeric_limits<double>::infinity();
    for(unsigned int iteration = 0; iteration < maxSteps; ++iteration)
    {
        result.steps = iteration + 1;
        const int count = static_cast<int>(active.size());
        Eigen::MatrixXd C(count, n), invHCt(n, count);
        for(int j = 0; j < count; ++j)
        {
            C.row(j) = B.row(active[static_cast<std::size_t>(j)]);
            invHCt.col(j) = cached_column(active[static_cast<std::size_t>(j)]);
        }
        const Eigen::VectorXd invHg = x + invHf;
        Eigen::VectorXd lambda(count);
        if(count > 0) lambda = (C * invHCt).ldlt().solve(C * invHg);
        Eigen::VectorXd dx = -invHg;
        if(count > 0) dx = invHCt * lambda - invHg;
        stepNorm = dx.norm();
        if(stepNorm <= tolerance)
        {
            std::vector<int> freeIndices;
            for(int j = 0; j < count; ++j)
                if(lambda(j) > 1e-6) freeIndices.push_back(j);
            if(!freeIndices.empty())
            {
                for(auto iterator = freeIndices.rbegin(); iterator != freeIndices.rend(); ++iterator)
                {
                    inactive.push_back(active[static_cast<std::size_t>(*iterator)]);
                    active.erase(active.begin() + *iterator);
                    ++result.activeSetChanges;
                }
            }
            else
            {
                // The frozen 1e-5 tolerance identifies a stable working set;
                // apply its remaining constrained Newton step rather than
                // returning the pre-step iterate as the legacy solver does.
                double alpha = 1.0;
                int blocker = -1;
                for(const int row : inactive)
                {
                    double movement;
                    double slack;
                    if(row < n) { movement = dx(row); slack = z(row) - x(row); }
                    else if(row < 2*n)
                    { const int v = row-n; movement = -dx(v); slack = z(row) + x(v); }
                    else { movement = B.row(row).dot(dx); slack = z(row)-B.row(row).dot(x); }
                    if(movement > 0.0)
                    {
                        const double ratio = slack / movement;
                        if(ratio < alpha) { alpha = ratio; blocker = row; }
                    }
                }
                x += alpha * dx;
                if(alpha < 1.0)
                {
                    active.push_back(blocker);
                    inactive.erase(std::remove(inactive.begin(), inactive.end(), blocker), inactive.end());
                    if(activationCount[static_cast<std::size_t>(blocker)]++ > 0)
                        ++result.repeatedActivations;
                    everActive[static_cast<std::size_t>(blocker)] = true;
                    ++result.activeSetChanges;
                    result.maximumActiveSetSize = std::max(
                        result.maximumActiveSetSize, static_cast<unsigned int>(active.size()));
                }
                else { result.converged = true; break; }
            }
        }
        else
        {
            double alpha = 1.0;
            int blocker = -1;
            for(const int row : inactive)
            {
                double movement;
                double slack;
                if(row < n) { movement = dx(row); slack = z(row) - x(row); }
                else if(row < 2*n)
                { const int v = row-n; movement = -dx(v); slack = z(row) + x(v); }
                else { movement = B.row(row).dot(dx); slack = z(row)-B.row(row).dot(x); }
                if(movement > 0.0)
                {
                    const double ratio = slack / movement;
                    if(ratio < alpha) { alpha = ratio; blocker = row; }
                }
            }
            x += alpha * dx;
            if(alpha < 1.0)
            {
                active.push_back(blocker);
                inactive.erase(std::remove(inactive.begin(), inactive.end(), blocker), inactive.end());
                if(activationCount[static_cast<std::size_t>(blocker)]++ > 0)
                    ++result.repeatedActivations;
                everActive[static_cast<std::size_t>(blocker)] = true;
                ++result.activeSetChanges;
                result.maximumActiveSetSize = std::max(
                    result.maximumActiveSetSize, static_cast<unsigned int>(active.size()));
            }
        }
    }
    result.finalStepSize = stepNorm;
    result.uniqueActiveConstraints = static_cast<unsigned int>(
        std::count(everActive.begin(), everActive.end(), true));
    result.solution = x;
    return result;
}

double percentile(std::vector<double> values, double p)
{
    std::sort(values.begin(), values.end());
    const double x = p * static_cast<double>(values.size()-1);
    const auto i = static_cast<std::size_t>(x);
    const double a = x - static_cast<double>(i);
    return values[i]*(1.0-a) + values[std::min(i+1,values.size()-1)]*a;
}

} // namespace

int main(int argc, char **argv)
{
    if(argc != 2) { std::cerr << "Usage: realtime_recovery_qp_replay QP_DIR\n"; return 2; }
    try
    {
        const std::filesystem::path directory(argv[1]);
        const Eigen::MatrixXd H = read_csv(directory/"H.csv");
        const Eigen::VectorXd f = read_csv(directory/"f.csv").col(0);
        Eigen::MatrixXd A = read_csv(directory/"A.csv");
        if(A.size() == 0) A.resize(0, H.cols());
        Eigen::VectorXd y = A.rows() ? read_csv(directory/"y.csv").col(0) : Eigen::VectorXd(0);
        const Eigen::MatrixXd B = read_csv(directory/"B.csv");
        const Eigen::VectorXd z = read_csv(directory/"z.csv").col(0);
        const Eigen::VectorXd seed = read_csv(directory/"seed.csv").col(0);

        SolverOptions<double> highOptions;
        highOptions.method = "active set";
        highOptions.stepSizeTolerance = 1e-12;
        highOptions.maxSteps = 2000;
        QPSolver<double> highSolver(highOptions);
        const Eigen::VectorXd high = highSolver.solve(H,f,A,y,B,z,seed);
        write_vector(directory/"high_accuracy_active_set_solution.csv", high);

        const auto candidate = solve_box_aware(H,f,B,z,seed,1e-5,50);
        write_vector(directory/"box_aware_candidate_solution.csv", candidate.solution);

        constexpr int repetitions = 50;
        std::vector<double> genericTimes, candidateTimes;
        genericTimes.reserve(repetitions); candidateTimes.reserve(repetitions);
        SolverOptions<double> productionOptions;
        productionOptions.method="active set"; productionOptions.stepSizeTolerance=1e-5;
        productionOptions.maxSteps=50;
        for(int i=-5; i<repetitions; ++i)
        {
            auto start=std::chrono::steady_clock::now();
            QPSolver<double> solver(productionOptions);
            const auto ignored=solver.solve(H,f,A,y,B,z,seed); (void)ignored;
            auto middle=std::chrono::steady_clock::now();
            const auto ignoredCandidate=solve_box_aware(H,f,B,z,seed,1e-5,50);
            (void)ignoredCandidate;
            auto end=std::chrono::steady_clock::now();
            if(i>=0)
            {
                genericTimes.push_back(std::chrono::duration<double,std::micro>(middle-start).count());
                candidateTimes.push_back(std::chrono::duration<double,std::micro>(end-middle).count());
            }
        }
        const auto mean=[](const std::vector<double>& v)
        { return std::accumulate(v.begin(),v.end(),0.0)/static_cast<double>(v.size()); };
        const auto highResults=highSolver.results();
        std::cout << std::setprecision(17) << "{"
                  << "\"high_accuracy_steps\":" << highResults.numberOfSteps << ','
                  << "\"high_accuracy_converged\":"
                  << (highResults.terminationReason==SolverTerminationReason::Converged) << ','
                  << "\"candidate_steps\":" << candidate.steps << ','
                  << "\"candidate_converged\":" << candidate.converged << ','
                  << "\"candidate_final_step_size\":" << candidate.finalStepSize << ','
                  << "\"candidate_active_set_changes\":" << candidate.activeSetChanges << ','
                  << "\"candidate_maximum_active_set_size\":" << candidate.maximumActiveSetSize << ','
                  << "\"candidate_unique_active_constraints\":" << candidate.uniqueActiveConstraints << ','
                  << "\"candidate_cache_hits\":" << candidate.cacheHits << ','
                  << "\"candidate_cache_misses\":" << candidate.cacheMisses << ','
                  << "\"candidate_repeated_activations\":" << candidate.repeatedActivations << ','
                  << "\"generic_timing_mean_us\":" << mean(genericTimes) << ','
                  << "\"generic_timing_p95_us\":" << percentile(genericTimes,.95) << ','
                  << "\"candidate_timing_mean_us\":" << mean(candidateTimes) << ','
                  << "\"candidate_timing_p95_us\":" << percentile(candidateTimes,.95)
                  << "}\n";
        return 0;
    }
    catch(const std::exception &error)
    { std::cerr << error.what() << '\n'; return 1; }
}
