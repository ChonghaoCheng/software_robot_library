#ifndef ROBOT_LIBRARY_BOX_AWARE_ACTIVE_SET_H
#define ROBOT_LIBRARY_BOX_AWARE_ACTIVE_SET_H

#include <Math/QPSolver.h>

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

enum class BoxAwareNumericalFailureReason
{
    None,
    ActiveWorkingSetRankDeficient,
    WorkingSetFactorizationFailed,
    WorkingSetSolutionNonfinite,
    WorkingSetResidualExceeded,
    InvalidBlockingStep
};

struct BoxAwareActiveSetResults
{
    Eigen::VectorXd solution;
    SolverResults<double> solver;
    unsigned int cacheHits{0};
    unsigned int cacheMisses{0};
    unsigned int repeatedActivations{0};
    std::vector<unsigned int> constraintActivationCounts;
    unsigned int maximumConstraintActivations{0};
    std::vector<int> releasedConstraintRows;
    unsigned int maximumConstraintsReleasedInOneIteration{0};
    std::vector<double> blockingStepFractions;
    double minimumBlockingStepFraction{1.0};
    double maximumWorkingSetRelativeResidual{0.0};
    bool numericalFailure{false};
    BoxAwareNumericalFailureReason numericalFailureReason{
        BoxAwareNumericalFailureReason::None};
    unsigned int numericalFailureIteration{0};
};

namespace box_aware_detail {

constexpr double activeWorkingSetRankTolerance=1e-10;
constexpr double workingSetRelativeResidualTolerance=1e-8;

struct BlockingStep
{
    double alpha{1.0};
    int blocker{-1};
};

inline BlockingStep blocking_step(
    const Eigen::MatrixXd &B, const Eigen::VectorXd &z,
    const Eigen::VectorXd &x, const Eigen::VectorXd &dx,
    const std::vector<int> &inactive, const int variableCount)
{
    BlockingStep output;
    for(const int row:inactive)
    {
        double movement,slack;
        if(row<variableCount){movement=dx(row);slack=z(row)-x(row);}
        else if(row<2*variableCount)
        {
            const int variable=row-variableCount;
            movement=-dx(variable);slack=z(row)+x(variable);
        }
        else{movement=B.row(row).dot(dx);slack=z(row)-B.row(row).dot(x);}
        if(movement>0.0)
        {
            const double rawRatio=slack/movement;
            if(!std::isfinite(rawRatio)) continue;
            const double ratio=std::clamp(rawRatio,0.0,1.0);
            if(ratio<output.alpha || (ratio==output.alpha
                                      && (output.blocker<0 || row<output.blocker)))
            {
                output.alpha=ratio;
                output.blocker=row;
            }
        }
    }
    return output;
}

inline double working_set_relative_residual(
    const Eigen::MatrixXd &matrix, const Eigen::VectorXd &solution,
    const Eigen::VectorXd &rightHandSide)
{
    if(!matrix.allFinite() || !solution.allFinite() || !rightHandSide.allFinite())
        return std::numeric_limits<double>::infinity();
    return (matrix*solution-rightHandSide).norm()
        / std::max(1.0,rightHandSide.norm());
}

inline void record_activation(const int row,
                              std::vector<unsigned int> &activationCounts,
                              BoxAwareActiveSetResults &output)
{
    if(row<0 || row>=static_cast<int>(activationCounts.size()))
        throw std::out_of_range("box-aware activation row out of range");
    if(activationCounts[static_cast<std::size_t>(row)]++>0)
        ++output.repeatedActivations;
    output.maximumConstraintActivations=std::max(
        output.maximumConstraintActivations,
        activationCounts[static_cast<std::size_t>(row)]);
}

} // namespace box_aware_detail

/** Active set specialized for [I;-I;dense] constraints with lazy H^-1 columns. */
inline BoxAwareActiveSetResults solve_box_aware_active_set(
    const Eigen::MatrixXd &H, const Eigen::VectorXd &f,
    const Eigen::MatrixXd &B, const Eigen::VectorXd &z,
    const Eigen::VectorXd &seed, const SolverOptions<double> &options)
{
    const int n=H.rows(), m=B.rows();
    if(H.cols()!=n || f.size()!=n || z.size()!=m || seed.size()!=n || m<2*n)
        throw std::invalid_argument("box-aware QP dimensions do not match");
    const int denseCount=m-2*n;
    for(int i=0;i<n;++i)
        if(std::abs(B(i,i)-1.0)>1e-15 || std::abs(B(n+i,i)+1.0)>1e-15
           || std::abs(B.row(i).squaredNorm()-1.0)>1e-15
           || std::abs(B.row(n+i).squaredNorm()-1.0)>1e-15)
            throw std::invalid_argument("box-aware QP requires leading +/- identity rows");
    const Eigen::LDLT<Eigen::MatrixXd> decomp=H.ldlt();
    if(decomp.info()!=Eigen::Success) throw std::runtime_error("box-aware QP H LDLT failed");
    const Eigen::VectorXd invHf=decomp.solve(f);
    Eigen::VectorXd x=seed;
    std::vector<int> active,inactive;
    std::vector<bool> seen(m,false);
    std::vector<unsigned int> activations(m,0);
    for(int row=0;row<m;++row)
    {
        const bool on=z(row)-B.row(row).dot(x)<=0.0;
        (on?active:inactive).push_back(row);
        if(on){seen[row]=true;++activations[row];}
    }
    if(static_cast<int>(active.size())>n) throw std::runtime_error("box-aware QP too many active rows");
    std::vector<std::optional<Eigen::VectorXd>> boxColumns(n),denseColumns(denseCount);
    BoxAwareActiveSetResults output;
    output.constraintActivationCounts=activations;
    output.maximumConstraintActivations=active.empty()?0U:1U;
    output.solver.terminationReason=SolverTerminationReason::MaxIterations;
    output.solver.maximumActiveSetSize=active.size();
    auto column=[&](int row) -> Eigen::VectorXd
    {
        if(row<2*n)
        {
            const int variable=row%n; const double sign=row<n?1.0:-1.0;
            auto &entry=boxColumns[variable];
            if(entry){++output.cacheHits;return sign * (*entry);}
            Eigen::VectorXd unit=Eigen::VectorXd::Zero(n);unit(variable)=1.0;
            entry=decomp.solve(unit);++output.cacheMisses;return sign * (*entry);
        }
        auto &entry=denseColumns[row-2*n];
        if(entry){++output.cacheHits;return *entry;}
        entry=decomp.solve(B.row(row).transpose());++output.cacheMisses;return *entry;
    };
    double stepNorm=0.0;
    for(unsigned int iteration=0;iteration<options.maxSteps;++iteration)
    {
        output.solver.numberOfSteps=iteration+1;
        const int count=active.size();
        Eigen::MatrixXd C(count,n),invHCt(n,count);
        for(int j=0;j<count;++j){C.row(j)=B.row(active[j]);invHCt.col(j)=column(active[j]);}
        if(count)
        {
            Eigen::FullPivLU<Eigen::MatrixXd> rankRevealing(C);
            rankRevealing.setThreshold(box_aware_detail::activeWorkingSetRankTolerance);
            if(rankRevealing.rank()!=count)
            {
                output.numericalFailure=true;
                output.numericalFailureReason=
                    BoxAwareNumericalFailureReason::ActiveWorkingSetRankDeficient;
                output.numericalFailureIteration=iteration+1;
                break;
            }
        }
        const Eigen::VectorXd invHg=x+invHf;
        Eigen::VectorXd lambda(count);
        if(count)
        {
            const Eigen::MatrixXd workingSet=C*invHCt;
            const Eigen::VectorXd rightHandSide=C*invHg;
            const Eigen::LDLT<Eigen::MatrixXd> workingSetLdlt=workingSet.ldlt();
            if(workingSetLdlt.info()!=Eigen::Success || !workingSetLdlt.isPositive())
            {
                output.numericalFailure=true;
                output.numericalFailureReason=
                    BoxAwareNumericalFailureReason::WorkingSetFactorizationFailed;
                output.numericalFailureIteration=iteration+1;
                break;
            }
            lambda=workingSetLdlt.solve(rightHandSide);
            if(workingSetLdlt.info()!=Eigen::Success || !lambda.allFinite())
            {
                output.numericalFailure=true;
                output.numericalFailureReason=
                    BoxAwareNumericalFailureReason::WorkingSetSolutionNonfinite;
                output.numericalFailureIteration=iteration+1;
                break;
            }
            const double relativeResidual=box_aware_detail::working_set_relative_residual(
                workingSet,lambda,rightHandSide);
            output.maximumWorkingSetRelativeResidual=std::max(
                output.maximumWorkingSetRelativeResidual,relativeResidual);
            if(relativeResidual>box_aware_detail::workingSetRelativeResidualTolerance)
            {
                output.numericalFailure=true;
                output.numericalFailureReason=
                    BoxAwareNumericalFailureReason::WorkingSetResidualExceeded;
                output.numericalFailureIteration=iteration+1;
                break;
            }
        }
        Eigen::VectorXd dx=-invHg;if(count)dx=invHCt*lambda-invHg;
        stepNorm=dx.norm();
        if(stepNorm<=options.stepSizeTolerance)
        {
            int freeIndex=-1;
            for(int j=0;j<count;++j)
            {
                if(lambda(j)<=1e-6) continue;
                if(freeIndex<0 || lambda(j)>lambda(freeIndex)
                   || (lambda(j)==lambda(freeIndex) && active[j]<active[freeIndex]))
                    freeIndex=j;
            }
            if(freeIndex>=0)
            {
                output.releasedConstraintRows.push_back(active[freeIndex]);
                output.maximumConstraintsReleasedInOneIteration=std::max(
                    output.maximumConstraintsReleasedInOneIteration,1U);
                inactive.push_back(active[freeIndex]);
                active.erase(active.begin()+freeIndex);
                ++output.solver.activeSetChanges;
                continue;
            }
            output.solver.terminationReason=SolverTerminationReason::Converged;
            break;
        }
        else
        {
            const auto blocking=box_aware_detail::blocking_step(B,z,x,dx,inactive,n);
            const double alpha=blocking.alpha;const int blocker=blocking.blocker;
            output.blockingStepFractions.push_back(alpha);
            output.minimumBlockingStepFraction=std::min(output.minimumBlockingStepFraction,alpha);
            if(!std::isfinite(alpha) || alpha<0.0 || alpha>1.0
               || (alpha<1.0 && blocker<0))
            {
                output.numericalFailure=true;
                output.numericalFailureReason=
                    BoxAwareNumericalFailureReason::InvalidBlockingStep;
                output.numericalFailureIteration=iteration+1;
                break;
            }
            x+=alpha*dx;
            if(alpha<1.0)
            {
                active.push_back(blocker);inactive.erase(std::remove(inactive.begin(),inactive.end(),blocker),inactive.end());
                box_aware_detail::record_activation(blocker,activations,output);seen[blocker]=true;
                ++output.solver.activeSetChanges;
            }
        }
        output.solver.maximumActiveSetSize=std::max(
            output.solver.maximumActiveSetSize,static_cast<unsigned int>(active.size()));
    }
    output.solver.finalStepSize=stepNorm;
    output.solver.objectiveFunction=x.dot(0.5*H*x+f);
    output.solver.solution=x;
    output.solver.uniqueActiveConstraints=std::count(seen.begin(),seen.end(),true);
    output.constraintActivationCounts=activations;
    output.solution=x;
    return output;
}

#endif
