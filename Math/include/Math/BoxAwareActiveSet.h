#ifndef ROBOT_LIBRARY_BOX_AWARE_ACTIVE_SET_H
#define ROBOT_LIBRARY_BOX_AWARE_ACTIVE_SET_H

#include <Math/QPSolver.h>

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <vector>

struct BoxAwareActiveSetResults
{
    Eigen::VectorXd solution;
    SolverResults<double> solver;
    unsigned int cacheHits{0};
    unsigned int cacheMisses{0};
    unsigned int repeatedActivations{0};
};

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
    auto blocking=[&](const Eigen::VectorXd &dx)
    {
        double alpha=1.0;int blocker=-1;
        for(const int row:inactive)
        {
            double movement,slack;
            if(row<n){movement=dx(row);slack=z(row)-x(row);}
            else if(row<2*n){const int v=row-n;movement=-dx(v);slack=z(row)+x(v);}
            else{movement=B.row(row).dot(dx);slack=z(row)-B.row(row).dot(x);}
            if(movement>0.0){const double ratio=slack/movement;if(ratio<alpha){alpha=ratio;blocker=row;}}
        }
        return std::pair<double,int>(alpha,blocker);
    };
    double stepNorm=0.0;
    for(unsigned int iteration=0;iteration<options.maxSteps;++iteration)
    {
        output.solver.numberOfSteps=iteration+1;
        const int count=active.size();
        Eigen::MatrixXd C(count,n),invHCt(n,count);
        for(int j=0;j<count;++j){C.row(j)=B.row(active[j]);invHCt.col(j)=column(active[j]);}
        const Eigen::VectorXd invHg=x+invHf;
        Eigen::VectorXd lambda(count);
        if(count)lambda=(C*invHCt).ldlt().solve(C*invHg);
        Eigen::VectorXd dx=-invHg;if(count)dx=invHCt*lambda-invHg;
        stepNorm=dx.norm();
        if(stepNorm<=options.stepSizeTolerance)
        {
            std::vector<int> freeIndices;
            for(int j=0;j<count;++j)if(lambda(j)>1e-6)freeIndices.push_back(j);
            if(!freeIndices.empty())
            {
                for(auto it=freeIndices.rbegin();it!=freeIndices.rend();++it)
                {inactive.push_back(active[*it]);active.erase(active.begin()+*it);++output.solver.activeSetChanges;}
                continue;
            }
            const auto [alpha,blocker]=blocking(dx);x+=alpha*dx;
            if(alpha>=1.0){output.solver.terminationReason=SolverTerminationReason::Converged;break;}
            active.push_back(blocker);inactive.erase(std::remove(inactive.begin(),inactive.end(),blocker),inactive.end());
            if(activations[blocker]++>0)++output.repeatedActivations;seen[blocker]=true;
            ++output.solver.activeSetChanges;
        }
        else
        {
            const auto [alpha,blocker]=blocking(dx);x+=alpha*dx;
            if(alpha<1.0)
            {
                active.push_back(blocker);inactive.erase(std::remove(inactive.begin(),inactive.end(),blocker),inactive.end());
                if(activations[blocker]++>0)++output.repeatedActivations;seen[blocker]=true;
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
    output.solution=x;
    return output;
}

#endif
