#include <Math/BoxAwareActiveSet.h>

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

bool check(bool condition,const std::string &message)
{
    if(!condition)std::cerr<<"FAIL: "<<message<<'\n';
    return condition;
}

SolverOptions<double> options()
{
    SolverOptions<double> value;
    value.method="active set";
    value.stepSizeTolerance=1e-5;
    value.maxSteps=50;
    return value;
}

void box(const int n,Eigen::MatrixXd &B,Eigen::VectorXd &z,const double limit=1.0)
{
    B.resize(2*n,n);B.topRows(n).setIdentity();B.bottomRows(n)=-Eigen::MatrixXd::Identity(n,n);
    z=Eigen::VectorXd::Constant(2*n,limit);
}

bool t0_release_one()
{
    Eigen::Matrix2d H=Eigen::Matrix2d::Identity();Eigen::Vector2d f=Eigen::Vector2d::Zero();
    Eigen::MatrixXd B;Eigen::VectorXd z;box(2,B,z);Eigen::Vector2d seed(1.0,1.0);
    const auto result=solve_box_aware_active_set(H,f,B,z,seed,options());
    return check(result.solver.terminationReason==SolverTerminationReason::Converged,"T0 converged")
        && check(result.releasedConstraintRows.size()==2,"T0 released exactly two rows overall")
        && check(result.releasedConstraintRows.front()==0,"T0 largest-lambda tie uses lower global row")
        && check(result.maximumConstraintsReleasedInOneIteration==1,"T0 releases one row per stationary iteration")
        && check(result.solution.norm()<=1e-12,"T0 optimum");
}

bool t1_alpha_lower_bound()
{
    Eigen::Matrix<double,2,1> B;B<<1.0,-1.0;
    Eigen::Vector2d z;z<<1.0,1.0;
    Eigen::Vector<double,1> x;x<<1.0+1e-12;
    Eigen::Vector<double,1> dx;dx<<1.0;
    const auto blocked=box_aware_detail::blocking_step(B,z,x,dx,{0,1},1);
    return check(blocked.alpha==0.0,"T1 negative raw ratio clamps to zero")
        && check(blocked.blocker==0,"T1 blocker remains valid")
        && check(blocked.alpha>=0.0 && blocked.alpha<=1.0,"T1 bounded alpha");
}

bool t2_rank_guard()
{
    Eigen::Matrix2d H=Eigen::Matrix2d::Identity();Eigen::Vector2d f=Eigen::Vector2d::Zero();
    Eigen::Matrix<double,5,2> B;B<<1,0, 0,1, -1,0, 0,-1, 1,0;
    Eigen::Matrix<double,5,1> z;z<<1,1,1,1,1;Eigen::Vector2d seed(1,0);
    const auto result=solve_box_aware_active_set(H,f,B,z,seed,options());
    return check(result.numericalFailure,"T2 fails closed")
        && check(result.numericalFailureReason==BoxAwareNumericalFailureReason::ActiveWorkingSetRankDeficient,
                 "T2 rank reason")
        && check(result.solver.terminationReason==SolverTerminationReason::MaxIterations,
                 "T2 cannot be accepted as converged");
}

bool t3_factorization_and_residual_guard()
{
    Eigen::Matrix<double,1,1> H;H<<-1.0;Eigen::Vector<double,1> f;f<<0.0;
    Eigen::MatrixXd B;Eigen::VectorXd z;box(1,B,z);Eigen::Vector<double,1> seed;seed<<1.0;
    const auto result=solve_box_aware_active_set(H,f,B,z,seed,options());
    const Eigen::Matrix<double,1,1> S=Eigen::Matrix<double,1,1>::Identity();
    Eigen::Vector<double,1> rhs;rhs<<1.0;Eigen::Vector<double,1> wrong;wrong<<2.0;
    Eigen::Vector<double,1> nonfinite;nonfinite<<std::numeric_limits<double>::infinity();
    return check(result.numericalFailure,"T3 factorization fails closed")
        && check(result.numericalFailureReason==BoxAwareNumericalFailureReason::WorkingSetFactorizationFailed,
                 "T3 factorization reason")
        && check(box_aware_detail::working_set_relative_residual(S,wrong,rhs)>1e-8,
                 "T3 excessive residual rejected by frozen threshold")
        && check(!std::isfinite(box_aware_detail::working_set_relative_residual(S,nonfinite,rhs)),
                 "T3 nonfinite solution rejected");
}

bool t4_repeated_activation_telemetry()
{
    BoxAwareActiveSetResults result;std::vector<unsigned int> counts(4,0);
    box_aware_detail::record_activation(2,counts,result);
    box_aware_detail::record_activation(2,counts,result);
    return check(result.repeatedActivations==1,"T4 repeat count")
        && check(counts[2]==2,"T4 per-constraint count")
        && check(result.maximumConstraintActivations==2,"T4 maximum activation count");
}

bool t5_easy_qp_no_regression()
{
    Eigen::Matrix2d H;H<<2.0,0.2,0.2,1.5;Eigen::Vector2d f(-2.0,-1.5);
    Eigen::Matrix<double,5,2> B;B<<1,0,0,1,-1,0,0,-1,1,1;
    Eigen::Matrix<double,5,1> z;z<<1,1,1,1,1.2;Eigen::Vector2d seed(0.5,0.5);
    const auto result=solve_box_aware_active_set(H,f,B,z,seed,options());
    const Eigen::Vector2d expected(0.6645161290322581,0.5354838709677419);
    return check(result.solver.terminationReason==SolverTerminationReason::Converged,"T5 converged")
        && check(!result.numericalFailure,"T5 no numerical guard")
        && check((result.solution-expected).lpNorm<Eigen::Infinity>()<1e-10,"T5 optimum")
        && check((B*result.solution-z).maxCoeff()<=1e-10,"T5 feasible");
}

} // namespace

int main()
{
    if(!t0_release_one())return 1;
    if(!t1_alpha_lower_bound())return 2;
    if(!t2_rank_guard())return 3;
    if(!t3_factorization_and_residual_guard())return 4;
    if(!t4_repeated_activation_telemetry())return 5;
    if(!t5_easy_qp_no_regression())return 6;
    std::cout<<"PASS: BoxAwareActiveSet hardening T0-T5\n";
    return 0;
}
