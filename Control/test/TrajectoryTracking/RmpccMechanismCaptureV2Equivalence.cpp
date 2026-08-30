#include <Control/TrajectoryTracking/SerialLinkRMPCC.h>
#include <Model/KinematicTree.h>
#include <Eigen/Geometry>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <system_error>
#include <vector>

namespace {
using RobotLibrary::Control::RmpccParameters;
using RobotLibrary::Control::SerialLinkParameters;
using RobotLibrary::Control::SerialLinkRMPCC;
using RobotLibrary::Model::Pose;
using RobotLibrary::Trajectory::CartesianSpline;
struct Fixture {
    Fixture() {
        path=std::filesystem::temp_directory_path()/("rmpcc_capture_v2_"+
            std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())+".urdf");
        std::ofstream s(path); s << R"(<robot name="capture"><link name="base"/><link name="arm"/><link name="tool"/>
<joint name="joint" type="revolute"><parent link="base"/><child link="arm"/><axis xyz="0 0 1"/>
<limit lower="-3.14" upper="3.14" effort="10" velocity="1"/></joint>
<joint name="tool_fixed" type="fixed"><parent link="arm"/><child link="tool"/><origin xyz="0.1 0 0" rpy="0 0 0"/></joint></robot>)";
    }
    ~Fixture(){std::error_code e;std::filesystem::remove(path,e);} std::filesystem::path path;
};
CartesianSpline make_path(){std::vector<Pose> p;std::vector<double> t;for(int i=0;i<9;++i){double u=double(i)/8.;p.emplace_back(Eigen::Vector3d(.1+.04*u,.01*u,0),Eigen::Quaterniond(Eigen::AngleAxisd(.35*u,Eigen::Vector3d::UnitY())));t.push_back(8*u);}return CartesianSpline(p,t,Eigen::Vector<double,6>::Zero());}
}
int main(int argc,char **argv){
    if(argc!=2)return 2; Fixture fixture;
    auto model=std::make_shared<RobotLibrary::Model::KinematicTree>(fixture.path.string());
    model->update_state(Eigen::VectorXd::Zero(1),Eigen::VectorXd::Zero(1));
    SerialLinkParameters serial;serial.qpsolver.maxSteps=50;serial.qpsolver.stepSizeTolerance=1e-5;
    RmpccParameters r;r.horizonSteps=20;r.referenceMotion=RobotLibrary::Control::RmpccReferenceMotion::FiniteStageExact;r.predictorGeometry=RobotLibrary::Control::RmpccPredictorGeometry::ExactSE3;
    const auto trajectory=make_path();
    unsetenv("ROBOT_LIBRARY_RMPCC_MECHANISM_CAPTURE_V2_DIR");
    SerialLinkRMPCC a(model,"tool",serial,r);a.set_trajectory(trajectory);a.step(.002,.2);const auto da=a.diagnostics();
    setenv("ROBOT_LIBRARY_RMPCC_MECHANISM_CAPTURE_V2_DIR",argv[1],1);setenv("ROBOT_LIBRARY_RMPCC_MECHANISM_CAPTURE_V2_STEPS","0",1);
    SerialLinkRMPCC b(model,"tool",serial,r);b.set_trajectory(trajectory);b.step(.002,.2);const auto db=b.diagnostics();
    const double u=(da.bodyTwist-db.bodyTwist).cwiseAbs().maxCoeff(),s=std::abs(da.progressRate-db.progressRate);
    const bool pass=u<=1e-12&&s<=1e-12&&da.qpConverged==db.qpConverged&&da.qpIterations==db.qpIterations&&da.hessianHash==db.hessianHash&&da.linearTermHash==db.linearTermHash;
    std::cout<<"u_inf="<<u<<" sdot="<<s<<" iterations_a="<<da.qpIterations<<" iterations_b="<<db.qpIterations<<" h_hash_equal="<<(da.hessianHash==db.hessianHash)<<" f_hash_equal="<<(da.linearTermHash==db.linearTermHash)<<'\n';return pass?0:1;
}
