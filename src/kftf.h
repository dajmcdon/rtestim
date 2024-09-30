#ifndef __KFTF_H
#define __KFTF_H

void f1step(double y,
            double c,
            double Z,
            double H,
            const Eigen::MatrixXd& A,
            double RQR,
            Eigen::VectorXd& a,
            Eigen::MatrixXd& P,
            double& vt,
            double& Ft,
            Eigen::VectorXd& Kt);
void df1step(double y,
             double Z,
             double H,
             const Eigen::MatrixXd& A,
             double RQR,
             Eigen::VectorXd& a,
             Eigen::MatrixXd& P,
             Eigen::MatrixXd& Pinf,
             int& rankp,
             double& vt,
             double& Ft,
             double& Finf,
             Eigen::VectorXd& Kt,
             Eigen::VectorXd& Kinf);
void kftf(const Eigen::VectorXd& y,
          int k,
          const Eigen::VectorXd& weights,
          const Eigen::VectorXd& x,
          const Eigen::VectorXd& c,
          double lambda,
          Eigen::VectorXd& theta);
Eigen::VectorXd kftf_test(Eigen::VectorXd y,
                          int k,
                          Eigen::VectorXd weights,
                          Eigen::VectorXd x,
                          Eigen::VectorXd c,
                          double lambda);
#endif
