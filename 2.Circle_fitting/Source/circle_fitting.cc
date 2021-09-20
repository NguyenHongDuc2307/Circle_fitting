#include <Eigen/Dense>
#include <iostream>
#include <random>
#include <numeric>
#include <cmath>
#include "matplotlibcpp.h"
#include <tuple>
#include <iostream>
#include <fstream>
#include <string>

using namespace Eigen;
namespace plt = matplotlibcpp;


Eigen::MatrixXd rodrigues_rot(Eigen::MatrixXd P, VectorXd n0, VectorXd n1){
  n0 = n0/n0.norm();
  n1= n1/n1.norm();  
  
  Eigen::MatrixXd P_rot = MatrixXd::Zero(P.rows(),3);


  Eigen::Matrix<double, 1, 3> n_0,n_1;
  Eigen::Matrix<double, 1, 3> k;
  Eigen::Matrix<double,1,3> temp;
  double theta;
  n_0 << n0(0), n0(1), n0(2);
  n_1 << n1(0), n1(1), n1(2);

  k = n_0.cross(n_1);
  k = k/k.norm();
  theta = acos(n0.dot(n1));


  for (int i=0; i< P.rows(); i++)
  {
     temp = P.block(i,0,1,3);
     P_rot.block(i,0,1,3) = temp*cos(theta)+ k.cross(temp)*sin(theta)+ k*k.dot(temp)*(1-cos(theta));
  }
  return P_rot;
}


Eigen::MatrixXd generate_circle_by_vectors(VectorXd t, MatrixXd C1, double r, VectorXd n, MatrixXd u){
  Eigen::Matrix<double, 1, 3> n1,u1;
  Eigen::MatrixXd P_circle(100,3);
  Eigen::VectorXd cos_t(100);
  Eigen::VectorXd sin_t(100);
  Eigen::MatrixXd concatC(100, 3);

  n1 << n(0), n(1), n(2);
  u1 << u(0,0), u(0,1), u(0,2);
  n1 = n1/n1.norm();
  u1 = u1/u1.norm();
 

  for (int i = 0; i < t.rows(); i++)
  {
    cos_t(i) = cos(t(i));
    sin_t(i) = sin(t(i));
  }

  for (int i = 0; i < concatC.rows(); i++) {
  for (int j = 0; j < concatC.cols(); j++) {
      concatC(i,j) = C1(0,j);
    }
  }

  P_circle = r*cos_t*u1 + r * sin_t*n1.cross(u1) + concatC;
  return P_circle;
}


std::tuple<double, double, double> fit_circle_2d(VectorXd x, VectorXd y){
  Eigen::MatrixXd A(x.rows(), 3);
  double xc,yc,r;

    for (int j = 0; j < A.cols(); j++) {
    for (int i = 0; i < A.rows(); i++) {
      if(j == 0){
        A(i,j) = x(i);
      } 
      else if(j == 1){
        A(i,j) = y(i);
      }
      else{
        A(i,j) = 1;
      }
    }
  }

  Eigen::VectorXd b(x.rows());
  for (int i = 0; i < b.rows(); ++i)
  {
    b(i) = pow(x(i),2)+ pow(y(i),2);
  }

  Eigen::VectorXd c(3);
  c = A.colPivHouseholderQr().solve(b); 
  xc = c(0)/2;
  yc = c(1)/2;
  r = sqrt(c(2)+pow(xc,2)+pow(yc,2));


  return {xc,yc,r};
}


int main(int argc, char *argv[]) {
  Eigen::VectorXd C(3);
  Eigen::VectorXd t(100);
  Eigen::MatrixXd P(100, 3);
  Eigen::MatrixXd m;


  std::fstream dataFile("data.txt",std::ios_base::in);
  double a;
  int row=0, col=0;
  while(dataFile >> a)
  {
    if (col <2)
    {
      P(row,col) = a;
      col++;
    }
    else
    {
      P(row,col) = a;
      row++;
      col = 0;
    }
  }

  Eigen::VectorXd P1 = P.block(0, 0, row, 1);
  Eigen::VectorXd P2 = P.block(0, 1, row, 1);
  Eigen::VectorXd P3 = P.block(0, 2, row, 1);

  
/*####################################################################################### 
  1) Fitting plane by SVD for the mean-centered data
# Eq. of plane is <p,n> + d = 0, where p is a point on plane and n is normal vector
####################################################################################### */
   
  Eigen::VectorXd P_mean(3);
  Eigen::MatrixXd P_centered(row,3);
  P_mean << P1.mean(), P2.mean(), P3.mean();


  for (int i = 0; i < P_centered.rows(); i++) {
    for (int j = 0; j < P_centered.cols(); j++) {
      P_centered(i,j) = P(i,j) - P_mean(j);
    }
  }

  Eigen::JacobiSVD<MatrixXd> svd(P_centered,ComputeThinU | ComputeThinV);

  Eigen::VectorXd normal(3);
  Eigen::VectorXd refVector(3);
  Eigen::MatrixXd V_transpose(3,3);
  double d;
  V_transpose=svd.matrixV().transpose();
  refVector << 0,0,1;

  normal << V_transpose(2,0), V_transpose(2,1), V_transpose(2,2);
  d = -P_mean.dot(normal);


/*####################################################################################### 
(2) Project points to coords X-Y in 2D plane
##################################################################################################*/

  Eigen::MatrixXd P_xy(row,3);
  P_xy =  rodrigues_rot(P_centered, normal, refVector);

  Eigen::VectorXd P_xy_1 = P_xy.block(0, 0, row, 1);
  Eigen::VectorXd P_xy_2 = P_xy.block(0, 1, row, 1);
  Eigen::VectorXd P_xy_3 = P_xy.block(0, 2, row, 1);


/*####################################################################################### 
(3) Fit circle in new 2D coords
##################################################################################################*/

  auto [xc,yc,rc] = fit_circle_2d(P_xy_1, P_xy_2);
  for (int i = 0; i < t.rows(); i++)
  {
    t(i) = i * (2.0 / 100 * M_PI);
  }
  Eigen::VectorXd xx(100);
  Eigen::VectorXd yy(100);
  Eigen::VectorXd xcc(1);
  Eigen::VectorXd xyy(1);
  xcc(0) = xc;
  xyy(0) = yc;
  for(int i =0; i< 100; i++)
  {
    xx(i) = xc+rc*cos(t(i));
    yy(i) = yc+rc*sin(t(i));
  }

/*####################################################################################### 
(4) Transform circle center back to 3D coords
##################################################################################################*/
  Eigen::VectorXd center(3);
  Eigen::MatrixXd C1(1,3);

  center << xc,yc,0;
  C1 = rodrigues_rot(center.transpose(), refVector, normal)+ P_mean.transpose();

  Eigen::MatrixXd u(1,3);
  u = P.block(0,0,1,3) - C1;

  Eigen::MatrixXd P_fitcircle(100,3);

  P_fitcircle = generate_circle_by_vectors(t, C1, rc, normal, u);


/*
##############################################################################
Generate points for fitting arc
*/

  Eigen::VectorXd P_fitcircle1 = P_fitcircle.block(0, 0, 100, 1);
  Eigen::VectorXd P_fitcircle2 = P_fitcircle.block(0, 1, 100, 1);
  Eigen::VectorXd P_fitcircle3 = P_fitcircle.block(0, 2, 100, 1);

  plt::figure();
  plt::subplot(1,3,1);
  plt::plot(P_fitcircle1,P_fitcircle2,"k-",{{"label","Fitting circle"}});
  plt::plot(P1,P2,"bo",{{"label","Cluster points P"}});
  plt::xlabel("x", {{"fontsize","18"}});
  plt::ylabel("y", {{"fontsize","18"}});
  plt::title("View X-Y");
  

  plt::subplot(1,3,2);
  plt::plot(P_fitcircle1,P_fitcircle3,"k-",{{"label","Fitting circle"}});
  plt::plot(P1,P3,"bo",{{"label","Cluster points P"}});
  plt::xlabel("x", {{"fontsize","18"}});
  plt::ylabel("z", {{"fontsize","18"}});
  plt::title("View X-Z");


  plt::subplot(1,3,3);
  plt::plot(P_fitcircle2,P_fitcircle3,"k-",{{"label","Fitting circle"}});
  plt::plot(P2,P3,"bo",{{"label","Cluster points P"}});
  plt::xlabel("y", {{"fontsize","18"}});
  plt::ylabel("z", {{"fontsize","18"}});
  plt::title("View Y-Z");


  plt::legend();
  plt::show();

/*
##############################################################################
Write fitting circle matrix to file and draw 3D figure in python.
Note: This is just a workaround because I can not find the way to plot 3D figures 
by matplotlibcpp 
*/

  std::ofstream file1("fitting_circle.txt");
  for (int j=0; j<P_fitcircle.cols(); j++)
  {
    for (int i = 0; i < P_fitcircle.rows(); ++i)
    {
      file1<<P_fitcircle(i,j);
      file1<<" ";
    }
    file1<<"\n";
  }

  file1.close();
  return 0;
}