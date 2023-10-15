///
/// Starter template for first baseball problem
/// Solve for the initial speed of the pitch given the initial parameters
/// xend : distance to home plate [18.5] m
/// z0 : height of release of ball [1.4] m
/// theta0 : angle of release above horizontal [1] degree
///
///  Do not change the interface for running the program
///  Fill in the value of vPitch in the print statement with your solution
///  at the end of main()
///

#include "RKn.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TF1.h"
#include "TCanvas.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>

using namespace std;

struct Params {
  double g;   // acceleration [m/s^2]
  double m;   // mass of object [kg], nb proj. In vacuum funcs do not depend on the mass
  double d;   // m diameter of ball
  double b;   // b,c params for air resistance
  double c;
};

double f_ri(double x, const vector<double> &y, void *params=0){ 
  (void) x;   // prevent unused variable warning
  return y[1];
}

/// \brief Change in velocity along  \f$\hat i\f$ axis
/// \param[in] x independent variable
/// \param[in] y dependent variables
double f_vi(double x, const vector<double> &y, void *params=0){ 
  (void) x;
  Params *p = (Params*)params;
  //F_x = -b*d*v_x - c*d*d*v*v_x
  //a_x = F_x/m
  double b = p->b;
  double c = p->c;
  double d = p->d;
  double m = p->m;
  double v=sqrt(y[1]*y[1] + y[3]*y[3]);
  return -(b*d*y[1] + c*d*d*v*y[1])/m;
  // return 0;  // if no air, no forces/acceleration along i direction in this problem
}

/// \brief Change in position along \f$\hat j\f$ axis
/// \param[in] x independent variable
/// \param[in] y dependent variables
///
/// Air resistance model: F= \f$k v^2\f$
///
double f_rk(double x, const vector<double> &y, void *params=0){  
  (void) x;   // prevent unused variable warning
  return y[3];
}

/// Change in velocity along  \f$\hat j\f$ axis
/// \param[in] x independent variable
/// \param[in] y dependent variables
double f_vk(double x, const vector<double> &y, void *params=0){  
  (void) x;
  Params *p = (Params*)params;
  //F_z = -b*d*v_z - c*d*d*v*v_z - mg
  //a_z = F_z/m
  double b = p->b;
  double c = p->c;
  double d = p->d;
  double m = p->m;
  double g = p->g;
  double v=sqrt(y[1]*y[1] + y[3]*y[3]);
  //printf("%lf \n", -(b*d*y[3] + c*d*d*v*y[3])/m - g);
  return -(b*d*y[3] + c*d*d*v*y[3])/m - g;
  //return -p->g;    // if no air constant acceleration along -j direction: F/m = -g
}

/// \brief Stopping condition
/// \param[in] x independent variable
/// \param[in] y dependent variables
///
/// Returns 0(1) to flag continuation(termination) of calculation
double f_stop(double x, const vector<double> &y, void *params=0){
  (void) x;
  if (y[2] < 0.9) return 1;
  //printf("%lf \n",y[2]);// stop calculation if z-height has dropped below 0.9, missing the strike zone
  return 0;  // continue calculation
}

double find_d(vector<double> y, void *params=0){
  Params *pars = (Params*)params;
  // setup default parameters
  /*
  Params pars;
  pars.g=9.81;
  pars.m=m;
  pars.air_k=0.1;
  void *p_par = (void*) &pars;

  double theta=pars.theta;   // initial angle degrees
  double v0=100;     // m/
  */

  // *** test 2: Use RK4SolveN to calculate simple projectile motion
  vector<pfunc_t> v_fun(4);   // 4 element vector of function pointers
  v_fun[0]=f_ri;
  v_fun[1]=f_vi;
  v_fun[2]=f_rk;
  v_fun[3]=f_vk;
  
  /*
  vector<double> y(4);
  // initial conditions are starting position, velocity and angle, equivalently ri,rj,vi,vj
  y[0]=0;   // init position on i-axis
  y[1]=v0*cos(theta*3.14159/180);  // init velocity along i axis
  y[2]=0;   // repeat for j-axis
  y[3]=v0*sin(theta*3.14159/180);
  //cout << "Vinit: " << v0 << " m/s" << endl;
  //cout << "Angle: " << theta << " deg" << endl;
  //cout << "(vx,vy) " << y[1] << " , "  <<  y[3] << " m/s" << endl;
  */

  
  double x=0;           // t0
  double xmax=20;  // tmax
  int nsteps=10000;
  // fixed step size algorithm
  auto tgN = RK4SolveN(v_fun, y, nsteps, x, xmax, pars, f_stop);
  // example of variable step algorithm, here the estimate accuracy is limited to 1e-4
  // in the plot you will see the change in step size thoughout the time interval
  //auto tgN = RK4SolveNA(v_fun, y, nsteps, x, xmax, pars, f_stop,1e-4);
  
  //cout << "Final velocity = " << sqrt(y[1]*y[1]+y[3]*y[3]) << endl;
  return y[0];
}


int main(int argc, char **argv){

  double xend=18.5;       // meters to plate
  double z0=1.4;             // height of release [m]
  double theta0=10;         // angle of velocity at release (degrees)
                                      // convert to radians before using!
  bool showPlot=true;    // keep this flag false by default
  
  //examples of paramaters
   Params pars;
  pars.g=9.81;
  pars.m=0.145;    
  pars.d=0.0075;   
  //pars.b=0; //no air resistance
  pars.b = 1.6e-4;  
  //pars.c=0; //no air resistance
  pars.c = 0.25;
  //pars.theta0 = theta0;
  //pars.z0 = z0;
  void *p_par = (void*) &pars;
  
  // allow changing the parameters from the command line
  int c;
  while ((c = getopt (argc, argv, "x:z:t:p")) != -1)
    switch (c) {
    case 'x':
      xend = atof(optarg);
      break;
    case 'z':
      z0 = atof(optarg);
      break;
    case 't':
      theta0 = atof(optarg);
      break;
    case 'p':
      showPlot=true;
      break;
    case '?':
      fprintf (stderr, "Unknown option `%c'.\n", optopt);
    }
  TApplication theApp("App", &argc, argv); // init ROOT App for displays


  double vPitch = 10;   // m/s of pitch needed to land in strike zone at 0.9 meters
  // write code to solve for vPitch here
  
  for (int i = 0; i < 400; i++) {
    //speed runs from 10m/s to 50m/s,
    vector<double> y(4);
    // initial conditions are starting position, velocity and angle, equivalently ri,rj,vi,vj
    y[0]=0;   // init position on i-axis
    y[1]=vPitch*cos(theta0*3.14159/180);  // init velocity along i axis
    y[2]=z0;   // repeat for k-axis; starting at pitcher's hand
    y[3]=vPitch*sin(theta0*3.14159/180);
  
    double traveled = find_d(y, p_par);
    printf("distance traveled: %lf \n", traveled);
    if (traveled > 18.5) break;
    else {
      vPitch += 0.1;
    }
  }
    // do not change these lines
  printf("********************************\n");
  printf("(xend,z0,theta0) = (%lf,%lf,%lf)\n",xend,z0,theta0);
  printf("v_pitch = %lf m/s\n",vPitch);
  printf("********************************\n");

  if (showPlot){
    cout << "Press ^c to exit" << endl;
    theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
     theApp.Run();
   }
  
  return 0;
}

