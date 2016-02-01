#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin);
void step(const cmplx* const psi0, cmplx* const psi1, const double dx,
	  const double dt, const int N, const double k, const double xmin);
void writeToFile(const cmplx* const v, const string s, const double dx,
         const int Nx, const double xmin, const double alpha,
         const double lambda, const double omega, const double t);
//-----------------------------------
int main(){

	const int Nx = 300;
	const double xmin = -40;
	const double xmax = 40;
	const double Tend = 10*M_PI;
	const double dx = ( xmax-xmin ) / ( Nx - 1 );
	const double dt = Tend / 1000;
	double t = 0;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);

	const double lambda = 10;
	const double omega = 0.2;
	const double k = omega * omega;
	const double alpha = sqrt(omega);
	

  stringstream strm;

	cmplx* psi0 = new cmplx[Nx];
	cmplx* psi1 = new cmplx[Nx];
	cmplx* h;

	init(psi0, alpha, lambda, dx, dt, Nx, xmin);

	writeToFile(psi0,"psi_0", dx,Nx,xmin, alpha, lambda, omega,t);



	for (int i = 1; i <= Na; i++) {
		for (int j = 1; j <= Nk-1; j++) {
		  step(psi0, psi1, dx, dt, Nx, k, xmin);
		  h = psi0;
		  psi0 = psi1;
		  psi1 = h;
		  t+=dt;
		}
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin, alpha, lambda, omega,t);
	}
  cout << "t = " << t << endl;
	
	delete[] psi0;
	delete[] psi1;
	return 0;
}
//-----------------------------------
void step(const cmplx* const psi0, cmplx* const psi1, const double dx, const double dt, const int N, const double k, const double xmin)
{
    cmplx a = cmplx( 0, -dt / (4 * dx * dx ) );
    cmplx* d = new cmplx[N];
    cmplx acc_up = cmplx( 0, dt / (4 * dx * dx ) );
    cmplx* acc_d = new cmplx[N];
    cmplx* dcc = new cmplx[N];    
    
    for(int j = 0; j<N; j++) d[j] = cmplx( 1, dt / ( 2 * dx * dx ) + dt / 2 * ( 0.5 * k * (xmin + j*dx)*(xmin + j*dx) ) );
    for(int j = 0; j<N; j++) dcc[j] = cmplx( 1, - dt / ( 2 * dx * dx ) - dt / 2 * ( 0.5 * k * (xmin + j*dx)*(xmin + j*dx) ) );
    for(int j = 0; j<N; j++) acc_d[j] = cmplx( 0, dt / (4 * dx * dx ) );
    
    for(int j = 1; j<N; j++) {
      d[j] -= 		a/d[j-1] * a;
      dcc[j] -= 	a/d[j-1] * acc_up;
      acc_d[j] -= 	a/d[j-1] * dcc[j-1];
    }
    
    psi1[N-1] = ( dcc[N-1] / d[N-1] ) * psi0[N-1];
    for(int j = N-2; j>=0; j--) psi1[j] = ( dcc[j]*psi0[j] + acc*psi0[j+1] - a*psi1[j+1] ) / d[j];
    
    delete[] d;
    delete[] dcc;
}
//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin, const double alpha,
                 const double lambda, const double omega, const double t)
{
	ofstream out(s.c_str());
  double x, xi, xil;
  double h1, h2, h3;
  cmplx ana;
	for(int i=0; i<Nx; i++){
		x = xmin + i * dx;
    xi = alpha * x;
    xil = alpha * lambda;
    h1 = -0.5 * pow(xi - xil*cos(omega*t),2 );
    h2 = omega*t/2 + xi * xil * sin(omega*t);
    h3 =  - 0.25 * xil*xil* sin(2*omega*t);
    ana = cmplx( h1 , h2 + h3  );
    ana = sqrt(alpha / sqrt(M_PI) ) * exp(ana);
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag()
         << "\t" << norm(ana) << "\t" << ana.real() << "\t" << ana.imag() <<  endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin)
{
	const double x0 = dx*Nx * 0.5; // why?
	for(int i=0;i<Nx; i++){
		double x = xmin + i*dx ;
		psi0[i] = sqrt(alpha/sqrt(M_PI)) * exp(- pow(alpha*(x-lambda),2)/2 );
	}
}
