/*
energy_calc.hpp
Jacob Wiebe & James Dolman
Reu1: Nov 2017
*/

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

template <class T> //probably don't need template
T psi_function (T x){

  //M: psiDerivative
  const double epsilon = 0.001;
  T e_psi = 1 / (2 * sqrt(x + epsilon));

  return e_psi;
}


//M: computePsidashFS_brox
template <class T>
mat generate_esmooth (T u, T v){
	double height = u.n_rows;
	double width = u.n_cols;
  mat e_smooth(2*height+1, 2*width+1, fill::zeros);
  mat operator(1, 2);
  mat operator.t(1, 2);

  operator(0) = 1;
  operator(1) = -1;

  mat u_dx = conv2(u, operator);
  mat v_dx = conv2(v, operator);

  mat u_dy = conv2(u, operator.t);
  mat v_dy = conv2(v, operator.t);

  operator(1) = 1;

  mat u_dx2 = conv2(u_dx, operator/2, "valid");
  mat v_dx2 = conv2(u_dx, operator/2, "valid");

  mat u_dy2 = conv2(u_dy, operator.t/2, "valid");
  mat v_dy2 = conv2(u_dy, operator.t/2, "valid");

  mat delta_ux = conv2(u_dy2, operator/2); //t
  mat u_paritaldx = pow(u_dx, 2) + pow(delta_ux, 2); //uxpd

  mat delta_uy = conv2(u_dy2, operator.t/2); //t
  mat u_paritaldy = pow(u_dx, 2) + pow(delta_uy, 2); //uypd

  mat delta_vx = conv2(v_dy2, operator/2); //t
  mat v_paritaldx = pow(v_dx, 2) + pow(delta_vx, 2); //vxpd

  mat delta_vy = conv2(v_dy2, operator.t/2); //t
  mat v_paritaldy = pow(v_dx, 2) + pow(delta_vy, 2); //vypd


	//fill odd rows / even columns ( uypd + vypd )
	//fill even rows / odd columns ( uxpd + vxpd )

  //FROM MATLAB:
  //psidashFS( 1:2:end, 2:2:end ) = psiDerivative( uypd + vypd ) ;
  //psidashFS( 2:2:end, 1:2:end ) = psiDerivative( uxpd + vxpd ) ;
  //where psidashFS is e_smooth
	//also e_smooth is huge

	return e_smooth;
}
