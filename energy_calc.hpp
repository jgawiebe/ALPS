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
  mat temp(1, 2);
  mat temp.t(1, 2);

  temp(0) = 1;
  temp(1) = -1;

  mat u_dx = conv2(u, temp);
  mat v_dx = conv2(v, temp);

  mat u_dy = conv2(u, temp.t);
  mat v_dy = conv2(v, temp.t);

  temp(1) = 1;

  mat u_dx2 = conv2(u_dx, temp/2, "valid");
  mat v_dx2 = conv2(u_dx, temp/2, "valid");

  mat u_dy2 = conv2(u_dy, temp.t/2, "valid");
  mat v_dy2 = conv2(u_dy, temp.t/2, "valid");

  mat delta_ux = conv2(u_dy2, temp/2); //t
  mat u_paritaldx = pow(u_dx, 2) + pow(delta_ux, 2); //uxpd

  mat delta_uy = conv2(u_dy2, temp.t/2); //t
  mat u_paritaldy = pow(u_dx, 2) + pow(delta_uy, 2); //uypd

  mat delta_vx = conv2(v_dy2, temp/2); //t
  mat v_paritaldx = pow(v_dx, 2) + pow(delta_vx, 2); //vxpd

  mat delta_vy = conv2(v_dy2, temp.t/2); //t
  mat v_paritaldy = pow(v_dx, 2) + pow(delta_vy, 2); //vypd


	//fill odd rows and even columns of esmooth with psifunc( uypd + vypd )
	//fill even rows and odd columns of esmooth with psifunc( uxpd + vxpd )

  //FROM MATLAB:
  //psidashFS( 1:2:end, 2:2:end ) = psiDerivative( uypd + vypd ) ;
  //psidashFS( 2:2:end, 1:2:end ) = psiDerivative( uxpd + vxpd ) ;
  //where psidashFS is e_smooth
	//also e_smooth is huge

	return e_smooth;
}
