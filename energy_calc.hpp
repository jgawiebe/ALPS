/*
energy_calc.hpp
Jacob Wiebe & James Dolman
Reu1: Nov 2017
Rev2: Jan 2018
*/

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

template <class T> //probably don't need template
T psi_function(T x)
{

  //M: psiDerivative
  const double epsilon = 0.001;
  T e_psi = 1 / (2 * sqrt(x + epsilon));

  return e_psi;
}

//M: computePsidashFS_brox
mat generate_esmooth(mat u, mat v)
{
  double height = u.n_rows;
  double width = u.n_cols;
  mat e_smooth(2 * height + 1, 2 * width + 1, fill::zeros);
  mat temp(1, 2);
  mat tr_temp(1, 2); //transposition of temp matrix

  temp(0) = 1;
  temp(1) = -1;

  mat u_dx = conv2(u, temp);
  mat v_dx = conv2(v, temp);

  mat u_dy = conv2(u, tr_temp);
  mat v_dy = conv2(v, tr_temp);

  temp(1) = 1;

  mat u_dx2 = conv2(u_dx, temp / 2, "valid");
  mat v_dx2 = conv2(u_dx, temp / 2, "valid");

  mat u_dy2 = conv2(u_dy, tr_temp / 2, "valid");
  mat v_dy2 = conv2(u_dy, tr_temp / 2, "valid");

  mat delta_ux = conv2(u_dy2, temp / 2);       //t
  mat u_pdx = pow(u_dx, 2) + pow(delta_ux, 2); //uxpd

  mat delta_uy = conv2(u_dy2, tr_temp / 2);     //t
  mat u_pdy = pow(u_dx, 2) + pow(delta_uy, 2); //uypd

  mat delta_vx = conv2(v_dy2, temp / 2);       //t
  mat v_pdx = pow(v_dx, 2) + pow(delta_vx, 2); //vxpd

  mat delta_vy = conv2(v_dy2, tr_temp / 2);     //t
  mat v_pdy = pow(v_dx, 2) + pow(delta_vy, 2); //vypd

  for (int col_ix = 0; col_ix < e_smooth.n_cols; col_ix++)
  {
    if (col_ix % 2 == 0)
    { //even column
      for (int row_ix = 1; row_ix < e_smooth.n_rows; row_ix += 2)
      { //increment through odd rows
        e_smooth(row_ix, col_ix) = psi_function(u_pdy + v_pdy);
      }
    }
    else
    { //odd column
      for (int row_ix = 0; row_ix < e_smooth.n_rows; row_ix += 2)
      { //increment through even rows
        e_smooth(row_ix, col_ix) = psi_function(u_pdx + v_pdx);
      }
    }
  }

  //FROM MATLAB:
  //psidashFS( 1:2:end, 2:2:end ) = psiDerivative( uypd + vypd ) ;
  //psidashFS( 2:2:end, 1:2:end ) = psiDerivative( uxpd + vxpd ) ;
  //where psidashFS is e_smooth
  //also e_smooth is huge

  return e_smooth;
}
