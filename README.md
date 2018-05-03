# DESCRIPTION OF PROJECT
Algorithm for Pixel Shiftmapping (ALPS). This is Jacob Wiebe and James Dolman's fourth year computer engineering project summited to the Royal Military College of Canada in April 2018. It uses [Brox's proposed method](http://www.mia.uni-saarland.de/Publications/brox-eccv04-of.pdf) for estimating a pixel shiftmap and [Visesh Chari's MATLAB implementation](https://www.mathworks.com/matlabcentral/fileexchange/17500-high-accuracy-optical-flow). It takes two images and generates the estimated pixel shiftmap between them in the *x* and *y* domains. The solving component has been implemented for optional parallel execution using the CUDA API.

## HOW TO USE PROJECT

### Executables
The Executables folder holds a series of compiled programs which use variations of the iterative solving method. They are as follows:
* **ALPS-SOR**: uses the Successive Over-Relaxation method for solving the system of linear equations in series.
* **ALPS-JACOBI**: uses the Jacobi method for solving the system of linear equations in series.
* **PALPS-v1**: uses the Jacobi method in parallel for solving the system of linear eqautions column-wise. In the source code this uses `jacobi_serial_v1`. This is only the solving portion of the program using test_sor.
* **PALPS-v2**: uses the Jacobi method in parallel for solving the system of linear eqautions row-wise (much more efficient). In the source code this uses `jacobi_serial_v2`. This is only the solving portion of the program using test_sor. NOTE: This is an incomplete version and does not converge.
These files were compiled on Windows 7 for a 64-bit architecture. They require the included libraries to run. The PALPS executables require a CUDA compatible GPU.

### Project Files
The `project_files` folder contains the settings, libraries and additional files used while writing this program. It was developed using Visual Studio 2017.

# NOTES ABOUT STYLE
All files (apart from CUDA) use the `std` and `arma` namespaces. This implies that *mat* is always an Armadillo-type matrix and OpenCV-type matrices use `cv::Mat`. Sometimes (in the case of the bilinear resize file) `arma::mat` is used for additional clarity.

# CONTACT INFO
Jacob Wiebe jgawiebe@gmail.com

