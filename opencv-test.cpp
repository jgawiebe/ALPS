//#include "stdafx.h"

#include <opencv2/opencv.hpp>

#include <iostream>
#include <string>

using namespace cv;
using namespace std;

int main(int argc, char** argv)
{
	// Read the image file
	Mat image = imread("C:/Users/s27508/this.jpg");

	Mat output;

	resize(image, output, Size(), 2, 2, 1);


	if (image.empty()) // Check for failure
	{
		cout << "Could not open or find the image" << endl;
		system("pause"); //wait for any key press
		return -1;
	}

	String window_name1 = "Initial Image 1";
	String window_name2 = "Initial Image 2";

	namedWindow(window_name1);
	namedWindow(window_name2);

	imshow(window_name1, image1);
	imshow(window_name2, image2);

	String windowName2 = "My HelloWorld Window 2"; //Name of the window

	namedWindow(windowName2); // Create a window

	imshow(windowName2, output); // Show our image inside the created window.

	waitKey(0); // Wait for any keystroke in the window

	destroyWindow(windowName2); //destroy the created window

	return 0;
}