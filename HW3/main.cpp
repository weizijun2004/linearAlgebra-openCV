#include<iostream>
#include<fstream>
#include<vector>
#include<opencv2/opencv.hpp>
using namespace std;
using namespace cv;
int main(int argc, char **argv)
{

	Mat inputImage = imread(argv[1], 1);
	cout << inputImage.size();
	imshow("test", inputImage);
	waitKey(0);


}