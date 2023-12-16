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
	vector<vector<Point>> contours;
	vector<Vec4i> point;
	Mat grayImage;
	// GaussianBlur(inputImage, inputImage, cv::Size(5, 5), 0);
	cvtColor(inputImage, grayImage, COLOR_BGR2GRAY);
	Mat binaryImage;
	threshold(grayImage, binaryImage, 130, 255, cv::THRESH_BINARY);

	imshow("test", binaryImage);
	waitKey(0);
}