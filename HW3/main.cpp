#include<iostream>
#include<fstream>
#include<vector>
#include<opencv2/opencv.hpp>
using namespace std;
using namespace cv;

bool isContourAllBlack(vector<Point>& contour, Mat& binaryImage) {
	for (Point& point : contour) {
		uchar pixelValue = binaryImage.at<uchar>(point.y, point.x);  // 注意这里的顺序

		if (pixelValue != 255) {
			// cout << "pixelValue : " << (int)pixelValue << endl;
			return false;
		}
	}
	return true;
}

int main(int argc, char **argv)
{

	Mat inputImage = imread(argv[1], 1);
	cout << inputImage.size();
	Mat grayImage;
	// to blur the image
	GaussianBlur(inputImage, inputImage, cv::Size(9, 9), 0);
	// turn image to gray
	cvtColor(inputImage, grayImage, COLOR_BGR2GRAY);
	Mat binaryImage;
	vector<vector<double>> temp;
	temp.resize(binaryImage.rows);
	for (int i = 0; i < temp.size(); ++i) temp[i].resize(binaryImage.cols);

	Mat noImage(1477, 1108, CV_8UC1, Scalar(255));
	// to Binarization the image, >100 -> 255, <100 -> 0;
	threshold(grayImage, binaryImage, 100, 255, cv::THRESH_BINARY);
	vector<vector<Point>> contours;
	// find all contours in iamge
	findContours(binaryImage, contours, cv::RETR_LIST, cv::CHAIN_APPROX_NONE);
	vector<vector<Point>> newContours;
	// turn image num type to uchar
	binaryImage.convertTo(binaryImage, CV_8UC1);
	for (int i = 0; i < contours.size(); i++) 
	{
		// find the contours that all black and also in the range of area and length 
		if (isContourAllBlack(contours[i], binaryImage) && contourArea(contours[i]) > 200 && arcLength(contours[i], true) < 200)
		{
			// cout << "all black ! " << endl;
			newContours.push_back(contours[i]);
		}
	}
	/*
	cout << "newContours size : " << newContours.size() << endl;
	for (int i = 0; i < newContours.size(); ++i) {
		cout << "area size : " << i << ' ' << contourArea(newContours[i]) << endl;
	}
	*/
	// four corners point and inisialize it to the middle
	Point upright(500, 500),
		upleft(500, 500),
		buttomleft(500, 500),
		buttomright(500, 500);
	for (int i = 0; i < newContours.size(); ++ i)
	{
		for (Point p : newContours[i])
		{
			// cout << p << endl;
			if (p.x < upleft.x && p.y < upleft.y) upleft = p;
			else if (p.x > upright.x && p.y < upright.y) upright = p;
			else if (p.x < buttomleft.x && p.y > buttomleft.y) buttomleft = p;
			else if (p.x > buttomright.x && p.y > buttomright.y) buttomright = p;
			// break;
		}
	}
	cout << newContours[0][0].x << ' ' << newContours[0][0].y << endl;
	cout << newContours[0][0] << endl;
	cout << newContours[0].size() << endl;
	cout << "----------------------------------" << endl;
	cout << upleft << endl;
	cout << upright << endl;
	cout << buttomleft << endl;
	cout << buttomright << endl;
	vector<vector<Point>> cornerContours;
	
	// find corners contours
	for (int i = 0; i < newContours.size(); ++i)
	{
		for (Point p : newContours[i])
		{
			if (p == upleft || p == upright || p == buttomleft || p == buttomright) cornerContours.push_back(newContours[i]);
		}
	}
	for (auto i : cornerContours[0]) cout << i << endl;
	drawContours(noImage, cornerContours, -1, Scalar(0, 0, 255), 1);

	// cout << "contours size : "  << contours.size() << endl;
	// cout << "square contours size : " << squareContours.size() << endl;
	imwrite(argv[2], noImage);
	waitKey(0);
}