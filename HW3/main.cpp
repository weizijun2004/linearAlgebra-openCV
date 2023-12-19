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
	Mat grayImage;
	GaussianBlur(inputImage, inputImage, cv::Size(5, 5), 0);
	cvtColor(inputImage, grayImage, COLOR_BGR2GRAY);
	Mat binaryImage;
	vector<vector<double>> temp;
	temp.resize(binaryImage.rows);
	for (int i = 0; i < temp.size(); ++i) temp[i].resize(binaryImage.cols);

	Mat noImage(1477, 1108, CV_8UC1, Scalar(255));
	// threshold(grayImage, noImage, 0, 255, cv::THRESH_BINARY);
	threshold(grayImage, binaryImage, 130, 255, cv::THRESH_BINARY);
	vector<vector<Point>> contours;
	vector<Vec4i> point;
	findContours(binaryImage, contours, cv::RETR_LIST, cv::CHAIN_APPROX_NONE);

	vector<vector<Point>> newContours;
	for (int i = 0; i < contours.size(); ++i) {
		if (contourArea(contours[i]) > 200 && arcLength(contours[i], false) < 200)
		{
			newContours.push_back(contours[i]);
			cout << i << " : " << contourArea(contours[i]) << endl;
		}
	}
	vector<Point> min_x_max_y, max_x_max_y, min_x_min_y, max_x_min_y;

	// ﹍て硂ㄇ跑计絋玂Τì镑翴秈︽ゑ耕
	min_x_max_y = newContours[0];
	min_x_min_y = newContours[0];
	max_x_min_y = newContours[0];
	max_x_max_y = newContours[0];

	// 筂菌┮Τcontoursт才兵ン翴
	for (const auto& contour : newContours) {
		for (const auto& point : contour) {
			// тx禸程y禸程翴
			if (point.x <= min_x_max_y[0].x && point.y >= min_x_max_y[0].y) {
				min_x_max_y = contour;
				cout << "0" << endl;
			}

			// тx禸程y禸程翴
			if (point.x >= max_x_max_y[0].x && point.y >= max_x_max_y[0].y) {
				max_x_max_y = contour;
				cout << "1" << endl;
			}

			// тx禸程y禸程翴
			if (point.x <= min_x_min_y[0].x && point.y <= min_x_min_y[0].y) {
				min_x_min_y = contour;
				cout << "2" << endl;
			}

			// тx禸程y禸程翴
			if (point.x >= max_x_min_y[0].x && point.y <= max_x_min_y[0].y) {
				max_x_min_y = contour;
				cout << "3" << endl;
			}
		}
	}
	vector<vector<Point>> squareContours;
	squareContours.push_back(max_x_max_y);
	squareContours.push_back(max_x_min_y);
	squareContours.push_back(min_x_max_y);
	squareContours.push_back(min_x_min_y);
	for(int i = 0;i < squareContours.size(); ++ i) cout << i << " : " << contourArea(squareContours[i]) << endl;
	cout << "size : " << noImage.size() << endl;
	// Mat resultImage;
	drawContours(noImage, squareContours, -1, Scalar(0, 0, 255), 1);
	// cout << "contours size : "  << contours.size() << endl;
	// cout << "square contours size : " << squareContours.size() << endl;
	imwrite(argv[2], noImage);
	waitKey(0);
}