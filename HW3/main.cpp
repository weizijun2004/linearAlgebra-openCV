#include<iostream>
#include<fstream>
#include<vector>
#include<algorithm>
#include<cmath>
#include<opencv2/opencv.hpp>
using namespace std;
using namespace cv;

bool isContourAllBlack(vector<Point>& contour, Mat& binaryImage) 
{

	for (Point& point : contour) {
		uchar pixelValue = binaryImage.at<uchar>(point.y, point.x);  // 注意这里的顺序

		if (pixelValue != 255) {
			// cout << "pixelValue : " << (int)pixelValue << endl;
			return false;
		}
	}
	return true;
}

///*
int whichChoose(vector<Point> ct)
{
	for (Point p : ct)
	{
		if (p.x == 275)return 1;
		else if (p.x == 340)return 2;
		else if (p.x == 400)return 3;
		else if (p.x == 460)return 4;
		else if (p.x == 520)return 5;
		else if (p.x == 580)return 6;
		else if (p.x == 640)return 7;
		else if (p.x == 700)return 8;
		else if (p.x == 765)return 9;
		else if (p.x == 825)return 10;
		else if (p.x == 890)return 11;
		else if (p.x == 950)return 12;
	}
}
//*/
	/*
	1 -> x = 275
	2 -> x = 340
	3 -> x = 400
	4 -> x = 460
	5 -> x = 520
	6 -> x = 580
	7 -> x = 640
	8 -> x = 700
	9 -> x = 765
	10 -> x = 825
	11 -> x = 890
	12 -> x = 950
	*/
bool inXRange(vector<Point>& contour)
{
	for (Point& point : contour) {
		if (point.x == 275 || point.x == 340 || point.x == 400 || point.x == 460 || point.x == 520 || point.x == 580
			|| point.x == 640 || point.x == 700 || point.x == 765 || point.x == 825 || point.x == 890 || point.x == 950)
		{
			return true;
		}
	}
	return false;
}

int main(int argc, char **argv)
{

	Mat inputImage = imread(argv[1], 1);
	cout << inputImage.size();
	Mat blurImage, grayImage;
	// to blur the image
	GaussianBlur(inputImage, blurImage, cv::Size(9, 9), 0);
	// turn image to gray
	cvtColor(blurImage, grayImage, COLOR_BGR2GRAY);
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
		if (contourArea(contours[i]) > 100 && arcLength(contours[i], true) < 100)
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
	Point upright(500, 500), upleft(500, 500), buttomleft(500, 500), buttomright(500, 500);
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
	cout << inputImage.size() << endl;
	// to optimize the corner point locate
	for (auto i : cornerContours)
	{
		double x, y;
		for (Point p : i)
		{
			if (p == upleft)
			{
				auto minElement = min_element(i.begin(), i.end(),
					[](const Point& p1, const Point& p2) {
						return (p1.x < p2.x) || ((p1.x == p2.x) && (p1.y < p2.y));
					});

				// get upleft point
				upleft = *minElement;
			}
			else if (p == upright)
			{
				auto maxElement = max_element(i.begin(), i.end(),
					[](const Point& p1, const Point& p2) {
						return (p1.x < p2.x) && (p1.y == p2.y) || (p1.y > p2.y);
					});

				// get upleft point
				upright = *maxElement;
			}
			else if (p == buttomleft)
			{
				auto maxElement = max_element(i.begin(), i.end(),
					[](const Point& p1, const Point& p2) {
						return (p1.x > p2.x) && (p1.y == p2.y) || (p1.y < p2.y);
					});

				// get upleft point
				buttomleft = *maxElement;
			}
			else if (p == buttomright)
			{
				auto maxElement = max_element(i.begin(), i.end(),
					[](const Point& p1, const Point& p2) {
						return (p1.x < p2.x) || ((p1.x == p2.x) && (p1.y < p2.y));
					});
				
				// get upleft point
				buttomright = *maxElement;
			}
		}
	}

	vector<Point2f> cornerPoint;
	cornerPoint.push_back(upleft);
	cornerPoint.push_back(upright);
	cornerPoint.push_back(buttomleft);
	cornerPoint.push_back(buttomright);
	vector<Point2f> detPoint;
	detPoint.push_back(Point2f(0, 0));
	detPoint.push_back(Point2f(1107, 0));
	detPoint.push_back(Point2f(0, 1476));
	detPoint.push_back(Point2f(1107, 1476));

	Mat rotaMat = getPerspectiveTransform(cornerPoint, detPoint);
	// for (auto i : cornerPoint) circle(inputImage, i, 2, Scalar(0, 0, 255), -1);

	warpPerspective(inputImage, inputImage, rotaMat, inputImage.size());
	contours.clear();
	newContours.clear();
	// drawContours(noImage, cornerContours, -1, Scalar(0, 0, 255), 1);
	GaussianBlur(inputImage, blurImage, cv::Size(13, 13), 0);
	cvtColor(blurImage, grayImage, COLOR_BGR2GRAY);
	threshold(grayImage, binaryImage, 120, 255, cv::THRESH_BINARY);
	findContours(binaryImage, contours, cv::RETR_LIST, cv::CHAIN_APPROX_NONE);
	circle(inputImage, Point(235, 540), 4, Scalar(0, 0, 255), -1);
	circle(inputImage, Point(970, 540), 4, Scalar(0, 0, 255), -1);
	circle(inputImage, Point(235, 1300), 4, Scalar(0, 0, 255), -1);
	circle(inputImage, Point(275, 615), 4, Scalar(0, 0, 255), -1);
	circle(inputImage, Point(340, 615), 4, Scalar(0, 0, 255), -1);
	circle(inputImage, Point(400, 615), 4, Scalar(0, 0, 255), -1);
	circle(inputImage, Point(460, 615), 4, Scalar(0, 0, 255), -1);
	circle(inputImage, Point(520, 615), 4, Scalar(0, 0, 255), -1);
	circle(inputImage, Point(580, 615), 4, Scalar(0, 0, 255), -1);
	circle(inputImage, Point(640, 615), 4, Scalar(0, 0, 255), -1);
	circle(inputImage, Point(700, 615), 4, Scalar(0, 0, 255), -1);
	circle(inputImage, Point(765, 615), 4, Scalar(0, 0, 255), -1);
	circle(inputImage, Point(825, 615), 4, Scalar(0, 0, 255), -1);
	circle(inputImage, Point(890, 615), 4, Scalar(0, 0, 255), -1);
	circle(inputImage, Point(950, 615), 4, Scalar(0, 0, 255), -1);
	/*
	1 -> x = 275
	2 -> x = 340
	3 -> x = 400
	4 -> x = 460
	5 -> x = 520
	6 -> x = 580
	7 -> x = 640
	8 -> x = 700
	9 -> x = 765
	10 -> x = 825
	11 -> x = 890
	12 -> x = 950
	*/
	vector<vector<Point>> edgeContours;
	for (auto ct : contours)
	{
		for (Point p : ct)
		{
			if (p.x == 105 && p.y > 550 && arcLength(ct, true) < 100)
			{
				edgeContours.push_back(ct);
				break;
			}
		}
	}
	
	vector<vector<Point>> ansContours;
	for (auto ct : contours)
	{
		for (Point p : ct)
		{
			if (p.x > 225 && p.y > 550 && p.y < 1350 && arcLength(ct, true) < 95
				&& contourArea(ct) > 55 && inXRange(ct) && isContourAllBlack(ct, binaryImage))
			{
				ansContours.push_back(ct);
				break;
			}
		}
	}
	cout << "ans size : " << ansContours.size() << endl;
	sort(edgeContours.begin(), edgeContours.end(), [](vector<Point> a, vector<Point> b) { return a[0].y < b[0].y; });
	vector<Point> maxPointArr, minPointArr;
	for (auto ct : edgeContours)
	{
		Point maxPoint(0, 0);
		Point minPoint(ct[0].x, ct[0].y);
		for (Point p : ct)
		{
			if (p.x > maxPoint.x) maxPoint.x = p.x;
			if (p.y > maxPoint.y) maxPoint.y = p.y;
			if (p.x < minPoint.x) minPoint.x = p.x;
			if (p.y < minPoint.y) minPoint.y = p.y;
		}
		maxPointArr.push_back(maxPoint);
		minPointArr.push_back(minPoint);
		cout << "min Y : " << minPoint.y << endl;
		cout << "max Y : " << maxPoint.y << endl;
		cout << endl;
	}
	cout << "edge size : " << edgeContours.size() << endl;
	vector<int> ans;
	vector<int> chooseArr;
	vector<Point> contoursShape;
	vector<vector<Point>> contoursShapeArr;
	// sort(ansContours.begin(), ansContours.end(), [](vector<Point> a, vector<Point> b) { return a[0].y < b[0].y; });
	for (int i = 0;i < edgeContours.size(); ++ i)
	{
		bool multiChoose = false;
		int temp = 0;
		for (auto ct : ansContours)
		{
			for (auto p : ct)
			{
				if (p.y > minPointArr[i].y && p.y < maxPointArr[i].y)
				{
					// if ((double)(min(p.y - minPointArr[i].y, maxPointArr[i].y - p.y)) / (double)(max(p.y - minPointArr[i].y, maxPointArr[i].y - p.y)) <= 0.8) continue;
					// cout << (double)(min(p.y - minPointArr[i].y, maxPointArr[i].y - p.y)) / (double)(max(p.y - minPointArr[i].y, maxPointArr[i].y - p.y)) << endl;
					approxPolyDP(ct, contoursShape, 2, true);
					contoursShapeArr.push_back(contoursShape);
					cout << "edge num : " << contoursShape.size() << endl;
					if (!isContourConvex(contoursShape) && contoursShape.size() > 2 ) break;
					temp++;
					if(!multiChoose) chooseArr.push_back(whichChoose(ct));
					multiChoose = true;
					break;
				}
			}
		}
		if (temp == 0) chooseArr.push_back(-1);
		cout << i + 1 << " : " << temp << endl;
		cout << "ans : " << chooseArr[chooseArr.size() - 1] << endl;
		ans.push_back(temp);
	}
	/*
	1 -> x = 275
	2 -> x = 340
	3 -> x = 400
	4 -> x = 460
	5 -> x = 520
	6 -> x = 580
	7 -> x = 640
	8 -> x = 700
	9 -> x = 765
	10 -> x = 825
	11 -> x = 890
	12 -> x = 950
	*/
	// for (int output : ans) cout << "ans : " << output << endl;

	cout << "choose arr size : " << chooseArr.size() << endl;
	cvtColor(binaryImage, binaryImage, COLOR_BayerBG2BGR);
	drawContours(binaryImage, ansContours, -1, Scalar(0, 255, 0), 1);
	// for(int i = 0, num = 600;i < 10; ++ i, num += 30) circle(inputImage, Point(100, num + i), 10, Scalar(0, 0, 255), -1);
	imwrite(argv[3], binaryImage);
	waitKey(0);
	fstream fs;
	fs.open(argv[2], ios::out);
	for (int i = 0; i < 24; ++ i)
	{
		if (ans[i] == 0) fs << 'X';
		else if (ans[i] > 1) fs << 'M';
		else
		{
			if (chooseArr[i] == 10) fs << 0;
			else if (chooseArr[i] == 11) fs << 'A';
			else if (chooseArr[i] == 12) fs << 'B';
			else fs << chooseArr[i];
		}
	}
	fs << "\n";
	fs.close();
}