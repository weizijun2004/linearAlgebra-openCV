#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include<opencv2/opencv.hpp>
using namespace std;
void readCSV(vector <vector <double>> &, string);
void splict(vector <double>&, string, char); 
void windowing(vector <vector<double>> &win, vector <vector<double>> &origin, double , double );
void turnCsvToMat(vector<vector<double>>& csv, cv::Mat& mat);
void createImage(vector<vector<double>> &, int length, string path);
cv::Mat mergeImage(vector<vector<double>>& b, vector<vector<double>>& g, vector<vector<double>>& r);


int main(int argc, char *argv[])
{

    string inputPath = argv[1];
    fstream fs;
    fs.open(inputPath, ios :: in);
    if(!fs.is_open()) cout << "file cannot open, in line 20" << endl;
    string inputContext[10];
    for(int i = 0;!fs.eof(); ++ i)
    {
        fs >> inputContext[i];
        // cout << inputContext[i] << endl;
    }
    fs.close();

    vector <vector <double>> origin;
    readCSV(origin, inputContext[1]);

    fs.open(inputContext[0], ios :: in);
    if(!fs.is_open()) cout << "file cannot open, in line 40" << endl;
    string temp;
    vector <string> arr(10);
    stringstream ss;
    double slope = 0, intercept = 0;
    while(!fs.eof())
    {
        getline(fs, temp);
        ss << temp;
        while(getline(ss, temp, ',')) 
        {
            arr.push_back(temp);
        }
        if (arr.size() > 0)
        {
            if (arr[0] == "Rescale Slope") slope = stod(arr[1]);
            else if (arr[0] == "Rescale Intercept") intercept = stod(arr[1]);
        }
        arr.clear();
        ss.clear();
    }
    fs.close();
    if (origin.size() == 513) origin.pop_back();

    for(int i = 0;i < origin.size(); ++ i)
    {
        for(int o = 0;o < origin[i].size(); ++ o) 
        {
            origin[i][o] = origin[i][o] * slope + intercept;
        }
    }

    vector <vector <double>> win1 = origin;
    double wl1 = stod(inputContext[3]);
    double ww1 = stod(inputContext[4]);
    windowing(win1, origin, wl1, ww1);
    /* // pass
    for(int i = 0;i < win.size(); ++ i)
    {
        for(int o = 0;o < win[i].size(); ++ o) 
        {
            cout << win[i][o] << ' ';
        }
        cout << endl;
    }  
    */
    createImage(win1, win1.size(), argv[2]);
    // start count mask
    vector <vector <double>> mask;
    readCSV(mask, inputContext[2]);
    if (mask.size() == 513) mask.pop_back();
    for (int i = 0; i < mask.size(); ++i)
    {
        for (int o = 0; o < mask[i].size(); ++o)
        {
            mask[i][o] = floor((double)(255) / (double)(16) * mask[i][o]);
        }
    }
    vector <vector <double>> win2 = origin;
    double wl2 = stod(inputContext[5]);
    double ww2 = stod(inputContext[6]);
    windowing(win2, origin, wl2, ww2);
    // createImage(win2, win2.size(), argv[3]);
    cv::Mat merged = mergeImage(mask, win1, win2);

    // deleted all to find falx *
    vector <vector <double>> falx = mask;
    for (int i = 0; i < falx.size(); ++i)
    {
        for (int o = 0; o < falx[i].size(); ++o)
        {
            // cout << newMask[i][o] << ',';
            if (falx[i][o] != 47) falx[i][o] = 0;
            else falx[i][o] = 255;
        }
    }
    cv::Mat falxImage(falx.size(), falx.size(), CV_64F);
    turnCsvToMat(falx, falxImage);
    // to find and draw contours
    vector<vector<cv::Point>> contours;
    vector<vector<cv::Vec4d>> hierarchy;
    falxImage.convertTo(falxImage, CV_8UC1);
    cv::findContours(falxImage, contours, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE);

    // find max and second biggest perimeter shape
    int maxPerimeterTarget = 0, secondPT = 0;
    double maxPerimeter = 0, secondP = 0;

    for (int t = 0; t < contours.size(); t++)
    {
        if (arcLength(contours[t], true) >= maxPerimeter)
        {
            if (maxPerimeterTarget == 0)
            {
                maxPerimeterTarget = t;
                maxPerimeter = arcLength(contours[t], true);
            }
            else
            {
                secondPT = maxPerimeterTarget;
                secondP = maxPerimeter;
                maxPerimeterTarget = t;
                maxPerimeter = arcLength(contours[t], true);
            }
            // cout << "t : " << t << endl;
            // cout << "maxPerimeter : " << maxPerimeter << endl;
        }
        else if (arcLength(contours[t], true) > secondP)
        {
            secondPT = t;
            secondP = arcLength(contours[t], true);
        }
        // cout << "the " << t << " length = " << arcLength(contours[t], true) << endl;
    }
    vector<vector<cv::Point>> targetContours;
    // cout << "maxPerimeter : " << maxPerimeter << endl;
    targetContours.push_back(contours[maxPerimeterTarget]);
    // if falx is complete, not record other noise
    if(maxPerimeter < 550) targetContours.push_back(contours[secondPT]);

    cv::Mat imageNoiseDeleted(falx.size(), falx.size(), CV_8UC1);
    cv::drawContours(imageNoiseDeleted, targetContours, -1, cv::Scalar::all(255), 1);

    for (int i = 0; i < falx.size(); ++i)
    {
        for (int o = 0; o < falx.size(); ++o)
        {
            if ((int)(imageNoiseDeleted.at<uchar>(i, o)) == 205) imageNoiseDeleted.at<uchar>(i, o) = 0;
        }
    }
    // cout << targetContours[0] << endl;
    vector<cv::Point> lineContours;
    for (int i = 0; i < targetContours.size(); ++i)
    {
        for (int o = 0; o < targetContours[i].size(); ++o) lineContours.push_back(targetContours[i][o]);
    }
    ///*
    cv::Vec4d linePoint;
    cv::fitLine(lineContours, linePoint, 1, 0, 0, 0);
    // cout << "linePoint : " << linePoint << endl;

    // find line
    cv::Point point0;
    point0.x = linePoint[2];
    point0.y = linePoint[3];
    double lineSlope = linePoint[1] / linePoint[0];
    cout << "lineSlope : " << lineSlope << endl;
    if (fabs(lineSlope) <= 1e-8) lineSlope = 0;

    double angle = (double)(atan(lineSlope)) * (double)(180) / (double)(CV_PI);
    // cout << "angle first : " << angle << endl;
    if (angle < 0) angle += 90;
    else angle -= 90;
    cout << "angle last : " << angle << endl;
    cv::Point point1(0, (lineSlope * (0 - point0.x) + point0.y)), point2(640, (lineSlope * (640 - point0.x) + point0.y));
    cv::line(falxImage, point1, point2, cv::Scalar(255, 0, 0), 2);
    // rotate image
    cv::Point center = cv::Point(falxImage.cols / 2, falxImage.rows / 2);
    cv::Mat rotMat = cv::getRotationMatrix2D(center, angle, 1);
    cv::Mat rotationedImage;
    cv::warpAffine(merged, rotationedImage, rotMat, falxImage.size());


    //*/
    bool success = cv::imwrite(argv[3], rotationedImage);
    if (success)
    {
        cout << "image has been save in " << argv[3] << endl;
    }
    else
    {
        cerr << "error in save " << argv[3] << "image" << endl;
    }

    fs.open(argv[4], ios::out);
    fs << fixed << setprecision(6) << angle;
    fs.close();
    /*
    for (int i = 0; i < line.size(); ++i)
    {
        for (int o = 0; o < line[i].size(); ++o)
        {
            // cout << newMask[i][o] << ',';
            cout << (int)(imagetest.at<uchar>(i, o)) << ',';
        }
        cout << endl;
    }
    */
    return 0;
}
///*

cv::Mat mergeImage(vector<vector<double>>&b, vector<vector<double>>&g, vector<vector<double>>&r)
{
    cv :: Mat bImage(b.size(), b.size(), CV_64F);
    cv :: Mat gImage(g.size(), g.size(), CV_64F);
    cv :: Mat rImage(r.size(), r.size(), CV_64F);
    turnCsvToMat(b, bImage);
    turnCsvToMat(g, gImage);
    turnCsvToMat(r, rImage);
    bImage.convertTo(bImage, CV_8U);
    gImage.convertTo(gImage, CV_8U);
    rImage.convertTo(rImage, CV_8U);
    vector<cv :: Mat> channels = { bImage, gImage, rImage};
    cv :: Mat mergedImage;
    merge(channels, mergedImage);
    /*
    bool success = cv::imwrite(path, mergedImage);
    if (success)
    {
        cout << "image has been save in " << path << endl;
    }
    else
    {
        cerr << "error in save " << path << "image" << endl;
    }
    */
    return mergedImage;
}

void createImage(vector<vector<double>> &csv, int length, string path)
{
    // cout << "win size : " << csv.size() << endl;
    // 創建一個Mat對象，類型為CV_64F（雙精度浮點數）
    cv::Mat image(length, length, CV_64F);

    turnCsvToMat(csv, image);
    image.convertTo(image, CV_8U);
    bool success = cv::imwrite(path, image);
    if (success) 
    {
        cout << "image has been save in " << path << endl;
    }
    else 
    {
        cerr << "error in save " << path << "image" << endl;
    }
}
//*/
void turnCsvToMat(vector<vector<double>>& csv, cv::Mat& mat)
{
    for (int i = 0; i < csv.size(); ++i)
    {
        for (int j = 0; j < csv[i].size(); j++)
        {
            mat.at<double>(i, j) = floor(csv[i][j]);
        }
        for (int o = csv[i].size(); o < csv.size(); ++o)
        {
            mat.at<double>(i, o) = 0;
        }

        // cout << "row : " << i << endl;
        // cout << endl;
    }
}
void splict(vector <double> &tempColumn, string temp, char symbol)
{
    stringstream ss(temp);
    while(getline(ss, temp, symbol))
    {
        //cout << temp << ',';
        tempColumn.push_back(stof(temp));
    }
}

void windowing(vector <vector<double>> &win, vector <vector<double>> &origin, double wl, double ww)
{
    for (int i = 0; i < win.size(); ++i)
    {
        for (int o = 0; o < win[i].size(); ++o)
        {
            if (win[i][o] <= wl - ww / 2) win[i][o] = 0;
            else if (win[i][o] > wl + ww / 2) win[i][o] = 255;
            else win[i][o] = 255 * ((origin[i][o] - wl) / ww + 0.5);
        }
    }
}

void readCSV(vector <vector <double>> &originalValue, string filePath)
{
    // cout << "in" << endl;
    fstream fs;
    fs.open(filePath, ios :: in);
    
    if(!fs.is_open()) cout << "file cannot open, in line 134" << endl;
    string temp;
    vector <double> tempColumn;
    int column = 0, row = 0;
    while(!fs.eof())
    {
        
        getline(fs, temp);
        
        splict(tempColumn, temp, ',');
        column ++ ;
        if(tempColumn.size() != 0) row = tempColumn.size();
        originalValue.push_back(tempColumn);
        tempColumn.clear();
    }
    // cout << column << ' ' << row << endl;
    // for(int i = 0;i < column; ++ i)
    // {
    //     for(int o = 0;o < originalValue[i].size(); ++ o) 
    //     {
    //         cout << originalValue[i][o] << ' ';
    //     }
    //     cout << endl;
    // }
    fs.close();
}