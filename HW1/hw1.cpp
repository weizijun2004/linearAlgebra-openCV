#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include<opencv2/opencv.hpp>
using namespace std;
void readCSV(vector <vector <double>> &, string);
void splict(vector <double> &, string, char);
// void rescaling();
// void windowing();


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
        cout << inputContext[i] << endl;
    }
    fs.close();
    vector <vector <double>> origin;
    readCSV(origin, inputContext[1]);
    // for(int i = 0;i < origin.size(); ++ i)
    // {
    //     for(int o = 0;o < origin[i].size(); ++ o) 
    //     {
    //         cout << origin[i][o] << ' ';
    //     }
    //     cout << endl;
    // }
    fs.open(inputContext[0], ios :: in);
    // cout << inputContext[0] << endl;
    if(!fs.is_open()) cout << "file cannot open, in line 40" << endl;
    string temp;
    vector <string> arr;
    stringstream ss;
    double slope, intercept;
    while(!fs.eof())
    {
        getline(fs, temp);
        //cout << "temp : " << temp << endl; 
        ss << temp;
        while(getline(ss, temp, ','))
        {
            // cout << temp << endl;
            arr.push_back(temp);
        }
        // cout << arr[0] << endl;
        if(arr[0] == "Rescale Slope") slope = stof(arr[1]);
        else if(arr[0] == "Rescale Intercept") intercept = stoi(arr[1]);
        arr.clear();
        ss.clear();
    }
    fs.close();
    // cout << "Rescale Slope : " << slope << endl;
    // cout << "Rescale Intercept : " << intercept << endl;
    for(int i = 0;i < origin.size(); ++ i)
    {
        for(int o = 0;o < origin[i].size(); ++ o) 
        {
            origin[i][o] = origin[i][o] * slope + intercept;
        }
    }
    // for(int i = 0;i < origin.size(); ++ i)
    // {
    //     for(int o = 0;o < origin[i].size(); ++ o) 
    //     {
    //         cout << origin[i][o] << ' ';
    //     }
    //     cout << endl;
    // }
    vector <vector <double>> win = origin;
    double wl = stof(inputContext[3]);
    double ww = stof(inputContext[4]);

    for(int i = 0;i < origin.size(); ++ i)
    {
        for(int o = 0;o < origin[i].size(); ++ o) 
        {
            if(origin[i][o] <= wl - ww/2) win[i][o] = 0;
            else if(origin[i][o] > wl + ww/2) win[i][o] = 255;
            else win[i][o] = 255 * ((origin[i][o] - wl) / ww + 0.5);
        }
    }

    for(int i = 0;i < win.size(); ++ i)
    {
        for(int o = 0;o < win[i].size(); ++ o) 
        {
            cout << win[i][o] << ' ';
        }
        cout << endl;
    }  
    cout << "wl : " << wl << endl;
    cout << "ww : " << ww << endl;

    fstream outputFile;
    outputFile.open(argv[2], ios :: out);
    cout << argv[2] << endl;
    for(int i = 0;i < win.size(); ++ i)
    {
        for(int o = 0;o < win[i].size(); ++ o) 
        {
            if(o != win[i].size() - 1) outputFile << win[i][o] << ',';
            else outputFile << win[i][o] << '\n';
        }
    }
    int rows = win.size();
    int cols = win[0].size();

    // 創建一個Mat對象，類型為CV_32F（單精度浮點數）
    cv::Mat image(rows, cols, CV_32F);

    // 將CSV數據轉換為Mat對象
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            image.at<float>(i, j) = win[i][j];
        }
    }
    cv::normalize(image, image, 0, 255, cv::NORM_MINMAX);

    // 將Mat對象轉換為8位無符號整數型別（CV_8U）
    cv::Mat outputImage;
    image.convertTo(outputImage, CV_8U);

    // 保存圖像為PNG檔案
    cv::imwrite("output.png", outputImage);

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


void readCSV(vector <vector <double>> &originalValue, string filePath)
{
    cout << "in" << endl;
    fstream fs;
    fs.open(filePath, ios :: in);
    
    if(!fs.is_open()) cout << "file cannot open, in line 134" << endl;
    string temp;
    vector <double> tempColumn;
    int column = 0, row;
    while(!fs.eof())
    {
        getline(fs, temp);
        splict(tempColumn, temp, ',');
        column ++ ;
        tempColumn.size() == 0 ? row = row : row = tempColumn.size();
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