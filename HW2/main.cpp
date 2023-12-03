#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <stack>
#include <cmath>
using namespace std;

class Matrix
{
public:
    Matrix(vector<vector<double>> matTemp) { matrix = matTemp; }
    Matrix(double matInput[4][4])
    {
        vector<double> vecTemp;
        for (int i = 0; i < 4; ++i)
        {
            for (int q = 0; q < 4; ++q)
            {
                vecTemp.push_back(matInput[i][q]);
            }
            matrix.push_back(vecTemp);
            vecTemp.clear();
        }
    }
    ~Matrix() { matrix.clear(); }
    int row() { return matrix[0].size(); }
    int col() { return matrix.size(); }
    void print()
    {
        for (int i = 0; i < matrix.size(); ++i)
        {
            for (int q = 0; q < matrix[i].size(); ++q)
                cout << fixed << setprecision(2) << matrix[i][q] << ' ';
            cout << endl;
        }
        cout << endl;
    }
    void writeFile(fstream &fs)
    {
        for(int m = 0;m < matrix.size(); ++ m)
        {
            for(int n = 0;n < matrix[m].size(); ++ n)
            {
                if(n == matrix[m].size() - 1)fs << fixed << setprecision(2) << matrix[m][n];
                else fs << fixed << setprecision(2) << matrix[m][n] << ' ';
            }
            fs << endl;
        }
    }
    vector<vector<double>> getMat(){ return matrix; }
    vector<double> operator*(vector<double> &vec)
    {
        try
        {
            if (vec.size() != this->row())
            {
                cout << "vec size : " << vec.size() << endl;
                cout << "matrix size : " << this->row() << endl;
                throw "not right size ! ";
            }
        }
        catch (string exp)
        {
            cerr << exp << endl;
            exit(1);
        }
        Matrix matTemp(matrix);
        vector<double> countTemp;
        double numTemp = 0;
        for (int i = 0; i < matTemp.matrix.size(); ++i)
        {
            for (int c = 0; c < vec.size(); ++c)
                numTemp += matTemp.matrix[i][c] * vec[c];
            countTemp.push_back(numTemp);
            numTemp = 0;
            // matTemp.matrix[i] = countTemp;
            // countTemp.clear();
        }
        return countTemp;
    }
    Matrix operator*(vector<vector<double>> &mat)
    {
        try
        {
            if (mat.size() != this->row())
                throw "not right size ! ";
        }
        catch (string exp)
        {
            cerr << exp << endl;
            exit(1);
        }
        Matrix matTemp(matrix);
        vector<double> countTemp;
        double numTemp = 0;
        for (int i = 0; i < matTemp.matrix.size(); ++i)
        {
            for (int c = 0; c < mat[0].size(); ++c)
            {
                for (int j = 0; j < mat.size(); ++j)
                    numTemp += matTemp.matrix[i][j] * mat[j][c];
                countTemp.push_back(numTemp);
                numTemp = 0;
            }
            matTemp.matrix[i] = countTemp;
            countTemp.clear();
        }
        return matTemp;
    }
    Matrix operator*(Matrix &mat)
    {
        try
        {
            if (mat.matrix.size() != this->row())
                throw "not right size ! ";
        }
        catch (string exp)
        {
            cerr << exp << endl;
            exit(1);
        }
        Matrix matTemp(matrix);
        vector<double> countTemp;
        double numTemp = 0;
        for (int i = 0; i < matTemp.matrix.size(); ++i)
        {
            for (int c = 0; c < mat.matrix[0].size(); ++c)
            {
                for (int j = 0; j < mat.matrix.size(); ++j)
                    numTemp += matTemp.matrix[i][j] * mat.matrix[j][c];
                countTemp.push_back(numTemp);
                numTemp = 0;
            }
            matTemp.matrix[i] = countTemp;
            countTemp.clear();
        }
        return matTemp;
    }

private:
    vector<vector<double>> matrix;
};

void Trans(double a, double b, double c, vector<Matrix> &matStorage)
{
    double matTemp[4][4] = {
            {1, 0, 0, a},
            {0, 1, 0, b},
            {0, 0, 1, c},
            {0, 0, 0, 1}
        };
    Matrix m(matTemp);
    matStorage.push_back(m);
}
void Scaling(double a, double b, double c, double x, double y, double z, vector<Matrix> &matStorage)
{
    double matTemp[4][4] = {
        {a, 0, 0, x},
        {0, b, 0, y},
        {0, 0, c, z},
        {0, 0, 0, 1}
    };
    // for (int i = 0; i < 3; ++i)
    //     matTemp[i][i] = stod(inputTemp[i + 4]);
    // for (int i = 0; i < 3; ++i)
    //     matTemp[i][3] = stod(inputTemp[i + 1]);
    Matrix m(matTemp);
    matStorage.push_back(m);
}
void Projection(string type, vector<Matrix> &matStorage)
{
    double matTemp[4][4] = {
            {1, 0, 0, 0},
            {0, 1, 0, 0},
            {0, 0, 1, 0},
            {0, 0, 0, 1}
        };
        if (type == "xy")
            matTemp[2][2] = 0;
        else if (type == "xz")
            matTemp[1][1] = 0;
        else if (type == "yz")
            matTemp[0][0] = 0;
        else
            cout << "Error, wrong P type of Matrix ! " << endl;
        Matrix m(matTemp);
        matStorage.push_back(m);
}
void Rotation(char type, double x, double y, double z, double angle_radians, vector<Matrix> &matStorage)
{
    // cout << "in R : " << type << endl;
    double c = cos(angle_radians), s = sin(angle_radians);
    if (type == 'x')
        {
            double matTemp[4][4] = {
                {1, 0, 0, x},
                {0, c, -1 * s, y},
                {0, s, c, z},
                {0, 0, 0, 1}
            };
            Matrix m(matTemp);
            matStorage.push_back(m);
        }
        else if (type == 'y')
        {
            // cout << "in Ry" << endl;
            double matTemp[4][4] = {
                {c, 0, s, x},
                {0, 1, 0, y},
                {s * -1, 0, c, z},
                {0, 0, 0, 1}
            };
            Matrix m(matTemp);
            matStorage.push_back(m);
        }
        else if (type == 'z')
        {
            double matTemp[4][4] = {
                {c, s * -1, 0, x},
                {s, c, 0, y},
                {0, 0, 1, z},
                {0, 0, 0, 1}
            };
            Matrix m(matTemp);
            matStorage.push_back(m);
        }
        else
            cout << "Error, wrong P type of Matrix ! " << endl;
}
void Shearing(char type, double x, double y, double z, double s, double t, vector<Matrix> &matStorage)
{
    if (type == 'x')
        {
            double matTemp[4][4] = {
                {1, s, t, x},
                {0, 1, 0, y},
                {0, 0, 1, z},
                {0, 0, 0, 1}
            };
            Matrix m(matTemp);
            matStorage.push_back(m);
        }
        else if (type == 'y')
        {
            double matTemp[4][4] = {
                {1, 0, 0, x},
                {s, 1, t, y},
                {0, 0, 1, z},
                {0, 0, 0, 1}
            };
            Matrix m(matTemp);
            matStorage.push_back(m);
        }
        else if (type == 'z')
        {
            double matTemp[4][4] = {
                {1, 0, 0, x},
                {0, 1, 0, y},
                {s, t, 1, z},
                {0, 0, 0, 1}
            };
            Matrix m(matTemp);
            matStorage.push_back(m);
        }
}

void readMatrix(string input, vector<Matrix> &matStorage)
{
    stringstream ss;
    // vector<vector<double>> mat;
    // vector<double> vecTemp;
    vector<string> inputTemp;
    ss << input;
    while (getline(ss, input, ' '))
        inputTemp.push_back(input);
    ss.clear();
    if (inputTemp[0][1] == 'T')
    {
        Trans(stod(inputTemp[1]), stod(inputTemp[2]), stod(inputTemp[3]), matStorage);
    }
    else if (inputTemp[0][1] == 'P')
    {
        string flat = "";
        flat += inputTemp[0][2];
        flat += inputTemp[0][3];
        Projection(flat, matStorage);
    }
    else if (inputTemp[0][1] == 'S')
    {
        Trans(stod(inputTemp[1]) * -1, stod(inputTemp[2]) * -1, stod(inputTemp[3]) * -1, matStorage);
        Scaling(
            stod(inputTemp[4]), stod(inputTemp[5]), stod(inputTemp[6]), 
            stod(inputTemp[1]), stod(inputTemp[2]), stod(inputTemp[3]),
            matStorage
            );
    }
    else if (inputTemp[0][1] == 'R')
    {
        // cout << "inside R" << endl;
        double a = stod(inputTemp[1]), b = stod(inputTemp[2]), c = stod(inputTemp[3]);
        double angle_degrees = stod(inputTemp[4]);
        double angle_radians = angle_degrees * M_PI / 180.0;
        Trans(a * -1, b * -1, c * -1, matStorage);
        Rotation(inputTemp[0][2], a, b, c, angle_radians, matStorage);
    }

    else if (inputTemp[0][1] == 'H')
    {
        double a = stod(inputTemp[1]), b = stod(inputTemp[2]), c = stod(inputTemp[3]);
        double s = stod(inputTemp[4]), t = stod(inputTemp[5]);
        Trans(a * -1, b * -1, c * -1, matStorage);
        Shearing(inputTemp[0][2], a, b, c, s, t, matStorage);
    }
    else
        cerr << "error Matrix type ! " << endl;
    
}


double det(Matrix matrix)
{
    double ans = 0;
    vector<vector<double>> mat(matrix.getMat());
    int size = mat.size();
    if(size == 2) return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];

    for(int i = 0;i < size; ++ i)
    {
        vector<vector<double>> subMat(size - 1, vector<double>(size - 1));
        for(int q = 1;q < size; ++ q)
        {
            for(int p = 0;p < size; ++ p)
            {
                if(p < i) subMat[q - 1][p] = mat[q][p];
                else if(p > i) subMat[q - 1][p - 1] = mat[q][p];
            }
        }
        ans += (i % 2 == 0 ? 1 : -1) * mat[0][i] * det(Matrix(subMat));
    }
    return ans;
}

vector<vector<double>> pointToVec(vector<vector<double>> vec)
{
    for(int i = 0;i < vec.size(); ++ i) vec[i].pop_back();
    vector<vector<double>> ans(vec.size() - 1, vector<double>(vec[0].size()));
    for(int i = 1;i < vec.size(); ++ i)
    {
        for(int q = 0;q < vec[i].size(); ++ q) ans[i - 1][q] = vec[i][q] - vec[0][q];
    }
    return ans;
}
void Gaus(vector<double> &a, double num, vector<double> &b)
{
    if(a.size() != b.size()) cout << "WTF ?" << endl;
    for(int i = 0;i < b.size(); ++ i)
    {
        b[i] += a[i] * num;
    }
}

vector<vector<double>> inverseMat(vector<vector<double>> mat)
{
    vector<vector<double>> ans = {
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    };
    // Matrix(ans).print();
    for(int i = 0;i < mat[0].size(); ++ i)
    {
        double devideTemp = mat[i][i];
        // cout << "i : " << i << endl;
        if(mat[i][i] != 1)
        {
            for(int q = 0;q < mat[i].size(); ++ q) mat[i][q] /= devideTemp;
            for(int q = 0;q < ans[i].size(); ++ q) ans[i][q] /= devideTemp;
        }
        for(int q = 0;q < mat.size(); ++ q)
        {
            if(q != i)
            {
                double countTemp = mat[q][i];
                if(countTemp != 0) 
                {
                    Gaus(mat[i], countTemp * -1, mat[q]);
                    Gaus(ans[i], countTemp * -1, ans[q]);
                }
                /*
                cout << "mat : " << devideTemp << endl;
                Matrix(mat).print();
                cout << "ans : " << endl;
                Matrix(ans).print();
                */
            }
        }
    }
    return ans;
}


int main(int argc, char **argv)
{
    fstream fs;
    // fs.open(argv[1], ios :: in);
    // fs.open("/home/zj/大二上-作業文件/code/linearAlgebra/HW2/case1/input1.txt", ios ::in);
    // read files
    fs.open("/home/weizijun/linearAlgebra-openCV/HW2/case3/input1.txt", ios :: in);
    string inputTemp;
    vector<string> input;
    while (getline(fs, inputTemp))
    {
        input.push_back(inputTemp);
    }
    fs.close();
    // turn string Vector to double and store in array
    vector<vector<double>> vecStorage;
    stringstream ss;
    vector<double> vecTemp;
    for (int i = 0; i < 5; ++i)
    {
        string temp = input[i];
        ss << temp;
        while (getline(ss, temp, ' ')) vecTemp.push_back(stod(temp));
        vecTemp.push_back(1);
        ss.clear();
        vecStorage.push_back(vecTemp);
        vecTemp.clear();
    }
    // get the vector that being reverse
    vector<double> reverseVec(vecStorage[vecStorage.size() - 1]);
    vecStorage.pop_back();
    // read Matrix type, turn the Matrix to double and store in matStorage
    vector<Matrix> matStorage;
    for (int i = 5; i < input.size(); ++i)
    {
        // cout << input[i] << endl;
        if (input[i][0] == '#' && input[i][1] != 'M')
            readMatrix(input[i], matStorage);
        else if (input[i][0] == '#' && input[i][1] == 'M')
        {
            // cout << "there have M" << endl;
            stringstream ss;
            vector<vector<double>> mat;
            vector<double> vecTemp;
            vector<string> inputTemp;
            string strTemp;
            for (int a = 1; a < 5; ++a)
            {
                strTemp = input[i + a];
                ss << strTemp;
                while (getline(ss, strTemp, ' '))
                    inputTemp.push_back(strTemp);
                ss.clear();
                for (int q = 0; q < inputTemp.size(); ++q)
                    vecTemp.push_back(stod(inputTemp[q]));
                mat.push_back(vecTemp);
                vecTemp.clear();
                inputTemp.clear();
            }
            matStorage.push_back(mat);
            mat.clear();
        }
    }
    ///* // test--pass
    // cout << "how many matrix : " << matStorage.size() << endl;
    //*/
    /*
    for(int i = 0;i < matStorage.size(); ++ i) 
    {
        cout << i << endl;
        matStorage[i].print();
    }
    */
    // for(Matrix output : matStorage) output.print();
    // (matStorage[matStorage.size() - 1] * matStorage[matStorage.size() - 2]).print();
    ///*
    // calculate T and 
    Matrix mat(matStorage[0]);
    for (int q = 1; q < matStorage.size(); ++ q)
    {
        // mat.print();
        mat = matStorage[q] * mat;
    }
    // cout << "-----------------------------------" << endl;
    // cout << "T : " << endl;
    // mat.print();
    // cout << "det : ";
    // get det of mat
    double detMat = det(mat);
    // cout << detMat << endl;

    // to compare old and new volume
    vector<vector<double>> vecAns;
    vector<double> ans;
    for (int i = 0; i < vecStorage.size(); ++i)
    {
        ans = mat * vecStorage[i];
        vecAns.push_back(ans);
    }
    vector<vector<double>> oldPoint = pointToVec(vecStorage);
    vector<vector<double>> newPoint = pointToVec(vecAns);
    double r = fabs((det(newPoint) / 6) / (det(oldPoint) / 6));

   //  cout << r << endl;
    // fs.open(argv[2], ios :: out);
    // write file
    fs.open("/home/weizijun/linearAlgebra-openCV/HW2/case3/outputTest1.txt", ios :: out);
    mat.writeFile(fs);
    Matrix(vecAns).writeFile(fs);
    fs << fixed << setprecision(2) << r << ' ' << detMat << endl;
    if(r == 0 && detMat == 0) fs << "zeros" << endl;
    else if(round(r * 100) == round(detMat * 100)) fs << "r==det(T)" << endl;
    else if(round(-1 * r * 100) == round(detMat * 100)) fs << "r==-det(T)" << endl;
    else fs << "others" << endl;
    if(detMat != 0)
    {
        vector<vector<double>> inMat = inverseMat(mat.getMat());
        vector<double> inVec = Matrix(inMat) * reverseVec;
        for(int i = 0;i < inVec.size(); ++ i)
        {
            if(i < inVec.size() - 1) fs << fixed << setprecision(2) << inVec[i] << ' ';
            else fs << fixed << setprecision(2) << inVec[i];
        }
    }
    else fs << "NaN";
    fs << endl;
    fs.close();


    // cout <<det(point) << endl;
    //*/
}


