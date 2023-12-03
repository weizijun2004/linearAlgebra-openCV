#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<cmath>
#include<opencv2/opencv.hpp>
using namespace std;
using namespace cv;

/*double[4][4] == double** */
void testMatrix(double T[4][4]); // test
void getTmpToAry(fstream&, string&, double[]);
void MatrixMultiply(double[4][4], double[4][4]);
void getParameter(double&, string&, stringstream&);
void Translation(double[4][4], stringstream&);
void OrthographicProjection(double[4][4], string);
void Scaling(double[4][4], stringstream&);
void Rotation(double[4][4], string, stringstream&);
void Shearing(double[4][4], string, stringstream&);
void Customize(double[4][4], fstream&);
void VecMatrixProduct(double[], double[4][4], double[]);
void caculateVector(double[], double[], double[]);
double caculateArea(double[], double[], double[]);
const int Tsize = 4;
double detT(double**,int);
void inverseMatrix(double[4][4], double[4][4]);

int main(int argc, char* argv[])
{
	string inputPath1 = argv[1], output1 = argv[2], inputPath2 = argv[3], output2 = argv[4];
	fstream Input(inputPath1, ios::in);

	double v1[4], v2[4], v3[4], v4[4];
	string tmp;
	int ct = 1;
	while (ct <= 4)
	{
		if (ct == 1)
			getTmpToAry(Input, tmp, v1);
		else if (ct == 2)
			getTmpToAry(Input, tmp, v2);
		else if (ct == 3)
			getTmpToAry(Input, tmp, v3);
		else
			getTmpToAry(Input, tmp, v4);
		ct++;
	}
	double u[4];
	getTmpToAry(Input, tmp, u);

	/*
	cout << "v:";
	for (int i = 0; i < 4; i++)
		cout << u[i] << " ";
	cout << "\n";
	system("pause");
	//*/

	int n;
	getline(Input, tmp);
	n = stoi(tmp);

	stringstream ss;
	string type, parameter;

	double T[4][4] = {
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1
	};
	ct = 0;
	tmp.clear();
	while (ct < n)
	{
		//find type
		getline(Input, tmp);
		ss << tmp;
		ss >> tmp;
		///*
		cout << ct << "\n";
		cout << "tmp[1]:" << tmp[1] << "\n";
		system("pause");
		//*/

		if (tmp[1] == 'T')
			Translation(T, ss);
		else if (tmp[1] == 'P') // Pxy, Pyz, Pxz
		{
			type += tmp[2];
			type += tmp[3];
			OrthographicProjection(T, type);
		}
		else if (tmp[1] == 'S')
			Scaling(T, ss);
		else if (tmp[1] == 'R') // Rx, Ry, Rz
		{
			type += tmp[2];
			Rotation(T, type, ss);
		}
		else if (tmp[1] == 'H') // Hx, Hy, Hz
		{
			type += tmp[2];
			Shearing(T, type, ss);
		}
		else if (tmp[1] == 'M')
			Customize(T, Input);
		type = "";
		ss.str("");
		ss.clear();
		ct++;
	}


	double u1[4] = {}, u2[4] = {}, u3[4] = {}, u4[4] = {};
	VecMatrixProduct(v1, T, u1);
	VecMatrixProduct(v2, T, u2);
	VecMatrixProduct(v3, T, u3);
	VecMatrixProduct(v4, T, u4);

	fstream fsout(output1, ios::out);
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			fsout << fixed << setprecision(2) << T[i][j];
			if (j < 3)
				fsout << " ";
		}
		fsout << endl;
	}

	///*
	cout << "T:\n";
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			cout << fixed << setprecision(2) << T[i][j] << " ";
		cout << endl;
	}
	system("pause");
	//*/
	for (int i = 0; i < 4; i++)
		fsout << fixed << setprecision(2) << u1[i] << " ";
	fsout << endl;
	for (int i = 0; i < 4; i++)
		fsout << fixed << setprecision(2) << u2[i] << " ";
	fsout << endl;
	for (int i = 0; i < 4; i++)
		fsout << fixed << setprecision(2) << u3[i] << " ";
	fsout << endl;
	for (int i = 0; i < 4; i++)
		fsout << fixed << setprecision(2) << u4[i] << " ";
	fsout << endl;

	///*
	//test
	cout << "u1:";
	for (int i = 0; i < 4; i++)
		cout << fixed << setprecision(2) << u1[i] << " ";
	cout << endl;
	cout << "u2:";
	for (int i = 0; i < 4; i++)
		cout << fixed << setprecision(2) << u2[i] << " ";
	cout << endl;
	cout << "u3:";
	for (int i = 0; i < 4; i++)
		cout << fixed << setprecision(2) << u3[i] << " ";
	cout << endl;
	cout << "u4:";
	for (int i = 0; i < 4; i++)
		cout << fixed << setprecision(2) << u4[i] << " ";
	cout << endl;
	system("pause");
	//*/

	
	double Vec1[3] = {}, Vec2[3] = {}, Vec3[3] = {};
	double vArea, uArea, r;

	caculateVector(v1, v2, Vec1);
	caculateVector(v1, v3, Vec2);
	caculateVector(v1, v4, Vec3);
	vArea = caculateArea(Vec1, Vec2, Vec3);

	caculateVector(u1, u2, Vec1);
	caculateVector(u1, u3, Vec2);
	caculateVector(u1, u4, Vec3);
	uArea = caculateArea(Vec1, Vec2, Vec3);

	r = round(fabs(uArea / vArea) * 100) / 100;
	fsout << r << " ";
	cout << r << " ";
	// cout << "r:" << r << "	";

	double det;

	double** tmpT = new double* [4];
	for (int i = 0; i < 4; i++)
		tmpT[i] = new double[4];

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			tmpT[i][j] = T[i][j];

	det = round(detT(tmpT, 4) * 100) / 100;
	fsout << det << endl;
	cout << det << endl;
	// cout << "det=" << det << "\n";

	///*
	//test
	if (r == det && r != 0 && det != 0)
		cout << "r==det(T)\n";
	else if (r == -det && r != 0 && det != 0)
		cout << "r==-det(T)\n";
	else if (r == 0 && det == 0)
		cout << "zero\n";
	else
		cout << "others\n";
	//*/
	if (r == det && r != 0 && det != 0)
		fsout << "r==det(T)\n";
	else if (r == -det && r != 0 && det != 0)
		fsout << "r==-det(T)\n";
	else if (r == 0 && det == 0)
		fsout << "zero\n";
	else
		fsout << "others\n";

	double inverse_T[4][4], v[4] = {};

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			inverse_T[i][j] = (i == j) ? 1 : 0;

	if (det == 0)
		fsout << "NaN";
	else {
		inverseMatrix(T, inverse_T);
		VecMatrixProduct(u, inverse_T, v);
		for (int i = 0; i < 4; i++)
		{
			fsout << fixed << setprecision(2) << v[i];
			if (i < 3)
				fsout << " ";
		}
	}
	///*
	// test
	if (det == 0)
		cout << "NaN\n";
	else {
		inverseMatrix(T, inverse_T);
		VecMatrixProduct(u, inverse_T, v);
		for (int i = 0; i < 4; i++)
			cout << fixed << setprecision(2) << v[i] << " ";
	}
	//*/

	fsout.close();
}

void testMatrix(double T[4][4])
{
	cout << "T:\n";
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
			cout << T[i][j] << " ";
		cout << endl;
	}
	system("pause");
}

void getTmpToAry(fstream& Input, string& tmp, double array[])
{
	stringstream ss;
	string tmpInteger;
	int index = 0;
	getline(Input, tmp);
	for (int i = 0; i < tmp.size(); i++)
	{
		if (tmp[i] != ' ')
			ss << tmp[i];
		else
		{
			ss >> tmpInteger;
			array[index] = stod(tmpInteger);
			index++;
			tmpInteger.clear();
			ss.str("");
			ss.clear();
		}

		if (i == tmp.size() - 1)
		{
			ss >> tmpInteger;
			array[index] = stod(tmpInteger);
			tmpInteger.clear();
			ss.str("");
			ss.clear();
		}
	}
	array[3] = 1;
}

void MatrixMultiply(double leftMatrix[4][4], double T[4][4])
{
	double tmpMatrix[4][4] = {};

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
				tmpMatrix[i][j] += leftMatrix[i][k] * T[k][j];
		}
	}

	//testMatrix(tmpMatrix);
	//testMatrix(leftMatrix);
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			T[i][j] = tmpMatrix[i][j];
}

void getParameter(double& value, string& Parameter, stringstream& tmpss)
{
	string tmpStr;

	tmpss << Parameter;
	tmpss >> tmpStr;
	value = stod(tmpStr);
}

void Translation(double T[4][4], stringstream& ss) // 平移T
{
	string Parameter;
	stringstream tmpss;

	//find type's parameter
	getline(ss, Parameter);
	double tx, ty, tz;
	getParameter(tx, Parameter, tmpss);
	getParameter(ty, Parameter, tmpss);
	getParameter(tz, Parameter, tmpss);

	double tmpMatrix[4][4] = {
		1, 0, 0, -tx,
		0, 1, 0, -ty,
		0, 0, 1, -tz,
		0, 0, 0,   1
	};
	MatrixMultiply(tmpMatrix, T);

	testMatrix(T);
}

void OrthographicProjection(double T[4][4], string type) //投影P
{
	string Parameter;

	//find type's parameter
	if (type == "xy")
	{
		double xyMatrix[4][4] = {
			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 0, 0,
			0, 0, 0, 1 
		};
		MatrixMultiply(xyMatrix, T);
	}
	else if (type == "yz")
	{

		double yzMatrix[4][4] = {
			0, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1 
		};
		MatrixMultiply(yzMatrix, T);
	}
	else if (type == "xz")
	{
		double xzMatrix[4][4] = {
			1, 0, 0, 0,
			0, 0, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1 
		};
		MatrixMultiply(xzMatrix, T);
	}
}

void Scaling(double T[4][4], stringstream& ss) // 縮放S
{
	string Parameter;
	stringstream tmpss;

	//find type's parameter
	getline(ss, Parameter);
	double cx, cy, cz, kx, ky, kz;
	getParameter(cx, Parameter, tmpss);
	getParameter(cy, Parameter, tmpss);
	getParameter(cz, Parameter, tmpss);
	getParameter(kx, Parameter, tmpss);
	getParameter(ky, Parameter, tmpss);
	getParameter(kz, Parameter, tmpss);

	double translationMatrix[4][4] = {
		1, 0, 0, -cx,
		0, 1, 0, -cy,
		0, 0, 1, -cz,
		0, 0, 0,   1
	};
	double tmpMatrix[4][4] = {
		kx,  0,  0, cx,
		 0, ky,  0, cy,
		 0,  0, kz, cz,
		 0,  0,  0,  1 
	};
	MatrixMultiply(translationMatrix, T);
	testMatrix(T);
	MatrixMultiply(tmpMatrix, T);
	testMatrix(T);
}

void Rotation(double T[4][4], string type, stringstream& ss) // 旋轉R
{
	string Parameter;
	stringstream tmpss;
	
	double PI = 3.14159265358979323846;
	//find type's parameter
	getline(ss, Parameter);
	double cx, cy, cz, angle;	
	getParameter(cx, Parameter, tmpss);
	getParameter(cy, Parameter, tmpss);
	getParameter(cz, Parameter, tmpss);
	getParameter(angle, Parameter, tmpss);
	double COS = cos(angle * PI / 180), SIN = sin(angle * PI / 180);
	///*
	cout << "COS:" << COS << "\n";
	cout << "SIN:" << SIN << "\n";
	cout << "-COS:" << -COS << "\n";
	cout << "-SIN:" << -SIN << "\n";
	//*/

	if (type == "x")
	{
		double translationMatrix[4][4] = {
			1, 0, 0, -cx,
			0, 1, 0, -cy,
			0, 0, 1, -cz,
			0, 0, 0,   1
		};
		double xRotationMatrix[4][4] = {
			1,   0,    0, cx,
			0, COS, -SIN, cy,
			0, SIN,  COS, cz,
			0,   0,    0,  1
		};
		MatrixMultiply(translationMatrix, T);
		testMatrix(T);
		MatrixMultiply(xRotationMatrix, T);
		testMatrix(T);
	}
	else if (type == "y")
	{
		double translationMatrix[4][4] = {
			1, 0, 0, -cx,
			0, 1, 0, -cy,
			0, 0, 1, -cz,
			0, 0, 0,   1
		};
		double yRotationMatrix[4][4] = {
			 COS,  0, SIN, cx,
			   0,  1,   0, cy,
			-SIN,  0, COS, cz,
			   0,  0,   0,  1
		};
		MatrixMultiply(translationMatrix, T);
		testMatrix(T);
		MatrixMultiply(yRotationMatrix, T);
		testMatrix(T);
	}
	else if (type == "z")
	{
		double translationMatrix[4][4] = {
			1, 0, 0, -cx,
			0, 1, 0, -cy,
			0, 0, 1, -cz,
			0, 0, 0,   1
		};
		double zRotationMatrix[4][4] = {
			COS, -SIN,   0, cx,
			SIN,  COS,   0, cy,
			  0,    0,   1, cz,
			  0,    0,   0,  1
		};

		MatrixMultiply(translationMatrix, T);
		testMatrix(T);
		MatrixMultiply(zRotationMatrix, T);
		testMatrix(T);
	}
}

void Shearing(double T[4][4], string type, stringstream& ss) // 推移H
{
	string Parameter;
	stringstream tmpss;

	//find type's parameter
	getline(ss, Parameter);
	double cx, cy, cz, s, t;
	getParameter(cx, Parameter, tmpss);
	getParameter(cy, Parameter, tmpss);
	getParameter(cz, Parameter, tmpss);
	getParameter(s, Parameter, tmpss);
	getParameter(t, Parameter, tmpss);

	if (type == "x")
	{
		double translationMatrix[4][4] = {
			1, 0, 0, -cx,
			0, 1, 0, -cy,
			0, 0, 1, -cz,
			0, 0, 0,  1
		};
		double xMatrix[4][4] = {
			1, s, t, cx,
			0, 1, 0, cy,
			0, 0, 1, cz,
			0, 0, 0,  1
		};
		MatrixMultiply(translationMatrix, T);
		testMatrix(T);
		MatrixMultiply(xMatrix, T);
		testMatrix(T);
	}
	else if (type == "y")
	{
		double translationMatrix[4][4] = {
			1, 0, 0, -cx,
			0, 1, 0, -cy,
			0, 0, 1, -cz,
			0, 0, 0,   1
		};
		double yMatrix[4][4] = {
			1, 0, 0, cx,
			s, 1, t, cy,
			0, 0, 1, cz,
			0, 0, 0,  1
		};
		MatrixMultiply(translationMatrix, T);
		testMatrix(T);
		MatrixMultiply(yMatrix, T);
		testMatrix(T);
	}
	else if (type == "z")
	{
		double translationMatrix[4][4] = {
			1, 0, 0, -cx,
			0, 1, 0, -cy,
			0, 0, 1, -cz,
			0, 0, 0,   1
		};
		double zMatrix[4][4] = {
			1, 0, 0, cx,
			0, 1, 0, cy,
			s, t, 1, cz,
			0, 0, 0,  1
		};
		MatrixMultiply(translationMatrix, T);
		testMatrix(T);
		MatrixMultiply(zMatrix, T);
		testMatrix(T);
	}
}

void Customize(double T[4][4], fstream& Input) // 自訂M
{
	double tmpMatrix[4][4];
	string Parameter;
	//find type's parameter
	for (int i = 0; i < 4; i++)
	{
		getline(Input, Parameter);
		/*
		cout << Parameter << "\n";
		system("pause");
		*/
		stringstream ss;
		string tmpInteger;
		int col = 0;

		for (int j = 0; j < Parameter.size(); j++)
		{
			if (Parameter[j] != ' ')
				ss << Parameter[j];
			else
			{
				ss >> tmpInteger;
				tmpMatrix[i][col] = stod(tmpInteger);
				col++;
				tmpInteger.clear();
				ss.str("");
				ss.clear();
			}

			if (j == Parameter.size() - 1)
			{
				ss >> tmpInteger;
				tmpMatrix[i][col] = stod(tmpInteger);
				tmpInteger.clear();
				ss.str("");
				ss.clear();
			}
		}

	}

	MatrixMultiply(tmpMatrix, T);
	//testMatrix(T);
		
}

void VecMatrixProduct(double v[], double T[4][4], double u[])
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			u[i] += T[i][j] * v[j];
		}
	}
}

void caculateVector(double P1[], double P2[], double Vec[])
{
	for (int i = 0; i < 3; i++)
		Vec[i] = P2[i] - P1[i];
}

double caculateArea(double V1[], double V2[], double V3[])
{
	/*
		v1[0], v2[0], v3[0],
		v1[1], v2[1], v3[1],
		v1[2], v2[2], v3[2],
	*/
	double Area = 0;
	Area = ((V1[0] * V2[1] * V3[2] + V2[0] * V3[1] * V1[2] + V3[0] * V1[1] * V2[2]) - (V3[0] * V2[1] * V1[2] + V3[1] * V2[2] * V1[0] + V3[2] * V2[0] * V1[1])) / 6;
	return Area;
}

double detT(double** T, int n)
{
	/*
	T[4][4] = {
		T11,T12,T13,T14
		T21,T22,T23,T24
		T31,T32,T33,T34
		T41,T42,T43,T44
	};
	//*/

	if (n == 2)
	{
		return T[0][0] * T[1][1] - T[0][1] * T[1][0];
	}
	else
	{
		double det = 0;
		for (int col = 0; col < n; col++)
		{
			double** subMatrix = new double* [n - 1];
			for (int i = 0; i < n -1; i++)
				subMatrix[i] = new double[n - 1];

			for (int i = 1; i < n; i++)
			{
				int subCol = 0;
				for (int j = 0; j < n; j++)
				{
					if (j != col) {
						subMatrix[i - 1][subCol] = T[i][j];
						subCol++;
					}
				}
			}
			

			det += (col % 2 == 0 ? 1 : -1) * T[0][col] * detT(subMatrix, n - 1);
		
			for (int i = 0; i < n - 1; i++)
				delete[] subMatrix[i];
			delete[] subMatrix;
		}
		return det;
	}
}

void inverseMatrix(double T[4][4], double inverse_T[4][4])
{
	for (int i = 0; i < 4; i++)
	{
		double division = T[i][i];
		for (int j = 0; j < 4; j++)
		{
			T[i][j]/= division;
			inverse_T[i][j] /= division;
		}

		for (int k = 0; k < 4; k++)
		{
			if (k != i)
			{
				double factor = T[k][i];
				for (int j = 0; j < 4; j++)
				{
					T[k][j] -= factor * T[i][j];
					inverse_T[k][j] -= factor * inverse_T[i][j] ;
				}
			}
		}
	}
}

