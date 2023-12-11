#define _USE_MATH_DEFINES // for M_PI
#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<cmath>
#include<math.h>
#include<vector>
using namespace std;

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

void Translation(double T[4][4], stringstream& ss) // ����T
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
		1, 0, 0, tx,
		0, 1, 0, ty,
		0, 0, 1, tz,
		0, 0, 0,   1
	};
	MatrixMultiply(tmpMatrix, T);
}

void OrthographicProjection(double T[4][4], string type) //��vP
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

void Scaling(double T[4][4], stringstream& ss) // �Y��S
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
	MatrixMultiply(tmpMatrix, T);
}

void Rotation(double T[4][4], string type, stringstream& ss) // ����R
{
	string Parameter;
	stringstream tmpss;

	//find type's parameter
	getline(ss, Parameter);
	double cx, cy, cz, angle;
	getParameter(cx, Parameter, tmpss);
	getParameter(cy, Parameter, tmpss);
	getParameter(cz, Parameter, tmpss);
	getParameter(angle, Parameter, tmpss);
	double COS = cos(angle * M_PI / 180), SIN = sin(angle * M_PI / 180);
	/*
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
		MatrixMultiply(xRotationMatrix, T);
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
		MatrixMultiply(yRotationMatrix, T);
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
		MatrixMultiply(zRotationMatrix, T);
	}
}

void Shearing(double T[4][4], string type, stringstream& ss) // ����H
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
		MatrixMultiply(xMatrix, T);
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
		MatrixMultiply(yMatrix, T);
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
		MatrixMultiply(zMatrix, T);
	}
}

void Customize(double T[4][4], fstream& Input) // �ۭqM
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

void caculateT(double T[4][4], fstream& Input, int n)
{
	stringstream ss;
	string type, tmp;
	int ct = 0;

	while (ct < n)
	{
		//find type
		getline(Input, tmp);
		ss << tmp;
		ss >> tmp;
		/*
		cout << ct << "\n";
		cout << "tmp[1]:" << tmp[1] << "\n";
		system("pause");
		*/

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
}

void VecMatrixProduct(double v[], double T[4][4], double u[])
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			//cout << "u : " << u[i] << endl;
			u[i] += T[i][j] * v[j];
			//cout << T[i][j] << "," << v[j] << "," << u[i] << endl;
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
			for (int i = 0; i < n - 1; i++)
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
			T[i][j] /= division;
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
					inverse_T[k][j] -= factor * inverse_T[i][j];
				}
			}
		}
	}
}

void findnearby8pt(double pt[], double C[8][4])
{
	// C000
	// C001
	// C010
	// C011
	// C100
	// C101
	// C110
	// C111
	int dx[] = { -1, -1, -1, -1, 1, 1, 1, 1 };
	int dy[] = { -1, -1, 1, 1, -1, -1, 1, 1 };
	int dz[] = { -1, 1, -1, 1, -1, 1, -1, 1 };

	int ax[] = { 0, 0, 0, 0, 1, 1, 1, 1 };
	int ay[] = { 0, 0, 1, 1, 0, 0, 1, 1 };
	int az[] = { 0, 1, 0, 1, 0, 1, 0, 1 };


	for (int i = 0; i < 8; ++i)
	{
		int x, y, z;
		if (floor(pt[0]) == ceil(pt[0]))
			x = (ax[i] == 0 ? floor(pt[0] + ax[i]) : ceil(pt[0] + ax[i]));
		else
			x = (dx[i] > 0 ? floor(pt[0] + dx[i]) : ceil(pt[0] + dx[i]));
		
		if (floor(pt[1]) == ceil(pt[1]))
			y = (ay[i] == 0 ? floor(pt[1] + ay[i]) : ceil(pt[1] + ay[i]));
		else
			y = (dy[i] > 0 ? floor(pt[1] + dy[i]) : ceil(pt[1] + dy[i]));

		if (floor(pt[2]) == ceil(pt[2]))
			z = (az[i] == 0 ? floor(pt[2] + az[i]) : ceil(pt[2] + az[i]));
		else 
			z = (dz[i] > 0 ? floor(pt[2] + dz[i]) : ceil(pt[2] + dz[i]));

		/*
		int x = (dx[i] > 0 ? floor(pt[0] + dx[i]) : ceil(pt[0] + dx[i])); //static_cast<int>(dx[i] > 0 ? floor(pt[0] + dx[i]) : ceil(pt[0] + dx[i]))
		int y = (dy[i] > 0 ? floor(pt[1] + dy[i]) : ceil(pt[1] + dy[i])); //^~~~~~~~~~~~~~~~ C++�j���૬ 
		int z = (dz[i] > 0 ? floor(pt[2] + dz[i]) : ceil(pt[2] + dz[i]));
		*/
		/*
		nx = max(nx, 0);
		ny = max(ny, 0);
		nz = max(nz, 0);
		*/

		C[i][0] = x;
		C[i][1] = y;
		C[i][2] = z;
		C[i][3] = 1;
		//cout << x << "," << y << "," << z << endl;
		//system("pause");
	}
}

void caculateOriginal(int l, int w, int h, double C[8][4], vector<vector<vector<double>>>& Original, vector<double>& tmp)
{
	// C000
	// C001
	// C010
	// C011
	// C100
	// C101
	// C110
	// C111
	//cout << Original[0].size() << "," << Original[0][0].size() << "," << Original.size() << endl;
	for (int i = 0; i < 8; i++)
	{
		//cout << "C:" << C[i][0] << "," << C[i][1] << "," << C[i][2] << ":" << endl;
		//system("pause");
		if (C[i][0] < 0 || C[i][0] >= l || C[i][1] < 0 || C[i][1] >= w || C[i][2] < 0 || C[i][2] >= h)
			tmp.push_back(0);
		else
			tmp.push_back(Original[C[i][2]][C[i][1]][C[i][0]]); //Original(z,y,x)
	}
}

double TrilinearInterpolation(double Anspt[], double C[8][4], vector<double>& Originaltmp)
{
	// 
	//   0    1    2    3    4    5    6    7
	// C000,C001,C010,C011,C100,C101,C110,C111
	// ^~~C00~~^,^~~C01~~^,^~~C10~~^,^~~C11~~^
	// ^~~~~~~~~C0~~~~~~~^ ^~~~~~~~~C1~~~~~~~^
	// ^~~~~~~~~~~~~~~~~~~C~~~~~~~~~~~~~~~~~~^
	double C00, C01, C10, C11;	
	double C0, C1;
	double Ansvalue;

	double x = Anspt[0], y = Anspt[1], z = Anspt[2];

	// to z-axis
	/*
	C00 = Originaltmp[0] * (C[1][2] - z) + Originaltmp[1] * (z - C[0][2]);
	C01 = Originaltmp[2] * (C[3][2] - z) + Originaltmp[3] * (z - C[2][2]);
	C10 = Originaltmp[4] * (C[5][2] - z) + Originaltmp[5] * (z - C[4][2]);
	C11 = Originaltmp[6] * (C[7][2] - z) + Originaltmp[7] * (z - C[6][2]);
	*/
	C00 = Originaltmp[0] * (z - C[0][2]) + Originaltmp[1] * (C[1][2] - z);
	C01 = Originaltmp[2] * (z - C[2][2]) + Originaltmp[3] * (C[3][2] - z);
	C10 = Originaltmp[4] * (z - C[4][2]) + Originaltmp[5] * (C[5][2] - z);
	C11 = Originaltmp[6] * (z - C[6][2]) + Originaltmp[7] * (C[7][2] - z);

	// to y-axis
	C0 = C00 * (C[2][1] - y) + C01 * (y - C[0][1]);
	C1 = C10 * (C[6][1] - y) + C11 * (y - C[4][1]);

	// to x-axis
	Ansvalue = C0 * (C[4][0] - x) + C1 * (x - C[0][0]);

	return Ansvalue;
}

int main(int argc, char* argv[])
{
	// string inputPath1 = argv[1], output1 = argv[2], inputPath2 = argv[3], output2 = argv[4];
	fstream Input1("/home/zj/大二上-作業文件/code/linearAlgebra/HW2/case1/inputTest.txt", ios::in);

	double v1[4], v2[4], v3[4], v4[4];
	string tmp;
	int ct = 1;
	while (ct <= 4)
	{
		if (ct == 1)
			getTmpToAry(Input1, tmp, v1);
		else if (ct == 2)
			getTmpToAry(Input1, tmp, v2);
		else if (ct == 3)
			getTmpToAry(Input1, tmp, v3);
		else
			getTmpToAry(Input1, tmp, v4);
		ct++;
	}
	double u[4];
	getTmpToAry(Input1, tmp, u);

	int n;
	getline(Input1, tmp);
	n = stoi(tmp);

	double T[4][4] = {
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1
	};
	//////////////////////
	caculateT(T, Input1, n);

	double u1[4] = {}, u2[4] = {}, u3[4] = {}, u4[4] = {};
	VecMatrixProduct(v1, T, u1);
	VecMatrixProduct(v2, T, u2);
	VecMatrixProduct(v3, T, u3);
	VecMatrixProduct(v4, T, u4);

	fstream fsout("/home/zj/大二上-作業文件/code/linearAlgebra/HW2/case1/outputTest1.txt", ios::out);
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

	/*
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
	{
		fsout << fixed << setprecision(2) << u1[i];
		if (i < 3)
			fsout << " ";
	}
	fsout << endl;
	for (int i = 0; i < 4; i++)
	{
		fsout << fixed << setprecision(2) << u2[i];
		if (i < 3)
			fsout << " ";
	}
	fsout << endl;
	for (int i = 0; i < 4; i++)
	{
		fsout << fixed << setprecision(2) << u3[i];
		if (i < 3)
			fsout << " ";
	}
	fsout << endl;
	for (int i = 0; i < 4; i++)
	{
		fsout << fixed << setprecision(2) << u4[i];
		if (i < 3)
			fsout << " ";
	}
	fsout << endl;

	/*
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

	double det;

	double** tmpT = new double* [4];
	for (int i = 0; i < 4; i++)
		tmpT[i] = new double[4];

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			tmpT[i][j] = T[i][j];

	det = round(detT(tmpT, 4) * 100) / 100;
	fsout << det << endl;

	/*
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
		fsout << "zeros\n";
	else
		fsout << "others\n";

	double inverse_T[4][4] = {}, v[4] = {};

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			inverse_T[i][j] = (i == j) ? 1.0 : 0.0;

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
	fsout << "\n";
	fsout.close();

	/*
	// test
	cout << "\n\n" << "inverse:\n";
	if (det == 0)
		cout << "NaN\n";
	else {
		inverseMatrix(T, inverse_T);
		VecMatrixProduct(u, inverse_T, v);
		for (int i = 0; i < 4; i++)
			cout << fixed << setprecision(2) << v[i] << " ";
	}
	//*/


	// Part2
	fstream Input2("/home/zj/大二上-作業文件/code/linearAlgebra/HW2/case1/input2.txt", ios::in);
	stringstream Matrixs, ss;
	string tmpstr;

	int l= 0, w = 0, h = 0;
	getline(Input2, tmp);
	ss << tmp;
	for (int i = 0; i < 3; i++)
	{
		ss >> tmpstr;
		if (i == 0)
			l = stoi(tmpstr);
		else if (i == 1)
			w = stoi(tmpstr);
		else if (i == 2)
			h = stoi(tmpstr);
		tmpstr.clear();
	}

	
	vector<vector<double>> M;
	vector<string> S;
	
	tmp.clear();
	tmpstr.clear();
	for (int i = 0;  i < l * h;  i++)
	{
		getline(Input2, tmp);
		S.push_back(tmp);
	}
	
	tmp.clear();
	ss.clear();
	ss.str("");
	ss.str("");

	vector<vector<vector<double>>> Originalimg;
	vector<vector<vector<double>>> newimg;

	for (int i = 0; i < S.size(); i++)
	{
		vector<double> tmpVec;
		stringstream tmpss;
		int imgIndex = 0;
		tmpstr = S[i];
		ss << tmpstr;
		for (int j = 0; j < w; j++)
		{
			ss >> tmp;
			tmpVec.push_back(stod(tmp));
			tmp.clear();
		}
		M.push_back(tmpVec);
		tmpstr.clear();
		ss.str("");
		ss.clear();
	
		if (i % 2 == 1 && i != 0)
		{
			/*
			cout << "test:\n";
			for (int j = 0; j < M.size(); j++)
			{
				for (int k = 0; k < M[0].size(); k++)
					cout << M[j][k] << " ";
				cout << endl;
			}
			//cout << M.size() << "," << M[0].size();
			system("pause");
			//*/
			vector<vector<double>> ndAns;
			for (int j = 0; j < M[0].size(); j++)
			{
				vector<double>stAns;
				for (int k = 0; k < M.size(); k++)
				{
					stAns.push_back(M[k][j]);
				}
				/*
				for (int p = 0; p < M[0].size(); p++)
					cout << stAns[p] << " ";
				system("pause");
				*/
				ndAns.push_back(stAns);
			}
			Originalimg.push_back(ndAns);
			imgIndex++;
			M.clear();
		}
	}

	/*
	cout << "original[0][0].size() = y: " << Originalimg[0][0].size() << "\n";
	cout << "original[0].size() = x: " << Originalimg[0].size() << "\n";
	cout << "original.size() = z: " << Originalimg.size() << "\n";
	///*
	for (int i = 0; i < Originalimg.size(); i++)
	{
		for (int j = 0; j < Originalimg[0].size(); j++)
		{
			for (int k = 0; k < Originalimg[0][0].size(); k++)
			{
				cout << Originalimg[i][j][k]<< " ";
			}
			cout << endl;
		}
	}
	return 0;
	//*/

	tmp.clear();
	getline(Input2, tmp);
	n = stoi(tmp);
	
	double T2[4][4] = {
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1
	};
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			inverse_T[i][j] = (i == j) ? 1.0 : 0.0;

	caculateT(T2, Input2, n);
	inverseMatrix(T2, inverse_T);

	double pt[4] = {}; // [x,y,z,1]
	double AnsValue;

	for (int i = 0; i < h; i++) // h
	{
		vector<vector<double>> tmpVec_nd;
		for (int j = 0; j < l; j++) // x
		{
			vector<double> tmpVec_st;
			for (int k = 0; k < w; k++) // y
			{
				pt[0] = j;
				pt[1] = k;
				pt[2] = i;
				pt[3] = 1;


				double Anspt[4] = {};
				VecMatrixProduct(pt, inverse_T, Anspt);
				/*
				for (int l = 0; l < 4; l++)
					cout << fixed << setprecision(2) << Anspt[l] << " ";
				cout << endl;
				system("pause");
				//*/

				double x = Anspt[0], y = Anspt[1], z = Anspt[2];
				//double x = max(Anspt[0], 0.0);
				//double y = max(Anspt[1], 0.0);
				//double z = max(Anspt[2], 0.0);

				//cout << fixed << setprecision(2) << x << "," << fixed << setprecision(2) << y << "," << fixed << setprecision(2) << z << endl;

				if (x == 0 || y == 0 || z == 0 || x >= l || y >= w || z >= h)
				{
					AnsValue = 0;
					/*
					cout << "AnsValue__0: " << AnsValue << endl;
					system("pause");
					//*/
				}
				else
				{
					vector<double> tmpVec;
					double C[8][4] = {};
					findnearby8pt(Anspt, C);
					caculateOriginal(l, w, h, C, Originalimg, tmpVec);
					AnsValue = TrilinearInterpolation(Anspt, C, tmpVec);
					/*
					cout << "AnsValue__: " << AnsValue << endl;
					system("pause");
					//*/
				}
				tmpVec_st.push_back(AnsValue);
			}
			tmpVec_nd.push_back(tmpVec_st);
			//cout << endl;
		}
		newimg.push_back(tmpVec_nd);
	}

	///*
	for (int i = 0; i < newimg.size(); i++)
	{
		for (int j = 0; j < newimg[0].size(); j++)
		{
			for (int k = 0; k < newimg[0][0].size(); k++)
			{
				cout << newimg[i][j][k] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
	//*/

}

/*
void testMatrix(double T[4][4]);
void getTmpToAry(fstream&, string&, double[]);
void caculateT(double[4][4], fstream&, int);
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
Mat MatrixToMat(vector<vector<double>>&);
void findnearby8pt(double[], double[8][4]);
void TrilinearInterpolation(double[], double[8][4], vector<vector<vector<double>>>, double);
*/