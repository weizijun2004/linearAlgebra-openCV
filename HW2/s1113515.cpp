#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<iomanip>
using namespace std;

/*double[4][4] == double** */
void getTmpToAry(fstream&, string&, double[]);
void MatrixMultiply(double[4][4], double[4][4]);
void getParameter(double&, string&, stringstream&);
void Translation(double[4][4], stringstream&);
void OrthographicProjection(double[4][4], string);
void Scaling(double[4][4], stringstream&);
void Rotation(double[4][4], string, stringstream&);
void Shearing(double[4][4], string, stringstream&);
void Customize(double[4][4], fstream&);

int main(int argc, char* argv[])
{
	// string inputPath1 = argv[1], output1 = argv[2], inputPath2 = argv[3], output2 = argv[4];
	fstream Input("/home/zj/大二上-作業文件/code/linearAlgebra/HW2/case1/input1.txt", ios::in);

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
	/*
	cout << "vec:";
	for (int i = 0; i < 4; i++)
		cout << v4[i] << " ";
	cout << "\n";
	return 0;
	//*/

	double u[4];
	getTmpToAry(Input, tmp, u);

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
		cout << ct << "\n";
		//find type
		getline(Input, tmp);
		ss << tmp;
		ss >> tmp;
		cout << "tmp[1]:" << tmp[1] << "\n";
		// system("pause");
		
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

	cout << "T:\n";
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << fixed << setprecision(2) << T[i][j] << " ";
		}
		cout << endl;
	}

}

void testMatrix(double T[4][4])
{
	cout << "T:\n";
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << T[i][j] << " ";
		}
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
    cout << "output : " << endl;
    for (int i = 0; i < 4; i++)
    {
		for (int j = 0; j < 4; j++)
			cout << T[i][j] << ' ';
        cout << endl;
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

	tx = -tx;
	ty = -ty;
	tz = -tz;

	double tmpMatrix[4][4] = {
		1, 0, 0, tx,
		0, 1, 0, ty,
		0, 0, 1, tz,
		0, 0, 0,  1
	};
	MatrixMultiply(tmpMatrix, T);

	//testMatrix(T);
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

	//testMatrix(T);
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
	MatrixMultiply(tmpMatrix, T);

	//testMatrix(tmpMatrix);
}

void Rotation(double T[4][4], string type, stringstream& ss) // 旋轉R
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
	*/

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
			 COS,   0,  SIN, cx,
			   0,   1,    0, cy,
			-SIN,   0,  COS, cz,
			   0,   0,    0,  1
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
		// testMatrix(T);
		MatrixMultiply(zMatrix, T);
	}

	//testMatrix(T);
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