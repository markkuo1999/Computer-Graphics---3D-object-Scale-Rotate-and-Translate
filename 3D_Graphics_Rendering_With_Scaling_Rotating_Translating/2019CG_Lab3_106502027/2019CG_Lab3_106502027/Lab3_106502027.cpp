#include<gl/GLUT.H>
#include<iostream>
#include "windows.h"
#include <math.h>
#include <fstream>
#include <string>
#include<sstream>
#include<vector>
using namespace std;

float x, y;
int obcount, vicount;

struct Node {
	float dot[2000][4];
	float plane[4000][4];
	int oneplanedot[4000];
	float dn, pn;
	
	struct Node* next;
};

struct Line {
	float dot1[20000][4];
	float dot2[20000][4];
	struct Line* next;
};

void DrawPoints()
{
	//cout << x << " " << y << endl;
	//glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(255.0, 255.0, 255.0);
	glBegin(GL_POINTS);
	glVertex2f(x, y);
	glEnd();
	glFlush();

}

void MulMatrix3x3(float a[][3], float Total[][3])
{
	int i, j, k;
	float b[3][3];


	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			b[i][j] = 0;
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++) {
			for (k = 0; k < 3; k++) {
				b[i][j] += (a[i][k] * Total[k][j]);
				//cout << i << " " << j << " " << (sum[i][k] * a[k][j]) << " " << b[i][j] << endl;;
			}
		}
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			Total[i][j] = b[i][j];
	/*for (i = 0; i < 4; i++) {
		cout << "T ";
		for (j = 0; j < 4; j++) {
			cout << Total[i][j] << " ";
		}
		cout << endl;
	}*/
	//cout << endl;
}

void DrawLines(bool reverse, float x1, float y1, float x2, float y2)
{
	float jx, jy, kx, ky, m;
	int qua;
	kx = ky = qua = 0;
	float tmp;

	if (x2 - x1 == 0)
	{
		if (y2 > y1)
			ky = 1;
		else
			ky = -1;
		if (ky == 1)
			while (y1 + ky <= y2)
			{
				y = y1 + ky;
				x = x1;
				if (reverse == true) {
					tmp = x;
					x = y;
					y = tmp;
					DrawPoints();
					tmp = x;
					x = y;
					y = tmp;
				}
				else
					DrawPoints();
				y1 += ky;
			}
		else
			while (y1 + ky >= y2)
			{
				y = y1 + ky;
				x = x1;
				if (reverse == true) {
					tmp = x;
					x = y;
					y = tmp;
					DrawPoints();
					tmp = x;
					x = y;
					y = tmp;
				}
				else
					DrawPoints();
				y1 += ky;
			}
	}
	else if (y2 - y1 == 0)
	{
		if (x2 > x1)
			kx = 1;
		else
			kx = -1;
		if (kx == 1)
			while (x1 + kx <= x2)
			{
				y = y1;
				x = x1 + kx;
				if (reverse == true) {
					tmp = x;
					x = y;
					y = tmp;
					DrawPoints();
					tmp = x;
					x = y;
					y = tmp;
				}
				else
					DrawPoints();
				x1 += kx;
			}
		else
			while (x1 + kx >= x2)
			{
				y = y1;
				x = x1 + kx;
				if (reverse == true) {
					tmp = x;
					x = y;
					y = tmp;
					DrawPoints();
					tmp = x;
					x = y;
					y = tmp;
				}
				else
					DrawPoints();
				x1 += kx;
			}
	}
	m = (y2 - y1) / (x2 - x1);
	//cout << "test " << m << endl;
	if (x2 > x1 && y2 > y1)
	{
		kx = ky = 1;
		qua = 1;
	}
	else if (x2 < x1 && y2 > y1)
	{
		kx = -1;
		ky = 1;
		qua = 2;
	}
	else if (x2 > x1 && y2 < y1)
	{
		kx = 1;
		ky = -1;
		qua = 4;
	}
	else if (x2 < x1 && y2 < y1) {
		kx = ky = -1;
		qua = 3;
	}

	jx = x1;
	jy = y1;
	//cout << x1 << " " << x2 << endl;
	if (jx > x2) {
		while (jx >= x2) {
			jx += kx;
			//cout << jx << endl;
			x = jx;
			if (m * (jx - x1) + y1 < (jy + jy + ky) / 2)
			{
				if (qua == 1 || qua == 2)
					y = jy;
				else
					y = jy + ky;
			}
			else
			{
				if (qua == 1 || qua == 2)
					y = jy + ky;
				else
					y = jy;
			}
			jy = y;
			if (reverse == true) {
				tmp = x;
				x = y;
				y = tmp;
				DrawPoints();
				tmp = x;
				x = y;
				y = tmp;
			}
			else
				DrawPoints();
			//cout << "test " << endl;
		}
	}
	else if (jx < x2) {
		while (jx <= x2) {
			jx += kx;
			//cout << jx << endl;
			x = jx;
			if (m * (jx - x1) + y1 < (jy + jy + ky) / 2)
			{
				if (qua == 1 || qua == 2)
					y = jy;
				else
					y = jy + ky;
			}
			else
			{
				if (qua == 1 || qua == 2)
					y = jy + ky;
				else
					y = jy;
			}
			jy = y;
			if (reverse == true) {
				tmp = x;
				x = y;
				y = tmp;
				DrawPoints();
				tmp = x;
				x = y;
				y = tmp;
			}
			else
				DrawPoints();
			//cout << "test " << endl;
		}
	}
	//cout << "test" << endl;
}

void CheckClipping(int &a, Line* pos, int &dcount, float x1, float y1, float z1, float w1, float x2, float y2, float z2, float w2)
{
	float c1, c2, t;

	//cout << x1 << " " << y1 << " " << z1 << " " << w1 << " " << x2 << " " << y2 << " " << z2 << " " << w2 << endl;

	if ((w1 + x1 < 0 && w2 + x2 < 0) || (w1 - x1 < 0 && w2 - x2 < 0) || (w1 + y1 < 0 && w2 + y2 < 0) || (w1 - y1 < 0 && w2 - y2 < 0) || (z1 < 0 && z2 < 0) || (w1 - z1 < 0 && w2 - z2 < 0))
		a = 1;

	if (a == 0) {
		c1 = w1 + x1;
		c2 = w2 + x2;
		if (c1 < 0) {		
			t = c1 / (c1 - c2);
			x1 = x1 + t * (x2 - x1);
			y1 = y1 + t * (y2 - y1);
			z1 = z1 + t * (z2 - z1);
			w1 = w1 + t * (w2 - w1);
		}
		else if(c2 < 0) {
			t = c1 / (c1 - c2);
			x2 = x1 + t * (x2 - x1);
			y2 = y1 + t * (y2 - y1);
			z2 = z1 + t * (z2 - z1);
			w2 = w1 + t * (w2 - w1);
		}

		if ((w1 + x1 < 0 && w2 + x2 < 0) || (w1 - x1 < 0 && w2 - x2 < 0) || (w1 + y1 < 0 && w2 + y2 < 0) || (w1 - y1 < 0 && w2 - y2 < 0) || (z1 < 0 && z2 < 0) || (w1 - z1 < 0 && w2 - z2 < 0))
			a = 1;

		c1 = w1 - x1;
		c2 = w2 - x2;
		if (c1 < 0) {
			t = c1 / (c1 - c2);
			x1 = x1 + t * (x2 - x1);
			y1 = y1 + t * (y2 - y1);
			z1 = z1 + t * (z2 - z1);
			w1 = w1 + t * (w2 - w1);
		}
		else if (c2 < 0) {
			t = c1 / (c1 - c2);
			x2 = x1 + t * (x2 - x1);
			y2 = y1 + t * (y2 - y1);
			z2 = z1 + t * (z2 - z1);
			w2 = w1 + t * (w2 - w1);
		}

		if ((w1 + x1 < 0 && w2 + x2 < 0) || (w1 - x1 < 0 && w2 - x2 < 0) || (w1 + y1 < 0 && w2 + y2 < 0) || (w1 - y1 < 0 && w2 - y2 < 0) || (z1 < 0 && z2 < 0) || (w1 - z1 < 0 && w2 - z2 < 0))
			a = 1;

		c1 = w1 + y1;
		c2 = w2 + y2;
		if (c1 < 0) {
			t = c1 / (c1 - c2);
			x1 = x1 + t * (x2 - x1);
			y1 = y1 + t * (y2 - y1);
			z1 = z1 + t * (z2 - z1);
			w1 = w1 + t * (w2 - w1);
		}
		else if (c2 < 0) {
			t = c1 / (c1 - c2);
			x2 = x1 + t * (x2 - x1);
			y2 = y1 + t * (y2 - y1);
			z2 = z1 + t * (z2 - z1);
			w2 = w1 + t * (w2 - w1);
		}

		if ((w1 + x1 < 0 && w2 + x2 < 0) || (w1 - x1 < 0 && w2 - x2 < 0) || (w1 + y1 < 0 && w2 + y2 < 0) || (w1 - y1 < 0 && w2 - y2 < 0) || (z1 < 0 && z2 < 0) || (w1 - z1 < 0 && w2 - z2 < 0))
			a = 1;

		c1 = w1 - y1;
		c2 = w2 - y2;
		if (c1 < 0) {
			t = c1 / (c1 - c2);
			x1 = x1 + t * (x2 - x1);
			y1 = y1 + t * (y2 - y1);
			z1 = z1 + t * (z2 - z1);
			w1 = w1 + t * (w2 - w1);
		}
		else if (c2 < 0) {
			t = c1 / (c1 - c2);
			x2 = x1 + t * (x2 - x1);
			y2 = y1 + t * (y2 - y1);
			z2 = z1 + t * (z2 - z1);
			w2 = w1 + t * (w2 - w1);
		}

		if ((w1 + x1 < 0 && w2 + x2 < 0) || (w1 - x1 < 0 && w2 - x2 < 0) || (w1 + y1 < 0 && w2 + y2 < 0) || (w1 - y1 < 0 && w2 - y2 < 0) || (z1 < 0 && z2 < 0) || (w1 - z1 < 0 && w2 - z2 < 0))
			a = 1;

		c1 = z1;
		c2 = z2;
		if (c1 < 0) {
			t = c1 / (c1 - c2);
			x1 = x1 + t * (x2 - x1);
			y1 = y1 + t * (y2 - y1);
			z1 = z1 + t * (z2 - z1);
			w1 = w1 + t * (w2 - w1);
		}
		else if (c2 < 0) {
			t = c1 / (c1 - c2);
			x2 = x1 + t * (x2 - x1);
			y2 = y1 + t * (y2 - y1);
			z2 = z1 + t * (z2 - z1);
			w2 = w1 + t * (w2 - w1);
		}

		if ((w1 + x1 < 0 && w2 + x2 < 0) || (w1 - x1 < 0 && w2 - x2 < 0) || (w1 + y1 < 0 && w2 + y2 < 0) || (w1 - y1 < 0 && w2 - y2 < 0) || (z1 < 0 && z2 < 0) || (w1 - z1 < 0 && w2 - z2 < 0))
			a = 1;

		c1 = w1 - z1;
		c2 = w2 - z2;
		if (c1 < 0) {
			t = c1 / (c1 - c2);
			x1 = x1 + t * (x2 - x1);
			y1 = y1 + t * (y2 - y1);
			z1 = z1 + t * (z2 - z1);
			w1 = w1 + t * (w2 - w1);
		}
		else if (c2 < 0) {
			t = c1 / (c1 - c2);
			x2 = x1 + t * (x2 - x1);
			y2 = y1 + t * (y2 - y1);
			z2 = z1 + t * (z2 - z1);
			w2 = w1 + t * (w2 - w1);
		}

		if ((w1 + x1 < 0 && w2 + x2 < 0) || (w1 - x1 < 0 && w2 - x2 < 0) || (w1 + y1 < 0 && w2 + y2 < 0) || (w1 - y1 < 0 && w2 - y2 < 0) || (z1 < 0 && z2 < 0) || (w1 - z1 < 0 && w2 - z2 < 0))
			a = 1;
	}

	if (a == 0) {
		pos->dot1[dcount][0] = x1;
		//cout << "pos->dot1[dcount][0]" << pos->dot1[dcount][0] << endl;
		pos->dot1[dcount][1] = y1;
		pos->dot1[dcount][2] = z1;
		pos->dot1[dcount][3] = w1;
		pos->dot2[dcount][0] = x2;
		pos->dot2[dcount][1] = y2;
		pos->dot2[dcount][2] = z2;
		pos->dot2[dcount][3] = w2;
		pos->dot1[dcount][0] /= pos->dot1[dcount][3];/*perspective devide*/
		pos->dot1[dcount][1] /= pos->dot1[dcount][3];
		pos->dot1[dcount][2] /= pos->dot1[dcount][3];
		pos->dot1[dcount][3] /= pos->dot1[dcount][3];
		pos->dot2[dcount][0] /= pos->dot2[dcount][3];
		pos->dot2[dcount][1] /= pos->dot2[dcount][3];
		pos->dot2[dcount][2] /= pos->dot2[dcount][3];
		pos->dot2[dcount][3] /= pos->dot2[dcount][3];
		dcount++;
	}
}

void CrossProduct(float v_A[], float v_B[], float c_P[]) {
	c_P[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
	c_P[1] = v_A[2] * v_B[0] - v_A[0] * v_B[2];
	c_P[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
}

void MulMatrix(float a[][4], float Total[][4])
{
	int i, j, k;
	float b[4][4];


	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			b[i][j] = 0;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++) {
			for (k = 0; k < 4; k++) {
				b[i][j] += (a[i][k] * Total[k][j]);
				//cout << i << " " << j << " " << (sum[i][k] * a[k][j]) << " " << b[i][j] << endl;;
			}
		}
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			Total[i][j] = b[i][j];
	/*for (i = 0; i < 4; i++) {
		cout << "T ";
		for (j = 0; j < 4; j++) {
			cout << Total[i][j] << " ";
		}
		cout << endl;
	}*/
	//cout << endl;

}

void Inverse(float T[][4])
{
	float inv[16], det, invOut[16];
	int i, j, k;

	inv[0] = T[1][1] * T[2][2] * T[3][3] -
		T[1][1] * T[2][3] * T[3][2] -
		T[2][1] * T[1][2] * T[3][3] +
		T[2][1] * T[1][3] * T[3][2] +
		T[3][1] * T[1][2] * T[2][3] -
		T[3][1] * T[1][3] * T[2][2];

	inv[4] = -T[1][0] * T[2][2] * T[3][3] +
		T[1][0] * T[2][3] * T[3][2] +
		T[2][0] * T[1][2] * T[3][3] -
		T[2][0] * T[1][3] * T[3][2] -
		T[3][0] * T[1][2] * T[2][3] +
		T[3][0] * T[1][3] * T[2][2];

	inv[8] = T[1][0] * T[2][1] * T[3][3] -
		T[1][0] * T[2][3] * T[3][1] -
		T[2][0] * T[1][1] * T[3][3] +
		T[2][0] * T[1][3] * T[3][1] +
		T[3][0] * T[1][1] * T[2][3] -
		T[3][0] * T[1][3] * T[2][1];

	inv[12] = -T[1][0] * T[2][1] * T[3][2] +
		T[1][0] * T[2][2] * T[3][1] +
		T[2][0] * T[1][1] * T[3][2] -
		T[2][0] * T[1][2] * T[3][1] -
		T[3][0] * T[1][1] * T[2][2] +
		T[3][0] * T[1][2] * T[2][1];

	inv[1] = -T[0][1] * T[2][2] * T[3][3] +
		T[0][1] * T[2][3] * T[3][2] +
		T[2][1] * T[0][2] * T[3][3] -
		T[2][1] * T[0][3] * T[3][2] -
		T[3][1] * T[0][2] * T[2][3] +
		T[3][1] * T[0][3] * T[2][2];

	inv[5] = T[0][0] * T[2][2] * T[3][3] -
		T[0][0] * T[2][3] * T[3][2] -
		T[2][0] * T[0][2] * T[3][3] +
		T[2][0] * T[0][3] * T[3][2] +
		T[3][0] * T[0][2] * T[2][3] -
		T[3][0] * T[0][3] * T[2][2];

	inv[9] = -T[0][0] * T[2][1] * T[3][3] +
		T[0][0] * T[2][3] * T[3][1] +
		T[2][0] * T[0][1] * T[3][3] -
		T[2][0] * T[0][3] * T[3][1] -
		T[3][0] * T[0][1] * T[2][3] +
		T[3][0] * T[0][3] * T[2][1];

	inv[13] = T[0][0] * T[2][1] * T[3][2] -
		T[0][0] * T[2][2] * T[3][1] -
		T[2][0] * T[0][1] * T[3][2] +
		T[2][0] * T[0][2] * T[3][1] +
		T[3][0] * T[0][1] * T[2][2] -
		T[3][0] * T[0][2] * T[2][1];

	inv[2] = T[0][1] * T[1][2] * T[3][3] -
		T[0][1] * T[1][3] * T[3][2] -
		T[1][1] * T[0][2] * T[3][3] +
		T[1][1] * T[0][3] * T[3][2] +
		T[3][1] * T[0][2] * T[1][3] -
		T[3][1] * T[0][3] * T[1][2];

	inv[6] = -T[0][0] * T[1][2] * T[3][3] +
		T[0][0] * T[1][3] * T[3][2] +
		T[1][0] * T[0][2] * T[3][3] -
		T[1][0] * T[0][3] * T[3][2] -
		T[3][0] * T[0][2] * T[1][3] +
		T[3][0] * T[0][3] * T[1][2];

	inv[10] = T[0][0] * T[1][1] * T[3][3] -
		T[0][0] * T[1][3] * T[3][1] -
		T[1][0] * T[0][1] * T[3][3] +
		T[1][0] * T[0][3] * T[3][1] +
		T[3][0] * T[0][1] * T[1][3] -
		T[3][0] * T[0][3] * T[1][1];

	inv[14] = -T[0][0] * T[1][1] * T[3][2] +
		T[0][0] * T[1][2] * T[3][1] +
		T[1][0] * T[0][1] * T[3][2] -
		T[1][0] * T[0][2] * T[3][1] -
		T[3][0] * T[0][1] * T[1][2] +
		T[3][0] * T[0][2] * T[1][1];

	inv[3] = -T[0][1] * T[1][2] * T[2][3] +
		T[0][1] * T[1][3] * T[2][2] +
		T[1][1] * T[0][2] * T[2][3] -
		T[1][1] * T[0][3] * T[2][2] -
		T[2][1] * T[0][2] * T[1][3] +
		T[2][1] * T[0][3] * T[1][2];

	inv[7] = T[0][0] * T[1][2] * T[2][3] -
		T[0][0] * T[1][3] * T[2][2] -
		T[1][0] * T[0][2] * T[2][3] +
		T[1][0] * T[0][3] * T[2][2] +
		T[2][0] * T[0][2] * T[1][3] -
		T[2][0] * T[0][3] * T[1][2];

	inv[11] = -T[0][0] * T[1][1] * T[2][3] +
		T[0][0] * T[1][3] * T[2][1] +
		T[1][0] * T[0][1] * T[2][3] -
		T[1][0] * T[0][3] * T[2][1] -
		T[2][0] * T[0][1] * T[1][3] +
		T[2][0] * T[0][3] * T[1][1];

	inv[15] = T[0][0] * T[1][1] * T[2][2] -
		T[0][0] * T[1][2] * T[2][1] -
		T[1][0] * T[0][1] * T[2][2] +
		T[1][0] * T[0][2] * T[2][1] +
		T[2][0] * T[0][1] * T[1][2] -
		T[2][0] * T[0][2] * T[1][1];

	det = T[0][0] * inv[0] + T[0][1] * inv[4] + T[0][2] * inv[8] + T[0][3] * inv[12];

	

	det = 1.0 / det;

	for (i = 0; i < 16; i++) 
		invOut[i] = inv[i] * det;

	T[0][0] = invOut[0];
	T[0][1] = invOut[1];
	T[0][2] = invOut[2];
	T[0][3] = invOut[3];
	T[1][0] = invOut[4];
	T[1][1] = invOut[5];
	T[1][2] = invOut[6];
	T[1][3] = invOut[7];
	T[2][0] = invOut[8];
	T[2][1] = invOut[9];
	T[2][2] = invOut[10];
	T[2][3] = invOut[11];
	T[3][0] = invOut[12];
	T[3][1] = invOut[13];
	T[3][2] = invOut[14];
	T[3][3] = invOut[15];
	
}

void DataMulMatrix(float &a, float &b, float &c, float &d, float Total[][4])
{
	float Data[4][1];
	float tmp[4][1], delta;
	int i, j, k;

	Data[0][0] = a;
	Data[1][0] = b;
	Data[2][0] = c;
	Data[3][0] = d;
	delta = 0.00001;
	for (i = 0; i < 4; i++)
			tmp[i][0] = 0;
	for (i = 0; i < 4; i++)	
		for (k = 0; k < 4; k++) {
			if (fabs(tmp[i][0] + Total[i][k] * Data[k][0]) <= delta) {
				//cout << "yes" << endl;
				tmp[i][0] = 0;
				//cout << tmp[i][0] << endl;
			}
			else
				tmp[i][0] += (Total[i][k] * Data[k][0]);

			//cout << "Total[" << i << "][" << k << "] = " << Total[i][k] << ", " << "Data[" << k << "][0] = " << Data[k][0] << ", " << "tmp[" << i << "][0] = " << tmp[i][0] << endl;
		}
	//cout << endl;
	a = tmp[0][0];
	b = tmp[1][0];
	c = tmp[2][0];
	d = tmp[3][0];
}

void Rotate(float Rx[][4], float Ry[][4], float Rz[][4], float x, float y, float z)
{
	float pi = 3.14159265359, c, s;

	if(x == 90)
		Rx[1][1] = Rx[2][2] = 0;
	else
		Rx[1][1] = Rx[2][2] = cosf((pi * x) / 180);
	if(y == 90)
		Ry[0][0] = Ry[2][2] = 0;
	else
		Ry[0][0] = Ry[2][2] = cosf((pi * y) / 180);
	if(z == 90)
		Rz[0][0] = Rz[1][1] = 0;
	else
		Rz[0][0] = Rz[1][1] = cosf((pi * z) / 180);

	if(x == 0)
		Rx[2][1] = 0;
	else
		Rx[2][1] = sinf((pi * x) / 180);
	Rx[1][2] = -Rx[2][1];
	if(y == 0)
		Ry[0][2] = 0;
	else
		Ry[0][2] = sinf((pi * y) / 180);
	Ry[2][0] = -Ry[0][2];
	if(z == 0)
		Rz[1][0] = 0;
	else
		Rz[1][0] = sinf((pi * z) / 180);
	Rz[0][1] = -Rz[1][0];
}

void DealParameters(vector<string> res, float P[], int len)
{
	int i, j, size;
	float mul;
	bool dot = false;

	for (i = 0; i < len; i++) 
		P[i] = 0;
	size = res.size();
	for (i = 1; i < size; i++) {
		mul = 1;
		for (j = 0; j < res[i].length() ;j++) {
			if (res[i][j] == '.') {
				j++;
				while (j < res[i].length()) {
					mul /= 10;
					j++;
				}
				break;
			}
		}
		for (j = res[i].length() - 1; j >= 0; j--) {
			if (res[i][j] >= '0' && res[i][j] <= '9') {				
				P[i - 1] += (res[i][j] - 48) * mul;			
				mul *= 10;
			}
			else if (res[i][j] == '-') 	
				P[i - 1] *= -1;
		}
	}
}

void Reset(float S[][4], float T[][4], float Rx[][4], float Ry[][4], float Rz[][4], float Total[][4])
{
	int i, j;

	for (i = 0; i < 4; i++) 
		for (j = 0; j < 4; j++) {
			if (i == j)
				S[i][j] = T[i][j] = Rx[i][j] = Ry[i][j] = Rz[i][j] = Total[i][j] = 1;
			else
				S[i][j] = T[i][j] = Rx[i][j] = Ry[i][j] = Rz[i][j] = Total[i][j] = 0;
		}
}

void Initial(float W, float H)//初始化函数 
{
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);//白色背景，前3个是RGB，最后是Alpha值，用来控制透明，1.0表示完全不透明
	glMatrixMode(GL_PROJECTION);//OpenGL按照三维方式来处理图像，所以需要一个投影变换将三维图形投影到显示器的二维空间中
	gluOrtho2D(0.0, W, 0.0, H);//指定使用正投影将一个x坐标在0~200，y坐标0~150范围内的矩形坐标区域投影到显示器窗口

}

void displayFunc(char* str)
{
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	glFlush();
	fstream file, ascfile;
	char buffer[100], ascbuffer[6000];
	float S[4][4], T[4][4], Rx[4][4], Ry[4][4], Rz[4][4], Total[4][4], PM[4][4], A[4][4], B[4][4], C[4][4], Width, Height, Hav, pi, a[3][3], b[3][3], c[3][3], EM[4][4];
	float Ex, Ey, Ez, COIx, COIy, COIz, Tilt, Hither, Yon;
	int i, dcount, pcount, j;
	bool nobackface, TMmul = false, EMmul = false;
	Node* head = new Node;
	Node* curr = head;
	head->next = NULL;
	curr->next = NULL;
	Line* start = new Line;
	Line* pos = start;
	start->next = NULL;
	pos->next = NULL;
	Reset(S, T, Rx, Ry, Rz, Total);
	pi = Width = Height = Hav = Ex = Ey = Ez = COIx = COIy = COIz = Tilt =  Hither = Yon = 0;
	nobackface = false;
	for(i = 0 ; i < 2000 ; i ++)
		for (j = 0; j < 4; j++) {
			curr->dot[i][j] = 0;
		}
	for (i = 0; i < 4000; i++)
		for (j = 0; j < 4; j++) {
			curr->plane[i][j] = 0;
		}
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			if (i == j)
				EM[i][j] = PM[i][j] = 1;
			else
				EM[i][j] = PM[i][j] = 0;
		}
	}
	file.open(str, ios::in);
	if (!file) {
		cout << "file can't open" << endl;
	}
	else {
		while (!file.eof()) {
			file.getline(buffer, sizeof(buffer));
			vector<string> res;
			string result;
			stringstream input(buffer);
			while (input >> result)
				res.push_back(result);
			if (buffer[0] >= '0' && buffer[0] <= '9') {
				Width = (buffer[0] - 48) * 100;
				Height = (buffer[4] - 48) * 100;
				cout << "Width = " << Width << endl;
				cout << "Height = " << Height << endl;
				glutInitWindowSize(Width, Height);//设定窗口的大小
				glutCreateWindow("3DpictureHW2");//创建一个窗口，参数是窗口标题名
				Initial(Width, Height);
				glClear(GL_COLOR_BUFFER_BIT);
				glutDisplayFunc(DrawPoints);
				//glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
			}
			else if (buffer[0] == 'r' && buffer[1] == 'e') {/*Reset*/
				Reset(S, T, Rx, Ry, Rz, Total);
				TMmul = EMmul = false;
			}			
			else if (buffer[0] == 's') {/*Scale*/
				float P[3];
				TMmul = true;
				DealParameters(res, P, 3);
				S[0][0] = P[0];
				S[1][1] = P[1];
				S[2][2] = P[2];
				/*cout << "Scale Matrix :" << endl;
				for (i = 0; i < 4; i++) {
					for (j = 0; j < 4; j++) {
						cout << S[i][j] << " ";
					}
					cout << endl;
				}*/
				MulMatrix(S, Total);
				/*cout << "Total Matrix :" << endl;
				for (i = 0; i < 4; i++) {
					for (j = 0; j < 4; j++) {
						cout << Total[i][j] << " ";
					}
					cout << endl;
				}*/
			}		
			else if (buffer[0] == 't') {/*Translate*/
				float P[3];
				TMmul = true;
				DealParameters(res, P, 3);
				T[0][3] = P[0];
				T[1][3] = P[1];
				T[2][3] = P[2];
				/*cout << "Translate Matrix :" << endl;
				for (i = 0; i < 4; i++) {
					for (j = 0; j < 4; j++) {
						cout << T[i][j] << " ";
					}
					cout << endl;
				}*/
				MulMatrix(T, Total);
				/*cout << "Total Matrix :" << endl;
				for (i = 0; i < 4; i++) {
					for (j = 0; j < 4; j++) {
						cout << Total[i][j] << " ";
					}
					cout << endl;
				}*/
			}					
			else if (buffer[0] == 'r' && buffer[1] == 'o') {/*Rotate*/
				float P[3];
				TMmul = true;
				DealParameters(res, P, 3);
				Rotate(Rx, Ry, Rz, P[0], P[1], P[2]);
				/*cout << "Rx Matrix :" << endl;
				for (i = 0; i < 4; i++) {
					for (j = 0; j < 4; j++) {
						cout << Rx[i][j] << " ";
					}
					cout << endl;
				}
				cout << "Ry Matrix :" << endl;
				for (i = 0; i < 4; i++) {
					for (j = 0; j < 4; j++) {
						cout << Ry[i][j] << " ";
					}
					cout << endl;
				}
				cout << "Rz Matrix :" << endl;
				for (i = 0; i < 4; i++) {
					for (j = 0; j < 4; j++) {
						cout << Rz[i][j] << " ";
					}
					cout << endl;
				}*/
				MulMatrix(Ry, Total);
				MulMatrix(Rx, Total);
				MulMatrix(Rz, Total);
				/*cout << "Total Matrix :" << endl;
				for (i = 0; i < 4; i++) {
					for (j = 0; j < 4; j++) {
						cout << Total[i][j] << " ";
					}
					cout << endl;
				}*/
			}
			else if (buffer[2] == 'j') {/*Object*/
				if (buffer[7] == 'c') {
					ascfile.open("cube.asc", ios::in);
				}
				else if (buffer[7] == 'b') {
					ascfile.open("bench.asc", ios::in);
				}
				else if (buffer[7] == 't') {
					ascfile.open("teapot.asc", ios::in);
				}
				if (!ascfile) {
					cout << "asc can't open" << endl;
				}
				else {
					if (buffer[7] == 'c') {
						curr->dn = 8;
						curr->pn = 6;
					}
					else if (buffer[7] == 'b') {
						curr->dn = 264;
						curr->pn = 244;
					}
					else if (buffer[7] == 't') {
						curr->dn = 1976;
						curr->pn = 3751;
					}
					dcount = pcount = 0;
					while (!ascfile.eof()) {
						ascfile.getline(ascbuffer, sizeof(ascbuffer));
						vector<string> ascres;
						string ascresult;
						stringstream ascinput(ascbuffer);
						while (ascinput >> ascresult)
							ascres.push_back(ascresult);
						if (ascres.size() == 3) {
							ascres.insert(ascres.begin(), "A");/*單純方便之後DealParameters的使用，A隨便選的*/
							float P[3];

							DealParameters(ascres, P, 3);
							curr->dot[dcount][0] = P[0];
							curr->dot[dcount][1] = P[1];
							curr->dot[dcount][2] = P[2];
							curr->dot[dcount][3] = 1;
							//cout << curr->dot[dcount][0] << " " << curr->dot[dcount][1] << " " << curr->dot[dcount][2] << endl;
							dcount++;
						}
						else if (ascres.size() == 4) {
							float P[3];

							DealParameters(ascres, P, 3);
							curr->plane[pcount][0] = P[0];
							curr->plane[pcount][1] = P[1];
							curr->plane[pcount][2] = P[2];
							curr->oneplanedot[pcount] = 3;
							//cout << curr->plane[pcount][0] << " " << curr->plane[pcount][1] << " " << curr->plane[pcount][2] << " " << curr->plane[pcount][3] << endl;
							pcount++;
						}
						else if (ascres.size() == 5) {
							float P[4];

							DealParameters(ascres, P, 4);
							curr->plane[pcount][0] = P[0];
							curr->plane[pcount][1] = P[1];
							curr->plane[pcount][2] = P[2];
							curr->plane[pcount][3] = P[3];
							curr->oneplanedot[pcount] = 4;
							//cout << curr->plane[pcount][0] << " " << curr->plane[pcount][1] << " " << curr->plane[pcount][2] << " " << curr->plane[pcount][3] << endl;
							pcount++;
						}
					}
					/*cout << "Total: " << endl;
					for (i = 0; i < 4; i++) {
						for (j = 0; j < 4; j++) {
							cout << Total[i][j] << " ";
						}
						cout << endl;
					}*/
					//curr = head;
					//while (curr->next != NULL) {
						for (i = 0; i < curr->dn; i++) {
							DataMulMatrix(curr->dot[i][0], curr->dot[i][1], curr->dot[i][2], curr->dot[i][3], Total);
							TMmul = true;
							//Totaltmp[i][0] = curr->dot[i][0];
							//Totaltmp[i][1] = curr->dot[i][1];
							//Totaltmp[i][2] = curr->dot[i][2];
							//Totaltmp[i][3] = curr->dot[i][3];
							//cout << "dot " << curr->dot[i][0] << " " << curr->dot[i][1] << " " << curr->dot[i][2] << " " << curr->dot[i][3] << endl;
						}
						//curr = curr->next;
					//}
					
					if (TMmul == true) {
						//curr = head;
						//while (curr->next != NULL) {
							for (i = 0; i < curr->dn; i++) {
								/*if (obcount >= 2) {
									Inverse(PM);
									DataMulMatrix(curr->dot[i][0], curr->dot[i][1], curr->dot[i][2], curr->dot[i][3], PM);
								}*/
								DataMulMatrix(curr->dot[i][0], curr->dot[i][1], curr->dot[i][2], curr->dot[i][3], EM);
								DataMulMatrix(curr->dot[i][0], curr->dot[i][1], curr->dot[i][2], curr->dot[i][3], PM);
								//curr->dot[i][0] = Totaltmp[i][0];
								//curr->dot[i][1] = Totaltmp[i][1];
								//curr->dot[i][2] = Totaltmp[i][2];
								//curr->dot[i][3] = Totaltmp[i][3];
								//cout << "EM dot " << curr->dot[i][0] << " " << curr->dot[i][1] << " " << curr->dot[i][2] << " " << curr->dot[i][3] << endl;
							}
							//curr = curr->next;
						//}
						//EMmul = true;
					}
						
					/*Inverse(Rz);/*plane的rotation要另外處理，所以先把rotate matrix拿掉*/
					/*Inverse(Rx);
					Inverse(Ry);
					MulMatrix(Rz, Total);
					MulMatrix(Rx, Total);
					MulMatrix(Ry, Total);
					Inverse(Total);
					Transpose(Total);*/
					Node* node = new Node;
					node->next = NULL;
					curr->next = node;
					curr = curr->next;
					curr->next = NULL;
					Line* lnode = new Line;
					lnode->next = NULL;
					pos->next = lnode;
					pos = pos->next;
					pos->next = NULL;
					if (curr == NULL)
						cout << "NULL" << endl;
					ascfile.close();
					//cout << endl;
					//file.getline(buffer, sizeof(buffer));
					//file.getline(buffer, sizeof(buffer));
				}

			}

			else if (buffer[0] == 'o' && buffer[2] == 's') {/*Observer*/
				
				float P[10];
				float Vz[3], Vt[3], V1[3], V2[3], V3[3], TiltM[4][4], EyeLoc[4][4];
				Ex = Ey = Ez = COIx = COIy = COIz = Tilt = Hither = Yon = Hav = 0;
				obcount++;
				curr = head;
				if (obcount >= 2) {
					Inverse(PM);
					Inverse(EM);
					while (curr->next != NULL) {
						for (i = 0; i < curr->dn; i++) {						
							DataMulMatrix(curr->dot[i][0], curr->dot[i][1], curr->dot[i][2], curr->dot[i][3], PM);
							DataMulMatrix(curr->dot[i][0], curr->dot[i][1], curr->dot[i][2], curr->dot[i][3], EM);
						}
						curr = curr->next;
					}
					Inverse(PM);
					//cout << "Inverse" << endl;
				}				
				P[0] = Ex;
				P[1] = Ey;
				P[2] = Ez;
				P[3] = COIx;
				P[4] = COIy;
				P[5] = COIz;
				P[6] = Tilt;
				P[7] = Hither;
				P[8] = Yon;
				P[9] = Hav;
				DealParameters(res, P, 10);
				Ex = P[0];
				Ey = P[1];
				Ez = P[2];
				COIx = P[3];
				COIy = P[4];
				COIz = P[5];
				Tilt = P[6];
				Hither = P[7];
				Yon = P[8];
				Hav = P[9];
				Vz[0] = V3[0] = COIx - Ex;
				Vz[1] = V3[1] = COIy - Ey;
				Vz[2] = V3[2] = COIz - Ez;
				Vt[0] = Vt[2] = 0;
				Vt[1] = 1; 
				pi = 3.14159265359;
				CrossProduct(Vt, Vz, V1);
				CrossProduct(V3, V1, V2);
				float tmp1[3], tmp2[3], tmp3[3];
				for (i = 0; i < 3; i++) {
					tmp1[i] = V1[i];
					tmp2[i] = V2[i];
					tmp3[i] = V3[i];
				}
				//cout << "V1 =" << V1[0] << " " << V1[1] << " " << V1[2] << endl;
				//cout << "V2 =" << V2[0] << " " << V2[1] << " " << V2[2] << endl;
				//cout << "V3 =" << V3[0] << " " << V3[1] << " " << V3[2] << endl;
				for (i = 0; i < 3; i++) {
					V1[i] /= sqrtf(tmp1[0] * tmp1[0] + tmp1[1] * tmp1[1] + tmp1[2] * tmp1[2]);
					V2[i] /= sqrtf(tmp2[0] * tmp2[0] + tmp2[1] * tmp2[1] + tmp2[2] * tmp2[2]);
					V3[i] /= sqrtf(tmp3[0] * tmp3[0] + tmp3[1] * tmp3[1] + tmp3[2] * tmp3[2]);
				}
				//cout << "V1 =" << V1[0] << " " << V1[1] << " " << V1[2] << endl;
				//cout << "V2 =" << V2[0] << " " << V2[1] << " " << V2[2] << endl;
				//cout << "V3 =" << V3[0] << " " << V3[1] << " " << V3[2] << endl;
				float GRM[4][4], Mirror[4][4];

				for (i = 0; i < 4; i++) {
					for (j = 0; j < 4; j++) {
						if (i == 0 && j != 3) {
							GRM[i][j] = V1[j];
						}
						else if (i == 1 && j != 3) {
							GRM[i][j] = V2[j];
						}
						else if (i == 2 && j != 3) {
							GRM[i][j] = V3[j];
						}
						else
							GRM[i][j] = 0;

						if (i == j)
							Mirror[i][j] = TiltM[i][j] = EyeLoc[i][j] = EM[i][j] = 1;
						else
							Mirror[i][j] = TiltM[i][j] = EyeLoc[i][j] = EM[i][j] = 0;
					}
				}
				GRM[3][3] = 1;
				Mirror[0][0] = -1;
				TiltM[0][0] = TiltM[1][1] = cosf((pi * Tilt) / 180);
				TiltM[0][1] = sinf((pi * Tilt) / 180);
				TiltM[1][0] = -TiltM[0][1];
				EyeLoc[0][3] = -Ex;
				EyeLoc[1][3] = -Ey;
				EyeLoc[2][3] = -Ez;
				/*cout << "Mirror Matrix :" << endl;
				for (i = 0; i < 4; i++) {
					for (j = 0; j < 4; j++) {
						cout << Mirror[i][j] << " ";
					}
					cout << endl;
				}
				cout << "EyeLoc Matrix :" << endl;
				for (i = 0; i < 4; i++) {
					for (j = 0; j < 4; j++) {
						cout << EyeLoc[i][j] << " ";
					}
					cout << endl;
				}
				cout << "GRM Matrix :" << endl;
				for (i = 0; i < 4; i++) {
					for (j = 0; j < 4; j++) {
						cout << GRM[i][j] << " ";
					}
					cout << endl;
				}
				cout << "TiltM Matrix :" << endl;
				for (i = 0; i < 4; i++) {
					for (j = 0; j < 4; j++) {
						cout << TiltM[i][j] << " ";
					}
					cout << endl;
				}*/
				MulMatrix(EyeLoc, EM);
				/*for (i = 0; i < 4; i++) {
					for (j = 0; j < 4; j++) {
						cout << "EM[" << i << "][" << j << "]" << EM[i][j] << endl;
					}
				}*/
				MulMatrix(GRM, EM);
				/*for (i = 0; i < 4; i++) {
					for (j = 0; j < 4; j++) {
						cout << "EM[" << i << "][" << j << "]" << EM[i][j] << endl;
					}
				}*/
				MulMatrix(Mirror, EM);
				/*for (i = 0; i < 4; i++) {
					for (j = 0; j < 4; j++) {
						cout << "EM[" << i << "][" << j << "]" << EM[i][j] << endl;
					}
				}*/
				MulMatrix(TiltM, EM);
				/*cout << "EM Matrix :" << endl;
				for (i = 0; i < 4; i++) {
					for (j = 0; j < 4; j++) {
						cout << EM[i][j] << " ";
					}
					cout << endl;
				}*/

				if (TMmul == true) {
					curr = head;
					while (curr->next != NULL) {
						for (i = 0; i < curr->dn; i++) {
							/*if (obcount >= 2) {
								Inverse(PM);
								DataMulMatrix(curr->dot[i][0], curr->dot[i][1], curr->dot[i][2], curr->dot[i][3], PM);
							}*/
							DataMulMatrix(curr->dot[i][0], curr->dot[i][1], curr->dot[i][2], curr->dot[i][3], EM);
							//DataMulMatrix(curr->dot[i][0], curr->dot[i][1], curr->dot[i][2], curr->dot[i][3], PM);
							//curr->dot[i][0] = Totaltmp[i][0];
							//curr->dot[i][1] = Totaltmp[i][1];
							//curr->dot[i][2] = Totaltmp[i][2];
							//curr->dot[i][3] = Totaltmp[i][3];
							//cout << "EM dot " << curr->dot[i][0] << " " << curr->dot[i][1] << " " << curr->dot[i][2] << " " << curr->dot[i][3] << endl;
						}
						curr = curr->next;
					}
					EMmul = true;
				}
				
				PM[3][3] = 0;
				PM[2][2] = (Yon * tanf((pi * Hav) / 180)) / (Yon - Hither);
				PM[2][3] = (Hither * Yon * tanf((pi * Hav) / 180)) / (Hither - Yon);
				PM[3][2] = tanf((pi * Hav) / 180);
				
			}
			else if (buffer[0] == 'v') {
				float P[4], vxl, vxr, vyb, vyt, Ar;
				vicount++;
				if (vicount >= 2) {
					Inverse(PM);
					while (curr->next != NULL) {
						for (i = 0; i < curr->dn; i++) {
							DataMulMatrix(curr->dot[i][0], curr->dot[i][1], curr->dot[i][2], curr->dot[i][3], PM);
							/*curr->dot[i][0] /= curr->dot[i][3];
							curr->dot[i][1] /= curr->dot[i][3];
							curr->dot[i][2] /= curr->dot[i][3];
							curr->dot[i][3] /= curr->dot[i][3];*/
							//cout << "PM dot " << curr->dot[i][0] << " " << curr->dot[i][1] << " " << curr->dot[i][2] << " " << curr->dot[i][3] << endl;
						}
						curr = curr->next;
					}
				}
				DealParameters(res, P, 4);
				vxl = P[0];
				vxr = P[1];
				vyb = P[2];
				vyt = P[3];
				Ar = (vxr - vxl) / (vyt - vyb);
				PM[1][1] = Ar;
				PM[3][3] = 0;
				PM[2][2] = (Yon * tanf((pi * Hav) / 180)) / (Yon - Hither);
				PM[2][3] = (Hither * Yon * tanf((pi * Hav) / 180)) / (Hither - Yon);
				PM[3][2] = tanf((pi * Hav) / 180);
				/*cout << "PM Matrix :" << endl;
				for (i = 0; i < 4; i++) {
					for (j = 0; j < 4; j++) {
						cout << PM[i][j] << " ";
					}
					cout << endl;
				}*/
				for (i = 0; i < 4; i++)
					for (j = 0; j < 4; j++) {
						if (i == j)
							A[i][j] = B[i][j] = C[i][j] = 1;
						else
							A[i][j] = B[i][j] = C[i][j] = 0;
					}
				for (i = 0; i < 3; i++)
					for (j = 0; j < 3; j++) {
						if (i == j)
							a[i][j] = b[i][j] = c[i][j] = 1;
						else
							a[i][j] = b[i][j] = c[i][j] = 0;
					}
				A[0][3] = (Width/2)*(1 + vxl);
				A[1][3] = (Height/2)*(1 + vyb);
				A[2][3] = 1;
				a[0][2] = (Width / 2) * (1 + vxl);
				a[1][2] = (Height / 2) * (1 + vyb);
				/*cout << "A Matrix :" << endl;
				for (i = 0; i < 4; i++) {
					for (j = 0; j < 4; j++) {
						cout << A[i][j] << " ";
					}
					cout << endl;
				}
				cout << "a Matrix :" << endl;
				for (i = 0; i < 3; i++) {
					for (j = 0; j < 3; j++) {
						cout << a[i][j] << " ";
					}
					cout << endl;
				}*/
				B[0][0] = ((Width / 2) * (1 + vxr) - (Width / 2) * (1 + vxl)) / (1 - (-1));
				b[0][0] = ((Width / 2) * (1 + vxr) - (Width / 2) * (1 + vxl)) / (1 - (-1));
				//cout << "B[0][0] = " << B[0][0] << endl;
				B[1][1] = ((Height / 2) * (1 + vyt) - (Height / 2) * (1 + vyb)) / (1 - (-1));
				b[1][1] = ((Height / 2) * (1 + vyt) - (Height / 2) * (1 + vyb)) / (1 - (-1));
				/*cout << "B Matrix :" << endl;
				for (i = 0; i < 4; i++) {
					for (j = 0; j < 4; j++) {
						cout << B[i][j] << " ";
					}
					cout << endl;
				}
				cout << "b Matrix :" << endl;
				for (i = 0; i < 3; i++) {
					for (j = 0; j < 3; j++) {
						cout << b[i][j] << " ";
					}
					cout << endl;
				}*/
				C[0][3] = -(-1);
				c[0][2] = -(-1);
				C[1][3] = -(-1);
				c[1][2] = -(-1);
				C[2][3] = -1;
				/*cout << "C Matrix :" << endl;
				for (i = 0; i < 4; i++) {
					for (j = 0; j < 4; j++) {
						cout << C[i][j] << " ";
					}
					cout << endl;
				}
				cout << "c Matrix :" << endl;
				for (i = 0; i < 3; i++) {
					for (j = 0; j < 3; j++) {
						cout << c[i][j] << " ";
					}
					cout << endl;
				}*/
				MulMatrix(B, C);
				MulMatrix3x3(b, c);
				MulMatrix(A, C);/*C為最終結果*/
				MulMatrix3x3(a, c);
				curr = head;
			}
			else if (buffer[0] == 'd') {
				glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
				glClear(GL_COLOR_BUFFER_BIT);
				glFlush();
				int a;
				bool reverse;

				reverse = false;
				curr = head;
				pos = start;
				dcount = a = 0;
				if (vicount >= 2) {

				}
				if (EMmul == true) {
					while (curr->next != NULL) {
						for (i = 0; i < curr->dn; i++) {
							
							DataMulMatrix(curr->dot[i][0], curr->dot[i][1], curr->dot[i][2], curr->dot[i][3], PM);
							/*curr->dot[i][0] /= curr->dot[i][3];
							curr->dot[i][1] /= curr->dot[i][3];
							curr->dot[i][2] /= curr->dot[i][3];
							curr->dot[i][3] /= curr->dot[i][3];*/
							//cout << "PM dot " << curr->dot[i][0] << " " << curr->dot[i][1] << " " << curr->dot[i][2] << " " << curr->dot[i][3] << endl;
						}
						curr = curr->next;
					}
				}			
				curr = head;
				pos = start;
				dcount = 0;		
				/*cout << "Window to Viewport Matrix :" << endl;
				for (i = 0; i < 4; i++) {
					for (j = 0; j < 4; j++) {
						cout << C[i][j] << " ";
					}
					cout << endl;
				}
				cout << endl;
				for (i = 0; i < 3; i++) {
					for (j = 0; j < 3; j++) {
						cout << c[i][j] << " ";
					}
					cout << endl;
				}*/
				float v1[3], v2[3], vtotal[3];
				
				vtotal[2] = -1;
				while (curr->next != NULL && pos->next != NULL) {
					for (i = 0; i < curr->pn; i++) {
						if (nobackface == true) {
							
							v1[0] = curr->dot[int(curr->plane[i][0]) - 1][0] / curr->dot[int(curr->plane[i][0]) - 1][3] - curr->dot[int(curr->plane[i][1]) - 1][0] / curr->dot[int(curr->plane[i][1]) - 1][3];
							v1[1] = curr->dot[int(curr->plane[i][0]) - 1][1] / curr->dot[int(curr->plane[i][0]) - 1][3] - curr->dot[int(curr->plane[i][1]) - 1][1] / curr->dot[int(curr->plane[i][1]) - 1][3];
							v1[2] = curr->dot[int(curr->plane[i][0]) - 1][2] / curr->dot[int(curr->plane[i][0]) - 1][3] - curr->dot[int(curr->plane[i][1]) - 1][2] / curr->dot[int(curr->plane[i][1]) - 1][3];
							
							v2[0] = curr->dot[int(curr->plane[i][2]) - 1][0] / curr->dot[int(curr->plane[i][2]) - 1][3] - curr->dot[int(curr->plane[i][1]) - 1][0] / curr->dot[int(curr->plane[i][1]) - 1][3];
							v2[1] = curr->dot[int(curr->plane[i][2]) - 1][1] / curr->dot[int(curr->plane[i][2]) - 1][3] - curr->dot[int(curr->plane[i][1]) - 1][1] / curr->dot[int(curr->plane[i][1]) - 1][3];
							v2[2] = curr->dot[int(curr->plane[i][2]) - 1][2] / curr->dot[int(curr->plane[i][2]) - 1][3] - curr->dot[int(curr->plane[i][1]) - 1][2] / curr->dot[int(curr->plane[i][1]) - 1][3];
							CrossProduct(v2, v1, vtotal);
							if (vtotal[2] >= 0)
								continue;
						}
						for (j = 0; j < curr->oneplanedot[i] - 1; j++) {
							//cout << "j = " << j << endl;
							/*取Plane的兩點做clipping*/
							/*cout << "curr->dot1[" << int(curr->plane[i][j]) - 1 << "][0] = " << curr->dot[int(curr->plane[i][j]) - 1][0] << endl;
							cout << "curr->dot1[" << int(curr->plane[i][j]) - 1 << "][1] = " << curr->dot[int(curr->plane[i][j]) - 1][1] << endl;
							cout << "curr->dot1[" << int(curr->plane[i][j]) - 1 << "][2] = " << curr->dot[int(curr->plane[i][j]) - 1][2] << endl;
							cout << "curr->dot1[" << int(curr->plane[i][j]) - 1 << "][3] = " << curr->dot[int(curr->plane[i][j]) - 1][3] << endl;
							cout << "curr->dot2[" << int(curr->plane[i][j + 1]) - 1 << "][0] = " << curr->dot[int(curr->plane[i][j + 1]) - 1][0] << endl;
							cout << "curr->dot2[" << int(curr->plane[i][j + 1]) - 1 << "][1] = " << curr->dot[int(curr->plane[i][j + 1]) - 1][1] << endl;
							cout << "curr->dot2[" << int(curr->plane[i][j + 1]) - 1 << "][2] = " << curr->dot[int(curr->plane[i][j + 1]) - 1][2] << endl;
							cout << "curr->dot2[" << int(curr->plane[i][j + 1]) - 1 << "][3] = " << curr->dot[int(curr->plane[i][j + 1]) - 1][3] << endl;*/
							//cout << endl;
							
							CheckClipping(a, pos, dcount, curr->dot[int(curr->plane[i][j]) - 1][0], curr->dot[int(curr->plane[i][j]) - 1][1], curr->dot[int(curr->plane[i][j]) - 1][2], curr->dot[int(curr->plane[i][j]) - 1][3], curr->dot[int(curr->plane[i][j + 1]) - 1][0], curr->dot[int(curr->plane[i][j + 1]) - 1][1], curr->dot[int(curr->plane[i][j + 1]) - 1][2], curr->dot[int(curr->plane[i][j + 1]) - 1][3]);														
							//cout << "hello " << a << endl;
							if (a == 0) {
								/*cout << "pos->dot1[" << dcount - 1 << "][0] = " << pos->dot1[dcount - 1][0] << endl;
								cout << "pos->dot1[" << dcount - 1 << "][1] = " << pos->dot1[dcount - 1][1] << endl;
								cout << "pos->dot1[" << dcount - 1 << "][2] = " << pos->dot1[dcount - 1][2] << endl;
								cout << "pos->dot1[" << dcount - 1 << "][3] = " << pos->dot1[dcount - 1][3] << endl;
								cout << "pos->dot2[" << dcount - 1 << "][0] = " << pos->dot2[dcount - 1][0] << endl;
								cout << "pos->dot2[" << dcount - 1 << "][1] = " << pos->dot2[dcount - 1][1] << endl;
								cout << "pos->dot2[" << dcount - 1 << "][2] = " << pos->dot2[dcount - 1][2] << endl;
								cout << "pos->dot2[" << dcount - 1 << "][3] = " << pos->dot2[dcount - 1][3] << endl;*/
								//WindowToViewport(c, pos->dot1[dcount - 1][0], pos->dot1[dcount - 1][1]);
								DataMulMatrix(pos->dot1[dcount - 1][0], pos->dot1[dcount - 1][1], pos->dot1[dcount - 1][2], pos->dot1[dcount - 1][3], C);
								/*cout << "cpos->dot1[" << dcount - 1 << "][0] = " << pos->dot1[dcount - 1][0] << endl;
								cout << "cpos->dot1[" << dcount - 1 << "][1] = " << pos->dot1[dcount - 1][1] << endl;
								cout << "cpos->dot1[" << dcount - 1 << "][2] = " << pos->dot1[dcount - 1][2] << endl;
								cout << "cpos->dot1[" << dcount - 1 << "][3] = " << pos->dot1[dcount - 1][3] << endl;
								cout << "cpos->dot2[" << dcount - 1 << "][0] = " << pos->dot2[dcount - 1][0] << endl;
								cout << "cpos->dot2[" << dcount - 1 << "][1] = " << pos->dot2[dcount - 1][1] << endl;
								cout << "cpos->dot2[" << dcount - 1 << "][2] = " << pos->dot2[dcount - 1][2] << endl;
								cout << "cpos->dot2[" << dcount - 1 << "][3] = " << pos->dot2[dcount - 1][3] << endl;*/
								
								DataMulMatrix(pos->dot2[dcount - 1][0], pos->dot2[dcount - 1][1], pos->dot2[dcount - 1][2], pos->dot2[dcount - 1][3], C);
								//WindowToViewport(c, pos->dot2[dcount - 1][0], pos->dot2[dcount - 1][1]);
								if ((pos->dot2[dcount - 1][1]  - pos->dot1[dcount - 1][1]) / (pos->dot2[dcount - 1][0] - pos->dot1[dcount - 1][0] ) <= 1 && (pos->dot2[dcount - 1][1] - pos->dot1[dcount - 1][1]) / (pos->dot2[dcount - 1][0] - pos->dot1[dcount - 1][0]) >= -1){
									reverse = false;
									DrawLines(reverse, pos->dot1[dcount - 1][0], pos->dot1[dcount - 1][1], pos->dot2[dcount - 1][0], pos->dot2[dcount - 1][1]);
									//system("pause");
								}

								else {
									reverse = true;
									DrawLines(reverse, pos->dot1[dcount - 1][1], pos->dot1[dcount - 1][0], pos->dot2[dcount - 1][1], pos->dot2[dcount - 1][0]);	
									
									//system("pause");
								}
								
							}
							a = 0;
							
						}	
						//cout << "done" << endl;
						if (vtotal[2] < 0) {
							
							CheckClipping(a, pos, dcount, curr->dot[int(curr->plane[i][curr->oneplanedot[i] - 1]) - 1][0], curr->dot[int(curr->plane[i][curr->oneplanedot[i] - 1]) - 1][1], curr->dot[int(curr->plane[i][curr->oneplanedot[i] - 1]) - 1][2], curr->dot[int(curr->plane[i][curr->oneplanedot[i] - 1]) - 1][3], curr->dot[int(curr->plane[i][0]) - 1][0], curr->dot[int(curr->plane[i][0]) - 1][1], curr->dot[int(curr->plane[i][0]) - 1][2], curr->dot[int(curr->plane[i][0]) - 1][3]);
							if (a == 0) {
								/*cout << "pos->dot1[" << dcount - 1 << "][0] = " << pos->dot1[dcount - 1][0] << endl;
								cout << "pos->dot1[" << dcount - 1 << "][1] = " << pos->dot1[dcount - 1][1] << endl;
								cout << "pos->dot1[" << dcount - 1 << "][2] = " << pos->dot1[dcount - 1][2] << endl;
								cout << "pos->dot1[" << dcount - 1 << "][3] = " << pos->dot1[dcount - 1][3] << endl;
								cout << "pos->dot2[" << dcount - 1 << "][0] = " << pos->dot2[dcount - 1][0] << endl;
								cout << "pos->dot2[" << dcount - 1 << "][1] = " << pos->dot2[dcount - 1][1] << endl;
								cout << "pos->dot2[" << dcount - 1 << "][2] = " << pos->dot2[dcount - 1][2] << endl;
								cout << "pos->dot2[" << dcount - 1 << "][3] = " << pos->dot2[dcount - 1][3] << endl;*/
								//WindowToViewport(c, pos->dot1[dcount - 1][0], pos->dot1[dcount - 1][1]);
								DataMulMatrix(pos->dot1[dcount - 1][0], pos->dot1[dcount - 1][1], pos->dot1[dcount - 1][2], pos->dot1[dcount - 1][3], C);
								/*cout << "cpos->dot1[" << dcount - 1 << "][0] = " << pos->dot1[dcount - 1][0] << endl;
								cout << "cpos->dot1[" << dcount - 1 << "][1] = " << pos->dot1[dcount - 1][1] << endl;
								cout << "cpos->dot1[" << dcount - 1 << "][2] = " << pos->dot1[dcount - 1][2] << endl;
								cout << "cpos->dot1[" << dcount - 1 << "][3] = " << pos->dot1[dcount - 1][3] << endl;
								cout << "cpos->dot2[" << dcount - 1 << "][0] = " << pos->dot2[dcount - 1][0] << endl;
								cout << "cpos->dot2[" << dcount - 1 << "][1] = " << pos->dot2[dcount - 1][1] << endl;
								cout << "cpos->dot2[" << dcount - 1 << "][2] = " << pos->dot2[dcount - 1][2] << endl;
								cout << "cpos->dot2[" << dcount - 1 << "][3] = " << pos->dot2[dcount - 1][3] << endl;*/
								//WindowToViewport(c, pos->dot2[dcount - 1][0], pos->dot2[dcount - 1][1]);
								DataMulMatrix(pos->dot2[dcount - 1][0], pos->dot2[dcount - 1][1], pos->dot2[dcount - 1][2], pos->dot2[dcount - 1][3], C);
								if ((pos->dot2[dcount - 1][1] - pos->dot1[dcount - 1][1]) / (pos->dot2[dcount - 1][0] - pos->dot1[dcount - 1][0]) <= 1 && (pos->dot2[dcount - 1][1] - pos->dot1[dcount - 1][1]) / (pos->dot2[dcount - 1][0] - pos->dot1[dcount - 1][0]) >= -1) {
									reverse = false;
									DrawLines(reverse, pos->dot1[dcount - 1][0], pos->dot1[dcount - 1][1], pos->dot2[dcount - 1][0], pos->dot2[dcount - 1][1]);
									//system("pause");
								}
								else {
									reverse = true;
									DrawLines(reverse, pos->dot1[dcount - 1][1], pos->dot1[dcount - 1][0], pos->dot2[dcount - 1][1], pos->dot2[dcount - 1][0]);	
									
								}

							}
							a = 0;
						}
					}
					/*for (i = 0; i < dcount; i++) {
						cout << i << endl;
						cout << "dot1 " << pos->dot1[i][0] << " " << pos->dot1[i][1] << " " << pos->dot1[i][2] << " " << pos->dot1[i][3] << endl;
						cout << "dot2 " << pos->dot2[i][0] << " " << pos->dot2[i][1] << " " << pos->dot2[i][2] << " " << pos->dot2[i][3] << endl;
					}*/

					curr = curr->next;
					pos = pos->next;
					dcount = 0;
					
				}
				system("pause");
				//cout << "HaHa" << endl;
			}
			else if (buffer[0] == 'n') {
				nobackface = true;
			}
			else if (buffer[0] == 'e') {
				exit(0);
			}
		}
	}
}

int main(int argc, char** argv)//这是使用glut库函数进行窗口管理
{
	system("pause");
	char* str;
	obcount = vicount = 0;
	str = argv[1];
	glutInit(&argc, argv);//使用glut库需要进行初始化
	glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);//设定窗口显示模式，颜色模型和缓存，这里是RGB颜色模型和单缓存
	glutInitWindowPosition(100, 100);//设定窗口的初始位置，屏幕左上角为原点，单位为像素
	displayFunc(str);
	glutMainLoop();//使窗口框架运行起来，使显示回调函数开始工作

	return 0;
}