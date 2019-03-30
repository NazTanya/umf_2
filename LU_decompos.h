#pragma once

#include <stdio.h>
#include <math.h>
#include <conio.h>

typedef double real;
typedef double reals;

#define inpform "%f "   //lf double //  f float
#define outform "%.8f\n"
#define outformP "%.2le\n"
#define outform_G "%.15lf \t %.15lf\n"

/*
massLocalMatrix[0][0] = coef0 * 0.09285714285714286 + coef1 * 0.04761904761904761 + coef2 * -7.1428571428571415e-3;
massLocalMatrix[0][1] = coef0 * 0.04761904761904761 + coef1 * 0.0380952380952381 + coef2 * -0.01904761904761905;
massLocalMatrix[0][2] = coef0 * -7.1428571428571415e-3 + coef1 * -0.01904761904761905 + coef2 * -7.1428571428571415e-3;

massLocalMatrix[1][0] = coef0 * 0.04761904761904761 + coef1 * 0.0380952380952381 + coef2 * -0.01904761904761905;
massLocalMatrix[1][1] = coef0 * 0.0380952380952381 + coef1 * 0.4571428571428571 + coef2 * 0.03809523809523809;
massLocalMatrix[1][2] = coef0 * -0.01904761904761905 + coef1 * 0.03809523809523809 + coef2 * 0.04761904761904761;

massLocalMatrix[2][0] = coef0 * -7.1428571428571415e-3 + coef1 * -0.01904761904761905 + coef2 * -7.1428571428571415e-3;
massLocalMatrix[2][1] = coef0 * -0.01904761904761905 + coef1 * 0.03809523809523809 + coef2 * 0.04761904761904761;
massLocalMatrix[2][2] = coef0 * -7.1428571428571415e-3 + coef1 * 0.04761904761904761 + coef2 * 0.09285714285714286;
*/

class Matrix {
public:
	//n - ���������� ��������� � ������ �������
	//m - ���������� �����
	int n, m;
	//������� � ���������� �������� �
	double x_now, x_previous;
	//������� � ���������� ��������� �������
	double t_now, t_previous;
	//������ � � �, � ������� D, b, �, � ����������� ��� �� �������
	real **matrix_L, *vector_D, **matrix_U, *vector_b, *vector_y, *vector_x;
	//���������� �������� q
	real *q_previous;
	//�����������
	Matrix();
	//�� ����������
	void LU_decomposition( );
	//������ ���
	void Ly_is_b( );
	//�������� ���
	void Ux_is_y( );
	//������������� �������� ����� ���������
	void vect( );
	//������ ������� � �����
	void enter_matrix(FILE* inp_file);
	//��������� ������� �
	void Local_A();
	//���������� ��������� ������ G
	void calc_G();
	//������ �������
	double lambda(double x);
	//���������� ��������� ������ �
	void calc_M();
	//����� �������
	double sigma(double x);
	//��������� ������ b
	void local_b();
	//������� (������ ������� �������)
	double Function(double x);
	//��������� ������� ��  ������
	real* Mult(Matrix M, real *&vector);
	//��������� ������� �� �����
	void Mult(Matrix M, real v);
	//����� ������ � � �
	void Sum(Matrix A, Matrix B, Matrix C);
};

