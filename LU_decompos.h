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
	//n - количество элементов в строке матрицы
	//m - полуширина ленты
	int n, m;
	//текущее и предыдущее значение х
	double x_now, x_previous;
	//текущее и предыдущее значениие времени
	double t_now, t_previous;
	//мтрицы Л и У, и вектора D, b, у, х необходимые для ЛУ решения
	real **matrix_L, *vector_D, **matrix_U, *vector_b, *vector_y, *vector_x;
	//предыдущее значение q
	real *q_previous;
	//конструктор
	Matrix();
	//ЛУ разложение
	void LU_decomposition( );
	//прямой ход
	void Ly_is_b( );
	//обратный ход
	void Ux_is_y( );
	//инициализация векторов перед решателем
	void vect( );
	//считка матрицы с файла
	void enter_matrix(FILE* inp_file);
	//локальная матрица А
	void Local_A();
	//вычисление локальной марицы G
	void calc_G();
	//лямбда функция
	double lambda(double x);
	//вычисление локальной марицы М
	void calc_M();
	//сигма функция
	double sigma(double x);
	//локальный вектор b
	void local_b();
	//функция (первое краевое условие)
	double Function(double x);
	//умножение матрицы на  вектор
	real* Mult(Matrix M, real *&vector);
	//умножение матрицы на число
	void Mult(Matrix M, real v);
	//сумма матриц В и С
	void Sum(Matrix A, Matrix B, Matrix C);
};

