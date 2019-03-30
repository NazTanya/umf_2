#include "LU_decompos.h"

//начально не изменяемая матрица G
Matrix G_;
//локальная матрица А
Matrix A;
//локальная матрица G
Matrix G;
//локальная матрица М
Matrix M;

Matrix::Matrix() {
	vector_D = new real[n];
	vector_b = new real[n];
	vector_x = new real[n];	
	vector_y = new real[n];
	matrix_L = new real*[m];
	for (int i = 0; i < m; i++)
		matrix_L[i] = new real[n];
	matrix_U = new real*[m];
	for (int i = 0; i < m; i++)
		matrix_U[i] = new real[n];
}

double Matrix::lambda(double x) {
	return 1;/////например
}

double Matrix::sigma(double x) {
	return 1;//например
}

double Matrix::Function(double x) {
	return 3 * x; //например
}

void Matrix:: LU_decomposition( ) {
	//LU разложение
	reals su, sl, sd;
	for (int i = 0; i < n; i++)
	{
		sd = 0;
		int j = i - m;
		for (int jL = 0; jL < m; jL++, j++)
		{
			su = 0;
			sl = 0;
			if (j < 0) continue;
			int ku = i - j; 
			for (int kl = 0; kl < jL; kl++,ku++)
			{
				if (ku < 0) continue;
				sl += matrix_L[i][kl] * matrix_U[j][ku]; 
				su += matrix_L[j][ku] * matrix_U[i][kl];
			}
			
			matrix_L[i][jL] -= sl;
			matrix_U[i][jL] = (matrix_U[i][jL] - su) / vector_D[j]; 
			sd += matrix_L[i][jL] * matrix_U[i][jL];  

		}
		vector_D[i] -= sd;
	}			

}

void Matrix:: Ly_is_b( ) {                
	//прямой ход	
	for (int i = 0; i<n; i++) {
		reals sum = 0;	
		int j = i - m;
		for (int jL = 0; jL < m; jL++, j++)
		{
			if (j < 0) continue;
			sum += matrix_L[i][jL] * vector_y[j];
		}
		vector_y[i] = (vector_b[i] - sum) / vector_D[i];
	}
}

void Matrix:: Ux_is_y( ) { 
	//обратный
	for (int i = n - 1; i >= 0; i--) 
	{		
		int j = i - m;		
		real xi = vector_y[i] / 1;
		for (int ju = 0; ju < m; ju++, j++)
		{
			if (j < 0) continue;
			vector_y[j] -= matrix_U[i][ju] * xi;
		}
		vector_x[i] = xi;
	}
}

void Matrix::enter_matrix(FILE* inp_file) {
	fscanf_s(inp_file, "%d %d", &n, &m); 
 //ввели размерность матрицы n и полуширину ленты m
 //количество чисел для записи матрицы в формате
	//m = n - 1;
	matrix_L = new real *[n];
	matrix_U = new real *[n];
	vector_D = new real[n];
	vector_b = new real[n];

	for (int i = 0; i < n; i++)
	{
		matrix_L[i] = new real[m];
		matrix_U[i] = new real[m];
	}
	
	for (int i = 0; i < n; i++)
	{
		for (int jL = 0; jL < m; jL++)
		{
			fscanf_s(inp_file, inpform, &matrix_L[i][jL]);
		}
	}

	for (int i = 0; i < n; i++)
	{
		fscanf_s(inp_file, inpform, &vector_D[i]);
	}

	for (int i = 0; i < n; i++)
	{
		for (int jL = 0; jL < m; jL++)
		{
			fscanf_s(inp_file, inpform, &matrix_U[i][jL]);
		}
	}
	for (int i = 0; i < n; i++)
		fscanf_s(inp_file, inpform, &vector_b[i]);
}

void Matrix::vect() {
	vector_y = vector_x = vector_b;
}

void Matrix::calc_G() {
	double h = x_now - x_previous;
	double c = lambda(x_now) / (3 * h);
	Mult(G, c);
}

void Matrix::calc_M() {
	double coef0, coef1, coef2;
	double c = (1 / (t_now - t_previous));
	double delta = x_now - x_previous;
	coef0 = delta * sigma(x_previous);
	coef1 = delta * sigma((x_previous + x_now) / 2);
	coef2 = delta * sigma(x_now);

	M.matrix_L[0][0] = 0;
	M.matrix_L[0][1] = 0;
	M.matrix_L[1][0] = 0;
	M.matrix_U[0][0] = 0;
	M.matrix_U[0][1] = 0;
	M.matrix_U[1][0] = 0;

	M.vector_D[0] = c * (coef0 * 0.09285714285714286 + coef1 * 0.04761904761904761 + coef2 * -7.1428571428571415e-3);

	//massLocalMatrix[0][1] = coef0 * 0.04761904761904761 + coef1 * 0.0380952380952381 + coef2 * -0.01904761904761905;
	M.matrix_U[1][1] = c * (coef0 * 0.04761904761904761 + coef1 * 0.0380952380952381 + coef2 * -0.01904761904761905);

	//massLocalMatrix[0][2] = coef0 * -7.1428571428571415e-3 + coef1 * -0.01904761904761905 + coef2 * -7.1428571428571415e-3;
	M.matrix_U[2][1] = c * (coef0 * -7.1428571428571415e-3 + coef1 * -0.01904761904761905 + coef2 * -7.1428571428571415e-3);

	//massLocalMatrix[1][0] = coef0 * 0.04761904761904761 + coef1 * 0.0380952380952381 + coef2 * -0.01904761904761905;
	M.matrix_L[1][1] = c * (coef0 * 0.04761904761904761 + coef1 * 0.0380952380952381 + coef2 * -0.01904761904761905);

	M.vector_D[1] = c * (coef0 * 0.0380952380952381 + coef1 * 0.4571428571428571 + coef2 * 0.03809523809523809);

	//massLocalMatrix[1][2] = coef0 * -0.01904761904761905 + coef1 * 0.03809523809523809 + coef2 * 0.04761904761904761;
	M.matrix_U[2][0] = c * (coef0 * -0.01904761904761905 + coef1 * 0.03809523809523809 + coef2 * 0.04761904761904761);

	//massLocalMatrix[2][0] = coef0 * -7.1428571428571415e-3 + coef1 * -0.01904761904761905 + coef2 * -7.1428571428571415e-3;
	M.matrix_L[2][1] = c * (coef0 * -7.1428571428571415e-3 + coef1 * -0.01904761904761905 + coef2 * -7.1428571428571415e-3);

	//massLocalMatrix[2][1] = coef0 * -0.01904761904761905 + coef1 * 0.03809523809523809 + coef2 * 0.04761904761904761;
	M.matrix_L[2][0] = c * (coef0 * -0.01904761904761905 + coef1 * 0.03809523809523809 + coef2 * 0.04761904761904761);

	M.vector_D[2] = c * (coef0 * -7.1428571428571415e-3 + coef1 * 0.04761904761904761 + coef2 * 0.09285714285714286);
}

void Matrix::Local_A() {
	G = G_;
	calc_G();
	calc_M();
	Sum(A, M, G);
}

void Matrix::local_b()
{
	calc_M();
	double delta = t_now - t_previous;
	double c = (1 / delta);
	Mult(M, c);
	double *u = Mult(M, q_previous);
	double *f = new double[3];
	f[0] = Function(x_previous);
	f[1] = Function((x_previous + x_now) / 2);
	f[2] = Function(x_now);
	for (int i = 0; i < 3; i++)
		vector_b[i] = u[i] + f[i];
}

real* Matrix:: Mult(Matrix M, real *&vector)
{
	real *ans;
	int i, j, jL;
	ans = new real[3];
	for (i = 0; i < 3; i++)
		ans[i] = 0;
	for (i = 0; i < 3; i++)
	{
		ans[i] = M.vector_D[i] * vector[i];
		j = i - 2;
		for (jL = 0; jL < 2; jL++, j++)
		{
			if (j < 0) continue;
			ans[i] += M.matrix_L[i][jL] * vector[j];
			ans[j] += M.matrix_U[i][jL] * vector[i];
		}
	}
	return ans;
}

void Matrix::Mult(Matrix M, real v)
{
	int i, j, jL;
	for (i = 0; i < 3; i++)
	{
		M.vector_D[i] *= v;
		j = i - 2;
		for (jL = 0; jL < 2; jL++, j++)
		{
			if (j < 0) continue;
			M.matrix_L[i][jL] *= v;
			M.matrix_U[i][jL] *= v;
		}
	}
}

void Matrix::Sum(Matrix A, Matrix B, Matrix C) {
	for (int i = 0; i < 3; i++)
		A.vector_D[i] = B.vector_D[i] + C.vector_D[i];
	for(int i = 0; i < 3; i++)
		for (int j = 0; j < 2; j++)
		{
			A.matrix_L[i][j] = B.matrix_L[i][j] + C.matrix_L[i][j];
			A.matrix_U[i][j] = B.matrix_U[i][j] + C.matrix_U[i][j];
		}
}