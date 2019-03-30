#include"LU_decompos.h"

void by_hand(); // ввод имени файла руками
//void Ak(); // матрицы Ak
//void Hylb(); // матрицы Гильберта
//void Gauss();

int main() {
	by_hand();
}

void by_hand() {

	Matrix A1;
	//начально не изменяемая матрица G
	Matrix G_;

	FILE* input_file_A;
	FILE* input_file_B;

	fopen_s(&input_file_A, "input.txt" , "r");
	G_.enter_matrix(input_file_A);
	fclose(input_file_A);

	A1.vect();
	A1.LU_decomposition( );
	A1.Ly_is_b( );
	A1.Ux_is_y( );


	fopen_s(&input_file_B, "answer.txt", "w");
	for (int i = 0; i < A1.n ; i++)
	{
		fprintf(input_file_B, outform, A1.vector_x[i]);
	}
	fprintf(input_file_B,"\n");
	for (int i = 0; i < A1.n; i++)
	{
		A1.vector_x[i] = (i + 1) - A1.vector_x[i];
		fprintf(input_file_B, outformP, A1.vector_x[i]);
	}
	fclose(input_file_B);
}


