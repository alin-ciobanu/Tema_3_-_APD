#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

#define NUM_COLORS 256

float absComplex (float re, float im) {

	return sqrtf(re * re + im * im);

}

void multiplyComplex (float x_re, float x_im, float y_re, float y_im, float* res_re, float* res_im) {

	float x = x_re * y_re - x_im * y_im;
	float y = x_re * y_im + x_im * y_re;
	*res_re = x;
	*res_im = y;

}

void addComplex (float x_re, float x_im, float y_re, float y_im, float* res_re, float* res_im) {

	float x = x_re + y_re;
	float y = x_im + y_im;
	*res_re = x;
	*res_im = y;

}

float absfloat (float x) {
	if (x >= 0)
		return x;
	else
		return -x;
}

int main(int argc, char **argv) {

	int rank, size;
	int tag = 1;

	MPI_Status status;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	FILE* in;
    FILE* out;

    int setType;
	float x_min, x_max, y_min, y_max;
	float res;
	int MAX_STEPS;
	float complexJuliaX, complexJuliaY;

    if (rank == 0) {

		in = fopen(argv[1], "r");
		
		fscanf(in, "%d", &setType);
		fscanf(in, "%f", &x_min);
		fscanf(in, "%f", &x_max);
		fscanf(in, "%f", &y_min);
		fscanf(in, "%f", &y_max);
		fscanf(in, "%f", &res);
		fscanf(in, "%d", &MAX_STEPS);
		if (setType == 1) {
			fscanf(in, "%f", &complexJuliaX);
			fscanf(in, "%f", &complexJuliaY);
		}

		int i;
		for (i = 1; i < size; i++) {
			MPI_Send(&x_min, 1, MPI_FLOAT, i, tag, MPI_COMM_WORLD);
			MPI_Send(&x_max, 1, MPI_FLOAT, i, tag, MPI_COMM_WORLD);
			MPI_Send(&y_min, 1, MPI_FLOAT, i, tag, MPI_COMM_WORLD);
			MPI_Send(&y_max, 1, MPI_FLOAT, i, tag, MPI_COMM_WORLD);
			MPI_Send(&res, 1, MPI_FLOAT, i, tag, MPI_COMM_WORLD);
			MPI_Send(&MAX_STEPS, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
			MPI_Send(&setType, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
			if (setType == 1) {
				MPI_Send(&complexJuliaX, 1, MPI_FLOAT, i, tag, MPI_COMM_WORLD);
				MPI_Send(&complexJuliaY, 1, MPI_FLOAT, i, tag, MPI_COMM_WORLD);
			}
		}

	}
	else {
		// rank != 0

		MPI_Recv(&x_min, 1, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &status);
		MPI_Recv(&x_max, 1, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &status);
		MPI_Recv(&y_min, 1, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &status);
		MPI_Recv(&y_max, 1, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &status);
		MPI_Recv(&res, 1, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &status);
		MPI_Recv(&MAX_STEPS, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
		MPI_Recv(&setType, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
		if (setType == 1) {
			MPI_Recv(&complexJuliaX, 1, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &status);
			MPI_Recv(&complexJuliaY, 1, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &status);
		}
	}

	int height = (int) floor((y_max - y_min) / res);
	int width = (int) floor((x_max - x_min) / res);
	float rankFloat = (float) rank;
	float sizeFloat = (float) size;

	if (setType == 0) {

		float cx, cy;
		int color;

		for (cx = x_min + rankFloat / sizeFloat * absfloat(x_max - x_min); 
			cx < x_min + (rankFloat + 1) / sizeFloat * absfloat(x_max - x_min); cx += res) {

			int col = (int) floor((cx - x_min) / res);
			int* colors = (int*) malloc ((height + 1) * sizeof(int));
			int i = 0;
			// la inceputul vectorului se pune coloana pe care se va scrie
			colors[i] = col;
			i++;


			for (cy = y_min; cy < y_max; cy += res) {

				float zx, zy;
				zx = zy = 0;
				int step = 0;

				while (absComplex(zx, zy) < 2.0 && step < MAX_STEPS) {

					float zx2, zy2;
					multiplyComplex(zx, zy, zx, zy, &zx2, &zy2);
					addComplex(zx2, zy2, cx, cy, &zx, &zy);
					step++;

				}

				color = step % NUM_COLORS;
				colors[i] = color;
				i++;

			}
			
			int next = 1;
			MPI_Send(&next, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
			MPI_Send(colors, height + 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
			free(colors);

		}
	}
	else if (setType == 1) {

		float cx, cy;
		int color;

		for (cx = x_min + rankFloat / sizeFloat * absfloat(x_max - x_min); 
			cx < x_min + (rankFloat + 1) / sizeFloat * absfloat(x_max - x_min); cx += res) {

			int col = (int) floor((cx - x_min) / res);
			int* colors = (int*) malloc ((height + 1) * sizeof(int));
			int i = 0;
			// la inceputul vectorului se pune coloana pe care se va scrie
			colors[i] = col;
			i++;


			for (cy = y_min; cy < y_max; cy += res) {

				float zx, zy;;
				zx = cx;
				zy = cy;
				int step = 0;

				while (absComplex(zx, zy) < 2.0 && step < MAX_STEPS) {

					float zx2, zy2;
					multiplyComplex(zx, zy, zx, zy, &zx2, &zy2);
					addComplex(zx2, zy2, complexJuliaX, complexJuliaY, &zx, &zy);
					step++;

				}

				color = step % NUM_COLORS;
				colors[i] = color;
				i++;

			}
			
			int next = 1;
			MPI_Send(&next, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
			MPI_Send(colors, height + 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
			free(colors);

		}
	}
	
	int next = 0;
	// semnalizeaza sfarsitul -- procesul nu mai trimite date
	MPI_Send(&next, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);

	int **image_matrix = (int**) malloc(height * sizeof(int*));
	int i;
	for (i = 0; i < width; i++) {
		image_matrix[i] = (int*) malloc (width * sizeof(int));
	}

	if (rank == 0) {

		for (i = 0; i < size; i++) {

			next = 2;
			MPI_Recv(&next, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &status);

			while (next == 1) {
				// has next
				int* colors = (int*) malloc ((height + 1) * sizeof(int));
				MPI_Recv(colors, height + 1, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
				int col = colors[0];
				int k;
				for (k = 1; k <= height; k++) {
					image_matrix[k - 1][col] = colors[k];
				}
				MPI_Recv(&next, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
				free(colors);
			}
		}

		out = fopen(argv[2], "w");
		fprintf(out, "P2\n");
		fprintf(out, "%d %d\n", width, height);
		fprintf(out, "%d\n", NUM_COLORS);
		int k;
		for (i = height - 1; i > 0; i--) {
			for (k = 0; k < width; k++) {
				fprintf(out, "%d ", image_matrix[i][k]);
			}
			if (i != height - 1)
				fprintf(out, "\n");
		}
		fclose(out);

	}

	MPI_Finalize();
	return 0;

}
