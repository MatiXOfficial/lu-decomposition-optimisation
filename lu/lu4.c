// 4: Loop unrolling (8 iterations)

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>

#include "configuration.h"

static double gtod_ref_time_sec = 0.0;

double dclock()
{
    double the_time, norm_sec;
    struct timeval tv;
    gettimeofday(&tv, NULL);
    if (gtod_ref_time_sec == 0.0)
        gtod_ref_time_sec = (double)tv.tv_sec;
    norm_sec = (double)tv.tv_sec - gtod_ref_time_sec;
    the_time = norm_sec + tv.tv_usec * 1.0e-6;
    return the_time;
}

int LUPDecompose(register double** A, register int N)
{
    register int i, j, k;
    register double A_i_i, A_j_i;

    for (i = 0; i < N; i++)
    {
        A_i_i = A[i][i];
        for (j = i + 1; j < N; j++)
        {
            A_j_i = A[j][i] / A_i_i;
            A[j][i] = A_j_i;

            for (k = i + 1; k < N;)
            {
                if (k + 8 < N)
                {
                    A[j][k] -= A_j_i * A[i][k];
                    A[j][k+1] -= A_j_i * A[i][k+1];
                    A[j][k+2] -= A_j_i * A[i][k+2];
                    A[j][k+3] -= A_j_i * A[i][k+3];
                    A[j][k+4] -= A_j_i * A[i][k+4];
                    A[j][k+5] -= A_j_i * A[i][k+5];
                    A[j][k+6] -= A_j_i * A[i][k+6];
                    A[j][k+7] -= A_j_i * A[i][k+7];

                    k += 8;
                }
                else
                {
                    A[j][k] -= A_j_i * A[i][k];
                    k++;
                }
            }
        }
    }

    return 1; //decomposition done
}

int main(int argc, const char* argv[])
{
    char fileName[128];
    if (argc == 1)
    {
        sprintf(fileName, "results/lu4_results.txt");
    }
    else
    {
        sprintf(fileName, "results/lu4_results_%s.txt", argv[1]);
    }
    FILE *fptr = fopen(fileName, "w");

    fprintf(fptr, "Size");
    for (int i_run = 0; i_run < N_RUNS; i_run++)
    {
        fprintf(fptr, ",Time%d", i_run);
    }
    fprintf(fptr, "\n");

    int size;
    for (size = SIZE_START; size <= SIZE_END; size += SIZE_STEP)
    {
        if (PRINT_LOGS)
            printf("=============================\n");
        printf("SIZE: %d\n", size);
        if (PRINT_LOGS)
            printf("=============================\n");
        fprintf(fptr, "%d", size);

        int i, j, retVal;
        double dtime;

        double** matrix = (double**)malloc(size * sizeof(double*));
        double* matrix_ = (double*)malloc(size * size * sizeof(double));
        for (int i = 0; i < size; i++)
        {
            matrix[i] = matrix_ + i * size;
        }

        srand(1);
        for (i = 0; i < size; i++)
        {
            for (j = 0; j < size; j++)
            {
                matrix[i][j] = rand();
            }
        }

        for (int i_run = 0; i_run < N_RUNS; i_run++)
        {
            if (PRINT_LOGS)
                printf("Run %d, call LUPDecompose", i_run);
            dtime = dclock();
            retVal = LUPDecompose(matrix, size);
            dtime = dclock() - dtime;
            if (retVal != 1)
            {
                printf("failure");
                fprintf(fptr, ",FAIL");
            }
            else
            {
                if (PRINT_LOGS)
                    printf("Time: %le s \n", dtime);
                fprintf(fptr, ",%le", dtime);
            }

            double check = 0.;
            for (i = 0; i < size; i++)
            {
                for (j = 0; j < size; j++)
                {
                    check += matrix[i][j];
                }
            }
            if (PRINT_LOGS)
                printf("Check: %le \n", check);
            fflush(stdout);
        }

        fprintf(fptr, "\n");
    }
    fclose(fptr);
}
