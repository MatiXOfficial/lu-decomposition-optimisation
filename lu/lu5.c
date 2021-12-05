// 5: Replacing a 2D matrix with a 1D counterpart

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>

#include "configuration.h"

#define IDX(i, j, n) (((j) + (i) * (n)))

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

int LUPDecompose(register double* A, register int N)
{
    register int i, j, k;
    register double A_i_i, A_j_i;

    for (i = 0; i < N; i++)
    {
        A_i_i = A[IDX(i, i, N)];
        for (j = i + 1; j < N; j++)
        {
            A_j_i = A[IDX(j, i, N)] / A_i_i;
            A[IDX(j, i, N)] = A_j_i;

            for (k = i + 1; k < N;)
            {
                if (k + 8 < N)
                {
                    A[IDX(j, k, N)] -= A_j_i * A[IDX(i, k, N)];
                    A[IDX(j, k+1, N)] -= A_j_i * A[IDX(i, k+1, N)];
                    A[IDX(j, k+2, N)] -= A_j_i * A[IDX(i, k+2, N)];
                    A[IDX(j, k+3, N)] -= A_j_i * A[IDX(i, k+3, N)];
                    A[IDX(j, k+4, N)] -= A_j_i * A[IDX(i, k+4, N)];
                    A[IDX(j, k+5, N)] -= A_j_i * A[IDX(i, k+5, N)];
                    A[IDX(j, k+6, N)] -= A_j_i * A[IDX(i, k+6, N)];
                    A[IDX(j, k+7, N)] -= A_j_i * A[IDX(i, k+7, N)];

                    k += 8;
                }
                else
                {
                    A[IDX(j, k, N)] -= A_j_i * A[IDX(i, k, N)];
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
        sprintf(fileName, "results/lu5_results.txt");
    }
    else
    {
        sprintf(fileName, "results/lu5_results_%s.txt", argv[1]);
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

        double* matrix = (double*)malloc(size * size * sizeof(double*));

        srand(1);
        for (i = 0; i < size; i++)
        {
            for (j = 0; j < size; j++)
            {
                matrix[IDX(i, j, size)] = rand();
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
                    check += matrix[IDX(i, j, size)];
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
