// 7: AVX

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>
#include <x86intrin.h>

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

    register __m256d mm_A_j_i;
    register __m256d tmp0, tmp1, tmp2, tmp3;

    for (i = 0; i < N; i++)
    {
        A_i_i = A[IDX(i, i, N)];
        for (j = i + 1; j < N; j++)
        {
            A_j_i = A[IDX(j, i, N)] / A_i_i;
            A[IDX(j, i, N)] = A_j_i;
            mm_A_j_i[0] = A_j_i;
            mm_A_j_i[1] = A_j_i;
            mm_A_j_i[2] = A_j_i;
            mm_A_j_i[3] = A_j_i;

            for (k = i + 1; k < N;)
            {
                if (k + 8 < N)
                {
                    tmp0 = _mm256_loadu_pd(A + IDX(j, k, N));
                    tmp1 = _mm256_loadu_pd(A + IDX(i, k, N));
                    tmp2 = _mm256_loadu_pd(A + IDX(j, k+4, N));
                    tmp3 = _mm256_loadu_pd(A + IDX(i, k+4, N));

                    tmp0 -= mm_A_j_i * tmp1;
                    tmp2 -= mm_A_j_i * tmp3;

                    _mm256_storeu_pd(A + IDX(j, k, N), tmp0);
                    _mm256_storeu_pd(A + IDX(j, k+4, N), tmp2);

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
        sprintf(fileName, "results/lu7_results.txt");
    }
    else
    {
        sprintf(fileName, "results/lu7_results_%s.txt", argv[1]);
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
        int* P = (int*)malloc((size + 1) * sizeof(int));

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
