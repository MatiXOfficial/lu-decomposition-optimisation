// 9: Blocked the matrix A; block size ABlockSize x N

#define min( i, j ) ( (i)<(j) ? (i): (j) )
#define max( i, j ) ( (i)>(j) ? (i): (j) )

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>
#include <x86intrin.h>

#include "configuration.h"

#define IDX(i, j, n) (((j) + (i) * (n)))

// Block sizes 
#define ABlockSize 200

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

void innerKernel(register double* A, register int N, register int ni, register double* ABlock, register int blockEnd)
{
    register int i, j, k;
    register double A_i_i, A_j_i;

    register __m256d mm_A_j_i;
    register __m256d tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;

    for (i = 0; i < ni; i++) // i will point to A
    {
        A_i_i = A[IDX(i, i, N)];
        for (j = max(0, i + 1 - ni + blockEnd); j < blockEnd; j++) // j will point to ABlock (j - only rows)
        {
            A_j_i = ABlock[IDX(j, i, N)] / A_i_i;
            ABlock[IDX(j, i, N)] = A_j_i;
            mm_A_j_i[0] = A_j_i;
            mm_A_j_i[1] = A_j_i;
            mm_A_j_i[2] = A_j_i;
            mm_A_j_i[3] = A_j_i;

            for (k = i + 1; k < N;) // k - only columns
            {
                if (k + 16 < N)
                {
                    tmp0 = _mm256_loadu_pd(ABlock + IDX(j, k, N));
                    tmp1 = _mm256_loadu_pd(A + IDX(i, k, N));
                    tmp2 = _mm256_loadu_pd(ABlock + IDX(j, k+4, N));
                    tmp3 = _mm256_loadu_pd(A + IDX(i, k+4, N));
                    tmp4 = _mm256_loadu_pd(ABlock + IDX(j, k+8, N));
                    tmp5 = _mm256_loadu_pd(A + IDX(i, k+8, N));
                    tmp6 = _mm256_loadu_pd(ABlock + IDX(j, k+12, N));
                    tmp7 = _mm256_loadu_pd(A + IDX(i, k+12, N));

                    tmp0 -= mm_A_j_i * tmp1;
                    tmp2 -= mm_A_j_i * tmp3;
                    tmp4 -= mm_A_j_i * tmp5;
                    tmp6 -= mm_A_j_i * tmp7;

                    _mm256_storeu_pd(ABlock + IDX(j, k, N), tmp0);
                    _mm256_storeu_pd(ABlock + IDX(j, k+4, N), tmp2);
                    _mm256_storeu_pd(ABlock + IDX(j, k+8, N), tmp4);
                    _mm256_storeu_pd(ABlock + IDX(j, k+12, N), tmp6);

                    k += 16;
                }
                else
                {
                    ABlock[IDX(j, k, N)] -= A_j_i * A[IDX(i, k, N)];
                    k++;
                }
            }
        }
    }
}

/* INPUT: A - array of pointers to rows of a square matrix having dimension N
 *        Tol - small tolerance number to detect failure when the matrix is near degenerate
 * OUTPUT: Matrix A is changed, it contains a copy of both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
 *        The permutation matrix is not stored as a matrix, but in an integer vector P of size N+1 
 *        containing column indexes where the permutation matrix has "1". The last element P[N]=S+N, 
 *        where S is the number of row exchanges needed for determinant computation, det(P)=(-1)^S    
 */
int LUPDecompose(register double* A, register int N)
{
    register int blockStart, blockEnd;

    for (blockStart = 0; blockStart < N; blockStart += ABlockSize)
    {
        blockEnd = min(N - blockStart, ABlockSize); // not to reach beyond matrix
        innerKernel(A, N, blockStart + blockEnd, A + IDX(blockStart, 0, N), blockEnd);
    }

    return 1; //decomposition done
}

int main(int argc, const char* argv[])
{
    char fileName[128];
    if (argc == 1)
    {
        sprintf(fileName, "results/lu9_results.txt");
    }
    else
    {
        sprintf(fileName, "results/lu9_results_%s.txt", argv[1]);
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
