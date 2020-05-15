#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MIN(a, b) ((a) < (b) ? (a) : (b))

int main(int argc, char *argv[])
{
    int count;           /* Local prime count */
    double elapsed_time; /* Parallel execution time */
    int first;           /* Index of first multiple */
    int global_count;    /* Global prime count */
    int high_value;      /* Highest value on this proc */
    int i;
    int id;    /* Process ID number */
    int index; /* Index of current prime */
    int base_index;
    int low_value; /* Lowest value on this proc */
    char *marked;  /* Portion of 2,...,'n' */
    char *base_marked;
    int n; /* Sieving from 2, ..., 'n' */
    int N;
    int p;          /* Number of processes */
    int proc0_size; /* Size of proc 0's subarray */
    int prime;      /* Current prime */
    int size;       /* Elements in 'marked' */
    int base_size;
    int base_low_value;
    int base_high_value;
    MPI_Init(&argc, &argv);

    /* Start the timer */

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    if (argc != 2)
    {
        if (!id)
            printf("Command line: %s <m>\n", argv[0]);
        MPI_Finalize();
        exit(1);
    }

    N = atoi(argv[1]);
    n = N % 2 ? N : N - 1;
    // printf("n= %d,N=%d\n",n,N);
    // n = N % 2 ? N - 1 : N;
    /* Figure out this process's share of the array, as
      well as the integers represented by the first and
      last array elements */

    low_value = 3 + (id * (n - 1) / 2 / p) * 2;
    high_value = 1 + ((id + 1) * (n - 1) / 2 / p) * 2;
    size = (high_value - low_value) / 2 + 1;
    base_low_value = 3;
    base_high_value = 1 + ((n - 1) / 2 / p) * 2;
    base_size = (int)sqrt((double)n);
    // printf("id: %d, low_value: %d, high_value: %d\n", id, low_value, high_value);
    // printf("id: %d, base_low_value: %d, base_high_value: %d\n", id, base_low_value, base_high_value);

    /* Bail out if all the primes used for sieving are
      not all held by process 0 */

    proc0_size = (n - 1) / p;

    if ((2 + proc0_size) < (int)sqrt((double)n))
    {
        if (!id)
            printf("Too many processes\n");
        MPI_Finalize();
        exit(1);
    }

    /* Allocate this process's share of the array. */

    marked = (char *)malloc(size);
    base_marked = (char *)malloc(base_size);
    if (marked == NULL || base_marked == NULL)
    {
        printf("Cannot allocate enough memory\n");
        MPI_Finalize();
        exit(1);
    }

    for (i = 0; i < size; i++)
        marked[i] = 0;
    for (i = 0; i < base_size; i++)
        base_marked[i] = 0;
    // if (!id)
    //     index = 0;
    index = 0;
    base_index = 0;
    prime = 3;
    do
    {
        // printf("prime: %d", prime);
        //find first index
        if (prime * prime > base_low_value)
            first = (prime * prime - base_low_value) / 2;
        else
        {
            if (!(base_low_value % prime))
                first = 0;
            else
            {
                int remain = base_low_value % prime;
                first = (prime * (1 + (remain + 1) % 2) - remain) / 2;
            }
        }
        // sieve
        for (i = first; i < base_size; i += prime)
            base_marked[i] = 1;
        while (base_marked[++base_index])
            ;
        prime = 2 * base_index + 3;
    } while (prime * prime <= n);
    prime = 3;
    base_index = 0;
    do
    {
        //find first index
        if (prime * prime > low_value)
            first = (prime * prime - low_value) / 2;
        else
        {
            if (!(low_value % prime))
                first = 0;
            else
            {
                int remain = low_value % prime;
                first = (prime * (1 + (remain + 1) % 2) - remain) / 2;
            }
        }
        // sieve
        for (i = first; i < size; i += prime)
            marked[i] = 1;
        while (base_marked[++base_index])
            ;
        prime = 2 * base_index + 3;
    } while (prime * prime <= n);
    /*  printf("id = %d\n", id);
    for (int j = 0; j < size; j++)
    {
        printf("%d ", marked[j]);
    }
    printf("\n"); */
    count = 0;
    for (i = 0; i < size; i++)
        if (!marked[i])
            count++;
    if (p > 1)
        MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM,
                   0, MPI_COMM_WORLD);

    /* Stop the timer */

    elapsed_time += MPI_Wtime();

    /* Print the results */

    if (!id)
    {
        if (p == 1)
        {
            global_count = count;
        }
        printf("There are %d primes less than or equal to %d\n",
               global_count + 1, N);
        printf("SIEVE (%d) %10.6f\n", p, elapsed_time);
    }
    MPI_Finalize();
    return 0;
}