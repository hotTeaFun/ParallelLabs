#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define llu unsigned long long
#define coreCount 14        //14 core
#define l3CacheSize 36700160 //35MB 35*2^20=36700160
#define cacheLine 64         //64B

struct aligned_llu
{
    llu value;
    char padding[cacheLine];
};
int main(int argc, char *argv[])
{
    llu count;                /* Local prime count */
    double elapsed_time;      /* Parallel execution time */
    struct aligned_llu first; /* Index of first multiple */
    llu global_count;         /* Global prime count */
    llu high_value;           /* Highest value on this proc */
    llu i;
    int id;    /* Process ID number */
    llu index; /* Index of current prime */
    llu base_index;
    struct aligned_llu b_first;
    struct aligned_llu b_last;
    struct aligned_llu b_high_value;
    struct aligned_llu b_low_value;
    struct aligned_llu prime;
    llu low_value; /* Lowest value on this proc */
    char *marked;  /* Portion of 2,...,'n' */
    char *base_marked;
    unsigned long long int n; /* Sieving from 2, ..., 'n' */
    unsigned long long int N;
    int p;          /* Number of processes */
    llu proc0_size; /* Size of proc 0's subarray */
    llu size;       /* Elements in 'marked' */
    llu base_size;
    llu base_low_value;
    llu base_high_value;
    llu count_cacheBlock;
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

    n = atoll(argv[1]);
    // printf("n= %d,N=%d\n",n,N);
    // n = N % 2 ? N - 1 : N;
    /* Figure out this process's share of the array, as
      well as the integers represented by the first and
      last array elements */
    base_low_value = 3;
    base_high_value = 1 + ((n - 1) / 2 / p) * 2;
    base_size = (int)sqrt((double)n);
    // printf("id: %d, low_value: %d, high_value: %d\n", id, low_value, high_value);
    // printf("id: %d, base_low_value: %d, base_high_value: %d\n", id, base_low_value, base_high_value);

    /* Bail out if all the primes used for sieving are
      not all held by process 0 */

    proc0_size = n - p;

    if ((2 + proc0_size) < (int)sqrt((double)n))
    {
        if (!id)
            printf("Too many processes\n");
        MPI_Finalize();
        exit(1);
    }

    /* Allocate this process's share of the array. */

    // marked = (char *)malloc(size);
    base_marked = (char *)malloc(base_size);
    if (base_marked == NULL)
    {
        printf("Cannot allocate enough memory\n");
        MPI_Finalize();
        exit(1);
    }

    for (i = 0; i < base_size; i++)
        base_marked[i] = 0;

    base_index = 0;
    prime.value = 3;
    N = (n - 1) >> 1;                                               
    low_value = (id * (N / p) + MIN(id, N % p)) * 2 + 3;            
    high_value = ((id + 1) * (N / p) + MIN(id + 1, N % p)) * 2 + 1; 
    size = ((high_value - low_value) >> 1) + 1;                     

    // cacheline数量
    int cacheline_count = l3CacheSize / coreCount / (cacheLine * p);
    // cache大小
    int cache_size = cacheline_count * 8;
    int block = cache_size - 1;
    //
    int line_need = size / cache_size;
    int line_rest = size % cache_size;
    int cache_index = 0;

    // allocate this process 's share of the array
    marked = (char *)malloc(cache_size);
    if (marked == NULL)
    {
        printf("Cannot allocate enough memory\n");
        MPI_Finalize();
        exit(1);
    }
    prime.value = 3;
    index = 0;
    while (cache_index <= line_need)
    {

        //cache更新；
        b_last.value = (cache_index + 1) * block + MIN(cache_index + 1, line_rest) - 1 + (id * (N / p) + MIN(id, N % p));
        b_first.value = cache_index * block + MIN(cache_index, line_rest) + (id * (N / p) + MIN(id, N % p));
        b_low_value.value = b_first.value * 2 + 3;
        if (cache_index == line_need)
        {
            b_high_value.value = high_value;
            b_last.value = (id + 1) * (N / p) + MIN(id + 1, N % p) - 1;
            cache_size = ((b_high_value.value - b_low_value.value) >> 1) + 1;
        }
        else
        {
            b_high_value.value = (b_last.value + 1) * 2 + 1;
        }

        index = 0;
        prime.value = 3;
        count_cacheBlock = 0;

        for (int i = 0; i < cache_size; i++)
            marked[i] = 0;

        // 块内素数
        do
        {
            if (prime.value * prime.value > b_low_value.value)
            {
                first.value = (prime.value * prime.value - b_low_value.value) >> 1;
            }
            else
            {
                if (!(b_low_value.value % prime.value))
                    first.value = 0;
                else if (b_low_value.value % prime.value % 2 == 0)
                    first.value = prime.value - ((b_low_value.value % prime.value) >> 1);
                else
                    first.value = (prime.value - (b_low_value.value % prime.value)) >> 1;
            }
            for (int i = first.value; i < cache_size; i += prime.value)
            {
                marked[i] = 1;
            }
            while (base_marked[++index])
                ;

            prime.value = index * 2 + 3;
        } while (prime.value * prime.value <= b_high_value.value);

        // 块内计数
        for (int i = 0; i < cache_size; i++)
        {
            if (marked[i] == 0)
            {
                count_cacheBlock++;
            }
        }

        count += count_cacheBlock;
        cache_index++;
    }
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
        printf("There are %llu primes less than or equal to %llu\n",
               global_count + 1, n);
        printf("SIEVE (%d) %10.6f\n", p, elapsed_time);
    }
    MPI_Finalize();
    return 0;
}