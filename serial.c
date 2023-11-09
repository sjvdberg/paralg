#include <stdio.h>
#include <stdbool.h>
#include <math.h>

void main()
{
    sieve(1000);
    sieve(10000);
    sieve(100000);
    sieve(1000000);
    sieve(10000000);
    sieve(100000000);
}

void sieve(long n)
{

    bool values[n]; 
    for(int i = 0; i < n; i++)
        values[i] = false;
    int strikecounter, ifcounter;
    strikecounter = 0; 
    ifcounter = 0;
    for(long long int j = 4; j < n; j+= 2)
        {
        strikecounter += 1;
        values[j] = true;
        }
    return;
    for(long long int i = 3; i < n; i++)
    {
        ifcounter += 1;
        if(values[i] != false)
            continue;
        ifcounter += 1;
        if(i*i >= n) 
            continue;
        for(long long int j = i*i; j < n; j+= 2*i)
            {
            strikecounter += 1;
            values[j] = true;
            }
    }
    printf("regular. n = %d.\n", n);
    printf("strike counter is %d \n", strikecounter);
    printf("if counter is %d \n", ifcounter);
    getchar();
}

