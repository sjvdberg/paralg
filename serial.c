#include <stdio.h>
#include <stdbool.h>
#include <math.h>


void main()
{

    long n = 1000000;

    bool values[1000000] = {false}; 
    int strikecounter, ifcounter;
    strikecounter = 0; 
    ifcounter = 0;
    for(long long int i = 2; i < n; i++)
    {
        ifcounter += 1;
        if(values[i] != false)
            continue;
        printf("%d \n", i);
        ifcounter += 1;
        if(i*i >= n) 
            continue;
        for(long long int j = i; j < n; j+=i)
            {
            strikecounter += 1;
            values[j] = true;
            }
    }
    printf("strike counter is %d \n", strikecounter);
    printf("if counter is %d \n", ifcounter);
    getchar();
}

