/*Pitfall of floating point arithmatic*/
/*calculate the result  for x and y*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
int main()
{
/* declaration of variable*/
double a = 1.0e17;
double b = -1.0e17;
double c = 1.0;

/*arithmatical process*/
double x = (a+b)+c;
double y = a+(b+c);	
printf("the value of x = %lf \n", x);
printf("the value of y = %lf \n", y);
return 0;	
}

/*by looking at the answer it's clear that x, in the calculation of y we deal with e17, where double only can save 15 decimal point number so any arithmatic operation happening below 15 decimal point is negelected for the chosen data type, b+c = b. and so the final answer in y will be zero, finally, which is wrong.*/
