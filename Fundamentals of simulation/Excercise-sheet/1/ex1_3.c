/*pitfalls of floating point representation*/
/*a)*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
long main()
{
float x = 0.01;
double y = x;
double z = 0.01;

int i = x * 10000;
int j = y * 10000;
int k = z * 10000;

printf("%d %d %d\n", i, j, k);

/*the value of integer i, j is different even though y intially declared to be equal to x, the reason can be explained by seeing the different datatype of both variable. As the both data type would have been saved in sigule-precision IEE-754 floating point, the float data type only considers the 7 point after decimal point where in double it consider 15 points after decimal point, so the number saved in binary in IEEE-754 is the same but y being double has differnt interpretation of same number than of x, being float.*/
/*b)*/

/*float a = y;*/
int kk = 0;
float c = 0.0;
long long int dx = 10;
/*long double m = y*dx;*/
long int b;
for(kk=0;kk<10;kk++)
/*while(y*dx-b>c)*/
	{
	dx*=10;
	long double abc = y*dx;
        long int b = labs(abc);
	printf("%li\n", b);
	}
	printf("the last value is numerator");
printf("the denominator is = %lli\n", dx);
return 0;
}


