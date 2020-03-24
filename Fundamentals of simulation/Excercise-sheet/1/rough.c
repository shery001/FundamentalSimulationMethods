#include<stdio.h>
#include<stdlib.h>
#include<math.h>
int main(){
long long int b = 10000000000000000000000000000;
float x = 0.01;
double y = x;
double z = y*b;
long int a = labs (z);
printf("%f\n",x);
printf("%.27g\n",y);
printf("%.27g\n",z);
printf("%li\n",a);
}


