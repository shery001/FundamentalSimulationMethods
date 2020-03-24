#include<stdio.h>
#include<conio.h>
#include<math.h>

void semincalculation(float x, float y){

	float a = (x + y**2)/y;
        return(a);

}

void main(){
	float x=3.0, y =5.0;
	float z;
	z = semincalculation(x,y);
	printf("the value of z is %f", & z);
}
