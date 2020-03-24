#include<stdio.h>
#include<stdlib.h>
#include<math.h>

float semincalculation(float x, float y){

	float a = (x + y);
        return(a);

}

void main(){
	float x=3.0, y =5.0;
	float z = semincalculation(x,y);
	printf("the value of z is %f \n ", z);
	
}
