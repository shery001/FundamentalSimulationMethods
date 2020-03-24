#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(){
/*declaration of variable*/
float pi1[2000], pi2[2000],f1[2000],f2[2000],f3[2000], f4[2000];
float m1 = 0.5, m2 = 1.0;
float l1 = 1.0 , l2= 2.0;
float g = 1.0;
float dt = 0.05;
float k1, k2, j1, j2;
int i=0;
f1[0]=0;
f2[0]=0;
pi1[0]=50;
pi2[0]=-120;
for(i==0;i<=2000;i++){

	f3[i] = (-m2*l1*l2*f1[i]*f2[i]*sin(pi1[i]-pi2[i]))-(l1*g*sin(pi1[i])*(m1+m2));
	f4[i] = m2*l1*l2*f1[i]*f2[i]*sin(pi1[i]-pi2[i])-l2*g*sin(pi2[i])*m2;
	k1 = f1[i];
	k2 = pi1[i] + k1*dt;
	pi1[i+1] = pi1[i]+((k1+k2)*dt)/(2) ; 
	j1 = f2[i];
	j2 = pi2[i] + j1*dt;
	pi2[i+1] = pi2[i]+((j1+j2)*dt)/(2) ;
	i=i+1;
	f1[i] = f3[i-1] - (l1*cos(pi1[i]-pi2[i])*f4[i-2])/(l2*(m1*l1*l1+m2*l1*l1-l1*l1*m2*cos(pi1[i]-pi2[i])*cos(pi1[i]-pi2[i])));
	f2[i] = (f4[i-1]*(m1*l1*l1+m2*l1*l1-f3[i-1]*m1 - 2*l1*l2*cos(pi1[i]-pi2[i])))/(m2*l2*(m1*l1*l1*l2*l2*cos(pi1[i]-pi2[i])));
}

}
