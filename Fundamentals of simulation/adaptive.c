#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>



double f(double t)
{
  return t*t;
}

/*Euler forward intergration function*/
/*...
 * ... to be filled in ...
 * ...*/

int main()
{
  double t = 0.0, x = 1.0;
  double dt = 0.5, tmax = 1.0;

  double max_error = 1.0e-5;
  double min_error = 0.1*max_error;
  double abserr, xlarge, xsmall;

  while(t < tmax)
    {
      xlarge = euler(x,t,dt);
      xsmall = euler(x,t,0.5*dt);
      xsmall = euler(xsmall,t+0.5*dt,0.5*dt);

      t += dt;
      x = xsmall;
      abserr = fabs(xlarge-xsmall);

      if(abserr/fabs(xsmall) < min_error)
        {
          dt *= 2.0;
        }
      if(abserr/fabs(xsmall) > max_error)
        {
          dt *= 0.5;
        }

      //Output
      /*...
       * ... to be filled in ...
       * ...*/

    }

  return 0;
}
