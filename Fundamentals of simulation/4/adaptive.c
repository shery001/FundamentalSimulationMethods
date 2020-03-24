#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

// define arrays for output
vector<double> time_st;
vector<double> x_st;
vector<double> delta_t_st;
vector<double> err_st;
vector<double> ana_st;


double f(double t)
{
  return t*t;
}

double euler(double x, double t, double delta_t)
{	
	return x + f(t)*delta_t;
}

void ExportData()
{
// Exports the data into a file for plotting
// the file's first line gives the number of data points
// then the lines alternate for time and temperature values
	ofstream file;
	file.open("data_ex4_1.dat");
	long n_pts = time_st.size();
	cout << "Writing "<<n_pts<<" data points to file..."<<endl;
	if(file.is_open()){
file<<"time, x, delta_t, error, analytical"<<endl;
		for(long i = 0; i < n_pts; ++i){
//&n_parts, &thetas, &n_nodes, &errors, &times_exact, &times_tree      
file<<time_st[i]<<", "<<x_st[i]<<", "<<delta_t_st[i]<<", "<<err_st[i]<<", "<<ana_st[i]<<endl;
		}
		file.close();
	}else{
		cout<<"could not write file!"<<endl;
	}
	file.close();
}

double get_analytical(double t)
{
    return  1./3. * t*t*t + 1;	
}

int main()
{
  double t = 0.0, x = 1.0;
  double dt = 0.5, tmax = 1.0;

  double max_error = 1.0e-5;
  double min_error = 0.1*max_error;
  double abserr, xlarge, xsmall;
  
    time_st.push_back(t);
    x_st.push_back(x);
    delta_t_st.push_back(dt);
    err_st.push_back(0);
    ana_st.push_back(get_analytical(t));  

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
	
	time_st.push_back(t);
	x_st.push_back(x);
	delta_t_st.push_back(dt);
	err_st.push_back(abserr);
	ana_st.push_back(get_analytical(t));

    }
	ExportData();
  return 0;
}
