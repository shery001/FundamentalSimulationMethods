#include <iostream>
#include <fstream>
#include <cmath> // for power command
#include <vector>

using namespace std;

double KB = 1.38e-23;
double L0 = 1e-35;
double T0 = 2e4;
double ALP = 10.0;
double BET = -.5;
double TINIT = 1e7;
double NH = 1e6;

double LambdaT(double temp){
// Compute the Lambda function of this temperature
	double frac = temp/T0;
	if (temp <= T0)
		return L0 * pow(frac, ALP);
	if (temp > T0)
		return L0 * pow(frac, BET);
}

double Rhs(double temp){
// Compute the RHS of the equation
	return -2./(3.*KB)*NH*LambdaT(temp);
}

double do_rk_step(double temp, double delta_t){
// Do the second order Runge Kutta integration step
	double k1 = Rhs(temp);
	double k2 = Rhs(temp + k1*delta_t);
	return temp + 0.5*(k1+k2)*delta_t;

}

double GetErrorStepHalf(double temp, double time_step, double &t_a_in, double &t_b_in){
// This function calculates the difference between the result using
// the current step size and the result using half the step size
// calculate the first prediction using current step size:
	double t_a;
        t_a = do_rk_step(temp, time_step);
// now calculate the same doing two steps of half the size
        double t_b;
	t_b = do_rk_step(temp, time_step/2.);
	t_b = do_rk_step(t_b, time_step/2.);
        
        t_a_in = t_a;
        t_b_in = t_b;
// now compare the two results to estimate the error
	return abs(t_a - t_b);    
}

void do_rk_step_adaptive(double &temp_in, double &time_step_in, double local_error){
// This will perform an Runge Kutta step using step refinemend to obey the 
// maximum local error bound.
        double t_a, t_b;
        double error = GetErrorStepHalf(temp_in, time_step_in, t_a, t_b);

        double temp = temp_in;
        double time_step = time_step_in;
//cout<<error<<endl;

// now check if we meet the error bound!
	if (error < local_error){
            
		if (error * pow(2,3) < local_error){
// here the error is much smaller. We dobule the step size
			temp = do_rk_step(t_a, time_step);
			time_step = 2*time_step;
		}else{
// here the error is smaller but not much smaller, so everything is good
			temp = t_b;
		}

	}else{
// Here we need to reduce the step size! We will try this recursively...
                time_step = time_step/2.;
		do_rk_step_adaptive(temp, time_step, local_error);
// Now the right time step and temperature are filled in temp and time_step
	}
// the only way to go here is a successfull pass through the if =>  set variables
	temp_in = temp;
	time_step_in = time_step;
	return;
}

void ExportData(vector<double> *time, vector<double> *temp, vector<double> *steps, string path_data, string path_steps){
// Exports the data into a file for plotting
// the file's first line gives the number of data points
// then the lines alternate for time and temperature values
	ofstream file;
        ofstream file_steps;
        
	file.open(path_data.c_str());
        file_steps.open(path_steps.c_str());
        
	long n_pts = time -> size();
	cout << "Writing "<<n_pts<<" data points to file..."<<endl;
	if(file.is_open() && file_steps.is_open()){
		file << n_pts<<endl;
                file_steps << n_pts<<endl;
		for(int i = 0; i < n_pts; ++i){
			file << time->at(i)<<endl;
			file << temp -> at(i)<<endl;
                        
                        file_steps << time -> at(i)<<endl;
			file_steps << steps -> at(i)<<endl;
		}
		file.close();
                file_steps.close();
	}else{
		cout<<"could not write file!"<<endl;
	}

}

int main(){

double final_temp = 6e3;
double curr_temp = TINIT;
double curr_time = 0;
// maximum local error that should be obeyed by the step refinement
double error_max = 50.;
// initial step size (if algorithm performes well the result should be
// independent of this value...
double delta_t = 1E10;

// these two vectos will store the time and temperature afeter each step
vector<double> *time_steps = new vector<double>;
time_steps -> push_back(curr_time);
vector<double> *temp_steps = new vector<double>;
temp_steps -> push_back(curr_temp);
vector<double> *delta_t_steps = new vector<double>;
delta_t_steps -> push_back(delta_t);

	while(curr_temp >= final_temp){
// in this function delta_t and curr_temp will get modified!	
		do_rk_step_adaptive(curr_temp, delta_t, error_max);
		curr_time += delta_t;

// Store the new values in the vectors

		time_steps -> push_back(curr_time);
		temp_steps -> push_back(curr_temp);
                delta_t_steps -> push_back(delta_t);
	}
	string exp_path_temp = "data_ex1_c.dat";
        string exp_path_time_steps = "steps_ex1_c.dat";
	ExportData(time_steps, temp_steps, delta_t_steps, exp_path_temp, exp_path_time_steps);
}

