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

void ExportData(vector<double> *time, vector<double> *temp){
// Exports the data into a file for plotting
// the file's first line gives the number of data points
// then the lines alternate for time and temperature values
	ofstream file;
	file.open("data_ex2_a.dat");
	long n_pts = time -> size();
	cout << "Writing "<<n_pts<<" data points to file..."<<endl;
	if(file.is_open()){
		file << n_pts<<endl;
		for(int i = 0; i < n_pts; ++i){
			file << time->at(i)<<endl;
			file << temp -> at(i)<<endl;
		}
		file.close();
	}else{
		cout<<"could not write file!"<<endl;
	}

}

int main(){

double final_temp = 6e3;
double curr_temp = TINIT;
double curr_time = 0;

double delta_t = 1E10;

// these two vectos will store the time and temperature afeter each step
vector<double> *time_steps = new vector<double>;
time_steps -> push_back(curr_time);
vector<double> *temp_steps = new vector<double>;
temp_steps -> push_back(curr_temp);

	while(curr_temp >= final_temp){
		
		curr_time += delta_t;
		curr_temp = do_rk_step(curr_temp, delta_t);
// Store the new values in the vectors

		time_steps -> push_back(curr_time);
		temp_steps -> push_back(curr_temp);
	}

	ExportData(time_steps, temp_steps);
}

