#include <iostream>
#include <fstream>
#include <cmath> // for power command
#include <vector>

using namespace std;

double PI = 3.14159265;
double L1 = 2.;
double L2 = 1.;
double M1 = 0.5;
double M2 = 1.;
double G = 1.;

double ToRad(double angle){
	return angle/360. * 2*PI;
}
double ToGrad(double angle){
	angle/(2*PI) * 360.;
}

// Functions to compute components of the RHS
double Phi1Dot(vector<double> *state){
	double l1cos = L1*cos(state->at(0)-state->at(1));
        double nummerator = state->at(2)*L2-l1cos*state->at(3);   
	double denominator = L2*((M1+M2)*L1*L1-M2*l1cos*l1cos);
	return nummerator/denominator;
}

double Phi2Dot(vector<double> *state){
	double m2l1l2cos = M2*L1*L2*cos(state->at(0)-state->at(1));
	double nummerator = (M1+M2)*L1*L1*state->at(3) - state->at(2)*m2l1l2cos;
	double denominator = M2*L2*(M1+M2)*L1*L1 - m2l1l2cos*m2l1l2cos;
	return nummerator / denominator;
}

double GetEnergy(vector<double>* state){
// This function computes the Energy of the system i.e. it evaluates the 
// Hamiltonian function for "state". E=T+V

        double phi_1 = state -> at(0);
        double phi_2 = state -> at(1);
        double phi_1_dot = Phi1Dot(state);
        double phi_2_dot = Phi2Dot(state);
        
        double e_kin_1 = M1*0.5*L1*L1*phi_1_dot*phi_1_dot;
        double e_kin_2_pure = M2*0.5*L2*L2*phi_2_dot*phi_2_dot;
        double e_kin_12 = 0.5*M2*(L1*L1*phi_1_dot*phi_1_dot + 2*L1*L2*phi_1_dot*phi_2_dot*cos(phi_1-phi_2));
        double e_kin_tot = e_kin_1 + e_kin_2_pure + e_kin_12;
        double e_pot_1 = M1*G*L1*(1-cos(phi_1));
        double e_pot_2 = M2*G*(L1*(1-cos(phi_1)) + L2*(1-cos(phi_2)));
        double e_pot_tot = e_pot_1 + e_pot_2;
        
        return e_kin_tot + e_pot_tot;    
}

vector<double> Rhs(vector<double>* state){
// Compute the RHS of the equation using 
// the last state vector y = (phi_1, phi_2, q_1, q_2)
	double f_1 = Phi1Dot(state);
	double f_2 = Phi2Dot(state);
	double f1f2m2l1l2sin = f_1*f_2*M2*L1*L2*sin(state->at(0)-state->at(1));
	
	vector<double> result;
	result.push_back(f_1);
	result.push_back(f_2);
	double rest_f3 = L1*G*sin(state->at(0))*(M1+M2);
	result.push_back(-f1f2m2l1l2sin - rest_f3);
	double rest_f4 = L2*G*sin(state ->at(1))*M2;
	result.push_back(f1f2m2l1l2sin - rest_f4);

	return result;
}

vector<double> do_rk2_step(vector<double> *curr_state, double delta_t){
// Do the second order Runge Kutta integration step for each element
// in the state vector
	int n_dim = curr_state -> size();
	vector<double> k1 = Rhs(curr_state);
// now add k1*\Delta t to the state vector
	vector<double> estimate;
	for(int i = 0; i < n_dim; ++i){
		estimate.push_back(curr_state->at(i) + k1.at(i)*delta_t);		
	}
// use the estimated value to compute f(y_est)	
	vector<double> k2 = Rhs(&estimate);
// now calculate the new value for y:
	vector<double> result;
	for(int i = 0; i < n_dim; ++i){
		double this_addition = 0.5*(k1.at(i)+k2.at(i))*delta_t;
		result.push_back(curr_state->at(i) + this_addition);
	}
	
	return result;	
}


void ExportData(vector<double> *time, vector<double> *phi_1, vector<double> *phi_2, vector<double>* energy_steps){
// Exports the data into a file for plotting
// the file's first line gives the number of data points
// Then the file contains time, phi_1, phi_2\\...
	ofstream file;
	file.open("data_ex3_c.dat");
	long n_pts = time -> size();
	cout << "Writing "<<n_pts<<" data points to file..."<<endl;
	if(file.is_open()){
//		file << n_pts<<endl;
		for(int i = 0; i < n_pts; ++i){
			file << time->at(i)<<", "<< phi_1 -> at(i)<<", "<< phi_2 -> at(i)<<", "<< energy_steps->at(i)<<endl;
		}
		file.close();
	}else{
		cout<<"could not write file!"<<endl;
	}
}

int main(){

// initial values for the angles
double phi_1_init = 50.;
double phi_2_init = -120.0;
// Initial values conjugate momenta
double q_1_init = 0.;
double q_2_init = 0.;

double t_start = 0.;
double t_end = 10.0;
int n_steps = 200;
double delta_t = (t_end-t_start)/(double) n_steps;

// define the vector holding the curent state of the system
vector<double> curr_state;
curr_state.push_back(ToRad(phi_1_init));
curr_state.push_back(ToRad(phi_2_init));
curr_state.push_back(q_1_init);
curr_state.push_back(q_2_init);

double energy_0 = GetEnergy(&curr_state);

// define the vectors holding the results vs. time
vector<double> time_steps;
time_steps.push_back(t_start);
vector<double> phi_1_steps;
phi_1_steps.push_back(phi_1_init);
vector<double> phi_2_steps;
phi_2_steps.push_back(phi_2_init);
vector<double> energy_steps;
energy_steps.push_back((GetEnergy(&curr_state)-energy_0)/energy_0);

for (int i = 0; i < n_steps; ++i){
cout<<"Step: "<<i<<endl;
	time_steps.push_back(time_steps.at(i)+delta_t);
cout<<"time = "<<time_steps.at(i+1)<<endl;
	curr_state = do_rk2_step(&curr_state, delta_t);
	phi_1_steps.push_back(curr_state.at(0));
	phi_2_steps.push_back(curr_state.at(1));
        energy_steps.push_back((GetEnergy(&curr_state)-energy_0)/energy_0);
cout<<"phi1 = "<<phi_1_steps.at(i+1)<<" phi2 = "<<phi_2_steps.at(i+1)<<endl;
	
}

ExportData(&time_steps, &phi_1_steps, &phi_2_steps, &energy_steps);

}

