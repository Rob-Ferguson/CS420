#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdio>
#include <vector>
#include <time.h>
#include <utility>
#include <cmath>
#include <sys/stat.h>

using namespace std;

typedef unsigned char uchar;

#define E (2.7182818284590452353602874713526624977572470937L )

//Class for Activator/Inhibitor Cellular Automaton Experiment Information
class AICA{
	public:
		string name;					//for simulation identification purposes
		double j1, j2;					//activation and inhibition parameters
		int h, r1, r2;					//command line parameters for calculations
		int id;							//simulation ID number
		int experiment;					//Overall Experiment 1, 2, or 3 as defined in the lab writeup
		double N2;						//number of cells in the 30x30 grid (900)
		vector <vector <int> > board;	//the grid of values (1 or -1)
		vector <double> sc;				//holds the calculated spatial correlation at each possible radius (0-14) by index
		vector <double> mi;				//holds the calculated mutual information at each possible radius (0-14) by index
		double H;						//holds the calculated entropy value of the entire system
		vector <double> Hj;				//holds the calculated joint entropy values at each possible radius (0-14) by index
		int lambda;						//the calculated correlation length

		AICA(double J1, double J2, int h_val, int R1, int R2, int sim_num, int exp_num, int run);	//initialize new sim
		void make_grid();				//randomly populate the grid with 1 or -1
		bool update_grid();				//update the randomized grid until the aica until it has stabilized
		int cell_distance(int i1, int j1, int i2, int j2);
		bool do_calculations();			//function that calls other calculation functions
		void calc_sc();					//calculate spatial correlation at every possible distance
		void calc_lambda();				//calculate the characteristic correlation length (lambda)
		void calc_entropy();			//calculate the overall entropy of the system
		void calc_joint_entropy();		//calculate the joint entropy at every possible distance
		void calc_mutual_info();		//calculate the mutual information at every possible distance
		void record_data();
		void make_pgm();
};

//each set of parameters are run multiple times and the results are averaged and stored in this class
class avg_AICA{
	public:
		vector <AICA *> simulations;
	
		vector <double> avg_sc;				//holds the calculated spatial correlation at each possible radius (0-14) by index
		vector <double> avg_mi;				//holds the calculated mutual information at each possible radius (0-14) by index
		double avg_H;						//holds the calculated entropy value of the entire system
		vector <double> avg_Hj;				//holds the calculated joint entropy values at each possible radius (0-14) by index
		double avg_lambda;					//the calculated correlation length

		avg_AICA(int runs);					//ensure all of the containers are set up before data is collected from each simulation
};


//possible parameters to test for a specified experiment (1-3)
struct Params{
	double J1, J2;
	vector <int> r1_vals;
	vector <int> r2_vals;
	vector <int> h_vals;
};

//random number generator function for random_shuffle calls
int rand_num_gen(int i) { return rand() % i;}

//for determining the parameters of a given experiment for automatic simulation generation
Params *get_parameters(int exp_num);

int main(int argc, char *argv[]){ 
	AICA *aica;			//instance of Activator/Inhibitor Cellular Automaton simulation class
	avg_AICA *avg_sims;	//instance of class that holds the average of a set of AICA sims for one set of parameters
	Params *args;		//struct that holds the system parameters defined by a particular experiment (1, 2, or 3)

	int sim_num = 1;						//current simulation number (for labeling)
	int exp_num;
	int runs = 3;
	double j1;								//the activation modifier within system
	double j2;								//the inhibition modifier within system 
	bool considered = true;

	//initialize random seed
	srand(unsigned (time(0)));

	//if experiment number (1-3) is all that is given, run simulations with varying parameters
	if(argc == 2){
		exp_num = atoi(argv[1]);
		
		//determine the parameters to use based on the experiment number
		args = get_parameters(exp_num);
		
		//vary r1
		for(int r1 = 0; r1 < args->r1_vals.size(); r1++){
			//vary r2
			for(int r2 = 0; r2 < args->r2_vals.size(); r2++){
				//vary h
				for(int h = 0; h < args->h_vals.size(); h++){

					//with these r1, r2, and h vals, run multiple times and record the average results
					avg_sims = new avg_AICA(runs);

					//run each set of parameters 3 times to find average results with random initial states
					for(int i = 0; i < runs; i++){
						//initialize values for calculations and AICA structure
						aica = new AICA(args->J1, args->J2, args->h_vals[h], args->r1_vals[r1], args->r2_vals[r2], sim_num, exp_num, i);

						considered = aica->do_calculations();

						if(considered){
							printf("Simulation %d_%c has finished\n", sim_num, 'a' + i);
							avg_sims->simulations.push_back(aica);
						}
						else{
							printf("Simulation %d_%c was discarded\n", sim_num, 'a' + i);
							delete aica;
						}
					}
					
					//write average results to file


					//free up memory for simulation instances and the avg_sims instance itself
					for(int j = 0; j < avg_sims->simulations.size(); j++) delete avg_sims->simulations[j];
					delete avg_sims;

					//increment to the next unique simulation
					sim_num++;

				}
			}
		}
		delete args;
	}
	//if enough parameters are given, run one simulation with those values
	else if(argc == 7){
		j1 = strtod(argv[1], NULL);
		j2 = strtod(argv[2], NULL);
		int h = atoi(argv[3]);
		int r1 = atoi(argv[4]);
		int r2 = atoi(argv[5]);
		int sim_num = atoi(argv[6]);

		//initialize values for calculations and AICA structure
		aica = new AICA(j1, j2, h, r1, r2, sim_num, 0, 0);

		considered = aica->do_calculations();

		if(considered) printf("Simulation %d has finished\n", sim_num);
		else printf("Simulation %d was discarded\n", sim_num);

		delete aica;
	}
	//if neither of the allowable formats are given, throw an error
	else{
		fprintf(stderr, "Please provide either 6 arguments for manual simulation or 1 for automatic\n");
		fprintf(stderr, "USAGE: ./aica j1 j2 h r1 r2 experiment_number\n");
		fprintf(stderr, "USAGE: ./aica experiment_number (where experiment_number is 1, 2, or 3)\n");
		exit(1);
	}

	return 0;
}


//averaged results container constructor
avg_AICA::avg_AICA(int runs){
	//resize vectors holding spatial correlation, mutual information, and joint entropy
	avg_sc.resize(15);
	avg_mi.resize(15);
	avg_Hj.resize(15);
}

//system constructor
AICA::AICA(double J1, double J2, int h_val, int R1, int R2, int sim_num, int exp_num, int run){
	char x = 'a' + run;					//for a set of simulations with same params, distinguish the run
	
	//primary experiment parameters
	j1 = J1;
	j2 = J2;
	r1 = R1;
	r2 = R2;
	h = h_val;
	id = sim_num;
	experiment = exp_num;
	name = "Simulation_" + to_string(experiment) + "_" + to_string(id) + "_"; 
	name.push_back(x);

	//resize vectors holding spatial correlation, mutual information, and joint entropy
	sc.resize(15);
	mi.resize(15);
	Hj.resize(15);

	//N2 is number of cells in grid (30^2)
	N2 = pow(30, 2);	
}


//primary function that performs all calculations
bool AICA::do_calculations(){
	bool varied = true;

	//populat the board with either 1 or -1
	make_grid();

	//update the AICA grid until the system stabilizes
	varied = update_grid();
	
	//if the system is varied enough to be interesting, perform calculations
	if(varied){
		//perform relevant calculations on system
		calc_sc();
		calc_lambda();
		calc_mutual_info();

		//save the data generated and create pictures of the system
		record_data();
		make_pgm();
		
		//the system was processed, so return true
		return true;
	} 
	//if more than 95% of the state converges to a single state, discard the experiment
	return false;
}

//for randomly populating a 30x30 grid with either -1 or 1
void AICA::make_grid(){
	int rand_num;
	
	//resize the board to be 30x30
	board.resize(30);
	for(int i = 0; i < 30; i++) board[i].resize(30);
	
	//populate it with 1 or -1
	for(int x = 0; x < 30; x++){
		for(int y = 0; y < 30; y++){
			if((rand() % 2) == 1) board[x][y] = 1;
			else board[x][y] = -1;
		}
	}
}


//continuously update the grid until the aica has fully stabalized
bool AICA::update_grid(){
	vector <vector <bool> > cells_updated;		//true at (i, j) if the cell at (i, j) has been updated already
	vector <int> indices;						//all possible grid indices if grid were a 1D array
	vector <pair <int, int> > coords;			//all possible indices as coordinate pairs
	int x, y;
	int d;
	int new_state;
	int psum = 0, nsum = 0;								//sum of positive and negative states after stabilization
	double near_sum, far_sum;
	double result;
	bool stabilized = false;

	//resize cells_updated 2D array to match size of grid (30x30)
	cells_updated.resize(30);
	for(int i = 0; i < 30; i++) cells_updated[i].resize(30);

	//resize indices vector to represent the grid as a 1D array
	indices.resize(30*30);
	for(int i = 0; i < indices.size(); i++) indices[i] = i;
		
	//to randomly choose a cell, shuffle indices vector and convert each value to (x, y) coords into grid
	random_shuffle(indices.begin(), indices.end(), rand_num_gen);
	for(int i = 0; i < indices.size(); i++){
		x = indices[i] / 30;
		y = indices[i] % 30;
		coords.push_back(make_pair(x, y));
	}

	//continue until cells are no longer changing (i.e. the AICA has stabilized)
	while(!stabilized){
		//initialized stabilized bool to true until a cell changes states
		stabilized = true;

		//reshuffle coords vector to pick a difference sequence of cell positions to update
		random_shuffle(coords.begin(), coords.end(), rand_num_gen);

	
		/*
		cout << endl;
		cout << "STARTING UPDATE PROCESS" << endl << endl;
		for(int i = 0; i < board.size(); i++){
			for(int j = 0; j < board[0].size(); j++){
				if(board[i][j] == 1) cout << "*";
				else if(board[i][j] == -1) cout << " ";
			}
			cout << endl;
		}
		*/
		
		


		//asychronously update each cell by traversing randomized coords vector 
		for(int z = 0; z < coords.size(); z++){
			//reinitialize variables for a new cell to be updated
			near_sum = 0;
			far_sum = 0;
		
			//compare distance of current cell at coords[i] to every other cell in grid
			for(int i = 0; i < board.size(); i++){
				for(int j = 0; j < board[0].size(); j++){
					d =	cell_distance(coords[z].first, coords[z].second, i, j);
				
					//if distance is zero, it's the same cell, so don't add to either sum
					if(d != 0){
						//if d < R1 add to near cell sum
						if(d < r1) near_sum += board[i][j];
						//if R1 <= d < R2, add to far away cell sum
						else if(d >= r1 && d < r2) far_sum += board[i][j];
					}
				}
			}
		
			//multiply near sum by j1, far sum by j2, sum the results, and also add in the h bias
			near_sum = near_sum * j1;
			far_sum = far_sum * j2;
			result = (h + near_sum + far_sum);		

			//if the result is negative, the cell at coords[z] will be -1 at the next time step
			if(result < 0){
				//only change state if cell is currently in state 1
				if(board[coords[z].first][coords[z].second] == 1){
					board[coords[z].first][coords[z].second] = -1;
					//because a state changed at this time step, the system has not stabilized 
					stabilized = false;
				}
			}
			//if the result is positive, the cell at coords[z] will be 1 at the next time step
			else if(result > 0){
				//only change state if cell is currently in state -1
				if(board[coords[z].first][coords[z].second] == -1){
					board[coords[z].first][coords[z].second] = 1;
					//because a state changed at this time step, the system has not stabilized 
					stabilized = false;
				}
			}
			//in the rare case that the transition function result equals zero, cell keeps same state

		}	//end of coords loop
		
	}

	//after the system has stabilized, ensure it does not converge to a single state (uninteresting)
	for(int x = 0; x < 30; x++){
		for(int y = 0; y < 30; y++){
			if(board[x][y] == 1) psum++;
			else if(board[x][y] == -1) nsum++;
		}
	}

	//if more than 95% of the space converges to one state or the other, discard the experiment
	if((psum > (0.95 * N2)) || (nsum > (0.95 * N2))) return 0;

	return 1;
}


//for calculating the distance between any two cells in the 30x30 grid
int AICA::cell_distance(int x1, int y1, int x2, int y2){
	int x_dist = abs(x1 - x2);
	int y_dist = abs(y1 - y2);

	//ensure that the indices given are within the 30x30 grid
	if((x1 < 0 || x1 > 29) || (y1 < 0 || y1 > 29) || (x2 < 0 || x2 > 29) || (y2 < 0 || y2 > 29)){
		fprintf(stderr, "The indices provided in cell_distance() cannot be negative or greater than 29\n");
		exit(1);
	}

	/* The space is a torus, so the indices technically wrap around on the top/bottom and left/right sides */
	if(x_dist > 15) x_dist = 30 - x_dist;
	if(y_dist > 15) y_dist = 30 - y_dist;
	
	//final distance between two cells with i and j being row and column index in grid, respectively
	return x_dist + y_dist;	
}


//function for calculation the spacial correlation of an AICA system
void AICA::calc_sc(){
	int cell_dist = 0;
	double cl;			
	double pl;					//spatial correlation at a given distance (l)
	double product_sum = 0.0;
	double cell_sum = 0.0;

	//iterate through all possible distances between cells (0-14)
	for(int l = 0; l < 15; l++){
		cell_sum = 0;
		product_sum = 0;
		pl = 0;
		cl = 4 * l;
					
		//loop through every cell in the grid
		for(int x1 = 0; x1 < 30; x1++){
			for(int y1 = 0; y1 < 30; y1++){
				//keep track of the sum of states for all cells in the grid
				cell_sum += board[x1][y1];
				
				//for each of these cells, iterate through all other cells
				for(int x2 = 0; x2 < 30; x2++){
					for(int y2 = 0; y2 < 30; y2++){
						//if((x1 != x2) || ((x1 == x2) && (y2 > y1))){
							cell_dist = cell_distance(x1, y1, x2, y2);

							//if current cell is within distance being considered, multiply states and add to sum
							if(cell_dist == l) product_sum += (board[x1][y1] * board[x2][y2]);
						
						//}
					}
				}	//end of inner row loop
			}
		}	//end of outer row loop
		
		//calculate the spatial correlation (pl) at the current distance
		if(l == 0){
			pl = abs(1 - pow(((1/N2) * cell_sum), 2));
		}else{
			pl = abs(((1/(N2 * cl)) * product_sum) - pow(((1/N2) * cell_sum), 2));
		}

		//add pl to sc vector at index l
		sc[l] = pl;
	}	//end of distance loop (l)
}


//function for calculation the characteristic correlation length of an AICA system
void AICA::calc_lambda(){
	double pl = sc[1]/E;
	double val = 0.0;
	double closest_val = 100.0;
	int closest_index = 0;

	//lambda is the index of the calculated spacial correlation value closest to sc[0]/e
	for(int i = 0; i < sc.size(); i++){
		val = abs(pl - sc[i]);
		if(val < closest_val){
			closest_val = val;
			closest_index = i;
		}
	}
	lambda = closest_index;
}


//function for calculating the entropy of an AICA system
void AICA::calc_entropy(){
	int bsum = 0;			//sum of converted binary states of all cells in grid
	double pospr = 0.0;		//probability of state value 1
	double negpr = 0.0;		//probability of state value -1
	double pos_term = 0.0;	//the state value 1 portion of the entropy equation
	double neg_term = 0.0;	//the state value -1 portion of the entropy equation

	//run through every cell and convert value from bipolar (-1,1) to binary (0,1)
	for(int x = 0; x < 30; x++){
		for(int y = 0; y < 30; y++){
			bsum += (1 + board[x][y]) / 2;
		}
	}

	//calculate probability of state value 1
	pospr = (1/N2) * bsum;
	//calculate probability of state value -1
	negpr = 1 - pospr;

	//using these probabilities, calculate the overall entropy
	if(pospr != 0) pos_term = pospr * log2(pospr);		
	else pos_term = 0;
	
	if(negpr != 0) neg_term = negpr * log2(negpr);		
	else neg_term = 0;
	
	H = (-1.0) * (pos_term + neg_term);
}


//function for calculating the joint entropy of an AICA system
void AICA::calc_joint_entropy(){
	int pbin_state1, pbin_state2;	//holds the positive state of two cells converted from bipolar to binary format
	int nbin_state1, nbin_state2;	//holds the negative state of two cells converted from bipolar to binary format
	double cl;
	double cell_dist;
	double jpos_sum;				//cumulative sum of the joint probability of two cells being in state 1
	double jneg_sum;				//cumulative sum of the joint probability of two cells being in state -1
	double pos_prob;				//overall probability of two cells both being positive (1) state
	double neg_prob;				//overall probability of two cells both being negative (-1) state
	double neut_prob;				//overall probability of two cells having different states (neutral)
	double pos_term;				//the state value 1 portion of the joint entropy equation
	double neg_term;				//the state value -1 portion of the joint entropy equation
	double neut_term;				//the mixed state value portion of the joint entropy equation
	double Hl = 0.0;				//joint entropy of system at distance l (0-14)

	//can't consider distance of 0, so joint entropy with l = 0 is just zero
	Hj[0] = 0;

	//iterate through all possible distances between cells (1-14)
	for(int l = 1; l < 15; l++){
		cl = 4 * l;
		jpos_sum = 0;
		jneg_sum = 0;

		//loop through every cell in the grid
		for(int x1 = 0; x1 < 30; x1++){
			for(int y1 = 0; y1 < 30; y1++){
				//convert first cell state value to binary equivalent
				pbin_state1 = (1 + board[x1][y1]) / 2; 
				nbin_state1 = (1 - board[x1][y1]) / 2; 
				
				//cout << "pbin1: " << pbin_state1 << ", nbin1: " << nbin_state1 << endl;

				//for each of these cells, iterate through all other cells
				for(int x2 = 0; x2 < 30; x2++){
					for(int y2 = 0; y2 < 30; y2++){
						//convert second cell state value to binary equivalent
						pbin_state2 = (1 + board[x2][y2]) / 2; 
						nbin_state2 = (1 - board[x2][y2]) / 2; 

						//condition to avoid double-counting cells
						//if((x1 != x2) || ((x1 == x2) && (y2 > y1))){
							cell_dist = cell_distance(x1, y1, x2, y2);

							//if current cell is within distance being considered, calculate state sums
							if(cell_dist == l){
								
								//calculate the product of states and add to positive sum
								jpos_sum += (pbin_state1) * (pbin_state2); 

								//calculate the product of states and add to negative sum
								jneg_sum += (nbin_state1) * (nbin_state2);
							}
						//}
					}
				}	//end of inner row loop
			}
		}	//end of outer row loop
		
		//calculate probability of two arbitrary cells both being positive
		pos_prob = (1/(N2 * cl)) * jpos_sum;

		//calculate probability of two arbitrary cells both being negative
		neg_prob = (1/(N2 * cl)) * jneg_sum;

		//calculate probability of two arbitrary cells having different state values
		neut_prob = 1 - pos_prob - neg_prob;
		if(neut_prob < 0){
			fprintf(stderr, "Invalid (negative) neutral probability in joint entropy function\n");
			exit(1);
		}
		
		//use these probabilities to calculate the joint entropy of the system at distance l
		if(pos_prob != 0) pos_term = pos_prob * log2(pos_prob);
		else pos_term = 0;
		if(neg_prob != 0) neg_term = neg_prob * log2(neg_prob);
		else neg_term = 0;
		if(neut_prob != 0) neut_term = neut_prob * log2(neut_prob);
		else neut_term = 0;
		
		Hl = (-1.0) * (pos_term + neg_term + neut_term);
		
		//add result to joint entropy vector at position l
		Hj[l] = Hl;
	}
}


//for calculating the mutual information at every possible distance
void AICA::calc_mutual_info(){

	//first, calculate the entropy of the system
	calc_entropy();

	//second, calculate the joint entropy of cell pairs at every possible distance (0-14)
	calc_joint_entropy();

	for(int i = 0; i < mi.size(); i++){
		mi[i] = (2 * H) - Hj[i];
	}
}


//after all the calculations have been finished, write the data to a .csv file
void AICA::record_data(){
	string datafilename = "MasterData_" + to_string(experiment) + ".csv";
	ofstream data;
	char outputline[1024];

	//append new data to the MasterData file associated with the experiment number being run
	data.open(datafilename, ios_base::app);
	if(!data.is_open()){
		fprintf(stderr, "Could not open current data vile (.csv)\n");
		exit(1);
	}

	//header info
	data << "Simulation# " << id << "\n";
	data << "Distance, Correlation, Lambda, Entropy, Joint Entropy, Mutual Information, J1, J2, H, R1, R2\n";
	
	//print out the data associated with each column specified above
	for(int i = 0; i < 15; i++){
		sprintf(outputline, "%d, %f, %d, %f, %f, %f, %f, %f, %d, %d, %d\n", i, sc[i], lambda, H, Hj[i], mi[i], j1, j2, h, r1, r2);
		data << outputline;
	}
	data << "\n";

	data.close();
}


//for creating a .pgm image of a given CA system
void AICA::make_pgm(){
	ofstream fout;
	//string image = "Simulation_" + to_string(experiment) + "_" + to_string(id) + ".pgm";
	string image = name + ".pgm";
	string dir_str = "./Simulation_Images/Experiment" + to_string(experiment);

	//if this is the first time the program has been run, make a new directory for images
	mkdir("./Simulation_Images", 0777);
	mkdir(dir_str.c_str(), 0777);

	fout.open(dir_str + "/" + image);
	
	if (!fout.is_open()){
		cout << "Can't open image output file: " << image << endl;
		exit(1);
	}

	// write the .pgm header
	fout << "P2\n"; 
	fout << 30 << " " << 30 << " \n255\n";

	// write the data
	for(int x = 0; x < 30; x++){
		for(int y = 0; y < 30; y++){
			if(board[x][y] == 1) fout << 0;
			else if(board[x][y] == -1) fout << 255;
			if(y < 29) fout << " ";
		}
		fout << "\n";
	}
	
	fout.close();
}

//for determining the parameters of a given experiment for automatic simulation generation
Params *get_parameters(int exp_num){
	Params *vals = new Params;

	if(exp_num == 1){
		vals->J1 = 1;
		vals->J2 = 0;
		vals->r1_vals.push_back(1);
		vals->r1_vals.push_back(3);
		vals->r1_vals.push_back(6);
		vals->r1_vals.push_back(9);
		vals->r1_vals.push_back(12);
		vals->r2_vals.push_back(6);
		vals->r2_vals.push_back(7);
		vals->r2_vals.push_back(12);
		vals->r2_vals.push_back(13);
		vals->h_vals.push_back(0);
		vals->h_vals.push_back(-1);
		vals->h_vals.push_back(-2);
		vals->h_vals.push_back(-3);
		vals->h_vals.push_back(1);
		vals->h_vals.push_back(2);
		vals->h_vals.push_back(3);
	}
	else if(exp_num == 2){
		vals->J1 = 0;
		vals->J2 = -0.1;
		vals->r1_vals.push_back(1);
		vals->r1_vals.push_back(4);
		vals->r1_vals.push_back(9);
		vals->r2_vals.push_back(2);
		vals->r2_vals.push_back(4);
		vals->r2_vals.push_back(5);
		vals->r2_vals.push_back(6);
		vals->r2_vals.push_back(7);
		vals->r2_vals.push_back(12);
		vals->r2_vals.push_back(13);
		vals->h_vals.push_back(0);
		vals->h_vals.push_back(-1);
		vals->h_vals.push_back(-2);
		vals->h_vals.push_back(-3);
		vals->h_vals.push_back(-5);
		vals->h_vals.push_back(-6);
	}
	else if(exp_num == 3){
		vals->J1 = 1;
		vals->J2 = -0.1;
		vals->r1_vals.push_back(1);
		vals->r1_vals.push_back(3);
		vals->r1_vals.push_back(7);
		vals->r1_vals.push_back(12);
		vals->r2_vals.push_back(2);
		vals->r2_vals.push_back(5);
		vals->r2_vals.push_back(9);
		vals->r2_vals.push_back(14);
		vals->h_vals.push_back(0);
		vals->h_vals.push_back(-1);
		vals->h_vals.push_back(-2);
		vals->h_vals.push_back(-3);
		vals->h_vals.push_back(-4);
		vals->h_vals.push_back(-6);
	}
	else{
		fprintf(stderr, "Please enter a valid experiment number (1-3)\n");
		exit(1);
	}
	return vals;
}




