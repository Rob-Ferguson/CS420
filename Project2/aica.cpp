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
		double j1, j2, h, r1, r2;		//command line parameters for calculations
		int id;							//experiment ID number
		double N2;						//number of cells in the 30x30 grid (900)
		//string ouputname;				//output file
		vector <vector <int> > board;	//the grid of values (1 or -1)
		vector <double> sc;				//holds the calculated spatial correlation at each possible radius (0-14) by index
		vector <double> mi;				//holds the calculated mutual information at each possible radius (0-14) by index
		double H;						//holds the calculated entropy value of the entire system
		vector <double> Hj;				//holds the calculated joint entropy values at each possible radius (0-14) by index
		double lambda;					//the calculated correlation length

		void make_grid();				//randomly populate the grid with 1 or -1
		void update_grid();				//update the randomized grid until the aica until it has stabilized
		int cell_distance(int i1, int j1, int i2, int j2);
		void calc_sc();					//calculate spatial correlation at every possible distance
		void calc_lambda();				//calculate the characteristic correlation length (lambda)
		void calc_entropy();			//calculate the overall entropy of the system
		void calc_joint_entropy();		//calculate the joint entropy at every possible distance
		void calc_mutual_info();		//calculate the mutual information at every possible distance
		void record_data();
		void make_pgm_webpage(string filename);
};

//random number generator function for random_shuffle calls
int rand_num_gen(int i) { return rand() % i;}

int main(int argc, char *argv[]){ 
	AICA aica;			//instance of Activator/Inhibitor Cellular Automaton experiment class

	srand(unsigned (time(0)));
	//srand(1);
	
	//check that the number of arguments is correct and assign parameters
	if(argc == 7){
		aica.j1 = strtod(argv[1], NULL);
		aica.j2 = strtod(argv[2], NULL);
		aica.h  = strtod(argv[3], NULL);
		aica.r1 = strtod(argv[4], NULL);
		aica.r2 = strtod(argv[5], NULL);
		aica.id = strtod(argv[6], NULL);
		//aica.ouputname = argv[7];
	} else{
		fprintf(stderr, "Please enter doubles j1, j2, h, r1, r2 for calculations, and an experiment ID number.\n");
		fprintf(stderr, "USAGE: ./aica j1 j2 h r1 r2 id\n");
		exit(1);
	}
	
	//resize vectors holding spatial correlation, mutual information, and joint entropy
	aica.sc.resize(15);
	aica.mi.resize(15);
	aica.Hj.resize(15);

	//N2 is number of cells in grid (30^2)
	aica.N2 = pow(30, 2);
	
	//populat the board with either 1 or -1
	aica.make_grid();
	
	//update the AICA grid until the system stabilizes
	aica.update_grid();

	//perform relevant calculations on system
	aica.calc_sc();
	aica.calc_lambda();
	aica.calc_mutual_info();

	//save the data generated
	aica.record_data();

	return 0;
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
void AICA::update_grid(){
	vector <vector <bool> > cells_updated;		//true at (i, j) if the cell at (i, j) has been updated already
	vector <int> indices;						//all possible grid indices if grid were a 1D array
	vector <pair <int, int> > coords;			//all possible indices as coordinate pairs
	int x, y;
	int d;
	int new_state;
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

		/*	
		//loop through every cell in the grid
		for(int x1 = 0; x1 < 30; x1++){
			for(int y1 = 0; y1 < 30; y1++){
				//keep track of the sum of states for all cells in the grid
				cell_sum += board[x1][y1];
				
				//for each of these cells, iterate through all other cells
				for(int x2 = x1; x2 < 30; x2++){
					for(int y2 = y1 + 1; y2 < 30; y2++){
						//condition to avoid double-counting cells
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
			pl = abs(1 - pow(((1.0/N2) * cell_sum), 2));
		}else{
			pl = abs(((2.0/(N2 * cl)) * product_sum) - pow(((1.0/N2) * cell_sum), 2));
		}
		*/	


					
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
	string experiment = "Experiment_" + to_string(id);
	ofstream data;
	char outputline[1024];

	//if this is the first time the program has been run, make a new directory for files
	mkdir("./Experiment_Data", 0777);

	data.open("./Experiment_Data/" +  experiment + ".csv");
	if(!data.is_open()){
		fprintf(stderr, "Could not open current data vile (.csv)\n");
		exit(1);
	}

	//header info
	data << "Experiment# " << id << "\n";
	data << "Distance, Correlation, Entropy, Joint Entropy, Mutual Information, J1, J2, H, R1, R2\n";
	for(int i = 0; i < 15; i++){
		sprintf(outputline, "%d, %f, %f, %f, %f, %f, %f, %lf, %lf, %lf\n", i, sc[i], H, Hj[i], mi[i], j1, j2, h, r1, r2);
		data << outputline;
	}

	data.close();

	
}


//for creating a .pgm image of a given CA system
void AICA::make_pgm_webpage(string filename){

	// define image dimensions
	unsigned short width   =  480;
	unsigned short height  =  640;

	// allocate memory for your data
	unsigned char *buff = new unsigned char[width*height*sizeof(uchar)];

	for (int i = 0; i < width*height; i++)
		buff[i] = rand() % 256;


	// output file streams
	ofstream fout (filename.c_str());

	if (!fout.is_open())
	{
		cout << "Can't open output file"  << filename << endl;
		exit(1);
	}

	// write the header
	fout << "P5\n" << width << " " << height << " 255\n";

	// write the data
	fout.write((char *)buff, width*height*sizeof(uchar));

	// close the stream
	fout.close();

	// free memory
	delete[] buff;
}









