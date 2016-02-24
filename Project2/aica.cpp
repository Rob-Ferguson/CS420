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

using namespace std;

typedef unsigned char uchar;

/* Calculate spatial correlation, mutual information, and characteristic correlation length lambda */

//Class for Activator/Inhibitor Cellular Automaton Experiment Information
class AICA{
	public:
		double j1, j2, h, r1, r2;		//command line parameters for calculations
		int id;							//experiment ID number
		string ouputname;				//output file
		vector <vector <int> > board;	//the grid of values (1 or -1)
		vector <double> sc;				//holds the calculated spatial correlation at each possible radius (0-14) by index
		vector <double> mi;				//holds the calculated mutual information at each possible radius (0-14) by index
		vector <double> lambda;			//holds the calculated correlation length (lambda) at each possible radius (0-14) by index

		void make_grid();				//randomly populate the grid with 1 or -1
		void update_grid();				//update the randomized grid until the aica until it has stabilized
		int cell_distance(int i1, int j1, int i2, int j2);
		int calc_sc();
		void make_pgm_webpage(string filename);
};

//random number generator function for random_shuffle calls
int rand_num_gen(int i) { return rand() % i;}

int main(int argc, char *argv[]){ 
	AICA aica;			//instance of Activator/Inhibitor Cellular Automaton experiment class

	srand(unsigned (time(0)));
	
	//check that the number of arguments is correct and assign parameters
	if(argc == 8){
		aica.j1 = strtod(argv[1], NULL);
		aica.j2 = strtod(argv[2], NULL);
		aica.h  = strtod(argv[3], NULL);
		aica.r1 = strtod(argv[4], NULL);
		aica.r2 = strtod(argv[5], NULL);
		aica.id = strtod(argv[6], NULL);
		aica.ouputname = argv[7];
	} else{
		fprintf(stderr, "Please enter doubles j1, j2, h, r1, r2 for calculations, an experiment ID number, and an output file name.\n");
		fprintf(stderr, "USAGE: ./aica j1 j2 h r1 r2 id outputimagename\n");
		exit(1);
	}
	
	//resize vectors holding spatial correlation, mutual information, and correlation length (lambda)
	aica.sc.resize(15);
	aica.mi.resize(15);
	aica.lambda.resize(15);

	//make_pgm_webpage(outputname);


	//populat the board with either 1 or -1
	aica.make_grid();
	aica.update_grid();

	aica.calc_sc();

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
	int cells_left = 30*30;						//number of cells that are left to be updated
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
		//while(cells_left > 0){
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
int AICA::calc_sc(){
	int cell_dist = 0;
	double N = 30, cl;			
	double pl;					//spatial correlation at a given distance (l)
	double product_sum = 0.0;
	double cell_sum = 0.0;

	//iterate through all possible distances between cells (0-14)
	for(int l = 0; l < 15; l++){
		cell_sum = 0;
		pl = 0;
		cl = 4 * l;
		
		//loop through every cell in the grid
		for(int x1 = 0; x1 < 30; x1++){
			for(int y1 = 0; y1 < 30; y1++){
				//printf("(%d,%d)\n", x1, y1);

				cell_sum += board[x1][y1];
				
				//for each of these cells, iterate through all other distinct pairs of cells
				for(int x2 = 0; x2 < 30; x2++){
					for(int y2 = 0; y2 < 30; y2++){
						if((x1 != x2) || ((x1 == x2) && (y2 > y1))){
							//printf("  (%d,%d)\n", x2, y2);
							cell_dist = cell_distance(x1, y1, x2, y2);

							if(cell_dist == l) product_sum += (board[x1][y1] * board[x2][y2]);
						
						}
					}
				}	//end of inner row loop
			}
		}	//end of outer row loop
		
		//calculate the spatial correlation (pl) at the current distance
		if(l == 0){
			pl = 1 - pow(((1/pow(N, 2)) * cell_sum), 2);
		}else{
			pl = ((2/(pow(N, 2) * cl)) * product_sum) - pow(((1/pow(N, 2)) * cell_sum), 2);
		}
	
		//add pl to sc vector at index l
		sc[l] = pl;
	}	//end of distance loop (l)

	for(int a = 0; a < sc.size(); a++) cout << a << ": " << sc[a] << endl;
	return 0;
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









