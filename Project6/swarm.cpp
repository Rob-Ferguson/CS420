#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstdio>
#include <vector>
#include <map>
#include <time.h>
#include <utility>
#include <cmath>
#include <random>

using namespace std;

const int WORLDWIDTH = 100; 			//the width of the world space [-50, 50]
const int WORLDHEIGHT = 100; 			//the height of the world space [-50, 50]
const int NUMEPOCHS = 1000; 			//the max number of iterations to consider if x and y error hasn't dropped below errorThreshold 
const double ERRORTHRESHOLD = 0.001; 	//the target error for both x and y coordinates in a given simulation
const double EPSILON = 0.001; 			//required distance of particle from global max to be considered "close"

//global random number engine and random_device for seeding
random_device rd;
mt19937 globalRNG;

//helper function prototypes
void runBatchSims(int problem); 							//runs simulations for all combinations of parameters in batch_parameters.txt 
string newDataFile(const string &currentTest, int problem); //for initializing data files for certain experiments
unsigned int randomize();									//randomize the generator seed (returns that seed)
char generateBitChar();										//generate random char of either '0' or '1'
bool generateBool(); 										//generate a random boolean value
double generateDouble(double lower, double upper);			//generate random double in range [lower, upper)
double generateIndex(int Size); 							//generate a random index into a container of length Size


//parameters from which to construct the current set of (identical) genetic simulation runs
typedef struct {
    unsigned int seed; 				//a specified seed value to use for RNG (seed 0 means generate a new seed)
	int problem; 					//the problem number to consider for this simulation
	int numParticles; 				//the number of different particles to include in the simulation 
	double inertia; 				//influences how particle velocities change over time
	double cognition; 				//guides particles toward their own local best values
	double social; 					//guides particles toward the global best values 
	double maxVelocity; 			//limits how much a particle can move per epoch (greater values will be scaled down) 
	string filename; 	            //stores the name of the csv file currently being written to (the current test)
} SimParams;

//calculated statistics for a single simulation run
typedef struct {
	//indexed by epoch
	vector <double> xErrors; 			//the x error of particle positions at that epoch
	vector <double> yErrors; 			//the y error of particle positions at that epoch
	vector <double> fractionNearMax; 	//the fraction of particles within EPSILON of global max 
	vector <double> avgPopFitness; 		//the average fitness of the population of particles (according to the problem being solved)
	vector <double> avgDistance; 	 	//the average Euclidean distance of all particles from the global maximum at each epoch
} Statistics;

//a given particle that is acting within the simulation
class Particle{
	public:
		double fitness; 				//this particle's current fitness value
		double xPos; 					//current x coordinate position
		double yPos; 					//current y coordinate position
		double xVelocity; 				//current x component of the velocity
		double yVelocity; 				//current y component of the velocity
		double bestFitness; 			//this particle's best fitness value thus far
		double xBest; 					//this particle's best x coordinate position (highest personal fitness)
		double yBest; 					//this particle's best y coordinate position (highest personal fitness)
		double distance; 				//current Euclidean distance from this particle to the global maximum

		Particle();
};

//the actual simulation constructed from one set of invariant parameters
class SwarmSim {
	public:
		unsigned int seed; 				//the seed used to generate random variables for this simulation
		int problem; 					//the lab writeup problem number (1 or 2) that is being solved with this simulation
		int numParticles; 				//the number of different particles to include in the simulation 
		double inertia; 				//influences how particle velocities change over time
		double cognition; 				//guides particles toward their own local best values
		double social; 					//guides particles toward the global best values 
		double maxDist; 				//the maximum distance across the world space (for quality functions)
		double maxVelocity; 			//limits how much a particle can move in a given epoch (greater values are scaled down)
		int convergenceEpoch; 			//the epoch at which this simulation's particle position errors converge to a value less than the error threshold
		bool convergedToLocal; 			//true if this simulation actually converged to the local max, rather than global max (Problem 2)
		string filename;                //stores the name of the csv file currently being written to (the current test)
		Particle *best; 				//handle to the current best particle (highest fitness)
		vector <Particle> allParticles; //stores all of the particles in this simulation
		Statistics stats; 				//holds data about the performance of this simulation

		SwarmSim(SimParams *params);
		int runSim(); 						//updates the simulation (returns true if x/y error < error threshold)	
		void calcVelocity(Particle &p); 	//calculate the new velocity of a particular particle
		void moveParticle(Particle &p); 	//move the current particle based on its current position and velocity
		double calcXError(); 				//calculates the total error in x position values at a given epoch
		double calcYError(); 				//calculates the total error in y position values at a given epoch
		double calcCloseParticles(); 		//calculates the fraction of particles that are within EPSILON distance of global maximum per epoch
		double posDist(const Particle &p);
		double negDist(const Particle &p);
		double Q1(const Particle &p); 		//quality function for Problem 1
		double Q2(const Particle &p); 		//quality function for Problem 2

		void recordResults();
};

int main(int argc, char **argv){
	SwarmSim *simulation;
	unsigned int seed = 0;
	int lineNum = 0;				//the current line number in parameters.txt
	int problem = 0; 				//will be set to Problem 1 or Problem 2 depending on which problem the simulations are solving
	string inputLine;				//a given line of input (containing network parameters) from parameters.txt
	string currentTest = "Generic";	//current experiment (i.e., varying mutation probability = "VariedMutationProb")
	string filename = "";

	//if a single argument is given (problem number), run a batch test on all of the combinations of all parameters listed batch_parameters.txt
	if(argc == 2){
		//determine the problem number
		problem = atoi(argv[1]);
		if(problem != 1 && problem != 2){
			fprintf(stderr, "A valid problem number was not specified.\n");
			fprintf(stderr, "To run simulations with all combinations of parameters specified in batch_parameters.txt,\n");
			fprintf(stderr, "  you must use the format: ./swarm problem\n");
			exit(1);
		}
		
		printf("Attempting to run a batch of simulations using parameters from batch_parameters.txt...\n");

		runBatchSims(problem);
	}
	
	//if 7 or 8 arguments are given, run one simulation with the given parameters 
	else if(argc == 8 || argc == 9){
		SimParams p;				//holds all of the parameter information
		SimParams *params = &p;		//points to that parameter struct for passing into SwarmSim constructor

		//grab the params
		p.problem = atoi(argv[1]);
		p.numParticles = atoi(argv[2]);
		p.inertia = atof(argv[3]);
		p.cognition = atof(argv[4]);
		p.social = atof(argv[5]);
		p.maxVelocity = atof(argv[6]);
		currentTest = argv[7];

		//set the filename that this data will be recorded to
		if(p.problem == 1) p.filename = "Results/Problem1/" + currentTest + ".csv";
		else if(p.problem == 2) p.filename = "Results/Problem2/" + currentTest + ".csv";
		else{
			fprintf(stderr, "A valid problem number was not specified. The problem number must be 1 or 2\n");
			exit(1);
		}

		//grab an actual seed value if it is given (otherwise, seed is 0, which will signal to generate new seed)
		if(argc == 9){
			stringstream buf(argv[8]);
			if(!(buf >> seed)) { fprintf(stderr, "Could not convert given seed to an unsigned int\n"); exit(1); }
		}
		p.seed = seed;

		printf("Attempting to run a single simulation with the specified parameters...\n");

		simulation = new SwarmSim(params);
		
		//run the simulation
		simulation->convergenceEpoch = simulation->runSim();
		if(simulation->convergenceEpoch < NUMEPOCHS) printf("The simulation successfully converged by %d epochs\n", simulation->convergenceEpoch);
		else printf("The simulation did not converge before %d epochs and was halted\n", NUMEPOCHS);

		//save the data for this simulation to a .csv file
		simulation->recordResults();

		//free allocated memory
		delete simulation;
	}

	//if no arguments are given, just accept the parameters from stdin (which will be piped in from parameters.txt) 
	else if(argc == 1){
		printf("Grabbing parameter values from parameters.txt...\n");

		//parse the input parameters line by line from stdin (piped from parameters.txt) for each new simulation's parameters
		while(getline(cin, inputLine)){
			lineNum++;

			istringstream parser(inputLine);

			//skip any blank lines
			if(inputLine.length() == 0 || inputLine.find_first_not_of(" \n\t") == string::npos){
				continue;
			}

			//skip lines that are "commented out"
			if(inputLine[0] == '/'){
				continue;
			}

			//grab lines that state which problem to use for subsequent simulations (Problem 1 or Problem 2)
			if(inputLine[0] == 'P' || inputLine[0] == 'p'){
				//Problem 1
				if(inputLine.find('1') != string::npos) problem = 1;
				else if (inputLine.find('2') != string::npos) problem = 2;
				else { 
					fprintf(stderr, "Line %d of parameters.txt started with character '%c' but did not list a valid problem number (1 or 2)\n", lineNum, inputLine[0]);
					fprintf(stderr, "The correct way to specify a problem number in parameters.txt is to set a line to say Problem X, where X is a valid problem number\n");
					fprintf(stderr, "After specifying the problem number in this way, all subsequent simulations will use that problem until a new problem number is specified\n");
					exit(1);
				}

				printf("\nProblem Number %d\n", problem);
				continue;
			}


			//grab strings representing a new experiment description/name, which are prepended by '*'
			if(inputLine[0] == '*'){
				currentTest = inputLine;
				currentTest.erase(currentTest.begin());
				printf("\n  Experiment: %s\n", currentTest.c_str());

				//create a new .csv file for this experiment
				filename = newDataFile(currentTest, problem);

				continue;
			}

			SimParams p;				//holds all of the parameter information
			SimParams *params = &p;		//points to that parameter struct for passing into GeneticSim constructor

			//otherwise, try to individually parse the arguments on each line
			if(parser >> p.numParticles >> p.inertia >> p.cognition >> p.social >> p.maxVelocity){
				if(currentTest == "Generic") {
					fprintf(stderr, "The test description was not specified in parameters.txt\n");
					fprintf(stderr, "Subsequent simulation data will be placed in Generic.csv until an experiment name is given\n");
				}
				
				//set the problem number
				if(problem != 0){
					p.problem = problem;
				}else{
					fprintf(stderr, "The problem number was not initialized in parameters.txt\n");
					exit(1);
				}
				
				//grab the current filename that this simulation's resulting data will be stored in
				p.filename = filename;
				
				//seed is 0 for batch tests (so constructor will create/record a new random seed for each simulation)
				p.seed = 0;

			}else{
				fprintf(stderr, "Incorrect arguments in parameters.txt on line %d\n", lineNum);
				fprintf(stderr, "Parameters.txt USAGE: numParticles inertia cognition social maxVelocity\n");
				exit(1);
			}

			//create a new instance of the genetic simulation for the current set of parameters
			simulation = new SwarmSim(params);

			//run the simulation
			simulation->convergenceEpoch = simulation->runSim();
			if(simulation->convergenceEpoch < NUMEPOCHS) printf("     The simulation successfully converged by %d epochs\n", simulation->convergenceEpoch);
			else printf("     The simulation did not converge before %d epochs and was halted\n", NUMEPOCHS);

			//save the data for this simulation to a .csv file
			simulation->recordResults();

			//free allocated memory
			delete simulation;
		}
	}
	else { 
		fprintf(stderr, "Please enter the correct arguments or pipe in parameters.txt in the correct format\n");
		fprintf(stderr, "USAGE: ./swarm problemNum numParticles inertia cognition social maxVelocity currentTestName (optional: seed)\n");
		fprintf(stderr, " 		If given, the seed value must be greater than 0\n");
		exit(1); 
	}

	return 0;
}


//constructor for a particle swarm optimization simulator
SwarmSim::SwarmSim(SimParams *params){
	//seed the random number generator with the given seed or a new one
	if(params->seed == 0) seed = randomize();
	else{
		seed = params->seed;
		globalRNG.seed(seed);
	}


	//set user-specified parameters for the current simulation
	problem = params->problem;
	numParticles = params->numParticles;
	inertia = params->inertia;
	cognition = params->cognition;
	social = params->social;
	maxVelocity = params->maxVelocity;
	filename = params->filename;

	//calculate the maximum distance across the world space
	maxDist = sqrt(pow(WORLDWIDTH, 2) + pow(WORLDHEIGHT, 2)) / 2.0;

	/*
	printf("  Seed: %u\n", seed);
	printf("  Problem: %d\n", problem);
	printf("  Particles: %d\n", numParticles);
	printf("  Inertia: %lg\n", inertia);
	printf("  Cognition: %lg\n", cognition);
	printf("  Social: %lg\n", social);
	printf("  MaxDist: %lg\n", maxDist);
	printf("  MaxVelocity: %lg\n", maxVelocity);
	printf("  filename: %s\n", filename.c_str());
	*/

	//fill the container of particles with new Particle instances
	allParticles.reserve(numParticles);
	for(int i = 0; i < numParticles; i++) allParticles.push_back(Particle());

	//initialize the best particle to simply be first particle
	best = &allParticles[0];

	//use the current problem's fitness function to calculate the initial fitness values for each particles
	for(auto &particle: allParticles){
		if(problem == 1) particle.fitness = Q1(particle); 
		else if(problem == 2) particle.fitness = Q2(particle);
		else { fprintf(stderr, "The problem number of %d is invalid. Must be 1 or 2\n", problem); exit(1); }
		
		//also initialize this particle's current best to be the fitness value that was just calculated
		particle.bestFitness = particle.fitness;

		//if a better particle is found point to it instead
		if(particle.bestFitness > best->bestFitness){
			best = &particle;
		}
	}

	//initialize the convergedToLocal bool to be false until proven otherwise
	convergedToLocal = false;
}


//runs the simulation until the max number of epochs have been reached or the x and y position errors are below the threshold
int SwarmSim::runSim(){
	int epoch = 0;
	double xError;
	double yError;
	double closeParticles;
	double avgFitness;
	double avgDistanceFromMax;

	//continuously update the position of all particles in the world space
	do{
		avgFitness = 0.0;
		avgDistanceFromMax = 0.0;

		//periodically increase the capacity of data containers depending on how many epochs occur
		if(epoch % 100 == 0){
			stats.xErrors.reserve(epoch + 100);
			stats.yErrors.reserve(epoch + 100);
			stats.fractionNearMax.reserve(epoch + 100);
			stats.avgPopFitness.reserve(epoch + 100);
	 		stats.avgDistance.reserve(epoch + 100);
		}
		
		//calculate the total error in the x and y position values compared to the current global best position
		xError = calcXError();
		yError = calcYError();
		stats.xErrors.push_back(xError);
		stats.yErrors.push_back(yError);

		//calculate the percentage of particles within a certain distance of the global maximum (problem "solution")
		closeParticles = calcCloseParticles();
		stats.fractionNearMax.push_back(closeParticles);

		//update the velocity and position of each particle for each epoch
		for(auto &particle: allParticles){
			//add the current particle's fitness to a cumulative sum before recalculating it
			avgFitness += particle.fitness;

			//add the current particle's distance from the global maximum to a cumulative sum
			avgDistanceFromMax += particle.distance;

			calcVelocity(particle);
			moveParticle(particle);
			
			//recalculate the current particle's fitness value
			if(problem == 1) particle.fitness = Q1(particle); 
			else particle.fitness = Q2(particle);

			//update the current particle's personal best values if it has a higher fitness value
			if(particle.fitness > particle.bestFitness){
				particle.bestFitness = particle.fitness;
				particle.xBest = particle.xPos;
				particle.yBest = particle.yPos;
			}

			//update the pointer to the global best particle if a new global best has been found
			if(particle.bestFitness > best->bestFitness){
				best = &particle;
			}
		}

		//record the average population fitness (before fitness values were updated this iteration)
		avgFitness /= (double)allParticles.size();
		stats.avgPopFitness.push_back(avgFitness);	

		//record the average population fitness (before fitness values were updated this iteration)
		avgDistanceFromMax /= (double)allParticles.size();
		stats.avgDistance.push_back(avgDistanceFromMax);	

		epoch++;
	} while(epoch < NUMEPOCHS && (xError > ERRORTHRESHOLD || yError > ERRORTHRESHOLD));

	//if the average distance to the global max is not very small after convergence, simulation converged to local max (or didn't converge)
	if(stats.avgDistance[stats.avgDistance.size() - 1] > 1.0) convergedToLocal = true; 

	return epoch;
}


void SwarmSim::calcVelocity(Particle &p){
	double r1 = generateDouble(0, 1);
	double r2 = generateDouble(0, 1);
	double xNewVelocity;
	double yNewVelocity;

	//calculate the x-component of the velocity
	p.xVelocity = (inertia * p.xVelocity) + (cognition * r1 * (p.xBest - p.xPos))
		+ (social * r2 * (best->xBest - p.xPos));
	xNewVelocity = p.xVelocity;

	//calculate the y-component of the velocity
	p.yVelocity = (inertia * p.yVelocity) + (cognition * r1 * (p.yBest - p.yPos))
		+ (social * r2 * (best->yBest - p.yPos));
	yNewVelocity = p.yVelocity;

	//if the new velocity values are too large, scale them down proportionally with respect to the max velocity
	if(sqrt(pow(p.xVelocity, 2) + pow(p.yVelocity, 2)) > pow(maxVelocity, 2)){
		p.xVelocity = (maxVelocity / sqrt(pow(xNewVelocity, 2) + pow(yNewVelocity, 2))) * xNewVelocity;
		p.yVelocity = (maxVelocity / sqrt(pow(xNewVelocity, 2) + pow(yNewVelocity, 2))) * yNewVelocity;
	}
}


//update a particles position based on its current position and its new velocity
void SwarmSim::moveParticle(Particle &p){
	p.xPos += p.xVelocity;
	p.yPos += p.yVelocity;

	//ensure that the position values will wrap around the space if they go out of the world space boundaries
	if(p.xPos < 0 - (WORLDWIDTH / 2.0)) p.xPos += WORLDWIDTH;
	if(p.yPos < 0 - (WORLDHEIGHT / 2.0)) p.yPos += WORLDHEIGHT;
	if(p.xPos > (WORLDWIDTH / 2.0)) p.xPos -= WORLDWIDTH;
	if(p.yPos > (WORLDHEIGHT / 2.0)) p.yPos -= WORLDHEIGHT;
}


//calculates the total x position error for all particles at the current epoch
double SwarmSim::calcXError(){
	double error = 0.0;

	//calculate the error for each particle and add to a cumulative sum
	for(auto &particle: allParticles){
		error += pow(particle.xPos - best->xBest, 2);
	}

	error = sqrt((1.0 / (2.0 * numParticles)) * error);

	return error;
}


//calculates the total y position error for all particles at the current epoch
double SwarmSim::calcYError(){
	double error = 0.0;

	//calculate the error for each particle and add to a cumulative sum
	for(auto &particle: allParticles){
		error += pow(particle.yPos - best->yBest, 2);
	}

	error = sqrt((1.0 / (2.0 * numParticles)) * error);
	
	return error;
}


//calculates the percentage of particles within EPSILON Euclidean distance of the global maximum 
double SwarmSim::calcCloseParticles(){
	int xMax = 20; 		//both problems have a global max at (20, 7)
	int yMax = 7; 		//both problems have a global max at (20, 7)
	int numClose = 0; 		

	//iterate through all of the particles and calculate each's current distance to global max 
	for(auto &particle: allParticles){
		particle.distance = sqrt(pow(particle.xPos - xMax, 2) + pow(particle.yPos - yMax, 2));
		//keep track of how many are within the EPSILON distance range
		if(particle.distance < EPSILON) numClose++;
	}
	
	return ((double)numClose / (double)allParticles.size());
}


//calculates the positive distance of a particular particle
double SwarmSim::posDist(const Particle &p){
	return sqrt(pow(p.xPos - 20, 2) + pow(p.yPos - 7, 2));
}


//calculates the negative distance of a particular particle
double SwarmSim::negDist(const Particle &p){
	return sqrt(pow(p.xPos + 20, 2) + pow(p.yPos + 7, 2));
}


//the quality function that rates a particle for Problem 1
double SwarmSim::Q1(const Particle &p){
	double pdist = posDist(p);
	double ndist = negDist(p);
	
	return 100 * (1 - (pdist / maxDist));
}


//the quality function that rates a particle for Problem 2
double SwarmSim::Q2(const Particle &p){
	double pdist = posDist(p);
	double ndist = negDist(p);
	double val = 10 - pow(pdist, 2); 
	if(val < 0) val = 0.0;

	return (9 * val) + (10 * (1 - (pdist / maxDist))) + (70 * (1 - ndist / maxDist));
}


//write data for the current simulation to the corresponding Experiment set data file
void SwarmSim::recordResults(){
	FILE *output;
	bool converged = false;

	//if the simulation did not converge, note that, but don't print its related data
	if(convergenceEpoch < NUMEPOCHS) converged = true;
	
	//append to the corresponding experiment's data file
	output = fopen(filename.c_str(), "a");

	//record data for the current simulation within this experimental set 
	if(output != NULL){
		fprintf(output, "SIMULATION INFORMATION\nSeed: %u\n", seed);
		
		//if it didn't converge, just print the parameter info and ending particle coordinates for this simulation
		if(!converged){
			fprintf(output, "Convergence Epoch:, Simulation did not converge before %d epochs\n", NUMEPOCHS);
			fprintf(output, "Particles:, %d\nInertia:, %lf\nCognition:, %lf\nSocial:, %lf\nMaximum Velocity:, %lf\n",
					numParticles, inertia, cognition, social, maxVelocity);
			
			//print particle position coordinates at time when simulation halted
			fprintf(output, "\nParticle Positions at epoch %d\nParticle Number:, ", NUMEPOCHS);
			for(int i = 0; i < allParticles.size(); i++) fprintf(output, "%d, ", i);
			fprintf(output, "\nX Coordinate:, ");
			for(auto &particle: allParticles) fprintf(output, "%lf, ", particle.xPos);	
			fprintf(output, "\nY Coordinate:, ");
			for(auto &particle: allParticles) fprintf(output, "%lf, ", particle.yPos);	

			fprintf(output, "\n\n");
			fclose(output);
			return;
		}
		
		//otherwise, also print out the collected data for this simulation
		fprintf(output, "Convergence Epoch:, %d\nParticles:, %d\nInertia:, %lf\nCognition:, %lf\nSocial:, %lf\nMaximum Velocity:, %lf\n",
				convergenceEpoch, numParticles, inertia, cognition, social, maxVelocity);

		//print final particle position coordinates upon simulation convergence
		fprintf(output, "\nParticle Positions upon Convergence\nParticle Number:, ");
		for(int i = 0; i < allParticles.size(); i++) fprintf(output, "%d, ", i);
		fprintf(output, "\nX Coordinate:, ");
		for(auto &particle: allParticles) fprintf(output, "%lf, ", particle.xPos);	
		fprintf(output, "\nY Coordinate:, ");
		for(auto &particle: allParticles) fprintf(output, "%lf, ", particle.yPos);	

		//record all of the data that relates to changes over time (epoch)
		fprintf(output, "\n\nEpoch:, X Error, Y Error, Avg Fitness, Avg Distance from Max, Fraction of Close Particles\n");
		for(int i = 0; i < stats.xErrors.size(); i++){
			fprintf(output, "%d, %lf, %lf, %lf, %lf, %lg\n", i, stats.xErrors[i], stats.yErrors[i], 
					stats.avgPopFitness[i], stats.avgDistance[i], stats.fractionNearMax[i]);
		}
		fprintf(output, "\n\n");

	} else {
		fprintf(stderr, "Could not open the current experiment's data file: %s\n", filename.c_str());
		exit(1);
	}
	fclose(output);
}




Particle::Particle(){
	double maxWidth = WORLDWIDTH / 2; 		//positive bound for width of world space
	double minWidth = maxWidth * (-1); 	//negative bound for width of world space
	double maxHeight = WORLDHEIGHT / 2; 	//positive bound for height of world space
	double minHeight = maxHeight * (-1); 	//negative bound for height of world space
	
	//initialize the known values for each particle
	xPos = generateDouble(minWidth, maxWidth);
	yPos = generateDouble(minHeight, maxHeight);
	xVelocity = 0.0;
	yVelocity = 0.0;
	xBest = xPos;
	yBest = yPos;

}





//for setting up new data files for a given experiment (set of simulations based on a varied parameter)
string newDataFile(const string &currentTest, int problem){
	FILE *output;
	string filename;

	if(problem == 1){
		//create the filename for a new set of experiments that are based on Problem 1
		filename = "Results/Problem1/" + currentTest + ".csv";
	}
	else if(problem == 2){
		//create the filename for a new set of experiments that are based on Problem 2
		filename = "Results/Problem2/" + currentTest + ".csv";
	}
	else{
		fprintf(stderr, "A valid problem number was not specified before an experiment name was declared in parameters.txt\n");
		fprintf(stderr, "Before defining experiments and parameter sets, set the problem number in parameters.txt\n");
		exit(1);
	}

	//if the file already exists, overwrite it for a new set of experiments
	output = fopen(filename.c_str(), "wb");
	
	//add header information about the experiment to the file
	if(output != NULL){
		fprintf(output, "Problem:, %d\nExperiment:, %s\n\n", problem, currentTest.c_str());
	} else {
		fprintf(stderr, "Could not create the data file for test: %s\n", currentTest.c_str());
		exit(1);
	}
	fclose(output);

	return filename;
}


//run many simulations with every combination of parameters in batch_parameters.txt, and rank them based on convergence rate
void runBatchSims(int problem){
	ifstream inputCombos ("batch_parameters.txt");
	SwarmSim *simulation;
	int runCount = 4;
	int totalSims = 0;
	int maxConvergence = 0;
	FILE *bestCombos;
	string outfilename;
	string parameter;
	string inputLine;
	vector <int> particleCountVals;
	vector <double> inertiaVals;
	vector <double> cognitionVals;
	vector <double> socialVals;
	vector <double> maxVelocityVals;
	multimap<int, SwarmSim *> bestSims;
	multimap<int, SwarmSim *>::iterator iter;
	

	//open the input file to start parsing all of the parameters to enumerate
	//inputCombos = fopen("batch_parameters.txt", "r");
	//if(inputCombos == NULL){
	if(!inputCombos){
		fprintf(stderr, "Could not open batch_parameters.txt");
		fprintf(stderr, "Please ensure that batch_parameters.txt exists in this directory and is formatted correctly\n");
		exit(1);
	}

	//parse the parameter file to start grabbing parameter values
	while(getline(inputCombos, inputLine)){
		int iVal = 0;
		double dVal = 0.0;
		istringstream parser(inputLine);
	
		//skip any blank lines
		if(inputLine.length() == 0 || inputLine.find_first_not_of(" \n\t") == string::npos){
			continue;
		}

		//skip lines that are "commented out"
		if(inputLine[0] == '/'){
			continue;
		}

		//determine the parameter being considered
		if(parser >> parameter){
			if(parameter == "PARTICLECOUNT"){
				while(parser >> iVal) particleCountVals.push_back(iVal);
			}
			else if(parameter == "INERTIA"){
				while(parser >> dVal) inertiaVals.push_back(dVal);
			}
			else if(parameter == "COGNITION"){
				while(parser >> dVal) cognitionVals.push_back(dVal);
			}
			else if(parameter == "SOCIAL"){
				while(parser >> dVal) socialVals.push_back(dVal);
			}
			else if(parameter == "MAXVELOCITY"){
				while(parser >> dVal) maxVelocityVals.push_back(dVal);
			}
		}
		else{
			fprintf(stderr, "Could not parse batch_parameters.txt\n");
			fprintf(stderr, "Ensure parameter lines start with PARTICLECOUNT, INERTIA, COGNITION,\n");
			fprintf(stderr, "   SOCIAL, or MAXVELOCITY before appending the respective parameter values\n");
			exit(1);
		}
	}
	//fclose(inputCombos);
	inputCombos.close();

	//if any of the parameter value lists are empty, there is a problem
	if(particleCountVals.size() == 0 || inertiaVals.size() == 0 || cognitionVals.size() == 0 ||
			socialVals.size() == 0 || maxVelocityVals.size() == 0){
		fprintf(stderr, "Could not parse batch_parameters.txt (one parameter was not assigned a set of values)\n");
		fprintf(stderr, "Ensure parameter lines start with PARTICLECOUNT, INERTIA, COGNITION,\n");
		fprintf(stderr, "   SOCIAL, or MAXVELOCITY before appending the respective parameter values\n");
		exit(1);
	}

	//determine the total number of different parameter combinations that were used
	totalSims = (particleCountVals.size() * inertiaVals.size() * cognitionVals.size() * socialVals.size() * maxVelocityVals.size());

	//enumerate all possible combinations of these parameters and run a simulation for each combination
	for(auto &numParticles: particleCountVals){
		for(auto &inertia: inertiaVals){
			for(auto &cognition: cognitionVals){
				for(auto &social: socialVals){
					for(auto &maxVelocity: maxVelocityVals){
						SimParams p;				//holds all of the parameter information
						SimParams *params = &p;		//points to that parameter struct for passing into SwarmSim constructor

						//grab the params
						p.problem = problem;
						p.numParticles = numParticles;
						p.inertia = inertia;
						p.cognition = cognition;
						p.social = social;
						p.maxVelocity = maxVelocity;
						p.seed = 0;
						p.filename = "Batch";

						//run the simulation multiple times (runCount) and average the epoch values at which it converges each time
						int requiredEpochs;
						double averageConvergence = 0;
						double localMaxFraction = 0.0;
						for(int i = 0; i < runCount; i++){
							simulation = new SwarmSim(params);

							//run the simulation and see how long it takes to converge
							simulation->convergenceEpoch = simulation->runSim();
							
							averageConvergence += simulation->convergenceEpoch;
							
							if(simulation->convergedToLocal) localMaxFraction += 1.0;

							//only save the last run of this set of parameters
							if(i < runCount - 1) delete simulation;
						}
						requiredEpochs = rint(averageConvergence / (double)runCount);
						localMaxFraction /= runCount;
						
						//if, on average, the simulation converged and fewer than a 3rd of the runs converged to the local max, consider it for the list of bestSims
						if(requiredEpochs < NUMEPOCHS && localMaxFraction < 0.3){
							//if there are fewer than 10 entries currently in the map, just go ahead and add this one
							if(bestSims.size() < 10) bestSims.insert(make_pair(requiredEpochs, simulation));
							
							//otherwise, check if this simulation is better than one currently in the map
							else{
								//grab the last simulation currently in the map
								SwarmSim *backSim = bestSims.rbegin()->second;

								//if this simulation converged faster than the slowest of the top 10 simulations, replace it 
								if(requiredEpochs < bestSims.rbegin()->first){
									//free the memory for the Sim being removed and then erase its entry in the map
									delete bestSims.rbegin()->second;
									iter = bestSims.end();
									--iter;
									bestSims.erase(iter);

									bestSims.insert(make_pair(requiredEpochs, simulation));
								}

								//otherwise, free the memory for the current simulation and move on
								else delete simulation;
							}
						}

						//ignore simulations that don't converge before NUMEPOCHS
						else{
							//free allocated memory for all simulations that do not converge
							delete simulation;
						}
					}
				}
			}
		}
	}

	/* Print the results to the console and save data to a .csv file */

	if(problem == 1){
		//create the filename for a new set of experiments that are based on Problem 1
		outfilename = "Results/Problem1/BatchRun.csv";
	}
	else if(problem == 2){
		//create the filename for a new set of experiments that are based on Problem 2
		outfilename = "Results/Problem2/BatchRun.csv";
	}

	//if the file already exists, overwrite it for a new batch test
	bestCombos = fopen(outfilename.c_str(), "wb");
	
	//add data for the best simulations of this batch test
	if(bestCombos != NULL){
		fprintf(bestCombos, "Batch Test Top %lu Simulation Parameters for Problem %d\n", bestSims.size(), problem);
		fprintf(bestCombos, "Total Parameter Combinations:, %d\n\n", totalSims);
		printf("\nBest Simulations (when averaged over %d runs):\n", runCount);

		//print out header info to the output file
		fprintf(bestCombos, "Best Simulations (when averaged over %d runs):\n", runCount);
		fprintf(bestCombos, "Rank (Sim #), Avg Convergence Epoch, Particles, Inertia, Cognition, Social, Maximum Velocity\n");

		//print the simulation data in order of average convergence rate and show the parameters used for each
		int num = 1;
		for(iter = bestSims.begin(); iter != bestSims.end(); iter++){
			printf("  [%d] Average Convergence Epoch: %d\n", num, iter->first);
			printf("       Particles: %d, Inertia: %lg, Cognition: %lg, Social: %lg, Maximum Velocity: %lg\n", iter->second->numParticles,
					iter->second->inertia, iter->second->cognition, iter->second->social, iter->second->maxVelocity);

			fprintf(bestCombos, "%d, %d, %d, %lf, %lf, %lf, %lf\n", num, iter->first, iter->second->numParticles,
					iter->second->inertia, iter->second->cognition, iter->second->social, iter->second->maxVelocity);
			
			num++;
		}

		//grab the largest converge epoch for the simulations contained in this list
		maxConvergence = bestSims.rbegin()->first; 

		//print out more header information for the error values associated with each epoch
		fprintf(bestCombos, "\nSim Number, X/Y Error, Epoch:,");
		for(int i = 0; i < maxConvergence; i++) fprintf(bestCombos, "%d,", i);
		fprintf(bestCombos, "\n");

		//print the error values (x and y) associated with each of the simulations listed above at each epoch 
		num = 1;
		for(iter = bestSims.begin(); iter != bestSims.end(); iter++){
			SwarmSim *sim = iter->second;
			Statistics *data = &(sim->stats);

			//x error for last simulation of the set of runs that use this combination of parameters
			fprintf(bestCombos, "%d, x,, ", num);
			for(int i = 0; i < data->xErrors.size(); i++) fprintf(bestCombos, "%lf, ", data->xErrors[i]);
			fprintf(bestCombos, "\n");
			
			//y error for last simulation of the set of runs that use this combination of parameters
			fprintf(bestCombos, ", y,, ");
			for(int i = 0; i < data->yErrors.size(); i++) fprintf(bestCombos, "%lf, ", data->yErrors[i]);
			fprintf(bestCombos, "\n");
			
			num++;
		}
	} else {
		fprintf(stderr, "Could not create the data file for the results of running a a batch test for Problem %d\n", problem);
		exit(1);
	}
	fclose(bestCombos);

	//free the memory that is still allocated for the best sims
	for(iter = bestSims.begin(); iter != bestSims.end(); iter++) delete iter->second;
}

/* Random number generator functions */

//randomize the seed for the global RNG
unsigned int randomize(){
	unsigned int seed = rd();
	globalRNG.seed(seed);
	return seed;
}

//randomly generate either '0' or '1'
char generateBitChar(){
	static uniform_int_distribution<char> c{};
	using parm_t = decltype(c)::param_type;
	return c(globalRNG, parm_t{'0','1'});
}

//randomly generate double in range [lower, upper)
double generateDouble(double lower, double upper){
	static uniform_real_distribution<double> d{};
	using parm_t = decltype(d)::param_type;
	return d(globalRNG, parm_t{lower, upper});
}

//randomly generate an index into a container of length Size
double generateIndex(int Size){
	static uniform_int_distribution<int> d{};
	using parm_t = decltype(d)::param_type;
	return d(globalRNG, parm_t{0, (Size - 1)});
}

//generate random boolean value
bool generateBool(){
	static bernoulli_distribution b{};
	using parm_t = decltype(b)::param_type;
	return b(globalRNG, parm_t{0.5});
}

