#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstdio>
#include <vector>
#include <time.h>
#include <utility>
#include <cmath>
#include <random>

using namespace std;

const int runCount = 10;				//the number of simulations to run for one set of invariant parameters (i.e., rerunning sim)

//global random number engine and random_device for seeding
random_device rd;
mt19937 globalRNG;

//helper function prototypes
void newDataFile(const string &filename, const string &currentTest); //for initializing data files for certain experiments
unsigned int randomize();							//randomize the generator seed (returns that seed)
char generateBitChar();								//generate random char of either '0' or '1'
bool generateBool(); 								//generate a random boolean value
double generateDouble(double lower, double upper);	//generate random double in range [lower, upper)
double generateIndex(int Size); 					//generate a random index into a container of length Size


//parameters from which to construct the current set of (identical) genetic simulation runs
typedef struct {
    unsigned int seed; 				//a specified seed value to use for RNG (seed 0 means generate a new seed)
    int geneLength;					//number of genes (bits) in each genetic string
    int popSize;					//size of the population
    int genCount;					//number of generations to simulate
    double mutationProb;			//probability of a mutation occuring
    double crossProb;				//probability of crossover occuring
    string filenameBase;            //stores the name of the csv file currently being written to (the current test)
} SimParams;

//holds info about the best individual in one generation of a particular run
typedef struct {
	double fitness; 				//fitness value of the best individual in some generation
	int numCorrectBits; 			//the number of correct bits ('1') in this individual's genetic string
	string genes; 					//the individual's actual genetic string
} Best;

//calculated statistics for a single simulation run
typedef struct {
	vector <Best *> bestIndividuals; 	//pointers to structs containing info about the best individual (indexed by generation)
	vector<int> numCorrectBits;			//the number of correct bits in the best individual in each generation (indexed by generation)
	vector <double> avgFitness;			//the average fitness of the population in each generation (indexed by generation)
	vector <double> bestFitness;		//the fitness of the best individual in each generation (indexed by generation)
} Statistics;


//the actual simulation constructed from one set of invariant parameters
class GeneticSim {
	public:
		unsigned int seed; 				//the seed used to generate random variables for this simulation
		int currentRun;					//the run number for this set of simulations based on the same initial parameters
		int currentGen; 				//the current generation (0-indexed)
		int geneLength;					//number of genes (bits) in each genetic string
		int popSize;					//size of the population
		int genCount;					//number of generations to simulate
		double mutationProb;			//probability of a mutation occuring
		double crossProb;				//probability of crossover occuring
		string filename;                //stores the name of the csv file currently being written to (the current test)

		vector <string> population;		//at each index is the genetic string of a given individual in the current generation
		vector <double> fitness;		//at each index is the calculated fitness of the respective individual at population[i]
		vector <double> normFitness;	//the normalized fitness values (indexed by individual) with respect to the total population fitness 
		vector <double> cumulativeNF; 	//the cumulative normalized fitness of all preceding individuals (sum of normFitness[0] up to normFitness[i])
		vector <Statistics * > data;	//[run number] -> Statistics struct for each run that contains important data from each generation per run

		GeneticSim(SimParams *p);
		~GeneticSim();
		void resetSim();
		void runSim();
		unsigned long int convertGeneticString(const string &genes);
		void calculatePopulationFitness();
		double fitnessFunction(const string &individual);
		void calculateStats();
		void createNextGeneration();
		void recordResults();
};

int main(int argc, char **argv){
	GeneticSim *simulation;
	unsigned int seed = 0;
	int lineNum = 0;				//the current line number in parameters.txt
	string inputLine;				//a given line of input (containing network parameters) from parameters.txt
	string currentTest = "Generic";	//current experiment (i.e., varying mutation probability = "VariedMutationProb")
	string filename = "";

	//if the correct number of arguments are given, run one simulation with the given parameters (for runCount repetitions)
	if(argc == 7 || argc == 8){
		SimParams p;				//holds all of the parameter information
		SimParams *params = &p;		//points to that parameter struct for passing into GeneticSim constructor

		//grab the params
		p.popSize = atoi(argv[1]);
		p.genCount = atoi(argv[2]);
		p.geneLength = atoi(argv[3]);
		p.mutationProb = atof(argv[4]);
		p.crossProb = atof(argv[5]);
		currentTest = argv[6];
		p.filenameBase = "Results/" + currentTest + ".csv";
	
		//grab the seed if it is given
		if(argc == 8){
			stringstream buf(argv[7]);
			if(!(buf >> seed)) { fprintf(stderr, "Could not convert given seed to an unsigned int\n"); exit(1); }
		}
		p.seed = seed;

		simulation = new GeneticSim(params);
		
		//run the simulation (with this set of parameters) for the specified number of different runs/iterations
		for(int run = 0; run < runCount; run++){
			simulation->runSim();
			simulation->resetSim();
		}
		
		//save the data for this simulation to a .csv file
		simulation->recordResults();

		//free allocated memory
		delete simulation;
	}

	//if no arguments are given, just accept the parameters from stdin (which will be piped in from parameters.txt) 
	else if(argc == 1){
		cout << "Grabbing parameter values from parameters.txt...\n";

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

			//grab strings representing a new experiment description/name, which are prepended by '*'
			else if(inputLine[0] == '*'){
				currentTest = inputLine;
				currentTest.erase(currentTest.begin());
				printf("\n  Experiment: %s\n", currentTest.c_str());

				//create the file name for this test
				filename = "Results/" + currentTest + ".csv";

				//create a new .csv file for this experiment
				newDataFile(filename, currentTest);

				continue;
			}

			SimParams p;				//holds all of the parameter information
			SimParams *params = &p;		//points to that parameter struct for passing into GeneticSim constructor

			//otherwise, try to individually parse the arguments on each line
			if(parser >> p.popSize >> p.genCount >> p.geneLength >> p.mutationProb
					>> p.crossProb){
				if(currentTest == "Generic") { fprintf(stderr, "The test description was not specified in parameters.txt\n"); exit(1); }
				
				p.filenameBase = filename;
				
				//seed is 0 for batch tests (so constructor will create/record a new random seed for each simulation)
				p.seed = 0;

			}else{
				cout << "Incorrect arguments in parameters.txt on line " << lineNum << endl;
				cout << "Parameters.txt USAGE: ./genetics popSize genCount geneLength mutationProb CrossProb\n";
				exit(1);
			}

			//create a new instance of the genetic simulation for the current set of parameters
			simulation = new GeneticSim(params);

			//run the simulation (with this set of parameters) for the specified number of different runs/iterations
			for(int run = 0; run < runCount; run++){
				simulation->runSim();
				simulation->resetSim();
			}

			//save the data for this simulation to a .csv file
			simulation->recordResults();

			//free allocated memory
			delete simulation;
		}
	}
	else { 
		fprintf(stderr, "Please enter the correct arguments or pipe in parameters.txt in the correct format\n");
		fprintf(stderr, "USAGE: popSize genCount geneLength mutationProb crossProb currentTestName (optional: seed)\n");
		fprintf(stderr, " 		If given, the seed value must be greater than 0\n");
		exit(1); 
	}

	return 0;
}


GeneticSim::GeneticSim(SimParams *p){
	//initialize the primary simulation to utilize the parameters in SimParams
	popSize = p->popSize;
	genCount = p->genCount;
	geneLength = p->geneLength;
	mutationProb = p->mutationProb;
	crossProb = p->crossProb;
	filename = p->filenameBase;
	
	//either create a new, random seed, or use one given on the command line
	if(p->seed == 0) seed = randomize();
	else{
		seed = p->seed;
		globalRNG.seed(seed);
	}
	
	printf("START SEED: %u\n", seed);

	//since this class instance was just instantiated, the currentRun and currentGen are 0 (0-based indexing)
	currentRun = 0;
	currentGen = 0;

	//resize all of the class containers to contain the correct number of elements
	population.resize(popSize);
	fitness.resize(popSize, 0.0);
	normFitness.resize(popSize, 0.0);
	cumulativeNF.resize(popSize, 0.0);
	data.resize(runCount);

	//initialize the starting poulation
	for(string &individual: population){
		individual.reserve(geneLength);
		for(int i = 0; i < geneLength; i++) individual.push_back(generateBitChar());
	}

	//data contains pointers to a Statistics struct for each run of this simulation
	for(auto &statPntr: data){
		statPntr = new Statistics;

		//resize the containers within the Statistics struct (each vector is indexed by generation)
		statPntr->avgFitness.resize(genCount, 0.0);
		statPntr->bestIndividuals.resize(genCount);

		//create new structs for the best individuals in each generation and initialize those values
		for(auto &best: statPntr->bestIndividuals){
			best = new Best;
			best->fitness = 0.0;
			best->numCorrectBits = 0;
			best->genes = "";
		}
	}
}


//free up all the memory that was allocated
GeneticSim::~GeneticSim(){
	//iterate through each of the Statistics * in the data vector
	for(auto &stats: data){
		//iterate through each of the Best * in stats->bestIndividuals and delete them
		for(auto &best: stats->bestIndividuals) delete best;
		
		//delete each Statistics *
		delete stats;
	}
}


//reinitialize the current simulation for a new run with the same parameters
void GeneticSim::resetSim(){
	//update the current run number
	currentRun++;
	
	//reset the current generation number
	currentGen = 0;

	//reset the initial population and generate a new gene string for each individual
	for(string &individual: population){
		for(int i = 0; i < geneLength; i++) individual[i] = generateBitChar();
	}
	
	//reset the current fitness and normalized fitness values
	fitness.assign(popSize, 0.0);
	normFitness.assign(popSize, 0.0);
	cumulativeNF.assign(popSize, 0.0);
}


//run the current simulation with the given parameters
void GeneticSim::runSim(){
	//continue for the specified number of generations
	while(currentGen < genCount){

		//calculate fitness values across the current population/generation
		calculatePopulationFitness();

		//calculate important statistical data for the current population/generation
		calculateStats();
		
		//update the population for the next generation
		createNextGeneration();
	}
}


//convert an individual's genetic string to an integer value (binary -> dec conversion)
unsigned long int GeneticSim::convertGeneticString(const string &bitString){
	unsigned long int bitSum = 0;
	string reversedBits = bitString;
	
	//reverse the string to invert the order of most significant and least significant bits
	reverse(reversedBits.begin(), reversedBits.end());

	for(int i = 0; i < geneLength; i++){
		if(reversedBits[i] == '1') bitSum += pow(2, i);
	}
	
	return bitSum;
}


//calculate the fitness and normalized fitness values for the current population
void GeneticSim::calculatePopulationFitness(){
	double fitnessSum = 0.0; 		//the sum of individual fitness values across the population
	
	//calculate the fitness level of each individual in the current population
	for(int i = 0; i < popSize; i++){
		fitness[i] = fitnessFunction(population[i]);
		fitnessSum += fitness[i];
	}

	//normalize each of these fitness values with respect to total population fitness
	for(int i = 0; i < popSize; i++){
		normFitness[i] = fitness[i] / fitnessSum;

		//for each individual, keep a cumulative sum of all normalized fitness values for all preceding individuals
		cumulativeNF[i]	= ((i == 0) ? (normFitness[i]) : (normFitness[i] + cumulativeNF[i - 1]));
	}
}


//determine an individual's fitness level according to the specified fitness function
double GeneticSim::fitnessFunction(const string &individual){
	return pow(((double)convertGeneticString(individual) / pow(2, (double)geneLength)), 10);
}


//calculate the statistical data for the current generation
void GeneticSim::calculateStats(){
	int bestIndividualIndex = 0;
	double bestFitness = 0.0;

	//grab a handle to the Statistics struct for this run
	Statistics *stats = data[currentRun];
	
	//grab a handle to the Best struct for the current generation in this run
	Best *topIndividual = stats->bestIndividuals[currentGen];
	
	//determine the average fitness of the whole population (and determine the most fit individual)
	for(int i = 0; i < popSize; i++){
		//running sum of all fitness values
		stats->avgFitness[currentGen] += fitness[i];
		
		//grab the current best fitness value and the population index for that individual
		if(bestFitness < fitness[i]){
			bestFitness = fitness[i];
			bestIndividualIndex = i;
		}
	}

	//actually calculate the average fitness
	stats->avgFitness[currentGen] /= ((double)popSize);

	//set the values of the best individual according to the bestIndividualIndex that was just found above
	topIndividual->fitness = bestFitness;
	topIndividual->genes = population[bestIndividualIndex];
	for(int i = 0; i < geneLength; i++){
		if(topIndividual->genes[i] == '1') topIndividual->numCorrectBits += 1;
	}
}


//produce the next generation from the current population
void GeneticSim::createNextGeneration(){
	int mother_index = 0; 					//index of one parent in the population vector
	int father_index = 0; 					//index of the other parent in the population vector
	int crossPoint; 						//the crossover point for mating
	double rand1, rand2; 					//random numbers to use for determining parents
	double tryCrossover, tryMutate; 		//random numbers to determine if crossover or mutation will take place
	vector<string> offspring; 				//the offspring that will become the population of the next generation
	
	//initialize the offspring container
	offspring.resize(popSize);
	for(string &individual: offspring){
		individual.resize(geneLength);
	}

	//update half of the current population with new offspring
	for(int x = 0; x < (popSize / 2); x++){
		//choose two individuals from the current generation to be parents (start with the first parent)
		rand1 = generateDouble(0, 1);

		//increment mother_index until the cumulative normalized fitness at that index is not less than rand1
		while(cumulativeNF[mother_index] < rand1) mother_index++;

		//now determine the second index (keep trying until mother_index != father_index)
		do{
			father_index = 0;
			rand2 = generateDouble(0, 1);
			
			//increment mother_index until the cumulative normalized fitness at that index is not less than rand1
			while(cumulativeNF[father_index] < rand2) father_index++;
		} while(father_index == mother_index);
		

		//now that parents are chosen, mate them and determine how to handle crossover (based on crossProb)
		tryCrossover = generateDouble(0, 1);
		
		//if the random number is less than the given crossover probability, commence crossover
		if(tryCrossover < crossProb){
			//choose a random crossover index
			crossPoint = generateIndex(geneLength);

			for(int g = 0; g < geneLength; g++){
				//before the crosspoint, child1 gets genes from mother, and child2 gets genes from father
				if(g < crossPoint){
					offspring[x * 2][g] = population[mother_index][g];
					offspring[x * 2 + 1][g] = population[father_index][g];
				}
				//after the crosspoint, child1 gets genes from father, and child2 gets genes from mother
				else{
					offspring[x * 2][g] = population[father_index][g];
					offspring[x * 2 + 1][g] = population[mother_index][g];
				}

			}
		}
		//the random number was not less than the crossover probability, so simply copy the parents as they were before
		else{
			offspring[x * 2] = population[mother_index];
			offspring[x * 2 + 1] = population[father_index];
		}
	}

	
	//now that crossover is done, loop through individuals in the offspring generation and handle random genetic mutations
	for(string &individual: offspring){
		for(char &bit: individual){
			tryMutate = generateDouble(0, 1);
			
			//if the random number is less than the mutation probability, flip the current bit in the current child's genetic string
			if(tryMutate < mutationProb){
				bit = ((bit == '0') ? ('1') : ('0'));
			}
		}
	}

	//update the current population
	population = offspring;

	//increase the generation counter to reflect the new population
	currentGen++;
}


void GeneticSim::recordResults(){
	int run = 0;
	ofstream output;
	Best *top;
	
	printf("Generating Output Data File: %s\n", filename.c_str());

	//open the output file and append data to the end of it
	output.open(filename.c_str(), ios_base::app);
	
	//write info about the parameters for the current simulation
	output << "Simulation Parameters:\n";
	output << "Seed:, " << seed << "\n";
	output << "Population Size:, " << popSize << "\n";
	output << "Generations:, " << genCount << "\n";
	output << "Number of Genes:, " << geneLength << "\n";
	output << "Mutation Probability:, " << mutationProb << "\n";
	output << "Crossover Probability:, " << crossProb << "\n";

	//print data contained in each Statistics struct for each run of this simulation
	for(auto &statPntr: data){
		output << "Run Number, Generation, Average Fitness, Best Fitness, Number of Correct Bits, Best Genetic String\n";
		
		printf(" Run: %d\n", run);

		//iterate through each generation
		for(int gen = 0; gen < genCount; gen++){
			top = statPntr->bestIndividuals[gen];
			
			//print the statistical data for each generation of the current run
			output << run << "," << gen << "," << statPntr->avgFitness[gen] << "," << top->fitness
				<< "," << top->numCorrectBits << ",'" << top->genes << "'\n";

			
			/*
			printf("   AvgFit: %lf\n", statPntr->avgFitness[gen]);
			printf("   Best: %s\n", top->genes.c_str());
			printf("    -> Fitness: %lf\n", top->fitness);
			printf("    -> CorrectBits: %d\n", top->numCorrectBits);
			*/
		}
		run++;
		output << "\n";
	}
}


//for setting up new data files for a given experiment (set of simulations based on a varied parameter)
void newDataFile(const string &filename, const string &currentTest){
	FILE *output;

	//if the file already exists, overwrite it for a new set
	output = fopen(filename.c_str(), "wb");
	
	//add header information about the experiment to the file
	if(output != NULL){
		fprintf(output, "Experiment:,%s\n\n", currentTest.c_str());
	} else {
		fprintf(stderr, "Could not create the data file for test: %s\n", currentTest.c_str());
		exit(1);
	}
	fclose(output);
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

