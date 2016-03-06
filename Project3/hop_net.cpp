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

typedef vector <vector <double> > basin_histogram;

//for a single instance of a hopfield network with 50 random patterns to imprint
class Hopfield{
	public:
		int run_number;
		int N;								//the number of neurons in the network
		vector <vector <int> > patterns;	//all of the patterns to be imprinted on network
		vector <vector <double> > weights;	//the weight between every (i, j) neuron pair for k imprinted patterns
		vector <vector <double> > *wm_pntr;	//pointer to the current weight matrix
		vector <int> network;				//the actual network of neurons
		vector <int> num_stable;			//index = current number of imprinted patterns, val = number of stable patterns
		vector <double> prob_unstable;		//the probability of unstable imprints for k (index) imprinted patterns	
		vector <int> bit_flips;
		vector <vector <double> > basin_size;	//the estimated basin of attraction size for each number of imprinted patterns

		Hopfield(int run);
		void build_wm(int k);				//build up a weight matrix from the first k patterns imprinted
		void test_networks(int k);			//calculate the fraction of unstable patterns
		bool determine_stability();			//determine the stability of a given pattern in the network
		int calc_basin(int pattern_index);	//estimate the basin size for an imprinted pattern among k imprinted patterns
		bool update_net(vector <int> net, const vector <int> &pat);	//update a network and see if it converges to the given pattern
		void record_data();
};

//the average values from running 50 different simulations of the hopfield network
class avg_Hopfield{
	public:
		int num_simulations;				//total number of simulations to run
		int num_patterns;					//total number of patterns to be imprinted in each simulation
		vector <int> avg_num_stable;		//the average number of stable imprinted patterns each number of total imprinted patterns
		vector <double> avg_prob_unstable;	//the average probability of an unstable pattern for each number of total imprinted patterns
		
		vector <vector <double> > all_basin_sizes;	//the estimated basin of attraction size for each number of imprinted patterns
		vector <basin_histogram *> b_histograms;

		avg_Hopfield();
		void record_avg_data();
};


int main(int argc, char **argv){
	//initialize random seed
	srand(time(NULL));

	avg_Hopfield main_test;						//to store average information about all 50 different hopfield simulations
	
	//start running simulations
	for(int x = 0; x < main_test.num_simulations; x++){
		Hopfield *h_net = new Hopfield(x + 1);	//information about one instance of the hopfield simulation

		for(int k = 0; k < h_net->patterns.size(); k++){
	
			//cout << "Number of imprinted patterns: " << k+1 << "\n\n";
			
			//build the associated weight matrix for the number of imprinted patterns k
			h_net->build_wm(k);

			//initialize the net to each of those patterns and determine how many are still stable
			h_net->test_networks(k);	
		}

		//save the resulting data for this simulation
		h_net->record_data();

		//grab that data to be averaged over 50 runs (add to cumulative sum for now)
		for(int i = 0; i < h_net->num_stable.size(); i++){
			main_test.avg_num_stable[i] += h_net->num_stable[i];
			main_test.avg_prob_unstable[i] += h_net->prob_unstable[i];
		}

		for(int i = 0; i < main_test.all_basin_sizes.size(); i++){
			for(int j = 0; j < main_test.all_basin_sizes[i].size(); j++){
				main_test.all_basin_sizes[i][j] += h_net->basin_size[i][j];
			}
		}
		
		
		//main_test.b_histograms[x] = new basin_histogram;
		//*(main_test.b_histograms[x]) = h_net->basin_size;

		delete h_net;
	}

	//calculate the actually average values for each number of imprinted patterns over all 50 simulations
	for(int i = 0; i < main_test.avg_num_stable.size(); i++){
		main_test.avg_num_stable[i] /= main_test.num_simulations;
		main_test.avg_prob_unstable[i] /= main_test.num_simulations;
	}


	for(int i = 0; i < main_test.all_basin_sizes.size(); i++){
		for(int j = 0; j < main_test.all_basin_sizes[i].size(); j++){
			main_test.all_basin_sizes[i][j] /= main_test.num_simulations;
			main_test.all_basin_sizes[i][j] /= (i + 1);
		}
	}
	/*
	   for(int i = 0; i < main_test.b_histograms.size(); i++){
		for(int x = 0; x < main_test.b_histograms[i])
	}
	*/

	//record the overall results
	main_test.record_avg_data();

	return 0;

}

/* avg_Hopfield class functions */

avg_Hopfield::avg_Hopfield(){
	//set the number of patterns to be considered in each simulation
	num_patterns = 50;

	/*
	 *	CHANGE THE NUMBER OF SIMULATIONS BACK TO 50 AFTER DEBUGGING
	 *
	 * */

	num_simulations = 5;

	//resize vectors appropriately
	avg_num_stable.resize(num_patterns, 0);
	avg_prob_unstable.resize(num_patterns, 0);

	all_basin_sizes.resize(num_patterns);
	for(int i = 0; i < all_basin_sizes.size(); i++) all_basin_sizes[i].resize(51, 0);

	b_histograms.resize(num_simulations);
}


void avg_Hopfield::record_avg_data(){
	string datafilename = "Average_Imprint_Stability.csv";
	ofstream data;
	char outputline[1024];

	//append new data to the MasterData file associated with the experiment number being run
	data.open(datafilename, ios_base::out);
	if(!data.is_open()){
		fprintf(stderr, "Could not open current data vile (.csv)\n");
		exit(1);
	}

	//header info
	data << "Number of Imprints, Average Fraction of Unstable Imprints, Average Number of Stable Imprints\n";

	//print out the data associated with each column specified above
	for(int i = 0; i < avg_prob_unstable.size(); i++){
		sprintf(outputline, "%d, %f, %d\n", i + 1, avg_prob_unstable[i], avg_num_stable[i]);
		data << outputline;
	}
	data << "\n\n";

	data << "Basin Sizes\nNumber of Imprints, ";
	for(int i = 0; i < all_basin_sizes[0].size(); i++) data << i << ", ";
	data << endl;
	for(int i = 0; i < all_basin_sizes.size(); i++){
		data << (i + 1) << ", ";
		for(int j = 0; j < all_basin_sizes[i].size(); j++){
			//data << (all_basin_sizes[i][j])/(i + 1) << ", ";	
			data << all_basin_sizes[i][j] << ", ";	
		}
		data << endl;
	}


	data.close();
}

/* Hopfield class functions */

Hopfield::Hopfield(int run){
	//save the current simulation number
	run_number = run;	
	
	//initialize the number of neurons in the hopfield network
	N = 100;

	//resize the vectors that hold each pattern (50 total patterns)
	patterns.resize(50);
	for(int i = 0; i < patterns.size(); i++){
		patterns[i].resize(N);
	}

	//randomly populate every pattern vector with 1 or -1
	for(int i = 0; i < patterns.size(); i++){
		for(int j = 0; j < patterns[i].size(); j++){
			if((rand() % 2) == 1) patterns[i][j] = 1;
    		else patterns[i][j] = -1;
		}
	}

	//resize the 100x100 vector of weights to hold each weight between pairs of neurons
	weights.resize(N);
	for(int i = 0; i < weights.size(); i++){
		weights[i].resize(N, 0);
	}

	//resize the vectors that hold number of stable imprints, probability of unstable patterns, and estimated basin size
	num_stable.resize(patterns.size(), 0);
	prob_unstable.resize(patterns.size(), 0);
	basin_size.resize(patterns.size());
	for(int i = 0; i < basin_size.size(); i++) basin_size[i].resize(51, 0);

	//fill bit_flips array with 0-99 (to be shuffled later)
	bit_flips.resize(100);
	for(int x = 0; x < 100; x++) bit_flips[x] = x;
}


//determine the stability of each pattern in the neural network
void Hopfield::test_networks(int k){
	bool stable;			//determine if a given pattern is stable in the network
	int count_stable = 0;	//reset the count of stable imprinted patterns for a new set
	int basin;				//the estimated basin size for a pattern among k imprinted patterns

	//set the neural net equal to each pattern and determine its stability
	for(int i = 0; i <= k; i++){
		network = patterns[i];
		stable = determine_stability();
		if(stable){
			count_stable++;
			basin = calc_basin(i);
			//cout << "Imprint " << i << "/" << k << ", basin size: " << basin << endl;
			basin_size[k][basin] += 1;
		}
		//the basin size for an unstable pattern is 0, so increase count of 0-sized basins for that number of imprinted patterns
		else basin_size[k][0] += 1;
	}

	//keep track of the number of stable patterns for the k imprinted patterns
	num_stable[k] = count_stable;
	//calculate the probability of an unstable pattern occuring with k imprinted patterns
	prob_unstable[k] = 1 - (count_stable / (double) (k + 1));	
}

bool Hopfield::determine_stability(){
	double h = 0;			//the local field value for a given neuron
	int new_state = 0;		//representing a new state after local field calculations

	//for each neuron i in the network, calculate its local field value
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			h += (weights[i][j] / N) * network[j];
		}
		//cout << "Cell " << i << " , h = " << h << endl;

		//if the local field is negative, the new state for neuron i would be -1
		if(h < 0) new_state = -1;
		//if the local field is positive, the new state for neuron i would be 1
		else if(h >= 0) new_state = 1;
		
		//if the new state is different than the current state, the pattern is not stable
		if(network[i] != new_state)	return false;

		h = 0;
	}
	return true;
}


//constructs a weight matrix for k imprinted patterns
void Hopfield::build_wm(int k){
	//for each new imprinted pattern, calculate the weights between each neuron pair (no self-coupling)
	//add that product of states for (i, j) to a cumulative sum at (i, j) in the weight matrix
	//when using a weight later, must divide by the number of neurons N
	for(int i = 0; i < patterns[k].size(); i++){
		for(int j = 0; j < patterns[k].size(); j++){
			if(i != j){
				weights[i][j] += (patterns[k][i] * patterns[k][j]);
			}
		}
	}
}


int Hopfield::calc_basin(int pattern_index){
	vector <int> net;
	int i;		
	int cur_basin;
	int avg_basin = 0;
	int num_iterations = 5;		//the number of times to recalculate the basin size for averaging 
	bool match;					//true if an altered network converges back to its original state after 10 iterations

	//calculate the estimated basin size for several bit_flip permutations and average them to get a general basin size
	for(int x = 0; x < num_iterations; x++){
		//reshuffle the indices that will be used to flip neuron state values 
		random_shuffle(bit_flips.begin(), bit_flips.end());

		//reinitialize a hopfield network to match the current pattern being considered
		net = patterns[pattern_index];

		for(i = 0; i < 50; i++){
			//flip the bits in the network that correspond to the first i bit indices in the bit_flips vector
			if(net[bit_flips[i]] == 1) net[bit_flips[i]] = -1;
			else if(net[bit_flips[i]] == -1) net[bit_flips[i]] = 1;
			else { fprintf(stderr, "Could not flip bits in network\n"); exit(1); }

			//update the neural network for 10 iterations and see if it hasn't converged back to the original imprinted pattern
			match = update_net(net, patterns[pattern_index]);
			if(!match || i == 49){
				cur_basin = i + 1;
				break;
			}
		}
		//cout << "   basin: " << cur_basin << endl;
		avg_basin += cur_basin;
	}
	avg_basin /= num_iterations;
	//cout << "   avg_basin: " << avg_basin << endl;

	return avg_basin;
}


bool Hopfield::update_net(vector <int> net, const vector <int> &pat){
	double h = 0;			//the local field value for a given neuron
	int new_state = 0;		//representing a new state after local field calculations

	//update the neural network for 10 iterations using the local fields of each neuron
	for(int x = 0; x < 10; x++){

		//if(net == pat) cout << "net and pat equal at x = " << x << endl;
		//else cout << "not equal at x = " << x << endl;
		
		//for each neuron i in the network, calculate its local field value
		for(int i = 0; i < N; i++){
			for(int j = 0; j < N; j++){
				h += (weights[i][j] / N) * net[j];
			}

			//if the local field is negative, the new state for neuron i would be -1
			if(h < 0) new_state = -1;
			//if the local field is positive, the new state for neuron i would be 1
			else if(h >= 0) new_state = 1;

			//set the new state for this neuron
			net[i] = new_state;

			h = 0;
		}
		//if the network has converged to the pattern, go ahead and return true
		if(net == pat) return true;
	}
	
	return false;
}


void Hopfield::record_data(){
	string datafilename = "Imprint_Stability.csv";
	ofstream data;
	char outputline[1024];

	//append new data to the MasterData file associated with the experiment number being run
	data.open(datafilename, ios_base::app);
	if(!data.is_open()){
		fprintf(stderr, "Could not open current data vile (.csv)\n");
		exit(1);
	}

	//header info
	data << "Run Number " << run_number << endl; 
	data << "Number of Imprints, Stable Imprints, Fraction of Unstable Imprints\n";

	//print out the data associated with each column specified above
	for(int i = 0; i < prob_unstable.size(); i++){
		sprintf(outputline, "%d, %d, %f\n", i + 1, num_stable[i], prob_unstable[i]);
		data << outputline;
	}
	data << "\n";

	data.close();
}
