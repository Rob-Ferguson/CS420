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

//variables to be used within the network instantiation
struct parameterData{
	int problem;						//the current problem number
	int hiddenLayers;					//the number of hidden layers in the network
	int origNeuronCount;				//the initial neuron count for hidden layers (to be adjusted for different tests)
	int inputCount;						//the number of input neurons (input variables)
	int outputCount;					//the number of output neurons (output variables)
	int epochs;
	double learnRate;					//the network's learning rate
	string trainingFile;				//full path name for the training file
	string validationFile;				//full path name for the validation file
	string testingFile;					//full path name for the testing file
    string filename;                    //stores the name of the csv file currently being written to

};

//a given node in the hidden layers or in the output layer
struct Node{
	double h;
	double sigma;
	double delta;
	double bias;
};

//a given backpropagation neural network to be constructed from one line of parameters.txt 
class Network{
	public:
		/* Network architecture */
		int problem;								//problem 1 or 2 from lab writeup
		int hiddenLayers;							//the number of hidden layers in the network
		int origNeuronCount;						//the initial neuron count for hidden layers (to be adjusted for different tests)
		int inputCount;								//the number of input neurons (input variables)
		int outputCount;							//the number of output neurons (output variables)
		vector <int> neurons;						//the adjusted number of neurons within each hidden layer 
		int epochs;									//the number of times to iterate through all of the input patterns and process them with the network (time steps)
		double learnRate;							//the network's learning rate
		string trainingFile;						//full path name for the training file
		string validationFile;						//full path name for the validation file
		string testingFile;							//full path name for the testing file
		string filename;                            //stores the name of the csv file currently being written to
        vector < vector < vector <double> > > hiddenWeights;	//the weight connections of all hidden nodes/output nodes
		vector < vector <double> > outputWeights;				//the weights (outputs) associated with nodes in the last hidden layer 
		vector < vector <Node *> > hiddenNodes;		//all of the hidden nodes indexed: [layer][node in layer]
		vector <Node *> outputNodes;				//all of the output nodes

		/* Training, validation, and testing input/output combos */
		vector < vector <double> > trainingData;	//contains the inputs followed by output value from the training file 
		vector < vector <double> > validationData;	//contains the inputs followed by output value from the validation file 
		vector < vector <double> > testingData;		//contains the inputs followed by output value from the testing file 

		/* Root mean square error calculations between expected output and network output (validation rmse done after each epoch, testing rmse done after all epochs) */
		vector < vector <double> > rmseValidation;	//indexed: [output node][epoch] where value there is the rmse of that output node after that epoch of training
		vector <double>	rmseTesting;				//indexed: [output node] where value there is the final rmse of that output node after training process

		/* Network member functions */
		Network(struct parameterData *params);		//construct a new network
		void resetNetwork();
		void resetNodeVals();
		void readFileData(int problem);
		void forwardPass(vector <double> &inputs);
		void backwardPass(vector <double> &inputs);
		void updateWeights(vector <double> &inputs);
		void printOutputData();
        void printCSV();

		/* Primary network functions */
		void trainNet();
		void validateNet(int curEpoch);
		void testNet();
};

//utility function prototypes
double calculateSigma(double h);
