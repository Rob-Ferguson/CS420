#include "networkdef.h"


int main(int argc, char **argv){
	int currentProblem = 0;			//valid problems are 1 and 2
	int lineNum = 0;				//the current line number in parameters.txt
	int networkNumber = 0;			//the current network number within a given set for an experiment
	string inputLine;				//a given line of input (containing network parameters) from parameters.txt
	string currentTest = "Generic";	//current experiment (i.e. effect of varying number of hidden layers = "VaryingHiddenLayers")
	string filename;
	int varyAmount = 0;				//the amount to increase the node count by in a given layer (only if the -v flag is given)
	string vpos = "";				//"e" or "o" for whether to vary neurons in even or odd layers


	//initialize random seed
	srand(time(NULL));
		
	//parse the input parameters line by line from stdin (piped from parameters.txt) for each new network's parameters
	while(getline(cin, inputLine)){
		lineNum++;

		istringstream parser(inputLine);
		
		//skip any blank lines
		if(inputLine.length() == 0 || inputLine.find_first_not_of(" \n\t") == string::npos){
			varyAmount = 0;
			vpos = "";
		
			continue;
		}

		//skip lines that are "commented out"
		if(inputLine[0] == '/'){
			varyAmount = 0;
			vpos = "";
			
			continue;
		}

		//see if we are looking at a new problem number
		else if(inputLine == "Problem 1"){
			currentProblem = 1;
			networkNumber = 0;
			varyAmount = 0;
			vpos = "";
			printf("\n%s\n", inputLine.c_str());

			continue;
		}
		else if(inputLine == "Problem 2"){
			currentProblem = 2;
			networkNumber = 0;
			printf("\n%s\n", inputLine.c_str());
			varyAmount = 0;
			vpos = "";
			
			continue;
		}

		//grab strings representing a new experiment description/name, which are prepended by '*'
		else if(inputLine[0] == '*'){
			currentTest = inputLine;
			currentTest.erase(currentTest.begin());
			printf("\n  Experiment: %s\n", currentTest.c_str());
			networkNumber = 0;
		
			//create the file name for this test
			filename = "Results/Problem" + to_string(currentProblem) + "/" + currentTest + ".csv";

			//create a new .csv file for this experiment
			newDataFile(filename, currentTest, currentProblem);
			
			varyAmount = 0;
			vpos = "";

			continue;
		}

		//also search for variants of the "-v" flag, which varies the neuron count in certain layers
		else if(inputLine[0] == '-' && inputLine[1] == 'v'){
			string flag;
			istringstream fparse(inputLine);

			//parse the -v flag line
			if(!(fparse >> flag)) { fprintf(stderr, "Could not parse out -v flag\n"); exit(1); }
			if(!(fparse >> vpos)) { fprintf(stderr, "Could not parse 'e' or 'o' within -v flag\n"); exit(1); }
			if(!(fparse >> varyAmount)){ 
				//if a value cannot be grabbed, just initialize it to 1
				varyAmount = 1;
			}
			
			continue;
		}

		struct parameterData p;				//holds all of the parameter information
		struct parameterData *params = &p;	//points to that parameter struct for passing into Network constructor
		
		//otherwise, try to individually parse the arguments on each line
		if(parser >> p.hiddenLayers >> p.origNeuronCount >> p.learnRate >> p.epochs
				>> p.trainingFile >> p.validationFile >> p.testingFile){
			
			//increment the network counter for this set
			networkNumber++;

			//set other variables that do not need to be read in each line (e.g., the current problem number, test name)
			if(currentProblem == 0){ fprintf(stderr, "The problem number was not correctly initialized from parameters.txt\n"); exit(1); }
			else if(currentTest == "Generic") { fprintf(stderr, "The test description was not specified in parameters.txt\n"); exit(1); }
			
			p.problem = currentProblem;
			p.filename = filename;
			p.networkNumber = networkNumber;
			p.vpos = vpos;
			p.varyAmount = varyAmount;
		}else{
			cout << "Incorrect arguments in parameters.txt on line " << lineNum << endl;
			cout << "USAGE: ./backprop hiddenlayers neurons learningrate epochs";
			cout << " trainingfilename validationfilename testingfilename" << endl; 
			exit(1);
		}

		printf("\n   Starting Network %d -> ", p.networkNumber);

		//create a new instance of the network (and initialize random weights/biases)
		Network n(params);

		//read in all of the training, validation, and testing data
		n.readFileData(n.problem);

		//print out important network parameters 
		printf("Hidden Layers: %d, Neurons/Layer: { ", n.hiddenLayers);
		for(int x = 0; x < n.neurons.size(); x++) printf("%d ", n.neurons[x]);
		printf("}, Learning Rate: %0.2lf, Epochs: %d\n", n.learnRate, n.epochs);
		
		//start processing network (run through inputs for the specified number of epochs)
		for(int i = 0; i < n.epochs; i++){
			n.trainNet();
			n.validateNet(i);
		}
		n.testNet();

		//print the results for this network to a .csv file
		n.printCSV();

		varyAmount++;
	}

	return 0;
}

//constructor for a completely new network
Network::Network(struct parameterData *params){
	//initialize random number generator	
	uniform_real_distribution <double> unif(-0.1, 0.1);
	random_device rand_dev;
	mt19937 rand_engine(rand_dev());
	
	//grab all of the parameter information from the struct
	problem = params->problem;
	networkNumber = params->networkNumber;
	hiddenLayers = params->hiddenLayers;
	origNeuronCount = params->origNeuronCount;
	epochs = params->epochs;
	learnRate = params->learnRate;
	trainingFile = params->trainingFile;
	validationFile = params->validationFile;
	testingFile = params->testingFile;
	filename = params->filename;	

	//set the number of function inputs/outputs based on the problem number
	if(problem == 1){
		inputCount = 2;
		outputCount = 1;
	}
	else if(problem == 2){
		inputCount = 3;
		outputCount = 1;
	}
	else{
		fprintf(stderr, "Could not properly determine inputs/outputs for problem number %d\n", problem);
		exit(1);
	}

	/*
	cout << "hiddenlayers: " << hiddenLayers << endl;
	cout << "neurons: " << origNeuronCount << endl;
	cout << "learnrate: " << learnRate << endl;
	cout << "epochs: " << epochs << endl;
	cout << "input_count: " << inputCount << endl;
	cout << "outputCount: " << outputCount << endl;
	cout << "trainingfile: " << trainingFile << endl;
	cout << "validationfile: " << validationFile << endl;
	cout << "testingfile: " << testingFile << endl;
	*/
	
	/* Resize all of the network vectors */
	neurons.resize(hiddenLayers, origNeuronCount);
	
	//if specified, vary the number of neurons in either even or odd hidden layers
	if(params->vpos == "e"){
		for(int i = 0; i < neurons.size(); i += 2) neurons[i] += params->varyAmount;
	}
	else if(params->vpos == "o"){
		for(int i = 1; i < neurons.size(); i += 2) neurons[i] += params->varyAmount;
	}


	outputWeights.resize(outputCount);

	//first dimmension of hiddenWeights -> the hidden layer number
	hiddenWeights.resize(hiddenLayers);
	
	//second dimension of hiddenWeights -> the neuron in that layer
	for(int i = 0; i < hiddenWeights.size(); i++){
		hiddenWeights[i].resize(neurons[i]);
		
		//third dimension of hiddenWeights -> the weight associated with a node from previous layer (input or hidden)
		if(i == 0){
			//the first hidden layer nodes are connected to the input nodes
			for(int j = 0; j < hiddenWeights[i].size(); j++){
				hiddenWeights[i][j].resize(inputCount);
				
				//initialize the weights to be random between -0.1 and 0.1
				for(int k = 0; k < hiddenWeights[i][j].size(); k++){
					hiddenWeights[i][j][k] = unif(rand_engine);
				}
			}
		} else{
			//the other hidden layer nodes are connected to the previous layer's nodes
			for(int j = 0; j < hiddenWeights[i].size(); j++){
				hiddenWeights[i][j].resize(neurons[i - 1]);
				
				//initialize the weights to be random between -0.1 and 0.1
				for(int k = 0; k < hiddenWeights[i][j].size(); k++){
					hiddenWeights[i][j][k] = unif(rand_engine);
				}
			}
		}
	}

	//the weights of output nodes to the last layer of hidden nodes
	for(int i = 0; i < outputWeights.size(); i++){
		outputWeights[i].resize(hiddenWeights[hiddenLayers - 1].size());

		//initialize the weights to be random between -0.1 and 0.1
		for(int j = 0; j < outputWeights[i].size(); j++){
			outputWeights[i][j] = unif(rand_engine);
		}

	}

	//resize the hiddenNodes [number of hidden layers][nodes in each hidden layer]
	hiddenNodes.resize(hiddenLayers);
	for(int i = 0; i < hiddenLayers; i++){
		hiddenNodes[i].resize(neurons[i]);
	}

	//resize the outputNodes
	outputNodes.resize(outputCount);
	
	/* Populate the Node vectors with the new Node structs */
	for(int i = 0; i < hiddenNodes.size(); i++){
		for(int j = 0; j < hiddenNodes[i].size(); j++){
			Node *aNode = new Node;
			aNode->bias = unif(rand_engine);
			aNode->h = 0;
			aNode->sigma = 0;
			aNode->delta = 0;
			hiddenNodes[i][j] = aNode;
		}
	}
	for(int i = 0; i < outputNodes.size(); i++){
		Node *aNode = new Node;
		aNode->bias = unif(rand_engine);
		aNode->h = 0;
		aNode->sigma = 0;
		aNode->delta = 0;
		outputNodes[i] = aNode;
	}

	/* Resize the vectors that contain the root mean square error information */
	rmseValidation.resize(outputCount);
	for(int i = 0; i < outputCount; i++) rmseValidation[i].resize(epochs);
	
	rmseTesting.resize(outputCount);
}


//for resetting the h, sigma, and delta values of each node
void Network::resetNodeVals(){
	//run through hidden layer nodes and reinitialize their values
	for(int i = 0; i < hiddenNodes.size(); i++){
		for(int j = 0; j < hiddenNodes[i].size(); j++){
			hiddenNodes[i][j]->h = 0;
			hiddenNodes[i][j]->sigma = 0;
			hiddenNodes[i][j]->delta = 0;
		}
	}
	//run through the output node(s) and reinitialize their values
	for(int i = 0; i < outputNodes.size(); i++){
		outputNodes[i]->h = 0;
		outputNodes[i]->sigma = 0;
		outputNodes[i]->delta = 0;
	}
}


//after each training epoch, the network is validated to measure accuracy gains/losses over time
void Network::validateNet(int curEpoch){
	vector <double> errorSums;					//the separate sum of expected/calculated output error for each output node (Problems 1/2 only have one output node)
	errorSums.resize(outputNodes.size(), 0);	//initialize each of these sums to 0
	
	//random check to ensure that numOutputs actually corresponds to the number of outputs in each line of the data file
	if(outputCount != (validationData[0].size() - inputCount)) { fprintf(stderr, "The specified number of outputs doesn't match the data files\n"); exit(1); }

	/* For each set inputs from the validation data file, initialize the network to use those inputs
		and then perform a forward pass, calculating outputs for every neuron */
	for(int i = 0; i < validationData.size(); i++){
		forwardPass(validationData[i]);

		//add to the total accumulated error (between expected outputs and network outputs) for each output node with respect to each new set of validation inputs 
		for(int o = 0; o < outputNodes.size(); o++){
			//inputs and expected output(s), respectively, are stored in same inputs vector(validationData[i]). Expected outputs start at position validationData[i][inputCount].
			//error is estimated using squared Euclidean distance between expected output and calculated output (an output node's sigma value)
			errorSums[o] += pow((validationData[i][inputCount + o] - outputNodes[o]->sigma), 2);
		}
	}

	//after all validation input patterns/sets have been evaluated, calculate the root mean square error for each output (Problems 1/2 will only have one output)
	for(int o = 0; o < outputNodes.size(); o++){
		rmseValidation[o][curEpoch] = sqrt((1.0 / (2.0 * validationData.size())) * errorSums[o]);	
	}

	//print the current validation RMSE to the screen every 1000 epochs
	if(curEpoch % 1000 == 0){
		printf("     Epoch %d\n", curEpoch);
		for(int o = 0; o < outputNodes.size(); o++) printf("       Output[%d] Validation RMSE: %.10lf\n", o, rmseValidation[o][curEpoch]);
	}
}


//after all of the training/validating has finished, test the final accuracy of the network
void Network::testNet(){
	vector <double> errorSums;					//the separate sum of expected/calculated output error for each output node (Problems 1/2 only have one output node)
	errorSums.resize(outputNodes.size(), 0);	//initialize each of these sums to 0
	
	//random check to ensure that numOutputs actually corresponds to the number of outputs in each line of the data file
	if(outputCount != (testingData[0].size() - inputCount)) { fprintf(stderr, "The specified number of outputs doesn't match the data files\n"); exit(1); }

	/* For each set inputs from the validation data file, initialize the network to use those inputs
		and then perform a forward pass, calculating outputs for every neuron */
	for(int i = 0; i < testingData.size(); i++){
		forwardPass(testingData[i]);

		//add to the total accumulated error (between expected outputs and network outputs) for each output node with respect to each new set of testing inputs 
		for(int o = 0; o < outputNodes.size(); o++){
			//inputs and expected output(s), respectively, are stored in same inputs vector(testingData[i]). Expected outputs start at position testingData[i][inputCount].
			//error is estimated using squared Euclidean distance between expected output and calculated output (an output node's sigma value)
			errorSums[o] += pow((testingData[i][inputCount + o] - outputNodes[o]->sigma), 2);
		}
	}

	printf("     Testing\n");
	//after all testing input patterns/sets have been evaluated, calculate the root mean square error for each output (Problems 1/2 will only have one output)
	for(int o = 0; o < outputNodes.size(); o++){
		rmseTesting[o] = sqrt((1.0 / (2.0 * testingData.size())) * errorSums[o]);	
		printf("      Output[%d] RMSE: %.10lf\n", o, rmseTesting[o]);
	}
}


//(forward pass) for computing the outputs of each successive node in the network for one pass
void Network::forwardPass(vector <double> &inputs){
	Node *curNeuron;
	
	//reset all of the output values within each node that were calculated from the previous set of training patterns
	resetNodeVals();

	//iterate through all the nodes and compute their output values hiddenNodes[layer][node]
	for(int l = 0; l < hiddenNodes.size(); l++){
		for(int n = 0; n < hiddenNodes[l].size(); n++){
			curNeuron = hiddenNodes[l][n];

			//if this is the first hidden layer, look back to the inputs
			if(l == 0){
				//iterate through all of the inputs in the input layer
				for(int i = 0; i < inputCount; i++){
					//update the h values (with respect to the current node and each input value)
					curNeuron->h += (hiddenWeights[l][n][i] * inputs[i]);
				}
			}
			//if not the first hidden layer, look back to the previous hidden layer
			else{
				//iterate through the nodes in the previous hidden layer
				for(int i = 0; i < hiddenNodes[l-1].size(); i++){
					//update the h values (with respect to the current node and each node in the previous layer)
					curNeuron->h += (hiddenWeights[l][n][i] * hiddenNodes[l-1][i]->sigma);
				}
			}

			//add the current node's bias weight to the h value that was just calculated
			curNeuron->h += curNeuron->bias;
		
			//determine the current node's sigma value based on its h value
			curNeuron->sigma = calculateSigma(curNeuron->h);	
		}
	}

	//iterate through the output node(s) and compute the output values
	for(int o = 0; o < outputNodes.size(); o++){
		//for each output, iterate through nodes in last hidden layer and grab their information
		for(int n = 0; n < hiddenNodes[hiddenLayers - 1].size(); n++){
			outputNodes[o]->h += (outputWeights[o][n] * hiddenNodes[hiddenLayers - 1][n]->sigma);
		}

		//add the current node's bias weight to the h value that was just calculated
		outputNodes[o]->h += outputNodes[o]->bias;

		//determine the current node's sigma value based on its h value
		outputNodes[o]->sigma = calculateSigma(outputNodes[o]->h);
	}
}


//after each forward pass and backward pass for a given pattern, all of the weights are adjusted
void Network::updateWeights(vector <double> &inputs){
	Node *curNeuron;
	
	//iterate through all the nodes and adjust their weight values hiddenNodes[layer][node]
	for(int l = 0; l < hiddenNodes.size(); l++){
		for(int n = 0; n < hiddenNodes[l].size(); n++){
			curNeuron = hiddenNodes[l][n];
			
			//if this is the first hidden layer, look back to the inputs
			if(l == 0){
				//iterate through all of the inputs in the input layer
				for(int i = 0; i < inputCount; i++){
					//update the weight connections from the current node to each of the input values
					hiddenWeights[l][n][i] += (learnRate * curNeuron->delta * inputs[i]);
				}
			}
			//if not the first hidden layer, look back to the previous hidden layer
			else{
				//iterate through the nodes in the previous hidden layer
				for(int i = 0; i < hiddenNodes[l-1].size(); i++){
					//update the weight connections from the current node to each of the previous layer's nodes
					hiddenWeights[l][n][i] += (learnRate * curNeuron->delta * hiddenNodes[l-1][i]->sigma);
				}
			}

			//update the current node's bias value
			curNeuron->bias += (learnRate * curNeuron->delta);
		}
	}

	//iterate through the output node(s) and update the weight values
	for(int o = 0; o < outputNodes.size(); o++){
		//for each output, iterate through nodes in last hidden layer and grab their information
		for(int n = 0; n < hiddenNodes[hiddenLayers - 1].size(); n++){
			outputWeights[o][n] += (learnRate * outputNodes[o]->delta * hiddenNodes[hiddenLayers - 1][n]->sigma);
		}
		
		//update the bias weight of output node(s)
		outputNodes[o]->bias += (learnRate * outputNodes[o]->delta);
	}
}


//simple logistic sigmoid function
double calculateSigma(double h){
	return (1.0 / (1.0 + exp(-1.0*h)));
}


//for creating a new experiment data file (.csv) and adding some header information to the file
void newDataFile(string filename, string currentTest, int currentProblem){
	FILE *outFile;
	//string filename = "Results/Problem" + to_string(currentProblem) + "/" + currentTest + ".csv";
	
	outFile = fopen(filename.c_str(), "wb");
	if(outFile != NULL){
		fprintf(outFile, "Problem %d\nExperiment:, %s\n\n", currentProblem, currentTest.c_str());
	} else{
		fprintf(stderr, "Could not create the data file for test: %s\n", currentTest.c_str());
		exit(1);
	}

	fclose(outFile);
}
