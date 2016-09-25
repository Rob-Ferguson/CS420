#include "networkdef.h"

void Network::trainNet(){
	/* For each set inputs from the training data file, initialize the network to use those inputs
		and then perform a forward pass, backward pass, and update the network weights */
	for(int i = 0; i < trainingData.size(); i++){
		//compute the outputs of all nodes within each successive layer of the network
		forwardPass(trainingData[i]);
		
		//traverse through layers in the network backwards and calculate delta values for each node
		backwardPass(trainingData[i]);

		//after all values have been calculated, update weight values between each pair of nodes
		updateWeights(trainingData[i]);
	}
}


void Network::readFileData(int problem){
    int i, j;
    ifstream file;
    string line;                                        //reads a line of input at a time
    stringstream buffer;                                //used for parsing the input line
    vector <double> lineValues(inputCount+outputCount);  //temporarily stores the values from each line of input to add to the main vector representing each input file 
    
    lineValues.resize(inputCount+outputCount);
    
    //opening the trainingFile
    file.open(trainingFile.c_str());
    if (!(file.is_open())){ 
        fprintf(stderr, "Training File did not open.\n");
        exit(1);
    }

    //reads every line in the training file
    while (getline(file, line)){
        
        //clears the error state and contents of the buffer
        buffer.clear();
        buffer.str("");

        //parses each line and stores the values in the trainingData vector
        buffer.str(line);    
        for(i = 0; i < inputCount+outputCount; i++){
            buffer >> lineValues[i];
        }
        trainingData.push_back(lineValues);
            
    }
    
    //closes the file so the next one can be opened
    file.close();
    
    //opening the validatingFile
    file.open(validationFile.c_str());
    if (!(file.is_open())){ 
        fprintf(stderr, "Validation File did not open.\n");
        exit(1);
    }

    //reads every line in the validating file
    while (getline(file, line)){
        
        //clears the error state and contents of the buffer
        buffer.clear();
        buffer.str("");

        //parses each line and stores the values in the validatingData vector
        buffer.str(line);
        for (i = 0; i < inputCount+outputCount; i++){
            buffer >> lineValues[i];
        }
        validationData.push_back(lineValues);

    }
    
    //closes the file so the next one can be opened
    file.close();

    //opening the testingFile
    file.open(testingFile.c_str());
    if (!(file.is_open())){ 
        fprintf(stderr, "Testing File did not open.\n");
        exit(1);
    }

    //reads every line in the testing file
    while (getline(file, line)){

        //clears the error state and contents of the buffer
        buffer.clear();
        buffer.str("");

        //parses each line and stores the values in the testingData vector
        buffer.str(line);
        for (i = 0; i < inputCount+outputCount; i++){
            buffer >> lineValues[i];
        }
        testingData.push_back(lineValues);

    }   
  
    //closes the file when all of the data has been read in
    file.close();
}

void Network::backwardPass(vector <double> &inputs){
    int l, n, k;
    double deltaoutput;
    double deltahidden;
    double sigma;
    vector <double> expectedoutputs(outputCount);
    
    //pushes back all of the outputs into an output vector
    for (l = inputCount; l < inputCount+outputCount; l++){
        expectedoutputs.push_back(inputs[l]);
    }
    
    //calculate the delta values for the output layer
    for (l = 0; l < outputNodes.size(); l++){
        sigma = outputNodes[l]->sigma;
        outputNodes[l]->delta = (2.0 * sigma * (1 - sigma) * (inputs[inputCount + l] - sigma));
    }

    //iterates backward through all of the nodes in the hidden layers
    for (l = hiddenLayers-1; l >= 0; l--){

        //iterates through individual nodes in the specific layer
        for (n = 0; n < hiddenNodes[l].size(); n++){
        

            //case for if currently in the last hidden layer
            if (l == hiddenLayers-1){
               
                for (k = 0; k < outputCount; k++){
                    sigma = hiddenNodes[l][n]->sigma;
                    hiddenNodes[l][n]->delta += (sigma * (1 - sigma) * outputNodes[k]->delta * outputWeights[k][n]); 
                }
            }
            else{

                for (k = 0; k < hiddenNodes[l+1].size(); k++){
                    sigma = hiddenNodes[l][n]->sigma;
                    hiddenNodes[l][n]->delta += (sigma * (1 - sigma) * hiddenNodes[l+1][k]->delta * hiddenWeights[l+1][k][n]);
                }
            }
        }
    }
}


void Network::printCSV(){

    int i, j;
    
    //opens the file for writing
    ofstream output(filename.c_str(), ios_base::app);
    
    output << "Network Number:" << "," << networkNumber << endl;
    output << "Hidden Layers:" << "," << hiddenLayers << endl;
    
    //prints the hidden layer numbers
    output << "Layer Number:" << ",";
    for (i = 0; i < hiddenLayers; i++){
        if (i != hiddenLayers-1){
            output << i << ",";
        }
        else{
            output << i << endl;
        }
    }
    
    //prints the neurons for a given layer
    output << "Neurons/Layer:" << ",";
    for (i = 0; i < hiddenLayers; i++){
        if (i != hiddenLayers-1){
            output << neurons[i] << ",";
        }
        else{
            output << neurons[i] << endl;
        }

    }
    
    output << "Learning Rate:" << "," << learnRate << endl;
    output << "Total Epochs:" << "," << epochs << endl;

    //prints the testing RMSE value
    output << "Testing RMSE:" << ",";
    for (i = 0; i < rmseTesting.size(); i++){
        output << rmseTesting[i] << endl;
    }
    
    output << "Epoch Number:" << "," << "Validation RMSE / Epoch:" << endl;
    for (i = 0; i < epochs; i++){
        output << i << ",";
        for (j = 0; j < rmseValidation.size(); j++){
            if (j != rmseValidation.size()-1){
                output << rmseValidation[j][i] << ",";
            }
            else{
                output << rmseValidation[j][i] << endl;
            }
        }
    }
    output << endl << endl;




    
    //prints the epoch numbers
   // output << "Epoch Number:" << ",";
    //for (i = 0; i < epochs; i++){
    //    if (i != epochs-1){
    //        output << i << ",";
    //    }
    //    else{
    //        output << i << endl;
    //    }

    //}

    //prints the validation RMSE values
    //output << "Validation RMSE / Epoch:" << ",";
    //for (i = 0; i < rmseValidation.size(); i++){

    //    for (j = 0; j < rmseValidation[i].size(); j++){
            
      //      if (j != rmseValidation[i].size()-1){
       //         output << rmseValidation[i][j] << ",";
        //    }
          //  else{
     //           output << rmseValidation[i][j];
          //  }

       // }
    //    output << endl;
   // }
    
    output << endl;

    output.close();

}
