#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdio>
#include <vector>
#include <time.h>

using namespace std;

typedef unsigned char uchar;

//for creating a .pgm image of a given CA system
void make_pgm(string filename);

//for randomly populating a 30x30 grid with either -1 or 1
void make_grid(vector <vector <int> > &board);


int main(int argc, char *argv[]){ 
	string ouputname = argv[1];
	vector <vector <int> > board;

	srand(time(NULL));
	//make_pgm(outputname);


	//populat the board with either 1 or -1
	make_grid(board);



	return 0;

}

//for creating a .pgm image of a given CA system
void make_pgm(string filename){

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

//for randomly populating a 30x30 grid with either -1 or 1
void make_grid(vector <vector <int> > &board){
	//resize the board to be 30x30
	board.resize(30);
	for(int i = 0; i < 30; i++) board[i].resize(30);
	
	//populate it
	int rand_num;
	for(int x = 0; x < 30; x++){
		for(int y = 0; y < 30; y++){
			if((rand() % 2) == 1) board[x][y] = 1;
			else board[x][y] = -1;
		}
	}
}
