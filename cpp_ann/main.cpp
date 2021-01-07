#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string.h>
#include "network.h"

using namespace std;

bool checkParameters(int); //Checks for parameter number
string returnFileName(char *argv[]); //Generates output file name
string returnFilePattern(int bS, int sS, char * ruta); //Returns input file pattern
string returnOutFile(int bS, int sS); //Returns output file pattern

int main(int argc, char *argv[])
{
	//Checking input parameters count
	bool continuar = checkParameters(argc);

	if (continuar == true) {

		//Generating output file name
		char file_out[256];
		strcpy(file_out, returnFileName(argv).c_str());

		//Assigning parameter values
	    int Neurons = atoi(argv[1]); //Number of nodes in the network
        int Degree = atoi(argv[2]); //K-neighbors per node
        double rewProb = atof(argv[3]); //Rewiring Probability: Omega parameter
        double sparseness = atof(argv[4]); //Sparseness (activity level) of the learning patterns
        int blocks = atoi(argv[5]); //Number of blocks in the patterns
        char th_fun = *argv[6]; //Threshold function
        double th_value = atof(argv[7]); //Threshold value: theta0
        double rho = atof(argv[8]); //Value of rho
        double np = atof(argv[9]); //Noise applied to initial states m0=1-np
        int time = atoi(argv[10]); //Time steps for network evolution
        int iniPat = atoi(argv[11]); //Initial pattern
        int patterns = atoi(argv[12]); //Final pattern
        int pat_int = atoi(argv[13]); //pattern interval
        int x_win = atoi(argv[14]); //x_win points (mxi_t): measures mesoscopic parameters in x-size windows
        int width = atoi(argv[15]); //pattern width
        int height = atoi(argv[16]); //pattern height
        char topology = *argv[19]; //network topology
        int subsetSize = atoi(argv[20]); //subnet size (K_b)
        int nNets = atoi(argv[21]);  // number of subnets, nNets x subsetSize = patterns

        FILE * oFile = fopen (file_out,"w");
		fclose(oFile);

        for (int ni=0; ni<nNets; ni++) {

        	//Generating small-world network int *ptr; ptr=new int[size];
    		Network Net(Neurons, Degree, rewProb, width, height, topology); //cout << "Red bien"; cin.get();
    		/*
    		Uncomment next line to printscreen the network topology
            Notice that N=widthxheigt, i.e. Use: N=6x6=36, K=8, width=6, height=6
            */
            // Net.toString(); cin.get();

            //Uncomment next line to printscreen adjacency matrix
            //Net.toAdjMat(); cin.get();

            /*
            w_file = true prints simulation results for every time step
            w_file = false prints simulation results for last time step
            */
    		bool w_file = false;

            int p = 0; //Learned patterns counter

            //Loop the set of patterns for learning and retrieve testing.
            //Each step a new pattern is learned by the network,
            //and right after initializing in a state with np noise
            //the network evolves for the given time steps.
            for (int il=0;il<subsetSize;il++) {

                for (int iil=6;iil<=pat_int;iil++) {

                    //Generating filename for learning patterns in path1
                    char file_in[256];
                    //strcpy(file_in, returnFilePattern(pS[il], iil, argv[17]).c_str());
                    strcpy(file_in, returnFilePattern(il+1+ni*subsetSize, iil, argv[17]).c_str());

                    //Read learning pattern file
                    Net.loadPatternFile(file_in);

                    //Hebb learning of file_in pattern
                    Net.hebbLearning();

                }

            }

            // Retrieval test for patterns
            for (int ir=1;ir<=patterns;ir++) {

                for (int iir=6;iir<=pat_int;iir++) {

                    //Generating filename for intial network states in path2
                    char file_in0[256];
                    strcpy(file_in0, returnFilePattern(ir, iir, argv[18]).c_str());

                    //Read intial state pattern
                    Net.loadPatternFile(file_in0);

                    //Applies a network initial condition with np noise
                    Net.networkInitialCodition(np);

                    /*
                    Perform network time update
                    for the given initial conditions and network parameters
                    */
                    vector<double> output_values = Net.updateNet(time, blocks, sparseness, th_fun,
                        th_value, patterns, file_out, w_file, x_win, rho);


                    oFile = fopen (file_out,"a");

    		        fprintf(oFile,"%d, %f, %d\n", ir, 
                        output_values[0],
                        (int)output_values[6]);

                    fclose(oFile);

                    p++; //Increase patterns counter

                }
            }

        }

	}
	return 0;
}

//Checks for parameter number and return some instructions if number is not OK
bool checkParameters(int argc) {
	if (argc < 20)
    {
        printf("Usage: \n");
	    printf("./sparsenet N K w a B T t rho np time p P pi x wth ht path1 path2 top subsetSize nNets\n");
	    printf("N:           number of neurons (nodes) in the network, N=widthxheigt\n");
	    printf("K:           number of neighbors per node\n");
	    printf("w:           rewiring probability, omega parameter: w=0 regular ring or grid, w=1 random net\n");
	    printf("a:           sparseness\n");
	    printf("B:           number of blocks, use always B=1 for fingerprints\n");
	    printf("T:           threshold function: r > rho, l > linear, s > square, t > sine, c > squarecut\n");
	    printf("t:           threshold value\n");
	    printf("rho:         rho value\n");
	    printf("np:          Noise value applied to initial state\n");
	    printf("time:        max simulation time\n");
	    printf("p:           initial pattern (fingerprint in path1 and path2)\n");
	    printf("P:           final pattern (fingerprint in path1 and path2)\n");
	    printf("pi:          pattern interval (digit after _ in files in path1 and path2)\n");
	    printf("x:           x points (mxi_t): measures mesoscopic parameters in x-size windows (Not used in fp)\n");
	    printf("wth:         pattern width\n");
	    printf("ht:          pattern height\n");
	    printf("path1:       path to learning patterns folder (i.e. rolled fingerprints)\n");
	    printf("path2:       path to initial state patterns folder (i.e. latent fingerprints)\n");
	    printf("top:         network topology: s > er-sym, a > er-asym, r > ring, c > Cross-Grid, x > X-Grid, l > l-side SquareGrid\n");
        printf("subsetSize:  subset size for each module\n");
        printf("nNets:       number of modules, nNets x subsetSize = patterns\n");
	    printf("Examples: \n");
        printf("single: ./sparsenet 89420 240 0.5 0.2258 1 r 0.656 0.7 0.0 100 1 10 6 100 263 340 patterns/ patterns/ c 10 1\n");
        printf("ensemble: ./sparsenet 89420 24 1 0.2258 1 r 0.656 0.7 0.0 100 1 100 6 100 263 340 patterns/ patterns/ c 10 10\n\n");

	    return false;
	}
	else {
		return true;
	}
}

/*
Returns the filename for the output text file using the input parameters
*/
string returnFileName(char *argv[]) {

    char * pch;
    char * ruta;
    char str0[256];
    strcpy (str0,argv[18]);
    pch = strtok (str0,"/");
    while (pch != NULL)
    {
        ruta = pch;
        pch = strtok (NULL, "/");
    }

    ostringstream outputFile;

	outputFile << "N"
	<< argv[1]
    << "K" << argv[2]
    << "w" << argv[3]
    << "a" << argv[4]
    << "b" << argv[5]
    << "T" << argv[6]
    << "t" << argv[7]
    << "rho" << argv[8]
    << "np" << argv[9]
    << "time" << argv[10]
    << "p" << argv[11]
    << "P" << argv[12]
    << "pi" << argv[13]
    << "x" << argv[14]
    << "w" << argv[15]
    << "h" << argv[16]
    << "TY" << argv[19]
    << "SNS" << argv[20]
    << "NN" << argv[21]
    << "PW" << ruta << ".txt";

    return outputFile.str();

}

//Returns input pattern filename in the specified path (ruta)
string returnFilePattern(int bS, int sS, char * ruta) {

        ostringstream inputFile;

        inputFile << ruta;
        inputFile << bS << '_' << sS;

        return inputFile.str();

}

//Returns output pattern filename in local directory
string returnOutFile(int bS, int sS) {

        ostringstream inputFile;

        inputFile << bS << '_' << sS;

        return inputFile.str();

}
