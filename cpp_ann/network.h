#ifndef NETWORK_H_
#define NETWORK_H_


#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <iostream>

using namespace std;

class Network {
private:
	int neurons; //number of neurons
	int neighbors; //number of neighbors (network degree)
	double rewiring; //small-world rewiring probability
    vector< vector<int> > C; //Adjacency matrix
	vector< vector<double> > W; //Weigth matrix
	vector<double> TH; //Threshold_i
	vector<bool> V_t; //Network state in time t
	vector<bool> V_o; //Network state for pattern hebb learning
	double THETA_0; //value of Theta_0 for all patterns

public:
	//Constructors
	Network(int nN, int nK, double rP, int width, int height, char topology);
	//Functions for topology matrix generation
	void swRingGenerator(int); //Generates a Small-world Ring Topology Matrix
    void erSymGenerator(int); //Generates a Erdos-Renyi Topology Matrix
	void swSquareGridGenerator(int, int, int, int, int); //Generates a Square-Grid Topology Matrix
	void swXGridGenerator(int, int); //Generates a Small-world X-Grid Topology Matrix
	void swCrossGridGenerator(int, int); //Generates a Small-world Cross-Grid Topology Matrix
	void swRewiring(int); //Rewires the connectivity adjacency list of the input node with the network rewiring probability

    //Reads pattern from file
	void loadPatternFile(char *);

	//Sets network initial condition
	void networkInitialCodition();

	//Sets network initial condition with the inpput noise value
	void networkInitialCodition(double);

	//Performs HEBB learning rule
	void hebbLearning();

	/*
	Perform network time update
    for the given initial conditions and network parameters
	*/
	vector<double> updateNet(int time, int blocks, double sparseness, char th_fun,
        double th_value, int pat, const char * file_name, bool w_filename, int x_win, double rho);

	/*
	Performs the calculation of the overlap between a network state and a learned pattern
	*/
	vector<double> mdCalculate(int, double, vector<bool> &, vector<bool> &);

	/*
	Performs the calculation of the intra-overlap between a network state and a learned pattern
	for a given windows size
	*/
	vector<double> mdCalculateWin(int, double, vector<bool> &, vector<bool> &);

	/*
	Compares overlaps and activity values between the network state and the pattern
	for the last two updating time steps
	Stopping criteria: If the values are the same then network update stops
	if not
	network continue updating until reaching the specified maximum updating time
	*/
	bool mdComparison(vector<double> &, vector<double> &);

	/*
	Searches a node value in the adjacency list for a given node in order to
	avoid duplicated connections
	*/
	bool searchValue(int, int);
	bool searchValue(vector<int> sV, int valtoSearch);

	//Print network state to a text file
	void printNetworkState(const char *, int);

	//Print main network parameters and the topology adjacency list to screen
    void toString();

    //Random seed initialization
    void seed();

    // return a uniform number in [0,1].
	double unifRand();

	//Threshold functions
    double cutlinearFunction(double, double, double, double);
    double sinFunction(double, double);
    double stepFunction(double, double);
    double stepCutFunction(double, double, double);
    double rhoFunction(double, double, double, double, double);

    //Calculates a theta0 value given the activity of each pattern
    double thetaZero();

    void thetaZero(double);

    /*
    Calculates the mean activity of the input vector
    Use to calculate the activity of patterns and network states
    */
    double vectorMean(vector<bool> &);

    //Print the time evolution of the overlap to a text file
    void mdtimeEvolution(int, double, vector<double> &, const char *);

    //Print to screen the topology adjacency matrix to screen
    void toAdjMat();

    //Random pattern with given sparseness
    void randomPattern(double sparseness, const char *);
    void randomPattern(double sparseness);

    //Generates a random subset from patterns set
    vector<int> randPatternSet(int setSize, int totalPat);

};

/*
Network Constructor.
nN: number of nodes (neurons)
nK: number of neighbors per node
rP: rewiring probability (omega parameter)
width: pattern width
height: pattern height
Network topology. r: Ring, x: X-Grid, c: Cross-Grid, l: Circle-Grid
*/
Network::Network(int nN, int nK, double rP, int width, int height, char topology) {

    //Random seed initialization
	seed();

	//Setting network parameters;
    neurons=nN;
    neighbors=nK;
    rewiring=rP;

    //Reserving Topology(C) and Weight(W) vectors capacity
	//C.reserve(neurons);
	//W.reserve(neurons);
	//V_t.reserve(neurons); //Network state in time t
	//V_o.reserve(neurons); //Network state in hebb learning phase
	//TH.reserve(neurons); //Neurons' threshold

	for (int n = 0; n < neurons; n++) {
	    //printf("Generating network......%d\r", n);
	    //C[n].reserve(neighbors);
	    //W[n].reserve(neighbors);
	    switch(topology) {
            case 'r':
                swRingGenerator(n);
            break;
            case 'x':
                swXGridGenerator(n, width);
            break;
            case 'c':
                swCrossGridGenerator(n, width);
            break;
            case 's':
                erSymGenerator(n);
            break;
            default:
            break;
        }

        V_o.push_back(0);
        V_t.push_back(0);
        TH.push_back(0);

	}

	if (topology == 'l') {
	    int lSide = (sqrt(neighbors+1)-1)/2;
        for (int i=0; i < width; i++) {
            for (int j=0; j < height; j++) {
                swSquareGridGenerator(i, j, width, height, lSide);
            }
        }
	}

}

//Builds small world SquareGrid topology matrix C
void Network::swSquareGridGenerator(int i, int j, int width, int height, int lSide) {

     vector<int> tmpC;
     vector<double> tmpW;

     int ni = i * width + j;

     int iI=i-lSide;
     int iJ=j-lSide;

     int fI=i+lSide;
     int fJ=j+lSide;

     //upper limit

     if (iI < 0) {
          iI = 0; fI = iI + 2*lSide;
     }

     if (iJ < 0) {
          iJ = 0; fJ = iJ + 2*lSide;
     }

     //bottom limit

     if (fI >= height) {
          fI = height-1; iI = fI - 2*lSide;
     }

     if (fJ >= width) {
          fJ = width-1; iJ = fJ - 2*lSide;
     }

     //Builiding connectivity matrix C
     int kC = 0; //K counter
     for (int ip=iI; ip <= fI; ip++) {
          for (int jp=iJ; jp <= fJ; jp++) {
               int kip = ip * width + jp;
                if (kip != ni) {
                     tmpC.push_back(kip);
                     tmpW.push_back(0.0);
                     kC++;
               }
          }
     }

	C.push_back(tmpC);
	W.push_back(tmpW);

	//Rewiring connectivity
	if (rewiring > 0.0) {
	    swRewiring(ni);
	}

}

//Builds erdos-renyi topology matrix C
void Network::erSymGenerator(int i) {
    vector<int> tmpC;
    vector<double> tmpW;
    int kn = neighbors;

    //Builiding connectivity matrix C
    for (int j = 0; j < kn; j++)
    {
        //Generating Right neighbors
        int rn = i + (j + 1);
        if (rn < neurons) {
            tmpC.push_back(rn);
        }
        else {
            tmpC.push_back(rn - neurons);
        }
        tmpW.push_back(0.0);

        //Generating Left neighbors
        int ln = i - (j+1);
        if (ln < 0) {
            tmpC.push_back(neurons + ln);
        }
        else {
            tmpC.push_back(ln);
        }
        tmpW.push_back(0.0);

    }

    C.push_back(tmpC);
    W.push_back(tmpW);

}

//Builds small world RING topology matrix C
void Network::swRingGenerator(int i) {
    vector<int> tmpC;
    vector<double> tmpW;
	int kn = neighbors/2;

    //Builiding connectivity matrix C
	for (int j = 0; j < kn; j++)
	{
		//Generating Right neighbors
		int rn = i + (j + 1);
		if (rn < neurons) {
			tmpC.push_back(rn);
		}
		else {
			tmpC.push_back(rn - neurons);
		}
		tmpW.push_back(0.0);

		//Generating Left neighbors
		int ln = i - (j+1);
		if (ln < 0) {
			tmpC.push_back(neurons + ln);
		}
		else {
			tmpC.push_back(ln);
		}
		tmpW.push_back(0.0);

	}

	C.push_back(tmpC);
	W.push_back(tmpW);

	//Rewiring connectivity
	if (rewiring > 0.0) {
	    swRewiring(i);
	}

}

//Builds small world XGRID topology matrix C
void Network::swXGridGenerator(int i, int width) {
    vector<int> tmpC;
    vector<double> tmpW;
	int kn = neighbors/4;

	//Builiding connectivity matrix C
	for (int j = 0; j < kn; j++)
	{
		//Generating Right-Down Neighbors
		int RDN = i + (width + 1)*(j+1);
		if (RDN < neurons) {
			tmpC.push_back(RDN);
		}
		else {
			tmpC.push_back(RDN - neurons);
		}
		tmpW.push_back(0.0);

		//Generating Left-Down Neighbors
		int LDN = i + (width - 1)*(j+1);
		if (LDN < neurons) {
			tmpC.push_back(LDN);
		}
		else {
		    tmpC.push_back(LDN - neurons);
        }
		tmpW.push_back(0.0);

		//Generating Right-Up Neighbors
		int RUN = i - (width - 1)*(j+1);
		if (RUN < 0) {
			tmpC.push_back(RUN + neurons);
		}
		else {
			tmpC.push_back(RUN);
		}
		tmpW.push_back(0.0);

		//Generating Left-Up Neighbors
		int LUN = i - (width + 1)*(j+1);
		if (LUN < 0) {
			tmpC.push_back(LUN + neurons);
		}
		else {
		    tmpC.push_back(LUN);
        }
		tmpW.push_back(0.0);

	}

    C.push_back(tmpC);
	W.push_back(tmpW);

	//Rewiring connectivity
	if (rewiring > 0.0) {
	    swRewiring(i);
	}

}

//Builds small world CROSS GRID neighborhood
void Network::swCrossGridGenerator(int i, int width) {
    vector<int> tmpC;
    vector<double> tmpW;
	int kn = neighbors/4;

	//Builiding connectivity matrix C
	for (int j = 0; j < kn; j++)
	{
		//Generating Right neighbors
		int rn = i + (j + 1);
		if (rn < neurons) {
			tmpC.push_back(rn);
		}
		else {
			tmpC.push_back(rn - neurons);
		}
		tmpW.push_back(0.0);

		//Generating Left neighbors
		int ln = i - (j+1);
		if (ln < 0) {
			tmpC.push_back(neurons + ln);
		}
		else {
			tmpC.push_back(ln);
		}
		tmpW.push_back(0.0);
	}

	for (int j = 0; j < kn; j++)
	{
		//Generating down neighbors
		int dn = i + (j+1)*width;
		if (dn < neurons) {
			tmpC.push_back(dn);
		}
		else {
			tmpC.push_back(dn - neurons);
		}
		tmpW.push_back(0.0);

		//Generating upper neighbors
		int un = i - (j+1)*width;
		if (un < 0) {
			tmpC.push_back(neurons + un);
		}
		else {
			tmpC.push_back(un);
		}
		tmpW.push_back(0.0);

	}

	C.push_back(tmpC);
	W.push_back(tmpW);

    //Rewiring connectivity
	if (rewiring > 0.0) {
	    swRewiring(i);
	}

}

//Rewires the connectivity adjacency list of the input node with the network rewiring probability
void Network::swRewiring(int i) {

    for (int j = 0; j < neighbors; j++) {
		double rg = unifRand(); //generates random value between 0 and 1
		bool found = true;
		if (rg  < rewiring) {
			do {
                int newNode = (int)floor(unifRand()*neurons); //Generates a new random node
                found = searchValue(newNode, i); //Test if new node is already a neighbor
                //Test for repeated nodes, self-connection, and if node is out of neuron's range
                if (found == false && newNode != i && newNode >= 0 && newNode < neurons)
                    C[i][j] = newNode; //Assings new random node to neighborhood
			} while (found == true);
		}
	}

}

//Reads pattern from file
void Network::loadPatternFile(char * file) {
	FILE * pFile;
  	pFile = fopen (file,"r");

  	char str[12];

  	for (int ni = 0; ni < neurons; ni++) {

  		fscanf(pFile, "%s", str);

        V_o[ni] = (short int)atoi(str);

  	}
  	fclose(pFile);
}

void Network::randomPattern(double sparseness, const char * file_out) {
    FILE * oFile = fopen (file_out,"w");

	for (int i = 0; i < neurons; i++) {
	    double rg = unifRand(); //generates random value between 0 and 1
	    if (rg < sparseness)  {
            //V_o[i] = 1;
            fprintf(oFile, "%d ",1);
        }
        else {
            //V_o[i] = 0;
            fprintf(oFile, "%d ",0);
        }
	}

	fclose(oFile);
}

void Network::randomPattern(double sparseness) {

	for (int i = 0; i < neurons; i++) {
	    double rg = unifRand(); //generates random value between 0 and 1
	    if (rg < sparseness)  {
            V_o[i] = 1;
        }
        else {
            V_o[i] = 0;
        }
	}
}

vector<int> Network::randPatternSet(int setS, int totalPat) {

    vector<int> pS;

    for (int pi = 0; pi < setS; pi++) {

        //int gP = ceil(unifRand()*totalPat); //cout << gP << " ";

		bool found = true;
		do {
            int gP = (int)floor(unifRand()*totalPat); //Generates a new random node
            found = searchValue(pS, gP); //Test if new node is already a neighbor
            //Test for repeated nodes, self-connection, and if node is out of neuron's range
            if (found == false && gP >= 1 && gP <= totalPat)
                pS.push_back(gP); //Assings new random node to neighborhood
        } while (found == true);

        //pS.push_back(gP);

	}

	//for (int i=0; i < pS.size(); i++)
    //    cout << pS[i] << " ";


    //cout << endl;

    return pS;

}

//Sets network initial condition
void Network::networkInitialCodition() {
    for (int i = 0; i < neurons; i++) {
        V_t[i] = V_o[i];
    }
}

//Sets network noisy initial condition with the input noise
void Network::networkInitialCodition(double noise) {
    double V_o_act = vectorMean(V_o);
    double rg = 0, rg1 = 0;
    for (int i = 0; i < neurons; i++) {
        rg = unifRand();
        if (rg < noise) {
            rg1 = unifRand();
            if (rg1 < V_o_act) {
                V_t[i] = 1;
            }
            else {
                V_t[i] = 0;
            }
        }
        else
        {
            V_t[i] = V_o[i];
        }
    }
}

//Performs hebb learning
void Network::hebbLearning() {
    double V_o_act = vectorMean(V_o); //Gets pattern global activtiy
    double W_std_factor = V_o_act * (1 - V_o_act); //Gets activity variance
	float tmphebb;
	for (int n = 0; n < neurons; n++)
	{
	    tmphebb = 0.0;
		for (int k = 0; k < neighbors; k++) {
		    //Performs hebb learning of pattern stored in V_o
			tmphebb = (V_o[n] - V_o_act) * (V_o[C[n][k]] - V_o_act) / ( W_std_factor );
			W[n][k] += tmphebb; //Update weight matrix
		}
	}
}

//Network update for every time step
vector<double> Network::updateNet(int s_time, int blocks, double sparseness, char th_fun,
    double th_value, int pat, const char * file_name, bool w_filename, int x_win, double rho) {

    //Initializing output vector
	vector<double> net_var;
	net_var.push_back(0); //value of m (overlap between pattern and network state)
	net_var.push_back(0); //value of d (block overlap)
 	net_var.push_back(0); //value of q_m (gloval activity)
	net_var.push_back(0); //value of q_d (block activity)
	net_var.push_back(0); //value of th_m (mean threshold value)
	net_var.push_back(0); //value of th_d (threshold std dev)

    /*Stores state of the network in t-1 step for parallel updating
    All neurons update their activity states simultaneously at discrete time steps
    The previous state in t-1 need to be stored to calculate the actual t state
    */
	vector<bool> V_tp;
    V_tp.reserve(neurons);

    //File variables values to store variables values at every time step to a text file
    char file_time[256];
    char file_win[256];

    //w_file = true: network update results are printed for every time step
    if ( w_filename == true) {
        //File name to store m (gobal overlap), d (block overlap) evolution in time
        strcpy(file_time, "mdqintime_");
        strcat(file_time, file_name);
        FILE * tFile = fopen (file_time,"w");
        fclose(tFile);

        //File name to store overlaps for the specified windows size for every time step
        strcpy(file_win, "xmi_");
        strcat(file_win, file_name);
        FILE * winFile = fopen (file_win,"w");
        fclose(winFile);

    }

    /*
    Calculating theta0 value for each pattern
    Comment next lines in order tu use a specified fixed threshold value
    not the calculated by thetaZero() fuction
    as well as the given sparseness as parameter
    */
    //th_value = thetaZero();
    //sparseness = vectorMean(V_o);

    //calculating slope for linear threshold function
    double slope = ((-2)*th_value) / (1 - 2 * sparseness);

    //Loop updates network for every time step
	for (int t = 0; t < s_time; t++) {
		//printf("Updating.................%d\r", t);

        //Network global activity
		double global_activity = 0.0;

        /*
        Storing previous state of the network
        and calculating global activiy of the network for every time step
		*/
		for (int n = 0; n < neurons; n++) {
		    global_activity += V_t[n];
            V_tp[n] = V_t[n];
		}

		global_activity /= neurons;

        //Calculates the percentage of bits changing every time step
		int hamm_dist = 0;

		/*
		Calculates overlap between initial network state and pattern
		If the initial overlap is lesser than 0.5 the rho1=1.0
		the threshold is more fleixible reducining the artifact
		of no retrieval for low alpha
		*/
		double rho1 = rho;
		//vector<double> net_var_0 = mdCalculateWin(blocks, sparseness, V_o, V_tp);
		if (t < 20)
            rho1 = 1.0/rho;
        else
            rho1 = rho;

        //Updating network node states
		for (int n = 0; n < neurons; n++) {

			double neural_field = 0.0; //Neural field calculated for each node n
			double local_activity = 0.0; //Local activity of node n neighborhood
			double varA; //Variance of local_activity

            //Calculating local activity of the k-neighbors of node n
            for (int k = 0; k < neighbors; k++) {
                local_activity += V_tp[C[n][k]];
            }

            local_activity /= neighbors;

            //Avoids division by zero when patterns are very sparse
            if (local_activity != 0.0) {

                varA = local_activity*(1.0-local_activity); //Variance of local_activity

                varA=sqrt(varA); //Std dev of local_activity

                //Calculating neural field of node n
                for (int k = 0; k < neighbors; k++) {

                    neural_field += W[n][k] * (V_tp[C[n][k]] - local_activity);

                }

                neural_field /= varA;

                neural_field /= neighbors;

                /*
                Calculates dynamic threshold TH[n] for every node at every time step
                */
                switch(th_fun) {
                    //linear threshold function
                    case 'l':
                    {
                        TH[n] = cutlinearFunction(sparseness, th_value, local_activity, slope);
                    }
                    break;
                    //rho threshold function
                    case 'r':
                    {
                        TH[n] = rhoFunction(sparseness, global_activity, local_activity, th_value, rho1);
                    }
                    break;
                    //step threshold function
                    case 's':
                    {
                        TH[n] = stepFunction(local_activity, th_value);
                    }
                    break;
                    //sine threshold function
                    case 't':
                    {
                        double th_v1 = th_value/rho1;
                        TH[n] = sinFunction(local_activity, th_v1);
                    }
                    break;
                    //Step-cut threshold function
                    case 'c':
                    {
                        TH[n] = stepCutFunction(local_activity, sparseness, th_value);
                    }
                    break;
                }

                neural_field -= TH[n];

                //Updating each node state V_t[n] at time t
                if (neural_field >= 0) {
                    V_t[n] = 1;
                }
                else {
                    V_t[n] = 0;
                }

			}

            //Calculating if state changed for hamming distance count
			if (V_t[n] != V_tp[n])
                hamm_dist++;

		}

        //Calculating overlap between net state and pattern for time t
	    vector<double> net_var_t = mdCalculate(blocks, sparseness, V_o, V_tp);

        //Comparing t network state with state at t-1 to test stop criterion
	    bool md_eq = mdComparison(net_var, net_var_t);

        //Calculating and printing results for every time step
	    if ( w_filename == true ) {

	        double h_d = hamm_dist*(100.0/neurons);

            if (t==0)
                puts("time, m, dm, q, dq, th, dth, hamming(%)");

	    	printf("%d, %f, %f, %f, %f, %f, %f %f\n", t, net_var_t[0],
                net_var_t[1], net_var_t[2], net_var_t[3], net_var_t[4], net_var_t[5], h_d);

            vector<double> win_output = mdCalculateWin(x_win, sparseness, V_o, V_tp);

            FILE * winFile = fopen (file_win,"a");

            //Printing mesoscopic results for the specified mesoscopic blocks in x_win
            for (int bi=0; bi < x_win; bi++) {
                fprintf(winFile, "%f, ", win_output[bi]);
            }
            fprintf(winFile, "\n");

            fclose(winFile);

            mdtimeEvolution(t, h_d, net_var_t, file_time); //Print results for each time t

	    }

        //if stop criterion is true update macroscopic results (m,d) and finish updating
		if (md_eq == true) {
	        net_var = net_var_t;
            net_var.push_back(t);
            break;
	    }
        else {
            net_var = net_var_t; //Save t-1 variable values to compare with t variable values
        }

        //If maximum time is reached store macroscopic results (m,d) and finish updating
		if (t == s_time-1) {

		    vector<double> net_var_t = mdCalculate(blocks, sparseness, V_o, V_tp);

		    net_var = net_var_t;
            net_var.push_back(t);

		}

	}

    //return macroscopic results (m,d)
    return net_var;

}

//Overlap calculation for mesoscopic blocks
vector<double> Network::mdCalculateWin(int bn, double sparseness, vector<bool> & V_in1, vector<bool> & V_in2) {

    int splitcut = neurons/bn; //Calculates block size

    vector<double> overlap_b(bn); //mesoscopic overlaps vector

    vector<double> q_b(bn); //mesoscopic pattern activity vector

    vector<double> q_net(bn); //mesoscopic network activity vector

    double q_std_factor = (neurons/bn);

	//Calculating mesoscopic overlaps for each block
	for (int b = 0; b < bn; b++) {
	    q_b[b] = 0;
	    q_net[b] = 0;
	    overlap_b[b] = 0;
        for (int i = b*splitcut; i < ((b+1)*splitcut); i++) {
             q_b[b] += V_in1[i];
             q_net[b] += V_in2[i];
        }
        q_b[b] /= q_std_factor; //Pattern activity in block b
        q_net[b] /= q_std_factor; //Network activity in block b
        //th_b[b] /= th_std_factor;

        for (int i = b*splitcut; i < ((b+1)*splitcut); i++) {
            overlap_b[b] += (V_in1[i]-q_b[b])*(V_in2[i]-q_net[b]);
        }

        overlap_b[b] /= ((neurons/bn)*(sqrt(q_b[b] * (1 - q_b[b])) * sqrt(q_net[b] * (1 - q_net[b])) ));

	}

	return overlap_b;

}

//Macroscopic overlap calculation
vector<double> Network::mdCalculate(int bn, double sparseness, vector<bool> & V_in1, vector<bool> & V_in2) {

    int splitcut = neurons/bn; //Calculates block size

    vector<double> overlap_b(bn); //mesoscopic overlaps vector

    vector<double> q_b(bn); //mesoscopic pattern activity vector

    vector<double> q_net(bn); //mesoscopic network activity vector

    vector<double> th_b(bn); //mesoscopic threshold vector

    double q_std_factor = (neurons/bn);
    double th_std_factor = (neurons/bn);

    //Calculating mesoscopic overlaps for each block
	for (int b = 0; b < bn; b++) {
	    q_b[b] = 0;
	    q_net[b] = 0;
	    th_b[b] = 0;
	    overlap_b[b] = 0;
        for (int i = b*splitcut; i < ((b+1)*splitcut); i++) {
             q_b[b] += V_in1[i];
             q_net[b] += V_in2[i];
             th_b[b] += TH[i];
        }

        q_b[b] /= q_std_factor;
        q_net[b] /= q_std_factor;
        th_b[b] /= th_std_factor;

        for (int i = b*splitcut; i < ((b+1)*splitcut); i++) {
            overlap_b[b] += (V_in1[i]-q_b[b])*(V_in2[i]-q_net[b]);
        }

        overlap_b[b] /= ((neurons/bn)*(sqrt(q_b[b] * (1 - q_b[b])) * sqrt(q_net[b] * (1 - q_net[b])) ));

	}

	double m = 0;
	double d_s = 0;

	double q_m = 0;
	double q_d_s = 0;

	double th_m = 0;
	double th_d_s = 0;

    //Calculating global overlaps and activity as mean and variance over the blocks
    for (int b=0; b < bn; b++) {
        m += overlap_b[b];
        d_s += pow(overlap_b[b], 2);

        q_m += q_net[b];
        q_d_s += pow(q_net[b],2);

        th_m += th_b[b];
        th_d_s += pow(th_b[b],2);
    }

    m /= bn;
    double d = sqrt(d_s/bn - pow(m, 2));

    q_m /= bn;
    double q_d = sqrt(q_d_s/bn - pow(q_m, 2));

    th_m /= bn;
    double th_d = sqrt(th_d_s/bn - pow(th_m, 2));

    vector<double> o_v;

    o_v.push_back(m);
    o_v.push_back(d);
    o_v.push_back(q_m);
    o_v.push_back(q_d);
    o_v.push_back(th_m);
    o_v.push_back(th_d);

    return o_v;

}

/*
Calculates the mean activity of the input vector
Use to calculate the activity of patterns and network states
*/
double Network::vectorMean(vector<bool> & V_in) {
    double sum_vec = 0.0;
    for (int ni = 0; ni < neurons; ni++)
        sum_vec += V_in[ni];
    return sum_vec / neurons;
}

/*
Compares overlaps and activity values between the network state and the pattern
for the last two updating time steps
Stopping criteria: If the values are the same then network update stops
if not
network continue updating until reaching the specified maximum updating time
*/
bool Network::mdComparison(vector<double> & v1, vector<double> & v2) {
    bool md_equality = false;
    int counter = 0;
    for (int i = 0; i < 4; i++) {
        if (v1[i] == v2[i])
            counter++;
    }
    if (counter == 4)
        md_equality = true;

    return md_equality;
}

// Reset the random number generator with the system clock.
void Network::seed()
{
    srand(time(NULL));
}

// return a uniform number in [0,1].
double Network::unifRand()
{
    return rand() / double(RAND_MAX + 1.0);
}

/*
Searches a node value in the adjacency list for a given node in order to
avoid duplicated connections
*/
bool Network::searchValue(int cij_value, int node) {
	bool found = false;
	for (int j = 0; j < neighbors; j++) {
		if (C[node][j] == cij_value) {
			found = true;
			break;
		}
	}
	return found;
}

bool Network::searchValue(vector<int> sV, int valtoSearch) {
	bool found = false;
	for (int j = 0; j < sV.size(); j++) {
		if (sV[j] == valtoSearch) {
			found = true;
			break;
		}
	}
	return found;
}

//Print network state to a text file
void Network::printNetworkState(const char * file, int width) {
    FILE * vsFile;
    vsFile = fopen (file,"w");

    //fprintf(vsFile, "%d",(int)V_tp[0]);

    for (int i = 0; i < neurons; i++) {
		fprintf(vsFile, "%d ",(int)V_t[i]);
		if ((i+1) % width == 0) {
		    fprintf(vsFile, "\n");
		}
	}

	fclose(vsFile);
}

//Cut linear threshold function
double Network::cutlinearFunction(double sparseness, double th_value, double local_act, double slope) {
    if (local_act < sparseness)
        return 0.0;
    else if (local_act > (1-sparseness) )
        return 0.0;
    else
        return slope*(local_act - sparseness) + th_value;
}

//Sine threshold function
double Network::sinFunction(double local_act, double pos_th) {
    if (local_act < 0.5)
        return pos_th*sin(2*M_PI*local_act);
    else
        return pos_th*sin(2*M_PI*local_act);
}

//Step threshold function
double Network::stepFunction(double local_act, double pos_th) {
    if (local_act == 0)
        return 0.0;
    else if (local_act < 0.5)
        return pos_th;
    else
        return pos_th*(-1.0);
}

//Rho-step threshold function
double Network::rhoFunction(double sparseness, double global_act, double local_act, double pos_th, double rhoTerm) {
    //double rhoTerm = 0.7;
    double oth = 0.0;
    if (global_act > (sparseness + 0.5) / 2  ) {
        oth = rhoTerm * pos_th;
    }
    else {
        oth = 1 / rhoTerm * pos_th;
    }

    if (local_act == 0)
        return 0.0;
    else if (local_act < 0.5)
        return oth;
    else
        return oth*(-1.0);
}

//Step-cut threshold function
double Network::stepCutFunction(double local_act, double sparseness, double pos_th) {
    if (local_act < sparseness)
        return 0.0;
    else if (local_act > (1-sparseness) )
        return 0.0;
    else if (local_act == 0)
        return 0.0;
    else if (local_act < 0.5)
        return pos_th;
    else
        return pos_th*(-1.0);
}

//Calculates a theta0 value given the activity (sparseness) of each pattern
double Network::thetaZero() {
    double theta0 = 0.0;
    double V_o_act = vectorMean(V_o);
    double W_std_factor = V_o_act * (1 - V_o_act);
    theta0 = (1-2*V_o_act)/(2*sqrt(W_std_factor));
    return theta0;
}

//Calculates a theta0 value given the activity (sparseness) of each pattern
void Network::thetaZero(double sparseness) {
    double std_factor = sparseness * (1 - sparseness);
    THETA_0 = (1-2*sparseness)/(2*sqrt(std_factor));
}

//Print the time evolution of the overlap to a text file
void Network::mdtimeEvolution(int t, double hamm_dist, vector<double> & V_in, const char * file_time) {
    FILE * tFile = fopen (file_time,"a");
    fprintf(tFile, "%d, %f, %f, %f, %f, %f, %f, %f\n", t, V_in[0],
        V_in[1], V_in[2], V_in[3], V_in[4], V_in[5],hamm_dist);
    fclose(tFile);
}

//Print main network parameters and the topology adjacency list to screen
void Network::toString() {

    cout << "N=" << neurons << ", K=" << neighbors
        << ", w=" << rewiring << endl;
    for (int i = 0; i < neurons; i++)
    {
        for (int j = 0; j < neighbors; j++) {
            //cout << "\"" << i << "\"" << "->" << "\"" << C[i][j] << "\"" << ", ";
            cout << i << "->" << C[i][j] << ", ";
        }
        cout << endl;
    }

}

//Print to screen the topology adjacency matrix to screen
void Network::toAdjMat() {

    for (int i = 0; i < neurons; i++)
    {
        for (int j = 0; j < neurons; j++) {

            bool fV=searchValue(i, j);
            if (fV == true)
                cout << 1 << " ";
            else
                cout << 0 << " ";
        }
        cout << endl;
    }

}

#endif /*NETWORK_H_*/
