// PID.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "mt19937ar.h"
#include "mymath.h"


//Evolutionary parameters
#define N 150					//Number of individuals
#define N_best (int)(N/10.0)	//Number of individuals in "hall of fame". Best individuals of all time
#define geneLength 24			//Number of bits per gene
#define nGenes 8				//Set this as the number of parameters in your model
int genomeLength = nGenes*geneLength; //Number of bits for all genes
int generations = 200; //Number of generations to run the simulation
double Tmin = 0.01; //Minimum temperature of the ladder. Preferrably between 0 and 1
double Tmax = 5; //Maximum temperature of the ladder. Preferrably around H_max
int seed = 8;//Seed for the random number generator

//Probabilities for exchange operator
double qm = 0.3; //The mutation probability vs crossover (1-qm for crossover)
double qma = 0.2; //smartCopy mutation rate (qm*(1-qma) for regular mutation)
double qmc = 0.2; //smartCopy crossover rate((1-qm)*(1-qmc) for regular crossover) 

//Probabilities for smart copy crossover
double p0 = 0.05;
double p1 = 0.15;
double p2 = 0.8;
double p_matrix[] = { p0*p0 + (1 - p0)*(1 - p0), 2 * p0*(1 - p0), p1*(1 - p2) + p2*(1 - p1), p1*p2 + (1 - p1)*(1 - p2) };

typedef struct {
	int genome[nGenes * geneLength]; //array 0's and 1's that constitute the genome. 3 parameters coded in 32 bits.
	double H;		//Fitness, or energy
	double parameter[nGenes];		//Proportional parameter	//Derivative parameter//Integral parameter
	double t; 		//temperature
	int stay_time;
	int gen; //generation
} personality;

typedef struct {
	//int n_params;
	double param_min[nGenes];
	double param_max[nGenes];
}properties;

typedef struct {
	int *loci;				//Array of points on a genome
	int store_counter;		
	int gen;				//generation
	double simulation_time;
	clock_t simulation_tic;
	clock_t simulation_toc;
	double evolution_time;
	clock_t evolution_tic;
	clock_t evolution_toc;
}arrays;

typedef struct datafiles {
	FILE *data;
	//FILE *surfdata;
	FILE *elitedata;
} datafiles;

double fitnessTest(double parameter[], int guynumber, arrays *arr) {
	double H;
	
	arr->simulation_tic = clock();//Start time counter


	//Write your fitness code here!!!!
	//H can be LSQ fitness or any other measure that needs to be minimized.
	//H should be  in the range [0, 30] approximately, because exp(H) is computed in the code. Large H leads to overflow.
	H = 3 + genrand_real1(); //Replace this line


	//I propose to use the function below to make H small overall, and sharp close to H = 0
	H = -1 / (log(H + 0.1) + 0.1) + log(H);

	//Record time spent
	arr->evolution_tic = arr->simulation_tic;
	arr->simulation_toc = clock();
	arr->simulation_time += (double)(arr->simulation_toc - arr->simulation_tic) / CLOCKS_PER_SEC;
	return H;
}
void setupParameters(properties *prop) {
	prop->param_min[0] = 0;
	prop->param_min[1] = 0;
	prop->param_min[2] = 0;
	prop->param_min[3] = 0;
	prop->param_min[4] = 0;
	prop->param_min[5] = 0;
	prop->param_min[6] = 0;
	prop->param_min[7] = 0;

	prop->param_max[0] = 200;
	prop->param_max[1] = 200;
	prop->param_max[2] = 200;
	prop->param_max[3] = 200;
	prop->param_max[4] = 200;
	prop->param_max[5] = 200;
	prop->param_max[6] = 200;
	prop->param_max[7] = 200;
}
void getFit(personality guys[], arrays *arr) {
	int i;
	for (i = 0; i < N; i++) {
		guys[i].H = fitnessTest(guys[i].parameter, i, arr);
		//printf("hej %.3f\n",guys[i].H);
	}
}
double bin2dec(int bin[], int binSize) {
	int i;
	//printf("%d\n",binSize);
	double dec = 0;
	for (i = 0; i < binSize; i++) {
		//printf("%d %d\n",i,bin[i]);
		if (bin[i] == 1) {
			dec += pow(2.0, binSize - i - 1);
		}
	}
	return dec;
}
double map2param(double param, int i, properties *prop) {
	double result;
	result = param / (pow(2.0, geneLength) - 1);
	result *= (prop->param_max[i] - prop->param_min[i]);
	result += prop->param_min[i];
	return result;
}
void getParameters4One(personality &guy, properties *prop) {
	int i, j;
	int gene[geneLength];
	//printf("%d\n",i);
	for (i = 0; i < nGenes; i++) {
		for (j = 0; j < geneLength; j++) {
			gene[j] = guy.genome[i*geneLength + j];
		}
		guy.parameter[i] = bin2dec(gene, sizeof(gene) / sizeof(gene[0]));
		guy.parameter[i] = map2param(guy.parameter[i], i, prop);
	}
}
void getParameters4All(personality guys[], properties *prop) {
	int i, j, k;
	int gene[geneLength];
	for (i = 0; i < N; i++) {
		for (j = 0; j < nGenes; j++) {
			for (k = 0; k < geneLength; k++) {
				gene[k] = guys[i].genome[j*geneLength + k];
			}
			guys[i].parameter[j] = bin2dec(gene, sizeof(gene) / sizeof(gene[0]));
			guys[i].parameter[j] = map2param(guys[i].parameter[j], j, prop);
		}
	}
}
void wakeUp(personality guys[], properties *prop) {
	int i, j;
	double T[N];
	//Initialize some temperature ladder, here logarithmic.
	logspace(T, Tmax, Tmin, N);
	//Initialize the genome sequence for all guys randomly
	for (i = 0; i < N; i++) {
		guys[i].t = T[i]; 		//Assign temperatures produced above
		guys[i].stay_time = 0;
		for (j = 0; j < genomeLength; j++) {
			guys[i].genome[j] = (int)nearbyint(genrand_real1());	//Distribute the genetic parameters
		}
	}
	getParameters4All(guys, prop);
}
void roulette_select(personality guys[], int selected[], double *Z, double s) {
	//Selected 0 is the guy with good fitness (lowest), selected 1 is random

	int i;
	//double total_H = 0;
	double roulette[N];
	double lucky_number;
	//Make a roulette wheel
	*Z = 0;
	for (i = 0; i < N; i++) {
		*Z += exp(-guys[i].H / s);
		roulette[i] = *Z; //Low fitness gives large area on roulette
	}
	lucky_number = *Z*genrand_real1();
	for (i = 0; i < N; i++) {
		if (lucky_number < roulette[i]) {
			selected[0] = i;
			break;
		}
	}
	selected[1] = selected[0];
	while (selected[1] == selected[0]) {
		selected[1] = (int)((N)*genrand_real2());
	}
}
void bitselector_smartCopy(personality guys[], personality newguys[], int selected[], int n_ab[]) {
	int i, j;
	for (i = 0; i < nGenes*geneLength; i++) {
		//If parents loci are the same
		if (guys[selected[0]].genome[i] == guys[selected[1]].genome[i]) {
			//Copy the values
			newguys[selected[0]].genome[i] = guys[selected[0]].genome[i];
			newguys[selected[1]].genome[i] = guys[selected[1]].genome[i];
			//Independently reverse with probability p0
			for (j = 0; j < 2; j++) {
				if (genrand_real1() < p0) {
					newguys[selected[j]].genome[i] = 1 - newguys[selected[j]].genome[i];
				}
			}
			if (newguys[selected[0]].genome[i] == newguys[selected[1]].genome[i]) {
				n_ab[0]++; //Kept the same}
			}
			else {
				n_ab[1]++;//Changed it
			}
		}
		//If parents loci are different
		else {
			//Copy the values
			newguys[selected[0]].genome[i] = guys[selected[0]].genome[i];
			newguys[selected[1]].genome[i] = guys[selected[1]].genome[i];
			//Reverse with probabilities p1 and p2 respectively
			if (genrand_real1() < p1) {
				newguys[selected[0]].genome[i] = 1 - newguys[selected[0]].genome[i];
			}
			if (genrand_real1() < p2) {
				newguys[selected[1]].genome[i] = 1 - newguys[selected[1]].genome[i];
			}
			if (newguys[selected[0]].genome[i] == newguys[selected[1]].genome[i]) {
				n_ab[2]++; //Kept the same}
			}
			else {
				n_ab[3]++;//Changed it
			}
		}
	}
}
void mutation(personality guys[], personality newguys[], arrays *arr, properties *prop) {
	int i, j;
	int mutantGenes;	//points to mutate
	int mutant;			//which guy to mutate
	for (i = 0; i < N; i++) {
		mutantGenes = (int)(genomeLength*genrand_real2());
		mutant = (int)(N*genrand_real2());
		rndChoice(arr->loci, mutantGenes, genomeLength);	//Choose which locus to mutate
		for (j = 0; j < mutantGenes; j++) {
			newguys[mutant].genome[arr->loci[j]] = 1 - newguys[mutant].genome[arr->loci[j]]; //reverse the gene(s)
		}
		getParameters4One(newguys[mutant], prop);
		newguys[mutant].H = fitnessTest(newguys[mutant].parameter, mutant, arr);
		//Perform metropolis
		double dH = newguys[mutant].H - guys[mutant].H;		//used to decide if we accept the new guy or not.
		j = mutant;
		if (dH < 0) {
			memcpy(&guys[j], &newguys[j], sizeof(guys[j]));
		}
		else if (exp(-dH / newguys[mutant].t) > genrand_real1()) {
			memcpy(&guys[j], &newguys[j], sizeof(guys[j]));
		}
		//Inform newguys of latest events... i.e sync them for the next round
		memcpy(&newguys[j], &guys[j], sizeof(newguys[j]));
	}
}
void mutation_elite(personality guys[], personality newguys[], personality bestguys[], arrays *arr, properties *prop) {
	int i,j;
	int mutantGenes;	//Number of points to mutate
	int mutant;			//which guy to mutate
	int elite_mutant;	//Which elite guy to receive wisdom from
	int count;
	int rand_int;
	for (i = 0; i < N; i++) {
		mutantGenes = 0;
		while (mutantGenes == 0) {
			mutantGenes = (int)(genomeLength*genrand_real2());
		}
		mutant = (int)(N*genrand_real2());
		elite_mutant = (int)(N_best*genrand_real2());
		//Fill loci with mutantGenes genome points to be mutated, from a triangular distribution for each parameter.
		count = 0;
		while (count < mutantGenes) {
			rand_int = tri_randint(0, geneLength) + (int)(genrand_real1()*geneLength*(nGenes - 1));
			arr->loci[count] = rand_int;
			count++;
		}
		
		//Copy an elite guy to a new guy to be a guinnea pig
		memcpy(&newguys[mutant].genome, &bestguys[elite_mutant].genome, sizeof(newguys[mutant].genome));

		for (j = 0; j < mutantGenes; j++) {
			newguys[mutant].genome[arr->loci[j]] = 1 - newguys[mutant].genome[arr->loci[j]]; //reverse the gene(s)
		}
		getParameters4One(newguys[mutant], prop);
		newguys[mutant].H = fitnessTest(newguys[mutant].parameter, mutant, arr);
		//Perform metropolis
		double dH = newguys[mutant].H - guys[mutant].H;		//used to decide if we accept the new guy or not.
		j = mutant; //Shorter notation...
		if (dH < 0) {
			memcpy(&guys[j], &newguys[j], sizeof(guys[j]));
			guys[j].stay_time = 0;
		}
		else if (exp(-dH / newguys[j].t) > genrand_real1()) {
			memcpy(&guys[j], &newguys[j], sizeof(guys[j]));
			guys[j].stay_time = 0;
		}
		//Inform newguys of latest events... i.e sync them for the next round
		memcpy(&newguys[j], &guys[j], sizeof(newguys[j]));
	}
}
void crossover(personality guys[], personality newguys[], arrays *arr, properties *prop) {
	//Start roulette selection
	//Should the fitnesses be ordered as H(xi) >=  H(xj)? Typ...
	int i, j, matings;
	//int 2 = 2;
	int nMatings = (int)(N / 5);
	int selected[2];
	int crossoverPoint;
	double s = guys[(int)(0.9*N)].t;
	//printf("s=%f\n",s);
	double rc, PXX, PYY, dHt0, dHt1;
	double expX[2];
	double expY[2];
	double Z;
	for (matings = 0; matings < nMatings; matings++) {
		roulette_select(guys, selected, &Z, s); //Selected 0 will be good, selected 1 random
												//expX[0] = exp(-guys[selected[0]].H / guys[selected[0]].t); //good guy
												//expX[1] = exp(-guys[selected[1]].H / guys[selected[1]].t); //bad guy
		expX[0] = exp(-guys[selected[0]].H / s); //good guy
		expX[1] = exp(-guys[selected[1]].H / s); //bad guy
		PXX = 1 / ((N - 1)*Z)*(expX[0] + expX[1]); //P((xi,xj) | x)
												   //Now mate the newguys to create offsprings
		crossoverPoint = (int)(genrand_real1()*genomeLength);
		for (i = crossoverPoint; i < genomeLength; i++) {
			newguys[selected[0]].genome[i] = guys[selected[1]].genome[i];
			newguys[selected[1]].genome[i] = guys[selected[0]].genome[i];
		}
		getParameters4One(newguys[selected[0]], prop);
		getParameters4One(newguys[selected[1]], prop);
		newguys[selected[0]].H = fitnessTest(newguys[selected[0]].parameter, selected[0], arr);
		newguys[selected[1]].H = fitnessTest(newguys[selected[1]].parameter, selected[1], arr);
		expY[0] = exp(-newguys[selected[0]].H / s);
		expY[1] = exp(-newguys[selected[1]].H / s);
		Z = 0;
		for (i = 0; i < N; i++) {
			Z += exp(-newguys[i].H / s);
		}
		PYY = 1 / ((N - 1)*Z)*(expY[0] + expY[1]); //P((xi,xj) | x)
		dHt0 = (newguys[selected[0]].H - guys[selected[0]].H) / guys[selected[0]].t; //good 
		dHt1 = (newguys[selected[1]].H - guys[selected[1]].H) / guys[selected[1]].t; //bad 

																					 //if (dHt1 + dHt2 < 0){ continue; }
		rc = exp(-dHt0 - dHt1)*PXX / PYY;
		//Accept or reject
		if (genrand_real1() < fmin(1, rc)) {
			for (i = 0; i < 2; i++) { //refresh newpos
				j = selected[i];
				memcpy(&guys[j], &newguys[j], sizeof(guys[j]));
			}
		}
		//Make sure newguys are up to speed on newest events
		for (i = 0; i < 2; i++) {
			j = selected[i];
			memcpy(&newguys[j], &guys[j], sizeof(newguys[j]));
		}
	}
}
void crossover_elite(personality guys[], personality newguys[], personality bestguys[], arrays *arr, properties *prop) {
	//Start roulette selection
	int i, j, matings;
	int nMatings = (int)(N_best / 5);
	int selected[2];
	int crossoverPoint;
	double rc, PXX, PYY, dHt1, dHt2;
	double expX[2];
	double expY[2];
	double Z;
	double s = guys[(int)(0.97*N)].t;
	int random_bestguy;
	for (matings = 0; matings < nMatings; matings++) {
		roulette_select(guys, selected, &Z, s);
		//Selected 0 is a boltzmann "just good" guy. Selected 1 is a random any guy.
		//Let a newguy selected 1 impersonate a random bestguy but keep the temperature. Then newguy gets superpowers
		random_bestguy = (int)(N_best*genrand_real2());
		memcpy(&newguys[selected[1]], &bestguys[random_bestguy], sizeof(newguys[selected[1]]));
		newguys[selected[1]].t = guys[selected[1]].t;
		newguys[selected[1]].gen = arr->gen;
		//Now selected[1] is some guy high on the temperature ladder with amazing bestguy-genes and fitness

		expX[1] = exp(-newguys[selected[1]].H / s); //good guy
		expX[0] = exp(-guys[selected[0]].H / s); //bad guy
		PXX = 1 / ((N - 1)*Z)*(expX[0] + expX[1]); //P((xi,xj) | x)
		//Now mate the newguys to create offsprings
		crossoverPoint = tri_inv_randint(0, geneLength) + (int)(genrand_real1()*geneLength*(nGenes - 1));
		//crossoverPoint = (int)(genrand_real1()*geneLength*(nGenes - 1));

		for (i = crossoverPoint; i < genomeLength; i++) {
			newguys[selected[0]].genome[i] = newguys[selected[1]].genome[i]; //Psssst this newguy 1 is the awesome one
			newguys[selected[1]].genome[i] = newguys[selected[0]].genome[i];
		}
		//From now on the newguys are offsprings. The parents are a good guy selected 1 and a bestguy random_bestguy. Adoptive father is guy selected 1
		getParameters4One(newguys[selected[0]], prop);
		getParameters4One(newguys[selected[1]], prop);
		newguys[selected[0]].H = fitnessTest(newguys[selected[0]].parameter, selected[0], arr);
		newguys[selected[1]].H = fitnessTest(newguys[selected[1]].parameter, selected[1], arr);
		expY[0] = exp(-newguys[selected[0]].H / s);
		expY[1] = exp(-newguys[selected[1]].H / s);
		Z = 0;
		for (i = 0; i < N; i++) {
			Z += exp(-newguys[i].H / s);
		}
		PYY = 1 / ((N - 1)*Z)*(expY[0] + expY[1]); //P((xi,xj) | x)
		dHt2 = (newguys[selected[0]].H - guys[selected[0]].H) / guys[selected[0]].t; //good 
		dHt1 = (newguys[selected[1]].H - guys[selected[1]].H) / guys[selected[1]].t; //bad 

		rc = exp(-dHt1 - dHt2)*PXX / PYY;
		//Accept or reject
		if (genrand_real1() < fmin(1, rc)) {
			for (i = 0; i < 2; i++) { //refresh newpos
				j = selected[i];
				memcpy(&guys[j], &newguys[j], sizeof(guys[j]));
			}
		}
		//Make sure newguys are up to speed on newest events
		for (i = 0; i < 2; i++) {
			j = selected[i];
			memcpy(&newguys[j], &guys[j], sizeof(newguys[j]));
		}
	}

}
void crossover_smartCopy(personality guys[], personality newguys[], arrays *arr, properties *prop) {
	//Start roulette selection																							
	int i, j, matings;
	int nMatings = (int)(N / 5);
	int selected[2];
	double s = guys[(int)(0.9*N)].t;
	double rc, PXX, PYY, dHt0, dHt1;
	double expX[2];
	int n_ab[4];
	double Z;
	for (matings = 0; matings < nMatings; matings++) {
		for (i = 0; i < 4; i++) {
			n_ab[i] = 0;
		}
		roulette_select(guys, selected, &Z, s); //Selected 0 will be good, selected 1 random
		expX[0] = exp(-guys[selected[0]].H / s); //good guy
		expX[1] = exp(-guys[selected[1]].H / s); //bad guy
		PXX = 1 / ((N - 1)*Z)*(expX[0] + expX[1]); //P((xi,xj) | x)
		 //Now mate the newguys to create offsprings
		bitselector_smartCopy(guys, newguys, selected, n_ab);
		PYY = 1;
		for (i = 0; i < 4; i++) {
			PYY *= pow(p_matrix[i], n_ab[i]);
		}
		getParameters4One(newguys[selected[0]], prop);
		getParameters4One(newguys[selected[1]], prop);
		newguys[selected[0]].H = fitnessTest(newguys[selected[0]].parameter, selected[0], arr);
		newguys[selected[1]].H = fitnessTest(newguys[selected[1]].parameter, selected[1], arr);
		dHt0 = (newguys[selected[0]].H - guys[selected[0]].H) / guys[selected[0]].t; //good 
		dHt1 = (newguys[selected[1]].H - guys[selected[1]].H) / guys[selected[1]].t; //bad 

		rc = exp(-dHt0 - dHt1)*PXX / PYY;
		//Accept or reject
		if (genrand_real1() < fmin(1, rc)) {
			for (i = 0; i < 2; i++) { //refresh newpos
				j = selected[i];
				memcpy(&guys[j], &newguys[j], sizeof(guys[j]));
			}
		}
		//Make sure newguys are up to speed on newest events
		for (i = 0; i < 2; i++) {
			j = selected[i];
			memcpy(&newguys[j], &guys[j], sizeof(newguys[j]));
		}
	}
}
void exchange(personality guys[], personality newguys[], arrays *arr) {
	int i, j, ex; //ex is the exchange iteration
	double re, dt_inv, dH;
	for (ex = 0; ex < N; ex++) {
		guys[ex].stay_time++;
		i = (int)(N*genrand_real2());
		if (i == 0) { j = 1; }
		else if (i == N - 1) { j = N - 2; }
		else { j = i + 2 * (int)(genrand_real1() + 0.5) - 1; }
		dH = guys[i].H - guys[j].H;
		dt_inv = 1 / guys[i].t - 1 / guys[j].t;
		re = exp(dH*dt_inv);
		if (genrand_real1() < fmin(1, re)) {
			guys[i].stay_time = 0;
			guys[j].stay_time = 0;
			//Copy guy[i] to newguy[j] and vice versa
			memcpy(&newguys[j], &guys[i], sizeof(newguys[i]));
			memcpy(&newguys[i], &guys[j], sizeof(newguys[j]));
			//But keep temperature
			newguys[j].t = guys[j].t;
			newguys[i].t = guys[i].t;
			//Now update the old guys
			memcpy(&guys[i], &newguys[i], sizeof(guys[i]));
			memcpy(&guys[j], &newguys[j], sizeof(guys[j]));
		}
	}
}
void insertguy(personality bestguys[], personality guys[], int from, int to) {
	int i;
	//Push out the worst bestguy and move everybody up until "to"
	for (i = 0; i < to; i++) {
		memcpy(&bestguys[i], &bestguys[i + 1], sizeof(bestguys[i]));
	}
	//Insert the exceptional guy into the illustrious group of best guys
	memcpy(&bestguys[to], &guys[from], sizeof(bestguys[to]));
}
void wakeUpBest(personality bestguys[], personality guys[]) {
	int i;
	int count = 0;
	//int best[N_best];
	int copied_guys[N_best];
	double lowest_H;
	int lowest_i;
	//printf("Hej\n");
	while (count < N_best) {
		lowest_H = 1e10;
		//Find the best guy yet among guys
		for (i = 0; i < N; i++) {
			//Check if i is in skip-list
			if (isvalueinarray(i, copied_guys, N_best) == 1) { continue; }
			if (guys[i].H < lowest_H) {
				lowest_H = guys[i].H;
				lowest_i = i;
				//printf("%d\n", lowest_i);
			}
		}
		//By now we should have a winner, copy him to bestguys and add him to skiplist
		memcpy(&bestguys[N_best - count - 1], &guys[lowest_i], sizeof(bestguys[N_best - count - 1]));
		copied_guys[count] = lowest_i;
		//printf("%d	%d	%f	%f\n", count, lowest_i, lowest_H, bestguys[N_best - count - 1].H);
		count++;

	}
}
void find_elite(personality guys[], personality bestguys[], arrays *arr) {
	int i;
	int count = 0;
	//int best[N_best];
	int tried_guys[N_best];
	double lowest_H;
	int lowest_i;
	int j = 0;
	int success = 1;
	while (success == 1) {
		lowest_H = 1e10;
		//Find the best guy yet among guys
		for (i = 0; i < N; i++) {
			//Check if i is in skip-list
			if (isvalueinarray(i, tried_guys, j) == 1) { continue; }
			if (guys[i].H < lowest_H) {
				lowest_H = guys[i].H;
				lowest_i = i;
			}
		}
		count++;
		//By now we should have a winner, check if he's better than any elite
		//We have already tried to match j guys, so we can start from j?
		for (i = 0; i < N_best; i++) {
			success = 0;
			if (lowest_H == bestguys[N_best - i - 1].H) {
				success = 1;
				tried_guys[j] = lowest_i;
				j++;
				return;//The guy is already here! Do not copy, just check next guy
			}
			if (lowest_H < bestguys[N_best - i - 1].H) {
				insertguy(bestguys, guys, lowest_i, N_best - i - 1);
				success = 1;
				tried_guys[j] = lowest_i;
				j++;
				break;
			}
		}

	}
}
void evolve(personality guys[], personality newguys[], personality bestguys[], arrays *arr, properties *prop) {
	//This function selects either mutation type operators or crossover type. 
	//Then does an exchange operator, and finally finds the best guys in the population
	int i;
	arr->evolution_tic = clock(); //Start timer for this generation
	//Update generation number
	for (i = 0; i < N; i++) {
		guys[i].gen = arr->gen;
		newguys[i].gen = arr->gen;
	}
	//Select mutation or crossover
	if (genrand_real1() < qm) {
		if (genrand_real2() < qma) {
			mutation_elite(guys, newguys, bestguys, arr, prop);
		}
		else {
			mutation(guys, newguys, arr, prop);
		}
	}
	else {
		if (genrand_real1() < qmc) {
			if (genrand_real1() < 0.5) {
				crossover_elite(guys, newguys, bestguys, arr, prop);
			}
			else {
				crossover_smartCopy(guys, newguys, arr, prop);
			}
		}
		else {
			crossover(guys, newguys, arr, prop);
		}
	}
	exchange(guys, newguys, arr);
	find_elite(guys, bestguys, arr);
	arr->evolution_toc = clock(); //End timer for this generation
	arr->evolution_time += (double)(arr->evolution_toc - arr->evolution_tic) / CLOCKS_PER_SEC; //Collect computation time
}

//Create necessary files and folders for data storage
void build_data_files(datafiles *data) {
	//Create folder for data storage
	char foldername[64];
	char subfoldername[64];
	sprintf(foldername, "mkdir data");
	sprintf(subfoldername, "mkdir data\\population");
	system(foldername);			//Make a directory to store data files
	system(subfoldername);			//Make a directory to store data files
	//Create filenames for datafiles

	//char filename_surf[64];
	char filename_dat[64];
	char filename_el[64];

	//sprintf(filename_surf, "data/surface.dat");
	sprintf(filename_dat, "data/parameters.dat");
	sprintf(filename_el, "data/elitedata.dat");

	data->data = fopen(filename_dat, "w+");
	//data->surfdata = fopen(filename_surf, "w+");
	data->elitedata = fopen(filename_el, "w+");
}
void print_to_file(personality guys[N], personality bestguys[N_best], properties *prop, datafiles *data) {
	int i, j;
	for (j = 0; j < nGenes; j++) {
		fprintf(data->data, "Parameter[%d]	", j);
		fprintf(data->elitedata, "Parameter[%d]	", j);
	}
	fprintf(data->data, "Fitness	Temperature	Generation\n");
	fprintf(data->elitedata, "Fitness	Temperature	Generation\n");
	for (i = 0; i < N; i++) {
		for (j = 0; j < nGenes; j++) {
			fprintf(data->data, "%8.5f	", guys[i].parameter[j]);
		}
		fprintf(data->data, "%8.5f	%8.5f	%d\n", guys[i].H, guys[i].t, guys[i].gen);
	}
	for (i = 0; i < N_best; i++) {
		for (j = 0; j < nGenes; j++) {
			fprintf(data->elitedata, "%8.3f	", bestguys[i].parameter[j]);
		}
		fprintf(data->elitedata, "%8.3f	%8.3f	%d\n", bestguys[i].H, bestguys[i].t, bestguys[i].gen);
	}
}
int main(void) {
	int k;
	//Start up some files and folder for saving data
	datafiles data;
	build_data_files(&data);
	//Start upp structures and arrays
	init_genrand(seed);
	personality guys[N]; //Make an array of N guys
	personality newguys[N]; //Make a temporary array of N guinneapigs
	personality bestguys[N_best]; //Make an array of N/10 good performers
	arrays arr;
	arr.loci = (int*)malloc(genomeLength*sizeof(int));
	arr.evolution_time = 0;
	arr.simulation_time = 0;
	arr.store_counter = 0;
	properties prop; 
	setupParameters(&prop);

	//Start algorithm

	wakeUp(guys, &prop);	//Let the guys get some initial random DNA, parameters and a temperature.
	getFit(guys, &arr); 	//Send the guys do some initial balancing to get a lousy fitness score H.
	wakeUpBest(bestguys, guys);
	memcpy(&newguys, &guys, sizeof(newguys)); //let the new guys start off as the old ones
	for (k = 0; k < generations; k++) {
		arr.gen = k+1;
		evolve(guys, newguys, bestguys, &arr, &prop); 	//Send the guys to training camp and let them evolve a bit
		printf("\rRunning... %.1f %%", 100.0*k / (generations - 1));
		//for (i = 0; i < N; i++) {
		//	for (j = 0; j < nGenes; j++) {
		//		fprintf(data.surfdata, "%8.5f	", guys[i].parameter[j]);
		//	}
		//	fprintf(data.surfdata, "%8.5f	%8.5f\n", bestguys[i].H, bestguys[i].t);
		//}
	}
	printf("\n");
	//Print data to files
	print_to_file(guys, bestguys, &prop, &data);
	//Print timing to console
	printf("Simulation time:	%.3f seconds\n", arr.simulation_time);
	printf("Evolution time:		%.3f seconds\n", arr.evolution_time);

	//Free allocated arrays
	free(arr.loci);
	fclose(data.data);
	//fclose(data.surfdata);
	fclose(data.elitedata);
	return 0;
}