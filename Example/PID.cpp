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
#define N 100
#define N_best (int)(N/10.0)
#define geneLength 24
#define nSelect 2
#define nGenes 6
int genomeLength = nGenes*geneLength; //Number of bits for all genes

int generations = 200;
double qm = 0.3; //The mutation probability vs crossover (1-qm for crossover)
double qma = 0.2; //Adaptive mutation rate (qm*(1-qma) for regular mutation)
double qmc = 0.2; //Adaptive crossover rate((1-qm)*(1-qmc) for regular crossover) 
double Tmin = 0.01;
double Tmax = 5;
int seed = 8;

//Rod parameters
double L = 0.5;		// Length of the rod
double l = L / 2; 	//Distance to center of mass
double m = 0.1;		// mass
double g = 9.82; //Gravitational constant
double I = m / 3.0*l*l;	//Inertial moment
double c = g / L;

//Cart parameters
double M = 0.4;		//cart mass
double Fmax = 30; //Around 5 Newtons?
double vmax = 10;
double friction = 0.1;

//Rod simulation parameters (for visualization and balancing, startup etc)
const double pi = 3.141592653589793238462;
double tmax = 3;			//how long to run balancing simulation
double dt = 0.002;			//time step length
int simLength = (int)ceil(tmax / dt);
int fps = 60;				//target framrate for visualization
int plot_freq = (int)(1.0 / dt / fps);//frequency for writing data to file
//double th_start = pi*1.1; 	//starting angle
double th_start;
double dth_start; 		//starting angular velocity
//double dth_start = 0.02; 		//starting angular velocity

double x_start = 0;			//starting x-position
double dx_start = 0;		//starting x-velocity
double Fx_start = 0;
double xmax = 1.0;			//Maximum x-displacement
double thmax = pi / 2;


//Probabilities for smart copy crossover
double p0 = 0.05;
double p1 = 0.15;
double p2 = 0.8;
double p_matrix[] = {p0*p0 + (1-p0)*(1-p0), 2*p0*(1-p0), p1*(1-p2)+p2*(1-p1), p1*p2 + (1-p1)*(1-p2)};



typedef struct{
	int genome[nGenes * geneLength]; //array 0's and 1's that constitute the genome. 3 parameters coded in 32 bits.
	double H;		//Fitness, or energy
	double parameter[nGenes];		//Proportional parameter	//Derivative parameter//Integral parameter
	double t; 		//temperature
	int stay_time;
} personality;

typedef struct{
	//char model;
	//int terms;
	int n_params;
	double param_min[nGenes];
	double param_max[nGenes];
}properties;

typedef struct{
	double *errorTh;
	double *errorX;
	double *speedX;
	double *errordTh;
	//double *errorTop; //Displacement of the rod top from reference value
	double *Fx;
	int *loci;
	bool store_flag;
	int store_counter;
	int gen;
	double simulation_time;
	clock_t simulation_tic;
	clock_t simulation_toc;
	double evolution_time;
	clock_t evolution_tic;
	clock_t evolution_toc;
}arrays;



typedef struct datafiles {
	FILE *data;
	FILE *surfdata;
	FILE *elitedata;
} datafiles;

double integral(double f[], int steps) {
	int i;
	double datIntegral = 0;
	for (i = 0; i < steps - 1; i++) {
		datIntegral += dt*0.5*(f[i] + f[i + 1]);
	}
	//printf("%f\n",datIntegral);
	//getchar();
	return datIntegral;
}
double areaIntegral(double f[], int steps) {
	int i;
	double datArea = 0;
	for (i = 0; i < steps - 1; i++) {
		datArea += dt*0.5*fabs(f[i] + f[i + 1]);
	}
	
	return datArea;
}

double acc_x(double *dx, double *th, double *dth, double ddth, double *Fx){
	//double ddx;
	//ddx = Fx/m + sin(th)*l*dth*dth - ddth*l*cos(th);
	//ddx = (Fx-m*L*ddth*cos(th)+m*L*dth*dth*sin(th))/(m+M);

	return  (*Fx - friction**dx - m*l*(ddth*cos(*th) - *dth**dth*sin(*th))) / (M + m);
	//return ddx; //return double derivative of x
}

double acc(double *th, double *dth, double *ddx){
	//double ddth;
	//ddth = -(ddx*l*cos(th) - sin(th)*(dx*dth - dth*dth*dx*l - g*l)) / (l*l + I / m);
	//ddth = -((g+l*cos(th)*dth*dth)*l*sin(th)+Fx/m*sin(th))/(l*l*sin(th)*sin(th)-I/m);
	return (-g*sin(*th) - *ddx*cos(*th)) / (l + I / m / l);
	//return ddth; //return the double derivative of theta
}

void rungekutta(double *dt, double *th, double *dth, double *dx, double *ddx, double *Fx){
	double a1, a2, a3, a4, b1, b2, b3, b4;
	double th_b1, dth_a1, th_b2, dth_a2, th_b3, dth_a3;
	//Runge-Kutta
	a1 = acc(th, dth, ddx)**dt;
	b1 = *dth**dt;
	th_b1 = *th + b1 / 2;
	dth_a1 = *dth + a1 / 2;

	a2 = acc(&th_b1, &dth_a1, ddx)**dt;
	b2 = (dth_a1)**dt;

	th_b2 = *th + b2 / 2;
	dth_a2 = *dth + a2 / 2;

	a3 = acc(&th_b2, &dth_a2, ddx)**dt;
	b3 = (dth_a2)**dt;

	th_b3 = *th + b3;
	dth_a3 = *dth + a3;

	a4 = acc(&th_b3, &dth_a3, ddx)**dt;
	b4 = (dth_a3)**dt;

	*dth += (1.0 / 6)*(a1 + 2 * (a2 + a3) + a4);
	*th += (1.0 / 6)*(b1 + 2 * (b2 + b3) + b4);
}

void Force(double *F, double *dx, double *dth, double *e, double *ex, double *e_int, double *ex_int, double parameter[]){
	//e and ex are the error of th and x, or deviation from the reference value of pi rad
	//-dth The derivative of the error of th
	//-dx The derivative of error of x

	//double FX = parameter[3]*ex[step]; //+ parameter[4]*dex + parameter[5]*integral(ex,step);
	//*F = parameter[3] * (-*dx) + parameter[0] * *e + parameter[1] * (-*dth) + parameter[2] * *e_int;
	*F = parameter[0] * *e + parameter[1] * (-*dth) + parameter[2] * *e_int; //

	//Activate the x-control only if angular velocity is low (i.e. if system is calm enough)
	//Activate it only if the rod is pointing towards origo
	if (fabs(-*dth) < 2*pi/60 && fabs(*e) < 20/ 360 * 2 * pi){
	//if (abs(*e) < 60 / 360 * 2 * pi){

		if (*e <= 0 && *ex < 0){
			*F += +1 * parameter[3] * (*ex) + 1 * parameter[4] * (-*dx) + 1 * parameter[5] * (*ex_int);
		}
		else if (*e  >= 0 && *ex > 0){
			*F += +1 * parameter[3] *  (*ex) + 1 * parameter[4] * (-*dx) + 1 * parameter[5] * (*ex_int);
		}
	}
	//*F += +1 * parameter[3] * *ex + 1 * parameter[4] * (-*dx) + 1 * parameter[5] * (*ex_int);

	if (*F > Fmax){
		*F = Fmax;
	}
	else if (*F < -Fmax){
		*F = -Fmax;
	}
	//printf("%f\n",*F);
}

double fitnessTest(double parameter[], int guynumber, arrays *arr){
	arr->simulation_tic = clock();
	double th, dth, errorTh, errorX;// , errorTop;
	double x, dx, ddx, ddx_dt;
	double t = 0;
	double Fx, H;
	double errorTh_integral, errorX_integral;// , errorTop_integral;
	double error_line_integral = 0;
	double errorTh_line_integral = 0;
	double equilibrium_time[50];
	equilibrium_time[0] = 0;
	int eq = 0;
	int i;
	int survival_steps = 0;
	th_start = 1.1*pi;// pi + bias + random
	dth_start = -0.1;
	th = th_start;
	x = x_start;
	dx = dx_start;
	dth = dth_start;
	ddx = 0;
	Fx = Fx_start;
	//Reference values
	errorTh = pi - th;
	errorX = 0 - x;
	//errorTop = pow(pow(-errorX + L*sin(errorTh), 2) + pow(L*(1 - cos(errorTh)), 2), 0.5);
	//Their integrals, if used in force calculation
	errorTh_integral = 0;
	errorX_integral = 0;
	//errorTop_integral = 0;
	//Otherwise just store in arrays
	arr->errorTh[0] = errorTh;
	arr->errordTh[0] = dth;
	arr->errorX[0] = errorX;
	arr->speedX[0] = dx;
	//arr->errorTop[0] = pow(pow(-errorX + L*sin(errorTh), 2) + pow(L*(1-cos(errorTh)),2),0.5);
	arr->Fx[0] = 0;

	//Store data if last generation
	//Try only scoring to reference value errorTop, not using it as feedback
	FILE *elitefile = NULL;
	if (arr->store_flag == true){
		char filename[100];
		sprintf(filename, "data/population/data%d.dat", arr->store_counter);
		elitefile = fopen(filename, "w+");
		fprintf(elitefile, "%8.5f	%8.5f	%8.5f	%8.5f	%8.5f	%8.5f	%8.5f	%d	%8.5f	%d\n", t, x, dx, th, dth, Fx, L, fps, arr->errorTh[0], arr->gen);
	}

	for (i = 1; i < simLength; i++){
		//Insert current values to arrays
		errorTh = pi - th;
		errorX = 0 - x;
		arr->errorTh[i] = errorTh;
		arr->errorX[i] = errorX;
		errorTh_integral += dt*0.5*(errorTh + arr->errorTh[i - 1]);
		errorX_integral += dt*0.5*(errorX + arr->errorX[i - 1]);
		Force(&Fx, &dx, &dth, &errorTh, &errorX, &errorTh_integral, &errorX_integral, parameter);
		arr->Fx[i] = Fx;
		//Perform the iteration
		t += dt;
		ddx = acc_x(&dx, &th, &dth, acc(&th, &dth, &ddx), &Fx);
		ddx = fmin(vmax, ddx); //Enforce maximum cart velocity
		ddx = fmax(-vmax, ddx);
		rungekutta(&dt, &th, &dth, &dx, &ddx, &Fx); //modifies th and dth
		ddx_dt = acc_x(&dx, &th, &dth, acc(&th, &dth, &ddx), &Fx);
		x += dx*dt + 0.5*ddx*dt*dt;
		dx += 0.5*(ddx_dt + ddx)*dt;
		arr->speedX[i] = dx;
		arr->errordTh[i] = dth;
		error_line_integral += dt*0.5*(pow(1 + dx*dx, 0.5) + pow(1 + arr->speedX[i - 1] * arr->speedX[i - 1], 0.5));
		errorTh_line_integral += dt*0.5*(pow(1 + dth*dth, 0.5) + pow(1 + arr->errordTh[i - 1] * arr->errordTh[i - 1], 0.5));
		//printf("%f %f %f %f %f %f %f %f\n", t, dx, ddx, ddx_dt, x, errorTh, errorX, Fx);
		if (isnan(errorX)){
			printf("Error! errorX: %f\n", errorX);
			//exit(1);
		}
		if (arr->store_flag == true){
			if (i % plot_freq == 0){
				fprintf(elitefile, "%8.5f	%8.5f	%8.5f	%8.5f	%8.5f	%8.5f	%8.5f	%d	%8.5f	%d\n", t, x, dx, th, dth, Fx, L, fps, arr->errorTh[i], arr->gen);
			}
		}
		if (get_sign(arr->errorTh[i - 1] * arr->errorTh[i]) < 0){
			eq++;
			equilibrium_time[eq] = t  - equilibrium_time[eq-1];
			//printf("Eq = %f\n", equilibrium_time[eq]);
		}

		//Fill the rest of the vector to shame bad performers
		//if (abs(errorTh) > thmax || abs(errorX) > xmax){
		//	for (j = survival_steps; j < simLength; j++){
		//		arr->errorTh[j] = errorTh / fmin(1, t);// errorTh;
		//		arr->errorX[j] = errorX / fmin(1, t);// errorX;
		//		arr->speedX[j] = dx/fmin(1,t);//
		//		//arr->errorTop[j] = pow(pow(-errorX + L*sin(errorTh), 2) + pow(L*(1 - cos(errorTh)), 2), 0.5);
		//		arr->Fx[j] = Fx;
		//	}
		//	break;

		//}
		survival_steps++;
	}
	if (arr->store_flag == true){
		fclose(elitefile);
	}
	//traveled distance! minimize it too!
	//double H1 = log((areaIntegral(arr->errorTh, i - 1))) / t;
	//double H2 = log((areaIntegral(arr->errorX, i - 1))) / t;
	//double H3 = log(areaIntegral(arr->speedX, i - 1)) / t;
	double H1 = areaIntegral(arr->errorTh, simLength);
	double H2 = areaIntegral(arr->errorX, simLength);
	//double H3 = fabs(integral(arr->speedX, simLength));
	//double H4 = areaIntegral(arr->Fx, simLength);
	//double H5 = fabs(dx);
	//double H6 = fabs(errorX);
	//double H7 = error_line_integral / t;
	double H8 = errorTh_line_integral / t;
	double H9;
	if (eq >= 1){
		H9 = mean(equilibrium_time, eq);
		//H9 = equilibrium_time[1];
	}
	else{ H9 = 5; }
	 
	//double H7 = areaIntegral(arr->errorTop, survival_steps);
	//double H8 = arr->errorTop[survival_steps];

	H = H1 +H2 + H8 + H9;// +H4*H4 + exp(1.0);
	if (isnan(H)){
		printf("Error H:	%.4f	%.4f	%.4f	%.4f	%d\n", H1, H2, H8, H9, survival_steps);
		exit(1);
	}



	H = -1 / log(H + 0.5);// +pow(H, 0.75);
	arr->evolution_tic = arr->simulation_tic;
	arr->simulation_toc = clock();
	arr->simulation_time += (double)(arr->simulation_toc - arr->simulation_tic) / CLOCKS_PER_SEC;
	return H;
}

void getFit(personality guys[], arrays *arr){
	int i;
	for (i = 0; i < N; i++){
		guys[i].H = fitnessTest(guys[i].parameter, i, arr);
		//printf("hej %.3f\n",guys[i].H);
	}

}
double bin2dec(int bin[], int binSize) {
	int i;
	//printf("%d\n",binSize);
	double dec = 0;
	for (i = 0; i < binSize; i++){
		//printf("%d %d\n",i,bin[i]);
		if (bin[i] == 1){
			dec += pow(2.0, binSize - i - 1);
		}
	}
	return dec;
}
double map2param(double param, int i, properties *prop){
	double result;
	result = param / (pow(2.0, geneLength) - 1);
	result *= (prop->param_max[i] - prop->param_min[i]);
	result += prop->param_min[i];
	return result;
}
void getParameters4One(personality &guy, properties *prop){
	int i, j;
	int gene[geneLength];
	//printf("%d\n",i);
	for (i = 0; i < nGenes; i++){
		for (j = 0; j < geneLength; j++){
			gene[j] = guy.genome[i*geneLength + j];
		}
		guy.parameter[i] = bin2dec(gene, sizeof(gene) / sizeof(gene[0]));
		guy.parameter[i] = map2param(guy.parameter[i], i, prop);
	}
}
void getParameters4All(personality guys[], properties *prop){
	int i, j, k;
	int gene[geneLength];
	for (i = 0; i < N; i++){
		for (j = 0; j < nGenes; j++){
			for (k = 0; k < geneLength; k++){
				gene[k] = guys[i].genome[j*geneLength + k];
			}
			guys[i].parameter[j] = bin2dec(gene, sizeof(gene) / sizeof(gene[0]));
			guys[i].parameter[j] = map2param(guys[i].parameter[j], j, prop);
		}
	}
}

void wakeUp(personality guys[], properties *prop){
	int i, j;
	double T[N];
	logspace(T, Tmax, Tmin, N);
	//Initialize some temperature ladder, here uniform.
	//for (i = 0; i < N; i++){
	//T[N - i - 1] = Tmin + (Tmax - Tmin)*i / (N - 1);
	//T[i] = Tmin + (Tmax-Tmin)*i/(N-1);
	//}
	//Initialize the genome sequence for all guys randomly
	for (i = 0; i < N; i++){
		guys[i].t = T[i]; 		//Assign temperatures produced above
		guys[i].stay_time = 0;
		for (j = 0; j < genomeLength; j++){
			guys[i].genome[j] = (int)nearbyint(genrand_real1());	//Distribute the genetic parameters
		}
	}
	getParameters4All(guys, prop);
}

double boltzmann(personality guys[]){
	int i;
	double Z = 0;
	double p[N]; //Boltzman selection probabilities for an individual
	double f;	//Boltzman distribution for the population
	for (i = 0; i < N; i++){
		Z += exp(-guys[i].H / guys[i].t);
	}
	for (i = 0; i < N; i++){
		p[i] = exp(-guys[i].H / guys[i].t) / Z;
	}
	f = p[0];
	for (i = 1; i < N; i++){
		f *= p[i];
	}
	return f;
}

void mutation(personality guys[], personality newguys[], arrays *arr, properties *prop){
	int i, j;
	int mutantGenes;	//points to mutate
	int mutant;			//which guy to mutate
	for (i = 0; i < N; i++){
		mutantGenes = (int)(genomeLength*genrand_real2());
		mutant = (int)(N*genrand_real2());
		rndChoice(arr->loci, mutantGenes, genomeLength);	//Choose which locus to mutate
		for (j = 0; j < mutantGenes; j++){
			//newguys[mutant].genome[arr->loci[j]] = (int)nearbyint(genrand_real1()); //reverse the gene(s)
			newguys[mutant].genome[arr->loci[j]] = 1 - newguys[mutant].genome[arr->loci[j]]; //reverse the gene(s)

		}
		getParameters4One(newguys[mutant], prop);
		newguys[mutant].H = fitnessTest(newguys[mutant].parameter, mutant, arr);
		//Perform metropolis
		double dH = newguys[mutant].H - guys[mutant].H;		//used to decide if we accept the new guy or not.
		j = mutant;
		if (dH < 0){
			memcpy(&guys[j], &newguys[j], sizeof(guys[j]));
		}
		else if (exp(-dH / newguys[mutant].t) > genrand_real1()){
			memcpy(&guys[j], &newguys[j], sizeof(guys[j]));
		}
		//Inform newguys of latest events... i.e sync them for the next round
		memcpy(&newguys[j], &guys[j], sizeof(newguys[j]));
	}
}

void elite_mutation(personality guys[], personality newguys[], personality bestguys[], arrays *arr, properties *prop){
	int i, j;
	int mutantGenes;	//Number of points to mutate
	int mutant;			//which guy to mutate
	int elite_mutant;	//Which elite guy to receive wisdom from
	for (i = 0; i < N; i++){
		mutantGenes = 0;
		while (mutantGenes == 0){
			mutantGenes = (int)(genomeLength*genrand_real2());
		}
		mutant = (int)(N*genrand_real2());
		elite_mutant = (int)(N_best*genrand_real2());
		//Fill loci with mutantGenes genome points to be mutated, from a triangular distribution for each parameter.
		int count = 0;
		int test_int;
		while (count < mutantGenes){
			test_int = tri_randint(0, geneLength) + (int)(genrand_real1()*geneLength*(nGenes - 1));
			arr->loci[count] = test_int;
			count++;
			
		}


		//Copy an elite guy to a new guy to be a guinnea pig
		//memcpy(&newguys[mutant], &bestguys[elite_mutant], sizeof(newguys[mutant]));
		memcpy(&newguys[mutant].genome, &bestguys[elite_mutant].genome, sizeof(newguys[mutant].genome));

		for (j = 0; j < mutantGenes; j++){
			newguys[mutant].genome[arr->loci[j]] = 1 - newguys[mutant].genome[arr->loci[j]]; //reverse the gene(s)
		}
		getParameters4One(newguys[mutant], prop);
		newguys[mutant].H = fitnessTest(newguys[mutant].parameter, mutant, arr);
		//Perform metropolis
		double dH = newguys[mutant].H - guys[mutant].H;		//used to decide if we accept the new guy or not.
		j = mutant; //Shorter notation...
		if (dH < 0){
			memcpy(&guys[j], &newguys[j], sizeof(guys[j]));
			guys[j].stay_time = 0;
		}
		else if (exp(-dH / newguys[j].t) > genrand_real1()){
			memcpy(&guys[j], &newguys[j], sizeof(guys[j]));
			guys[j].stay_time = 0;

		}
		//Inform newguys of latest events... i.e sync them for the next round
		memcpy(&newguys[j], &guys[j], sizeof(newguys[j]));
	}
}

void roulette_select(personality guys[], int selected[], double *Z, double s){
	//Selected 0 is the guy with good fitness (lowest), selected 1 is random

	int i;
	//double total_H = 0;
	double roulette[N];
	double lucky_number;
	//Make a roulette wheel
	*Z = 0;
	for (i = 0; i < N; i++){
		*Z += exp(-guys[i].H / s);
		roulette[i] = *Z; //Low fitness gives large area on roulette
	}
	lucky_number = *Z*genrand_real1();
	for (i = 0; i < N; i++){
		if (lucky_number < roulette[i]){
			selected[0] = i;
			break;
		}
	}
	selected[1] = selected[0];
	while (selected[1] == selected[0]){
		selected[1] = (int)((N)*genrand_real2());
	}
}

void elite_roulette_select(personality guys[], personality bestguys[], int selected[], double *Z, double s){
	//Selected 0 is some random bestguy. Selected 1 is a regular guy selected with boltzmann weight
	int i;
	//double total_H = 0;
	double roulette[N];
	double lucky_number;
	//Make a roulette wheel
	*Z = 0;
	for (i = 0; i < N; i++){
		//if (exp(-guys[i].H / guys[i].t) > 1e6){
		//	printf("%f %d %f %f\n", *Z, i, guys[i].H,exp(-guys[i].H / guys[i].t));
		//}
		*Z += exp(-guys[i].H / s);
		roulette[i] = *Z; //Low fitness gives large area on roulette
	}
	lucky_number = *Z*genrand_real1();
	for (i = 0; i < N; i++){
		if (lucky_number < roulette[i]){
			selected[1] = i;
			break;
		}
	}
	selected[0] = selected[1];
	while (selected[0] == selected[0]){
		selected[1] = (int)((N_best)*genrand_real2());
	}
}

void elite_crossover(personality guys[], personality newguys[], personality bestguys[], int gen, arrays *arr, properties *prop){
	//Start roulette selection
	//Should the fitnesses be ordeed as H(xi) >=  H(xk)?
	int i, j, matings;
	//int nSelect = 2;
	int nMatings = (int)(N_best / 5);
	int selected[nSelect];
	int crossoverPoint;
	double rc, PXX, PYY, dHt1, dHt2;
	double expX[nSelect];
	double expY[nSelect];
	double Z;
	double s = guys[(int)(0.97*N)].t;
	//int random_gene;
	int random_bestguy;
	for (matings = 0; matings < nMatings; matings++){
		roulette_select(guys, selected, &Z, s);
		//Selected 0 is a boltzmann "just good" guy. Selected 1 is a random any guy.
		//Let a newguy selected 1 impersonate a random bestguy but keep the temperature. Then newguy gets superpowers
		random_bestguy = (int)(N_best*genrand_real2());
		memcpy(&newguys[selected[1]], &bestguys[random_bestguy], sizeof(newguys[selected[1]]));
		newguys[selected[1]].t = guys[selected[1]].t;
		//Now selected[1] is some guy high on the temperature ladder with amazing bestguy-genes and fitness

		expX[1] = exp(-newguys[selected[1]].H / s); //good guy
		expX[0] = exp(-guys[selected[0]].H / s); //bad guy
		PXX = 1 / ((N - 1)*Z)*(expX[0] + expX[1]); //P((xi,xj) | x)
		//Now mate the newguys to create offsprings
		//crossoverPoint = (int)(genrand_real1()*genomeLength);
		//random_gene = geneLength * (int)(nGenes*genrand_real2()); //Chroose crossover point with triangular distribution.
		crossoverPoint = tri_inv_randint(0, geneLength) + (int)(genrand_real1()*geneLength*(nGenes - 1));
		//crossoverPoint = (int)(genrand_real1()*geneLength*(nGenes - 1));

		for (i = crossoverPoint; i < genomeLength; i++){
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
		for (i = 0; i < N; i++){
			Z += exp(-newguys[i].H / s);
		}
		PYY = 1 / ((N - 1)*Z)*(expY[0] + expY[1]); //P((xi,xj) | x)
		dHt2 = (newguys[selected[0]].H - guys[selected[0]].H) / guys[selected[0]].t; //good 
		dHt1 = (newguys[selected[1]].H - guys[selected[1]].H) / guys[selected[1]].t; //bad 

		//if (dHt1 + dHt2 < 0){ continue; }
		rc = exp(-dHt1 - dHt2)*PXX / PYY;
		//Accept or reject
		j = 0;
		if (genrand_real1() < fmin(1, rc)){
			for (i = 0; i < nSelect; i++){ //refresh newpos
				j = selected[i];
				memcpy(&guys[j], &newguys[j], sizeof(guys[j]));
			}
		}
		//Make sure newguys are up to speed on newest events
		for (i = 0; i < nSelect; i++){
			j = selected[i];
			memcpy(&newguys[j], &guys[j], sizeof(newguys[j]));
		}
	}

}

void crossover(personality guys[], personality newguys[], int gen, arrays *arr, properties *prop){
	//Start roulette selection
	//Should the fitnesses be ordered as H(xi) >=  H(xj)? Typ...
	int i, j, matings;
	//int nSelect = 2;
	int nMatings = (int)(N / 5);
	int selected[nSelect];
	int crossoverPoint;
	double s = guys[(int)(0.9*N)].t;
	//printf("s=%f\n",s);
	double rc, PXX, PYY, dHt0, dHt1;
	double expX[nSelect];
	double expY[nSelect];
	double Z;
	//double maxexp = 50; //protect from overflow
	//double minexp = -50; //protect from overflow
	for (matings = 0; matings < nMatings; matings++){
		roulette_select(guys, selected, &Z, s); //Selected 0 will be good, selected 1 random
		//expX[0] = exp(-guys[selected[0]].H / guys[selected[0]].t); //good guy
		//expX[1] = exp(-guys[selected[1]].H / guys[selected[1]].t); //bad guy
		expX[0] = exp(-guys[selected[0]].H / s); //good guy
		expX[1] = exp(-guys[selected[1]].H / s); //bad guy
		PXX = 1 / ((N - 1)*Z)*(expX[0] + expX[1]); //P((xi,xj) | x)
		//Now mate the newguys to create offsprings
		crossoverPoint = (int)(genrand_real1()*genomeLength);
		for (i = crossoverPoint; i < genomeLength; i++){
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
		for (i = 0; i < N; i++){
			Z += exp(-newguys[i].H / s);
		}
		PYY = 1 / ((N - 1)*Z)*(expY[0] + expY[1]); //P((xi,xj) | x)
		dHt0 = (newguys[selected[0]].H - guys[selected[0]].H) / guys[selected[0]].t; //good 
		dHt1 = (newguys[selected[1]].H - guys[selected[1]].H) / guys[selected[1]].t; //bad 

		//if (dHt1 + dHt2 < 0){ continue; }
		rc = exp(-dHt0 - dHt1)*PXX / PYY;
		//Accept or reject
		//j = 0;
		if (genrand_real1() < fmin(1, rc)){
			for (i = 0; i < nSelect; i++){ //refresh newpos
				j = selected[i];
				memcpy(&guys[j], &newguys[j], sizeof(guys[j]));
			}
		}
		//Make sure newguys are up to speed on newest events
		for (i = 0; i < nSelect; i++){
			j = selected[i];
			memcpy(&newguys[j], &guys[j], sizeof(newguys[j]));
		}
	}
}


void adaptive_bitselector(personality guys[], personality newguys[], int selected[], int n_ab[] ){
	int i, j;
	for (i = 0; i < nGenes*geneLength; i++){
		//If parents loci are the same
		if (guys[selected[0]].genome[i] == guys[selected[1]].genome[i]){
			//Copy the values
			newguys[selected[0]].genome[i] = guys[selected[0]].genome[i];
			newguys[selected[1]].genome[i] = guys[selected[1]].genome[i];
			//Independently reverse with probability p0
			for (j = 0; j < nSelect; j++){
				if (genrand_real1() < p0){
					newguys[selected[j]].genome[i] = 1 - newguys[selected[j]].genome[i];
				}
			}
			if (newguys[selected[0]].genome[i] == newguys[selected[1]].genome[i]){
				n_ab[0]++; //Kept the same}
			}
			else{
				n_ab[1]++;//Changed it
			}
		}
		//If parents loci are different
		else {
			//Copy the values
			newguys[selected[0]].genome[i] = guys[selected[0]].genome[i];
			newguys[selected[1]].genome[i] = guys[selected[1]].genome[i];
			//Reverse with probabilities p1 and p2 respectively
			if (genrand_real1() < p1){
				newguys[selected[0]].genome[i] = 1 - newguys[selected[0]].genome[i];
			}
			if (genrand_real1() < p2){
				newguys[selected[1]].genome[i] = 1 - newguys[selected[1]].genome[i];
			}
			if (newguys[selected[0]].genome[i] == newguys[selected[1]].genome[i]){
				n_ab[2]++; //Kept the same}
			}
			else{
				n_ab[3]++;//Changed it
			}
		}
	}
}
void adaptive_crossover(personality guys[], personality newguys[], int gen, arrays *arr, properties *prop) { //Not working yet
	//Start roulette selection
	//Should the fitnesses be ordered as H(xi) >=  H(xj)? Typ...
	int i, j, matings;
	//int nSelect = 2;
	int nMatings = (int)(N / 5);
	int selected[nSelect];
	//int crossoverPoint;
	double s = guys[(int)(0.9*N)].t;
	//printf("s=%f\n",s);
	double rc, PXX, PYY, dHt0, dHt1;
	double expX[nSelect];
//	double expY[nSelect];
	int n_ab[4];
	double Z;

	//double maxexp = 50; //protect from overflow
	//double minexp = -50; //protect from overflow
	for (matings = 0; matings < nMatings; matings++){
		for (i = 0; i < 4; i++){
			n_ab[i] = 0;
		}
		roulette_select(guys, selected, &Z, s); //Selected 0 will be good, selected 1 random
		//expX[0] = exp(-guys[selected[0]].H / guys[selected[0]].t); //good guy
		//expX[1] = exp(-guys[selected[1]].H / guys[selected[1]].t); //bad guy
		expX[0] = exp(-guys[selected[0]].H / s); //good guy
		expX[1] = exp(-guys[selected[1]].H / s); //bad guy
		PXX = 1 / ((N - 1)*Z)*(expX[0] + expX[1]); //P((xi,xj) | x)
		//Now mate the newguys to create offsprings
		//crossoverPoint = (int)(genrand_real1()*genomeLength);
		//for (i = crossoverPoint; i < genomeLength; i++){
		//	newguys[selected[0]].genome[i] = guys[selected[1]].genome[i];
		//	newguys[selected[1]].genome[i] = guys[selected[0]].genome[i];
		//}
		adaptive_bitselector(guys, newguys, selected, n_ab);
		PYY = 1;
		for (i = 0; i < 4; i++){
			PYY *= pow(p_matrix[i], n_ab[i]);
		}

		getParameters4One(newguys[selected[0]], prop);
		getParameters4One(newguys[selected[1]], prop);
		newguys[selected[0]].H = fitnessTest(newguys[selected[0]].parameter, selected[0], arr);
		newguys[selected[1]].H = fitnessTest(newguys[selected[1]].parameter, selected[1], arr);
		//expY[0] = exp(-newguys[selected[0]].H / s);
		//expY[1] = exp(-newguys[selected[1]].H / s);
		//Z = 0;
		//for (i = 0; i < N; i++){
		//	Z += exp(-newguys[i].H / s);
		//}
		//PYY = 1 / ((N - 1)*Z)*(expY[0] + expY[1]); //P((xi,xj) | x)
		dHt0 = (newguys[selected[0]].H - guys[selected[0]].H) / guys[selected[0]].t; //good 
		dHt1 = (newguys[selected[1]].H - guys[selected[1]].H) / guys[selected[1]].t; //bad 

		//if (dHt1 + dHt2 < 0){ continue; }
		rc = exp(-dHt0 - dHt1)*PXX / PYY;
		//Accept or reject
		//j = 0;
		if (genrand_real1() < fmin(1, rc)){
			for (i = 0; i < nSelect; i++){ //refresh newpos
				j = selected[i];
				memcpy(&guys[j], &newguys[j], sizeof(guys[j]));
			}
		}
		//Make sure newguys are up to speed on newest events
		for (i = 0; i < nSelect; i++){
			j = selected[i];
			memcpy(&newguys[j], &guys[j], sizeof(newguys[j]));
		}
	}
}

void exchange(personality guys[], personality newguys[], int gen){
	int i, j, ex; //ex is the exchange iteration
	double re, dt_inv, dH;
	for (ex = 0; ex < N; ex++){
		guys[ex].stay_time++;
		i = (int)(N*genrand_real2());
		if (i == 0){ j = 1; }
		else if (i == N - 1){ j = N - 2; }
		else { j = i + 2 * (int)(genrand_real1() + 0.5) - 1; }
		dH = guys[i].H - guys[j].H;
		dt_inv = 1 / guys[i].t - 1 / guys[j].t;
		re = exp(dH*dt_inv);
		if (genrand_real1() < fmin(1, re)){
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

/* bool isvalueinarray(int val, int *arr, int size){
	int i;
	for (i = 0; i < size; i++) {
		if (arr[i] == val)
			return true;
	}
	return false;
} */
void insertguy(personality bestguys[], personality guys[], int from, int to){
	int i;
	//Push out the worst bestguy and move everybody up until "to"
	for (i = 0; i < to; i++){
		memcpy(&bestguys[i], &bestguys[i + 1], sizeof(bestguys[i]));
	}
	//Insert the exceptional guy into the illustrious group of best guys
	memcpy(&bestguys[to], &guys[from], sizeof(bestguys[to]));
}
void wakeUpBest(personality bestguys[], personality guys[]){
	int i;
	int count = 0;
	//int best[N_best];
	int copied_guys[N_best];
	double lowest_H;
	int lowest_i;
	//printf("Hej\n");
	while (count < N_best){
		lowest_H = 1e10;
		//Find the best guy yet among guys
		for (i = 0; i < N; i++){
			//printf("%f\n",guys[i].H);
			//Check if i is in skip-list
			if (isvalueinarray(i, copied_guys, N_best)){ continue; }
			if (guys[i].H < lowest_H){
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


void find_elite(personality guys[], personality bestguys[]){
	int i;
	int count = 0;
	//int best[N_best];
	int tried_guys[N_best];
	double lowest_H;
	int lowest_i;
	int j = 0;
	bool success = true;
	while (success == true){
		lowest_H = 1e10;
		//Find the best guy yet among guys
		for (i = 0; i < N; i++){
			//Check if i is in skip-list
			if (isvalueinarray(i, tried_guys, j)){ continue; }
			if (guys[i].H < lowest_H){
				lowest_H = guys[i].H;
				lowest_i = i;
			}
		}
		count++;
		//By now we should have a winner, check if he's better than any elite
		//We have already tried to match j guys, so we can start from j?
		for (i = 0; i < N_best; i++){
			success = false;
			if (lowest_H == bestguys[N_best - i - 1].H){
				success = true;
				tried_guys[j] = lowest_i;
				j++;
				return;//The guy is already here! Do not copy, just check next guy
			}
			if (lowest_H < bestguys[N_best - i - 1].H){
				insertguy(bestguys, guys, lowest_i, N_best - i - 1);
				success = true;
				tried_guys[j] = lowest_i;
				j++;
				//copied_guys[count] = lowest_i;
				break;
			}
		}

	}
}

void evolve(personality guys[], personality newguys[], personality bestguys[], int gen, arrays *arr, properties *prop){

	arr->evolution_tic = clock();
	//Select mutation or crossover
	if (genrand_real1() < qm){
		if (genrand_real2() < qma){
			elite_mutation(guys, newguys, bestguys, arr, prop);
		}
		else{
			mutation(guys, newguys, arr, prop);
		}
	}
	else{
		if (genrand_real1() < qmc){
			if (genrand_real1() < 0.5){
				elite_crossover(guys, newguys, bestguys, gen, arr, prop);
			}
			else{
				adaptive_crossover(guys, newguys, gen, arr, prop);
			}
		}
		else{
			crossover(guys, newguys, gen, arr, prop);
		}
	}
	exchange(guys, newguys, gen);

	find_elite(guys, bestguys);
	//find_elite();
	//arr->evolution_toc;
	arr->evolution_toc = clock();
	arr->evolution_time += (double)(arr->evolution_toc - arr->evolution_tic) / CLOCKS_PER_SEC;
}

//Create necessary files and folders for data storage
void build_data_files(datafiles *data){
	//Create folder for data storage
	char foldername[64];
	char subfoldername[64];
	sprintf(foldername, "mkdir data");
	sprintf(subfoldername, "mkdir data\\population");
	system(foldername);			//Make a directory to store data files
	system(subfoldername);			//Make a directory to store data files
	//Create filenames for datafiles

	char filename_surf[64];
	char filename_dat[64];
	char filename_el[64];

	sprintf(filename_surf, "data/surface.dat");
	sprintf(filename_dat, "data/parameters.dat");
	sprintf(filename_el, "data/elitedata.dat");

	data->data = fopen(filename_dat, "w+");
	data->surfdata = fopen(filename_surf, "w+");
	data->elitedata = fopen(filename_el, "w+");
}

int main(void){
	int i, j, k, who;
	//Start up some files and folder for saving data
	datafiles data;
	build_data_files(&data);
	//Start upp structures and arrays
	init_genrand(seed);
	personality guys[N]; //Make an array of N guys
	personality newguys[N]; //Make a temporary array of N guinneapigs
	personality bestguys[N_best]; //Make an array of N/10 good performers

	arrays arr;
	arr.errorTh = (double*)malloc(simLength * sizeof(double));
	arr.errordTh = (double*)malloc(simLength * sizeof(double));
	arr.errorX = (double*)malloc(simLength * sizeof(double));
	arr.speedX = (double*)malloc(simLength * sizeof(double));

	arr.Fx = (double*)malloc(simLength * sizeof(double));
	arr.loci = (int*)malloc(genomeLength*sizeof(int));
	arr.store_flag = false;
	arr.evolution_time = 0;
	arr.simulation_time = 0;
	arr.store_counter = 0;
	properties prop; //
	prop.param_min[0] = 0;
	prop.param_min[1] = 0;
	prop.param_min[2] = 0;
	prop.param_min[3] = 0;
	prop.param_min[4] = 0;
	prop.param_min[5] = 0;
	prop.param_max[0] = 200;
	prop.param_max[1] = 200;
	prop.param_max[2] = 200;
	prop.param_max[3] = 200;
	prop.param_max[4] = 200;
	prop.param_max[5] = 200;
	prop.n_params = nGenes;

	//Start algorithm

	wakeUp(guys, &prop);	//Let the guys get some initial random DNA, parameters and a temperature.
	getFit(guys, &arr); 	//Send the guys do some initial balancing to get a lousy fitness score H.
	wakeUpBest(bestguys, guys);
	memcpy(&newguys, &guys, sizeof(newguys)); //let the new guys start off as the old ones
	for (k = 0; k < generations; k++){
		evolve(guys, newguys, bestguys, k, &arr, &prop); 	//Send the guys to training camp and let them evolve a bit
		//memcpy(&newguys, &guys, sizeof(newguys)); //let the new guys start off as the old ones
		printf("\rRunning... %.1f %%", 100.0*k / (generations - 1));
		for (i = 0; i < N; i++){
			fprintf(data.surfdata, "%.5f	%.5f	%.5f	%.5f	%.5f\n", guys[i].parameter[0], guys[i].parameter[1], guys[i].parameter[2], guys[i].parameter[3], guys[i].H);
		}
		if ( mod(k,200) == 0 || k == 1 || k == 10 || k == 20|| k == generations - 1){
			arr.store_flag = true;
			if (k < 5) {
				arr.gen = k + 1;
				who = (int)(0.3*(N - 1));
				guys[who].H = fitnessTest(guys[who].parameter, i, &arr);
			}
			else {
				arr.gen = k;
				who = (int)(0.3*(N_best - 1));
				if (k == generations - 1) { who = N_best - 1; arr.gen = k + 1; }
				bestguys[i].H = fitnessTest(bestguys[who].parameter, i, &arr);
			}
			arr.store_flag = false;
			arr.store_counter++;

		}
	}
	printf("\n");
	//Save last generation to file
	//for (i = 0; i < N; i++){
	//	fitnessTest2File(guys[i].parameter, i, &arr);
	//}


	for (i = 0; i < N; i++){
		for (j = 0; j < prop.n_params; j++){
			fprintf(data.data, "%8.5f	", guys[i].parameter[j]);
		}
		fprintf(data.data, "%8.5f	%8.5f	%d\n", guys[i].H, guys[i].t, i);
	}
	for (i = 0; i < N_best; i++){
		for (j = 0; j < prop.n_params; j++){
			fprintf(data.elitedata, "%8.5f	", bestguys[i].parameter[j]);
		}
		fprintf(data.elitedata, "%8.5f	%8.5f	%d\n", bestguys[i].H, bestguys[i].t, i);
	}
	printf("Simulation time:	%.3f seconds\n", arr.simulation_time);
	printf("Evolution time:		%.3f seconds\n", arr.evolution_time);

	free(arr.errorTh);
	free(arr.errordTh);
	free(arr.errorX);
	free(arr.speedX);
	//free(arr.errorTop);
	free(arr.loci);
	fclose(data.data);
	fclose(data.surfdata);
	fclose(data.elitedata);
	return 0;
}