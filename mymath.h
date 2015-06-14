#define pi  3.141592653589793238462

void linspace(double array[], double xi, double xf, int n){
	int i;
	double dx, range = xf - xi;
	dx = range / (n - 1);
	array[0] = xi;
	for (i = 1; i < n; i++){
		array[i] = array[i - 1] + dx;
	}
}
void logspace(double array[], double xi, double xf, int n){
	int i;
	double s = log(xi);
	double f = log(xf);
	linspace(array, s, f, n);
	for (i = 0; i < n; i++){
		array[i] = exp(array[i]);
	}
}

int get_sign(double value)
{
	if (value == 0)
		return 0;
	else if (value > 0)
		return 1;
	else
		return -1;
}
void rndChoice(int randomArray[], int nn, int NN){
	//Choose nn integers from the set [0,NN-1] without repetition. Used the Knuth algorithm found on the interwebz (stack overflow). Ordo(Npmove).
	int in, im;
	int rn, rm;
	im = 0;
	for (in = 0; in < NN && im < nn; in++){
		rn = NN - in;
		rm = nn - im;
		if (genrand_real1()*rn < rm){
			randomArray[im] = in;
			im += 1;
		}
	}
}
int tri_randint(int a, int b){
	//Returns number between a and b with triangular distribution
	return (int)(a + pow(genrand_real1(), 0.5)*(b - a));
}
int tri_inv_randint(int a, int b){
	//Returns number between a and b with triangular distribution, sloping downwards
	return (int)(b - pow(1 - genrand_real1(), 0.5)*(b - a));
}
double mean(double array[], int len){
	double sum = 0;
	int i;
	for (i = 0; i < len; i++){
		sum += array[i];
	}
	return sum / len;
}
int mod(int a, int b){
	if (b < 0) 
		return mod(-a, -b);
	
	int ret = a % b;
	if (ret < 0) 
		ret += b;
		return ret;
	
}
int isvalueinarray(int val, int *arr, int size){
	int i;
	for (i = 0; i < size; i++) {
		if (arr[i] == val)
			return 1;
	}
	return 0;
}


int heaviside(double x){
	if (x > 0){
		return 1;
	}
	else{
		return 0;
	}
}


double var(double a[], int n) {
	if (n == 0){
		return 0.0;
	}
	double sum = 0;
	for (int i = 0; i < n; i++){
		sum += a[i];
	}
	double mean = sum / n;
	double sq_diff_sum = 0;
	for (int i = 0; i < n; ++i) {
		double diff = a[i] - mean;
		sq_diff_sum += diff * diff;
	}
	double variance = sq_diff_sum / n;
	return variance;
}
//Autocorrelation. You must free acf after calling it
double * autocorr(double a[], int len, int maxlag){
	int i, k;
	double *acf= (double*)malloc(maxlag * sizeof(double));
	double term1;
	double term21;
	double term22;
	double norm;

	for (k = 0; k < maxlag ; k++){
		term1 = 0;
		term21 = 0;
		term22 = 0;
		for (i = 0; i < len - k; i++){
			term1  += a[i]*a[i + k];
			term21 += a[i];
			term22 += a[i + k];
		}
		term1  /= i;
		term21 /= i;
		term22 /= i;
		acf[k] = term1 - term21*term22;
	}
	norm = acf[0];
	for (k = 0; k < maxlag; k++) {
		acf[k] = acf[k] /norm;
	}
	return acf;
}

double tau(double a[], int len, int maxlag){
	int i;
	double sum = 0;
	double *acf = autocorr(a,len,maxlag);

	for (i = 0; i < maxlag; i++){
		sum += acf[i]* (1 - i/len);
	}
	free(acf);
	return 0.5 + sum;
}

int find_index(int a[], int n, int value) {
	int i;
	for (i = 0; i < n; i++)
	{
		if (a[i] == value)
		{
			return i;  /* it was found */
		}
	}
	return(-1);  /* if it was not found */
}

