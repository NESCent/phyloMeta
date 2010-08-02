// Statistics.h
#ifndef STATISTICS_H
#define STATISTICS_H

#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30
#define PI 3.14159265

struct BIAS_CI
{
	double L_CI;
	double U_CI;
};


struct T_TEST
{
	double t_value;
	double p_value;
	int df;
};

struct MOMENT
{
	int N;
	double mean;
	double SD;
	double AD; //average deviation
	double variance;
	double skewness;
	double kurtosis;
};

struct REGRESSION
{
	double correlation_coeficient;
	double r_square;
	double F_value;
	int df;
	double p_value;
	double alpha;
	double alpha_SD;
	double beta;
	double beta_SD;
	int N;
	double PHYLO;
};

struct EFFECTSIZE
{
	Matrix covarianceMatrix;
	double effectSize;
	double effect;
	int N;
	int Nc;
	int Ne;
	int k;
	double variance;
	double L95CI;
	double U95CI;
	double Q;
	int Q_df;
	double Q_p;
	double residual;
};


bool isPositiveDefinate(const Matrix & aMatrix)
{
	bool posDef = true;
	Matrix temp; temp = aMatrix;
	while(!temp.IsZero()) {
		if(temp.Determinant() < 0) posDef = false;
		temp = temp.SubMatrix(1, temp.Nrows()-1, 1, temp.Ncols() -1);
	}
	return posDef;
}




/*****************************************************************************
Functions used to estimate probabilities from F, t and ChiSquare distributions
*****************************************************************************/

float gammln(float xx) //Returns the value ln[gamma(xx)] for xx > 0.
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;

	return -tmp+log(2.5066282746310005*ser/x);
}

void gcf(float *gammcf, float a, float x, float *gln) // Returns the incomplete gamma function Q(a, x) evaluated by its continued fraction representation as gammcf. Also returns lnÃ(a) as gln.
{
	int i;
	float an,b,c,d,del,h;
	*gln=gammln(a);
	b=x+1.0-a; //Set up for evaluating continued fraction by modified Lentz’s method (§5.2) with b0 = 0.
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1; i<= ITMAX;i++) { //Iterate to convergence.
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) nrerror("a too large, ITMAX too small in gcf");
	*gammcf=exp(-x+a*log(x)-(*gln))*h; //Put factors in front.
}

void gser(float *gamser, float a, float x, float *gln) // Returns the incomplete gamma function P(a, x) evaluated by its series representation as gamser. // Also returns ln Ã(a) as gln.
{
	int n;
	float sum,del,ap;
	*gln=gammln(a);
	if (x <= 0.0) {
	if (x < 0.0) nrerror("x less than 0 in routine gser");
		*gamser=0.0;
		return;
	}
	else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		nrerror("a too large, ITMAX too small in routine gser");
		return;
	}
}

float gammp(float a, float x) //Returns the incomplete gamma function P(a, x).
{
	float gamser,gammcf,gln;
	if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammp");
	if (x < (a+1.0)) { //Use the series representation.
		gser(&gamser,a,x,&gln);
		return gamser;
	}
	else { //Use the continued fraction representation
		gcf(&gammcf,a,x,&gln);
		return 1.0-gammcf; //and take its complement.
	}
}

float gammq(float a, float x) //Returns the incomplete gamma function Q(a, x) = 1 - P(a, x).
{
	float gamser,gammcf,gln;
	if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammq");
	if (x < (a+1.0)) { //Use the series representation
		gser(&gamser,a,x,&gln);
		return 1.0-gamser; //and take its complement.
	}
	else { //Use the continued fraction representation.
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}

float betacf(float a, float b, float x)
{
	int m,m2;
	float aa,c,d,del,h,qab,qam,qap;
	qab=a+b;
	qap=a+1.0;
	qam=a-1.0;
	c=1.0;
	d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN) d=FPMIN;
	d=1.0/d;
	h=d;
	for (m=1;m<=ITMAX;m++) {
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}

	if (m > ITMAX) nrerror("a or b too big, or MAXIT too small in betacf function");

	return h;
}


float betai(float a, float b, float x)
{
	float bt;
	if (x < 0.0 || x > 1.0) nrerror("Bad x in routine betai function");
	if (x == 0.0 || x == 1.0) bt=0.0;
	else bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
	if (x < (a+1.0)/(a+b+2.0)) return bt*betacf(a,b,x)/a;
	else return 1.0-bt*betacf(b,a,1.0-x)/b;
}


/***********************************************************************
Functions to calculate p-values of F, t and ChiSquare values
************************************************************************/

float F_Prob(int df1, int df2, float f)
{
	return betai(df2 * 0.5, df1 * 0.5, df2 / (df2 + df1 * f));
}

float CHIsq_Prob(double chsq, double df)
{
	if(chsq < 0.0) chsq = 0.0;
	else chsq = chsq/2.0;
	return gammq(df/2.0,chsq);
}

double t_Prob(double t, int dof)
#define ONEBYPI ((double) 0.3183098861837906715377675) /* Abramowitz & Stegun, of course. */
{
	/* Arguments:
	t : wanna evaluate integral from t to +infinity. dof: degrees of freedom of this t distribution.
	author http://lib.stat.cmu.edu/general/gaut.c
	This functional interface is not exactly what AS 3 used to be.
	*/
	double d_dof, s, c, f, a, b, probability;
	int fk, ks, im2, ioe, k;

	if (dof < 1) nrerror("tprob function is mishandled");
	d_dof = (double) dof; /* d_dof is F of fortran code. */

	a = t / sqrt(d_dof);
	b = d_dof/(d_dof + (t*t));
	im2 = dof - 2;
	ioe = dof % 2;
	s = c = f = 1.0;
	fk = ks = 2 + ioe;
	if (im2 > 2)
	for (k=ks; k<=im2; k+=2) {
		c = c*b*(fk - 1.0)/fk;
		s += c;
		if (s == f) break; /* == ? */
		f = s;
		fk += 2.0;
	}
	if (ioe != 1) { /* Label 20 of fortran code. */
		probability = 0.5 + (0.5*a*sqrt(b)*s);
		return probability;
	}
	else { /* Label 30 of fortran code. */
	if (dof == 1) s = 0.0;
		probability = 0.5 + ((a*b*s + atan(a))*ONEBYPI);
		return probability;
	}
}

//not
double Ztest(double Zval)
{
	return (1 - (1/sqrt(2*PI))*exp(-(Zval*Zval)/2)) * 2;
}



/***********************************************************************
General moments for summerizing data, mean, SE, etc.
************************************************************************/

//skewness and kurtosis wrong
MOMENT moment(std::vector<double> aVector)
{
	MOMENT results;
	results.mean = results.AD = results.SD = results.variance = results.skewness = results.kurtosis = 0.0;
	results.N = aVector.size();
	double ep = 0.0, s = 0.0, p = 0.0;

	if (results.N <= 1) nrerror("n must be at least 2 in moment");

	for (int j = 0; j < results.N; j++) s += aVector[j];
	results.mean = s/results.N;

	for (j = 0; j < results.N; j++) {
		results.AD += fabs(s = aVector[j] - results.mean);
		ep += s;
		results.variance += (p = s * s);
		results.skewness += (p *= s);
		results.kurtosis += (p *= s);
	}
	results.AD /= results.N;
	results.variance = (results.variance - ep * (ep/results.N))/(results.N - 1);
	results.SD = sqrt(results.variance);
	if (results.variance) {
		results.skewness /= (results.N * results.variance * results.SD);
		results.kurtosis = (results.kurtosis)/(results.N * results.variance * results.variance) - 3.0;
	}
	else nrerror("No skew/kurtosis when variance = 0 (in moment)");

	std::cout << "mean = " << results.mean << ", SD = " << results.SD << ", variance = " << results.variance
	<< ", skewness = " << results.skewness << ", kurtosis = " << results.kurtosis
	<< ", N = " << results.N << std::endl;

	return results;
}

double vectorMean(std::vector<double> aVector)
{
	double mean = 0;
	for(int i = 0; i < aVector.size(); i++) mean += aVector[i];
	return mean/aVector.size();
}

double vectorVariance(std::vector<double> aVector)
{
	double variance = 0, mean = vectorMean(aVector);
	for(int i = 0; i < aVector.size(); i++) variance += ((aVector[i] - mean)*(aVector[i] - mean));
	return variance / (aVector.size() - 1);
}

double vectorSD(std::vector<double> aVector)
{
	double variance = 0, mean = vectorMean(aVector);
	for(int i = 0; i < aVector.size(); i++) variance += ((aVector[i] - mean)*(aVector[i] - mean));
	return variance / (aVector.size() - 1);


}


/***********************************************************************
General statistics, t-test, z-test, regression, etc.
************************************************************************/

void pT_TEST(T_TEST aResult)
{
	cout << "t = " << aResult.t_value << ", d.f. = " << aResult.df << ", p = " << aResult.p_value << endl;
}

T_TEST ttest(std::vector<double> aVector1, std::vector<double> aVector2)
{
	double ave1, ave2, var1, var2, svar;
	T_TEST temp;

	ave1 = vectorMean(aVector1);
	ave2 = vectorMean(aVector2);
	var1 = vectorVariance(aVector1);
	var2 = vectorVariance(aVector2);

	temp.df = aVector1.size()+aVector2.size()-2;
	svar = ((aVector1.size()-1)*var1+(aVector2.size()-1)*var2)/temp.df;
	temp.t_value = (ave1-ave2)/sqrt(svar*(1.0/aVector1.size()+1.0/aVector2.size()));
	temp.p_value = betai(0.5*temp.df,0.5,temp.df/(temp.df+(temp.t_value)*(temp.t_value)));

	//pT_TEST(temp);

	return temp;
}

//checked with results obtained with systat;
REGRESSION simpleLinearRegression(double xObs[], double yObs[], int numObs)
{
	ColumnVector Ones(numObs), Y(numObs), Y2, beta(2), SSres(1), SStotal(1), SSmodel(1); Y2=0.0; Y = 0.0; beta = 0.0; Ones = 1.0;
	Matrix X(numObs,2), X2(numObs,2), varBeta, b, num, num2; X = 0.0; varBeta = 0.0;
	double rSQR = 0.0, pearsonProduct = 0.0 ;
	REGRESSION results;

	for(int j = 1; j < numObs + 1; j++) {
		//xObs[j-1] = 1/xObs[j-1]; //only for inverse polinomial regression
		Y.Row(j)<< yObs[j-1];
		X.Row(j) << 1 << xObs[j-1];
		X2.Row(j) << 1.0 << xObs[j-1];
	}

	beta = (X.t() * X).i() * (X.t() * Y);
	SSmodel = beta.t() * X.t() * Y - (Y.Sum() * Y.Sum()) / numObs;
	SSres = (Y - X * beta).t() * (Y - X * beta); // also Y.t() * (I - X * (X.t() * X).i() * X.t()) *Y
	SStotal = Y.t() * Y - (Y.Sum() * Y.Sum()) / numObs;
	varBeta = (X.t() * X).i() * (SSres(1)/(numObs-2));
	rSQR = SSmodel(1) / SStotal(1);

	cout << "correction = " << (Y.Sum() * Y.Sum()) / numObs << " " << Y.t() * Ones << endl;//(Y.t()* Y * Ones * Ones)/numObs << endl;
	//pearson product moment coeficient
	Y2 = Y - (Y.Sum())/numObs;
	X2.Column(2) = X2.Column(2) - ((X2.Column(2).Sum())/numObs);
	b = X2.Column(2).t() * Y2;
	pearsonProduct = b(1,1) / (sqrt(X2.Column(2).SumSquare() * Y2.SumSquare()));
	num = ((Y.t() * Y).i());
	num2 = Y.t() * X * (X.t() * X).i() * X.t() * Y;
	// cout << "-> " << num << " " << num2 << " " << num2(1,1) / num(1,1) << " " << Y.t()*(I-X*(X.t()*X).i()*X.t())*Y << endl;
	/*cout << "_____________________________________________ANOVA__" << endl << endl;
	cout << "Source " << setw(4) << "df" << setw(12) << "SS" << setw(12) << "MS" << setw(12) << "F (P)" << endl;
	cout << "____________________________________________________" << endl << endl;
	cout << "Regression" << setw(4) << 1 << setw(12) << SSmodel(1) << setw(12) << SSmodel(1) << setw(12) << SSmodel(1)/(SSres(1)/(numObs-2)) << endl;
	cout << "Error " << setw(4) << numObs-2 << setw(12) << SSres(1) << setw(12) << SSres(1)/(numObs-2) << setw(4) << "(" << setw(8) << F_Prob(1, numObs-2, SSmodel(1)/(SSres(1)/(numObs-2))) << ")" << endl;
	cout << "Total " << setw(4) << numObs-1 << setw(12) << SStotal(1) << endl;
	cout << "____________________________________________________" << endl;
	cout << "N=" << numObs << ", Rp=" << pearsonProduct << ", R=" << sqrt(rSQR) << ", R\375=" << rSQR << ", R\375(adj)=" << (((1 - (SSres(1)/(numObs-2))/(SStotal(1)/(numObs-1))) <= 0.0) ? 0 : (1 - (SSres(1)/(numObs-2))/(SStotal(1)/(numObs-1)))) << endl;
	cout << "y = " << beta(1) << "(" << sqrt(varBeta(1,1)) << ") + " << beta(2) << "(" <<sqrt(varBeta(2,2)) << ")x" << endl << endl;
	//cout << "shit = " << ((pearsonProduct*pearsonProduct)/1)/((1-pearsonProduct*pearsonProduct)/(numObs-1-1)) << endl;
*/
	results.correlation_coeficient = pearsonProduct;
	results.r_square = rSQR;
	results.F_value = SSres(1)/(numObs-2);
	results.p_value = F_Prob(1, numObs-2, SSmodel(1)/(SSres(1)/(numObs-2)));
	results.df = numObs-1;
	results.alpha = beta(1);
	results.alpha_SD = sqrt(varBeta(1,1));
	results.beta = beta(2);
	results.beta_SD = sqrt(varBeta(2,2));
	results.N = numObs;

	return results;
}

REGRESSION simpleLinearRegression(std::vector<double> xObs, std::vector<double> yObs, int numObs)
{
	ColumnVector Ones(numObs), Y(numObs), Y2, beta(1), SSres(1), SStotal(1), SSmodel(1); Y2=0.0; Y = 0.0; beta = 0.0; Ones = 1.0;
	Matrix X(numObs,2), X2(numObs,1), varBeta, b, num, num2; X = 0.0; varBeta = 0.0;
	double rSQR = 0.0, pearsonProduct = 0.0 ;
	REGRESSION results;

	for(int j = 1; j < numObs + 1; j++) {
		//xObs[j-1] = 1/xObs[j-1]; //only for inverse polinomial regression
		Y.Row(j)<< yObs[j-1];
		X.Row(j) << 1.0 << xObs[j-1];
		//X2.Row(j) << 1.0 << xObs[j-1];
	}

	beta = (X.t() * X).i() * (X.t() * Y);

	// cout << beta << endl;

	SSmodel = beta.t() * X.t() * Y - (Y.Sum() * Y.Sum()) / numObs;
	SSres = (Y - X * beta).t() * (Y - X * beta); // also Y.t() * (I - X * (X.t() * X).i() * X.t()) *Y
	SStotal = Y.t() * Y - (Y.Sum() * Y.Sum()) / numObs;
	varBeta = (X.t() * X).i() * (SSres(1)/(numObs-2));
	rSQR = SSmodel(1) / SStotal(1);

	// cout << SSmodel << SSres << SStotal << varBeta << rSQR << endl;

	//cout << "correction = " << (Y.Sum() * Y.Sum()) / numObs << " " << Y.t() * Ones << endl;//(Y.t()* Y * Ones * Ones)/numObs << endl;
	//pearson product moment coeficient
	// Y2 = Y - (Y.Sum())/numObs;
	// X2.Column(2) = X2.Column(2) - ((X2.Column(2).Sum())/numObs);
	// b = X2.Column(2).t() * Y2;
	// pearsonProduct = b(1,1) / (sqrt(X2.Column(2).SumSquare() * Y2.SumSquare()));
	// num = ((Y.t() * Y).i());
	// num2 = Y.t() * X * (X.t() * X).i() * X.t() * Y;
	// cout << "-> " << num << " " << num2 << " " << num2(1,1) / num(1,1) << " " << Y.t()*(I-X*(X.t()*X).i()*X.t())*Y << endl;
/*	cout << "_____________________________________________ANOVA__" << endl << endl;
	cout << "Source " << setw(4) << "df" << setw(12) << "SS" << setw(12) << "MS" << setw(12) << "F (P)" << endl;
	cout << "____________________________________________________" << endl << endl;
	cout << "Regression" << setw(4) << 1 << setw(12) << abs(SSmodel(1)) << setw(12) << abs(SSmodel(1)) << setw(12) << abs(SSmodel(1))/(abs(SSres(1))/(numObs-2)) << endl;
	cout << "Error " << setw(4) << numObs-2 << setw(12) << SSres(1) << setw(12) << SSres(1)/(numObs-2) << setw(4) << "(" << setw(8) << F_Prob(1, numObs-2, SSmodel(1)/(SSres(1)/(numObs-2))) << ")" << endl;
	cout << "Total " << setw(4) << numObs-1 << setw(12) << SStotal(1) << endl;
	cout << "____________________________________________________" << endl;
	cout << "N=" << numObs << ", Rp=" << pearsonProduct << ", R=" << sqrt(rSQR) << ", R\375=" << rSQR << ", R\375(adj)=" << (((1 - (SSres(1)/(numObs-2))/(SStotal(1)/(numObs-1))) <= 0.0) ? 0 : (1 - (SSres(1)/(numObs-2))/(SStotal(1)/(numObs-1)))) << endl;
	cout << "y = " << beta(1); //<< "(" << sqrt(varBeta(1)) << ")x" << endl << endl;
	cout << "shit = " << ((pearsonProduct*pearsonProduct)/1)/((1-pearsonProduct*pearsonProduct)/(numObs-1-1)) << endl;
*/
	results.correlation_coeficient = ((rSQR == 0) ? 0.0: pearsonProduct);
	results.r_square = rSQR;
	results.F_value = SSmodel(1)/(SSres(1)/(numObs-2));
	results.p_value = F_Prob(1, numObs-2, SSmodel(1)/(SSres(1)/(numObs-2))); 
	results.df = numObs-1;
	results.alpha = beta(1);
	results.alpha_SD = sqrt(varBeta(1,1));
	results.beta = beta(2);
	results.beta_SD = sqrt(varBeta(2,2));
	results.N = numObs;
	
	return results;
}


double stdnormal_cdf(double u)
{

/*
 * The standard normal CDF, for one random variable.
 *
 *   Author:  W. J. Cody
 *   URL:   http://www.netlib.org/specfun/erf
 *
 * This is the erfc() routine only, adapted by the
 * transform stdnormal_cdf(u)=(erfc(-u/sqrt(2))/2;
 */
 const double a[5] = {
  1.161110663653770e-002,3.951404679838207e-001,2.846603853776254e+001,
  1.887426188426510e+002,3.209377589138469e+003
 };
 const double b[5] = {
  1.767766952966369e-001,8.344316438579620e+000,1.725514762600375e+002,
  1.813893686502485e+003,8.044716608901563e+003
 };
 const double c[9] = {
  2.15311535474403846e-8,5.64188496988670089e-1,8.88314979438837594e00,
  6.61191906371416295e01,2.98635138197400131e02,8.81952221241769090e02,
  1.71204761263407058e03,2.05107837782607147e03,1.23033935479799725E03
 };
 const double d[9] = {
  1.00000000000000000e00,1.57449261107098347e01,1.17693950891312499e02,
  5.37181101862009858e02,1.62138957456669019e03,3.29079923573345963e03,
  4.36261909014324716e03,3.43936767414372164e03,1.23033935480374942e03
 };
 const double p[6] = {
  1.63153871373020978e-2,3.05326634961232344e-1,3.60344899949804439e-1,
  1.25781726111229246e-1,1.60837851487422766e-2,6.58749161529837803e-4
 };
 const double q[6] = {
  1.00000000000000000e00,2.56852019228982242e00,1.87295284992346047e00,
  5.27905102951428412e-1,6.05183413124413191e-2,2.33520497626869185e-3
 };
 register double y, z;

 if (_isnan(u))
  return _Nan._D;
 if (!_finite(u))
  return (u < 0 ? 0.0 : 1.0);
 y = fabs(u);
    if (y <= 0.46875*M_SQRT2) {
  /* evaluate erf() for |u| <= sqrt(2)*0.46875 */
  z = y*y;
  y = u*((((a[0]*z+a[1])*z+a[2])*z+a[3])*z+a[4])
       /((((b[0]*z+b[1])*z+b[2])*z+b[3])*z+b[4]);
  return 0.5+y;
 }
 z = exp(-y*y/2)/2;
 if (y <= 4.0) {
  /* evaluate erfc() for sqrt(2)*0.46875 <= |u| <= sqrt(2)*4.0 */
  y = y/M_SQRT2;
  y =
((((((((c[0]*y+c[1])*y+c[2])*y+c[3])*y+c[4])*y+c[5])*y+c[6])*y+c[7])*y+c[8])


/((((((((d[0]*y+d[1])*y+d[2])*y+d[3])*y+d[4])*y+d[5])*y+d[6])*y+d[7])*y+d[8]);

  y = z*y;
    } else {
  /* evaluate erfc() for |u| > sqrt(2)*4.0 */
  z = z*M_SQRT2/y;
  y = 2/(y*y);
        y = y*(((((p[0]*y+p[1])*y+p[2])*y+p[3])*y+p[4])*y+p[5])
    /(((((q[0]*y+q[1])*y+q[2])*y+q[3])*y+q[4])*y+q[5]);
        y = z*(M_1_SQRTPI-y);
    }
 return (u < 0.0 ? y : 1-y);
};

double cndev(double u)
{
	/* returns the inverse of cumulative normal distribution function
	Reference> The Full Monte, by Boris Moro, Union Bank of Switzerland
							RISK 1995(2)*/

	static double a[4]={2.50662823884,
											-18.61500062529,
											41.39119773534,
											-25.44106049637};
	static double b[4]={-8.47351093090,
											23.08336743743,
											-21.06224101826,
											3.13082909833};
	static double c[9]={0.3374754822726147,
											0.9761690190917186,
											0.1607979714918209,
											0.0276438810333863,
											0.0038405729373609,
											0.0003951896511919,
											0.0000321767881768,
											0.0000002888167364,
											0.0000003960315187};
	double x,r;
	x=u-0.5;
	if (fabs(x)<0.42)
	{
	r=x*x;
	r=x*(((a[3]*r+a[2])*r+a[1])*r+a[0])/((((b[3]*r+b[2])*r+b[1])*r+b[0])*r+1.0);
	return(r);
	}
	r=u;
	if(x>0.0) r=1.0-u;
			r=log(-log(r));
			r=c[0]+r*(c[1]+r*(c[2]+r*(c[3]+r*(c[4]+r*(c[5]+r*(c[6]+r*(c[7]+r*c[8])))))));
			if(x<0.0) r=-r;
	return(r);
}


BIAS_CI biasCorrectedBootstrap(std::vector<double> aVector, double alpha) 
{
	/*bias corrected CI's, count = # (reps < real value)
	zscore = cndev(count/straps);
	pl = stdnormal_cdf(2*zscore-1.96); pu = stdnormal_cdf(2*zscore+1.96);
	sort2(allData, straps);
	U = ceil((straps+1)*pl); L = floor((straps+1)*pu);
	RRCI[0] = allData[U]; RRCI[1] = allData[L];*/

	double zL = cndev(alpha/2), zU = cndev(1-alpha/2);
	
	BIAS_CI ci;
	double mean = 0, zscore = 0, pl = 0, pu = 0, count = 0, U = 0, L = 0;


	for(std::vector<double>::iterator i = aVector.begin(); i != aVector.end(); i++) mean+= *i;
	mean = mean/aVector.size();
	//cout << "mean:" << mean << endl;

	std::sort(aVector.begin(), aVector.end());

	
	//use the median rather than mean for centering the 95%CI
	/*double temp = ceil((aVector.size()-1)/2);
	temp = aVector[temp]; */

	for(i = aVector.begin(); i != aVector.end(); i++) 
		if(*i < mean) count++;// cout << *i << endl;}

	//cout << count << endl;
	
	//cout << "inZ " << (count/aVector.size()) << endl;
	zscore = cndev(count/aVector.size());
	//cout << "Z " << zscore << " " << zL <<  " " << zU << endl;
	pl = stdnormal_cdf(zL + 2*zscore); pu = stdnormal_cdf(zU + 2*zscore);
	//cout << pl << "-" << pu << endl;
	//cout << aVector.size()*pl << "-" << aVector.size()*pu << endl;
	
	//for(i = aVector.begin(); i != aVector.end(); i++) cout << *i << " ";
	//cout << endl;
	
	//for( i = aVector.begin(); i != aVector.end(); i++) cout << *i << " ";
	//cout << endl;
	//cout << zscore << " " << pl << " " << pu << endl;



	U = ceil((aVector.size()-1)*pu); L = floor((aVector.size()-1)*pl);

	ci.L_CI = aVector[L]; ci.U_CI = aVector[U];

	//cout << mean << " " << L << "-" << U << "->" << ci.L_CI << " " << ci.U_CI << endl;
	return ci;
}


double Zscore(double Zval)
{
	return (2.0 * (1.0 - stdnormal_cdf(Zval)));
}

#endif