#ifndef COMPARATIVEANALYSIS_H
#define COMPARATIVEANALYSIS_H

#include "C:\Documents and Settings\Work\Desktop\Biology_TOOLS\phylogeny.h"

//checked with results obtained with Martins' compare 4.6, with alpha = 0;
REGRESSION phyloLinearRegression(std::vector<double> xObs, std::vector<double> yObs, int numObs, Phylogeny aPhylogeny)
{
	ColumnVector Ones(numObs), Y(numObs), Y2(numObs), beta, SSres(1), SStotal(1), SSmodel(1); Y = 0.0; beta = 0.0; Ones = 1.0;
	Matrix X(numObs,2), X2(numObs,2), temp(4,4), varBeta,a,b,c; X2 = 0.0; X = 0.0; varBeta = 0.0, temp;
	IdentityMatrix I(numObs);
	double rSQR,pearsonProduct=0.0;
	REGRESSION results;

	//aPhylogeny = aPhylogeny.transform_phylogenyConstant();

	//aPhylogeny.transform_phylogenyPagel();
	Matrix V = aPhylogeny.format_phylogeny_to_matrix();

	//std::cout << std::endl << "Variance-covariance matrix:" << std::endl << std::endl;
	//V = aPhylogeny.matrix_standardize_UNQUAL(V);
	//V = aPhylogeny.matrix_polarize(V);
	//V = aPhylogeny.
	//std::cout << "standardized" << std::endl;
	//cout << setprecision(4) << V << endl;

	//V = I;

	/* float x = V(1,1), y = V(2,2), z = V(3,3), v = V(4,4);


	temp << x << 0 << 0 << 0 <<
	0 << y << 0 << 0 <<
	0 << 0 << z << 0 <<
	0 << 0 << 0 << v;

	V = V*temp.i();
	cout << V << endl;*/
	//V << 1 << 0 << 0 << 0 <<
	// 0 << 1 << 0.4714 << 0.4714 <<
	// 0 << 0.4714 << 1 << 0.83333 <<
	// 0 << 0.4714 << 0.83333 << 1;

	// V << 1 << 0 << 0 << 0
	//<< 0 << 1 << 0.4714 << 0.4714
	//<< 0 << 0.4714 << 1 << 0.8333
	//<< 0 << 0.4714 << 0.8333 << 1;

	/* V <<1.0000 << 148.4132 << 8886110.5205 << 1096.6332
	<< 148.4132 << 1.0000 << 3269017.3725 << 403.4288
	<< 8886110.5205 << 3269017.3725 << 1.0000 << 59874.1417
	<< 1096.6332 << 403.4288 << 59874.1417 << 1.0000;*/
	//cout << V << endl;

	for(int j = 1; j < numObs + 1; j++) {
		Y.Row(j) << yObs[j-1];
		X.Row(j) << 1 << xObs[j-1];
	}

	//cout << (X.t() * V.i() * X).i() << endl;

	//cout << X << V << Y << endl;

	beta = (X.t() * V.i() * X).i() * X.t() * V.i() * Y;


	//cout << beta << endl;

	//cout << beta << endl;

	//cout << "beta = " << beta << endl;
	SSmodel = beta.t() * X.t() * V.i() * Y - (Y.t() * V.i() * Ones * Y.t() * V.i() * Ones)/numObs;//- (Y.Sum() * Y.Sum()) / numObs;//Y.t() * V.i() * X * beta;// - (Y.Sum() * Y.Sum()) / numObs;
	SStotal = Y.t() * V.i() * Y - (Y.t() * V.i()* Ones * Y.t() * V.i() * Ones)/numObs;//- (Y.Sum() * Y.Sum()) / numObs;
	SSres = Y.t() * V.i() * Y - beta.t() * X.t() * V.i() * Y;
	varBeta = (X.t() * V.i() * X).i() * (SSres(1)/(numObs - 2)); // Martins & Hansen (1997) does not multiply with (SSres(1)/(numObs-2)), but all models SimpleLS & GgenerakLS show that you should

	//cout << "turd" << SStotal << " " << SSmodel(1)+SSres(1) << endl;

	rSQR = SSmodel(1)/SStotal(1);//1 - SSres(1)/SStotal(1);//

/*	cout << setprecision(5) << "__________________________(phylo_corrected)__ANOVA__" << endl << endl;
	cout << "Source " << setw(4) << "df" << setw(12) << "SS" << setw(12) << "MS" << setw(12) << "F (P)" << endl;
	cout << "____________________________________________________" << endl << endl;
	cout << "Regression" << setw(4) << 1 << setw(12) << SSmodel(1) << setw(12) << SSmodel(1) << setw(12) << SSmodel(1)/(SSres(1)/(numObs-2)) << endl;
	cout << "Error " << setw(4) << numObs-2 << setw(12) << SSres(1) << setw(12) << SSres(1)/(numObs-2) << setw(4) << "(" << setw(8) << F_Prob(1, numObs-2, SSmodel(1)/(SSres(1)/(numObs-2))) << ")" << endl;
	cout << "Total " << setw(4) << numObs-1 << setw(12) << SStotal(1) << endl;
	cout << "____________________________________________________" << endl;
	cout << "N = " << numObs << ", R = " << ((rSQR > 0) ? sqrt(rSQR) : -(sqrt(-rSQR)) ) << ", R\375 = " << ((rSQR > 0) ? rSQR : rSQR) << ", R\375(adj) = " << (1 - (SSres(1)/(numObs-2))/(SStotal(1)/(numObs-1))) << endl;
	cout << "y = " << beta(1) << "(" << sqrt(varBeta(1,1)) << ") + " << beta(2) << "(" << sqrt(varBeta(2,2)) << ")x" << endl;
	double AICcorrection = (2 * (numObs-2) * ((numObs-2) + 1)) / (numObs - (numObs - 2) - 1);
	cout << "AIC = " << (2*(numObs-2) + numObs * log(SSres(1)/numObs) + AICcorrection) << ", log-likelihood = " << (-numObs/2)*(1+log(2*PI)+log(SSres(1)/numObs)) << endl << endl;
	cout << "_____________(conservative phylo_corrected)__ANOVA__" << endl << endl;
	cout << "Number of polytomies in tree: " << aPhylogeny.stats_numberOfPolytomies() << endl;
	cout << "New d.f. adjusted for the # of polytomies: " << (numObs - aPhylogeny.stats_numberOfPolytomies() - 2) << endl;
	cout << "F = " << SSmodel(1)/(SSres(1)/(numObs-2)) << ", d.f. = (1," << (numObs - aPhylogeny.stats_numberOfPolytomies() - 2) << "), P = " << F_Prob(1, numObs-aPhylogeny.stats_numberOfPolytomies()-2, SSmodel(1)/(SSres(1)/(numObs-2))) << endl;
	cout << "____________________________________________________" << endl;
*/

	//cout << "AIC = " << numObs * (log((2 * PI * SSres(1))/numObs)+1);
	results.correlation_coeficient = ((rSQR == 0) ? 0.0 : pearsonProduct);
	results.r_square = rSQR;
	results.F_value = SSmodel(1)/(SSres(1)/(numObs-2));
	results.p_value = F_Prob(1, numObs-1, SSmodel(1)/(SSres(1)/(numObs-2)));
	results.df = numObs-1;
	results.alpha = beta(1);
	results.alpha_SD = sqrt(varBeta(1,1));
	results.beta = beta(2);
	results.beta_SD = sqrt(varBeta(2,2));
	results.N = numObs;
	results.PHYLO = (SSres(1)/(numObs-2))/SStotal(1);

	return results;
}

//checked with results obtained with Martins' compare 4.6, with alpha = 0;
REGRESSION linearRegression(std::vector<double> xObs, std::vector<double> yObs, int numObs)
{
	ColumnVector Ones(numObs), Y(numObs), Y2(numObs), beta, SSres(1), SStotal(1), SSmodel(1); Y = 0.0; beta = 0.0; Ones = 1.0;
	Matrix X(numObs,2), X2(numObs,2), temp(4,4), varBeta,a,b,c; X2 = 0.0; X = 0.0; varBeta = 0.0, temp;
	IdentityMatrix I(numObs);
	double rSQR,pearsonProduct=0.0;
	REGRESSION results;

	//aPhylogeny = aPhylogeny.transform_phylogenyConstant();

	//aPhylogeny.transform_phylogenyPagel();
	Matrix V;

	//std::cout << std::endl << "Variance-covariance matrix:" << std::endl << std::endl;
	V = I ;


	//V = aPhylogeny.matrix_polarize(V);
	//V = aPhylogeny.
	//std::cout << "standardized" << std::endl;
	//cout << V << endl;



	//V = I;

	/* float x = V(1,1), y = V(2,2), z = V(3,3), v = V(4,4);


	temp << x << 0 << 0 << 0 <<
	0 << y << 0 << 0 <<
	0 << 0 << z << 0 <<
	0 << 0 << 0 << v;

	V = V*temp.i();
	cout << V << endl;*/
	//V << 1 << 0 << 0 << 0 <<
	// 0 << 1 << 0.4714 << 0.4714 <<
	// 0 << 0.4714 << 1 << 0.83333 <<
	// 0 << 0.4714 << 0.83333 << 1;

	// V << 1 << 0 << 0 << 0
	//<< 0 << 1 << 0.4714 << 0.4714
	//<< 0 << 0.4714 << 1 << 0.8333
	//<< 0 << 0.4714 << 0.8333 << 1;

	/* V <<1.0000 << 148.4132 << 8886110.5205 << 1096.6332
	<< 148.4132 << 1.0000 << 3269017.3725 << 403.4288
	<< 8886110.5205 << 3269017.3725 << 1.0000 << 59874.1417
	<< 1096.6332 << 403.4288 << 59874.1417 << 1.0000;*/
	//cout << V << endl;

	for(int j = 1; j < numObs + 1; j++) {
		Y.Row(j) << yObs[j-1];
		X.Row(j) << 1 << xObs[j-1];
	}

	//cout << V << X << Y << endl;
	//cout << (X.t() * V.i() * X).i() << endl;

	//cout << V << Y << endl;

	beta = (X.t() * V.i() * X).i() * (X.t() * V.i() * Y);

	// cout << beta << endl;

	//cout << "beta = " << beta << endl;
	SSmodel = beta.t() * X.t() * V.i() * Y - (Y.t() * V.i() * Ones * Y.t() * V.i() * Ones)/numObs;//- (Y.Sum() * Y.Sum()) / numObs;//Y.t() * V.i() * X * beta;// - (Y.Sum() * Y.Sum()) / numObs;
	SStotal = Y.t() * V.i() * Y - (Y.t() * V.i()* Ones * Y.t() * V.i() * Ones)/numObs;//- (Y.Sum() * Y.Sum()) / numObs;
	SSres = Y.t() * V.i() * Y - beta.t() * X.t() * V.i() * Y;
	varBeta = (X.t() * V.i() * X).i() * (SSres(1)/(numObs - 2)); // Martins & Hansen (1997) does not multiply with (SSres(1)/(numObs-2)), but all models SimpleLS & GgenerakLS show that you should

	//cout << "turd" << SStotal << " " << SSmodel(1)+SSres(1) << endl;

	/* X.Column(2) = X.Column(2) - ((X.Column(2).Sum())/numObs);
	b = X.Column(2).t() * Y;
	pearsonProduct = b(1,1) / (sqrt(X.Column(2).SumSquare() * Y.SumSquare()));
	cout << pearsonProduct << endl;
	*/
	rSQR = SSmodel(1)/SStotal(1);//1 - SSres(1)/SStotal(1);//

	/*cout << "_____________________________________________ANOVA__" << endl << endl;
	cout << "Source " << setw(4) << "df" << setw(12) << "SS" << setw(12) << "MS" << setw(12) << "F (P)" << endl;
	cout << "____________________________________________________" << endl << endl;
	cout << "Regression" << setw(4) << 1 << setw(12) << SSmodel(1) << setw(12) << SSmodel(1) << setw(12) << SSmodel(1)/(SSres(1)/(numObs-2)) << endl;
	cout << "Error " << setw(4) << numObs-2 << setw(12) << SSres(1) << setw(12) << SSres(1)/(numObs-2) << setw(4) << "(" << setw(8) << F_Prob(1, numObs-2, SSmodel(1)/(SSres(1)/(numObs-2))) << ")" << endl;
	cout << "Total " << setw(4) << numObs-1 << setw(12) << SStotal(1) << endl;
	cout << "____________________________________________________" << endl;
	cout << "N = " << numObs << ", R = " << ((rSQR > 0) ? sqrt(rSQR) : -(sqrt(-rSQR)) ) << ", R\375 = " << ((rSQR > 0) ? rSQR : rSQR) << ", R\375(adj) = " << (1 - (SSres(1)/(numObs-2))/(SStotal(1)/(numObs-1))) << endl;
	cout << "y = " << beta(1) << "(" << sqrt(varBeta(1,1)) << ") + " << beta(2) << "(" << sqrt(varBeta(2,2)) << ")x" << endl;
	double AICcorrection = (2 * (numObs-2) * ((numObs-2) + 1)) / (numObs - (numObs - 2) - 1);
	cout << "AIC = " << (2*(numObs-2) + numObs * log(SSres(1)/numObs) + AICcorrection) << ", log-likelihood = " << (-numObs/2)*(1+log(2*PI)+log(SSres(1)/numObs)) << endl << endl;
	//cout << "AIC = " << numObs * (log((2 * PI * SSres(1))/numObs)+1);
 */
	results.correlation_coeficient = ((rSQR == 0) ? 0.0 : pearsonProduct);
	results.r_square = rSQR;
	results.F_value = SSmodel(1)/(SSres(1)/(numObs-2));
	results.p_value = F_Prob(1, numObs-1, SSmodel(1)/(SSres(1)/(numObs-2)));
	results.df = numObs-1;
	results.alpha = beta(1);
	results.alpha_SD = sqrt(varBeta(1,1));
	results.beta = beta(2);
	results.beta_SD = sqrt(varBeta(2,2));
	results.N = numObs;
	results.PHYLO = (SSres(1)/(numObs-2))/SStotal(1);

	return results;
}


REGRESSION oneWAYANOVA(std::vector< std::string > xObs, std::vector<double> yObs, int numObs, std::vector < std::string > theFactors)
{
	ColumnVector Ones(numObs), Y(numObs), Y2(numObs), beta, SSres(1), SStotal(1), SSmodel(1); Y = 0.0; beta = 0.0; Ones = 1.0;
	Matrix X, varBeta, design; X = 0.0; varBeta = 0.0;
	IdentityMatrix I(numObs);
	double rSQR,pearsonProduct=0.0;
	REGRESSION results;

	Matrix V;

	V = I ;
	int count = 0;
	for(int j = 1; j < numObs + 1; j++) Y.Row(j) << yObs[j-1];

	//design matrix

	std::string tempFactor = xObs[0];
	ColumnVector temp(numObs); temp = 1.0;
	design = temp;

	for(j = 0; j < theFactors.size(); j++) {
		for(int i = 1; i < numObs+1; i++) {
			if(xObs[i-1] == theFactors[j]) temp.Row(i) << 1.0;
			else temp.Row(i) << 0.0;
		}
		design = design | temp;
	}

	//cout << design << endl;

	X = design;

	beta = (X.t() * V.i() * X).i() * (X.t() * V.i() * Y);
	SSmodel = beta.t() * X.t() * V.i() * Y - (Y.t() * V.i() * Ones * Y.t() * V.i() * Ones)/numObs;//- (Y.Sum() * Y.Sum()) / numObs;//Y.t() * V.i() * X * beta;// - (Y.Sum() * Y.Sum()) / numObs;
	SStotal = Y.t() * V.i() * Y - (Y.t() * V.i()* Ones * Y.t() * V.i() * Ones)/numObs;//- (Y.Sum() * Y.Sum()) / numObs;
	SSres = Y.t() * V.i() * Y - beta.t() * X.t() * V.i() * Y;
	varBeta = (X.t() * V.i() * X).i() * (SSres(1)/(numObs - 2)); // Martins & Hansen (1997) does not multiply with (SSres(1)/(numObs-2)), but all models SimpleLS & GgenerakLS show that you should

	rSQR = SSmodel(1)/SStotal(1);//1 - SSres(1)/SStotal(1);//

	cout << "_______________________________________1-way_ANOVA__" << endl << endl;
	cout << "Source " << setw(4) << "df" << setw(12) << "SS" << setw(12) << "MS" << setw(12) << "F (P)" << endl;
	cout << "____________________________________________________" << endl << endl;
	cout << "Model " << setw(4) << (theFactors.size()-1) << setw(12) << SSmodel(1) << setw(12) << SSmodel(1)/(theFactors.size()-1) << setw(12) << (SSmodel(1)/(theFactors.size()-1))/(SSres(1)/(numObs-(theFactors.size()))) << endl;
	cout << "Error " << setw(4) << numObs-theFactors.size() << setw(12) << SSres(1) << setw(12) << SSres(1)/(numObs-theFactors.size()) << setw(4) << "(" << setw(8) << F_Prob(theFactors.size()-1, numObs-theFactors.size(), (SSmodel(1)/(theFactors.size()-1))/(SSres(1)/(numObs-(theFactors.size())))) << ")" << endl;
	cout << "Total " << setw(4) << numObs-1 << setw(12) << SStotal(1) << endl;
	cout << "____________________________________________________" << endl;
	cout << "N = " << numObs << ", R = " << ((rSQR > 0) ? sqrt(rSQR) : -(sqrt(-rSQR)) ) << ", R\375 = " << ((rSQR > 0) ? rSQR : rSQR) << ", R\375(adj) = " << (1 - (SSres(1)/(numObs-(theFactors.size())))/(SStotal(1)/(numObs-1))) << endl;
	// cout << "y = " << beta(1) << "(" << sqrt(varBeta(1,1)) << ") + " << beta(2) << "(" << sqrt(varBeta(2,2)) << ")x" << endl;
	// double AICcorrection = (2 * (numObs-2) * ((numObs-2) + 1)) / (numObs - (numObs - 2) - 1);
	// cout << "AIC = " << (2*(numObs-2) + numObs * log(SSres(1)/numObs) + AICcorrection) << ", log-likelihood = " << (-numObs/2)*(1+log(2*PI)+log(SSres(1)/numObs)) << endl << endl;

	results.correlation_coeficient = ((rSQR == 0) ? 0.0 : pearsonProduct);
	results.r_square = rSQR;
	results.F_value = (SSmodel(1)/(theFactors.size()-1))/(SSres(1)/(numObs-(theFactors.size())));
	results.p_value = F_Prob(theFactors.size()-1, numObs-theFactors.size(), (SSmodel(1)/(theFactors.size()-1))/(SSres(1)/(numObs-(theFactors.size()))));
	results.df = numObs-theFactors.size();
	results.alpha = beta(1);
	results.alpha_SD = sqrt(varBeta(1,1));
	results.beta = beta(2);
	results.beta_SD = sqrt(varBeta(2,2));
	results.N = numObs;
	results.PHYLO = (SSres(1)/(numObs-theFactors.size()))/SStotal(1);

	return results;
}

REGRESSION phylo_oneWAYANOVA(std::vector< std::string > xObs, std::vector<double> yObs, int numObs, std::vector < std::string > theFactors, Phylogeny aPhylogeny)
{
	ColumnVector Ones(numObs), Y(numObs), Y2(numObs), beta, SSres(1), SStotal(1), SSmodel(1); Y = 0.0; beta = 0.0; Ones = 1.0;
	Matrix X, varBeta, design; X = 0.0; varBeta = 0.0;
	IdentityMatrix I(numObs);
	double rSQR,pearsonProduct=0.0;
	REGRESSION results;

	Matrix V = aPhylogeny.format_phylogeny_to_matrix();
	V = aPhylogeny.matrix_standardize_UNQUAL(V);

	int count = 0;
	for(int j = 1; j < numObs + 1; j++) Y.Row(j) << yObs[j-1];

	//design matrix

	std::string tempFactor = xObs[0];
	ColumnVector temp(numObs); temp = 1.0;
	design = temp;

	for(j = 0; j < theFactors.size(); j++) {
		for(int i = 1; i < numObs+1; i++) {
			if(xObs[i-1] == theFactors[j]) temp.Row(i) << 1.0;
			else temp.Row(i) << 0.0;
		}
		design = design | temp;
	}

	//cout << design << endl;

	X = design;

	//cout << X << endl << V << endl;

	beta = (X.t() * V.i() * X).i() * (X.t() * V.i() * Y);
	SSmodel = beta.t() * X.t() * V.i() * Y - (Y.t() * V.i() * Ones * Y.t() * V.i() * Ones)/numObs;//- (Y.Sum() * Y.Sum()) / numObs;//Y.t() * V.i() * X * beta;// - (Y.Sum() * Y.Sum()) / numObs;
	SStotal = Y.t() * V.i() * Y - (Y.t() * V.i()* Ones * Y.t() * V.i() * Ones)/numObs;//- (Y.Sum() * Y.Sum()) / numObs;
	SSres = Y.t() * V.i() * Y - beta.t() * X.t() * V.i() * Y;
	varBeta = (X.t() * V.i() * X).i() * (SSres(1)/(numObs - 2)); // Martins & Hansen (1997) does not multiply with (SSres(1)/(numObs-2)), but all models SimpleLS & GgenerakLS show that you should

	rSQR = SSmodel(1)/SStotal(1);//1 - SSres(1)/SStotal(1);//

	cout << "_________________________________phylo_1-way_ANOVA__" << endl << endl;
	cout << "Source " << setw(4) << "df" << setw(12) << "SS" << setw(12) << "MS" << setw(12) << "F (P)" << endl;
	cout << "____________________________________________________" << endl << endl;
	cout << "Model " << setw(4) << (theFactors.size()-1) << setw(12) << SSmodel(1) << setw(12) << SSmodel(1)/(theFactors.size()-1) << setw(12) << (SSmodel(1)/(theFactors.size()-1))/(SSres(1)/(numObs-(theFactors.size()))) << endl;
	cout << "Error " << setw(4) << numObs-theFactors.size() << setw(12) << SSres(1) << setw(12) << SSres(1)/(numObs-theFactors.size()) << setw(4) << "(" << setw(8) << F_Prob(theFactors.size()-1, numObs-theFactors.size(), (SSmodel(1)/(theFactors.size()-1))/(SSres(1)/(numObs-(theFactors.size())))) << ")" << endl;
	cout << "Total " << setw(4) << numObs-1 << setw(12) << SStotal(1) << endl;
	cout << "____________________________________________________" << endl;
	cout << "N = " << numObs << ", R = " << ((rSQR > 0) ? sqrt(rSQR) : -(sqrt(-rSQR)) ) << ", R\375 = " << ((rSQR > 0) ? rSQR : rSQR) << ", R\375(adj) = " << (1 - (SSres(1)/(numObs-(theFactors.size())))/(SStotal(1)/(numObs-1))) << endl;
	// cout << "y = " << beta(1) << "(" << sqrt(varBeta(1,1)) << ") + " << beta(2) << "(" << sqrt(varBeta(2,2)) << ")x" << endl;
	// double AICcorrection = (2 * (numObs-2) * ((numObs-2) + 1)) / (numObs - (numObs - 2) - 1);
	// cout << "AIC = " << (2*(numObs-2) + numObs * log(SSres(1)/numObs) + AICcorrection) << ", log-likelihood = " << (-numObs/2)*(1+log(2*PI)+log(SSres(1)/numObs)) << endl << endl;

	results.correlation_coeficient = ((rSQR == 0) ? 0.0 : pearsonProduct);
	results.r_square = rSQR;
	results.F_value = (SSmodel(1)/(theFactors.size()-1))/(SSres(1)/(numObs-(theFactors.size())));
	results.p_value = F_Prob(theFactors.size()-1, numObs-theFactors.size(), (SSmodel(1)/(theFactors.size()-1))/(SSres(1)/(numObs-(theFactors.size()))));
	results.df = numObs-theFactors.size();
	results.alpha = beta(1);
	results.alpha_SD = sqrt(varBeta(1,1));
	results.beta = beta(2);
	results.beta_SD = sqrt(varBeta(2,2));
	results.N = numObs;
	results.PHYLO = (SSres(1)/(numObs-theFactors.size()))/SStotal(1);

	return results;
}

double theVariance(EFFECTSIZE d)
{
	return ((d.Ne + d.Nc)/(d.Ne * d.Nc)) + (d.effectSize * d.effectSize)/(2*(d.Ne + d.Nc));
}

double theCoVariance(EFFECTSIZE d1, EFFECTSIZE d2, double P)
{
	return ((d1.Ne + d1.Nc)/(d1.Ne * d1.Nc)) * P + (d1.effectSize * d2.effectSize * P * P)/(2*(d1.Ne + d1.Nc));
}

Matrix varianceMatrix(std::vector< EFFECTSIZE > d)
{
	Matrix theCovarianceMatrix(d.size(),d.size()); theCovarianceMatrix = 0.0;
	for(int i = 1; i <= theCovarianceMatrix.Ncols(); i++) theCovarianceMatrix(i,i) = theVariance(d[i-1]);
	return theCovarianceMatrix;
}


Matrix covarianceMatrix(std::vector< EFFECTSIZE > d, Phylogeny aPhylo, bool weighted)
{
	Matrix theCovarianceMatrix(d.size(),d.size()), P; theCovarianceMatrix = 0.0;
	P = aPhylo.matrix_standardize_UNQUAL(aPhylo.format_phylogeny_to_matrix());

	// matrixPRINT(P, "phylogeny matrix");

	for(int i = 1; i <= theCovarianceMatrix.Ncols(); i++)
	for(int j = i; j <= theCovarianceMatrix.Nrows(); j++) {
		if(i == j) {
			if(weighted == true)
			theCovarianceMatrix(i,j) = theVariance(d[i-1]);
			else theCovarianceMatrix(i,j) = 1.0;
		}
		else theCovarianceMatrix(i,j) = theCovarianceMatrix(j,i) = theCoVariance(d[i-1], d[j-1], P(i,j));
	}

	return theCovarianceMatrix;
}



//checked with results obtained with Martins' compare 4.6, with alpha = 0;
EFFECTSIZE pool_effects_phylo(std::vector< EFFECTSIZE > yObs, int numObs, Phylogeny aPhylogeny, bool weighted, bool phyloCorrect)
{
	ColumnVector Ones(numObs), X(numObs), Y(numObs), pooledEffect, Q, SSmodel, SStotal, SSres; Y = 0.0; pooledEffect = 0.0; X = 1.0; Ones = 1.0;
	Matrix varBeta, V, QE; varBeta = 0.0; V = 0.0;
	IdentityMatrix I(numObs);

	double pearsonProduct=0.0;
	EFFECTSIZE results;

	if(weighted == true) {
		if (phyloCorrect == true) V = covarianceMatrix(yObs, aPhylogeny, true);
		else V = varianceMatrix(yObs);
	}
	else {
		if (phyloCorrect == true) V = covarianceMatrix(yObs, aPhylogeny, false);
		else V = I;
	}

	// matrixPRINT(V, "covariance matrix");

	for(int j = 1; j < numObs + 1; j++) Y.Row(j) << (yObs[j-1]).effectSize;

	pooledEffect = (X.t() * V.i() * X).i() * X.t() * V.i() * Y;
	varBeta = (X.t() * V.i() * X).i();
	Q = Y.t() * V.i() * Y - (pooledEffect.t() * (X.t() * V.i() * X).i() * pooledEffect);

	SStotal = Y.t() * V.i() * Y - (Y.t() * V.i()* Ones * Y.t() * V.i() * Ones)/numObs;//- (Y.Sum() * Y.Sum()) / numObs;
	SSres = Y.t() * V.i() * Y - pooledEffect.t() * X.t() * V.i() * Y;
	SSmodel = pooledEffect.t() * X.t() * V.i() * Y - (Y.t() * V.i() * Ones * Y.t() * V.i() * Ones)/numObs;


	//QE = pooledEffect.t() * (X.t() * V.i() * X).i() * pooledEffect;
	//cout << "crap Q= " << QE(1,1) << " " << CHIsq_Prob(QE(1,1), 1) << endl;


	results.effectSize = pooledEffect(1);
	results.variance = varBeta(1,1);
	results.L95CI = results.effectSize - (1.96 * sqrt(results.variance));
	results.U95CI = results.effectSize + (1.96 * sqrt(results.variance));
	results.k = yObs.size();
	results.Q = Q(1);
	results.Q_df = results.k - 1;
	results.Q_p = CHIsq_Prob(results.Q, results.Q_df);
	results.covarianceMatrix = V;

	results.residual = SSmodel(1)/SStotal(1);
	cout << SSmodel(1) << "->" << SStotal(1) << "->" << SSres(1) << " " <<results.k << " " << results.residual << endl;

	cout << setprecision(4) << char(230) << " = " << results.effectSize << ", 95%CI(" << results.L95CI << ", " << results.U95CI << "), Q = " << results.Q << ", d.f. = " << results.Q_df << ", p = " << results.Q_p << endl;
	return results;
}

#endif