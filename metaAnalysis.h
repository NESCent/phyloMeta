// metaAnalysis.h

#ifndef METAANALYSIS_H
#define METAANALYSIS_H


/*struct META_D
{
	double d;
	double var;
	double CI_L95;
	double CI_U95;
	double SE;
	double SD;
	double 
	double Q;
	double Q_df;
	double Q_p;
	double k;
}

*/


struct META_REGR            
{
	double r;
	double CI_L95;
	double CI_U95;
	double SE;
	double SD;
	double rp;
	double p;
	double Q;
	double Q_df;
	double Q_p;
	double k;
};

struct META_RESULTS            
{
	double e;
	double var;
	double k;
	double CI_L95;
	double CI_U95;
	
	double Qzero;
	double Qzero_df;
	double Qzero_p;

	double Q;
	double Q_df;
	double Q_p;
};

double Z_transform(double r) 
{
	return 0.5 * log((1.0+r)/(1.0-r));
}


META_REGR meta_correlation_coefficients_Hunter_Schmidt(std::vector <double> rList, std::vector<int> NList, bool print)
{
	META_REGR result;
	double sumN = 0.0, sumRbyN = 0.0, sumA = 0.0;

	result.k = rList.size(); result.Q = 0.0;

	for(int i = 0; i < result.k; i++) {
		sumN += NList[i];
		sumRbyN += (rList[i] * NList[i]);
	}

	result.r = sumRbyN/sumN;

	for(i = 0; i < result.k; i++) {
		sumA += (NList[i] * (rList[i] - result.r) * (rList[i] - result.r));
		result.Q += ((NList[i] - 1) * (rList[i] - result.r) * (rList[i] - result.r))/((1 - result.r * result.r) * (1 - result.r * result.r));
	}

	result.SD = sqrt(sumA/sumN);
	result.SE = result.SD/sqrt(result.k);

	result.CI_L95 = result.r - 1.96/(sqrt(sumN - 3 * result.k));
	result.CI_U95 = result.r + 1.96/(sqrt(sumN - 3 * result.k));

	result.Q_df = result.k - 1;

	result.Q_p = CHIsq_Prob(result.Q, result.Q_df);

/*	result.rp = result.r/result.SE;
	result.p = Ztest(result.rp); 
	std::cout << "Z test =" << result.rp << ", p = " << result.p << std::endl;*/
	
	if(print) {
		std::cout << "___________________________________________________________Hunter-Schmidt Method_____" << std::endl << std::endl;
		std::cout << "mean r = " << result.r << ", k = " << result.k << ", SE = " << result.SE << ", SD = " << result.SD << std::endl;
		std::cout << "95% CI lower = " << result.CI_L95 << ", 95% CI upper = " << result.CI_U95 << std::endl;
		std::cout << "homogeneity test: Q = " << result.Q << ", d.f. = " << result.Q_df <<", p = " << result.Q_p << ((result.Q_p <= 0.05) ? " | Significant heterogeneity!" : " | The data is homogeneous.") << std::endl;
		std::cout << "_____________________________________________________________________________________" << std::endl;
	}

	return result;
}


META_REGR meta_correlation_coefficients_Hedges_Olkin_fixed(std::vector <double> rList, std::vector<int> NList, bool print, bool simplePrint)
{
	META_REGR result;
	std::vector<double> zList, wList;

	double sumN = 0.0, sumW = 0.0, sumZbyN = 0.0;

	result.k = rList.size(); result.Q = 0.0;

	for(int i = 0; i < result.k; i++) {
		zList.push_back(Z_transform(rList[i]));
		wList.push_back(NList[i] - 3);
	}

	for(i = 0; i < result.k; i++) {
		sumN += NList[i];
		sumW += wList[i];
		sumZbyN += (zList[i] * wList[i]);
	}
	
	result.r = sumZbyN/sumW;

	for(i = 0; i < result.k; i++) result.Q += wList[i] * ((zList[i] - result.r) * (zList[i] - result.r));

	result.SE = 1/sqrt(sumW);

	result.CI_L95 = result.r - 1.96/(sqrt(sumN - 3 * result.k));
	result.CI_U95 = result.r + 1.96/(sqrt(sumN - 3 * result.k));

	result.Q_df = result.k - 1;
	//std::cout << result.Q << "->" << result.Q_df << std::endl;
	result.Q_p = CHIsq_Prob(result.Q, result.Q_df);

/*	result.rp = result.r/result.SE;
	result.p = Ztest(result.rp); 
	std::cout << "Z test =" << result.rp << ", p = " << result.p << std::endl;*/

	if(print) {
		if(simplePrint) {
			std::cout << result.k << " " << result.r << " " << result.CI_L95 << " " << result.CI_U95 << " " << result.Q << " " << result.Q_p << " fixed - Hedges" << std::endl;
		}
		else {
			std::cout << "______________________________________________________Hedges-Olkin Fixed Effects_____" << std::endl << std::endl;
			std::cout << "mean r = " << result.r << ", k = " << result.k << ", SE = " << result.SE << std::endl;
			std::cout << "95% CI lower = " << result.CI_L95 << ", 95% CI upper = " << result.CI_U95 << std::endl;
			std::cout << "homogeneity test: Q = " << result.Q << ", d.f. = " << result.Q_df <<", p = " << result.Q_p << ((result.Q_p <= 0.05) ? " | Significant heterogeneity!" : " | The data is homogeneous.") << std::endl;
			std::cout << "_____________________________________________________________________________________" << std::endl;
		}
	}
	return result;
}

META_REGR meta_correlation_coefficients_Rosenthal_Rubin_random(std::vector <double> rList, std::vector<int> NList, bool print, bool simplePrint)
{
	META_REGR result, resultFixed = meta_correlation_coefficients_Hedges_Olkin_fixed(rList, NList, false, false);
	std::vector<double> zList, wList, w2List, wstarList, wstarbyZlist;

	double sumN = 0.0, sumW = 0.0, sumW2 = 0.0, sumZbyN = 0.0, c = 0.0, sumWstar = 0.0, sumWstarByZ = 0.0, t2 = 0.0;

	result.k = rList.size(); result.Q = 0.0;

	for(int i = 0; i < result.k; i++) {
		zList.push_back(Z_transform(rList[i]));
		wList.push_back(NList[i] - 3);
		w2List.push_back((NList[i] - 3)*(NList[i] - 3));
	}

	for(i = 0; i < result.k; i++) {
		sumN += NList[i];
		sumW += wList[i];
		sumW2 += w2List[i];
		sumZbyN += (zList[i] * wList[i]);
	}

	c = sumW - (sumW2/sumW);
	t2 = (resultFixed.Q - (result.k - 1))/c;

	for(i = 0; i < result.k; i++) wstarList.push_back(1/(1/wList[i] + t2));
	for(i = 0; i < result.k; i++) wstarbyZlist.push_back(zList[i] * wstarList[i]);
	
	for(i = 0; i < result.k; i++) {
		sumWstar += wstarList[i];
		sumWstarByZ += wstarbyZlist[i];
	}

	result.r = sumWstarByZ/sumWstar;

	for(i = 0; i < result.k; i++) result.Q += wstarList[i] * ((zList[i] - result.r) * (zList[i] - result.r));

	result.SE = 1/sqrt(sumWstar);

	result.CI_L95 = result.r - 1.96/(sqrt(sumN - 3 * result.k));
	result.CI_U95 = result.r + 1.96/(sqrt(sumN - 3 * result.k));

	result.Q_df = result.k - 1;
	//std::cout << result.Q << " " << result.Q_df << std::endl;
	result.Q_p = CHIsq_Prob(abs(result.Q), result.Q_df);


/*	result.rp = result.r/result.SE;
	result.p = Ztest(result.rp); 
	std::cout << "Z test =" << result.rp << ", p = " << result.p << std::endl;*/

	if(print) {
		if(simplePrint) {
			std::cout << result.k << " " << result.r << " " << result.CI_L95 << " " << result.CI_U95 << " " << result.Q << " " << result.Q_p << " random" << std::endl;
		}
		else {
			std::cout << "__________________________________________________Rosenthal-Rubin Random Effects_____" << std::endl << std::endl;
			std::cout << "mean r = " << result.r << ", k = " << result.k << ", SE = " << result.SE << std::endl;
			std::cout << "95% CI lower = " << result.CI_L95 << ", 95% CI upper = " << result.CI_U95 << std::endl;
			std::cout << "homogeneity test: Q = " << result.Q << ", d.f. = " << result.Q_df <<", p = " << result.Q_p << ((result.Q_p <= 0.05) ? " | Significant heterogeneity!" : " | The data is homogeneous.") << std::endl;
			std::cout << "_____________________________________________________________________________________" << std::endl;
		}
	}


	return result;
}




#endif