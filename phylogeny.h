// PHYLOGENY.h
// Decleration of PHYLOGENY class and member functions.

#ifndef PHYLOGENY_H
#define PHYLOGENY_H

#define MAXRESAMPLES 5000
//#define PI 3.14159265

struct PHYLODIV
{
	double PD;
	double PDvar;
	int HR;
};

struct TIP            
{
	bool node;								// identidies whether node is a true node or species
	int id;									// unique identifier of node
	double BRL;								// branch-length (x-axis BL)              
	double BRL_y;							// y axis BL used for drawing
	std::string species;					// species name, when a true node is "node"
	bool ignore;                            //used in tree drawing
	double DATA_x;
	double var;
	double covar;
	int optima;
};

class Phylogeny {
public:

	Phylogeny() {tree<TIP> aPhylogeny;};								// constructor
	Phylogeny(tree<TIP> aTree) {aPhylogeny = aTree;};							// conversion constructor a phylogeny with a tree<TIP> class member
	~Phylogeny() {aPhylogeny.clear();};	//destructor

		
	Matrix metaCovariance(bool, bool, double);
	void wtgMatrix(const ColumnVector &);

	void getCoordinates();
	void showCoordinates();

	void inputTaxaXdata(std::vector < std::string >, std::vector < double>, std::vector < int >, std::vector<double>);
	void inputTaxaXdataCOVAR(std::vector < std::string >, std::vector < double>, std::vector < int >, std::vector<double>, std::vector<double>);

	void file_Open_NEWICKPhylogeny(char []);								// loads phylogeny with info on file (file data must be in parenthical NEWICK)
	std::string file_Save_NEWICKPhylogeny(char [], bool, bool);					// saves phylogeny to file in NEWICK format with (true) or without (false) branch-lengths
	std::string file_Save_NEWICKPhylogeny_STRING(std::string , bool , bool );
	void file_Save_VARCOVAR_Matrix();

	std::string format_phylogeny_to_NEWICK(bool);								// prints tree in parenthical format, take bool for true:with BRL or false:without BRL, returns string
	Matrix format_phylogeny_to_matrix(void);
	Matrix format_phylo_to_matrix(void);
	ColumnVector format_phylo_to_DATACOLUMN(void);
	ColumnVector format_phylo_to_DATACOLUMN_COVAR(void);
	ColumnVector format_phylo_to_DATACOLUMN_VAR(void);
	std::vector< std::string > format_phylo_to_DATACOLUMN_moderators(void);
	
	void print_phylogeny();														// prints phylogeny on screen as text art
	void print_allSpecies();
	void print_allGenera(std::string);
	void print_taxaData();

	bool isSpeciesOnTree(std::string);									// checks if species is on tree
	std::vector< std::string > getAllSpecies();
	std::vector< std::string > getAllGenera(std::string);

	Phylogeny transform_NEWICK_to_Phylogeny(std::string);
	void transform_phylogenyPagel();
	Phylogeny transform_phylogenyPagel(Phylogeny);
	Phylogeny transform_phylogenyConstant();
	void transform_trim_BL_Tips(double);
	Phylogeny transform_make_Ultrametric(); 
	void transform_scale_Tree(double theScale);

	//Phylogeny statistics
	int stats_numberOfSpecies();
	int stats_numberOfNodes();
	int stats_size() {return aPhylogeny.size();};								// returns number of nodes on tree
	bool stats_anyPolytomies();
	int stats_numberOfPolytomies();
	bool stats_isUltrametric();
	bool stats_isSpeciational();
	int stats_maxTreeDepth();
	int stats_minTreeDepth();
	std::vector < std::string > stats_mostDerivedSpecies();
	std::vector < std::string > stats_mostBasalSpecies();
	double stats_maxTreeBRL();
	double stats_minTreeBRL();
	double stats_averageSpeciesBRL();
	double stats_totalBLtree();
	int stats_speciesNodeDepth(std::string);
	double stats_speciesTotalBRL(std::string);
	float stats_Nbar();
	float stats_Nvar();
	float stats_C();
	float stats_B1();
	float stats_meanIprime();
	double stats_NTI_Webb(std::vector< std::string >, bool); 
	double stats_NRI_Webb(std::vector< std::string >, bool); 
	double stats_NTIspecies_Webb(std::vector< std::string >, std::string, bool);
	double stats_NRIspecies_Webb(std::vector< std::string >, std::string, bool);
	PHYLODIV stats_phylodiversity(bool, bool, std::vector<std::string>);
	int stats_howManyKtaxaCanBeResolvedIntoFullBifurcatingTrees(int);
	long unsigned stats_numberOfResolutions(int );
	long unsigned stats_numberOfResolutionsForPolytomies();

	double stats_Pagel_lambda(const ColumnVector & , bool , bool);
	double stats_Bulter_selection(const ColumnVector & , bool );
	int stats_maxOptima(void);

	//tree manipulations	
	Phylogeny treeSpeciesSubset(std::vector< std::string > );			// makes a tree based on a subset of species
	Phylogeny treeSpeciesAntiSubset(std::vector< std::string > );
	Phylogeny generate_randomBifurcatedTree(int); 
	Phylogeny randomlyResolvePolytomies(double);
	Phylogeny randomPhylogeny(int, double);
	Phylogeny addSpeciesPolitomy(std::string, int);

	Phylogeny changeSpeciesNames(std::vector < std::string >);

	//Matrix functions
	Matrix GLS_tranformation(ColumnVector, Matrix, bool, double);
	Matrix GLS_tranformation_2(Matrix , Matrix );
	LowerTriangularMatrix matrix_phylogenyGLM(const Matrix &, int);
	Matrix matrix_phylogenyTransformOU(const Matrix & , int , float);
	Matrix matrix_phylogenyTransformOU(double );

	Matrix matrix_standardize(const Matrix & );
	Matrix matrix_standardize_UNQUAL(const Matrix &);
	Matrix matrix_polarize(const Matrix & );
	bool matrix_isPositiveDefinate(const Matrix &);
	Matrix matrix_squareRoot(const Matrix &);
	bool matrix_isSymmetric(const Matrix & );
	Matrix matrix_transformLambda(double );

	void getVECTORS_BUTLER_KING();

private:
	tree<TIP> aPhylogeny;

};



//-------------------------------------------------------------------------------------
//used by createPhylogenyFromFile() to double check wether paratheses are balanced (if not something wrong with tree on file)
bool parenOK(std::string phylo)
{
	int countL = 0, countR = 0;
	for(int i = 0; i < phylo.size(); i++) {
		if(phylo[i] == '(') countL++;
		if (phylo[i] == ')') countR++;
	}
	if (!countL || !countR) return 0;
	return countL == countR;
}

void Phylogeny::file_Open_NEWICKPhylogeny(char fileName[])
{
	char phylogenyBit[15000];
	std::string phylogeny, sp, aFile = fileName; 
	std::vector<std::string> testPhylo;
	int countNode = 0;
	float BRL = 0.0, BRL2 = 0.0;
		
	ifstream phyloFile(fileName, ios::in); if(!phyloFile) errorACTION("Could not open: " + aFile, "Exiting system", true);
	while(phyloFile >> phylogenyBit) phylogeny += phylogenyBit;
	phyloFile.close();	

	if(!parenOK(phylogeny)) nrerror("phylogy format incorrect, unbalanced paraenthese...");	//write function that removes all spaces from phylogeny string!

	//puts file contents into a vector
	for(int i = 0; i < phylogeny.size(); i++) {
		if(phylogeny[i] == '(') testPhylo.push_back("(");
		if((phylogeny[i] != '(') && (phylogeny[i] != ')') && (phylogeny[i] != ',') && (phylogeny[i] != ';') && (phylogeny[i] != ':')) 
			sp += phylogeny[i];
		if((phylogeny[i] == ',') || (phylogeny[i] == ')') || (phylogeny[i] == ':')) {
				if(sp.size()) testPhylo.push_back(sp);
				sp = "";
				if(phylogeny[i] == ')') testPhylo.push_back(")");
				if(phylogeny[i] == ':') testPhylo.push_back(":");
		}
	}
	testPhylo.push_back(":"); testPhylo.push_back("0.0");

	//printStringVector(testPhylo);

	tree<TIP> treeA;
	tree<TIP>::pre_order_iterator rootI;
	std::vector<tree<TIP>::pre_order_iterator> pointerList;	
	
	TIP root, tp; 
	root.node = true; root.id = 0; root.DATA_x = 0.0; root.ignore = true; root.BRL=0.0; root.BRL_y = -5.0; root.species = "node"; root.optima = 0; //root.parent = treeA.begin();
	tp.node = false; tp.id = 1; tp.DATA_x = 0.0; tp.ignore = true; tp.BRL=0.0; tp.BRL_y = -5.0; tp.species = "tip"; tp.optima = 0; //tp.parent = treeA.begin();
	
	//construct phylogeny tree with branchlengths
	for(i = 0; i < testPhylo.size(); i++) {	
		if(testPhylo[i] == "(") {
			root.id = countNode++; // root id always zero
			//root.parent = rootI;
			if(!treeA.size()) rootI = treeA.insert(treeA.begin(), root);
			else rootI = treeA.append_child(rootI,root);
			pointerList.push_back(rootI);
		}
		else if(testPhylo[i] == ":" && testPhylo[i-1] != ")") {
			std::istringstream sin(testPhylo[i+1]); sin >> BRL;	
			tp.BRL = BRL; tp.id = countNode++; tp.species = testPhylo[i-1];
			//tp.parent = rootI;
			treeA.append_child(rootI,tp);
		} 
		else if(testPhylo[i] == ")") {
			if(testPhylo[i+1] == ":") {
				std::istringstream sin(testPhylo[i+2]); sin >> BRL2;
				(*rootI).BRL = BRL2;	
			}
			pointerList.erase(pointerList.end()-1,pointerList.end()); 
			rootI = pointerList[pointerList.size()-1];
		}
	}
	


	aPhylogeny = treeA;
}

std::string Phylogeny::file_Save_NEWICKPhylogeny(char fileName[], bool withBRL, bool append)
{
	Phylogeny aPhylo = aPhylogeny;
	std::string aNEWIKphylogeny = aPhylo.format_phylogeny_to_NEWICK(withBRL), aFile = fileName;
	const char *charNEWIKphylogeny = aNEWIKphylogeny.c_str();
	return saveFile(fileName, aNEWIKphylogeny, append);
}

std::string Phylogeny::file_Save_NEWICKPhylogeny_STRING(std::string fileName, bool withBRL, bool append)
{
	Phylogeny aPhylo = aPhylogeny;
	std::string aNEWIKphylogeny = aPhylo.format_phylogeny_to_NEWICK(withBRL);
	const char *charNEWIKphylogeny = aNEWIKphylogeny.c_str();
	return saveFile_STRING(fileName, aNEWIKphylogeny, append);
}

Matrix getRelationMatrix(const ColumnVector &sp, int numSpecies, bool printMatrices)
{
	Matrix a, b, c(numSpecies, numSpecies);
	a = sp; b = sp.t();
	
	for(int z = 0; z < numSpecies-1; z++) 
		{a = a | sp; b = b & sp.t();}
	
	for(int i = 1; i < numSpecies+1; i++)
		for(int j = 1; j < numSpecies+1; j++) 
			c(i,j) = a(i,j) && b(i,j);
	
	if(printMatrices) {
		bool f;
		for(int q = 1; q < numSpecies+1; q++) cout << sp(q) <<" ";
		cout << "sp. with shared branches" << endl << endl;
		for(i = 1; i < numSpecies+1; i++) 
			{for(int j = 1; j < numSpecies+1; j++) {f = c(i,j); cout << f << " ";} cout << endl;}
		cout << endl;	
	}

	return c;
}


int numSpeciesFromNEWICK(std::string phylo)
{
	int count=0;
	for(int i = 0; i < phylo.size(); i++) 
		if(phylo[i] == ',') count++;
	return count + 1;
}

Matrix Phylogeny::format_phylogeny_to_matrix(void)
{
	//char phylogenyBit[5000];
	Phylogeny aPhylo = aPhylogeny;
	std::string sp, phylogeny = aPhylo.format_phylogeny_to_NEWICK(true); 
	std::vector<std::string> speciesNames, testPhylo;
	int count = 0, speciesALL = 0, para = 0, rPara = 0, lPara = 0;
	float BRL = 0.0;
	bool first = true;
	
	if(!parenOK(phylogeny)) nrerror("phylogy format incorrect...");
	//write function that removes all spaces from phylogeny string!
	speciesALL = numSpeciesFromNEWICK(phylogeny);
	
	//converts NEWIK into a vector
	for(int i = 0; i < phylogeny.size(); i++) {
		if(phylogeny[i] == '(') testPhylo.push_back("(");
		if((phylogeny[i] != '(') && (phylogeny[i] != ')') && (phylogeny[i] != ',') && (phylogeny[i] != ';') && (phylogeny[i] != ':')) 
			sp += phylogeny[i];
		if((phylogeny[i] == ',') || (phylogeny[i] == ')') || (phylogeny[i] == ':')) {
				if(sp.size()) testPhylo.push_back(sp);
				sp = "";
				if(phylogeny[i] == ')') testPhylo.push_back(")");
				if(phylogeny[i] == ':') testPhylo.push_back(":");
		}
	}
	testPhylo.push_back(":"); testPhylo.push_back("0");

	Matrix varCov(speciesALL, speciesALL); varCov = 0.0; 
	ColumnVector speciesRelations(speciesALL); speciesRelations = 0.0;

	//fills diagonal with species branch lengths;
	for(int j = 0; j < testPhylo.size(); j++) {
		if(testPhylo[j] == ":" && testPhylo[j-1] != ")") {
			std::istringstream sin(testPhylo[j+1]); sin >> BRL;	
			varCov.element(count,count) = BRL;
			count++;
			speciesNames.push_back(testPhylo[j-1]);
		}
	}
	count = 0;

	//print species names & order
	//for(int lk=0; lk < speciesNames.size(); lk++) std::cout << speciesNames[lk] << std::endl;

	//fills relatedness matrix with species branch lengths;
//	std::cout << "constructing var-covar matrix ";
	for(j = testPhylo.size()-1; j >= 0; j--) {
		if(testPhylo[j] == ":" && testPhylo[j-1] == ")") {
			std::istringstream sin(testPhylo[j+1]); sin >> BRL;	
			for(int k = j+1; k >= 0; k--) {
				if(testPhylo[k] == ")") {
					para++; 
					if(first) {rPara = k; first = false;}
				}
				if(testPhylo[k] == "(") {
					para--; 
					if(!para) {lPara = k+1;	break;}
				}
			}
			//get species relation vector
			for(int i = lPara; i < rPara; i++) 
				for(int z = 0; z < speciesALL; z++) 
					if(speciesNames[z] == testPhylo[i])	speciesRelations(z+1) = 1; 

			//use species relation vector to make relationship matrices
			varCov = varCov + BRL * getRelationMatrix(speciesRelations, speciesALL, false);
			speciesRelations = 0;
			para = 0; first = true; lPara = 0; rPara = 0;
		}
	//	std::cout << ".";
	}
//	std::cout << " DONE!" << std::endl;

	return varCov;
}

//-------------------------------------------------------------------------------------

void Phylogeny::print_phylogeny()
{
	tree<TIP>::pre_order_iterator it = aPhylogeny.begin(), begin = aPhylogeny.begin(), end = aPhylogeny.end(), toParent, toParent2;
	int maxDepth = 0, m = 0, sp = 1, line = 0;
	std::vector<std::string> speciesNames;

	while(it != end) {
		if ((*it).node == false) {
			if (aPhylogeny.depth(it) > maxDepth) maxDepth = aPhylogeny.depth(it); 
			speciesNames.push_back((*it).species);
		}
		++it;
	}
	
	Matrix treeMatrix(aPhylogeny.size(), maxDepth+2); treeMatrix = 0.0;


	it = aPhylogeny.begin();

	while(it != end) {
		if((*it).node == true) {
			treeMatrix(sp,aPhylogeny.depth(it)+1) = 5;
			treeMatrix(sp,aPhylogeny.depth(it)+2) = 7;
		}
		else {
			treeMatrix(sp,aPhylogeny.depth(it)+1) = 5;
			treeMatrix(sp,aPhylogeny.depth(it)+2) = -1;
		}
		sp++;
		++it;
	}


	it = aPhylogeny.begin(); sp = 1;

	while(it != end) {
			if(sp != 1) {
				line = sp-1;
				while(treeMatrix(line, aPhylogeny.depth(it)+1) != 7 && treeMatrix(line, aPhylogeny.depth(it)+1) != 5) {
					treeMatrix(line, aPhylogeny.depth(it)+1) = 1;
					--line;
				}
				if(treeMatrix(line, aPhylogeny.depth(it)+1) == 5) treeMatrix(line, aPhylogeny.depth(it)+1) = 3;
				
			}
		sp++;
		++it;
	}

	sp=0;

	for(int x = 1; x < aPhylogeny.size()+1; x++) {
		for(int y = 1; y < maxDepth+3; y++) {
			if(treeMatrix(x,y) == -1) {
				std::cout << (char)196 << speciesNames[sp];
				sp++;
			}
			else if (treeMatrix(x,y) == 1) std::cout << (char)179; //(char)179;//|
			else if (treeMatrix(x,y) == 2) std::cout << (char)180; //7
			else if (treeMatrix(x,y) == 3) std::cout << (char)195;//L
			else if (treeMatrix(x,y) == 4) std::cout << (char)197;//T
			else if (treeMatrix(x,y) == 5) std::cout << (char)192;//_|_
			else if (treeMatrix(x,y) == 7) std::cout << (char)191;//_|_
			else if (treeMatrix(x,y) == 6) std::cout << (char)218;//_|_
			else if (treeMatrix(x,y) == 8) std::cout << (char)92;//_|_
			else std::cout << " ";
		}
		std::cout << std::endl;
	}
}

std::vector< std::string > speciesVector(const tree<TIP>& tr)
{
	tree<TIP>::pre_order_iterator it = tr.begin(), end = tr.end();
	std::vector< std::string > speciesList;
	while(it != end) {
		if((*it).node == false) speciesList.push_back((*it).species);
		++it;
    }
	return speciesList;
}


std::vector< std::string > Phylogeny::getAllSpecies()
{
	return speciesVector(aPhylogeny);
}


Phylogeny Phylogeny::changeSpeciesNames(std::vector < std::string > theSpecies)
{
	tree<TIP> temp = aPhylogeny;
	tree<TIP>::pre_order_iterator it = temp.begin(), end = temp.end();

	while(it != end) {
		if((*it).node == false) 
			(*it).species = theSpecies[stringToInt((*it).species) - 1];
		++it;
    }
	
	Phylogeny theTree = temp;
	return theTree;
}

std::string onlyGenus(std::string aSpecies, std::string theDivider) 
{
	std::string temp;
	for(int i = 0; i < aSpecies.size(); i++) {
		if(aSpecies[i] == theDivider[0]) return temp;
		temp += aSpecies[i];
	}
	return temp;
}

std::vector< std::string > Phylogeny::getAllGenera(std::string theDivider)
{
	std::vector< std::string > allSpecies = speciesVector(aPhylogeny), temp;
	for(int i = 0; i < allSpecies.size(); i++) 
		temp.push_back(onlyGenus(allSpecies[i], theDivider));
	return cleanStringVector(temp);
}

bool Phylogeny::isSpeciesOnTree(std::string aSpecies)
{
	std::vector< std::string > aSpeciesList = speciesVector(aPhylogeny);
	for(int i = 0; i < aSpeciesList.size(); i++) if(aSpecies == aSpeciesList[i]) return true;
	return false;
}

void Phylogeny::print_allSpecies()
{
	std::vector< std::string > aSpeciesList = speciesVector(aPhylogeny);
	for(int i = 0; i < aSpeciesList.size(); i++) std::cout << aSpeciesList[i] << std::endl;
}

void Phylogeny::print_allGenera(std::string theDivider)
{
	Phylogeny tempPhylo = aPhylogeny;
	std::vector< std::string > aGenusList = tempPhylo.getAllGenera(theDivider);
	for(int i = 0; i < aGenusList.size(); i++) std::cout << aGenusList[i] << std::endl;	
}



//-------------------------------------------------------------------------------------


std::string parenthicalHelper(const tree<TIP>& tr, bool BL) 
{
	tree<TIP>::iterator sib = tr.begin();
	std::string temp;
	for(tree<TIP>::sibling_iterator sib2 = tr.begin(sib); sib2 != tr.end(sib); ++sib2) {
		if((*sib2).node == false) {
			temp += (((*sib2).species) + (!BL ? "" : (":" + floatToString((*sib2).BRL))));
		}
		else temp += ("(" + parenthicalHelper(tr.subtree(sib2), BL) + ")" + (!BL ? "" : (":" + floatToString((*sib2).BRL))));
		if(tr.next_sibling(sib2) != tr.end(sib)) temp += ",";
	} 
	return temp;
}


std::string Phylogeny::format_phylogeny_to_NEWICK(bool BL)
{
	return "(" + parenthicalHelper(aPhylogeny, BL) + ");";
}

//START treeSpeciesSubset(std::vector< std::string > aSpeciesList) functions ----------------------------------------------

bool speciesOnList(std::vector< std::string > aList, std::string aSpecies) 
{
	for(int i = 0; i < aList.size(); i++) if(aList[i] == aSpecies) return true;
	return false;
}

bool anyDanglingNodes(const tree<TIP> aTree) 
{
	tree<TIP>::pre_order_iterator it = aTree.begin(), end = aTree.end();
	while(it != end) {
		if(((*it).node == true) && (aTree.number_of_children(it) == 0)) return true;
		++it;
	}
	return false;
}

tree<TIP> killAllDanglingNodes(tree<TIP> aTree) 
{
	tree<TIP>::pre_order_iterator it = aTree.begin(), end = aTree.end();
	std::vector< tree<TIP>::pre_order_iterator > aIteratorList;

	while(anyDanglingNodes(aTree) == true) {
		while(it != end) {
			if( ((*it).node == true) && (aTree.number_of_children(it) == 0) ) aIteratorList.push_back(it);
			++it;
		}
		for(int i = 0; i < aIteratorList.size(); i++) aTree.erase(aIteratorList[i]);
		it = aTree.begin(), end = aTree.end(); aIteratorList.clear();
	}
	return aTree;
}

bool anyExtendedNodes(const tree<TIP> aTree) 
{
	tree<TIP>::pre_order_iterator it = aTree.begin(), end = aTree.end();
	while(it != end) {
		if(((*it).node == true) && (aTree.number_of_children(it) == 1)) return true; //removed ((*it).id != 0) &&
		++it;
	}
	return false;
}

int numberOfExtendedNode(const tree<TIP> aTree) 
{
	int count = 0;
	tree<TIP>::pre_order_iterator it = aTree.begin(), end = aTree.end();
	while(it != end) {
		if(((*it).node == true) && (aTree.number_of_children(it) == 1)) count++;
		++it;
	}
	return count;
}

tree<TIP> killAllExtendedNodesHelper(tree<TIP> aTree) 
{
 	tree<TIP> temp;
	tree<TIP>::pre_order_iterator it = aTree.begin(), end = aTree.end(), tempIt, tempParent;
	std::vector< tree<TIP>::pre_order_iterator > aIteratorList;
	double BRLtemp = 0.0;
	while(it != end) {
		if((*it).id != 0) {
			if(((*it).node == true) && (aTree.number_of_children(it) == 1)) {
				temp = aTree.subtree(it); tempIt = temp.begin(); 
				BRLtemp = (*tempIt).BRL; 
				++tempIt;
				temp = aTree.subtree(tempIt); tempIt = temp.begin();
				(*tempIt).BRL += BRLtemp;
				tempParent = aTree.parent(it);
				aTree.append_child(tempParent, tempIt);
				aTree.erase(it);
				return aTree;
			}
		} 
		else {
			if(((*it).node == true) && (aTree.number_of_children(it) == 1)) {
				temp = aTree.subtree(it); tempIt = temp.begin(); ++tempIt;
				if(temp.number_of_children(tempIt) > 2) std::cout << "poly" <<std::endl;
				temp = aTree.subtree(tempIt);
				tempIt = temp.begin();
				(*tempIt).node = true;	(*tempIt).BRL = 0.0; (*tempIt).id = 0;
				return temp;
			}
		}
		++it;
	}
	return aTree;
}

tree<TIP> killAllExtendedNodes(tree<TIP> aTree) 
{
	tree<TIP> temp = aTree;
	//std::cout << numberOfExtendedNode(aTree) << " need cleaning...";
	while(anyExtendedNodes(temp)) {
		temp = killAllExtendedNodesHelper(temp);
	}
	return temp;
}

Phylogeny Phylogeny::treeSpeciesSubset(std::vector< std::string > aSpeciesList) 
{
	tree<TIP> aTempTree;
	aTempTree = aPhylogeny;
	std::string tempo;
	
	//printStringVector(theSpecies, 5);
//	std::cout << "---->" << aTempTree.size() << std::endl;
	Phylogeny aPhylo = aPhylogeny;
	if(aSpeciesList.size() == 1) {
		tempo = "(" + aSpeciesList[0] + ':' + floatToString(aPhylo.stats_speciesTotalBRL(aSpeciesList[0])) + ");";
		return aPhylo.transform_NEWICK_to_Phylogeny(tempo);
	}

	tree<TIP>::pre_order_iterator it = aTempTree.begin(), end = aTempTree.end();
	std::vector< tree<TIP>::pre_order_iterator > aIteratorList;

	//erase all species
	while(it != end) {
		if(((*it).node != true) && (!speciesOnList(aSpeciesList, (*it).species))) aIteratorList.push_back(it);
		++it;
	}
	for(int i = 0; i < aIteratorList.size(); i++) aTempTree.erase(aIteratorList[i]);

	//clean up tree
	aTempTree = killAllDanglingNodes(aTempTree);
	aTempTree = killAllExtendedNodes(aTempTree);

//	std::cout << "---->" << aTempTree.size() << std::endl;
	
	Phylogeny aTree = aTempTree;
	return aTree;
}




Phylogeny Phylogeny::treeSpeciesAntiSubset(std::vector< std::string > aSpeciesList) 
{
	Phylogeny aPhylo = aPhylogeny, nullPhylo;
	std::vector< std::string > allSpecies = aPhylo.getAllSpecies(), theSpecies;
	std::string tempo;	
	bool exclude = false;
	for(int i = 0; i < allSpecies.size(); i++) {
		for(int j = 0; j < aSpeciesList.size(); j++) {
			//std::cout << allSpecies[i] << " " << aSpeciesList[j] << std::endl;
			if(allSpecies[i] == aSpeciesList[j]) exclude = true;
		}
		if(exclude == false) theSpecies.push_back(allSpecies[i]);
		exclude = false;
	}

	if(theSpecies.size() == 1) {
		tempo = "(" + theSpecies[0] + ':' + floatToString(aPhylo.stats_speciesTotalBRL(theSpecies[0])) + ");";
		return aPhylo.transform_NEWICK_to_Phylogeny(tempo);
	}
	//printStringVector(theSpecies, 5);

	//printStringVector(theSpecies);

	return aPhylo.treeSpeciesSubset(theSpecies);
}

//END treeSpeciesSubset(std::vector< std::string > aSpeciesList) functions ----------------------------------------------

bool Phylogeny::matrix_isPositiveDefinate(const Matrix & aMatrix)
{
	bool posDef = true;
	Matrix temp; temp = aMatrix;
	while(!temp.IsZero()) {
		if(temp.Determinant() < 0) posDef = false;
		temp = temp.SubMatrix(1, temp.Nrows()-1, 1, temp.Ncols() -1);
	}
	return posDef;
}

/*Matrix Phylogeny::matrix_squareRoot(const Matrix & temp)
{
	SymmetricMatrix S(temp.Ncols()); 
	for(int i = 1; i < temp.Ncols()+1; i++) 
		for(int j = 1; j < temp.Ncols()+1; j++) S(i,j) = temp(i,j);
	LowerTriangularMatrix L = Cholesky(S);
	DiagonalMatrix D;
	Matrix V, U, result;
	SVD(L, D, U, V);
	return U * D * U.t();
}*/

bool Phylogeny::matrix_isSymmetric(const Matrix & aMatrix)
{
	Matrix temp = aMatrix.t();
	for(int i = 1; i <= aMatrix.Nrows(); i++)
		for(int j = 1; j <= aMatrix.Ncols(); j++)
			if(aMatrix(i,j) != temp(i,j)) {
				cout << i << " " << j << " " << aMatrix(i,j) << " " << temp(i,j) << endl;
					return false;
			}
	return true;
}


Matrix Phylogeny::matrix_phylogenyTransformOU(const Matrix & phylo, int numObs, float alpha)
{
	
	Matrix phylo2; phylo2 = phylo;
	float maxBL = phylo(1,1); 
	for(int i = 1; i < numObs+1; i++)
		for(int j = 1; j < numObs+1; j++)
			phylo2(i,j) = exp(-alpha * 2 * (maxBL - phylo(i,j)));
	return phylo2;
}


Matrix Phylogeny::matrix_phylogenyTransformOU(double selection)
{
	Phylogeny thePhylo = aPhylogeny;
	Matrix phyloMatrix = thePhylo.format_phylogeny_to_matrix(), matrixOU = phyloMatrix;
	matrixOU = 0.0;
	
	for(int i = 1; i < phyloMatrix.Ncols()+1; i++)
		for(int j = 1; j < phyloMatrix.Nrows()+1; j++) {
			//if(i == j) matrixOU(i,j) = (1.0 - exp(-2.0 * selection * phyloMatrix(i,j)))/(2.0 * selection);
			matrixOU(i,j) = (exp(-2.0 * selection * (phyloMatrix(i,i) - phyloMatrix(i,j))) * (1.0 - exp(-2.0 * selection * phyloMatrix(i,j))))/(2.0*selection);
		}
	return matrixOU;
}

/*Matrix Phylogeny::matrix_designMatrixOUOptima(double selection, int numOptima, Matrix optima)
{
	Phylogeny thePhylo = aPhylogeny;
	Matrix phyloMatrix = thePhylo.format_phylogeny_to_matrix(), designMatrix (1, phyloMatrix.Nrows());
	ColumnVector optimaVector(phyloMatrix.Nrows()); optimaVector = 0.0;
	
	for(int i = 1; i < optimaVector.Nrows()+1; i++)	designMatrix(1, i) = exp(-1.0 * selection * phyloMatrix(i,i));
}
*/

LowerTriangularMatrix Phylogeny::matrix_phylogenyGLM(const Matrix &V, int N) 
{
	Matrix varCov = matrix_phylogenyTransformOU(V, N, 0.1);
	varCov = V;
	varCov = varCov.i();
	SymmetricMatrix S(N); for(int i = 1; i < N+1; i++) for(int j = 1; j < N+1; j++) S(i,j) = varCov(i,j);
	LowerTriangularMatrix L = Cholesky(S);
	return L;
}


Matrix Phylogeny::matrix_standardize_UNQUAL(const Matrix & phylo)
{
	double temp = 0;
	Matrix phyloTemp; phyloTemp = phylo;
	std::vector<double> totalBRL;
	
	for(int k = 1; k < phylo.Ncols()+1; k++) totalBRL.push_back(phylo(k,k));

	for(int i = 1; i < phylo.Ncols()+1; i++) 
		for(int j = 1; j < phylo.Nrows()+1; j++) {
			temp = phyloTemp(i,j);
			phyloTemp(i,j) = temp/((totalBRL[i-1]+totalBRL[j-1])/2);
		}
		
	return phyloTemp;

}


Matrix Phylogeny::matrix_standardize(const Matrix & phylo)
{
	double temp = 0;
	Matrix phyloTemp;
	for(int i = 1; i < phylo.Ncols()+1; i++) {
		for(int j = 1; j < phylo.Nrows()+1; j++) {
			if(phylo(i,j) > temp) temp = phylo(i,j);
		}
	} 
	//cout << "temp===" << 1/temp << endl;
	return (phylo/temp);
}

Matrix Phylogeny::matrix_polarize(const Matrix & phylo)
{
	Matrix phyloTemp; phyloTemp = phylo;
	for(int i = 1; i < phylo.Ncols()+1; i++) {
		for(int j = 1; j < phylo.Nrows()+1; j++) {
			phyloTemp(i,j) = 1 - phylo(i,j);
		}
	} 
	return phyloTemp;
}

std::vector<tree<TIP>::pre_order_iterator> speciesPointerVector(const tree<TIP>& tr)
{
	tree<TIP>::pre_order_iterator it = tr.begin(), end = tr.end();
	std::vector<tree<TIP>::pre_order_iterator> speciesListPointers;
	while(it != end) {
		if((*it).node == false) speciesListPointers.push_back(it);
		++it;
    }
	return speciesListPointers;
}



void Phylogeny::transform_trim_BL_Tips(double valToTrim)
{
	std::vector<tree<TIP>::pre_order_iterator> speciesListPointers = speciesPointerVector(aPhylogeny);
	for(int i = 0; i < speciesListPointers.size(); i++) 
		(*speciesListPointers[i]).BRL = (*speciesListPointers[i]).BRL - valToTrim;
}


//-------------------------------------------------------------------------------------

int Phylogeny::stats_numberOfSpecies()
{
	tree<TIP> tr = aPhylogeny;
	tree<TIP>::pre_order_iterator it = tr.begin(), end = tr.end();
	int i = 0;
	while(it != end) {
		if((*it).node == false) ++i;
		++it;
    }
    return i;
}

int Phylogeny::stats_numberOfNodes()
{
	tree<TIP> tr = aPhylogeny;
	tree<TIP>::pre_order_iterator it = tr.begin(), end = tr.end();
	int i = 0;
	while(it != end) {
		if((*it).node == true) ++i;
		++it;
    }
    return i;
}

int Phylogeny::stats_maxTreeDepth()
{
	tree<TIP> tr = aPhylogeny;
	tree<TIP>::pre_order_iterator it = tr.begin(), end = tr.end();
	int maxDepth = 0;
	while(it != end) {
		if ((*it).node == false) if (tr.depth(it) > maxDepth) maxDepth = tr.depth(it); 
		++it;
	}
	return maxDepth;
}

int Phylogeny::stats_minTreeDepth()
{
	Phylogeny aPhylo = aPhylogeny;
	tree<TIP> tr = aPhylogeny;
	tree<TIP>::pre_order_iterator it = tr.begin(), end = tr.end();
	int minDepth = aPhylo.stats_maxTreeDepth();
	while(it != end) {
		if ((*it).node == false) if (tr.depth(it) < minDepth) minDepth = tr.depth(it); 
		++it;
	}
	return minDepth;
}

bool Phylogeny::stats_isSpeciational()
{
	Phylogeny aPhylo = aPhylogeny, tempPhylo = aPhylogeny;
	tempPhylo = tempPhylo.transform_phylogenyConstant();
	Matrix phyloA, phyloB;
	phyloA = aPhylo.format_phylogeny_to_matrix();
	phyloB = tempPhylo.format_phylogeny_to_matrix();

	//cout << phyloA << phyloB << endl;

	if(phyloA == phyloB) return true;
	else return false;
}


std::vector < std::string > Phylogeny::stats_mostDerivedSpecies()
{
	Phylogeny aPhylo = aPhylogeny;
	std::vector < std::string > speciesList = aPhylo.getAllSpecies(), theSpeciesVector;
	std::string theSpecies;
	int depth = 0;

	for(int i = 0; i < speciesList.size(); i++) {
		if(depth < aPhylo.stats_speciesNodeDepth(speciesList[i])) {
			theSpecies = speciesList[i];
			depth = aPhylo.stats_speciesNodeDepth(speciesList[i]);
		}
	}

	for(i = 0; i < speciesList.size(); i++) {
		if(aPhylo.stats_speciesNodeDepth(theSpecies) == aPhylo.stats_speciesNodeDepth(speciesList[i])) {
			theSpeciesVector.push_back(speciesList[i]);
		}
	}
	return theSpeciesVector;
}


std::vector < std::string > Phylogeny::stats_mostBasalSpecies()
{
	Phylogeny aPhylo = aPhylogeny;
	std::vector< std::string > speciesList = aPhylo.getAllSpecies(), theSpeciesVector;
	std::string theSpecies;
	int depth = aPhylo.stats_maxTreeDepth();

	for(int i = 0; i < speciesList.size(); i++) {
		if(aPhylo.stats_speciesNodeDepth(speciesList[i]) < depth) {
			theSpecies = speciesList[i];
			depth = aPhylo.stats_speciesNodeDepth(speciesList[i]);
		}
	}

	//collect all sp with same max depth
	for(i = 0; i < speciesList.size(); i++) {
		if(aPhylo.stats_speciesNodeDepth(theSpecies) == aPhylo.stats_speciesNodeDepth(speciesList[i])) {
			theSpeciesVector.push_back(speciesList[i]);
		}
	}
	return theSpeciesVector;
}

float maxTreeDepthofTree(tree<TIP> tr)
{
	tree<TIP>::pre_order_iterator it = tr.begin(), end = tr.end();
	int maxDepth = 0;
	while(it != end) {
		if ((*it).node == false) if (tr.depth(it) > maxDepth) maxDepth = tr.depth(it); 
		++it;
	}
	return maxDepth;
}



int Phylogeny::stats_speciesNodeDepth(std::string aSpecies) 
{
	tree<TIP> tr = aPhylogeny;
	tree<TIP>::pre_order_iterator it = tr.begin(), end = tr.end();
	int maxDepth = 0; 
	bool found = false;
	while((it != end) && (found != true)) {
		if ((*it).species == aSpecies) {
			maxDepth = tr.depth(it);
			found = true;
		}
		++it;
	}
	return maxDepth;	
}

double Phylogeny::stats_speciesTotalBRL(std::string aSpecies) 
{
	tree<TIP> tr = aPhylogeny;
	tree<TIP>::pre_order_iterator it = tr.begin(), end = tr.end(), theSpecies;
	double maxDepth = 0.0; 
	bool found = false;
	while((it != end) && (found != true)) {
		if ((*it).species == aSpecies) {
			found = true;
			theSpecies = it;
		}
		++it;
	}

	while((*theSpecies).id != 0) {
		maxDepth += (*theSpecies).BRL;
		theSpecies = aPhylogeny.parent(theSpecies);
	}
	return maxDepth;	
}


double Phylogeny::stats_maxTreeBRL()
{
	std::vector<tree<TIP>::pre_order_iterator> allTips = speciesPointerVector(aPhylogeny);
	tree<TIP>::pre_order_iterator temp;
	double maxBRL = 0.0, tempBRL = 0.0;
	for(int i = 0; i < allTips.size(); i++) {
		temp = allTips[i];
		while((*temp).id != 0) {
			tempBRL += (*temp).BRL;
			temp = aPhylogeny.parent(temp);
		}
		if(tempBRL > maxBRL) maxBRL = tempBRL;
		tempBRL = 0.0;
	}
	return maxBRL;
}

double Phylogeny::stats_minTreeBRL()
{
	Phylogeny thePhylo = aPhylogeny;
	std::vector< tree<TIP>::pre_order_iterator > allTips = speciesPointerVector(aPhylogeny);
	tree<TIP>::pre_order_iterator temp;
	double minBRL = thePhylo.stats_maxTreeBRL(), tempBRL = 0.0;
	for(int i = 0; i < allTips.size(); i++) {
		temp = allTips[i];
		while((*temp).id != 0) {
			tempBRL += (*temp).BRL;
			temp = aPhylogeny.parent(temp);
		}
		if(tempBRL < minBRL) minBRL = tempBRL;
		tempBRL = 0.0;
	}
	return minBRL;
}

double Phylogeny::stats_totalBLtree()
{
	Phylogeny thePhylo = aPhylogeny;
	std::vector< std::string > allSpecies = thePhylo.getAllSpecies();
	double BRL = 0.0;

	for(int i = 0; i < allSpecies.size(); i++) BRL += thePhylo.stats_speciesTotalBRL(allSpecies[i]);

	return BRL;
}

double Phylogeny::stats_averageSpeciesBRL()
{
	Phylogeny thePhylo = aPhylogeny;
	//std::vector< std::string > allSpecies = thePhylo.getAllSpecies();
	return (thePhylo.stats_totalBLtree() / thePhylo.stats_numberOfSpecies());
}

tree<TIP>::pre_order_iterator findSpeciesPointer(std::vector<tree<TIP>::pre_order_iterator> allSpeciesPointer, std::string sp) 
{
	int i = 0;
	while((*allSpeciesPointer[i]).species != sp) {
		i++;
		//std::cout << i << " " << sp << std::endl;
		if(i == allSpeciesPointer.size()) {
			std::cout << sp << " not on tree!!!" << std::endl;
			nrerror("check file & species list...");
		}
	}
	return allSpeciesPointer[i];
}

std::vector<tree<TIP>::pre_order_iterator> getSpeciesPointerFromString(std::vector<tree<TIP>::pre_order_iterator> pointerVector, std::vector<std::string> speciesVector) 
{
	tree<TIP>::pre_order_iterator it;
	std::vector<tree<TIP>::pre_order_iterator> speciesListPointers;
	for(int i = 0; i < speciesVector.size(); i++) 
		speciesListPointers.push_back(findSpeciesPointer(pointerVector, speciesVector[i]));
	return speciesListPointers;
}

tree<TIP>::pre_order_iterator shareNode(const tree<TIP>& tr, std::string speciesA, std::string speciesB)
{
	std::vector<tree<TIP>::pre_order_iterator> speciesListPointers;
	tree<TIP>::pre_order_iterator pointA, pointB, start = tr.begin(), sharedNode = tr.begin(), temp;
	speciesListPointers = speciesPointerVector(tr);
	bool found = false;

	for(int i = 0; i < speciesListPointers.size(); i++) {
		if(speciesA == (*speciesListPointers[i]).species) pointA = speciesListPointers[i];
		if(speciesB == (*speciesListPointers[i]).species) pointB = speciesListPointers[i];
	}

	temp = pointB;
	while(found != true) {
		while(found != true) {
			if((*pointA).id == (*pointB).id) {sharedNode = pointB; found = true;}
			if((*pointB).id == 0) break;
			else pointB = tr.parent(pointB);
		}
		if((*pointA).id == 0) break;
		else pointA = tr.parent(pointA);
		pointB = temp;
    }
	return sharedNode;
}


double speciesToSpeciesDistance(const tree<TIP>& tr, std::string speciesA, std::string speciesB, bool sharedNodeDistance)
{
	tree<TIP>::pre_order_iterator end = tr.begin(), pointA, pointB;
	std::vector<tree<TIP>::pre_order_iterator> speciesListPointers;
	speciesListPointers = speciesPointerVector(tr);
	double A = 0, B = 0;
	if(sharedNodeDistance == true) {
		for(int i = 0; i < speciesListPointers.size(); i++) {
			if(speciesA == (*speciesListPointers[i]).species) pointA = speciesListPointers[i];
			if(speciesB == (*speciesListPointers[i]).species) pointB = speciesListPointers[i];
		}
		end = shareNode(tr, (*pointA).species, (*pointB).species);
	}		
	
	if(pointA == pointB) {
		if ((*pointA).id == (*pointB).id) A = ((*pointA).BRL)/2;
		else A = (*pointA).BRL;
	}
	else {
		while (pointA != end) {
			A += (*pointA).BRL;
			pointA = tr.parent(pointA);
		}
	}

	return A;
} 



PHYLODIV Phylogeny::stats_phylodiversity(bool withBRL, bool useAllSpecies, std::vector<std::string> spVector)
{
	Phylogeny tempTree = aPhylogeny;
	tree<TIP>::pre_order_iterator it;
	std::vector<tree<TIP>::pre_order_iterator> speciesListPointers, speciesListPointers2;
	double PD = 0.0, N = 0, sumSpecies = 0.0, wij = 0.0, w = 0.0, varPD = 0.0;

	speciesListPointers2 = speciesPointerVector(aPhylogeny);
	if(useAllSpecies == true) speciesListPointers = speciesListPointers2;
	else speciesListPointers = getSpeciesPointerFromString(speciesListPointers2, spVector);
	N = speciesListPointers.size();
	PHYLODIV theResult;
	//cout << "Max tree depth = " << tempTree.stats_maxTreeBRL() << endl;

	if(N == 0 || N == 1) {
		cout << N << " " << 0 << " " << 0 << endl; 
		theResult.PD = 0.0; theResult.PDvar = 0.0;
		return theResult;
	}

	for(int i = 0; i < speciesListPointers.size(); i++) {
		for(int j = 0; j < speciesListPointers.size(); j++) {
			if (i != j) {
				wij = speciesToSpeciesDistance(aPhylogeny,  (*speciesListPointers[i]).species, (*speciesListPointers[j]).species, true);
				sumSpecies += wij;
				w += (wij * wij);
			}
		}
	}

	PD = (1/(N*(N-1))) * sumSpecies;
	varPD = (w/(N*(N-1))) - (PD * PD);

	theResult.PD = PD;
	theResult.PDvar = varPD;
	theResult.HR = N;
	cout << N << " " << PD << " " << varPD << endl;
    return theResult;
}


double alpha(const ColumnVector &y, const Matrix &V)
{
	ColumnVector X(y.Nrows()); X = 1.0;
	Matrix alpha = (X.t() * V.i() * X).i() * (X.t() * V.i() * y);
	return alpha(1,1);
}

double alphaError(const ColumnVector &y, const Matrix &V, double alpha)
{
	ColumnVector X(y.Nrows()); X = 1.0;
	//outcomes varies a bit with Pagel's bayestraits depending on wether y.nrows is substracted by 1.0;
	//in with -1.0 is absent ML estimate = same, lambda = same, alpha = same, but var-alpha differs
	//chose to include -1.0 becuase lambda= same, and alpha & alpha-var same as bayestraits
	Matrix error = (1.0 / (y.Nrows() - 1.0)) * ((y - alpha * X).t() * V.i() * (y - alpha * X)); 
	//Matrix error = alpha * (X.t()*V.i()*X) * alpha;
	//Matrix error = (X.t()*V.i()*X).i();
	return error(1,1);
}

Matrix transformMatrixWithLambda(const Matrix &covariance, double lambda, bool asButler)
{
	Matrix temp = covariance;

	if(asButler == true) lambda = lambda * lambda;

	for(int i = 1; i < temp.Ncols() + 1; i++)
		for(int j = 1; j < temp.Ncols() + 1; j++) 
			if (i != j) temp(i,j) = lambda * temp(i,j);
	
	return temp;
}

Matrix Phylogeny::matrix_transformLambda(double lambda)
{
	Phylogeny tempPhylo = aPhylogeny;
	
	Matrix temp = tempPhylo.format_phylo_to_matrix();
	//temp = tempPhylo.matrix_standardize_UNQUAL(temp);

	for(int i = 1; i < temp.Ncols() + 1; i++)
		for(int j = 1; j < temp.Ncols() + 1; j++) //temp(i,j) = lambda * temp(i,j);
			if (i != j) temp(i,j) = lambda * temp(i,j);
	
	return temp;
}

double MaxLik(const ColumnVector &y, const Matrix &V, double alpha, double error)
{
	ColumnVector X(y.Nrows()); X = 1.0;
	Matrix lik2;
	
	//using Butler & King 2004 Am. Nat. 164's definition
	lik2 = (y - alpha * X).t() * V.i() * (y - alpha * X);
	double lik = (1.0 / sqrt(pow(2.0 * PI * error, y.Nrows()) * V.Determinant())) * exp(lik2(1,1)/(-2.0 * error));
	
	//result same as above but using Freckleton et al (2002) AM. Nat. 160's definition
	//lik2 = (-1.0 / (2.0 * error)) * ((y - alpha * X).t() * V.i() * (y - alpha * X));
	//double lik = (1.0 / pow(2.0 * PI * error, y.Nrows()/2.0)) * (1.0/sqrt(V.Determinant())) * exp(lik2(1,1));

	return log(lik);		
}





double likLambda(const ColumnVector &x, const Matrix &V, double lambda, bool asButler) 
{
	double a = 0.0, e = 0.0, ML = 0.0;
	Matrix temp = transformMatrixWithLambda(V, lambda, asButler);
	a = alpha(x, temp); 	
	e = alphaError(x, temp, a);	
	return ML = MaxLik(x, temp, a, e); 
}


Matrix matrixCovariance(const Matrix &correlations, const tree<TIP>& aTree)
{
//	Matrix temp = correlations;
	
//	temp(5,4) = 0.0; temp(4,5) = 0.0;
//	temp(5,3) = 0.0; temp(3,5) = 0.0;
//	temp(4,3) = 0.0; temp(3,4) = 0.0;

//	cout << correlations << endl << temp << endl;
	std::vector<tree<TIP>::pre_order_iterator> allTips = speciesPointerVector(aTree);
	
	Matrix covariance(correlations.Ncols(), correlations.Ncols()); 
	Matrix SD(correlations.Ncols(), correlations.Ncols()); SD = 0.0;

	for(int i = 1; i < covariance.Ncols()+1; i++)
		for(int j = 1; j < covariance.Ncols()+1; j++) {
			if (i == j) SD(i,j) = sqrt((*allTips[i-1]).var);
			//if(i == j) covariance(i,j) = (*allTips[i-1]).var;// * correlations(i,j);
			//else covariance(i,j) = sqrt((*allTips[i-1]).var) * correlations(i,j)  * sqrt((*allTips[j-1]).var);
		}

		//cout << setprecision(3) << covariance << endl;
		//cout << correlations << endl << endl;
		//cout << SD * correlations * SD << endl;

	return SD * correlations * SD;
}

Matrix Phylogeny::metaCovariance(bool withCorrelations, bool lambda, double theLambda)
{
	Phylogeny tempTree = aPhylogeny;
	std::vector<tree<TIP>::pre_order_iterator> allTips = speciesPointerVector(aPhylogeny);
	Matrix theCov(allTips.size(), allTips.size()); theCov = 0.0;

	if(withCorrelations == false) {
		for(int i = 1; i < allTips.size()+1; i++) theCov(i,i) = (*allTips[i-1]).var;
		//cout << theCov.i() << endl;
		return theCov;
	}
	else {
		Matrix A = tempTree.format_phylogeny_to_matrix();
		A = tempTree.matrix_standardize_UNQUAL(A);
		if(lambda == true) A = tempTree.matrix_transformLambda(theLambda);
		return matrixCovariance(A, aPhylogeny);
	}
}

double lambda(const ColumnVector & x, const Matrix & V, double MAX, double MIN, bool showML, bool asButler, const tree<TIP>& aTree) 
{
	double MLMzero = 0.0, MLMone = 0.0, a = 0.0, am= 0.0, em = 0.0, MLm = 0.0, b = 0.0, e = 0.0, ML = 0.0, MLMAX = 0.0, start = 0.0,  tempLambda = MIN;
	Matrix lambdaV;
		
	MLMzero = likLambda(x, V, 0.0, asButler);
	MLMone = likLambda(x, V, 1.0, asButler);
	std::vector <double> MLvector;
	double progress = MIN;

	cout << "Estmating ML estimate of lambda:" << endl;
	for(double lambda = MIN; lambda <= MAX + MIN; lambda += MIN) {
		lambdaV = transformMatrixWithLambda(V, lambda, asButler);
		//cout << lambdaV << endl;
		 //ambdaV = matrixCovariance(lambdaV, aTree);
		a = alpha(x, lambdaV); 	
		e = alphaError(x, lambdaV, a);	
		ML = MaxLik(x, lambdaV, a, e);  

		if(lambda == MIN) {MLMAX = ML; tempLambda = lambda; start = ML;}
		if(ML >= MLMAX) {MLMAX = ML; tempLambda = lambda;}

		//cout << lambda << " " << progress << endl;
		if(lambda > progress) {
			printProgress(lambda, MAX);
			progress += (250.0 * 0.0001);
		}

		if(showML == true) cout << lambda << " | " << ML << " (" << ML - MLMzero << ", " << ML - MLMone << ") ->" << MLMAX << endl;
	}

	cout << setprecision(4) << "lambda = " << tempLambda << ", " << (char)224 << " = " << alpha(x, transformMatrixWithLambda(V, tempLambda, asButler)) << ", " << (char)229 << (char)253 << "(" << (char)224 << ") = " << alphaError(x, transformMatrixWithLambda(V, tempLambda, asButler), alpha(x, transformMatrixWithLambda(V, tempLambda, asButler))) << ", ln(L) = " << MLMAX << endl;
	cout << "  " << char(192) << char(196) << char(196) << char(194) << char(196) << char(170) << " Equal to 1? Chi-square = " << -2.0*(MLMzero - MLMAX) << ", d.f. = 1, p = " << CHIsq_Prob(-2.0*(MLMzero - MLMAX), 1) << ", ln(L) = " << MLMzero << endl;
	cout << "     " << char(192) << char(196) << char(170) << " Equal to 0? Chi-square = " << -2.0*(MLMone - MLMAX) << ", d.f. = 1, p = " << CHIsq_Prob(-2.0*(MLMone - MLMAX), 1) << ", ln(L) = " << MLMone << endl;
	
	//cout << endl << "    a) lambda = one? MLE = " << MLMzero << " | Chi-square = " << -2.0*(MLMzero - MLMAX) << ", p = " << CHIsq_Prob(-2.0*(MLMzero - MLMAX), 1) << endl;
	//cout << "    lambda differs from one? MLE = " << MLMone << " | Chi-square = " << -2.0*(MLMone - MLMAX) << ", p = " << CHIsq_Prob(-2.0*(MLMone - MLMAX), 1) << endl;

	return tempLambda;
}


double Phylogeny::stats_Pagel_lambda(const ColumnVector & y, bool showDetails, bool asButler)
{
	Phylogeny tempTree = aPhylogeny;
	Matrix A = tempTree.format_phylogeny_to_matrix(), temp;
	//IdentityMatrix I(A.Ncols());
	//temp = I;
	A = tempTree.matrix_standardize_UNQUAL(A);
	//temp = matrixCovariance(I, aPhylogeny);
	//cout << alpha(y, temp) << " " << alphaError(y, temp, alpha(y, temp)) << endl;

	return lambda(y, A, 2, 0.0001, showDetails, asButler, aPhylogeny);
}

Matrix phylogenyTransformOU(const Matrix &V, double theSelection)
{
	Matrix phyloMatrix = V, matrixOU = phyloMatrix;	matrixOU = 0.0;
	for(int i = 1; i < phyloMatrix.Ncols()+1; i++)
		for(int j = 1; j < phyloMatrix.Nrows()+1; j++) 
			matrixOU(i,j) = (exp(-2.0 * theSelection * (phyloMatrix(i,i) - phyloMatrix(i,j))) * (1.0 - exp(-2.0 * theSelection * phyloMatrix(i,j))))/(2.0 * theSelection);
	return matrixOU;
}

double nodalDistanceToRoot(const tree<TIP>& aTree, tree<TIP>::pre_order_iterator aNode) 
{
	double maxDepth = 0.0;
	while((*aNode).id != 0) {
		maxDepth += (*aNode).BRL;
		aNode = aTree.parent(aNode);
	}
	return maxDepth;
}


Matrix weightMatrix(const tree<TIP>& aTree, bool withoutSelection, double theSelection)
{
	Phylogeny thePhylo = aTree;
	std::vector<tree<TIP>::pre_order_iterator> allTips = speciesPointerVector(aTree);
	Matrix W(thePhylo.stats_numberOfSpecies(),2); W = 0.0;
	double s = 0.0;

	for(int i = 1; i < W.Nrows()+1; i++) {
		//cout << thePhylo.stats_speciesTotalBRL((*allTips[i-1]).species) << endl;

		//goodW(i,2) = exp(-1.0 * theSelection * thePhylo.stats_speciesTotalBRL((*allTips[i-1]).species));
		
		//W(i,3) = exp(-1.0*theSelection*(thePhylo.stats_speciesTotalBRL((*allTips[i-1]).species)));

		if ((*allTips[i-1]).optima == 0) {
			//goodW(i,3) = 1.0 - exp(-1.0 * theSelection * thePhylo.stats_speciesTotalBRL((*allTips[i-1]).species));
			//W(i,3) = 0.000001;

			//W(i,2) = exp(-1.0 * theSelection * thePhylo.stats_speciesTotalBRL((*allTips[i-1]).species));

			//W(i,2) = 0;
			//good2W(i,1) = 1.0;
			//W(i,2) = W(i,1) * exp(-1.0 * theSelection * thePhylo.stats_speciesTotalBRL((*allTips[i-1]).species));

			//s = thePhylo.stats_speciesTotalBRL((*allTips[i-1]).species);
			//W(i,2) = exp(-1.0 * theSelection * 0.0) - exp(-1.0 * theSelection * s);
			//W(i,1) = 1 - exp(-1.0 * theSelection * 0.0);



			// USE when optima zero is unkown ancestral, not another optima
			if(withoutSelection == 1) {W(i,1) = 0.0; W(i,2) = 1.0;}
			else {
				W(i,1) = exp(-1.0 * theSelection * thePhylo.stats_speciesTotalBRL((*allTips[i-1]).species));
				W(i,2) = 1 - W(i,1);
			}

			/* 
			// when optima 1 is derived from 0, w/ no ancestor
			W(i,1) = 0;
			W(i,2) = 1;
			*/
			  
			

		}
		else if ((*allTips[i-1]).optima == 1) {
			
			if(withoutSelection == 1) {W(i,1) = 1.0; W(i,2) = 0.0;}
			else {
				s = nodalDistanceToRoot(aTree, aTree.parent(allTips[i-1]));
				W(i,2) = exp(-1.0 * theSelection * s) - exp(-1.0 * theSelection * thePhylo.stats_speciesTotalBRL((*allTips[i-1]).species));
				W(i,1) = 1 - exp(-1.0 * theSelection * s);
			}
			
			
			 /* 
			// when optima 1 is derived from 0, w/ no ancestor
			s = nodalDistanceToRoot(aTree, aTree.parent(allTips[i-1]));
			W(i,2) = exp(-1.0 * theSelection * (thePhylo.stats_speciesTotalBRL((*allTips[i-1]).species) - s));
			W(i,1) = 1.0 - exp(-1.0 * theSelection * (thePhylo.stats_speciesTotalBRL((*allTips[i-1]).species) - s));
			*/
			
			
			//s = nodalDistanceToRoot(aTree, aTree.parent(allTips[i-1]));
			//goodW(i,4) = (exp(-1.0 * theSelection * (thePhylo.stats_speciesTotalBRL((*allTips[i-1]).species) - s))) * (1.0 - exp(-1.0 * theSelection * s));
			//goodW(i,3) = 1.0 - exp(-1.0 * theSelection * (thePhylo.stats_speciesTotalBRL((*allTips[i-1]).species) - s));
			//good2W(i,2) = (1.0 - exp(-1.0 * theSelection * (thePhylo.stats_speciesTotalBRL((*allTips[i-1]).species) - s)));
			//W(i,2) = exp(-1.0 * theSelection * s);
			//good2W(i,1) = (exp(-1.0 * theSelection * (thePhylo.stats_speciesTotalBRL((*allTips[i-1]).species) - s)));

			
			
			//W(i,3) = 1.0;
			//W(i,2) = -1.0;
		}
		else {
			//W(i,2) = 1.0;
			cout << "error" << endl;
		}
		//W(i,2) = exp(-1.0*theSelection*s) - exp(-1.0*theSelection*thePhylo.stats_speciesTotalBRL((*allTips[i-1]).species));
		//W(i,1) = 1 - exp(-1.0*theSelection*s);
	}

	return W;
}






/*
Matrix weightMatrix(const tree<TIP>& aTree, double theSelection)
{
	Phylogeny thePhylo = aTree;
	std::vector<tree<TIP>::pre_order_iterator> allTips = speciesPointerVector(aTree);
	Matrix W(thePhylo.stats_numberOfSpecies(),4); W = 0.0;
	double s = 0.0;


//	std::vector< int> theOptima;
//	for(int j = 0; j < allTips.size(); j++) 
//		theOptima.push_back((*allTips[j]).optima);
//	theOptima = cleanIntVector(theOptima);
//	printIntVector(theOptima);

	for(int i = 1; i < W.Nrows()+1; i++) {
		W(i,1) = exp(-1.0 * theSelection * thePhylo.stats_speciesTotalBRL((*allTips[i-1]).species));
		if ((*allTips[i-1]).optima == 0) {
			//W(i,2) = 1.0 - exp(-1.0 * theSelection * thePhylo.stats_speciesTotalBRL((*allTips[i-1]).species));
			W(i,2) = W(i,1) * exp(-1.0 * theSelection * thePhylo.stats_speciesTotalBRL((*allTips[i-1]).species));
		}
		else if ((*allTips[i-1]).optima == 1) {
			s = nodalDistanceToRoot(aTree, aTree.parent(allTips[i-1]));
			//cout << s << endl;
			W(i,2) = (exp(-1.0 * theSelection * (thePhylo.stats_speciesTotalBRL((*allTips[i-1]).species) - s))) * (1.0 - exp(-1.0 * theSelection * s));
			W(i,3) = 1.0 - exp(-1.0 * theSelection * (thePhylo.stats_speciesTotalBRL((*allTips[i-1]).species) - s));
		}
		else if ((*allTips[i-1]).optima == 2) {
			s = nodalDistanceToRoot(aTree, aTree.parent(allTips[i-1]));
		//	//cout << s << endl;
			W(i,2) = (exp(-1.0 * theSelection * (thePhylo.stats_speciesTotalBRL((*allTips[i-1]).species) - s))) * (1.0 - exp(-1.0 * theSelection * s));
			W(i,4) = 1.0 - exp(-1.0 * theSelection * (thePhylo.stats_speciesTotalBRL((*allTips[i-1]).species) - s));
		}
		else {
			cout << "error" << endl;
		}
		//cout << i << endl;
	}

//	cout << W << endl;
	return W;
}
*/




void Phylogeny::wtgMatrix(const ColumnVector &zzzz)
{
	
	Phylogeny thePhylo;
	thePhylo.file_Open_NEWICKPhylogeny("Butler_King_2004_fig5.phy");
	thePhylo.print_phylogeny();

	double yVal;
	ColumnVector y(thePhylo.stats_numberOfSpecies());
	int j = 1;
	ifstream theFile("Butler_King_2004_fig5.dat");
	while(theFile >> yVal) {
		y(j) = log(yVal);
		j++;
	}
	theFile.close();
	
		
	Matrix V = thePhylo.format_phylogeny_to_matrix(), W;
//	IdentityMatrix I(V.Ncols()); 
//	V = I;
	
	ColumnVector te(V.Nrows()), te2(V.Nrows()); 

	
	//for(double a = 0.01; a < 15; a+= 0.01) {
	//W = weightMatrix(aPhylogeny, true, 1);
//	W = (1.0/W(1,1)) * W;

	//ColumnVector W2(13); W2 = 0;
	//W = W2;
	//W(1,1) = 1.0; W(2,1) = 1.0; W(6,1) = 1.0; W(10,1) = 1.0; W(11,1) = 1.0; W(12,1) = 1.0;
	
	//V = phylogenyTransformOU(V, 1);

	//double lambda = tempTree.stats_Pagel_lambda(y, false, false);

	//V = tempTree.matrix_transformLambda(lambda);

	//V = tempTree.metaCovariance(true,true,1);

	//V = tempTree.matrix_standardize_UNQUAL(V);
//	V = tempTree.metaCovariance(true,false,0.01);
	Matrix mean, var, qr, qh;
	
	cout << W << V << endl;
		
	/*	mean = (W.t() * V.i() * W).i() * (W.t() * V.i() * y);
		var = ((W.t() * V.i() * W).i()).t();
		qr = mean.t() * W.t() * V.i() * W * mean;
		qh = y.t() * V.i() * y - qr;

		cout << mean << var << qr << qh << endl;

		cout << "mean 1 = " << mean(1,1) << ", (" << (mean(1,1) - 1.96 * sqrt(var(1,1))) << ", " << (mean(1,1) + 1.96 * sqrt(var(1,1))) << ")" << endl;
		cout << "mean 2 = " << mean(2,1) << ", (" << (mean(2,1) - 1.96 * sqrt(var(2,2))) << ", " << (mean(2,1) + 1.96 * sqrt(var(2,2))) << ")" << endl;
		//cout << "variance = " << var(1,1) << endl;
		cout << "Is non zero? Q = " << qr(1,1) << ", d.f. = " << y.Nrows() - 1 << ", p = " << CHIsq_Prob(qr(1,1), y.Nrows() - 1) << endl;
		cout << "Q = " << qh(1,1) << ", d.f. = " << 1 << ", p = " << CHIsq_Prob(qh(1,1), 1) << endl;

*/

	
	//Matrix lik2 = (1.0/error(1,1)) * (y - W * mean).t() * V.i() * (y - W * mean);
	//double lik = y.Nrows()*log(2*PI*error(1,1)) + log(V.Determinant()) + lik2(1,1);
	//cout << "AIC = " << 2+2*log(lik) << endl;
	
	double tempU = 500, min = 500;
	Matrix tempV = V;
/*	for(double i = 0.001; i < 5.0; i += 0.001) {
		te = exp(-i * thePhylo.stats_maxTreeBRL());
		te2 = 1 - exp(-i * thePhylo.stats_maxTreeBRL());
		W = te;
		W = W | te2;
		//cout << W << endl;
		//cout << i << endl;
		//V = (i * i) * tempV;
		//V = tempTree.matrix_standardize_UNQUAL((i * i) * tempV);	

		V = phylogenyTransformOU(tempV, i);
		//V = transformMatrixWithLambda(V, 1, false);
	//cout << V << endl;
		if((W.t() * V.i() * W).Determinant() != 0) { 
			mean = (W.t() * V.i() * W).i() * (W.t() * V.i() * y);
			//cout << mean << endl;
			Matrix error = (1.0 / y.Nrows()) * ((y - W * mean).t() * V.i() * (y - W * mean));	
			double U = (y.Nrows() * (1 + log(2 * PI * error(1,1))) + log(V.Determinant()));
			//cout << mean.t() << ", error = " << error(1,1) << ", U = " << U << endl;
			if(U < tempU) {tempU = U; min = i;}
		}
	}
*/
	cout << min << endl;

	thePhylo.getVECTORS_BUTLER_KING();
		
		//cout << (1.0/V.Ncols())*(y - W * (W.t() * V.i() * W).i() * (W.t() * V.i() * y)).t()*V.i()*(y - W * (W.t() * V.i() * W).i() * (W.t() * V.i() * y));
		//cout << (W.t() * V.i() * W).i() << endl;
	//}
}


/*Matrix alphaSelection(const ColumnVector &y, const Matrix &V, const Matrix &W)
{
	Matrix temp = (W.t() * V.i() * W).i() * (W.t() * V.i() * y);
	return temp;
}

Matrix alphaErrorSelection(const ColumnVector &y, const Matrix &V, const Matrix &W, const Matrix &alpha)
{
	Matrix e = (1.0 / V.Nrows()) * (y - alpha * W).t() * V.i() * (y - alpha * W); 
	return e;
}

Matrix powerMatrix (const Matrix & y, double power)
{
	Matrix temp = y;
	for(int i = 1; i <= y.Ncols(); i++)
		for(int j = 1; j <= y.Nrows(); j++)
			temp(i,j) = pow(y(i,j), power);
	return temp;
}

Matrix exponentMatrix (const Matrix & y)
{
	Matrix temp = y;
	for(int i = 1; i <= y.Ncols(); i++)
		for(int j = 1; j <= y.Nrows(); j++)
			temp(i,j) = exp(y(i,j));
	return temp;
}

Matrix MaxLikSelection(const ColumnVector & y, const Matrix & V, const Matrix & W, const Matrix & alpha, const Matrix & error)
{
	Matrix lik2, lik3, power;
	lik2 = (-0.5 * error.i()) * ((y - alpha * W).t() * V.i() * (y - alpha * W));
	power = powerMatrix(2.0 * PI * error, y.Nrows()/2.0);
	lik3 = power.i() * (1.0/sqrt(V.Determinant())) * exponentMatrix(lik2);
	return lik3;
}

double selection(const ColumnVector & x, const Matrix & V, double MAX, double MIN, bool showML) 
{
	double MLMzero = 0.0, MLMone = 0.0, MLm = 0.0, b = 0.0, MLMAX = 0.0, start = 0.0,  tempselection = MIN;
	Matrix selectionV, selectionW, a, e, ML;
	ColumnVector y(x.Nrows()); y = 1.0;
		
	std::vector <double> MLvector;
	double progress = MIN;

	cout << "Estmating ML  of selection:" << endl;
	for(double theselection = MIN; theselection <= MAX + 0.001; theselection += 0.001) {
		selectionV = phylogenyTransformOU(V, theselection);
		selectionW = weightMatrix(aPhylogeny, theselection);

		a = alphaSelection(x, selectionV, selectionW); 	
		//e = alphaErrorSelection(x, selectionV, selectionW, a);	
		//ML = MaxLikSelection(x, selectionV, selectionW, a, e);  

		cout << a << e << ML << endl;
	
	//	if(selection == MIN) {MLMAX = ML; tempselection = selection; start = ML;}
	//	if(ML >= MLMAX) {MLMAX = ML; tempselection = selection;}

		//cout << lambda << " " << progress << endl;
		if(theselection > progress) {
			printProgress(theselection, MAX);
			progress += (250.0 * 0.001);
		}

		if(showML == true) cout << theselection << " | " << ML << " (" << ML - MLMzero << ", " << ML - MLMone << ") ->" << MLMAX << endl;
	}

	//cout << setprecision(4) << "selection = " << tempselection << ", " << (char)224 << " = " << alpha(x, phylogenyTransformOU(V, tempselection)) << ", " << (char)229 << (char)253 << "(" << (char)224 << ") = " << alphaError(x, phylogenyTransformOU(V, tempselection), alpha(x, phylogenyTransformOU(V, tempselection))) << ", ln(L) = " << MLMAX << endl;
	
	
	//cout << "  " << char(192) << char(196) << char(196) << char(194) << char(196) << char(170) << " Equal to 1? Chi-square = " << -2.0*(MLMzero - MLMAX) << ", d.f. = 1, p = " << CHIsq_Prob(-2.0*(MLMzero - MLMAX), 1) << ", ln(L) = " << MLMzero << endl;
	//cout << "     " << char(192) << char(196) << char(170) << " Equal to 0? Chi-square = " << -2.0*(MLMone - MLMAX) << ", d.f. = 1, p = " << CHIsq_Prob(-2.0*(MLMone - MLMAX), 1) << ", ln(L) = " << MLMone << endl;
	
	//cout << endl << "    a) lambda = one? MLE = " << MLMzero << " | Chi-square = " << -2.0*(MLMzero - MLMAX) << ", p = " << CHIsq_Prob(-2.0*(MLMzero - MLMAX), 1) << endl;
	//cout << "    lambda differs from one? MLE = " << MLMone << " | Chi-square = " << -2.0*(MLMone - MLMAX) << ", p = " << CHIsq_Prob(-2.0*(MLMone - MLMAX), 1) << endl;

	return tempselection;
}





double Phylogeny::stats_Bulter_selection(const ColumnVector & y, bool showDetails)
{
	bool asButler = false;
	Phylogeny tempTree = aPhylogeny;
	Matrix A = tempTree.format_phylogeny_to_matrix();
	//double val = lambda(y, A, 2.0, 0.00001, showDetails, asButler);
	selection(y, transformMatrixWithLambda(A, 1.0, asButler), 15, 0.001, showDetails);
	
	return 0.0;
	//return maxOU(y, A, W, 1.1, 0.0, showDetails);
}


*/


//-------------------------------------------------------------------------------------
float BRL(const tree<TIP>& tr, tree<TIP>::pre_order_iterator it)
{
	tree<TIP>::pre_order_iterator start = tr.begin();
	float BRL = 0;
	while (it != start) {
		BRL += (*it).BRL;
		it = tr.parent(it);
	}
	return BRL;
}

void Phylogeny::transform_phylogenyPagel() 
{
	Phylogeny temp = aPhylogeny;
	temp.transform_phylogenyConstant();
	tree<TIP> tr = temp.aPhylogeny;
	float maxDepth = maxTreeDepthofTree(aPhylogeny), subtreeDepth = 0, nodeDepth = 0, totalDepth = 0;
	tree<TIP>::pre_order_iterator it = tr.begin(), end = tr.end();
	++it; //skip root node
	while(it != end) {
		subtreeDepth = maxTreeDepthofTree(tr.subtree(it));
		nodeDepth = BRL(tr,it);
		totalDepth = subtreeDepth + nodeDepth - 1;
		//if(subtreeDepth = !nodeDepth)  (*it).BRL = maxDepth - nodeDepth;
		(*it).BRL = (maxDepth - totalDepth);
		++it;
	}
	aPhylogeny = tr;
}

Phylogeny Phylogeny::transform_phylogenyConstant() 
{
	tree<TIP> tr = aPhylogeny; 
	tree<TIP>::pre_order_iterator it = tr.begin(), end = tr.end();
	while(it != end) {
		(*it).BRL = 1.0;
		++it;
	}
	aPhylogeny = tr;
	return tr;
}

void Phylogeny::transform_scale_Tree(double theScale) 
{
	Phylogeny temp = aPhylogeny;
	double MAXBL = temp.stats_maxTreeBRL();

	tree<TIP> tr = aPhylogeny; 
	tree<TIP>::pre_order_iterator it = tr.begin(), end = tr.end();
	while(it != end) {
		(*it).BRL = (theScale / MAXBL) * ((*it).BRL);
		++it;
	}
	aPhylogeny = tr;
}



Phylogeny Phylogeny::transform_phylogenyPagel(Phylogeny aPhylo) 
{
	aPhylo.transform_phylogenyConstant();
	tree<TIP> tr = aPhylo.aPhylogeny;
	double maxDepth = maxTreeDepthofTree(tr), subtreeDepth = 0, nodeDepth = 0, totalDepth = 0;
	tree<TIP>::pre_order_iterator it = tr.begin(), end = tr.end();
	++it; //skip root node
	while(it != end) {
		subtreeDepth = maxTreeDepthofTree(tr.subtree(it));
		nodeDepth = BRL(tr,it);
		totalDepth = subtreeDepth + nodeDepth - 1.0;
		if(subtreeDepth == nodeDepth)  (*it).BRL = maxDepth - nodeDepth;
		else (*it).BRL = (maxDepth - totalDepth);
		++it;
	}
	Phylogeny aPhy = tr;
	return aPhy;
}


//-------------------------------------------------------------------------------------


bool Phylogeny::stats_isUltrametric()
{
	Phylogeny aPhylo = aPhylogeny;
	std::vector< std::string > theSpecies = aPhylo.getAllSpecies();
	double temp = aPhylo.stats_speciesTotalBRL(theSpecies[0]);
	for(int i = 1; i < theSpecies.size(); i++) if(temp != aPhylo.stats_speciesTotalBRL(theSpecies[i])) return false;
	return true;
}

bool Phylogeny::stats_anyPolytomies()
{
	Phylogeny aPhylo = aPhylogeny;
	return ((aPhylo.stats_numberOfSpecies()*2-1) > aPhylo.stats_size());
}

int Phylogeny::stats_numberOfPolytomies() 
{
	tree<TIP>::pre_order_iterator it = aPhylogeny.begin(), end = aPhylogeny.end();
	int numPoly = 0;

	while(it != end) {
		if((*it).node && (aPhylogeny.number_of_children(it) > 2)) numPoly++;
		++it;
	}
	return numPoly;
}

float Phylogeny::stats_Nbar()
{
	Phylogeny tempPhylo = aPhylogeny;
	tree<TIP> tr = aPhylogeny;
	tree<TIP>::pre_order_iterator it = tr.begin(), end = tr.end();
	int numSpecies = 0;
	float N = 0;
	
	if(tempPhylo.stats_anyPolytomies()) errorACTION("Phylogeny contains polytomies", "Nbar cannot be calculated, exited to system", true); 

	while(it != end) {
		if(!(*it).node) {
			N += tr.depth(it);
			numSpecies++;
		}
		++it;
	}
	
	return N/numSpecies;
}

float Phylogeny::stats_Nvar()
{
	Phylogeny tempPhylo = aPhylogeny;
	tree<TIP> tr = aPhylogeny;
	tree<TIP>::pre_order_iterator it = tr.begin(), end = tr.end();
	int numSpecies = 0;
	float N = 0.0, varN = 0.0;
	
	if(tempPhylo.stats_anyPolytomies()) errorACTION("Phylogeny contains polytomies", "Nvar cannot be calculated, exited to system", true); 

	while(it != end) {
		if(!(*it).node) {
			N += tr.depth(it);
			numSpecies++;
		}
		++it;
	}
	
	N = N/numSpecies;

	it = tr.begin();
	while(it!=end) {
		if(!(*it).node) varN += ((tr.depth(it)-N)*(tr.depth(it)-N));
		++it;
	}
	
	return varN/numSpecies;
}

int numberOfSpeciesonTREE(const tree<TIP>& tr)
{
	tree<TIP>::pre_order_iterator it = tr.begin(), end = tr.end();
	int i = 0;
	while(it != end) {
		if((*it).node == false) ++i;
		++it;
    }
    return i;
}

float Phylogeny::stats_C()
{
	Phylogeny tempPhylo = aPhylogeny;
	tree<TIP> tr = aPhylogeny, sub1, sub2;
	tree<TIP>::pre_order_iterator it = tr.begin(), end = tr.end();
	std::vector<tree<TIP>::pre_order_iterator> nodes, nodesTMP, nodes_mem;
	float C = 0.0;
	
	if(tempPhylo.stats_anyPolytomies()) errorACTION("Phylogeny contains polytomies", "C cannot be calculated, exited to system", true); 

	while(it != end) {
		if((*it).node) 
			if((*it).id != 0) 
				nodes.push_back(it); 
		++it;
	}
	
	nodes_mem = nodes; // used for I prime calculation

	bool stillPairs = false;
	while(!nodes.empty()) {
		stillPairs = false;
		for(int z = 0; z < nodes.size(); z++) {//checks for pairs
			for(int x = 0; x < nodes.size(); x++) 
				if((z != x) && (tr.parent(nodes[z]) == tr.parent(nodes[x]))) 
					{stillPairs = true; break;}
			if(stillPairs == true) break;
		}
		
		if(stillPairs == true) {
			for(int j = 0; j < nodes.size(); j++) {
				for(int z = 0; z < nodes.size(); z++) {
					if(z != j) {
						if(tr.parent(nodes[z]) == tr.parent(nodes[j])) {
							sub1 = tr.subtree(nodes[z]); sub2 = tr.subtree(nodes[j]);
							if(numberOfSpeciesonTREE(sub1) >= numberOfSpeciesonTREE(sub2)) 
								C += (numberOfSpeciesonTREE(sub1) - numberOfSpeciesonTREE(sub2));
							else 
								C += (numberOfSpeciesonTREE(sub2) - numberOfSpeciesonTREE(sub1));
							for(int k = 0; k < nodes.size(); k++) 
								if((k != j) && (k != z)) 
									nodesTMP.push_back(nodes[k]); 
							nodes = nodesTMP; nodesTMP.clear();
							break; 
						}
					}
				}
			}
		}
		else {
			for(int j = 0; j < nodes.size(); j++)
				C += (numberOfSpeciesonTREE(tr.subtree(nodes[j])) - 1);
			break;
		}
	}

	return (2 / (numberOfSpeciesonTREE(tr) * (numberOfSpeciesonTREE(tr) - 3) + 2)) * C;
}


float Phylogeny::stats_B1()
{
	Phylogeny tempPhylo = aPhylogeny;
	tree<TIP> tr = aPhylogeny;
	tree<TIP>::pre_order_iterator it = tr.begin(), end = tr.end(), it3 = tr.begin();
	std::vector< tree<TIP> > treeVector;
	float B1 = 0.0, MAX = 0.0;

	if(tempPhylo.stats_anyPolytomies()) errorACTION("Phylogeny contains polytomies", "C cannot be calculated, exited to system", true); 

	while(it != end) {
		if((*it).node) {
			if((*it).id != 0)
				treeVector.push_back(tr.subtree(it)); 
		}
		++it;
	}

	for(int p = 0; p < treeVector.size(); p++) {
		it3 = (treeVector[p]).begin();
		MAX = 0.0;
		while(it3 != (treeVector[p]).end()) {
			if(!(*it3).node) 
				if((treeVector[p]).depth(it3) > MAX) MAX = (treeVector[p]).depth(it3);
			++it3;
		}
		B1 += (1/MAX);
	}

	return B1;
}


float Phylogeny::stats_meanIprime()
{
	Phylogeny tempPhylo = aPhylogeny;
	tree<TIP> tr = aPhylogeny;

	std::vector<tree<TIP>::pre_order_iterator> nodes, nodes_mem;
	tree<TIP>::sibling_iterator biggest;
	tree<TIP>::pre_order_iterator it = tr.begin(), end = tr.end(), it3 = tr.begin();

	float Iprime = 0.0, IprimeM = 0.0;
	int S = 0, m = 0, B = 0;

	while(it != end) {
		if((*it).node) 
			if((*it).id != 0) 
				nodes.push_back(it); 
		++it;
	}
	
	for(int i = 0; i < nodes.size(); i++)
		if(numberOfSpeciesonTREE(tr.subtree(nodes[i])) >= 4) nodes_mem.push_back(nodes[i]);

	nodes = nodes_mem;

	for(i = 0; i < nodes.size(); i++) {
		S = numberOfSpeciesonTREE(tr.subtree(nodes[i]));
		if(fmod(S,2) == 0) m = S / 2; else m = (S/2)+1;

		biggest = tr.child(nodes[i], 0);
		while(tr.index(biggest) != tr.number_of_children(nodes[i])-1) {
			if(numberOfSpeciesonTREE(tr.subtree(biggest)) > B) B = numberOfSpeciesonTREE(tr.subtree(biggest));
			++biggest;
		}
		Iprime = (B - m) / (S - m - 1);
		if(isEven(S)) Iprime = (Iprime * (S-1)) / S;
		IprimeM += Iprime;
		Iprime = 0;	B = 0; S = 0; m = 0;
	}

	return IprimeM/nodes.size();
}

long unsigned Phylogeny::stats_numberOfResolutions(int aNumber)
{
	long unsigned solves;
	if(aNumber <= 2) solves = 0;
	else solves = factorial((2*aNumber-3))/(pow(2,aNumber-2)*factorial(aNumber-2));
	return solves;
}

long unsigned Phylogeny::stats_numberOfResolutionsForPolytomies()
{
	tree<TIP>::pre_order_iterator it = aPhylogeny.begin(), end = aPhylogeny.end();
	long unsigned resolutions = 1, solves = 0;
	Phylogeny temp, tempB = aPhylogeny;
	
	if(tempB.stats_anyPolytomies() == false) return 0;

	while(it != end) {
		if((*it).node && (aPhylogeny.number_of_children(it) > 2)) {
			solves = temp.stats_numberOfResolutions(aPhylogeny.number_of_children(it));
			resolutions *= solves;
		}
		solves = 0;
		++it;
	}
	return resolutions;
}

//*********** Webb 2000 Am. Nat. | START

tree<TIP>::pre_order_iterator getSpeciesPointer(const tree<TIP>& tr, std::string aSpecies) 
{
	tree<TIP>::pre_order_iterator it = tr.begin(), end = tr.end();
	std::vector<tree<TIP>::pre_order_iterator> speciesListPointers;
	while(it != end) {
		if((*it).species == aSpecies) return it;
		++it;
    }
	return tr.begin();	
}

/*std::vector<tree<TIP>::pre_order_iterator> getSpeciesPointerFromString(std::vector<tree<TIP>::pre_order_iterator> pointerVector, std::vector<std::string> speciesVector) 
{
	tree<TIP>::pre_order_iterator it;
	std::vector<tree<TIP>::pre_order_iterator> speciesListPointers;
	for(int i = 0; i < speciesVector.size(); i++) speciesListPointers.push_back(findSpeciesPointer(pointerVector, speciesVector[i]));
	return speciesListPointers;
}*/

std::vector<std::string> random_species_list_vector(std::vector<std::string> masterList, int hostRange)
{
	std::vector<std::string> randomList;
	int N = masterList.size();
	for(int i = 0 ; i < hostRange; i++) 
		randomList.push_back(masterList[rand()%(N)]);
	return randomList;
}

double nodalDistanceToNode(const tree<TIP>& aTree, tree<TIP>::pre_order_iterator sp, tree<TIP>::pre_order_iterator aNode) 
{
	double NDTN = 0.0;
	while(sp != aNode) {
		NDTN += (*sp).BRL;
		sp = aTree.parent(sp);
	}
	return NDTN;
}



tree<TIP>::pre_order_iterator webb_sharedNode(const tree<TIP>& aTree, tree<TIP>::pre_order_iterator sp1, tree<TIP>::pre_order_iterator sp2) 
{
	tree<TIP>::pre_order_iterator temp = sp2, root = aTree.begin();
	while(sp1 != root) {
		while(sp2 != root) {
			if(sp1 == sp2) return sp1;
			sp2 = aTree.parent(sp2);
		}
		sp1 = aTree.parent(sp1);
		sp2 = temp;
	}
	return root;
}

double webb_PairwiseNodalDistance(const tree<TIP>& aTree, tree<TIP>::pre_order_iterator sp1, tree<TIP>::pre_order_iterator sp2) 
{
	tree<TIP>::pre_order_iterator sharedNode;
	sharedNode = webb_sharedNode(aTree, sp1, sp2);
	double PNDsp1 = 0.0, PNDsp2 = 0.0;
	PNDsp1 = nodalDistanceToNode(aTree, sp1, sharedNode);
	PNDsp2 = nodalDistanceToNode(aTree, sp2, sharedNode);
	return (PNDsp1 + PNDsp2)/2;
}

double webb_meanPairwiseNodalDistance(const tree<TIP>& aTree, std::vector< std::string > aSpeciesList, bool withBRL) 
{
	std::vector<tree<TIP>::pre_order_iterator> speciesListPointers, allspeciesListPointers;
	allspeciesListPointers = speciesPointerVector(aTree);
	speciesListPointers = getSpeciesPointerFromString(allspeciesListPointers, aSpeciesList);
	int temp = 1, N = aSpeciesList.size();
	double Nc = (N * (N - 1))/2, NRI = 0.0;
	for(int i = 0; i < N; i++) {
		for(int j = temp; j < N; j++) NRI += webb_PairwiseNodalDistance(aTree, speciesListPointers[i], speciesListPointers[j]);
		temp++;
	}
	return NRI/Nc;
}

double webb_MAX_PairwiseNodalDistance(const tree<TIP>& aTree, std::vector< std::string > aSpeciesList) 
{
	std::vector< std::string > randomList, speciesList = speciesVector(aTree);
	int N = aSpeciesList.size();
	double max = 0.0, temp = 0.0;
	for(int i = 0; i < MAXRESAMPLES; i++) 
	{
		randomList = random_species_list_vector(speciesList, N);		
		temp = webb_meanPairwiseNodalDistance(aTree, randomList, true);
		if(temp > max) max = temp;
	}
	return max;
}

double webb_closestRelative(const tree<TIP>& aTree, std::vector< std::string > aSpeciesList) 
{
	double temp = 10000000000.0, nodeDist = 0.0, sum = 0.0;
	int N = aSpeciesList.size();
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			if(i != j) {
				nodeDist = webb_PairwiseNodalDistance(aTree, getSpeciesPointer(aTree, aSpeciesList[i]), getSpeciesPointer(aTree, aSpeciesList[j]));
				if (nodeDist < temp) temp = nodeDist;
			}
		}
		sum += temp; temp = 1000000.0; 
	}
	return sum/N;
}

double webb_MAX_closestRelative(const tree<TIP>& aTree, std::vector< std::string > aSpeciesList) 
{
	std::vector< std::string > randomList, speciesList = speciesVector(aTree);
	int N = aSpeciesList.size();
	double max = 0.0, temp = 0.0;
	for(int i = 0; i < MAXRESAMPLES; i++) 
	{
		randomList = random_species_list_vector(speciesList, N);		
		temp = webb_closestRelative(aTree, randomList);
		if(temp > max) max = temp;
	}
	return max;
}

double webb_meanPairwiseNodalDistance_species(const tree<TIP>& aTree, std::vector< std::string > aSpeciesList, std::string aSpecies) 
{
	std::vector<tree<TIP>::pre_order_iterator> speciesListPointers, allspeciesListPointers;
	allspeciesListPointers = speciesPointerVector(aTree);
	speciesListPointers = getSpeciesPointerFromString(allspeciesListPointers, aSpeciesList);
	tree<TIP>::pre_order_iterator aSpeciesPtr = getSpeciesPointer(aTree, aSpecies);
	int N = aSpeciesList.size();
	double Nc = (N * (N - 1))/2, NRI = 0.0;
	for(int i = 0; i < N; i++) {
		NRI += webb_PairwiseNodalDistance(aTree, aSpeciesPtr, speciesListPointers[i]);
	}
	return NRI/Nc;
}

double webb_MAX_PairwiseNodalDistance_species(const tree<TIP>& aTree, std::vector< std::string > aSpeciesList, std::string aSpecies) 
{
	std::vector< std::string > randomList, speciesList = speciesVector(aTree);
	int N = aSpeciesList.size();
	double max = 0.0, temp = 0.0;
	for(int i = 0; i < MAXRESAMPLES; i++) 
	{
		randomList = random_species_list_vector(speciesList, N);		
		temp = webb_meanPairwiseNodalDistance_species(aTree, randomList, aSpecies);
		if(temp > max) max = temp;
	}
	return max;
}

double webb_closestRelative_species(const tree<TIP>& aTree, std::vector< std::string > aSpeciesList, std::string aSpecies) 
{
	double temp = 100000000000000000.0, nodeDist = 0.0;
	int N = aSpeciesList.size();
	tree<TIP>::pre_order_iterator aSpeciesPtr = getSpeciesPointer(aTree, aSpecies);
	for(int i = 0; i < N; i++) {
		nodeDist = webb_PairwiseNodalDistance(aTree, aSpeciesPtr, getSpeciesPointer(aTree, aSpeciesList[i]));
		if (nodeDist < temp) temp = nodeDist;
	}
	return temp/N;
}

double webb_MAX_closestRelative_species(const tree<TIP>& aTree, std::vector< std::string > aSpeciesList, std::string aSpecies) 
{
	std::vector< std::string > randomList, speciesList = speciesVector(aTree);
	int N = aSpeciesList.size();
	double max = 0.0, temp = 0.0;
	for(int i = 0; i < MAXRESAMPLES; i++) 
	{
		randomList = random_species_list_vector(speciesList, N);		
		temp = webb_closestRelative_species(aTree, randomList, aSpecies);
		if(temp > max) max = temp;
	}
	return max;
}

double Phylogeny::stats_NRIspecies_Webb(std::vector< std::string > aSpeciesList, std::string aSpecies, bool resample) 
{
	Phylogeny thePhylo = aPhylogeny;
	if(resample == false) return (webb_meanPairwiseNodalDistance_species(aPhylogeny, aSpeciesList, aSpecies)/thePhylo.stats_maxTreeBRL());

	double MPND = webb_meanPairwiseNodalDistance_species(aPhylogeny, aSpeciesList, aSpecies), maxPND =  webb_MAX_PairwiseNodalDistance_species(aPhylogeny, aSpeciesList, aSpecies);
	return (1.0 - MPND/maxPND);
}

double Phylogeny::stats_NTIspecies_Webb(std::vector< std::string > aSpeciesList, std::string aSpecies, bool resample) 
{
	Phylogeny thePhylo = aPhylogeny;
	if(resample == false) return (webb_closestRelative_species(aPhylogeny, aSpeciesList, aSpecies)/thePhylo.stats_maxTreeBRL());

	double CR = webb_closestRelative_species(aPhylogeny, aSpeciesList, aSpecies), maxCR = webb_MAX_closestRelative_species(aPhylogeny, aSpeciesList, aSpecies);
	return (1.0 - CR/maxCR);
}

double Phylogeny::stats_NRI_Webb(std::vector< std::string > aSpeciesList, bool withBRL) 
{
	double MPND = webb_meanPairwiseNodalDistance(aPhylogeny, aSpeciesList, true), maxPND =  webb_MAX_PairwiseNodalDistance(aPhylogeny, aSpeciesList);
	return (1.0 - MPND/maxPND);
}

double Phylogeny::stats_NTI_Webb(std::vector< std::string > aSpeciesList, bool withBRL) 
{
	double CR = webb_closestRelative(aPhylogeny, aSpeciesList), maxCR = webb_MAX_closestRelative(aPhylogeny, aSpeciesList);
	return (1.0 - CR/maxCR);
}

//*********** Webb 2000 Am. Nat. | END

/////////////tree generation functions
Phylogeny Phylogeny::generate_randomBifurcatedTree(int numberOfSpecies)
{
	Phylogeny newTree;
	std::vector< std::string > speciesList;
	for(int i = 0; i < numberOfSpecies; i++) speciesList.push_back(intToString(i));
	return newTree;
}

Phylogeny Phylogeny::transform_NEWICK_to_Phylogeny(std::string phylogeny)
{
	std::string sp; 
	std::vector<std::string> testPhylo;
	int countNode = 0;
	float BRL = 0.0, BRL2 = 0.0;
	
	if(!parenOK(phylogeny)) nrerror("phylogy format incorrect, unbalanced paraenthese...");	//write function that removes all spaces from phylogeny string!

	//puts file contents into a vector
	for(int i = 0; i < phylogeny.size(); i++) {
		if(phylogeny[i] == '(') testPhylo.push_back("(");
		if((phylogeny[i] != '(') && (phylogeny[i] != ')') && (phylogeny[i] != ',') && (phylogeny[i] != ';') && (phylogeny[i] != ':')) 
			sp += phylogeny[i];
		if((phylogeny[i] == ',') || (phylogeny[i] == ')') || (phylogeny[i] == ':')) {
				if(sp.size()) testPhylo.push_back(sp);
				sp = "";
				if(phylogeny[i] == ')') testPhylo.push_back(")");
				if(phylogeny[i] == ':') testPhylo.push_back(":");
		}
	}
	testPhylo.push_back(":"); testPhylo.push_back("0.0");

	tree<TIP> treeA;
	tree<TIP>::pre_order_iterator rootI;
	std::vector<tree<TIP>::pre_order_iterator> pointerList;	
	
	TIP root, tp; 
	root.node = true; root.id = 0; root.DATA_x = 0.0; root.ignore = true; root.BRL = 0; root.BRL_y = -5.0; root.species = "node"; root.optima = 0; //root.parent = treeA.begin();
	tp.node = false; tp.id = 1; tp.DATA_x = 0.0; tp.ignore = true; tp.BRL = 0.0; tp.BRL_y = -5.0; tp.species = "tip"; tp.optima = 0; //tp.parent = treeA.begin();

	//construct phylogeny tree with branchlengths
	for(i = 0; i < testPhylo.size(); i++) {	
		if(testPhylo[i] == "(") {
			root.id = countNode++; // root id always zero
			//root.parent = rootI;
			if(!treeA.size()) rootI = treeA.insert(treeA.begin(), root);
			else rootI = treeA.append_child(rootI,root);
			pointerList.push_back(rootI);
		}
		else if(testPhylo[i] == ":" && testPhylo[i-1] != ")") {
			std::istringstream sin(testPhylo[i+1]); sin >> BRL;	
			tp.BRL = BRL; tp.id = countNode++; tp.species = testPhylo[i-1];
			//tp.parent = rootI;
			treeA.append_child(rootI,tp);
		} 
		else if(testPhylo[i] == ")") {
			if(testPhylo[i+1] == ":") {
				std::istringstream sin(testPhylo[i+2]); sin >> BRL2;
				(*rootI).BRL = BRL2;	
			}
			pointerList.erase(pointerList.end()-1,pointerList.end()); 
			rootI = pointerList[pointerList.size()-1];
		}
	}
	
	Phylogeny final = treeA;
	return final;
}

Phylogeny Phylogeny::addSpeciesPolitomy(std::string aSpecies, int number)
{
	Phylogeny aPhylo = aPhylogeny;
	std::string temp = aPhylo.format_phylogeny_to_NEWICK(true), theString, aName, aNother, final;
	std::vector < std::string > aVectorString;

	theString += "(";
	for(int i = 0; i < number; i++) {
		theString += aSpecies + ":0.000001";
		if(i != number-1) theString += ",";
	}
	theString += ")";

	for(int j = 0; j <= temp.size(); j++) {
		aName += temp[j];
		if(temp[j] == '(' || temp[j] == ':' || temp[j] == ',' || temp[j] == ')') {
			aName = aName.erase(aName.size()-1);
			aVectorString.push_back(aName);
			aNother = temp[j];
			aVectorString.push_back(aNother);
			aName = "";
		}
	}

	for(int k = 0; k < aVectorString.size(); k++) {
		if (aVectorString[k] == aSpecies) final += theString; 
		else final += aVectorString[k];
	}
	
	aPhylo = aPhylo.transform_NEWICK_to_Phylogeny(final + ";");
	return aPhylo;
}

Matrix GLS_tranformation(ColumnVector theDatas, Matrix theMatrix, bool OU, double alpha)
{
	//ColumnVector theDatas(theData.size()); for(int j = 1; j < theData.size()+1; j++) theDatas.Row(j) = theData[j-1];
	//Matrix theMatrix = aPhylo.format_phylogeny_to_matrix(), a(3,3);

	//theMatrix = aPhylo.matrix_standardize_UNQUAL(theMatrix);
		/* Desman & Beavers 1976 example
	a << 120000 << 230 << 10 << 230 << 1000 << 1 << 10 << 1 << 0.5;
	SymmetricMatrix S(3); for(int i = 1; i < 3+1; i++) for(int j = 1; j < 3+1; j++) S(i,j) = a(i,j);
	*/
	//cout << setprecision(3) << theMatrix << endl;


	SymmetricMatrix S(theDatas.Nrows()); 
	for(int i = 1; i < theDatas.Nrows()+1; i++) 
		for(int j = 1; j < theDatas.Nrows()+1; j++) S(i,j) = theMatrix(i,j);
	LowerTriangularMatrix L = Cholesky(S);
	DiagonalMatrix D;
	Matrix V, U, result;


	//method 1
	//SVD(L, D, U, V);
	//result = U * D * U.t();
	//return result.i() * theDatas;
	

	//cout << setprecision(3) << theMatrix << endl << result << (theMatrix == (result * result.t())) << endl;

	

//	SymmetricMatrix S(theData.size()); for(int i = 1; i < theData.size()+1; i++) for(int j = 1; j < theData.size()+1; j++) S(i,j) = theMatrix(i,j);

//	Matrix V, U, result;
//	DiagonalMatrix D;

//method 3
//	EigenValues(S, D, V);
//	for(i = 1; i < D.Ncols() + 1; i++) D(i) = sqrt(D(i));
//	result = V * D * V.t();

	//method 2 ADAMS'
	SVD(S, D, U, V);
	for(i = 1; i < D.Ncols() + 1; i++) D(i) = sqrt(D(i));
	result = U * D * V.t();
	return result.i() * theDatas;
}

Matrix GLS_tranformation_2(Matrix theDatas, Matrix theMatrix)
{
	SymmetricMatrix S(theDatas.Nrows()); 
	for(int i = 1; i < theDatas.Nrows()+1; i++) 
		for(int j = 1; j < theDatas.Nrows()+1; j++) S(i,j) = theMatrix(i,j);
	LowerTriangularMatrix L = Cholesky(S);
	DiagonalMatrix D;
	Matrix V, U, result;

	SVD(S, D, U, V);
	for(i = 1; i < D.Ncols() + 1; i++) D(i) = sqrt(D(i));
	result = U * D * V.t();
	return result.i() * theDatas;
}

Phylogeny Phylogeny::randomlyResolvePolytomies(double internodeBRL)
{
	tree<TIP> temp = aPhylogeny;
	
	tree<TIP>::pre_order_iterator it = temp.begin(), end = temp.end(), it2 = temp.begin(), end2 = temp.end(), tempIt, nodeA, nodeB;
	std::vector < tree<TIP>::pre_order_iterator > theTHEY;

	TIP newNode;
	newNode.BRL = internodeBRL; newNode.id = 666; newNode.node = true; newNode.species = "node";
	
	Phylogeny somePhylo = temp;

	while (somePhylo.stats_numberOfPolytomies() > 0) {
		it = temp.begin();	end = temp.end();
		it2 = temp.begin();	end2 = temp.end();
		//cout << "--------->" << somePhylo.stats_numberOfPolytomies() << endl;
		while(it != end) {
			if((*it).node && (temp.number_of_children(it) > 2)) {
				while(it2 != end2) { //collects all children at the it polytomy
					if(it == temp.parent(it2)) theTHEY.push_back(it2);
					++it2;
				}
				it2 = temp.begin();
				end2 = temp.end();
				//std::cout << temp.number_of_children(it) << " " << theTHEY.size() << std::endl;
				//for(int j = 0; j < theTHEY.size(); j++) std::cout << "               " << (*theTHEY[j]).species << std::endl;
				
				//randomly pick two children to resolve

				int aNUM = rand()%(theTHEY.size()), tempNUM;
				tempNUM = aNUM;

				nodeA = theTHEY[tempNUM];
				while(tempNUM == aNUM) aNUM = rand()%(theTHEY.size());
				nodeB = theTHEY[aNUM];

				tempIt = temp.append_child(it,newNode);
				temp.append_child(tempIt, nodeA);
				temp.append_child(tempIt, nodeB);
				temp.erase(nodeA);
				temp.erase(nodeB);

				theTHEY.clear();
			}
			++it;
		}
	somePhylo = temp;
	//somePhylo.print_phylogeny();
	//std::cout << somePhylo.format_phylogeny_to_NEWICK(true) << std::endl;
	}		
	
	//temp.sort(temp.begin(), temp.end(), true);
	//somePhylo = temp;
	return somePhylo;
}


std::vector<tree<TIP>::pre_order_iterator> killElementOfVector(std::vector<tree<TIP>::pre_order_iterator> theVector, int element)
{
	std::vector<tree<TIP>::pre_order_iterator>::iterator it;
	int count = 0;
	for(it = theVector.begin(); it != theVector.end(); it++) 	{
		if(count == element) {
			theVector.erase(it);
			return theVector;
		}
		count++;
	}
	return theVector;
}


//ISSUES biased speciation towards species belonging early in IT iteration

Phylogeny Phylogeny::randomPhylogeny(int numSpecies, double speciationRate)
{
	tree<TIP> newTree;
	tree<TIP>::pre_order_iterator it, end, rootIt;

	Phylogeny finalTree;

	TIP root, tp; 
	int ID = 1, species = 0, randomNumber = 0;
	double rate = 1.0;
	root.node = true; root.id = 0; root.BRL = 0.0; root.species = "root";
	tp.node = false; tp.id = 1; tp.BRL = 0.0; tp.species = intToString(ID); 
	
	rootIt = newTree.insert(newTree.begin(), root);
	
	it = newTree.begin(); end = newTree.end();

	newTree.append_child(it,tp);
	ID++;
	tp.id = ID; tp.species = intToString(2);
	newTree.append_child(it,tp);
	
	finalTree = newTree;
	species = finalTree.stats_numberOfSpecies();
	
	bool exit = false;

	while(exit != true){
		it = newTree.begin(); end = newTree.end();
		while(it != end && exit != true) {
			if((!((*it).node)) && ((randomNumber = rand()%100) <= speciationRate*100)) {
				if((*it).BRL == 0 && (*it).species != "root") {
					(*it).BRL = rate; tp.BRL = -(rate);
				}
				else tp.BRL = 0.0;
				tp.species = (*it).species;	ID++; tp.id = ID; tp.species += intToString(ID); tp.node = false; 
				newTree.append_child(it,tp);
				tp.species = (*it).species;	ID++; tp.id = ID; tp.species += intToString(ID); tp.node = false; 
				newTree.append_child(it,tp);
				(*it).node = true; (*it).species = "node";
			}
			++it;
			finalTree = newTree;
			if(finalTree.stats_numberOfSpecies() == numSpecies) exit = true;
		}

		it = newTree.begin(); end = newTree.end();
			
		while(it != end) {
			if(!(*it).node) (*it).BRL += rate;
			++it;
				finalTree = newTree;
		}
	}

	bool stop = false;
	while(stop != true) {
		it = newTree.begin(); end = newTree.end();
		while(it != end) {
			if(!((*it).node)) {
				(*it).BRL += rate;
				if((randomNumber = rand()%100) <= speciationRate*100) stop = true;
			}
			++it;
		}
	}

	finalTree = newTree;
	return finalTree;
}

Phylogeny Phylogeny::transform_make_Ultrametric() 
{
	tree<TIP> somePhylo = aPhylogeny;
	Phylogeny aPhylo = aPhylogeny;
	double maxBRL = aPhylo.stats_maxTreeBRL(), spBRL = 0.0;
	tree<TIP>::pre_order_iterator theSpecies;
	std::vector< std::string > allSpecies = aPhylo.getAllSpecies();
	for(int i = 0; i < allSpecies.size(); i++) {
		spBRL = aPhylo.stats_speciesTotalBRL(allSpecies[i]);
		if((maxBRL-spBRL) != 0.0) {
			theSpecies = getSpeciesPointer(somePhylo, allSpecies[i]);
			(*theSpecies).BRL += (maxBRL-spBRL);
		}
	}
	Phylogeny thePhylo = somePhylo;
	return thePhylo;
}


//tree drawing functions


std::vector<tree<TIP>::pre_order_iterator> ignorePointerVector(const tree<TIP>& tr)
{
	tree<TIP>::pre_order_iterator it = tr.begin(), end = tr.end();
	std::vector<tree<TIP>::pre_order_iterator> ignoreListPointers;
	while(it != end) {
		if((*it).ignore == false && (*it).id != 0) ignoreListPointers.push_back(it);
		++it;
    }
	return ignoreListPointers;
}


tree<TIP> takeBRL_y(tree<TIP> & theTree) 
{
	std::vector<tree<TIP>::pre_order_iterator> allPTRX = ignorePointerVector(theTree);
	tree<TIP>::pre_order_iterator tempPTR, parentA, parentB;

//	cout << allPTRX.size() << endl;

	double y_distance = 0.0;
	for(int k = 0; k < allPTRX.size(); k++) {
		for(int j = 0; j < allPTRX.size(); j++) {
			if(k != j) {
				parentA = theTree.parent(allPTRX[k]); parentB = theTree.parent(allPTRX[j]);
				if((*parentA).id == (*parentB).id) {
					//cout << (*parentA).id << " " << (*parentB).id << endl;
					tempPTR = theTree.parent(allPTRX[k]);
					(*tempPTR).BRL_y = ((*allPTRX[k]).BRL_y + (*allPTRX[j]).BRL_y)/2.0;
					(*tempPTR).ignore = false;
					(*allPTRX[k]).ignore = true; (*allPTRX[j]).ignore = true;
					return theTree;
				}	
			}
		}
	}
	return theTree;
}

bool anyBRL_yNULL( tree<TIP> theTree) 
{
	tree<TIP>::pre_order_iterator it = theTree.begin(), itend = theTree.end();
	
	while(it != itend) {
		if ((*it).BRL_y == -5.0) return true;
		++it;
	}	
	return false;
}

bool anyPairs(tree<TIP> theTree) 
{
	tree<TIP>::pre_order_iterator it = theTree.begin(), itend = theTree.end();
	int count = 0;
	while(it != itend) {
		if (((*it).ignore == false) && ((*it).id != 0)) {
			//std::cout << "id " << (*it).id << " " << (*it).species << std::endl;
			count++;
		}
		++it;
	}	
	//cout << "ggg = " << count << endl;
	return ((count > 1) ? true: false);
}

void Phylogeny::getCoordinates() 
{
	tree<TIP> somePhylo = aPhylogeny;
	Phylogeny tempPhylo = aPhylogeny;
	
	std::vector<tree<TIP>::pre_order_iterator> allSpPointers = speciesPointerVector(somePhylo);
	tree<TIP>::pre_order_iterator tempPTR;

	double maxBRL = tempPhylo.stats_maxTreeBRL()/tempPhylo.stats_numberOfSpecies();


	double y_distance = 0.0;

	for(int k = 0; k < allSpPointers.size(); k++) {
		(*allSpPointers[k]).BRL_y = y_distance;
		(*allSpPointers[k]).ignore = false;
		y_distance += maxBRL;
	}

	tree<TIP>::pre_order_iterator it = somePhylo.begin(), itend = somePhylo.end();
	
	while(anyBRL_yNULL(somePhylo) == true) {
		it = somePhylo.begin(); itend = somePhylo.end();
		somePhylo = takeBRL_y(somePhylo);
	}
	aPhylogeny = somePhylo;
}


tree<TIP> findPairs_y(tree<TIP> & theTree) 
{
	std::vector<tree<TIP>::pre_order_iterator> allPTRX = ignorePointerVector(theTree);
	tree<TIP>::pre_order_iterator tempPTR, parentA, parentB, tempPA, tempPB, tempPC, it = theTree.begin(), itend = theTree.end();
	std::vector < tree<TIP>::pre_order_iterator > theThey;
	double max = 0.0, min = 0.0;
	int maxPTR, minPTR;

	double y_distance = 0.0, PA = 0.0, PB = 0.0;
	for(int k = 0; k < allPTRX.size(); k++) {
		for(int j = 0; j < allPTRX.size(); j++) {
			if(k != j) {
				//cout << k << " " << j << endl;
				if((*allPTRX[k]).id != 0) parentA = theTree.parent(allPTRX[k]); 
				else parentA = allPTRX[k];
				if((*allPTRX[j]).id != 0) parentB = theTree.parent(allPTRX[j]); 
				else parentB = allPTRX[j];
				
				if((*parentA).id == (*parentB).id) {
					//tempPC = theTree.parent(parentA);
					if(theTree.number_of_children(parentA) > 2) {
						while(it != itend) { //collects all children at the it polytomy
							if(parentA == theTree.parent(it)) {
								//std::cout << (*it).species << std::endl;
								theThey.push_back(it);
							}
							++it;
						}
						//pick the two most distant
						max = (*theThey[0]).BRL_y; maxPTR = 0;
						min = (*theThey[0]).BRL_y; minPTR = 0;
						(*theThey[0]).ignore = true;
						for(int u = 1; u < theThey.size(); u++) {
							if((*theThey[u]).BRL_y > max) {
								max = (*theThey[u]).BRL_y;
								maxPTR = u;
							}
							if((*theThey[u]).BRL_y < min) {
								min = (*theThey[u]).BRL_y;
								minPTR = u;
							}
							if((*theThey[u]).id != 0) (*theThey[u]).ignore = true;
						}
						tempPA = theThey[minPTR]; tempPB = theThey[maxPTR];
						while((*tempPA).id != 0) {
							PA += (*tempPA).BRL;
							tempPA = theTree.parent(tempPA);
						}
						while((*tempPB).id != 0) {
							PB += (*tempPB).BRL;
							tempPB = theTree.parent(tempPB);
						}
						cout << (PA - (*theThey[minPTR]).BRL) << " " << (*theThey[minPTR]).BRL_y << " " << (PB - (*theThey[maxPTR]).BRL) << " " << (*theThey[maxPTR]).BRL_y << " node" << endl;
					}
					else {
						tempPA = allPTRX[k]; tempPB = allPTRX[j];
						while((*tempPA).id != 0) {
							PA += (*tempPA).BRL;
							tempPA = theTree.parent(tempPA);
						}
						while((*tempPB).id != 0) {
							PB += (*tempPB).BRL;
							tempPB = theTree.parent(tempPB);
						}
						cout << (PA - (*allPTRX[k]).BRL) << " " << (*allPTRX[k]).BRL_y << " " << (PB - (*allPTRX[j]).BRL) << " " << (*allPTRX[j]).BRL_y << " node" << endl;
						if((*allPTRX[k]).id != 0) (*allPTRX[k]).ignore = true;
						if((*allPTRX[j]).id != 0) (*allPTRX[j]).ignore = true;
					}
				 
					return theTree;
				}	
				PA = 0.0;
				PB = 0.0;
			}
		}
	}
	return theTree;
}


void Phylogeny::showCoordinates()
{
	tree<TIP> somePhylo = aPhylogeny;	
	tree<TIP>::pre_order_iterator it = somePhylo.begin(), itend = somePhylo.end(), tempIt, tempIt2;
	double x = 0.0, y= 0.0, x2= 0.0, y2= 0.0;
	
	while(it != itend) {
		tempIt = it;
		while((*tempIt).id != 0) {
			x2 += (*tempIt).BRL;
			tempIt = aPhylogeny.parent(tempIt);
		}
		std::cout << (x2 - (*it).BRL) << " " << (*it).BRL_y << " " << x2 << " " << (*it).BRL_y << " " << (*it).species << std::endl;
		
		++it;
		x2 = 0;
	}	

	it = somePhylo.begin();

	while(it != itend) {
		(*it).ignore = false;
		++it;
	}
	//cout <<"bitch" << endl;
	while(anyPairs(somePhylo) == true) {
		
		somePhylo = findPairs_y(somePhylo);
	}

}


void Phylogeny::inputTaxaXdata(std::vector < std::string> theTAXA, std::vector < double> theDATA, std::vector < int > theOptima, std::vector< double > theVar)
{
	std::vector<tree<TIP>::pre_order_iterator> allTips = speciesPointerVector(aPhylogeny);
	//printStringVector(theTAXA);
	for(int i = 0; i < theTAXA.size(); i++) 
		for(int j = 0; j < theTAXA.size(); j++)
			if(theTAXA[i] == (*allTips[j]).species) {
				//std::cout << (*allTips[j]).species << std::endl;
				(*allTips[j]).DATA_x = theDATA[i];
				(*allTips[j]).optima = theOptima[i];
				(*allTips[j]).var = theVar[i];
				(*allTips[j]).covar = -1.0;
			}
}


void Phylogeny::inputTaxaXdataCOVAR(std::vector < std::string> theTAXA, std::vector < double> theDATA, std::vector < int > theOptima, std::vector< double > theVar, std::vector< double > theCovar)
{
	std::vector<tree<TIP>::pre_order_iterator> allTips = speciesPointerVector(aPhylogeny);
	for(int i = 0; i < theTAXA.size(); i++) 
		for(int j = 0; j < theTAXA.size(); j++)
			if(theTAXA[i] == (*allTips[j]).species) {
				(*allTips[j]).DATA_x = theDATA[i];
				(*allTips[j]).optima = theOptima[i];
				(*allTips[j]).var = theVar[i];
				(*allTips[j]).covar = theCovar[i];
			}
}
				

void Phylogeny::print_taxaData()
{
	std::vector<tree<TIP>::pre_order_iterator> allTips = speciesPointerVector(aPhylogeny);
	std::cout << "species" << " " << "data" << " " << "optima" << " " << "var" << " " << "covar" <<std::endl;	
	for(int i = 0; i < allTips.size(); i++) 
		std::cout << (*allTips[i]).species << " " << (*allTips[i]).DATA_x << " " << (*allTips[i]).optima << " " << (*allTips[i]).var << " " << (*allTips[i]).covar << std::endl;	
}

Matrix Phylogeny::format_phylo_to_matrix(void)
{
	Phylogeny tempPhylo = aPhylogeny;
	std::vector<tree<TIP>::pre_order_iterator> allTips = speciesPointerVector(aPhylogeny);
	tree<TIP>::pre_order_iterator temp;

	Matrix theMatrix(tempPhylo.stats_numberOfSpecies(), tempPhylo.stats_numberOfSpecies());
	double tempBRL = 0.0; theMatrix = 0.0;

	for(int i = 1; i < theMatrix.Ncols()+1; i++)
		for(int j = 1; j < theMatrix.Nrows()+1; j++) {
			if(i == j) theMatrix(i,j) = tempPhylo.stats_speciesTotalBRL((*allTips[i-1]).species);
			else {
				temp = shareNode(aPhylogeny, (*allTips[i-1]).species, (*allTips[j-1]).species);
				while((*temp).id != 0) {
					tempBRL += (*temp).BRL;
					temp = aPhylogeny.parent(temp);
				}
				theMatrix(i,j) = tempBRL;
				tempBRL = 0.0;
			}
		}

	return theMatrix;
}

ColumnVector Phylogeny::format_phylo_to_DATACOLUMN(void)
{
	Phylogeny tempPhylo = aPhylogeny;
	std::vector<tree<TIP>::pre_order_iterator> allTips = speciesPointerVector(aPhylogeny);

	ColumnVector theDATA(tempPhylo.stats_numberOfSpecies()); 
	theDATA = 0.0; 
	
	for(int i = 1; i < theDATA.Nrows()+1; i++) theDATA(i) = (*allTips[i-1]).DATA_x;

	return theDATA;
}

ColumnVector Phylogeny::format_phylo_to_DATACOLUMN_VAR(void)
{
	Phylogeny tempPhylo = aPhylogeny;
	std::vector<tree<TIP>::pre_order_iterator> allTips = speciesPointerVector(aPhylogeny);

	ColumnVector theDATA(tempPhylo.stats_numberOfSpecies()); 
	theDATA = 0.0; 
	
	for(int i = 1; i < theDATA.Nrows()+1; i++) theDATA(i) = (*allTips[i-1]).var;

	return theDATA;
}

std::vector< std::string > Phylogeny::format_phylo_to_DATACOLUMN_moderators(void)
{
	Phylogeny tempPhylo = aPhylogeny;
	std::vector<tree<TIP>::pre_order_iterator> allTips = speciesPointerVector(aPhylogeny);
	std::vector< std::string > theMods;
	for(int i = 1; i < tempPhylo.stats_numberOfSpecies()+1; i++) theMods.push_back(intToString((*allTips[i-1]).optima));
	return theMods;
}


ColumnVector Phylogeny::format_phylo_to_DATACOLUMN_COVAR(void)
{
	Phylogeny tempPhylo = aPhylogeny;
	std::vector<tree<TIP>::pre_order_iterator> allTips = speciesPointerVector(aPhylogeny);

	ColumnVector theDATA(tempPhylo.stats_numberOfSpecies()); 
	theDATA = 0.0; 
	
	for(int i = 1; i < theDATA.Nrows()+1; i++) theDATA(i) = (*allTips[i-1]).covar;

	return theDATA;
}

/*int getMaxOptima(void)
{
	std::vector<tree<TIP>::pre_order_iterator> allTips = speciesPointerVector(aPhylogeny);
	int MAX = 0;
	for(int i = 0; i < allTips.size(); i++) if((*allTips[i]).optima > MAX) MAX = (*allTips[i]).optima;
	return MAX;
}
*/

void Phylogeny::getVECTORS_BUTLER_KING() {
	tree<TIP>::pre_order_iterator it = aPhylogeny.begin(), end = aPhylogeny.end(), start = aPhylogeny.begin(), ancestor;

	while(it != end) {
		if((*it).species == "node") std::cout << "NA, ";
		else std::cout << '"' << (*it).species << '"' << ", ";
		++it;
    }

	std::cout << std::endl << std::endl;
	it = start;

	while(it != end) {
		std::cout << '"' << ((*it).id+1) << '"' << ", ";
		++it;
    }

	std::cout << std::endl << std::endl;
	it = start;

	while(it != end) {
		if ((*it).id == 0) std::cout << "NA, ";
		//else if((*aPhylogeny.parent(it)).id == 0) std::cout << "NA, "; 
		else std::cout << '"' << ((*aPhylogeny.parent(it)).id + 1) << '"' << ", ";
		++it;
    }


		std::cout << std::endl << std::endl;
	it = start;

	while(it != end) {
		std::cout << '"' << (*it).id << '"' << ", ";
		++it;
    }

	std::cout << std::endl << std::endl;
	it = start;

	while(it != end) {
		if ((*it).species != "node") std::cout << '"' << (*it).DATA_x << '"' << ", ";
		else std::cout << "NA, ";		
		++it;
    }

	std::cout << std::endl << std::endl;
	it = start;

	while(it != end) {
		std::cout << '"' << nodalDistanceToRoot(aPhylogeny, it) << '"' << ", ";
		++it;
    }

		std::cout << std::endl << std::endl;
	it = start;

	while(it != end) {
		std::cout << '"' << (*it).optima << '"' << ", ";
		++it;
    }

	std::cout << std::endl << std::endl;
	it = start;

	while(it != end) {
		if((*it).optima == 0) std::cout << '"' << 1 << '"' << ", ";
		else std::cout << '"' << 0 << '"' << ", ";
		
		++it;
    }

}




#endif


