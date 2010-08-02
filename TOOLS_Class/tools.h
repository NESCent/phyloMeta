#ifndef TOOLS_H
#define TOOLS_H

#include <sys/stat.h>
#include <conio.h>


void exit() 
{
	std::cout << std::endl << "Press any key to continue...";
     cin.ignore(0,'\n');
     getch();
}

bool isMatrixPositiveDefinate(const Matrix & aMatrix)
{
	bool posDef = true;
	Matrix temp; temp = aMatrix;
	while(!temp.IsZero()) {
		if(temp.Determinant() < 0) posDef = false;
		temp = temp.SubMatrix(1, temp.Nrows()-1, 1, temp.Ncols() -1);
	}
	return posDef;
}

Matrix MatrixSquareRoot(const Matrix & temp)
{
	SymmetricMatrix S(temp.Ncols()); 
	for(int i = 1; i < temp.Ncols()+1; i++) 
		for(int j = 1; j < temp.Ncols()+1; j++) S(i,j) = temp(i,j);
	LowerTriangularMatrix L = Cholesky(S);
	DiagonalMatrix D;
	Matrix V, U, result;
	SVD(L, D, U, V);
	return U * D * U.t();
}

bool isMatrixSymmetric(const Matrix & aMatrix)
{
	return aMatrix == aMatrix.t();
}

// DATA CONVERSIONS

unsigned long factorial(unsigned long aNum)
{
	unsigned long temp = aNum;
	for(int i = aNum-1; i > 1; i--) temp *= i;
	return temp;
}

char* stringToChar(std::string theString)
{
	char theChar[20];
	strcpy(theChar, theString.c_str());
	return theChar;
}

std::string floatToString(double aDouble)
{
	std::ostringstream ostr ;
	ostr << aDouble;
	std::string s ( ostr.str () ) ;
	return s;
}

std::string floatToString(double aDouble, int precision)
{
	std::ostringstream ostr ;
	ostr << std::setprecision ( precision ) << aDouble ;
	std::string s ( ostr.str () ) ;
	return s;
}


int stringToInt(std::string aString)
{
	std::istringstream istr(aString); int anInt; istr >> anInt;
	return anInt;
}

double stringToDouble(std::string aString)
{
	std::istringstream istr(aString); double anInt; istr >> anInt;
	return anInt;
}

std::string intToString(int anInt)
{
	std::stringstream sstr; sstr << anInt; std::string aString = sstr.str();
	return aString;
}


// ERROR FUNCTIONS

void nrerror(char error_text[])
{
	std::cerr << "Run-time error..." << std::endl
	<< error_text << std::endl
	<< "...now exiting to system..." << std::endl;
	exit(1);
}

void errorACTION(std::string error_text, std::string error_action, bool exit_system)
{
	const char *error = error_text.c_str(), *action = error_action.c_str();
	std::cerr << "ERROR: " << error << " ACTION: " << action << std::endl;
	if(exit_system) exit(1);
}

// FILE ACCESS FUNCTIONS

int numberOfRowsInFile(const char* FileName)
{
	char c, d; int lines = 0;
	ifstream theFile(FileName);
	while(theFile.get(c)) {
		if(c == '\n') lines++;
		d = c;
	}
	if(d != '\n') lines++; //cout << d << endl;
	theFile.close();
	return (lines);
}

int numberOfColumnsInFile(const char* FileName)
{
	char c; int lines = 0, columns = 0;
	ifstream theFile(FileName);
	while(theFile.get(c)) {
		if((c == '\t') || (c == ' ')) columns++;
		if(c == '\n') break;
	}
	theFile.close();
	return (columns+1);
}

bool FileExist(const char* FileName)
{
	struct stat my_stat;
	return (stat(FileName, &my_stat) == 0);
}

bool IsDirectory(const char* FileName)
{
	struct stat my_stat;
	if (stat(FileName, &my_stat) != 0) return false;
	return ((my_stat.st_mode & S_IFDIR) != 0);
}

std::string getFileName(std::string aFileName)
{
	int size = aFileName.size(), step = size - 1;
	std::string aNum, aChar, file, temp;
	temp = aFileName[aFileName.size() - 1];
	if(temp != ")") return (aFileName + "(1)");

	while(aChar != "(") {
		aChar = aFileName[step];
		if(aChar != "(") step--;
		if(step == 1) return aFileName;
	}

	while(aChar != ")") {
		aChar = aFileName[step];
		if((aChar != "(") && (aChar != ")")) aNum += aChar;
		step++;
		if(step > aFileName.size()) return aFileName;
	}

	aChar = ""; step = 0;
	while(aChar != "(") {
		aChar = aFileName[step];
		if(aChar != "(") file += aChar;
		step++;
	}

	int i = stringToInt(aNum); i++;

	return file + "(" + intToString(i) + ")";
}

std::string saveFile(const char *aFileName, std::string aStringToSave, bool append)
{
	std::string fileName = aFileName, newName = aFileName;
	const char *aString = aStringToSave.c_str(), *aName = fileName.c_str();

	if(!append) {
		while(FileExist(aName)) {
			newName = getFileName(newName);
			aName = newName.c_str();
		}
		if(newName != fileName) errorACTION(fileName + " already exists", "Saved file to " + newName, false);
	}

	ofstream aFile(aName, ((append == true) ? (ios::out | ios::app) : ios::out));
	if(!aFile) errorACTION("Could not open File", "Exiting system", true);
	aFile << aString << endl;
	aFile.close();
	return aName;
}

std::string saveFile_STRING(std::string aFileName, std::string aStringToSave, bool append)
{
	std::string fileName = aFileName, newName = aFileName;
	const char *aString = aStringToSave.c_str(), *aName = fileName.c_str();

	if(!append) {
		while(FileExist(aName)) {
			newName = getFileName(newName);
			aName = newName.c_str();
		}
		if(newName != fileName) errorACTION(fileName + " already exists", "Saved file to " + newName, false);
	}

	ofstream aFile(aName, ((append == true) ? (ios::out | ios::app) : ios::out));
	if(!aFile) errorACTION("Could not open File", "Exiting system", true);
	aFile << aString << endl;
	aFile.close();
	return aName;
}

std::string saveFileOK(char aFileName[])
{
	std::string fileName = aFileName, newName = aFileName;
	const char *aName = fileName.c_str();

	while(FileExist(aName)) {
		newName = getFileName(newName);
		aName = newName.c_str();
	}
	if(newName != fileName) errorACTION(fileName + " already exists", "Saved file to " + newName, false);
	
	return aName;
}


//printing functions

void printStringVector(std::vector< std::string > aStringVector)
{
	for(int i = 0; i < aStringVector.size(); i++) std::cout << aStringVector[i] << " ";
	std::cout << std::endl;
}

void printStringVector(std::vector< std::string > aStringVector, int spaces)
{
	for(int i = 0; i < aStringVector.size(); i++) std::cout << std::setw(spaces) << aStringVector[i];
	std::cout << std::endl;
}


void printDoubleVector(std::vector< double > aStringVector)
{
	int count = 0;
	for(int i = 0; i < aStringVector.size(); i++) {
		std::cout << std::setprecision(6) << aStringVector[i] << " ";
		if(count == 10) {std::cout << std::endl; count = -1;}
		count++;
	}
	std::cout << std::endl;
}

void printIntVector(std::vector< int > aStringVector)
{
	for(int i = 0; i < aStringVector.size(); i++) std::cout << aStringVector[i] << " ";
	std::cout << std::endl;
}

double absolute(double Nbr)
{
	return (Nbr >= 0) ? Nbr : -Nbr;
}

int isEven(int number)
{
	if(number == 2 || number == 0) return true;
	else if (number == 1) return false;
	isEven(number-2);
}

bool allCONTENTSsame(std::vector < std::string > aStringVector, std::string aSameString)
{
	if(aStringVector[0] != aSameString) return false;
	for(int i = 0; i < aStringVector.size(); i++)
		if(aStringVector[0] != aStringVector[i]) return false;
	return true;
}

bool isInVectorINT(int aNum, std::vector<int> theList)
{
	for(int i = 0; i < theList.size(); i++)
		if(aNum == theList[i]) return true;
	return false;
}

bool isInVectorSTRING(std::string aString, std::vector< std::string > theList)
{
	for(int i = 0; i < theList.size(); i++)
		if(aString == theList[i]) return true;
	return false;
}

bool isInVectorDOUBLE(int aNum, std::vector<int> theList)
{
	for(int i = 0; i < theList.size(); i++)
		if(aNum == theList[i]) return true;
	return false;
}

void myBanner (std::string Name, std::string Version, std::string Date)
{
	int sizeName = Name.size(), sizeVersion = Version.size(), first = 0, last = 0;
	bool even = isEven(sizeName);
	sizeName = (even ? sizeName/2 : (sizeName/2 + 1));
	first = 24 - sizeName;
	last = 24 + sizeName;
	std::string spacesF, spacesL, theName;

	for(int i = 0; i < first; i++) spacesF += " ";
	for(i = last; i <= (even ? 47 : 48); i++) spacesL += " ";

	theName = spacesF + Name + spacesL;

	std::cout << "               "<< (char)218 << (char)196<< (char)196<<(char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196 << (char)191 << std::endl;
	std::cout << "               "<< (char)179 << "                                                " << (char)179 << std::endl;
	std::cout << "               "<< (char)179 << theName << (char)179 << std::endl;
	std::cout << "               "<< (char)179 << "            version " << Version << " (" << Date << ")            " << (char)179 << std::endl;
	std::cout << "               "<< (char)179 << "                                                " << (char)179 << std::endl;
	std::cout << "               "<< (char)179 << "                " << (char)196<<(char)196<<(char)196<< (char)196<<(char)196<<(char)196<< (char)196<<(char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< "                " << (char)179 << std::endl;
	std::cout << "               "<< (char)179 << "                                                " << (char)179 << std::endl;
	std::cout << "               "<< (char)179 << "               Marc J. Lajeunesse               " << (char)179 << std::endl;
	std::cout << "               "<< (char)179 << "           marc.lajeunesse@NESCent.org          " << (char)179 << std::endl;
	std::cout << "               "<< (char)179 << "                                                " << (char)179 << std::endl;
	std::cout << "               "<< (char)179 << "   The National Evolutionary Synthesis Center   " << (char)179 << std::endl;
	std::cout << "               "<< (char)179 << "         2024 W. Main Street, Suite A200        " << (char)179 << std::endl;
	std::cout << "               "<< (char)179 << "             Durham, NC 27705-4667              " << (char)179 << std::endl;
	std::cout << "               "<< (char)179 << "                                                " << (char)179 << std::endl;
	std::cout << "               "<< (char)192 << (char)196<< (char)196<<(char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196<< (char)196 << (char)217 << std::endl << std::endl << std::endl;

}


//return vector without duplicate strings
std::vector < std::string > cleanStringVector(std::vector < std::string > theData) 
{
	std::vector < std::string > theSpecies;
	bool found = false;
	for(int i = 0; i < theData.size(); i++) { 
		if(theSpecies.size() == 0 ) theSpecies.push_back(theData[i]);
		else {
			for(int j = 0; j < theSpecies.size(); j++) if(theData[i] == theSpecies[j]) found = true;
			if(found == false) theSpecies.push_back(theData[i]);
			found = false;
		}
	}
	return theSpecies;
}

std::vector < int > cleanIntVector(std::vector < int > theData) 
{
	std::vector < int > theSpecies;
	bool found = false;
	for(int i = 0; i < theData.size(); i++) { 
		if(theSpecies.size() == 0 ) theSpecies.push_back(theData[i+1]);
		else {
			for(int j = 0; j < theSpecies.size(); j++) if(theData[i] == theSpecies[j]) found = true;
			if(found == false) theSpecies.push_back(theData[i]);
			found = false;
		}
	}
	return theSpecies;
}



std::vector<std::string> deleteVectorRow(std::vector<std::string> theVector, int location)
{
	if(location > theVector.size()) return theVector;

	std::vector<std::string> temp;
	for (int i = 0; i < theVector.size(); i++)
		if(location != i) temp.push_back(theVector[i]);
	return temp;
}

/*void matrixPRINT(Matrix theMatrix, char aString[] )
{
int middle = theMatrix.Nrows()/2;
std::string theString = aString;

for(int i = 1; i <= theMatrix.Nrows(); i++) {
if(i == middle) cout << setprecision(3) << aString << " = " << theMatrix.Row(i);
else {
for(int j = 0; j < (theString.size() + 3); j++) cout << " ";
cout << setprecision(3) << theMatrix.Row(i);
}
}
cout << endl;

}*/

void printProgress(double value, double MAX) 
{
		cout << " " << setprecision(3 ) << (value/MAX)*100 << "%        \r";
		cout.flush();
}

Matrix matrixFromFile(const char* FileName, bool printMatrix)
{
	int theRows = numberOfRowsInFile(FileName), theCols = numberOfColumnsInFile(FileName), i = 1, j = 1;
	double theData = 0.0;
	Matrix theMatrix(theRows, theCols); theMatrix = 0.0;

	ifstream theFile(FileName);
	while(theFile >> theData) {
		theMatrix(i,j) = theData;
		if(j == theCols) {
			j = 1; 
			i++;
		}	
		else j++;
	}
	theFile.close();

	if(printMatrix) cout << theMatrix << endl;
	
	return theMatrix;
}


#endif