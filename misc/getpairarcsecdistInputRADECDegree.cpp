/*
 * Sungryong Hong, UT Austin
 *
 *
 *
 *  The code will ONLY look for input files in the running directory!!!
 *
 *	gsl compile option : gcc -Wall -I/usr/local/include -L/usr/local/lib -lgsl -lgslcblas -lm hellogsl.c -o hellogsl.bin
 *
 *	Compile > g++ -Wall -I/usr/local/include -L/usr/local/lib -lgsl -lgslcblas -lm -O3 getpairarcsecdistInputRADECDegree.cpp -o getpairarcsecdistInputRADECDegree.bin
*/

#include <iostream>
#include <iomanip> 
#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

// c library
#include <ctime>  //for random generator
#include <cmath> // for spherical trigonometry
#include <gsl/gsl_rng.h> // gsl random number generator 


using namespace std;

/**************************************************
Collection of Most Used Constants in My calculation
**************************************************/
class MyMathConstant {
	public: 
		MyMathConstant (void); 
		double pi,twopi,halfpi;
};

MyMathConstant::MyMathConstant (void){
	pi=3.141592653589793238462643383279502884;
	twopi = 2.0 * pi;
	halfpi = 0.5 * pi;
}


/**************************************************
Collection of Spherical Trigonometry Tools
**************************************************/
class SphericalTrigonometry {    // all quantities in radian scale
	public: 
		SphericalTrigonometry (void);
		double dumpra, dumpdec,positionErrorCap;
		

		double sphdist(double sara, double sadec, double sbra, double sbdec); //return sphdist in radian scale
		void dumpWhereWhenPlusThis(double inra, double indec, double dirangle, double dirlength); // at inra,indec.. go to dir vector 
		void printRaDecCircle(double inra, double indec, double inradius, int numpolygon, string ofile); // print ra-dec circle at (inra,indec)

		double hourtorad,degtorad,radtoarcsec,radtohour,radtodeg; // ra hour to radian; dec deg to radian
		double halfpi,pi,twopi;
};


SphericalTrigonometry::SphericalTrigonometry (void) {    // null constructor
	hourtorad = 3.141592653589793238462643383279502884 * 2.0 / 24.0; 
	degtorad =  3.141592653589793238462643383279502884 * 2.0 / 360.0;
	radtoarcsec = 360.0 * 3600.0 / (3.141592653589793238462643383279502884 * 2.0); 
	radtohour = 24.0 / (3.141592653589793238462643383279502884 * 2.0); 
	radtodeg = 360.0 / (3.141592653589793238462643383279502884 * 2.0); 
	pi = 3.141592653589793238462643383279502884;
	halfpi = pi /2.0;
	twopi = pi * 2.0;
	dumpra = -1.0; // buffer;; temporary storage for various methods within this class
	dumpdec = -1.0;
	positionErrorCap =  0.1 / 3600.0 * degtorad; // 0.2 arcsec is a tolerance error for position coordinate
}

double SphericalTrigonometry::sphdist (double sara,double sadec, double sbra, double sbdec){  //return spherical distance
	double dtmp;
	double ax,ay,az,bx,by,bz;

	dtmp = halfpi - sadec;
	ax = sin(dtmp) * cos(sara);
	ay = sin(dtmp) * sin(sara);
	az = cos(dtmp);
	
	dtmp = halfpi - sbdec;
	bx = sin(dtmp) * cos(sbra);
	by = sin(dtmp) * sin(sbra);
	bz = cos(dtmp);

	dtmp = acos(ax*bx + ay*by + az*bz);

	return dtmp;
}

// severe problem when the position is near dec = +-90  
void SphericalTrigonometry::dumpWhereWhenPlusThis(double inra, double indec, double dirangle, double dirlength) { // at inra,indec.. go to dir vector .. 
	double optdec, optra, tmpd, tmpcosdec, tmpra, tmpdec, zerora,
		zerocosdec,steplength, tmplength, stepdirnow,stepdirprevious ; // cosdec = cos (pi/2 - dec) ; all radian scale
	int iwalk; // number of root finding walks 
	//int interrupt;

	zerora = inra;
	zerocosdec = cos(halfpi - indec);

	// initalize the root finding iteration
	steplength = dirlength/2.0; //initial steplength
	stepdirnow = 1.0; // positive direction for the first step
	stepdirprevious = 1.0; //
	tmplength = dirlength; // initial position length
	tmpra = tmplength * cos(dirangle);
	tmpcosdec = tmplength * sin(dirangle);
	tmpdec = halfpi - acos(tmpcosdec); // dec radian


	// root finding start ;; step back and forth til the step size < positionError
	iwalk = 0; //counting the number of steps
	while (steplength > positionErrorCap) {
			
		tmpd = sphdist(inra, indec,inra+tmpra,indec+tmpdec); // get the sphdist for the first step
		if ( tmpd < dirlength) {stepdirnow = 1.0;} 
		else {stepdirnow = -1.0;}


//		cout << "SphericalTrigonometry::dumpWhereWhenPlusThis : inra indec iwalk tmpra(arcsec) tmpdec(arcsec) sphdist inlength dirPrev dirNow "
//			<< inra <<" "<< indec <<" "<< iwalk <<" "<< tmpra*radtoarcsec <<" "<< tmpdec*radtoarcsec <<" "<< tmpd*radtoarcsec <<" "<< dirlength*radtoarcsec 
//			<<" "<<stepdirprevious <<" "<< stepdirnow << endl; 
//		cout << "steplength before " << steplength; 
		if ( stepdirnow*stepdirprevious < 0.0) { // opposite direction 		
			steplength *= 0.5; //half-size the step
		}
//		cout << " steplength after " << steplength << endl; 
		
		// assign the next walk 
//		cout << "tmplength before " << tmplength; 
		tmplength = tmplength + stepdirnow * steplength;
//		cout << " tmplength after " << tmplength << endl; 

		stepdirprevious = stepdirnow; 
		tmpra = tmplength * cos(dirangle);
		tmpcosdec = tmplength * sin(dirangle);
		tmpdec = halfpi - acos(tmpcosdec); // dec radian
		iwalk++;

		// interrupt for debug
		//cin >> interrupt;

		// print out if not converging! 
		if ( iwalk > 1000 ){ cout << "SphericalTrigonometry::dumpWhereWhenPlusThis : Too many steps to find roots : iwalk = " << iwalk << endl;}

	}

	optdec = tmpdec;
	optra = tmpra;

	// write the results to dumps
	dumpra = inra + optra;
	dumpdec = indec + optdec; 
}


void SphericalTrigonometry::printRaDecCircle(double inra, double indec, double inradius, int numpolygon, string ofile){  // using gsl random num generator 
	double rastart, decstart,tmpincr,tmpdirvalue;
	ofstream of;
  	of.open (ofile);

	tmpdirvalue = 0.0; 
	tmpincr = (double) numpolygon;
	tmpincr = twopi/tmpincr;

	dumpWhereWhenPlusThis(inra, indec,tmpdirvalue,inradius);
	rastart = dumpra;
	decstart = dumpdec;
  	of << rastart << " " << decstart << endl;
	for (int i=0; i < numpolygon; i++){
		tmpdirvalue += tmpincr;
		dumpWhereWhenPlusThis(inra, indec,tmpdirvalue,inradius);
  		of << dumpra << " " << dumpdec << endl;
	}
  	of.close();
}



/*************************************************
SkyPolygon in All Shapes 
and 
Tools to (re)shape itself. 
*************************************************/
class SkyPolygon {  // ra dec array forming a polyong on sky
	public:
		int size;
		int sizeplusone; // n_gons need 1 more to match start vertex = end vertex ;; making loop polygon
		std::vector<double> ragon, decgon; //all in radian scale

		SkyPolygon (int size);		
//		assignCircle (double cenra, double cendec, double cradius); // make a circle with the size of "sizeplusone"
};

SkyPolygon::SkyPolygon (int s) {
	size = s;
	sizeplusone = s+1;
	ragon.resize(sizeplusone);
	decgon.resize(sizeplusone);
}	



/*************************************************
GalaxyPopulation classs
*************************************************/
class GalaxyPopulation {  // ra dec array forming a polyong on sky
	public:
		int size;
		SphericalTrigonometry stool;
		std::vector<double> ra, dec; //all in radian scale

		GalaxyPopulation (void);		
		GalaxyPopulation (int size);
		void addGalaxy(double ra, double dec);
		double getra(int idx);
		double getdec(int idx);
		void printPairWiseAngularDistance(string ofilename);
};
//void contructor
GalaxyPopulation::GalaxyPopulation (void) {
	size = 0;
}

GalaxyPopulation::GalaxyPopulation (int s) {
	size = s;
	ra.resize(size);
	dec.resize(size);
}
void GalaxyPopulation::addGalaxy (double inra, double indec) {
	size += 1;ra.push_back(inra);dec.push_back(indec);
}
double GalaxyPopulation::getra (int i) {
	if (i < size) { return ra[i]; } else { return -1.0;} //return negative ra, when out-of-index range
}
double GalaxyPopulation::getdec (int i) {
	if (i < size) { return dec[i]; } else { return -1.0;} //return negative dec, when out-of-index range
}
void GalaxyPopulation::printPairWiseAngularDistance(string ofile){  // using gsl random num generator 
	double dij;
	int ix,iy ; //index for i and j .. pair dist is dij

	ofstream of;
  	of.open (ofile);

	dij = 0.0; 

	for (ix=0; ix < size; ix++){
		for (iy=(ix+1); iy < size; iy++){
			dij = stool.sphdist(ra[ix] *stool.degtorad , dec[ix]*stool.degtorad , ra[iy]*stool.degtorad , dec[iy]*stool.degtorad );
			dij *= stool.radtoarcsec;
  			//cout << ix << " " << iy << " " << dij <<  endl; // debug : display the pair indices 
			of << fixed << setprecision(5); 
  			of << ix << " " << iy << " " << dij << endl;
		}
	}
  	of.close();
}


int main(int argc, char* argv[])
{
	int i;
	string line,rastr,decstr;
	ifstream infile;
	GalaxyPopulation g;

	// parsing main arguments 
	cout << "argc = " << argc << endl; 
	for(int i = 0; i < argc; i++) 
		cout << "argv[" << i << "] = " << argv[i] << endl; 
	if (argc!=3){
		cout << "Usage : <command> <in_text_file; ra dec in degree> <out_text_file; i j dij(arcsec)> "<< endl;
		exit(EXIT_FAILURE);
	}
	infile.open(argv[1]);
	if (infile.is_open()){
		while (getline(infile,line)){
			cout << line << endl;
			stringstream ss(line);
			ss >> rastr >> decstr;
			//cout << stod(rastr) << endl;
			//cout << stod(decstr) << endl;
			g.addGalaxy(stod(rastr),stod(decstr));
		}
	}
	else{
		cout << "Unable to open input-file.";
		exit(EXIT_FAILURE);
	}
//	cout << fixed << setprecision(5); 
//	for (i=0; i < g.size+3; ++i){
//		cout << i << " " << g.getra(i) << " " << g.getdec(i) << endl;
//	}
	
	g.printPairWiseAngularDistance(argv[2]);


	// stool instance
	//SphericalTrigonometry stool;
	//==test log
	//stool.printRaDecCircle(3.4,0.9,3.0*stool.degtorad,100,"tmpcircle.txt");

	infile.close();

    return EXIT_SUCCESS;
}
