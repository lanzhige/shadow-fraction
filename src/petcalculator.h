#ifndef CALCULATEMRT_PETCALCULATOR_H_
#define CALCULATEMRT_PETCALCULATOR_H_

#include <sstream> 

#include "jsonstruct.h"

namespace calculate {
	class PETCalculator {
	public:
		PETCalculator(double airTemperature, double mrt, double vaporPressure, double windVelocity, const file::JsonStruct *json);
		void calculate(std::stringstream *ss);
	private:
		double acl = 0;
		double adu;
		double aeff;
		double age;
		double c[11];
		double cair;
		double cb;
		double cbare;
		double cclo;
		double csum;
		double di;
		double ed;
		double emcl;
		double emsk;
		double enbal;
		double enbal2;
		double ere;
		double erel;
		double eres;
		double esw;
		double eswdif;
		double eswphy;
		double eswpot;
		double eta;
		double evap;
		double facl;
		double fcl;
		double fec;
		double feff;
		double food;
		double h;
		double hc;
		double he;
		double ht;
		double htcl;
		double icl;
		double mbody;
		double met;
		double metbf;
		double metbm;
		double metf;
		double metm;
		double p;
		double po;
		double r1;
		double r2;
		double rbare;
		double rcl;
		double rclo;
		double rclo2;
		double rdcl;
		double rdsk;
		double rob;
		double rsum;
		double rtv;
		double sigm;
		double sw;
		double swf;
		double swm;
		double ta;
		double tbody;
		double tcl;
		double tcore[8];
		double tex;
		double tmrt;
		double tsk;
		double tx;
		double v;
		double vb;
		double vb1;
		double vb2;
		double vpa;
		double vpex;
		double vpts;
		double wetsk;
		double wd;
		double work;
		double wr;
		double ws;
		double wsum;
		double xx;

		double contr;
		int count1;
		int count2;
		int count3;
		int j;
		double pos;
		int sex;
		double esc;

		double y;
		double rsum_temp;
		double csum_temp;
		double ere_temp;

		std::stringstream *ss;

		void label_90();
		void label_20();
		void label_22();
		void label_24();
		void label_30();
		void label_40();
		void label_50();
		void label_60();
		void label_70();
		void label_80();

		void label_150();
		void label_160();
		void label_170();

		void subroutine_PET();
	};
}

#endif // !CALCULATEMRT_PETCALCULATOR_H_
