#include "petcalculator.h"

namespace calculate {
	/**
		@brief Initialize the PETCalculator, assign input variables
		@param airTemperature The current air temperature
		@param mrt The Mean Radiant Temperature
		@param vaporPressure The current vapor pressure
		@param windVelocity The current wind velocity
		@param json The configuration settings
	*/
	PETCalculator::PETCalculator(double airTemperature, double mrt, double vaporPressure, double windVelocity, const file::JsonStruct *json) {
		ta = airTemperature;
		tmrt = mrt;
		vpa = vaporPressure;
		v = windVelocity;
		
		po = 1013.25;
		p = 1013.25;
		rob = 1.06;
		cb = 3.64 * 1000.0;
		food = 0.0;
		emsk = 0.99;
		emcl = 0.95;
		evap = 2.42 * pow(10, 6);
		sigm = 5.67 * pow(10, -8);

		age = json->age;
		mbody = json->weight;
		ht = json->height;
		work = json->activity;
		eta = 0.0;
		icl = json->clothing;
		fcl = 1.15;
		pos = json->standing;
		sex = json->male;
	}
	
	/**
		@brief Calculate PET results
		@param ss The stringstream for output buffer
	*/
	void PETCalculator::calculate(std::stringstream *ss) {
		this->ss = ss;
		metbf = 3.19 * pow(mbody, 0.75) * (1. + 0.004 *(30. - age) + 0.018 *((ht*100. / (pow(mbody, (1.0 / 3.0)))) - 42.1));
		metbm = 3.45 * pow(mbody, 0.75) * (1. + 0.004 *(30. - age) + 0.010 *((ht*100. / (pow(mbody, (1.0 / 3.0)))) - 43.4));
		metm = work + metbm;
		metf = work + metbf;

		if (sex == 1)	met = metm;
		if (sex == 2)	met = metf;
		h = met * (1.0 - eta);

		cair = 1.01 * 1000.;
		tex = 0.47 * ta + 21.0;
		rtv = 1.44 * pow(10, -6) * met;
		eres = cair * (ta - tex) * rtv;

		vpex = 6.11 * pow(10, (7.45 * tex / (235. + tex)));
		erel = 0.623 * evap / p * (vpa - vpex) * rtv;

		ere = eres + erel;

		wetsk = 0.0;
		adu = 0.203 * pow(mbody, 0.425) * pow(ht, 0.725);

		hc = 2.67 + (6.5 * pow(v, 0.67));
		hc = hc * pow((p / po), 0.55);
		feff = 0.725;
		facl = (-2.36 + 173.51 * icl - 100.76 * icl * icl + 19.28 * (pow(icl, 3.0))) / 100.0;

		if (facl > 1.0)   facl = 1.0;
		rcl = (icl / 6.45) / facl;
		
		if (icl>2.0) y = 1.0;
		if ((icl > 0.6) && (icl < 2.0))  y = (ht - 0.2) / ht;
		if ((icl <= 0.6) && (icl < 0.3)) y = 0.5;
		if ((icl <= 0.3) && (icl > 0.))  y = 0.1;

		r2 = adu * (fcl - 1. + facl) / (2. * 3.14 * ht * y);
		r1 = facl * adu / (2. * 3.14 * ht * y);

		di = r2 - r1;

		j = 0;
		label_90();

		ws = sw * 3600.0 * 1000.0;
		if (ws > 2000.0) ws = 2000.0;
		wd = ed / evap * 3600. * (-1000.);
		wr = erel / evap * 3600. * (-1000.);

		wsum = ws + wr + wd;

		return subroutine_PET();
	}

	void PETCalculator::label_90() {
		j++;
		tsk = 34.0;
		count1 = 0;
		tcl = (ta + tmrt + tsk) / 3.0;
		count3 = 1;
		enbal2 = 0.0;

		label_20();
	}

	void PETCalculator::label_20() {
		acl = adu * facl + adu * (fcl - 1.0);
		rclo2 = emcl * sigm * (pow((tcl + 273.2), 4.0) - pow((tmrt + 273.2), 4.0)) * feff;
		htcl = 6.28 * ht * y * di / (rcl * log(r2 / r1) * acl);
		tsk = 1. / htcl * (hc * (tcl - ta) + rclo2) + tcl;

		//      STRAHLUNGSSALDO
		aeff = adu * feff;
		rbare = aeff * (1. - facl) * emsk * sigm * (pow((tmrt + 273.2), 4.0) - pow((tsk + 273.2), 4.0));
		rclo = feff * acl * emcl * sigm * (pow((tmrt + 273.2), 4.0) - pow((tcl + 273.2), 4.0));
		rsum = rbare + rclo;

		//       KONVEKTION
		cbare = hc * (ta - tsk) * adu * (1. - facl);
		cclo = hc * (ta - tcl) * acl;
		csum = cbare + cclo;

		//       KERNTEMPERATUR
		c[0] = h + ere;
		c[1] = adu * rob * cb;
		c[2] = 18. - 0.5 * tsk;
		c[3] = 5.28 * adu * c[2];
		c[4] = 0.0208 * c[1];
		c[5] = 0.76075 * c[1];
		c[6] = c[3] - c[5] - tsk * c[4];
		c[7] = -c[0] * c[2] - tsk * c[3] + tsk * c[5];
		c[8] = c[6] * c[6] - 4.0 * c[4] * c[7];
		c[9] = 5.28 * adu - c[5] - c[4] * tsk;
		c[10] = c[9] * c[9] - 4.0 * c[4] * (c[5] * tsk - c[0] - 5.28 * adu * tsk);

		if (tsk == 36.0) tsk = 36.01;
		tcore[7] = c[0] / (5.28 * adu + c[1] * 6.3 / 3600.) + tsk;
		tcore[3] = c[0] / (5.28 * adu + (c[1] * 6.3 / 3600.) / (1 + 0.5 * (34. - tsk))) + tsk;


		if (c[10] < 0.0){
			label_22();
		} else{
			tcore[6] = (-c[9] - pow(c[10], 0.5)) / (2. * c[4]);
			tcore[1] = (-c[9] + pow(c[10], 0.5)) / (2. * c[4]);
			label_22();
		}

		tbody = 0.1 * tsk + 0.9 * tcore[j];
		swm = 304.94 * (tbody - 36.6) * adu / 3600000.0;
		vpts = 6.11 * pow(10., (7.45 * tsk / (235. + tsk)));

		if (tbody <= 36.6) swm = 0.0;
		swf = 0.7 * swm;

		if (sex == 1) sw = swm;
		if (sex == 2) sw = swf;
		eswphy = -sw * evap;
		he = 0.633 * hc / (p * cair);
		fec = 1. / (1. + 0.92 * hc * rcl);
		eswpot = he * (vpa - vpts) * adu * evap * fec;
		wetsk = eswphy / eswpot;
		if (wetsk > 1.0) wetsk = 1.0;

		eswdif = eswphy - eswpot;

		if (eswdif <= 0.0) esw = eswpot;
		if (eswdif > 0.0) esw = eswphy;
		if (esw > 0.0) esw = 0.0;

		rdsk = 0.79 * pow(10.0, 7.0);
		rdcl = 0.0;
		ed = evap / (rdsk + rdcl) * adu * (1 - wetsk) * (vpa - vpts);

		vb1 = 34. - tsk;
		vb2 = tcore[j] - 36.6;

		if (vb2 < 0.0) vb2 = 0.0;
		if (vb1 < 0.0) vb1 = 0.0;
		vb = (6.3 + 75. * (vb2)) / (1. + 0.5 * vb1);

		enbal = h + ed + ere + esw + csum + rsum + food;

		if (count1 == 0) xx = 1.0;
		if (count1 == 1) xx = 0.1;
		if (count1 == 2) xx = 0.01;
		if (count1 == 3) xx = 0.001;

		if (enbal > 0.0) tcl = tcl + xx;
		if (enbal < 0.0) tcl = tcl - xx;
		if ((enbal <= 0.0) && (enbal2 > 0.0)) {
			label_30();
		}
		else if ((enbal >= 0.0) && (enbal2 < 0.0)) {
			label_30();
		} else {
			enbal2 = enbal;
			count3 = count3 + 1;
			if (count3 > 200) {
				label_30();
			} else{
				label_20();
			}
		}

	}

	void PETCalculator::label_22() {
		if (c[8] < 0.0) {
			label_24();
		}else{
			tcore[2] = (-c[6] + pow(abs(c[8]), 0.5)) / (2. * c[4]);
			tcore[5] = (-c[6] - pow(abs(c[8]), 0.5)) / (2. * c[4]);
			label_24();
		}
	}

	void PETCalculator::label_24() {
		tcore[4] = c[0] / (5.28 * adu + c[1] * 1.0 / 40.0) + tsk;
	}

	void PETCalculator::label_30() {
		if ((count1 == 0) || (count1 == 1) || (count1 == 2)) {
			count1 = count1 + 1;
			enbal2 = 0.0;
			label_20();
		}
		else if (count1 == 3) {
			if ((j == 2) || (j == 5)) {
				label_40();
			}
			else if ((j == 6) || (j == 1)) {
				label_50();
			}
			else if (j == 3) {
				label_60();
			}
			else if (j == 7) {
				label_70();
			}
			else if (j == 4) {
				label_80();
			}
		}
		else {
			label_40();
		}
	}

	void PETCalculator::label_40() {
		if (c[8] < 0.0) {
			label_90();
		}
		else if ((tcore[j] >= 36.6) && (tsk <= 34.050)) {
			label_80();
		}else{
			label_90();
		}
	}
	void PETCalculator::label_50() {
		if (c[10] < 0.0) {
			label_90();
		}
		else if ((tcore[j] >= 36.6) && (tsk > 33.850)) {
			label_80();
		}else{
			label_90();
		}
	}
	void PETCalculator::label_60() {
		if ((tcore[j] < 36.6) && (tsk <= 34.000)) {
			label_80();
		}else{
			label_90();
		}
	}
	void PETCalculator::label_70() {
		if ((tcore[j] < 36.6) && (tsk > 34.000)) {
			label_80();
		}else{
			label_90();
		}
	}
	void PETCalculator::label_80() {
		if ((j != 4) && (vb >= 91.0)) {
			label_90();
		}
		else if
			((j == 4) && (vb < 89.0)) {
			label_90();
		}else if(vb > 90.0){
			vb = 90.0;
		}
	}

	void PETCalculator::label_150() {
		hc = 2.67 + 6.5 * pow(0.1, 0.67);
		hc = hc * pow((p / po), 0.55);

		aeff = adu * feff;
		rbare = aeff * (1. - facl) * emsk * sigm *	(pow((tx + 273.2), 4.0) - pow((tsk + 273.2), 4.0));
		rclo = feff * acl * emcl * sigm * (pow((tx + 273.2), 4.0) - pow((tcl + 273.2), 4.0));
		rsum_temp = rbare + rclo;

		cbare = hc * (tx - tsk) * adu * (1.0 - facl);
		cclo = hc * (tx - tcl) * acl;
		csum_temp = cbare + cclo;

		ed = evap / (rdsk + rdcl) * adu * (1. - wetsk) * (12. - vpts);

		tex = 0.47 * tx + 21.0;
		eres = cair * (tx - tex) * rtv;
		vpex = 6.11 * pow(10.0, (7.45 * tex / (235. + tex)));
		erel = 0.623 * evap / p * (12. - vpex) * rtv;
		ere_temp = eres + erel;

		enbal = h + ed + ere_temp + esw + csum_temp + rsum_temp;

		if (count1 == 0)  xx = 1.0;
		if (count1 == 1)  xx = 0.1;
		if (count1 == 2)  xx = 0.01;
		if (count1 == 3)  xx = 0.001;
		if (enbal > 0.0)   tx = tx - xx;
		if (enbal < 0.0)   tx = tx + xx;
		if ((enbal <= 0.0) && (enbal2 > 0.0)) {
			return label_160();
		}
		else if ((enbal >= 0.0) && (enbal2 < 0.0)) {
			return label_160();
		}else{
			enbal2 = enbal;
			return label_150();
		}
	}

	void PETCalculator::label_160() {
		count1 = count1 + 1;
		if (count1 == 4) {
			return label_170();
		}else{
			return label_150();
		}
	}

	void PETCalculator::label_170() {
		double airTemperature = ta;
		
		double radiationTemperature = tmrt;

		double steamPressure = vpa;

		double windSpeed = v;

		double coreTemperature = tcore[j];

		double skinTemperature = tsk;

		double totalWaterLoss = wsum;

		double skinWetting = wetsk;

		double internalHeat = h;

		double radiationBalance = rsum;
		
		double convection = csum;
		
		double waterVaporDiffusion = ed;
		
		double weldingEvaporation = esw;
		
		double Respiration = ere;
		
		double pet = tx;

		(*ss) << pet << "," << coreTemperature << "," << skinTemperature << "," << totalWaterLoss << "," << skinWetting << "," << internalHeat << "," << radiationBalance << "," << convection << "," << waterVaporDiffusion << "," << weldingEvaporation << "," << Respiration;
	}

	void PETCalculator::subroutine_PET() {
		tx = ta;
		enbal2 = 0.0;
		count1 = 0;
		rsum_temp = rsum;
		csum_temp = csum;
		ere_temp = ere;
		return label_150();
	}
}