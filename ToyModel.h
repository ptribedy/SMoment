/*
  Class to generate number of particles for a given value of efficiency
  Copyright (C) Prithwish Tribedy - April 5,2015 
  -pls- share your comments to <ptribedy@vecc.gov.in>,<ptribedy@bnl.gov>,<prithwish2005@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

//version 1.0 modified April 20,2015 -Actual Momemnt expressions added
//version 1.1 -current version
//---------------------------------------------------------
//Example (1)
/* 
   Macro:
	#include "ToyModel.h"; 
	vector<double> eff;
	vector<double> Np;
	vector<double> Npinc;
	eff.push_back(1.);
	DeltaDist dt(10);
	dt.GetEvnt(1,eff,Np,Npinc);
*/
//
//
#ifndef ToyModel_h
#define ToyModel_h

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "TRandom.h"

#include<iostream>
using namespace std;

class ToyModel
{

	protected:
		int mdim;
		int ndim;

		vector<double> meff;
		vector<double> neff;
		vector<int> Np;
		vector<int> Npinc;

	public:
		ToyModel();
		ToyModel(const int mdim_)
		{
			mdim=mdim_;
		}

		ToyModel(const int mdim_, const int ndim_)
		{
			mdim=mdim_;
			ndim=ndim_;
		}

		virtual ~ToyModel()
		{
			vector<double>().swap(meff);
			vector<double>().swap(neff);
			vector<int>().swap(Np);
			vector<int>().swap(Npinc);
		}

		virtual vector<int> getNp() const {return Np;}
		virtual	vector<int> getNpinc() const {return Npinc;}
};


//Class DeltaDist generates particles according to a delta function ∂(Np, <Np>, \eff)
//<Np> : The mean of the distribution, \eff : The efficiency of detection, Np : The rnadom variate/multiplicity generated
//
class DeltaDist : public ToyModel{

	protected:
		vector<int> Mean;
	public:

		DeltaDist(const vector<int> Mean_) : ToyModel(0,0)
		{
			Mean=Mean_;
			cout<<"Multiplicity generation from Delta function"<<endl;
		}


		~DeltaDist()
		{
			vector<int>().swap(Mean);
		}

		void GetEvnt(const int mdim_, vector<double> meff_, vector<int>& Np_, vector<int>& Npinc_)
		{
			mdim=mdim_;
			if(meff_.size()<mdim){
			       	cout<<"ERROR: Efficiency input is not correct"<<endl;
				abort();
			}
			meff=meff_;

			Np.clear(); Npinc.clear();

			gRandom->SetSeed(0);
			for(int i=0; i<mdim; i++)
			{

				double TEMP_Mult=0;

				TEMP_Mult=gRandom->Binomial(Mean.at(i),meff.at(i));

			//	for(int j=0; j<Mean.at(i); j++)
			//	{
			//		double TEMP_eff = gRandom->Uniform(0.,1.);
			//		//cout<<i<<" "<<j<<" "<<TEMP_eff<<" "<<TEMP_Mult<<endl;
			//		if(TEMP_eff<meff.at(i)) TEMP_Mult++;
			//	}

				//cout<<"∂(Np, <Np>="<<Mean.at(i)<<", eff="<<meff.at(i)<<") ="<<TEMP_Mult<<endl;
				//
				Np.push_back(TEMP_Mult);   Npinc.push_back(Mean.at(i));
			}

			Np_=Np; Npinc_=Npinc;
		}

		long long ActualFMoment(int norder_)
		{
			long Mean__=0L;
			for(int i=0; i<mdim; i++) Mean__ += Mean.at(i);
	
			long long TEMP_fmom=1LL;
			if(norder_>0)
				for(int i=0; i<norder_; i++) TEMP_fmom *= (Mean__-i);
			return TEMP_fmom;
		}
};



//Class PoissonDist generates particles according to a Poisson distribution P(Np, <Np>, \eff)
//
//
class PoissonDist : public ToyModel{

	protected:
		vector<double> Mean;
	public:

		PoissonDist(const vector<double> Mean_) : ToyModel(0,0)
		{
			Mean=Mean_;
			cout<<"Multiplicity generation from Poisson distribution"<<endl;
		}


		~PoissonDist()
		{
			vector<double>().swap(Mean);
		}

		void GetEvnt(const int mdim_, vector<double> meff_, vector<int>& Np_, vector<int>& Npinc_)
		{
			mdim=mdim_;
			if(meff_.size()<mdim){
			       	cout<<"ERROR: Efficiency input is not correct"<<endl;
				abort();
			}
			meff=meff_;

			Np.clear(); Npinc.clear();

			gRandom->SetSeed(0);
			for(int i=0; i<mdim; i++)
			{
				//long

				double TEMP_Mult=0;
				double TEMP_Poisson=gRandom->Poisson(Mean.at(i));

				TEMP_Mult = gRandom->Binomial(TEMP_Poisson,meff.at(i));

				//for(int j=0; j<TEMP_Poisson; j++)
				//{
				//	double TEMP_eff = gRandom->Uniform(0.,1.);
				//	//cout<<i<<" "<<j<<" "<<TEMP_eff<<" "<<TEMP_Mult<<endl;
				//	if(TEMP_eff<meff.at(i)) TEMP_Mult++;
				//}

				//cout<<"P(Np, <Np>="<<Mean.at(i)<<", eff="<<meff.at(i)<<") ="<<TEMP_Mult<<endl;
				//
				Np.push_back(TEMP_Mult);   Npinc.push_back(TEMP_Poisson);
			}

			Np_=Np; Npinc_=Npinc;
		}

		long long ActualFMoment(int norder_)
		{
			long Mean__=0L;
			for(int i=0; i<mdim; i++) Mean__ += Mean.at(i);
			long long TEMP_fmom=1LL;
			if(norder_>=0) TEMP_fmom = pow(Mean__,norder_);
			return TEMP_fmom;
		}

};


//Class GaussianDist generates particles according to a Gaussian distribution P(Np, <Np>, \eff)
//
//
class GaussianDist : public ToyModel{

	protected:
		vector<int> Mean;
		vector<int> Sigma;
	public:

		GaussianDist(const vector<int> Mean_, const vector<int> Sigma_) : ToyModel(0,0)
		{
			Mean=Mean_;
			Sigma=Sigma_;
			cout<<"Multiplicity generation from Gaussian distribution"<<endl;
		}


		~GaussianDist()
		{
			vector<int>().swap(Mean);
			vector<int>().swap(Sigma);
		}

		void GetEvnt(const int mdim_, vector<double> meff_, vector<int>& Np_, vector<int>& Npinc_)
		{
			mdim=mdim_;
			if(meff_.size()<mdim){
			       	cout<<"ERROR: Efficiency input is not correct"<<endl;
				abort();
			}
			meff=meff_;

			Np.clear(); Npinc.clear();

			gRandom->SetSeed(0);
			for(int i=0; i<mdim; i++)
			{

				double TEMP_Mult=0;
				int TEMP_Gaussian=gRandom->Gaus(Mean.at(i),Sigma.at(i));

				TEMP_Mult = gRandom->Binomial(TEMP_Gaussian,meff.at(i));

				//for(int j=0; j<TEMP_Gaussian; j++)
				//{
				//	double TEMP_eff = gRandom->Uniform(0.,1.);
				//	//cout<<i<<" "<<j<<" "<<TEMP_eff<<" "<<TEMP_Mult<<endl;
				//	if(TEMP_eff<meff.at(i)) TEMP_Mult++;
				//}

				//cout<<"P(Np, <Np>="<<Mean.at(i)<<", eff="<<meff.at(i)<<") ="<<TEMP_Mult<<endl;
				//
				Np.push_back(TEMP_Mult);   Npinc.push_back(TEMP_Gaussian);
			}

			Np_=Np; Npinc_=Npinc;
		}

		long long ActualFMoment(int norder_)
		{
			long Mean__=0L;
			long Sigma__=0L;
			for(int i=0; i<mdim; i++) {Mean__ += Mean.at(i); Sigma__ += pow(Sigma.at(i),2);}
			Sigma__ = sqrt(Sigma__);

			long long TEMP_fmom=1LL;

			if(norder_==0) TEMP_fmom = 1.;

			if(norder_==1) TEMP_fmom = Mean__;

			if(norder_==2) TEMP_fmom = -Mean__ + pow(Mean__,2) + pow(Sigma__,2);
			if(norder_==3) TEMP_fmom = 2*Mean__ - 3*(pow(Mean__,2) + pow(Sigma__,2)) +  Mean__*(pow(Mean__,2) + 3*pow(Sigma__,2));

			if(norder_==4) TEMP_fmom = -6*Mean__ + pow(Mean__,4) + 6*pow(Mean__,2)*pow(Sigma__,2) + 3*pow(Sigma__,4) + 11*(pow(Mean__,2) + pow(Sigma__,2)) - 6*Mean__*(pow(Mean__,2) + 3*pow(Sigma__,2));


			if(norder_==5) TEMP_fmom = 24*Mean__ - 50*(pow(Mean__,2) + pow(Sigma__,2)) + 35*Mean__*(pow(Mean__,2) + 3*pow(Sigma__,2)) - 10*(pow(Mean__,4) + 6*pow(Mean__,2)*pow(Sigma__,2) + 3*pow(Sigma__,4)) + Mean__*(pow(Mean__,4) + 10*pow(Mean__,2)*pow(Sigma__,2) + 15*pow(Sigma__,4));

			if(norder_>5) cout<<"ERROR: can't calculate moment >5"<<endl;

			return TEMP_fmom;
		}


};



#endif
