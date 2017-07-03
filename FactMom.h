/*
  Class to calculate factorial moments of different orders
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

//Version 1.0 -modified April 9,2015 // Stirling series added
//version 1.1 -modified April 11,2015 // CentVec class added
//version 1.2 -modified April 16,2015 // Truncation problem reduced
//version 1.3 -modified April 17,2015 //
//version 1.4 -modified April 27,2015 // CentVecN class added
//version 1.5 -modified April 28,2015 // MomFactN upto 4th order added
//version 1.6 -modified April 29,2015 // FactMom for any order added
//version 1.7 -modified May 2,2015 // Bug removed for array dim ->0
//version 1.8 -modified May 5,2015 // Correlation to Central moment conversion added
//version 1.9 -modified May 9,2015 // Calculation of E-by-E central moments added
//version 2.0 -modified May 30,2015 // Bug removed for correlation variable & CalcFactVec added
//version 2.1 -modified June 20,2015 // Error calculation for central moments added
//version 2.2 -modified July 4,2015 // Notation for printD() n>9 modified & CalcFactVec for Error added
//version 2.3 -modified July 12,2015 // Cumulant vector class modified 
//version 2.4 -modified July 14,2015 // Cumulant vector class modified for conversion to other moments
//version 2.5 -modified July 18,2015 // Cumulant class modified and CalcFactVec added
//version 2.6 -modified July 19,2015 // Error of Cumulant function added
//version 2.7 -modified Aug 4,2015 //AddFactVec function added
//version 2.8 -modified Aug 15,2015 //MomVecN added
//version 2.9 -modified Sept 15,2015 //SIC functionality implemented
//version 3.0 -modified Sept 20,2015 //Cumulant Ratio class added
//version 3.1 -modified Dec 10,2015
//version 3.2 -modified Dec 21,2015 //FactVec assignment with diff expo modified 
//version 3.3 -modified Dec 28,2015 //CorrErr class modified
//version 3.4 -modified July 1,2017 //KochRatio class added
//version 3.5 -modified July 3,2017 //KochRatio class modified
//version 3.6 -current version
#ifndef FactMom_h
#define FactMom_h

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <iterator>
#include<iostream>

using namespace std;

class FactVec;
class FactorialMoment;


int Mymin (int a_, int b_) {return (a_>b_) ? b_:a_;}
int Mymax (int a_, int b_) {return (a_>b_) ? a_:b_;}

class FactVec
{
	private:
		int ndim;
		vector<int> fmnab;
		//int fexp;
		float fexp;
	public:
		FactVec(const int ndim_=0)
		{
			ndim=ndim_;
			if(ndim!=0){
				for (int i=0; i<ndim; i++) {fmnab.push_back(0);}
			}
			fexp=1;
		}

		~FactVec()
		{
			vector<int>().swap(fmnab);
		}

		void setV(const vector<int> V_){fmnab = V_; ndim = static_cast<int>(V_.size());}

		void setD(const int i_, const int D_){if(i_>ndim-1) cout<<"ERROR: out side FactVec size"<<endl; fmnab.at(i_)=D_;}
		void addD(const int i_, const int D_){if(i_>ndim-1) cout<<"ERROR: out side FactVec size"<<endl; fmnab.at(i_) += D_;}
		

		//void setE(const int E_){fexp=E_;}
		//void addE(const int E_){fexp += E_;}

		void setE(const float E_){fexp=E_;}
		void addE(const float E_){fexp += E_;}




		int  sumD() const {int sum__=0; for(int i_=0; i_<ndim; i_++) sum__ += fmnab.at(i_); return sum__;}
		int  getD(const int i_) const {return fmnab.at(i_);}

		int  getE() const {return fexp;}

		int  size() const {return ndim;}
		//int  expo() const {return fexp;}
		float  expo() const {return fexp;}


		void clear(){fmnab.clear(); ndim=0; fexp=1;};

		void resize(const int ndim_=0)
		{
			this->clear();
			ndim=ndim_;
			if(ndim!=0){
				for (int i=0; i<ndim; i++) {fmnab.push_back(0);}
			}
		}


		void printD()
		{
			for(int i__=0; i__<static_cast<int>(fmnab.size()); i__++) {
				if(fmnab.at(i__)>9)
				{
					if(i__==0)
					{
						cout<<fmnab.at(i__)<<",";
					}else if(i__ != (static_cast<int>(fmnab.size())-1)){
						cout<<","<<fmnab.at(i__)<<",";
					}else{
						cout<<","<<fmnab.at(i__);
					}
				}else{
					cout<<fmnab.at(i__);
				}
			}
		}

		FactVec& operator = (const FactVec& a)
		{
			if(a.size() != ndim) {
				fprintf(stderr," ERROR: Setting vectors of differnt dimensions");
			}else if(a.expo() != fexp) {
			//	fprintf(stderr," WARNING: Setting vectors of differnt exponents");

			//	cout<<" ";
				//this->printD();
			//	cout<<"^"<<this->expo()<<"≠";
				//a.printD();
			//	cout<<"^"<<a.expo()<<" ";

				for(int i=0; i<ndim; i++){
					fmnab.at(i) = a.getD(i);
				}
				fexp = a.expo();

			}else{
				for(int i=0; i<ndim; i++){
					fmnab.at(i) = a.getD(i);
				}
			}
			return *this;
		}

		FactVec& operator += (FactVec& a)
		{
			if(a.size() != ndim) {
				fprintf(stderr," ERROR: Adding vectors of differnt dimensions");
			}else if(a.expo() != fexp) {
				fprintf(stderr," ERROR: Adding vectors of differnt exponents");
			}else{
				for(int i=0; i<ndim; i++){
					fmnab.at(i) += a.getD(i);
				}
			}
			return *this;
		}

		FactVec& operator *= (FactVec& a)
		{
			if(a.size() != ndim) {
				fprintf(stderr,"ERROR: Multiplying vectors of differnt dimensions");
			}else if(a.expo() != fexp) {
				fprintf(stderr,"ERROR: Multiplying vectors of differnt exponents");
			}else{
				for(int i=0; i<ndim; i++){
					fmnab.at(i) += a.getD(ndim-i-1);
				}
			}
			return *this;
		}


		FactVec& operator *= (const double a)
		{
			for(int i=0; i<ndim; i++){
				fmnab.at(i) *= a;
			}

			return *this;
		}


		bool operator == (FactVec& a)
		{
			if(a.size() != ndim) return false;
			for(int i=0; i<a.size(); i++)
				if(fmnab.at(i) != a.getD(i)) return false;
		//	if(a.expo() != fexp) return false;
			
			return true;
		}

		bool operator != (FactVec& a)
		{
			if(a.size() != ndim) return true;
				for(int i=0; i<a.size(); i++)
					if(fmnab.at(i) != a.getD(i)) return true;
		//	if(a.expo() != fexp) return true;

			return false;
		}

};


      bool sortsumD (FactVec fvv_, FactVec _fvv) {
	      bool __GE=0; 
	      if(fvv_.sumD()>_fvv.sumD()){
		      __GE=1;
	      }else if(fvv_.sumD()==_fvv.sumD()){ 
		      __GE=(fvv_.getD(0)>_fvv.getD(0));
	      }else{
		      __GE=0;
	      }
	      return __GE;
      }


class FactorialMoment
{

	private:
		int morder;
		int norder;

		vector<int> Mval;
		vector<int> Nval;
		vector<int> mpow;
		vector<int> npow;
		vector<double> meff;
		vector<double> neff;

		long double TFmnInc;
		long long TFmnObs;

		string VERBOSE;

	public:
		FactorialMoment(const int morder_=0, const int norder_=0)
		{
			morder=morder_;
			norder=norder_;

			if(morder!=0){
				for (int i=0; i<morder; i++) {Mval.push_back(0); mpow.push_back(0); meff.push_back(1.);}
			}

			if(norder!=0){
				for (int i=0; i<norder; i++) {Nval.push_back(0); npow.push_back(0); neff.push_back(1.);}
			}

		}

		~FactorialMoment()
		{
			vector<int>().swap(Mval);
			vector<int>().swap(Nval);
			vector<int>().swap(mpow);
			vector<int>().swap(npow);
			vector<double>().swap(meff);
			vector<double>().swap(neff);
		}

		void setVERBOSE(const string verbose_){VERBOSE=verbose_;};

		long double getFmnInc() const {return TFmnInc;}
		long long getFmnObs() const {return TFmnObs;}

		long long calcfact(int Nval_, int nval_)
		{
			long long TEMP_FACT=1LL;

			for(int i=0; i<nval_; i++) TEMP_FACT *= (Nval_-i);
			return TEMP_FACT;
		}


		void add(const vector<int> Mval_, const vector<int> Nval_)
		{
			Mval=Mval_; Nval=Nval_;
			fprintf(stderr,"Multiplicities for %d bins are =",morder+norder);
			if(morder!=0) for (int i=0; i<morder; i++) {fprintf(stderr,"%d,", Mval.at(i));}
			if(norder!=0) for (int i=0; i<norder-1; i++) {fprintf(stderr,"%d,", Nval.at(i));}
			if(norder!=0) fprintf(stderr,"%d \n", Nval.at(norder-1));
		}

		void add(const vector<int> Mval_, const vector<int> Nval_, const vector<double> meff_, const vector<double> neff_)
		{

			if(VERBOSE=="LOUD"){
				cout<<" I am here "<<endl;

				cout<<"Sizes Np, Npbar_, Npeff, Npbareff_ = "<<Mval_.size()<<" "<<Nval_.size()<<" "<<meff_.size()<<" "<<neff_.size()<<endl;

				cout<<"morder= "<<morder<<endl;
				cout<<"norder= "<<norder<<endl;
			}

			Mval=Mval_; Nval=Nval_; meff=meff_; neff=neff_;

			if(VERBOSE=="LOUD"){
				fprintf(stderr,"\n Multiplicities for %d bins are =",morder+norder);
				if(morder!=0) for (int i=0; i<morder; i++) {fprintf(stderr,"%d,", Mval.at(i));}
				if(norder!=0) for (int i=0; i<norder-1; i++) {fprintf(stderr,"%d,", Nval.at(i));}
				if(norder!=0) fprintf(stderr,"%d \n", Nval.at(norder-1));

				fprintf(stderr,"\n Efficiencies for %d bins are =",morder+norder);
				if(morder!=0) for (int i=0; i<morder; i++) {fprintf(stderr,"%g,", meff.at(i));}
				if(norder!=0) for (int i=0; i<norder-1; i++) {fprintf(stderr,"%g,", neff.at(i));}
				if(norder!=0) fprintf(stderr,"%g \n", neff.at(norder-1));
			}
		}

		void calcmom(const vector<int> mpow_, const vector<int> npow_)
		{
			if(VERBOSE=="LOUD"){
				cout<<"morder= "<<morder<<endl;
				cout<<"norder= "<<norder<<endl;

				cout<<endl;
				cout<<"mpow= ";
				for(size_t k_=0; k_<mpow_.size(); k_++) cout<<mpow_.at(k_);
				cout<<endl;

				cout<<endl;
				cout<<"npow= ";
				for(size_t k_=0; k_<npow_.size(); k_++) cout<<npow_.at(k_);
				cout<<endl;
			}

			mpow=mpow_; npow=npow_;

			long long TEMPFmnObs=1LL;
			long double TEMPFmnInc=1.0L;

			for (int i=0; i<morder; i++) 
			{
				TEMPFmnObs *= calcfact(Mval.at(i),mpow.at(i));
				TEMPFmnInc *= static_cast<long double>(calcfact(Mval.at(i),mpow.at(i)))/pow(meff.at(i),mpow.at(i));
			}


			for (int i=0; i<norder; i++) 
			{

				TEMPFmnObs *= calcfact(Nval.at(i),npow.at(i));
				TEMPFmnInc *= static_cast<long double>(calcfact(Nval.at(i),npow.at(i)))/pow(neff.at(i),npow.at(i));
			}

			if(VERBOSE=="LOUD"){
				cout<<"f_{";
				if(morder!=0) {for (int i=0; i<morder; i++) cout<<mpow.at(i);}
				if(norder!=0) {for (int i=0; i<norder; i++) cout<<npow.at(i);}
				cout<<"}="<<TEMPFmnObs<<" ";

				cout<<"F_{";
				if(morder!=0) {for (int i=0; i<morder; i++) cout<<mpow.at(i);}
				if(norder!=0) {for (int i=0; i<norder; i++) cout<<npow.at(i);}
				cout<<"}="<<TEMPFmnInc<<"("<<static_cast<long long>(round(TEMPFmnInc))<<")"<<endl;			
			}

			TFmnInc = TEMPFmnInc;
			TFmnObs = TEMPFmnObs;
		}




};

class FactMom
{
	private:
		int mth;
		int nth;
		int abin;
		int bbin;

	protected:
		int Norder;
		int Nelements;
		int Nsize;
		string VERBOSE;

		vector<FactVec> Fmna;
		vector<FactVec> Fmnb;
		vector<FactVec> Fmnab;

	public:

		FactMom(string VERBOSE_){VERBOSE=VERBOSE_;}

		//FactMom(const int mth_, const int nth_, const int abin_, const int bbin_);
		FactMom(const int mth_, const int nth_, const int abin_, const int bbin_, string VERBOSE_="");


		//FactMom::FactMom(const int mth_, const int nth_, const int abin_, const int bbin_, vector<int> Np, vector<int> Npbar,  vector<double> Npeff, vector<double> Npbareff)

		virtual ~FactMom()
		{
			vector<FactVec>().swap(Fmna);
			vector<FactVec>().swap(Fmnb);
			vector<FactVec>().swap(Fmnab);
		}

		virtual void setVERBOSE(string verbose_){VERBOSE=verbose_;};

		void SetFmnab();

		void GetMoment(vector<int> Np, vector<int> Npbar, vector<double> Npeff, vector<double> Npbareff, long long& ObsMom_, long double& IncMom_);
	//	void GetMoment(vector<int> Np,  vector<double> Npeff, long long& ObsMom_, long double& IncMom_);

		int get_mth() const {return mth;}
		int get_nth() const {return nth;}
		int get_abin() const {return abin;}
		int get_bbin() const {return bbin;}

};




//FactMom::FactMom(const int mth_, const int nth_, const int abin_, const int bbin_, vector<int> Np, vector<int> Npbar,  vector<double> Npeff, vector<double> Npbareff)



FactMom::FactMom(const int mth_, const int nth_, const int abin_, const int bbin_, string VERBOSE_)
{

	        mth=mth_;
		nth=nth_;
		abin=abin_;
		bbin=bbin_;
        	
	        VERBOSE=VERBOSE_;
        	
	        if(abin==0){
			abin=bbin;
			bbin=0;
		}
        	
	        if(mth==0){
			mth=nth;
			nth=0;
		}
        	
        
	        Norder=(abin + bbin); 
		Nelements=(pow(abin,mth)*pow(bbin,nth)); 
        	
	        fprintf(stderr,"FactMom class set F_{m=%d,n=%d} , Bins (%d x %d), Total Order = %d, No. of elements = %d \n",mth,nth,abin,bbin,Norder,Nelements);
        	
	        SetFmnab();
        	
	        if(VERBOSE=="LOUD"){
			fprintf(stderr,"F_{%d %d %d %d} = ",mth,nth,abin,bbin);
        	
	        	for(int k=0; k<static_cast<int>(Fmnab.size()); k++){
				fprintf(stderr,"F_{");
				for(int l=0; l<Fmnab.at(k).size(); l++) fprintf(stderr,"%d",Fmnab.at(k).getD(l));
        	
	        		if(k==static_cast<int>(Fmnab.size())-1) {fprintf(stderr,"} \n");
				}else{
					fprintf(stderr,"} + ");
				}
			}
		}
        	
	        /*
		   long double IncMom = 0.0L;
		   long long ObsMom = 0LL;
        	
	           GetMoment(Np, Npbar,Npeff, Npbareff, ObsMom, IncMom);
        	
	           cout<<"_________________________________________________________________"<<endl;
		   cout<<"f_{"<<mth<<nth<<abin<<bbin<<"}= "<<ObsMom<<" ";
		   cout<<"F_{"<<mth<<nth<<abin<<bbin<<"}= "<<IncMom<<"("<<static_cast<long long>(round(IncMom))<<")"<<endl;
		   */
        	
}

void FactMom::GetMoment(vector<int> Np, vector<int> Npbar, vector<double> Npeff, vector<double> Npbareff, long long& ObsMom_, long double& IncMom_)
{

		FactorialMoment fmn(abin,bbin); //= new FactorialMoment(abin,bbin);
        
		fmn.add(Np,Npbar,Npeff,Npbareff);
		fmn.setVERBOSE(VERBOSE);
        
		vector<int> m; vector<int> n;
		long double IncMom = 0.0L;
		long long ObsMom = 0LL;
        
		for(int k=0; k<static_cast<int>(Fmnab.size()); k++){
        
			m.clear();
			n.clear();
        
			for(int i=0; i<abin; i++){m.push_back(Fmnab.at(k).getD(i));}
			for(int j=0; j<bbin; j++){n.push_back(Fmnab.at(k).getD(j+abin));}
        
			fmn.calcmom(m,n);
        
			ObsMom += fmn.getFmnObs();
			IncMom += fmn.getFmnInc();
		}
        
		//cout<<"_________________________________________________________________"<<endl;
		//cout<<"f_{"<<mth<<nth<<abin<<bbin<<"}= "<<ObsMom<<" ";
		//cout<<"F_{"<<mth<<nth<<abin<<bbin<<"}= "<<IncMom<<"("<<static_cast<long long>(round(IncMom))<<")"<<endl;
        
		ObsMom_=ObsMom; IncMom_=IncMom;
		vector<int>().swap(m); vector<int>().swap(n);
}


void FactMom::SetFmnab()
{

	for(int i=0; i<abin; i++)
	{
		FactVec TEMP_FA(Norder);
		TEMP_FA.setD(i,1);
		Fmna.push_back(TEMP_FA);
		//		fprintf(stderr,"Should be called 1 \n");
	}

	for(int j=0; j<bbin; j++)
	{
		FactVec TEMP_FB(Norder);
		TEMP_FB.setD(j+abin,1);
		Fmnb.push_back(TEMP_FB);
		//		fprintf(stderr,"Should be called 2 \n");
	}


	vector<FactVec> Fmnab_;

	//Set up a new vector with Nelements (= a^m x b^n) dimesions 

	//Step 1: Set Fmnab[a] = Fmna[a] 


	for(int j=0; j<static_cast<int>(Fmna.size()); j++){

		FactVec TEMP_FAB(Norder);
		TEMP_FAB += Fmna.at(j);
		Fmnab_.push_back(TEMP_FAB);
		//		fprintf(stderr,"Should be called 3 \n");
	}


	//Step 2: Multiply Fmnab[a] with Fmna[a] for m-1 times

	int OLD_SIZE=0;
	int NEW_SIZE=0;

	for(int i=0; i<mth-1; i++){


		NEW_SIZE=static_cast<int>(Fmnab_.size());

		for(int j=OLD_SIZE; j<NEW_SIZE; j++){

			for(int k=0; k<static_cast<int>(Fmna.size()); k++){


				FactVec TEMP_FAB(Norder);

				TEMP_FAB += Fmnab_.at(j);

				TEMP_FAB += Fmna.at(k);

				Fmnab_.push_back(TEMP_FAB);

				//				fprintf(stderr,"\n %d (%d) \n",Fmnab_.size(),Nelements);
			}
		}

		OLD_SIZE=NEW_SIZE;

	}



	//Step 3: Multiply Fmnab[axa] with Fmnb[b] for n times


	for(int i=0; i<nth; i++){


		NEW_SIZE=static_cast<int>(Fmnab_.size());

		for(int j=OLD_SIZE; j<NEW_SIZE; j++){

			for(int k=0; k<static_cast<int>(Fmnb.size()); k++){


				FactVec TEMP_FAB(Norder);

				TEMP_FAB += Fmnab_.at(j);

				TEMP_FAB += Fmnb.at(k);

				Fmnab_.push_back(TEMP_FAB);

				//				fprintf(stderr,"\n %d (%d) \n",Fmnab_.size(),Nelements);
			}
		}

		OLD_SIZE=NEW_SIZE;

	}

	Nsize=static_cast<int>(Fmnab_.size())-Nelements;

	for(int i=0; i<static_cast<int>(Fmnab_.size()); i++)
	{
		if(i>=Nsize){
			Fmnab.push_back(Fmnab_.at(i));
		}
	}

	vector<FactVec>().swap(Fmnab_);
}



class FactMomN : public FactMom
{
	private:

		vector<int> mthN;
		vector<int> abinN;
		vector<vector<FactVec> > FmnN;

	public:

		FactMomN(string VERBOSE_="") : FactMom(VERBOSE_)
		{

		}

		FactMomN(const vector<int> mthN_, const vector<int> abinN_, string VERBOSE_="") : FactMom(VERBOSE_)
		{	
			if(mthN_.size() != abinN_.size()) cout<<"ERROR : Input not set correctly"<<endl;

			for(size_t k_=0; k_<abinN_.size();k_++){
				mthN.push_back(mthN_.at(k_));
				abinN.push_back(abinN_.at(k_));
			}

			Norder=0;
			Nelements =1;
			for(size_t k_=0; k_<abinN.size(); k_++) 
			{
				Norder += abinN.at(k_);
				Nelements *= (pow(abinN.at(k_),mthN.at(k_))); 
			}

			if(VERBOSE=="LOUD"){
				cout<<"FactMom class set F_{";
				for(size_t k_=0; k_<mthN.size(); k_++) cout<<mthN.at(k_); 
				cout<<"} , Bins (";
				for(size_t k_=0; k_<abinN.size(); k_++) cout<<abinN.at(k_); 
				cout<<"), Total Order = ";
				cout<<Norder;
				cout<<", No. of elements = ";
				cout<<Nelements<<endl;
			}

			FmnN.resize(abinN.size()-1);

			SetFmnabN();
		}


		~FactMomN()
		{
			vector<int>().swap(mthN);
			vector<int>().swap(abinN);
		}


		void SetFmnabN();


		void SetFmnabN(const vector<int> mthN_, const vector<int> abinN_)
		{
			if(mthN_.size() != abinN_.size()) cout<<"ERROR : Input not set correctly"<<endl;

			abinN.clear(); mthN.clear();

			for(size_t k_=0; k_<abinN_.size();k_++){
				mthN.push_back(mthN_.at(k_));
				abinN.push_back(abinN_.at(k_));
			}

			Norder=0;
			Nelements =1;
			for(size_t k_=0; k_<abinN.size(); k_++) 
			{
				Norder += abinN.at(k_);
				Nelements *= (pow(abinN.at(k_),mthN.at(k_))); 
			}

			if(VERBOSE=="LOUD"){
				cout<<"FactMom class set F_{";
				for(size_t k_=0; k_<mthN.size(); k_++) cout<<mthN.at(k_); 
				cout<<"} , Bins (";
				for(size_t k_=0; k_<abinN.size(); k_++) cout<<abinN.at(k_); 
				cout<<"), Total Order = ";
				cout<<Norder;
				cout<<", No. of elements = ";
				cout<<Nelements<<endl;
			}

			Fmnab.clear();
			Fmna.clear();
			FmnN.clear();

			FmnN.resize(abinN.size()-1);

			SetFmnabN();
		}




//		void GetMoment(vector<vector<int> > Np, vector<vector<double> > Npeff, long long& ObsMom_, long double& IncMom_);

		void GetMomentN(const vector<int> Np, const vector<double> Npeff, long long& ObsMom_, long double& IncMom_);


		void printD()
		{
			if(mthN.size()==0 || abinN.size()==0){
				cout<<"ERROR : Empty bins found"<<endl; 
			}else{
				if(VERBOSE=="LOUD")
				{
					cout<<endl;
					cout<<"F_{";
					for(size_t k_=0; k_<mthN.size(); k_++) cout<<mthN.at(k_)<<" ";  
					cout<<", ";
					for(size_t k_=0; k_<abinN.size()-1; k_++) cout<<abinN.at(k_)<<" ";
					cout<<abinN.at(abinN.size()-1);
					cout<<"} =";
				}
				for(size_t k_=0; k_<Fmnab.size(); k_++){
					cout<<"F_{";
					for(int l=0; l<Fmnab.at(k_).size(); l++) cout<<Fmnab.at(k_).getD(l);

					if(k_==(Fmnab.size()-1)) {
						cout<<"}";
						if(VERBOSE=="LOUD") cout<<endl;
					}else{
						cout<<"} + ";
					}
				}

				if(VERBOSE=="LOUD") cout<<endl;
			}
		}

		void printD(string FORMAT_)
		{
			if(mthN.size()==0 || abinN.size()==0){
				cout<<"ERROR : Empty bins found"<<endl; 
			}else{
				if(VERBOSE=="LOUD")
				{
					cout<<endl;
					cout<<"F_{";
					for(size_t k_=0; k_<mthN.size(); k_++) cout<<mthN.at(k_)<<" ";  
					cout<<", ";
					for(size_t k_=0; k_<abinN.size()-1; k_++) cout<<abinN.at(k_)<<" ";
					cout<<abinN.at(abinN.size()-1);
					cout<<"} =";
				}

				int counter=0;
				for(size_t k_=0; k_<Fmnab.size(); k_++){

					if(FORMAT_=="latex"){cout<<"$frac{{$rm f}_{";}else{cout<<"f_{";}

					for(int l=0; l<Fmnab.at(k_).size(); l++) cout<<Fmnab.at(k_).getD(l);

					if(FORMAT_=="latex")cout<<"}}{";

					if(k_==(Fmnab.size()-1)) {

						if(FORMAT_=="eff")
						{       
							cout<<"}";
							for(int l_=0; l_<Fmnab.at(k_).size(); l_++)
							{
								if(Fmnab.at(k_).getD(l_)!=0){
									cout<<"/€";
									cout<<l_+1;
									if(Fmnab.at(k_).getD(l_)>1) cout<<"^"<<Fmnab.at(k_).getD(l_);
								}
							}
						}else if(FORMAT_=="latex"){
							for(int l_=0; l_<Fmnab.at(k_).size(); l_++)
							{
								if(Fmnab.at(k_).getD(l_)!=0){
									cout<<"$varepsilon_"<<l_+1;
									if(Fmnab.at(k_).getD(l_)>1) cout<<"^"<<Fmnab.at(k_).getD(l_);
								}
							}
							cout<<"}";
							counter ++;
						}




						if(VERBOSE=="LOUD") cout<<endl;
					}else{

						if(FORMAT_=="eff")
						{       
							cout<<"}";
							for(int l_=0; l_<Fmnab.at(k_).size(); l_++)
							{
								if(Fmnab.at(k_).getD(l_)!=0){
									cout<<"/€";
									cout<<l_+1;
									if(Fmnab.at(k_).getD(l_)>1) cout<<"^"<<Fmnab.at(k_).getD(l_);
								}
							}
							cout<<"	+ ";

						}else if(FORMAT_=="latex"){
							//	cout<<"}";
							for(int l_=0; l_<Fmnab.at(k_).size(); l_++)
							{
								if(Fmnab.at(k_).getD(l_)!=0){
									cout<<"$epsilon_"<<l_+1;
//									cout<<"~";
									if(Fmnab.at(k_).getD(l_)>1) cout<<"^"<<Fmnab.at(k_).getD(l_);
								}
							}
							//cout<<"}";
							cout<<"	} + ";
							counter ++;
						}

					}

					if(FORMAT_=="latex"){
						//if(k_!=0 && k_%2==0){
						if(counter%3==0){
							cout<<" $$";
							cout<<endl;
							cout<<" && ";
						}
					}
				}

				if(VERBOSE=="LOUD") cout<<endl;
			}
		}

};


		void FactMomN::SetFmnabN()
		{
			//cout<<"abin.size()= "<<abinN.size()<<endl;
			//cout<<"abinN.at(0).size()= "<<abinN.at(0)<<endl;

			//cout<<"("<<Fmna.size()<<") Fmna ="<<endl;

			for(int i__=0; i__<abinN.at(0); i__++)
			{
				//cout<<"i= "<<i__<<endl;
				FactVec TEMP_FA_(Norder);
				bool TEMP_D = mthN.at(0);
				TEMP_FA_.setD(i__,TEMP_D);
				Fmna.push_back(TEMP_FA_);
				//cout<<endl;
				//cout<<"TEMP_D= "<<TEMP_D<<endl;
				//TEMP_FA_.printD();
				//cout<<endl;
			}

			//cout<<endl;
			//cout<<"("<<Fmna.size()<<") Fmna ="<<endl;


			//for(size_t k_=0; k_<Fmna.size(); k_++)
			//{
			//	cout<<endl;		
			//	Fmna.at(k_).printD();
			//	cout<<endl;		
			//}

			int TEMP_POS = Fmna.size();

			//cout<<"TEMP_POS= "<<TEMP_POS<<endl;

			for(size_t k_=1; k_<abinN.size(); k_++)
				for(int j=0; j<abinN.at(k_); j++)
				{
					bool TEMP_D = mthN.at(k_);
					//
					FactVec TEMP_FB(Norder);
					TEMP_FB.setD(TEMP_POS,TEMP_D);
					FmnN.at(k_-1).push_back(TEMP_FB);
					TEMP_POS ++;
					//FmnN.at(k_-1).printD();
				}

			//cout<<endl;
			//cout<<"FmnN ="<<endl;
			//for(size_t k_=0; k_<FmnN.size(); k_++)
			//	for(size_t j_=0; j_<FmnN.at(k_).size(); j_++){cout<<"k_="<<k_<<" j_="<<j_<<" "; FmnN.at(k_).at(j_).printD(); cout<<endl;}


			vector<FactVec> Fmnab_;

			for(size_t j=0; j<Fmna.size(); j++)
			{
				FactVec TEMP_FAB(Norder);
				TEMP_FAB += Fmna.at(j);
				Fmnab_.push_back(TEMP_FAB);
			}



			size_t OLD_SIZE=0;
			size_t NEW_SIZE=0;

			for(int i=0; i<mthN.at(0)-1; i++){

				NEW_SIZE=Fmnab_.size();

				for(size_t j=OLD_SIZE; j<NEW_SIZE; j++){

					for(size_t k=0; k<Fmna.size(); k++){


						FactVec TEMP_FAB(Norder);

						TEMP_FAB += Fmnab_.at(j);

						TEMP_FAB += Fmna.at(k);

						Fmnab_.push_back(TEMP_FAB);

					}
				}

				OLD_SIZE=NEW_SIZE;

			}



			for(size_t l__=1; l__<mthN.size(); l__++){

				for(int i=0; i<mthN.at(l__); i++){


					NEW_SIZE=Fmnab_.size();

					for(size_t j=OLD_SIZE; j<NEW_SIZE; j++){

						for(size_t k=0; k<FmnN.at(l__-1).size(); k++){


							FactVec TEMP_FAB(Norder);

							TEMP_FAB += Fmnab_.at(j);

							TEMP_FAB += FmnN.at(l__-1).at(k);

							Fmnab_.push_back(TEMP_FAB);
						}
					}

					OLD_SIZE=NEW_SIZE;
				}

			}


			//cout<<"Fmnab_ ="<<endl;
			//for(size_t k_=0; k_<Fmnab_.size(); k_++)
			//{
			//	cout<<endl;		
			//	Fmnab_.at(k_).printD();
			//	cout<<endl;		
			//}



			Nsize=static_cast<int>(Fmnab_.size())-Nelements;

			for(size_t i__=0; i__<Fmnab_.size(); i__++)
			{
				if(i__>=static_cast<unsigned int>(Nsize))
					Fmnab.push_back(Fmnab_.at(i__));
			}

			vector<FactVec>().swap(Fmnab_);
		}


		void FactMomN::GetMomentN(const vector<int> Np,  const vector<double> Npeff, long long& ObsMom_, long double& IncMom_)
		{

			int isize =0;

			if(VERBOSE=="LOUD")
				cout<<endl;
				for(size_t k_=0; k_<abinN.size(); k_++) {isize += abinN.at(k_); if(VERBOSE=="LOUD") cout<<" "<<abinN.at(k_)<<" ";}
			if(VERBOSE=="LOUD")	cout<<endl;
			if(VERBOSE=="LOUD")	cout<<" isize= "<<isize <<endl;
		
		
			if(Np.size() != static_cast<unsigned int>(isize)){cout<<endl; cout<<"ERROR: Input multiplicity array dimension mismath" <<endl;}
			if(Npeff.size() != static_cast<unsigned int>(isize)){cout<<endl; cout<<"ERROR: Input efficiency array dimension mismath" <<endl;}

			FactorialMoment fmn(isize); //= new FactorialMoment(abin,bbin);

			vector<int> Npbar_;
			vector<double> Npbareff_;

//			for(size_t k_=0; k_<abinN.size(); k_++) {Npbar_.push_back(0); Npbareff_.push_back(1);}
			
//			cout<<"Sizes abinN= Np, Npbar_, Npeff, Npbareff_ = "<<abinN.size()<<" "<<Np.size()<<" "<<Npbar_.size()<<" "<<Npeff.size()<<" "<<Npbareff_.size()<<endl;

			fmn.setVERBOSE(VERBOSE);

			fmn.add(Np,Npbar_,Npeff,Npbareff_);

			vector<int> m; vector<int> n;
			long double IncMom = 0.0L;
			long long ObsMom = 0LL;

			for(size_t k_=0; k_<Fmnab.size(); k_++){

				m.clear();
//				n.clear();
				
//				for(size_t l_=0; l_<abinN.size(); l_++)
//				for(int i__=0; i__<abinN.at(l_); i__++) 
//
				for(int i__=0; i__<Fmnab.at(k_).size(); i__++) 
					m.push_back(Fmnab.at(k_).getD(i__));

//				for(int l=0; l<Fmnab.at(k_).size(); l++) cout<<Fmnab.at(k_).getD(l);
//
//				cout<<endl;
//				for(size_t j_=0; j_<m.size(); j_++) cout<<m.at(j_);
//				cout<<endl;

				fmn.calcmom(m,n);

				ObsMom += fmn.getFmnObs();
				IncMom += fmn.getFmnInc();

			}

			if(VERBOSE=="LOUD"){
			      cout<<endl;
			      cout<<"_________________________________________________________________"<<endl;
			      cout<<"f_{<<mth<<nth<<abin<<bbin<<}= "<<ObsMom<<" ";
			      cout<<"F_{<<mth<<nth<<abin<<bbin<<}= "<<IncMom<<"("<<static_cast<long long>(round(IncMom))<<")"<<endl;
			}

			ObsMom_=ObsMom; IncMom_=IncMom;
			vector<int>().swap(m); vector<int>().swap(n);
			vector<int>().swap(Npbar_);
			vector<double>().swap(Npbareff_);
		}



class CentVec
{

	protected:
		int ndim;
		int mdim;
		int nsize;
		vector<int> f10pow;
		vector<int> f01pow;


		vector<long> cmn;
		vector<FactVec> cfmn;

		void  setcmn(const long cmn_){cmn.push_back(cmn_);}
		void  setf10pow(const int f10pow_){f10pow.push_back(f10pow_);}
		void  setf01pow(const int f01pow_){f01pow.push_back(f01pow_);}


		void  setcfmn(const FactVec cfmn_){cfmn.push_back(cfmn_);}

	public:
		CentVec(const int mdim_=0, const int ndim_=0)
		{
			mdim=mdim_;
			ndim=ndim_;
			nsize=0;
		}

		~CentVec()
		{

			vector<int>().swap(f10pow);
			vector<int>().swap(f01pow);
			vector<long>().swap(cmn);
			vector<FactVec>().swap(cfmn);
		}

		void setvec(const int j_, const int i_, const long cmnval_){

			cmn.push_back(cmnval_);
			FactVec TEMP_cfmn(2);
			TEMP_cfmn.setD(0,Mymin(mdim-j_,i_+1));
			cfmn.push_back(TEMP_cfmn);
			f10pow.push_back(j_);
			nsize++;
		}

		void setvec2(const int j_, const int k_, FactVec f1_, FactVec f2_, const long cmnval_){

			cmn.push_back(cmnval_);

			FactVec TEMP_cfmn(2);

			TEMP_cfmn += f1_;
			TEMP_cfmn *= f2_;

			cfmn.push_back(TEMP_cfmn);

			f10pow.push_back(j_);
			f01pow.push_back(k_);

			nsize++;
		}



		long     getcmn(const int i_){return cmn.at(i_);}
		int      getf10pow(const int i_){return f10pow.at(i_);}
		int      getf01pow(const int i_){return f01pow.at(i_);}
		FactVec  getcfmn(const int i_){return cfmn.at(i_);}

		int  size() const {return nsize;}
		int  morder() const {return mdim;}
		int  norder() const {return ndim;}

};




class CentVecN
{
	protected:
		vector<FactVec> CCmn;
		vector<int> vvmn;
                string MOMENT;
		vector<vector<FactVec> >nthErr;
		vector<long> Stirling(const int order);
		
	private:
		size_t mdim;
		int nsize;

		vector<int> abinNN;
		vector<int> mthNN;

		vector<int> ndim;

		vector<vector<int> > fpow;
		vector<long> cmn;
		vector<FactVec> cfmn;
		


		void  setcmn(const long cmn_){cmn.push_back(cmn_);}
		void  setfpow(const int fpow_, const int i__){fpow.at(i__).push_back(fpow_);}
		void  setcfmn(const FactVec cfmn_){cfmn.push_back(cfmn_);}

		CentVec * Cent2Fact(const int order);

		void setvec(const vector<int> j_, vector<FactVec> fj_, long cmnval_){

			cmn.push_back(cmnval_);

			FactVec TEMP_cfmn(static_cast<int>(mdim));

			for(size_t i_=0; i_<fj_.size(); i_++)
			{       
				TEMP_cfmn += fj_.at(i_);
			}

			cfmn.push_back(TEMP_cfmn);

			for(size_t k_=0; k_<j_.size(); k_++) 
				fpow.at(k_).push_back(j_.at(k_));

			nsize++;
		}

		void clear(){
			mdim=0;
			ndim.clear();
			for(size_t i_=0;i_<fpow.size();i_++) fpow.at(i_).clear();
			fpow.clear();
			cmn.clear();
			cfmn.clear();
			ndim.clear();
			nsize=0;
		}

	public:
		CentVecN(const size_t mdim_=0, string MOMENT_="")
		{
			mdim=mdim_;
			fpow.resize(mdim);
			MOMENT=MOMENT_;
			nsize=0;
		}

		~CentVecN()
		{
			for(size_t i_=0;i_<fpow.size();i_++) vector<int>().swap(fpow.at(i_));
			vector<vector<int> >().swap(fpow);
			vector<long>().swap(cmn);
			vector<FactVec>().swap(cfmn);
			vector<int>().swap(ndim);
			vector<FactVec>().swap(CCmn);
			vector<int>().swap(vvmn);
		}


		long calccomb(int Nval_, int nval_);

		long     getcmn(const int j_){return cmn.at(j_);}
		int      getfpow(const int i__, const int j_){return fpow.at(i__).at(j_);}
		vector<int> getfpow(const int i__){return fpow.at(i__);}
		FactVec getcfmn(const int j_) const {return cfmn.at(j_);}
		vector<FactVec> getcfmn() const {return cfmn;}


		//vector<FactVec> getCCmn() const {return CCmn;};
		vector<int> getvvmn() const {return vvmn;}
		vector<FactVec> getCCmn() const {return CCmn;}

		int  size() const {return nsize;}
		size_t  norder() const {return mdim;}

		void setndim(const vector<int> ndim_){
			ndim=ndim_;
		}

		
		void Corr2CentN(const vector<int> mth__, string FORMAT="", string FORMAT2="");
		void Corr2FactN(const vector<int> mth__, string FORMAT="", string VERBOSE_="LOUD");
		void Corr2FactN(const vector<int> mth_, const vector<int> abin__, string FORMAT="",string VERBOSE_="LOUD");


		void Cent2FactN(const vector<int> morder);

		void SetFmnabN(const vector<int> abinN_) {abinNN =abinN_;}

		void GetCentMomentN(const vector<int> mth__, const vector<int> abin__,const vector<int> Np, const vector<double> Npeff, long long& ObsMom_, long double& IncMom_, string FORMAT="");


		vector<FactVec> AddFactVec(vector<FactVec> fvv1, vector<FactVec> fvv2)
		{
			vector<FactVec> fvv;
			vector<FactVec> fvv__;

			for(size_t i__=0; i__<fvv1.size(); i__++) fvv.push_back(fvv1.at(i__));
			for(size_t i__=0; i__<fvv2.size(); i__++) fvv.push_back(fvv2.at(i__));

			for(size_t j_=0; j_<fvv.size(); j_++){
				for(size_t k_=0; k_<fvv.size(); k_++) {

					if(fvv.at(k_).sumD()==0) continue;
					if(fvv.at(j_)==fvv.at(k_)) {
						fvv__.push_back(fvv.at(k_));
						fvv.at(k_) *= 0;
					}

				}
			}

			cout<<"old vec 1 :<";
			for(size_t j_=0; j_<fvv1.size(); j_++) {cout<<" F_"; fvv1.at(j_).printD();cout<<" ";}
			cout<<">"<<endl;

			cout<<"old vec 2 :<";
			for(size_t j_=0; j_<fvv2.size(); j_++) {cout<<" F_"; fvv2.at(j_).printD(); cout<<" ";}
			cout<<">"<<endl;


			cout<<"new vec :<";
			for(size_t j_=0; j_<fvv__.size(); j_++) {cout<<" F_"; fvv__.at(j_).printD(); cout<<" ";}
			cout<<">"<<endl;

			vector<FactVec>().swap(fvv);
			return fvv__;
		}


		vector<FactVec> CalcFactVec(string VERBOSE_="LOUD")//const CentVecN * cnt__)
		{

			vector<FactVec> fvv;

			if(MOMENT=="Nudyn"){
				FactVec F10(2);
				F10.setD(0,1);
				FactVec F01(2);
				F01.setD(1,1);
				fvv.push_back(F10);
				fvv.push_back(F01);
			}

			if(MOMENT=="Error"){

				for(size_t j__=0; j__<nthErr.size(); j__++)
				{
					bool check=1;

					for(size_t k_=0; k_<nthErr.at(j__).size(); k_++)
					{
						if(nthErr.at(j__).at(k_).sumD() <=1) check*=0;
					}

					if(check){
						for(size_t k_=0; k_<nthErr.at(j__).size(); k_++)
						{
							vector<int> mthE_;
							for(int l__=0; l__<nthErr.at(j__).at(k_).size(); l__++)
								mthE_.push_back(nthErr.at(j__).at(k_).getD(l__));

							CentVecN * cnt__ =new CentVecN();
							cnt__->Corr2CentN(mthE_,"NOEXPS");  

							for(size_t l__=0; l__<cnt__->getCCmn().size(); l__++){ 
								vector<int> norder__;
								for(int i__=0; i__<cnt__->getCCmn().at(l__).size(); i__++){
									norder__.push_back(cnt__->getCCmn().at(l__).getD(i__));
								}

								CentVecN * cnt_ =new CentVecN();
								cnt_->Cent2FactN(norder__);

								for(size_t k__=0; k__<cnt_->getcfmn().size(); k__++){
									fvv.push_back(cnt_->getcfmn().at(k__));
								}						
							}
						}
					}
				}

			}else{

				for(size_t l__=0; l__<this->getCCmn().size(); l__++){ 

					vector<int> norder__;
					for(int i__=0; i__<this->getCCmn().at(l__).size(); i__++){
						norder__.push_back(this->getCCmn().at(l__).getD(i__));
					}

					CentVecN * cnt_ =new CentVecN();
					cnt_->Cent2FactN(norder__);


					for(size_t k__=0; k__<cnt_->getcfmn().size(); k__++){

						if(MOMENT=="Nudyn"){
							if(cnt_->getcfmn().at(k__).sumD() == this->getCCmn().at(0).sumD()) fvv.push_back(cnt_->getcfmn().at(k__));
						}else{
							fvv.push_back(cnt_->getcfmn().at(k__));
						}
					}
				}
			}
			cout<<endl;

			for(size_t j_=0; j_<fvv.size(); j_++){
				for(size_t k_=j_+1; k_<fvv.size(); k_++) {
					if(fvv.at(j_) == fvv.at(k_)) { fvv.erase(fvv.begin() + k_); }
				}
			}

			if(MOMENT=="Error"){
				for(size_t j_=0; j_<fvv.size(); j_++){
					for(size_t k_=j_+1; k_<fvv.size(); k_++) {
						if(fvv.at(j_) == fvv.at(k_)) { fvv.erase(fvv.begin() + k_); }
					}
				}
				//fvv.erase(fvv.begin() + fvv.size()-1);
			}

			if(VERBOSE_=="LOUD"){
				cout<<"Factorial Moments to be calculated are :"<<endl;

				for(size_t j_=0; j_<fvv.size(); j_++) {cout<<"F_"; fvv.at(j_).printD(); cout<<endl;}

				cout<<endl;
				cout<<endl;
				cout<<endl;
				cout<<endl;
			}

			return fvv;
		}


		void printD(string FORMAT=""){

			if(FORMAT != "NOEXPS"){
				cout<<endl;
				cout<<endl;
				cout<<"C_";
				for(size_t k_=0; k_<ndim.size(); k_++) cout<<ndim.at(k_);
				cout<<" = ";
			}

			vector<vector<FactVec> > temp_nthC;
			vector<long > temp_coeffC;
			this->Squeeze(temp_nthC,temp_coeffC);

			for(size_t j__=0; j__<temp_nthC.size(); j__++)
			{
				if(temp_coeffC.at(j__)==1){ if(j__!=0) cout<<" + ";
				}else if(temp_coeffC.at(j__)==-1){ cout<<" - ";
				}else if(temp_coeffC.at(j__)<0){
					cout<<" - "<<fabs(temp_coeffC.at(j__));
				}else{
					cout<<" + "<<fabs(temp_coeffC.at(j__));
				}

				for(size_t k_=0; k_<temp_nthC.at(j__).size(); k_++)
				{
					//	cout<<" f_{";
					cout<<" f";
					temp_nthC.at(j__).at(k_).printD();
					//cout<<"}";
					if(temp_nthC.at(j__).at(k_).getE()>1)cout<<"^"<<temp_nthC.at(j__).at(k_).getE();
				}
			}


/*
			for(int i=0;i<this->size(); i++) {

				//cout<<"(";
				bool TEMP=1;

				if(this->getcmn(i)==1){ if(i!=0) cout<<" + ";
				}else if(this->getcmn(i)==-1){ cout<<" - ";
				}else if(this->getcmn(i)<0){
					cout<<" - "<<fabs(this->getcmn(i));
				}else{
					cout<<" + "<<fabs(this->getcmn(i));
				}

				cout<<" ";

				for(size_t j=0; j<this->norder(); j++) 
				{
					int ADDPOW=0;
					if(this->getfpow(j,i)==0) continue;

					FactVec temp(this->norder());
					temp.setD(j,1);

					FactVec temp2(this->norder());
					for(size_t k=0; k<this->norder(); k++) temp2.setD(k,this->getcfmn(i).getD(k));

					cout<<"f";
					temp.printD();

					ADDPOW += this->getfpow(j,i);

					if(temp == temp2) {ADDPOW += 1; TEMP=0;}

					if(ADDPOW>1)cout<<"^"<<ADDPOW;
					cout<<" ";
				}

				if(i<this->size()-1 && TEMP==1)
				{
					cout<<"f";
					this->getcfmn(i).printD();
				}
				//cout<<")";

			}
*/

			if(FORMAT!="NOEXPS")cout<<endl;
			if(FORMAT!="NOEXPS")cout<<endl;
		}


		void printD(const vector<int> abinN__, string FORMAT=""){

			if(abinN__.size()!=this->norder()){
				this->printD();
			}else{
				vector<int> mthN_;
				vector<int> abinN_;

				FactMomN * fact = new FactMomN();

				if(FORMAT!="NOEXPS"){
					cout<<endl;
					cout<<endl;
					cout<<"C_";
					for(size_t k_=0; k_<ndim.size(); k_++) cout<<ndim.at(k_);
					cout<<" = ";
				}

				vector<vector<FactVec> > temp_nthC;
				vector<long > temp_coeffC;
				this->Squeeze(temp_nthC,temp_coeffC);

				for(size_t j__=0; j__<temp_nthC.size(); j__++)
				{
					if(temp_coeffC.at(j__)==1){ if(j__!=0) cout<<" + ";
					}else if(temp_coeffC.at(j__)==-1){ cout<<" - ";
					}else if(temp_coeffC.at(j__)<0){
						cout<<" - "<<fabs(temp_coeffC.at(j__));
					}else{
						cout<<" + "<<fabs(temp_coeffC.at(j__));
					}

					for(size_t k_=0; k_<temp_nthC.at(j__).size(); k_++)
					{
						//	cout<<" f_{";
						//cout<<" f";
						//temp_nthC.at(j__).at(k_).printD();

						cout<<"(";
						mthN_.clear(); abinN_.clear();

						FactVec __temp(temp_nthC.at(j__).at(k_));

						for(int i__=0; i__<__temp.size(); i__++) {
							mthN_.push_back(__temp.getD(i__)); abinN_.push_back(abinN__.at(i__));
						}

						fact->SetFmnabN(mthN_,abinN_);
						fact->printD("eff");
						cout<<")";

						//cout<<"}";
						if(temp_nthC.at(j__).at(k_).getE()>1)cout<<"^"<<temp_nthC.at(j__).at(k_).getE();
					}
				}



				/*
				   for(int i=0; i<this->size(); i++) {

				   bool TEMP=1;

				   if(this->getcmn(i)==1){ if(i!=0) cout<<" + ";
				   }else if(this->getcmn(i)==-1){ cout<<" - ";
				   }else if(this->getcmn(i)<0){
				   cout<<" - "<<fabs(this->getcmn(i));
				   }else{
				   cout<<" + "<<fabs(this->getcmn(i));
				   }

				   cout<<" ";

				   for(size_t j=0; j<this->norder(); j++) 
				   {
				   int ADDPOW=0;
				   if(this->getfpow(j,i)==0) continue;

				   FactVec temp(this->norder());
				   temp.setD(j,1);

				   FactVec temp2(this->norder());
				   for(size_t k=0; k<this->norder(); k++) temp2.setD(k,this->getcfmn(i).getD(k));

				//cout<<"f";
				//temp.printD();

				cout<<"(";
				mthN_.clear(); abinN_.clear();

				for(int i__=0; i__<temp.size(); i__++) {
				mthN_.push_back(temp.getD(i__)); abinN_.push_back(abinN__.at(i__));
				}

				fact->SetFmnabN(mthN_,abinN_);
				fact->printD();
				cout<<")";

				ADDPOW += this->getfpow(j,i);

				if(temp == temp2) {ADDPOW += 1; TEMP=0;}

				if(ADDPOW>1)cout<<"^"<<ADDPOW;
				cout<<" ";
				}

				if(i<this->size()-1 && TEMP==1)
				{
				//cout<<"f";
				//this->getcfmn(i).printD();
				cout<<"(";
				mthN_.clear(); abinN_.clear();
				for(int i__=0; i__<this->getcfmn(i).size(); i__++) {
				mthN_.push_back(this->getcfmn(i).getD(i__)); abinN_.push_back(abinN__.at(i__));
				}
				fact->SetFmnabN(mthN_,abinN_);
				fact->printD();
				cout<<")";
				} 
				}


				 */
			}
		}




		void CalcCentMomentN(const vector<int> Np, const vector<double> Npeff, long long& ObsMom_, long double& IncMom_)
		{

			size_t isize =0;
			for(size_t k_=0; k_<abinNN.size(); k_++) {isize += abinNN.at(k_);}

			if(Np.size()!=isize){
				cout<<"ERROR: Np.size() != CentVecN::norder()"<<endl;
			}else if(Npeff.size()!=isize){
				cout<<"ERROR: Npeff.size() != CentVecN::norder()"<<endl;
			}else{

				FactMomN * fact = new FactMomN();

				vector<long double> IncF00__;
				vector<long long> ObsF00__;

				   for(size_t j=0; j<this->norder(); j++) 
				   {
					   FactVec temp(this->norder());
					   temp.setD(j,1);
					   mthNN.clear();

					   for(int i__=0; i__<temp.size(); i__++) {
						   mthNN.push_back(temp.getD(i__));
					   }

					   fact->SetFmnabN(mthNN,abinNN);
					   fact->GetMomentN(Np,Npeff, ObsMom_, IncMom_);
					   IncF00__.push_back(IncMom_);
					   ObsF00__.push_back(ObsMom_);

					   cout<<" |"<<ObsMom_<<","<<IncMom_<<"| ";
				   }

				   long double _IncMom__ = 0.0L;
				   long long _ObsMom__ = 0LL;

				   for(int i=0; i<this->size(); i++)
				   {
					long double IncMom__ = 1.0L;
					long long ObsMom__ = 1LL;

					bool TEMP=1;


					if(this->getcmn(i)==1){ 
						if(i!=0) cout<<" + ";
					}else if(this->getcmn(i)==-1){ cout<<" - ";
					}else if(this->getcmn(i)<0){
						cout<<" - "<<fabs(this->getcmn(i));
					}else{
						cout<<" + "<<fabs(this->getcmn(i));
					}
					cout<<" ";
					
					for(size_t j=0; j<this->norder(); j++) 
					{

						int ADDPOW=0;
						FactVec temp(this->norder());
						temp.setD(j,1);

						FactVec temp2(this->norder());
						for(size_t k=0; k<this->norder(); k++) temp2.setD(k,this->getcfmn(i).getD(k));

						ADDPOW += this->getfpow(j,i);

						if(temp == temp2) {ADDPOW += 1; TEMP=0;}

						IncMom__ *= pow(IncF00__.at(j),ADDPOW);	
						ObsMom__ *= pow(ObsF00__.at(j),ADDPOW);
						
						if(j<this->norder()-1) cout<<" {"<<IncMom__<<"} ";
					}

					if(i<this->size()-1 && TEMP==1)
					{
						mthNN.clear();
						for(int i__=0; i__<this->getcfmn(i).size(); i__++) {
							mthNN.push_back(this->getcfmn(i).getD(i__));
						}
						fact->SetFmnabN(mthNN,abinNN);
						fact->GetMomentN(Np,Npeff, ObsMom_, IncMom_);
						cout<<" ["<<ObsMom_<<"] ";

						ObsMom__ *= ObsMom_; 
						IncMom__ *= IncMom_;
					} 

					IncMom__ *= (this->getcmn(i));
					ObsMom__ *= (this->getcmn(i));

					_IncMom__ += IncMom__;
					_ObsMom__ += ObsMom__;
				   }


				   ObsMom_ = _ObsMom__; IncMom_ = _IncMom__;

				   cout<<"("<<IncMom_<<")"<<endl;

				   vector<long double>().swap(IncF00__);
				   vector<long long>().swap(ObsF00__);
			}
		}


//		void Squeeze()//vector<vector<FactVec> >& __temp_nthC, vector<long >& __temp_coeffC)
		void Squeeze(vector<vector<FactVec> >& __temp_nthC, vector<long >& __temp_coeffC)
                {
//			cout<<"Inside Squeeze: "<<this->size()<<" "<<this->getCCmn().size()<<" "<<this->getcfmn().size()<<endl;
			vector<vector<FactVec> > temp_nthC;
			vector<long > temp_coeffC;

			for(int i=0;i<this->size();i++){
				bool TEMP=1;

				//if(this->getcmn(i)==1){ if(i!=0) cout<<" + ";
				//}else if(this->getcmn(i)==-1){ cout<<" - ";
				//}else if(this->getcmn(i)<0){
				//	cout<<" - "<<fabs(this->getcmn(i));
				//}else{
				//	cout<<" + "<<fabs(this->getcmn(i));
				//}

				//cout<<" ";
				temp_coeffC.push_back(this->getcmn(i));

				vector<FactVec>  temp_temp_nthC;

				for(size_t j=0; j<this->norder(); j++) 
				{
					int ADDPOW=0;
					if(this->getfpow(j,i)==0) continue;

					FactVec temp(this->norder());
					temp.setD(j,1);

					FactVec temp2(this->norder());
					for(size_t k=0; k<this->norder(); k++) temp2.setD(k,this->getcfmn(i).getD(k));

					//cout<<"f";
					//temp.printD();

					ADDPOW += this->getfpow(j,i);

					if(temp == temp2) {ADDPOW += 1; TEMP=0;}

					//if(ADDPOW>1)cout<<"^"<<ADDPOW;
					//cout<<" ";

					FactVec temp_temp_temp_nthC(temp);
					temp_temp_temp_nthC.setE(ADDPOW);
					temp_temp_nthC.push_back(temp_temp_temp_nthC);

				}


				if(i<this->size()-1 && TEMP==1)
				{
					//cout<<"f";
					//this->getcfmn(i).printD();

					temp_temp_nthC.push_back(this->getcfmn(i));

				}

				temp_nthC.push_back(temp_temp_nthC);
				vector<FactVec>().swap(temp_temp_nthC);

			}

//			cout<<endl;

			//re-arrange
			for(size_t j__=0; j__<temp_nthC.size(); j__++)
			{

				vector<FactVec> fvv_;

				for(size_t k_=0; k_<temp_nthC.at(j__).size(); k_++)
				{
					//FactVec temp_fvv_(temp_nthC.at(j__).at(k_));
					//temp_fvv_.setE(temp_nthC.at(j__).at(k_).expo());
					//fvv_.push_back(temp_fvv_);//temp_nthC.at(j__).at(k_));
					fvv_.push_back(temp_nthC.at(j__).at(k_));
				}

				std::sort(fvv_.begin(),fvv_.end(),sortsumD);

				temp_nthC.at(j__).swap(fvv_);

			}

	//		cout<<endl;

			//add-exponents 
			//for(size_t j__=0; j__<temp_nthC.size(); j__++)
			//{
			//	vector<FactVec> fvv_;
			//	for(size_t k_=0; k_<temp_nthC.at(j__).size(); k_++) fvv_.push_back(temp_nthC.at(j__).at(k_));

			//	for(size_t i__=0; i__<fvv_.size(); i__++)
			//	{
			//		if(fvv_.at(i__).expo()==0) continue;
			//		for(size_t k__=i__+1; k__<fvv_.size(); k__++)
			//		{
			//			if(fvv_.at(i__) == fvv_.at(k__))
			//			{
			//				fvv_.at(i__).addE(fvv_.at(k__).expo());
			//				fvv_.at(k__).setE(0);
			//			}
			//		}
			//	}

			//	vector<FactVec> fvv__;
			//	for(size_t i__=0; i__<fvv_.size(); i__++) if(fvv_.at(i__).expo()!=0)fvv__.push_back(fvv_.at(i__));
			//	temp_nthC.at(j__).swap(fvv__);

			//	vector<FactVec>().swap(fvv_);
			//}


			vector<FactVec> __fvv;
			for(size_t i__=0; i__<temp_nthC.size(); i__++)
			{
				vector<int> __v;
				for(size_t _j=0; _j<temp_nthC.at(i__).size(); _j++)
				{
					for(int _k=0; _k<temp_nthC.at(i__).at(_j).size(); _k++)
						__v.push_back(temp_nthC.at(i__).at(_j).getD(_k));

					__v.push_back(temp_nthC.at(i__).at(_j).getE());
				}

				FactVec temp__fvv;
				temp__fvv.setV(__v);
				__fvv.push_back(temp__fvv);

			}

		//	cout<<endl;
		//	for(size_t j__=0; j__<__fvv.size(); j__++){cout<<" "; __fvv.at(j__).printD(); cout<<" ";}
		//	cout<<endl;

			for(size_t j__=0; j__<__fvv.size(); j__++)
			{
				for(size_t i__=j__+1; i__<__fvv.size(); i__++)
				{
					if(__fvv.at(j__)==__fvv.at(i__)){
						temp_coeffC.at(j__) += temp_coeffC.at(i__);
						temp_coeffC.at(i__) *=0;
						temp_nthC.at(i__).at(0) *=0;
					}
				}
			}
			vector<FactVec>().swap(__fvv);

			vector<vector<FactVec> > temp_nthC__;
			vector<long > temp_coeffC__;

			for(size_t j__=0; j__<temp_nthC.size(); j__++) 
			{
				//bool check=1;

				//for(size_t k_=0; k_<temp_nthC.at(j__).size(); k_++)
				//{
				//	if(temp_nthC.at(j__).at(k_).sumD() <=1) check*=0;
				//}

				if(temp_coeffC.at(j__)!=0){
					temp_nthC__.push_back(temp_nthC.at(j__));
					temp_coeffC__.push_back(temp_coeffC.at(j__));
				}
			}

			temp_nthC.swap(temp_nthC__); temp_coeffC.swap(temp_coeffC__);

/*
			for(size_t j__=0; j__<temp_nthC.size(); j__++)
			{
				if(temp_coeffC.at(j__)==1){ if(j__!=0) cout<<" + ";
				}else if(temp_coeffC.at(j__)==-1){ cout<<" - ";
				}else if(temp_coeffC.at(j__)<0){
					cout<<" - "<<fabs(temp_coeffC.at(j__));
				}else{
					cout<<" + "<<fabs(temp_coeffC.at(j__));
				}

				for(size_t k_=0; k_<temp_nthC.at(j__).size(); k_++)
				{
					//	cout<<" f_{";
					cout<<" f";
					temp_nthC.at(j__).at(k_).printD();
					//cout<<"}";
					if(temp_nthC.at(j__).at(k_).getE()>1)cout<<"^"<<temp_nthC.at(j__).at(k_).getE();
				}
			}
*/
			//cout<<endl;


			__temp_nthC.swap(temp_nthC);
			__temp_coeffC.swap(temp_coeffC);


		}



};



long CentVecN::calccomb(int Nval_, int nval_)
{
	long NUM_FACT=1;
	long DEN_FACT=1;
	long comb;

	for(int i=0; i<nval_; i++) NUM_FACT *= (Nval_-i);
	for(int j=1; j<=nval_; j++) DEN_FACT *= j;

	comb = NUM_FACT/DEN_FACT;

	return comb;
}	


void CentVecN::Corr2CentN(const vector<int> mth__,string FORMAT, string FORMAT2)
{
	size_t isize=mth__.size()*2;

	vector<vector<FactVec> > Cmn;
	vector<vector<int> > vcmn;

	Cmn.resize(mth__.size());
	vcmn.resize(mth__.size());

	int ipos=0;

	for(size_t k_=0; k_<vcmn.size(); k_++)
	{
		for(int i_=0; i_<=mth__.at(k_); i_++)
		{
			vcmn.at(k_).push_back(this->calccomb(mth__.at(k_),i_)*pow(-1,i_));

			FactVec TEMP_Cmn(isize);

			TEMP_Cmn.setD(0+(k_*2),mth__.at(k_)-i_);
			TEMP_Cmn.setD(1+(k_*2),i_);

			Cmn.at(k_).push_back(TEMP_Cmn);
			ipos++;
		}
	}



	vector<int> vcmn2;
	vector<FactVec> Cmn2;

	for(size_t i__=1; i__<vcmn.size(); i__++){
		for(size_t j__=0; j__<vcmn.at(i__).size(); j__++)
		{
			for(size_t k__=0; k__<vcmn.at(0).size(); k__++)
			{
				vcmn2.push_back(vcmn.at(i__).at(j__)*vcmn.at(0).at(k__));
				FactVec TEMP_Cmn_(Cmn.at(0).at(k__));
				TEMP_Cmn_ += Cmn.at(i__).at(j__);
				Cmn2.push_back(TEMP_Cmn_);
			}

		}

		vcmn.at(0) = vcmn2;
		Cmn.at(0) = Cmn2;
		vcmn2.clear();

		Cmn2.clear();
	}

	CCmn = Cmn.at(0);
	vvmn = vcmn.at(0);


	if(FORMAT!="NOEXPS"){
		if(FORMAT2!="NOEXPS"){
			cout<<"Ç_";
			for(size_t k_=0; k_<mth__.size(); k_++) cout<<mth__.at(k_);
			cout<<"=";
		}
		for(size_t l__=0; l__<vvmn.size(); l__++){ 

			if(vvmn.at(l__)==1){ if(l__!=0) cout<<" + ";
			}else if(vvmn.at(l__)==-1){ cout<<" - ";
			}else if(vvmn.at(l__)<0){
				cout<<" - "<<fabs(vvmn.at(l__));
			}else{
				cout<<" + "<<fabs(vvmn.at(l__));
			}

			cout<<" C_"; 
			CCmn.at(l__).printD();
		}
		if(FORMAT2!="NOEXPS")
			cout<<endl;
	}

}




void CentVecN::Corr2FactN(const vector<int> mth__, string FORMAT, string VERBOSE_)
{

	CCmn.clear();
	vvmn.clear();

	this->Corr2CentN(mth__,FORMAT);

	if(FORMAT!="NOEXPS")	
	{
		cout<<"Ç_";
		for(size_t k_=0; k_<mth__.size(); k_++) cout<<mth__.at(k_);
		cout<<"=";
	}

	for(size_t l__=0; l__<vvmn.size(); l__++){ 

		if(VERBOSE_=="LOUD")
		{
			if(vvmn.at(l__)==1){ if(l__!=0) cout<<" + ";
			}else if(vvmn.at(l__)==-1){ cout<<" - ";
			}else if(vvmn.at(l__)<0){
				cout<<" - "<<fabs(vvmn.at(l__));
			}else{
				cout<<" + "<<fabs(vvmn.at(l__));
			}

			cout<<" ("; 

		}
		vector<int> norder__;
		for(int i__=0; i__<CCmn.at(l__).size(); i__++){
			norder__.push_back(CCmn.at(l__).getD(i__));
		}

		CentVecN * cnt__ =new CentVecN();
		cnt__->Cent2FactN(norder__);
		if(VERBOSE_=="LOUD")cnt__->printD("NOEXPS");
		if(VERBOSE_=="LOUD")cout<<") ";
	}

	if(FORMAT!="NOEXPS")	cout<<endl;
}

void CentVecN::Corr2FactN(const vector<int> mth__, const vector<int> abin__, string FORMAT, string VERBOSE_)
{

	if(abin__.size()==0 || abin__.size()==1)
	{
		this->Corr2FactN(mth__,FORMAT);
	}else{
		CCmn.clear();
		vvmn.clear();

		this->Corr2CentN(mth__,FORMAT);

		if(FORMAT!="NOEXPS"){
			cout<<"Ç_";
			for(size_t k_=0; k_<mth__.size(); k_++) cout<<mth__.at(k_);
			cout<<"=";
		}

		for(size_t l__=0; l__<vvmn.size(); l__++){ 

			if(VERBOSE_=="LOUD"){
				if(vvmn.at(l__)==1){ if(l__!=0) cout<<" + ";
				}else if(vvmn.at(l__)==-1){ cout<<" - ";
				}else if(vvmn.at(l__)<0){
					cout<<" - "<<fabs(vvmn.at(l__));
				}else{
					cout<<" + "<<fabs(vvmn.at(l__));
				}

				cout<<" ("; 
			}

			vector<int> norder__;
			for(int i__=0; i__<CCmn.at(l__).size(); i__++){
				norder__.push_back(CCmn.at(l__).getD(i__));
			}

			CentVecN * cnt__ =new CentVecN();
			cnt__->Cent2FactN(norder__);
			if(VERBOSE_=="LOUD")	cnt__->printD(abin__,"NOEXPS");
			if(VERBOSE_=="LOUD")		cout<<") ";
		}

		if(FORMAT!="NOEXPS") cout<<endl;
	}
}



void CentVecN::GetCentMomentN(const vector<int> mth__, const vector<int> abin__, const vector<int> Np, const vector<double> Npeff, long long& ObsMom_, long double& IncMom_, string FORMAT)
{

	long double _IncMom__ = 0.0L;
	long long _ObsMom__ = 0LL;

	if(abin__.size()==0 || abin__.size()==1)
	{
		this->Corr2FactN(mth__,FORMAT);
	}else{
		CCmn.clear();
		vvmn.clear();

		this->Corr2CentN(mth__,FORMAT);

		cout<<"Ç_";
		for(size_t k_=0; k_<mth__.size(); k_++) cout<<mth__.at(k_);
		cout<<"=";


		for(size_t l__=0; l__<vvmn.size(); l__++){ 

			long double IncMom__ = 1.0L;
			long long ObsMom__ = 1LL;

			if(vvmn.at(l__)==1){ if(l__!=0) cout<<" + ";
			}else if(vvmn.at(l__)==-1){ cout<<" - ";
			}else if(vvmn.at(l__)<0){
				cout<<" - "<<fabs(vvmn.at(l__));
			}else{
				cout<<" + "<<fabs(vvmn.at(l__));
			}

			cout<<" ("; 

			IncMom__ *= (vvmn.at(l__));
			ObsMom__ *= (vvmn.at(l__));

			vector<int> norder__;
			for(int i__=0; i__<CCmn.at(l__).size(); i__++){
				norder__.push_back(CCmn.at(l__).getD(i__));
			}

			CentVecN * cnt__ =new CentVecN();
			cnt__->Cent2FactN(norder__);

			//	cnt__->printD(abin__,"NOEXPS");
			//	long double IncMom = 0.0L;
			//	long long ObsMom = 0LL;
			//
			cnt__->SetFmnabN(abin__);
			//	cnt__->CalcCentMomentN();//Np,eff,ObsMom,IncMom);
			cnt__->CalcCentMomentN(Np, Npeff, ObsMom_, IncMom_);
			//
			IncMom__ *= IncMom_;
			ObsMom__ *= ObsMom_;

			_IncMom__ += IncMom__;
			_ObsMom__ += ObsMom__;

			cout<<") ";
		}

		cout<<endl;
	}

	IncMom_ = _IncMom__; 
	ObsMom_ = _ObsMom__;

	cout<<"("<<_IncMom__<<","<<_ObsMom__<<")"<<endl;
}




vector<long> CentVecN::Stirling(const int order){

	//Function to return Stirling series for a given order

	vector<long> anm1;

	if(order==0){
		anm1.push_back(1);
	}else{

		long Ninit=1;
		anm1.push_back(Ninit);

		for(int i=1; i<=order; i++)
		{

			vector<long> an;

			for(int k=0;k<i;k++) an.push_back(0);

			for(int j=0; j<i; j++) 
			{
				long an_m=0, an_m_1=0;


				size_t dim1=i-j-1;
				size_t dim2=i-j-1-1;

				if(anm1.size()>dim1) { an_m=anm1.at(dim1);}
				if(anm1.size()>dim2) { an_m_1=anm1.at(dim2);}

				long new_an = (i-j)*an_m + an_m_1;

				an.at(i-j-1)=new_an;
			}

			anm1.clear();
			anm1=an;
		}

	}
	return anm1;
}



CentVec * CentVecN::Cent2Fact(const int order)
{
	CentVec *cmnab = new CentVec(order);

	for(int j=0;j<=order;j++){

		vector<long> amn;
		amn=Stirling(order-j);
		int stsize=amn.size();

		for(int i=0;i<stsize;i++)
		{       
			cmnab->setvec(j,i,int(pow(-1,j)*calccomb(order,j)*amn.at(i)));
		}

	}

	return cmnab;
}


void CentVecN::Cent2FactN(const vector<int> morder)
{

	this->clear();

	mdim=morder.size();
	fpow.resize(morder.size());

	CentVec * cmn__[morder.size()];

	for(size_t i_=0; i_<morder.size(); i_++) cmn__[i_] = Cent2Fact(morder.at(i_));

	int isize=0;

	if(morder.size()==1)
	{
		for(int i=0; i<cmn__[0]->size(); i++)
		{
			vector<int> fpow_;
			vector<FactVec> cfmn_;

			FactVec TEMP_CFMN(morder.size());
			TEMP_CFMN.setD(0,cmn__[0]->getcfmn(i).getD(0));

			fpow_.push_back(cmn__[0]->getf10pow(i));
			cfmn_.push_back(TEMP_CFMN);
			long cmnval_ = cmn__[0]->getcmn(i);
			this->setvec(fpow_,cfmn_,cmnval_);
		}

	}else{
		for(int i=0; i<cmn__[0]->size(); i++)
		{
			for(int k=0; k<cmn__[1]->size(); k++) 
			{
				vector<int> fpow_;
				vector<FactVec> cfmn_;

				FactVec TEMP_CFMN(morder.size());
				TEMP_CFMN.setD(0,cmn__[0]->getcfmn(i).getD(0));

				fpow_.push_back(cmn__[0]->getf10pow(i));
				cfmn_.push_back(TEMP_CFMN);
				long cmnval_ = cmn__[0]->getcmn(i);


				FactVec TEMP_CFMN_(morder.size());
				TEMP_CFMN_.setD(1,cmn__[1]->getcfmn(k).getD(0));

				fpow_.push_back(cmn__[1]->getf10pow(k));
				cfmn_.push_back(TEMP_CFMN_);

				cmnval_ *= cmn__[1]->getcmn(k);

				this->setvec(fpow_,cfmn_,cmnval_);

				isize++;
			}
		}


		for(size_t l__=2; l__<morder.size(); l__++)
		{

			CentVecN * cmnab_ = new CentVecN(morder.size());


			for(int i=0; i<this->size(); i++)
			{
				for(int k=0; k<cmn__[l__]->size(); k++) 
				{
					vector<int> fpow_;
					vector<FactVec> cfmn_;
					long cmnval_ =1L;

					cmnval_ *= this->getcmn(i);

					for(size_t j=0;j<this->norder();j++)
					{
						if(this->getfpow(j).size()!=0) fpow_.push_back(this->getfpow(j,i));
					}

					cfmn_.push_back(this->getcfmn(i));


					FactVec TEMP_CFMN__(morder.size());
					TEMP_CFMN__.setD(l__,cmn__[l__]->getcfmn(k).getD(0));

					fpow_.push_back(cmn__[l__]->getf10pow(k));
					cfmn_.push_back(TEMP_CFMN__);

					cmnval_ *= cmn__[l__]->getcmn(k);

					cmnab_->setvec(fpow_,cfmn_,cmnval_);
				}
			}
			*this = *cmnab_;
		}
	}
	this->setndim(morder);
//	this->Squeeze();
}



class MomErrorN
{

	protected:
		vector<int> mthE;
		vector<vector<FactVec> > nthE;
		vector<int> coeff;
		size_t fsize;
	public:
		MomErrorN(const vector<int> mth_)
		{
			mthE=mth_;
			fsize=(int(mthE.size()/2));
		}	

		MomErrorN()
		{
			fsize=0;
		}

		~MomErrorN()
		{
			vector<vector<FactVec> >().swap(nthE);
			vector<int>().swap(coeff);
		}

		void setmthE(const vector<int> mth_)
		{
			mthE=mth_;
			fsize=(int(mthE.size()/2));
		}

		vector<vector<FactVec> > getnthE(){return nthE;};
		vector<int> getcoeff(){return coeff;}

		size_t size() const {return fsize;}

		void clear()
		{
			mthE.clear();
			nthE.clear();
			coeff.clear();
			fsize=0;
		}


		void CalcMomError(const vector<int> mth_, string FORMAT="LOUD") //
		{
			this->clear();
			mthE=mth_;
			fsize=(int(mthE.size()/2));
			CalcMomError(FORMAT); //
		}

		void CalcMomError(string FORMAT="LOUD") //
		{

			//Following are examples of the expressions for first few orders of moments
			////
			//(order-1) 5-terms:
			//∑_{i | j} = Cov(µ_{i}, µ_{j}) = 
			//
			//μ_{i+j} − μ_{i} μ_{j} 
			//
			//
			//(order-2) 10-terms:
			//∑_{i,j | k,l} = Cov(µ_{i,j}, µ_{k,l}) = 
			//
			//μ_{i+k,j+l} − μ_{i,j} μ_{k,l} 
			//
			//
			//(order-3) 17-terms:
			//∑_{i,j,k | l,m,n} = Cov(µ_{i,j,k}, µ_{l,m,n}) = 
			//
			//μ_{i+l,j+m,k+n} − μ_{i,j,k} μ_{l,m,n}
			//


			vector<FactVec> TEMP_nthE;
			FactVec TEMP_TEMP_nthE (static_cast<int>(fsize));

			//
			//step-I (symmetric terms-1, total 1) ;
			//
			for(size_t j_=0; j_<fsize; j_++)
			{
				int index_=0;
				for(size_t k_=j_; k_<mthE.size(); k_+=fsize) index_ +=mthE.at(k_); 
				TEMP_TEMP_nthE.setD(j_,index_);
			}

			TEMP_nthE.push_back(TEMP_TEMP_nthE);  
			nthE.push_back(TEMP_nthE); coeff.push_back(1);


			//step-II (symmetric terms-2, total 1) ;
			//
			TEMP_nthE.clear();

			for(size_t j_=0; j_<fsize; j_++) TEMP_TEMP_nthE.setD(j_,mthE.at(j_));
			TEMP_nthE.push_back(TEMP_TEMP_nthE);  


			for(size_t j_=fsize; j_<2*fsize; j_++) TEMP_TEMP_nthE.setD((j_-fsize),mthE.at(j_));
			TEMP_nthE.push_back(TEMP_TEMP_nthE);  

			nthE.push_back(TEMP_nthE); coeff.push_back(-1);


			//

			if(FORMAT=="LOUD")
			{
				cout<<endl;
				cout<<endl;
				cout<<endl;
				cout<<endl;
				cout<<endl;
				cout<<"∑_{";
				nthE.at(1).at(0).printD();
				cout<<",";
				nthE.at(1).at(1).printD();
				cout<<"}=";


				for(size_t j__=0; j__<nthE.size(); j__++)
				{
					bool check=1;

					for(size_t k_=0; k_<nthE.at(j__).size(); k_++)
					{
						if(nthE.at(j__).at(k_).sumD() <=1) check*=0;
					}

					if(check){
						if(coeff.at(j__)==1){ if(j__!=0) cout<<" + ";
						}else if(coeff.at(j__)==-1){ cout<<" - ";
						}else if(coeff.at(j__)<0){
							cout<<" - "<<fabs(coeff.at(j__));
						}else{
							cout<<" + "<<fabs(coeff.at(j__));
						}

						for(size_t k_=0; k_<nthE.at(j__).size(); k_++)
						{
							cout<<" µ_{";
							nthE.at(j__).at(k_).printD();
							cout<<"} ";
						}
					}
				}

				cout<<endl;
				cout<<endl;
				cout<<endl;
				cout<<endl;
			}

		}


};

class MomVecN : public CentVecN, MomErrorN 
{
	private:
		size_t mdim;
		int nsize;
		bool CALCHECK;
		vector<int> morder;
	public:
		MomVecN(const size_t mdim_=0) : CentVecN(0,"Moment") , MomErrorN()
	{
		mdim=mdim_;
		nsize=0;
		CALCHECK=0;
	}

		~MomVecN()
		{
			vector<int>().swap(morder);
		}

		void clear(){
			mdim=0;
			nsize=0;
			CALCHECK=0;
			morder.clear();
		}

		void CalcVariance(const vector<int> morder_)
		{

			vector<int> merr_;
			for(size_t i_=0; i_<morder_.size(); i_++) merr_.push_back(morder_.at(i_));
			for(size_t i_=0; i_<morder_.size(); i_++) merr_.push_back(morder_.at(i_));

			CalcMomError(merr_);
		}


		void Error(const vector<int> morder_)
		{

			vector<int> merr_;
			for(size_t i_=0; i_<morder_.size(); i_++) merr_.push_back(morder_.at(i_));
			for(size_t i_=0; i_<morder_.size(); i_++) merr_.push_back(morder_.at(i_));

			CalcMomError(merr_);
		}

		void Mom2FactN(const vector<int> morder_, string FORMAT="")
		{

			cout<<"Func to convert moments µ_m =<N^m> to factorial moments f_m = <N(N-1)..(N-m+1)> "<<endl; 

			this->clear();
			morder=morder_;
			mdim=morder.size();

			vector<vector<FactVec> > CCmn_;
			vector<vector<int> > vvmn_;

			vector<FactVec> temp_CCmn_;
			vector<int> temp_vvmn_;

			for(size_t i_=0; i_<morder.size(); i_++)
			{
				temp_CCmn_.clear();
				temp_vvmn_.clear();

				vector<long> amn =Stirling(morder.at(i_));

				for(size_t j_=0; j_<amn.size(); j_++)
				{       
					FactVec temp_amn(static_cast<int>(mdim));

					if(morder.at(i_)>0){	
						temp_amn.setD(i_,j_+1);
					}else{
						temp_amn.setD(i_,0);
					}

					temp_CCmn_.push_back(temp_amn);
					temp_vvmn_.push_back(amn.at(j_));
				}

				CCmn_.push_back(temp_CCmn_);
				vvmn_.push_back(temp_vvmn_);
			}

			temp_CCmn_.clear();
			temp_vvmn_.clear();

			for(size_t i__=1; i__<morder.size(); i__++)
			{
				for(size_t j__=0; j__<vvmn_.at(0).size(); j__++){
					for(size_t k__=0; k__<vvmn_.at(i__).size(); k__++){

						int temp_c =1;
						temp_c *=  vvmn_.at(0).at(j__);
						temp_c *=  vvmn_.at(i__).at(k__); 

						FactVec temp_f(static_cast<int>(mdim));
						temp_f += CCmn_.at(0).at(j__);
						temp_f += CCmn_.at(i__).at(k__);

						temp_vvmn_.push_back(temp_c);
						temp_CCmn_.push_back(temp_f);
					}
				}

				CCmn_.at(0) = temp_CCmn_;
				vvmn_.at(0) = temp_vvmn_;

				temp_CCmn_.clear();
				temp_vvmn_.clear();
			}

			temp_CCmn_ = CCmn_.at(0);
			temp_vvmn_ = vvmn_.at(0);


			CCmn =temp_CCmn_;
			vvmn =temp_vvmn_;

			vector<vector<FactVec> >().swap(CCmn_);
			vector<vector<int> >().swap(vvmn_);
			vector<FactVec>().swap(temp_CCmn_);
			vector<int>().swap(temp_vvmn_);

			CALCHECK=1;

			if(FORMAT!="NOEXPS") this->printD(FORMAT);
		}

		void printD(string FORMAT="")
		{

			if(FORMAT!="NOEXPS"){
				cout<<"µ_";
				for(size_t __i=0; __i<morder.size(); __i++) cout<<morder.at(__i);
				cout<<"=";
			}

			for(size_t k_=0;k_<vvmn.size();k_++)
			{
				if(vvmn.at(k_)==1){ if(k_!=0) cout<<" + ";
				}else if(vvmn.at(k_)==-1){ cout<<" - ";
				}else if(vvmn.at(k_)<0){
					cout<<" - "<<fabs(vvmn.at(k_));
				}else{
					cout<<" + "<<fabs(vvmn.at(k_));
				}

				if(FORMAT!="NOEXPS"){
					cout<<" f";
					CCmn.at(k_).printD();
				}else if(FORMAT=="latex"){
					cout<<"{$rm f}_{";
					CCmn.at(k_).printD();
					cout<<"}";
				}else{
					cout<<"f";
					CCmn.at(k_).printD();
				}
			}

			if(FORMAT!="NOEXPS")cout<<endl;
		}

		void Mom2FactN(const vector<int> morder_, const vector<int> abinN__, string FORMAT="")
		{
			this->Mom2FactN(morder_,"NOEXPS");
			if(FORMAT!="NOEXPS") this->printD(abinN__,FORMAT);
		}


		void printD(const vector<int> abinN__, string FORMAT="")
		{
			if(FORMAT!="NOEXPS"){
				cout<<"µ_";
				for(size_t __i=0; __i<morder.size(); __i++) cout<<morder.at(__i);
				cout<<"=";
			}

			for(size_t k_=0;k_<vvmn.size();k_++)
			{
				if(vvmn.at(k_)==1){ if(k_!=0) cout<<" + ";
				}else if(vvmn.at(k_)==-1){ cout<<" - ";
				}else if(vvmn.at(k_)<0){
					cout<<" - "<<fabs(vvmn.at(k_));
				}else{
					cout<<" + "<<fabs(vvmn.at(k_));
				}

				if(FORMAT!="NOEXPS"){
					cout<<" f";
					CCmn.at(k_).printD();
				}else if(FORMAT=="latex"){
					cout<<"{$rm f}_{";
					CCmn.at(k_).printD();
					cout<<"}";
				}else{
					vector<int> mthN__;
					for(int i__=0; i__<CCmn.at(k_).size(); i__++)
					{
						mthN__.push_back(CCmn.at(k_).getD(i__));
					}

					FactMomN  fact__;

					fact__.SetFmnabN(mthN__,abinN__);
					fact__.printD("eff");
					vector<int>().swap(mthN__);
				}
			}

			if(FORMAT!="NOEXPS")cout<<endl;
		}


		vector<FactVec> CalcFactVec(const vector<int> morder)
		{

			if(CALCHECK!=1) this->Mom2FactN(morder);

			vector<FactVec> fvv;

			for(size_t k_=0; k_<CCmn.size(); k_++)
			{
				fvv.push_back(CCmn.at(k_));
			}

			return fvv;
		}


};





class NudynVecN : public CentVecN 
{
	private:
		size_t mdim;
		int nsize;
	public:
		NudynVecN(const size_t mdim_=0) : CentVecN(0,"Nudyn")
	{
		mdim=mdim_;
		nsize=0;
	}

		~NudynVecN()
		{
		}


		void Nudyn2FactN(const int nth__)
		{
			vector<int> mth__;
			mth__.push_back(nth__);
			Corr2CentN(mth__,"NOEXPS");

			//FactVec F10(2);
			//F10.setD(0,1);
			//FactVec F01(2);
			//F01.setD(1,1);

			cout<<"Nu_{dyn}(";
			for(size_t k_=0; k_<mth__.size(); k_++) cout<<mth__.at(k_);
			cout<<")=";

			for(size_t l__=0; l__<vvmn.size(); l__++){ 

				if(vvmn.at(l__)==1){ if(l__!=0) cout<<" + ";
				}else if(vvmn.at(l__)==-1){ cout<<" - ";
				}else if(vvmn.at(l__)<0){
					cout<<" - "<<fabs(vvmn.at(l__));
				}else{
					cout<<" + "<<fabs(vvmn.at(l__));
				}

				cout<<" F_"; 
				CCmn.at(l__).printD();

				cout<<"/";
				cout<<"( ";
				if(CCmn.at(l__).getD(0) !=0) cout<<"F_10^"<<CCmn.at(l__).getD(0);
				cout<<" ";
				if(CCmn.at(l__).getD(1) !=0) cout<<"F_01^"<<CCmn.at(l__).getD(1);
				cout<<")";
			}
			cout<<endl;

		}

};



class CentErrorN : public CentVecN 
{

	protected:
		vector<int> mthE;
		vector<vector<FactVec> > nthE;
		vector<int> coeff;
		size_t fsize;
	public:
		CentErrorN(const vector<int> mth_) : CentVecN(0,"Error")
		{
			mthE=mth_;
			fsize=(int(mthE.size()/2));
		}

		~CentErrorN()
		{
			vector<vector<FactVec> >().swap(nthE);
			vector<int>().swap(coeff);
		}


		vector<vector<FactVec> > getnthE(){return nthE;};
		vector<int> getcoeff(){return coeff;}

		size_t size() const {return fsize;}

		void CalcCentError(string FORMAT="LOUD") //
		{

			//Following are examples of the expressions for first few orders of moments
			////
			//(order-1) 5-terms:
			//∑_{i | j} = Cov(µ_{i}, µ_{j}) = 
			//
			//μ_{i+j} − μ_{i} μ_{j} 
			//
			//− i μ_{i−1} μ_{j+1} − j μ_{i+1} μ_{j−1} 
			//
			//+ ij μ_{i−1} μ_{j−1} μ_{2} 
			//
			//
			//(order-2) 10-terms:
			//∑_{i,j | k,l} = Cov(µ_{i,j}, µ_{k,l}) = 
			//
			//μ_{i+k,j+l} − μ_{i,j} μ_{k,l} 
			//
			//− i μ_{i−1,j} μ_{k+1,l} − j μ_{i,j-1} μ_{k,l+1} - k µ_{i+1,j} µ_{k-1,l} - l µ_{i,j+1} µ_{k,l-1}  
			//
			//+ ik μ_{i−1,j} μ_{k−1,l} μ_{2,0} + jl µ_{i,j-1} µ_{k,l-1} µ_{0,2} 
			//+ il µ_{i-1,j} µ_{k,l-1} µ_{1,1} + jk µ_{i,j-1} µ_{k-1,l} µ_{1,1}
			//
			//
			//(order-3) 17-terms:
			//∑_{i,j,k | l,m,n} = Cov(µ_{i,j,k}, µ_{l,m,n}) = 
			//
			//μ_{i+l,j+m,k+n} − μ_{i,j,k} μ_{l,m,n}
			//
			//- i µ_{i-1,j,k} µ_{l+1,m,n} - j µ_{i,j-1,k} µ_{l,m+1,n} - k µ_{i,j,k-1} µ_{l,m,n+1} 
			//- l µ_{i+1,j,k} µ_{l-1,m,n} - m µ_{i,j+1,k} µ_{l,m-1,n} - n µ_{i,j,k+1} µ_{l,m,n-1}
			// 
			//+ il µ_{i-1,j,k} µ_{l-1,m,n} µ_{2,0,0} + jm µ_{i,j-1,k} µ_{l,m-1,n} µ_{0,2,0} 
			//+ kn µ_{i,j,k-1} µ_{l,m,n-1} µ_{0,0,2} + im µ_{i-1,j,k} µ_{l,m-1,n} µ_{1,1,0} 
			//+ in µ_{i-1,j,k} µ_{l,m,n-1} µ_{1,0,1} + jn µ_{i,j-1,k} µ_{l,m,n-1} µ_{0,1,1}
			//+ jl µ_{i,j-1,k} µ_{l-1,m,n} µ_{1,1,0} + kl µ_{i,j,k-1} µ_{l-1,m,n} µ_{1,0,1}
			//+ km µ_{i,j,k-1} µ_{i,m-1,n} µ_{0,1,1}
			//



			vector<FactVec> TEMP_nthE;
			FactVec TEMP_TEMP_nthE (static_cast<int>(fsize));

			//
			//step-I (symmetric terms-1, total 1) ;
			//
			for(size_t j_=0; j_<fsize; j_++)
			{
				int index_=0;
				for(size_t k_=j_; k_<mthE.size(); k_+=fsize) index_ +=mthE.at(k_); 
				TEMP_TEMP_nthE.setD(j_,index_);
			}

			TEMP_nthE.push_back(TEMP_TEMP_nthE);  
			nthE.push_back(TEMP_nthE); coeff.push_back(1);


			//step-II (symmetric terms-2, total 1) ;
			//
			TEMP_nthE.clear();

			for(size_t j_=0; j_<fsize; j_++) TEMP_TEMP_nthE.setD(j_,mthE.at(j_));
			TEMP_nthE.push_back(TEMP_TEMP_nthE);  


			for(size_t j_=fsize; j_<2*fsize; j_++) TEMP_TEMP_nthE.setD((j_-fsize),mthE.at(j_));
			TEMP_nthE.push_back(TEMP_TEMP_nthE);  

			nthE.push_back(TEMP_nthE); coeff.push_back(-1);


			//step-III (single asymmetric terms-3, total mthE.size()) ;
			//
			for(size_t i__=0; i__<fsize; i__++){

				TEMP_nthE.clear();
				for(size_t j_=0; j_<fsize; j_++)
				{
					int index_=mthE.at(j_);
					if(j_==i__) index_ -= 1;
					TEMP_TEMP_nthE.setD(j_,index_);
				}

				TEMP_nthE.push_back(TEMP_TEMP_nthE);

				for(size_t j_=fsize; j_<2*fsize; j_++)
				{
					int index_=mthE.at(j_);
					if(j_==i__+fsize) index_ += 1;
					TEMP_TEMP_nthE.setD(j_-fsize,index_);
				}

				TEMP_nthE.push_back(TEMP_TEMP_nthE);

				nthE.push_back(TEMP_nthE); coeff.push_back(-mthE.at(i__));
			}


			for(size_t i__=0; i__<fsize; i__++){
				TEMP_nthE.clear();
				for(size_t j_=0; j_<fsize; j_++)
				{
					int index_=mthE.at(j_);
					if(j_==i__) index_ += 1;
					TEMP_TEMP_nthE.setD(j_,index_);
				}

				TEMP_nthE.push_back(TEMP_TEMP_nthE);

				for(size_t j_=fsize; j_<2*fsize; j_++)
				{
					int index_=mthE.at(j_);
					if(j_==i__+fsize) index_ -= 1;
					TEMP_TEMP_nthE.setD(j_-fsize,index_);
				}

				TEMP_nthE.push_back(TEMP_TEMP_nthE);
				nthE.push_back(TEMP_nthE); coeff.push_back(-mthE.at(i__+fsize));
			}

			//step-IV (double asymmetric terms-4, total mthE.size()^2/4) ;
			//
			for(size_t i__=0; i__<fsize; i__++){
				for(size_t j__=0; j__<fsize; j__++){

					TEMP_nthE.clear();

					int index_=0;
					for(size_t k__=0; k__<fsize; k__++)
					{
						index_=mthE.at(k__);
						if(k__==i__) index_ -= 1;
						TEMP_TEMP_nthE.setD(k__,index_);
					}
					TEMP_nthE.push_back(TEMP_TEMP_nthE);

					index_=0;
					for(size_t k__=fsize; k__<2*fsize; k__++)
					{
						index_=mthE.at(k__);
						if(k__==j__+fsize) index_ -= 1;
						TEMP_TEMP_nthE.setD(k__-fsize,index_);
					}
					TEMP_nthE.push_back(TEMP_TEMP_nthE);

					TEMP_TEMP_nthE.clear();
					TEMP_TEMP_nthE.resize(static_cast<int>(fsize));

					for(size_t k__=0; k__<fsize; k__++)
					{
						if(k__==i__)TEMP_TEMP_nthE.addD(k__,1);
						if(k__==j__)TEMP_TEMP_nthE.addD(k__,1);
					}

					TEMP_nthE.push_back(TEMP_TEMP_nthE);

					coeff.push_back(mthE.at(i__)*mthE.at(j__+fsize));

					nthE.push_back(TEMP_nthE);

				}

			}

			//

			//this->Squeeze();


			if(FORMAT=="LOUD")
			{
				cout<<endl;
				cout<<endl;
				cout<<endl;
				cout<<endl;
				cout<<endl;
				cout<<"∑_{";
				nthE.at(1).at(0).printD();
				cout<<",";
				nthE.at(1).at(1).printD();
				cout<<"}=";


				for(size_t j__=0; j__<nthE.size(); j__++)
				{
					bool check=1;

					for(size_t k_=0; k_<nthE.at(j__).size(); k_++)
					{
						if(nthE.at(j__).at(k_).sumD() <=1) check*=0;
					}

					if(check){
						if(coeff.at(j__)==1){ if(j__!=0) cout<<" + ";
						}else if(coeff.at(j__)==-1){ cout<<" - ";
						}else if(coeff.at(j__)<0){
							cout<<" - "<<fabs(coeff.at(j__));
						}else{
							cout<<" + "<<fabs(coeff.at(j__));
						}

						for(size_t k_=0; k_<nthE.at(j__).size(); k_++)
						{
							cout<<" Ç_{";
							nthE.at(j__).at(k_).printD();
							cout<<"} ";
						}
					}
				}

				cout<<endl;
				cout<<endl;
				cout<<endl;
				cout<<endl;
			}

			this->Squeeze();


			nthErr=nthE;


		}


		//void printD()
		//{

		//	//					this->Squeeze();
		//	for(size_t j__=0; j__<nthE.size(); j__++)
		//	{
		//		bool check=1;

		//		for(size_t k_=0; k_<nthE.at(j__).size(); k_++)
		//		{
		//			if(nthE.at(j__).at(k_).sumD() <=1) check*=0;
		//		}

		//		if(check){
		//			if(coeff.at(j__)==1){ if(j__!=0) cout<<" + ";
		//			}else if(coeff.at(j__)==-1){ cout<<" - ";
		//			}else if(coeff.at(j__)<0){
		//				cout<<" - "<<fabs(coeff.at(j__));
		//			}else{
		//				cout<<" + "<<fabs(coeff.at(j__));
		//			}

		//			for(size_t k_=0; k_<nthE.at(j__).size(); k_++)
		//			{
		//				cout<<" Ç_";
		//				nthE.at(j__).at(k_).printD();

		//				//CentVecN __cnt;
		//				//vector<int> mthE_;
		//				//for(int l__=0; l__<nthE.at(j__).at(k_).size(); l__++)
		//				//	mthE_.push_back(nthE.at(j__).at(k_).getD(l__));

		//				//cout<<"{";
		//				//__cnt.Corr2FactN(mthE_,"NOEXPS");  
		//				//cout<<"}";
		//				if(nthE.at(j__).at(k_).expo()>1) cout<<"^"<<nthE.at(j__).at(k_).getE();
		//				cout<<" ";
		//			}

		//		}
		//	}

		//}



		void printD()
		{
			cout<<endl;
			cout<<endl;
			cout<<endl;
			cout<<endl;
			cout<<endl;


			cout<<"∑_{";
			//nthE.at(1).at(0).printD();
			//cout<<",";
			//nthE.at(1).at(1).printD();
			for(size_t i_=0; i_<mthE.size(); i_ ++) //if(i_<mthE.size()-2) 
				cout<<mthE.at(i_);//<<" "<<mthE.at(i_+1)<<endl;
			cout<<"}=";


			for(size_t j__=0; j__<nthE.size(); j__++)
			{
				bool check=1;

				for(size_t k_=0; k_<nthE.at(j__).size(); k_++)
				{
					if(nthE.at(j__).at(k_).sumD() <=1) check*=0;
				}

				if(check){
					if(coeff.at(j__)==1){ if(j__!=0) cout<<" + ";
					}else if(coeff.at(j__)==-1){ cout<<" - ";
					}else if(coeff.at(j__)<0){
						cout<<" - "<<fabs(coeff.at(j__));
					}else{
						cout<<" + "<<fabs(coeff.at(j__));
					}

					cout<<"[";
					for(size_t k_=0; k_<nthE.at(j__).size(); k_++)
					{
						CentVecN __cnt;
						vector<int> mthE_;
						for(int l__=0; l__<nthE.at(j__).at(k_).size(); l__++)
							mthE_.push_back(nthE.at(j__).at(k_).getD(l__));

						cout<<"{";
						__cnt.Corr2FactN(mthE_,"NOEXPS");  
						cout<<"}";
					}
					cout<<"]";
				}
			}


			cout<<endl;
			cout<<endl;
			cout<<endl;
			cout<<endl;

		}

		void printD(vector<int> __abin)
		{
			cout<<endl;
			cout<<endl;
			cout<<endl;
			cout<<endl;
			cout<<endl;


			cout<<"∑_{";
			//nthE.at(1).at(0).printD();
			//cout<<",";
			//nthE.at(1).at(1).printD();
			for(size_t i_=0; i_<mthE.size(); i_ ++) //if(i_<mthE.size()-2) 
				cout<<mthE.at(i_);//<<" "<<mthE.at(i_+1)<<endl;
			cout<<"}=";


			for(size_t j__=0; j__<nthE.size(); j__++)
			{
				bool check=1;

				for(size_t k_=0; k_<nthE.at(j__).size(); k_++)
				{
					if(nthE.at(j__).at(k_).sumD() <=1) check*=0;
				}

				if(check){
					if(coeff.at(j__)==1){ if(j__!=0) cout<<" + ";
					}else if(coeff.at(j__)==-1){ cout<<" - ";
					}else if(coeff.at(j__)<0){
						cout<<" - "<<fabs(coeff.at(j__));
					}else{
						cout<<" + "<<fabs(coeff.at(j__));
					}

					cout<<"[";
					for(size_t k_=0; k_<nthE.at(j__).size(); k_++)
					{
						CentVecN __cnt;
						vector<int> mthE_;
						for(int l__=0; l__<nthE.at(j__).at(k_).size(); l__++)
							mthE_.push_back(nthE.at(j__).at(k_).getD(l__));

						cout<<"{";
						Corr2FactN(mthE_,__abin,"NOEXPS");  
						//__cnt.Corr2FactN(mthE_,abin__,"NOEXPS","");  
						//__cnt.printD(abin__,"NOEXPS");
						cout<<"}";
					}
					cout<<"]";
				}
			}


			cout<<endl;
			cout<<endl;
			cout<<endl;
			cout<<endl;

		}


	//	void printD(const vector<int> abinN__)
	//	{
	//		cout<<"∑_{";
	//		//nthE.at(1).at(0).printD();
	//		//cout<<",";
	//		//nthE.at(1).at(1).printD();
	//		for(size_t i_=0; i_<mthE.size(); i_ ++) //if(i_<mthE.size()-2) 
	//			cout<<mthE.at(i_);//<<" "<<mthE.at(i_+1)<<endl;
	//		cout<<"}=";


	//		for(size_t j__=0; j__<nthE.size(); j__++)
	//		{
	//			bool check=1;

	//			for(size_t k_=0; k_<nthE.at(j__).size(); k_++)
	//			{
	//				if(nthE.at(j__).at(k_).sumD() <=1) check*=0;
	//			}

	//			if(check){
	//				if(coeff.at(j__)==1){ if(j__!=0) cout<<" + ";
	//				}else if(coeff.at(j__)==-1){ cout<<" - ";
	//				}else if(coeff.at(j__)<0){
	//					cout<<" - "<<fabs(coeff.at(j__));
	//				}else{
	//					cout<<" + "<<fabs(coeff.at(j__));
	//				}

	//				cout<<"[";
	//				for(size_t k_=0; k_<nthE.at(j__).size(); k_++)
	//				{
	//					vector<int> mthE_;
	//					for(int l__=0; l__<nthE.at(j__).at(k_).size(); l__++)
	//						mthE_.push_back(nthE.at(j__).at(k_).getD(l__));

	//					cout<<"{";

	//					//		CentVecN cnt__;
	//					//		cnt__.Corr2FactN(mthE_,abinN,"NOEXPS");  
	//					Corr2FactN(mthE_,abinN__,"NOEXPS");  
	//					cout<<"}";
	//				}
	//				cout<<"]";
	//			}
	//		}


	//		cout<<endl;
	//		cout<<endl;
	//		cout<<endl;
	//		cout<<endl;
	
	//}

		void Squeeze()
		{

			//re-arrange
			for(size_t j__=0; j__<nthE.size(); j__++)
			{
				//	bool check=1;

				vector<FactVec> fvv_;
				//for(size_t k_=0; k_<nthE.at(j__).size(); k_++)
				//{
				//	if(nthE.at(j__).at(k_).sumD() <=1) check*=0;
				//}

				//						if(check) 
				for(size_t k_=0; k_<nthE.at(j__).size(); k_++) fvv_.push_back(nthE.at(j__).at(k_));

				//cout<<"before:"<<endl;
				//for(size_t _i=0; _i<fvv_.size();_i++){cout<<"Ç_";fvv_.at(_i).printD();cout<<" ";}
				//cout<<endl;

				std::sort(fvv_.begin(),fvv_.end(),sortsumD);

				//cout<<"after:"<<endl;
				//for(size_t _i=0; _i<fvv_.size();_i++){cout<<"Ç_";fvv_.at(_i).printD();cout<<" ";}
				//cout<<endl;

				nthE.at(j__).swap(fvv_);

			}


			for(size_t j__=0; j__<nthE.size(); j__++)
			{
				vector<FactVec> fvv_;
				for(size_t k_=0; k_<nthE.at(j__).size(); k_++) fvv_.push_back(nthE.at(j__).at(k_));

				for(size_t i__=0; i__<fvv_.size(); i__++)
				{
					if(fvv_.at(i__).expo()==0) continue;
					for(size_t k__=i__+1; k__<fvv_.size(); k__++)
					{
						if(fvv_.at(i__) == fvv_.at(k__))
						{
							fvv_.at(i__).addE(fvv_.at(k__).expo());
							fvv_.at(k__).setE(0);
						}
					}
				}

				vector<FactVec> fvv__;
				for(size_t i__=0; i__<fvv_.size(); i__++) if(fvv_.at(i__).expo()!=0)fvv__.push_back(fvv_.at(i__));
				nthE.at(j__).swap(fvv__);

				vector<FactVec>().swap(fvv_);
			}



			vector<FactVec> __fvv;
			for(size_t i__=0; i__<nthE.size(); i__++)
			{
				vector<int> __v;
				for(size_t _j=0; _j<nthE.at(i__).size(); _j++)
				{
					for(int _k=0; _k<nthE.at(i__).at(_j).size(); _k++)
						__v.push_back(nthE.at(i__).at(_j).getD(_k));

					__v.push_back(nthE.at(i__).at(_j).getE());
				}

				FactVec temp__fvv;
				temp__fvv.setV(__v);
				__fvv.push_back(temp__fvv);

			}

			//cout<<endl;
			//for(size_t j__=0; j__<__fvv.size(); j__++){cout<<" "; __fvv.at(j__).printD(); cout<<" ";}
			//cout<<endl;

			for(size_t j__=0; j__<__fvv.size(); j__++)
			{
				for(size_t i__=j__+1; i__<__fvv.size(); i__++)
				{
					if(__fvv.at(j__)==__fvv.at(i__)){
						coeff.at(j__) += coeff.at(i__);
						coeff.at(i__) *=0;
						nthE.at(i__).at(0) *=0;
					}
				}
			}

			vector<FactVec>().swap(__fvv);

			vector<vector<FactVec> > nthE__;
			vector<int> coeff__;

			for(size_t j__=0; j__<nthE.size(); j__++) 
			{
				bool check=1;

				for(size_t k_=0; k_<nthE.at(j__).size(); k_++)
				{
					if(nthE.at(j__).at(k_).sumD() <=1) check*=0;
				}

				if(check==1 && coeff.at(j__)!=0){
					nthE__.push_back(nthE.at(j__));
					coeff__.push_back(coeff.at(j__));
				}
			}

			nthE.swap(nthE__); coeff.swap(coeff__);

		}


		void printE()
		{

			//					this->Squeeze();
			for(size_t j__=0; j__<nthE.size(); j__++)
			{
				bool check=1;

				for(size_t k_=0; k_<nthE.at(j__).size(); k_++)
				{
					if(nthE.at(j__).at(k_).sumD() <=1) check*=0;
				}

				if(check){
					if(coeff.at(j__)==1){ if(j__!=0) cout<<" + ";
					}else if(coeff.at(j__)==-1){ cout<<" - ";
					}else if(coeff.at(j__)<0){
						cout<<" - "<<fabs(coeff.at(j__));
					}else{
						cout<<" + "<<fabs(coeff.at(j__));
					}

					for(size_t k_=0; k_<nthE.at(j__).size(); k_++)
					{
						cout<<" Ç_";
						nthE.at(j__).at(k_).printD();
						if(nthE.at(j__).at(k_).expo()>1) cout<<"^"<<nthE.at(j__).at(k_).getE();
						cout<<" ";
					}

				}
			}

		}



};



class CumulantVec : public CentVecN
{
	private:
		int mdim;
		int nsize;

		vector<vector<FactVec> > Kappa;
		vector<vector<FactVec> > Mumn;
		vector<vector<long> > cmn;

		vector<vector<vector<FactVec> > >nthC;
		vector<vector<long> > coeffC;

		bool CALCHECK;
	public:
		CumulantVec(const size_t mdim_=0) : CentVecN()
	{
		mdim=mdim_;
		nsize=0;
		CALCHECK=0;
		Kappa.resize(mdim);
		Mumn.resize(mdim);
		cmn.resize(mdim);
		MOMENT="Cumulant";
	}

		~CumulantVec()
		{
			vector<vector<FactVec> >().swap(Kappa);
			vector<vector<FactVec> >().swap(Mumn);
			vector<vector<long> >().swap(cmn);

			vector<vector<vector<FactVec> > >().swap(nthC);
			vector<vector<long> >().swap(coeffC);
		}


		vector<vector<vector<FactVec> > > getnthC(){return nthC;}
		vector<vector<long> > getcoeffC() {return coeffC;}


		void clearnthC(){nthC.clear();}
		void clearcoeffC() {coeffC.clear();}

		bool calcheck() const {return CALCHECK;}

		void Cumulant2Cumulant(const int mdim__)
		{
			cout<<"Func to convert cumulant K_m to cumulant K_m"<<endl; 

			for(int k_=1; k_<=mdim__; k_++)
			{

				vector<FactVec> Kappa_;
				vector<FactVec> Mumn_;
				vector<long> cmn_;

				for(int i_=1; i_<=k_-1; i_++)
				{
					FactVec TEMP_Kappa_(1);
					TEMP_Kappa_.setD(0,i_);
					FactVec TEMP_Mumn_(1);
					TEMP_Mumn_.setD(0,k_-i_);

					cmn_.push_back(-calccomb(k_-1,i_-1));
					Kappa_.push_back(TEMP_Kappa_);
					Mumn_.push_back(TEMP_Mumn_);
				}


				//if(k_==mdim__){


				cout<<"K_"<<k_<<"= µ_"<<k_;
				for(size_t l__=0; l__<cmn_.size(); l__++){

					if(cmn_.at(l__)==1){ if(l__!=0) cout<<" + ";
					}else if(cmn_.at(l__)==-1){ cout<<" - ";
					}else if(cmn_.at(l__)<0){
						cout<<" - "<<fabs(cmn_.at(l__));
					}else{
						cout<<" + "<<fabs(cmn_.at(l__));
					}

					cout<<" K_";
					Kappa_.at(l__).printD();
					cout<<" µ_";
					Mumn_.at(l__).printD();
					cout<<" ";
				}
				cout<<endl;

				//}

				cmn_.clear();
				Kappa_.clear();
				Mumn_.clear();

			}
		}




		void Cumulant2Moment(const int mdim__, string FORMAT="")
		{

			//					cout<<"(Inside Cumulant2Moment) mdim ="<<mdim__<<endl;

			//cout<<"-------------------------------------------------------------------------------"<<endl;
			//cout<<"K_1= Mu_1"<<endl;
			//cout<<"K_2= Mu_2 -  K_1 Mu_1"<<endl;
			//cout<<"K_3= Mu_3 -  K_1 Mu_2  - 2 K_2 Mu_1"<<endl; 
			//cout<<"K_4= Mu_4 -  K_1 Mu_3  - 3 K_2 Mu_2  - 3 K_3 Mu_1"<<endl; 
			//cout<<"K_5= Mu_5 -  K_1 Mu_4  - 4 K_2 Mu_3  - 6 K_3 Mu_2  - 4 K_4 Mu_1"<<endl; 
			//cout<<"K_6= Mu_6 -  K_1 Mu_5  - 5 K_2 Mu_4  - 10 K_3 Mu_3  - 10 K_4 Mu_2  - 5 K_5 Mu_1"<<endl; 
			//cout<<"-------------------------------------------------------------------------------"<<endl;


			if(FORMAT!="NOEXPS")cout<<"Func to convert cumulant K_m to moments µ_m = <∆N^m>"<<endl; 

			for(int k_=1; k_<=mdim__; k_++)
			{

				vector<FactVec> Kappa_;
				vector<FactVec> Mumn_;
				vector<long> cmn_;

				for(int i_=1; i_<=k_-1; i_++)
				{
					FactVec TEMP_Kappa_(1);
					TEMP_Kappa_.setD(0,i_);
					FactVec TEMP_Mumn_(1);
					TEMP_Mumn_.setD(0,k_-i_);

					cmn_.push_back(-calccomb(k_-1,i_-1));
					Kappa_.push_back(TEMP_Kappa_);
					Mumn_.push_back(TEMP_Mumn_);
				}

				vector<vector<FactVec> > temp_nthC;
				vector<long > temp_coeffC;


				FactVec TEMP_Mumn__(1);
				TEMP_Mumn__.setD(0,k_);

				vector<FactVec > TEMP_TEMP_Mumn__;
				TEMP_TEMP_Mumn__.push_back(TEMP_Mumn__);
				temp_nthC.push_back(TEMP_TEMP_Mumn__);

				temp_coeffC.push_back(1.);


				for(size_t i__=0; i__<Mumn_.size(); i__++)
				{

					for(size_t i=0; i<nthC.at(Kappa_.at(i__).getD(0)-1).size(); i++)
					{
						vector<FactVec> temp_temp_nthC;

						for(size_t ll_=0; ll_<nthC.at(Kappa_.at(i__).getD(0)-1).at(i).size(); ll_++)
						{
							temp_temp_nthC.push_back(nthC.at(Kappa_.at(i__).getD(0)-1).at(i).at(ll_));
						}
						temp_temp_nthC.push_back(Mumn_.at(i__));

						temp_nthC.push_back(temp_temp_nthC);
						temp_coeffC.push_back(cmn_.at(i__)*coeffC.at(Kappa_.at(i__).getD(0)-1).at(i));

					}
				}

				nthC.push_back(temp_nthC);
				coeffC.push_back(temp_coeffC);

				vector<long>().swap(cmn_);
				vector<FactVec>().swap(Kappa_);
				vector<FactVec>().swap(Mumn_);

				vector<vector<FactVec> >().swap(temp_nthC);
				vector<long >().swap(temp_coeffC);
				vector<FactVec >().swap(TEMP_TEMP_Mumn__);
			}


			for(size_t i__=0; i__<nthC.size(); i__++)
			{
				vector<FactVec > Kexps;
				vector<long > Cexps;
				if(FORMAT!="NOEXPS") cout<<" K_"<<i__+1<<"= ";
				for(size_t j__=0; j__<nthC.at(i__).size(); j__++)
				{

					long prefact = coeffC.at(i__).at(j__);

					vector<short > addpow;
					addpow.resize(i__+2);

					for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
					{
						addpow.at(nthC.at(i__).at(j__).at(k__).getD(0)) += 1;
					}


					vector<int> vkexps;
					for(size_t l__=0; l__<addpow.size(); l__++)
					{
						if(addpow.at(l__)>0){
							vkexps.push_back(l__);
							vkexps.push_back(addpow.at(l__));
						}
					}

					FactVec temp_Kexps;
					temp_Kexps.setV(vkexps);


					Kexps.push_back(temp_Kexps);
					Cexps.push_back(prefact);

					vector<int>().swap(vkexps);
					vector<short >().swap(addpow);
				}

				vector<vector<FactVec> > temp_nthC;
				vector<long > temp_coeffC;

				for(size_t __k=0; __k<Kexps.size(); __k++) 
				{
					if(Kexps.at(__k).sumD()==0) continue;
					double tcoeff=Cexps.at(__k);

					for(size_t __j=__k+1; __j<Kexps.size(); __j++) 
					{
						if(Kexps.at(__k) == Kexps.at(__j)) 
						{ 
							tcoeff += Cexps.at(__j);	
							Kexps.at(__j) *= 0;
						}
					}

					if(FORMAT!="NOEXPS"){ 
						if(tcoeff==1){ if(__k!=0) cout<<" + ";
						}else if(tcoeff==-1){ cout<<" - ";
						}else if(tcoeff<0){
							cout<<" - "<<fabs(tcoeff);
						}else{
							cout<<" + "<<fabs(tcoeff);
						}
					}

					vector<FactVec> temp_temp_nthC;

					for(int __i=0; __i<Kexps.at(__k).size()-1; __i += 2) 
					{
						if(FORMAT!="NOEXPS"){ 
							cout<<" µ_";
							cout<<Kexps.at(__k).getD(__i);

							if(Kexps.at(__k).getD(__i+1)>1)
							{
								cout<<"^";
								cout<<Kexps.at(__k).getD(__i+1);
							}
						}

						FactVec temp_temp_temp_nthC(1);

						temp_temp_temp_nthC.setD(0,Kexps.at(__k).getD(__i));
						temp_temp_temp_nthC.setE(Kexps.at(__k).getD(__i+1));

						temp_temp_nthC.push_back(temp_temp_temp_nthC);
					}


					temp_nthC.push_back(temp_temp_nthC);
					temp_coeffC.push_back(tcoeff);

					vector<FactVec>().swap(temp_temp_nthC);
				}

				if(FORMAT!="NOEXPS")cout<<endl;

				Kexps.clear();
				Cexps.clear();

				vector<FactVec >().swap(Kexps);
				vector<long >().swap(Cexps);

				nthC.at(i__).swap(temp_nthC);
				coeffC.at(i__).swap(temp_coeffC);

			}
			if(FORMAT!="NOEXPS")cout<<endl;

			CALCHECK=1;
		}





		void Cumulant2CorrMoment(const int mdim__, string FORMAT="")
		{
			//					cout<<"(Inside Cumulant2CorrMoment) mdim ="<<mdim__<<endl;

			if(FORMAT!="NOEXPS")cout<<"Func to convert cumulant K_m to correlation variables Ç_m = (∆N-<∆N>)^m or central moments w.r.to ∆N"<<endl; 
			if(CALCHECK!=1) this->Cumulant2Moment(mdim__,"NOEXPS");

			size_t i__=mdim__-1;
			//for(size_t i__=mdim__-1; i__<mdim__; i__++)
			{
				if(FORMAT!="NOEXPS")cout<<" K_"<<i__+1<<"= ";

				for(size_t j__=0; j__<nthC.at(i__).size(); j__++)
				{
					bool CHECK=1;

					for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
					{
						if(nthC.at(i__).at(j__).at(k__).getD(0)==1) CHECK *=0;
					}

					if(CHECK!=0){

						long prefact = coeffC.at(i__).at(j__);

						//if(FORMAT!="NOEXPS"){
						if(prefact==1){ if(j__!=0) cout<<" + ";
						}else if(prefact==-1){ cout<<" - ";
						}else if(prefact<0){
							cout<<" - "<<fabs(prefact);
						}else{
							cout<<" + "<<fabs(prefact);
						}			
						//}

						for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
						{
							//if(FORMAT!="NOEXPS"){
							cout<<" Ç_";
							nthC.at(i__).at(j__).at(k__).printD();
							if(nthC.at(i__).at(j__).at(k__).expo()>1) cout<<"^"<<nthC.at(i__).at(j__).at(k__).expo();
							//}

						}

					}

				}
				if(FORMAT!="NOEXPS")cout<<endl;
			}
			if(FORMAT!="NOEXPS")cout<<endl;

		}





		void Cumulant2CentralMoment(const int mdim__, string FORMAT="")
		{
			if(FORMAT!="NOEXPS")cout<<"Func to convert cumulant K_m to central moments C_m = <(N-<N>)^m>"<<endl; 

			if(CALCHECK!=1) this->Cumulant2Moment(mdim__,"NOEXPS");

			size_t i__=mdim__-1;

//			for(size_t i__=0; i__<nthC.size(); i__++)
			{
				if(FORMAT!="NOEXPS")cout<<" K_"<<i__+1<<"= ";

				for(size_t j__=0; j__<nthC.at(i__).size(); j__++)
				{
					bool CHECK=1;

					for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
					{
						if(nthC.at(i__).at(j__).at(k__).getD(0)==1) CHECK *=0;
					}

					if(CHECK!=0){

						long prefact = coeffC.at(i__).at(j__);

						if(FORMAT!="NOEXPS"){
							if(prefact==1){ if(j__!=0) cout<<" + ";
							}else if(prefact==-1){ cout<<" - ";
							}else if(prefact<0){
								cout<<" - "<<fabs(prefact);
							}else{
								cout<<" + "<<fabs(prefact);
							}			
						}

						for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
						{
							CentVecN cnt_;
							vector<int> mthC;
							mthC.push_back(nthC.at(i__).at(j__).at(k__).getD(0));

							if(FORMAT!="NOEXPS"){
								//cout<<" C_";
								//nthC.at(i__).at(j__).at(k__).printD();
								cout<<"(";
								cnt_.Corr2CentN(mthC,"","NOEXPS");
								cout<<" )";
								if(nthC.at(i__).at(j__).at(k__).expo()>1) cout<<"^"<<nthC.at(i__).at(j__).at(k__).expo();
							}

						}

					}


				}
				if(FORMAT!="NOEXPS")cout<<endl;
			}
			if(FORMAT!="NOEXPS")cout<<endl;

		}

		void Cumulant2FactMoment(const int mdim__, string FORMAT="")
		{
			if(FORMAT!="NOEXPS") cout<<"Func to convert cumulant K_m to factorial moments f_m = <N(N-1)..(N-m+1)>"<<endl; 

			if(CALCHECK!=1) this->Cumulant2Moment(mdim__,"NOEXPS");

			size_t i__=mdim__-1;

			//for(size_t i__=0; i__<nthC.size(); i__++)
			{
				if(FORMAT!="NOEXPS")cout<<" K_"<<i__+1<<"= ";
				if(i__==0){
					cout<<"f10 - f01"<<endl;
					cout<<endl;
				//	continue;
				}else{

				for(size_t j__=0; j__<nthC.at(i__).size(); j__++)
				{


					bool CHECK=1;

					for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
					{
						if(nthC.at(i__).at(j__).at(k__).getD(0)==1) CHECK *=0;
					}

					if(CHECK!=0){

						long prefact = coeffC.at(i__).at(j__);

						if(FORMAT!="NOEXPS"){
							if(prefact==1){ if(j__!=0) cout<<" + ";
							}else if(prefact==-1){ cout<<" - ";
							}else if(prefact<0){
								cout<<" - "<<fabs(prefact);
							}else{
								cout<<" + "<<fabs(prefact);
							}			
						}

						for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
						{
							CentVecN cnt_;
							vector<int> mthC;
							mthC.push_back(nthC.at(i__).at(j__).at(k__).getD(0));

							if(FORMAT!="NOEXPS"){
								//cout<<" C_";
								//nthC.at(i__).at(j__).at(k__).printD();
								cout<<"(";
								cnt_.Corr2FactN(mthC,"NOEXPS");
								cout<<" )";
								if(nthC.at(i__).at(j__).at(k__).expo()>1) cout<<"^"<<nthC.at(i__).at(j__).at(k__).expo();
							}

						}

					}


				}
				if(FORMAT!="NOEXPS")cout<<endl;
				if(FORMAT!="NOEXPS")cout<<endl;
			}
			}
			if(FORMAT!="NOEXPS")cout<<endl;

		}

//here
		void Cumulant2FactMoment(const int mdim__,vector<int> abin__, string FORMAT="")
		{
			if(FORMAT!="NOEXPS") cout<<"Func to convert cumulant K_m to factorial moments f_m = <N(N-1)..(N-m+1)>"<<endl; 

			if(abin__.size()!=2){
				cout << "\033[1;35mWARNING !! Input efficiency vector size ≠2 detected for bi-variate cumulant \033[0m\n"<<" "<<endl;
				abin__.resize(2);
			}

			if(CALCHECK!=1) this->Cumulant2Moment(mdim__,"NOEXPS");

			size_t i__=mdim__-1;
			//for(size_t i__=0; i__<nthC.size(); i__++)
			{
				if(FORMAT!="NOEXPS")cout<<" K_"<<i__+1<<"= ";

				if(i__==0){
					FactMomN fact_;
					vector<int> mthC;
					mthC.push_back(1);
					mthC.push_back(0);

					if(FORMAT!="NOEXPS"){
						cout<<"(";
						fact_.SetFmnabN(mthC,abin__);
						fact_.printD("eff");
						cout<<" )";
					}

					mthC.clear();
					mthC.push_back(0);
					mthC.push_back(1);

					if(FORMAT!="NOEXPS"){
						cout<<"-(";
						fact_.SetFmnabN(mthC,abin__);
						fact_.printD("eff");
						cout<<" )";
					}

					cout<<endl;
					cout<<endl;
					cout<<endl;
					//continue;
				}else{

				for(size_t j__=0; j__<nthC.at(i__).size(); j__++)
				{


					bool CHECK=1;

					for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
					{
						if(nthC.at(i__).at(j__).at(k__).getD(0)==1) CHECK *=0;
					}

					if(CHECK!=0){

						long prefact = coeffC.at(i__).at(j__);

						if(FORMAT!="NOEXPS"){
							if(prefact==1){ if(j__!=0) cout<<" + ";
							}else if(prefact==-1){ cout<<" - ";
							}else if(prefact<0){
								cout<<" - "<<fabs(prefact);
							}else{
								cout<<" + "<<fabs(prefact);
							}			
						}

						for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
						{
							CentVecN cnt_;
							vector<int> mthC;
							mthC.push_back(nthC.at(i__).at(j__).at(k__).getD(0));

							if(FORMAT!="NOEXPS"){
								//cout<<" C_";
								//nthC.at(i__).at(j__).at(k__).printD();
								cout<<"(";
								cnt_.Corr2FactN(mthC,abin__,"NOEXPS");
								cout<<" )";
								if(nthC.at(i__).at(j__).at(k__).expo()>1) cout<<"^"<<nthC.at(i__).at(j__).at(k__).expo();
							}

						}

					}


				}
				if(FORMAT!="NOEXPS")cout<<endl;
				if(FORMAT!="NOEXPS")cout<<endl;
			}
			}
			if(FORMAT!="NOEXPS")cout<<endl;

		}





		//Overriding the following function from the base class
		vector<FactVec> CalcFactVec(const int __mdim, string WError="")
		{
			const int mdim__ = ((WError=="NOERROR") ? __mdim : 2*__mdim);

			//this->CalcVariance(mdim__);

			cout<<"order : "<<__mdim<<" "<<mdim__<<endl;
			vector<FactVec> fvv;

			this->clearnthC();
			this->clearcoeffC();

			//if(CALCHECK!=1) 
			this->Cumulant2Moment(mdim__,"NOEXPS");

			size_t i__=nthC.size()-1;

			if(i__==0){
				//cout<<"f10 - f01"<<endl;
				FactVec F10(2);
				F10.setD(0,1);
				FactVec F01(2);
				F01.setD(1,1);
				fvv.push_back(F10);
				fvv.push_back(F01);
				//continue;
			}

			for(size_t j__=0; j__<nthC.at(i__).size(); j__++)
			{

				vector<FactVec> temp_fvv;

				bool CHECK=1;

				for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
				{
					if(nthC.at(i__).at(j__).at(k__).getD(0)==1) CHECK *=0;
				}

				if(CHECK!=0){
					//long prefact = coeffC.at(i__).at(j__);

					for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
					{
						CentVecN cnt_;
						vector<int> mthC;
						mthC.push_back(nthC.at(i__).at(j__).at(k__).getD(0));

						cnt_.Corr2FactN(mthC,"NOEXPS","");
						temp_fvv=cnt_.CalcFactVec("");

					}

				}

				fvv.insert(fvv.end(), temp_fvv.begin(), temp_fvv.end());
			}

			for(size_t j_=0; j_<fvv.size(); j_++){
				for(size_t k_=j_+1; k_<fvv.size(); k_++) {
					if(fvv.at(j_) == fvv.at(k_)) { fvv.erase(fvv.begin() + k_); }
				}
			}

			cout<<"Factorial Moments to be calculated are :"<<endl;

			for(size_t j_=0; j_<fvv.size(); j_++) {cout<<"F_"; fvv.at(j_).printD(); cout<<endl;}

			return fvv;
		}



		void CalcVariance(const int mdim__)//, vector<vector<FactVec > >& Kdelt__, vector<vector<long > >& Cdelt__, vector<FactVec>& fvv__)
		{
			vector<FactVec> fvv;

			if(CALCHECK!=1) this->Cumulant2Moment(mdim__,"NOEXPS");

			size_t i__=nthC.size()-1;


			for(size_t j__=0; j__<nthC.at(i__).size(); j__++)
			{

				bool CHECK=1;

				for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
				{
					if(nthC.at(i__).at(j__).at(k__).getD(0)==1) CHECK *=0;
					//cout<<" ";
					//nthC.at(i__).at(j__).at(k__).printD();
					//cout<<" ";

				}

				if(CHECK!=0){

					for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
					{
						fvv.push_back(nthC.at(i__).at(j__).at(k__));
					}
				}
			}

			for(size_t j_=0; j_<fvv.size(); j_++){
				for(size_t k_=j_+1; k_<fvv.size(); k_++) {
					if(fvv.at(j_) == fvv.at(k_)) { fvv.at(k_) *= 0;}
				}
			}

			vector<FactVec> temp_fvv;

			for(size_t j_=0; j_<fvv.size(); j_++) {if(fvv.at(j_).sumD() != 0) temp_fvv.push_back(fvv.at(j_));}

			fvv.swap(temp_fvv);

			//cout<<"Cumulants to be calculated are :"<<endl;

			//for(size_t j_=0; j_<fvv.size(); j_++) {cout<<"Ç_"; fvv.at(j_).printD(); cout<<endl;}


			vector<vector<FactVec > >Kdelt;
			vector<vector<long > >Cdelt;

			for(size_t __j=0; __j<fvv.size(); __j++)
			{

				//cout<<"∂K_"<<i__+1<<"/∂Ç_";
				//fvv.at(__j).printD();
				//cout<<" = ";
				for(size_t j__=0; j__<nthC.at(i__).size(); j__++)
				{

					bool CHECK=1;

					for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
					{
						if(nthC.at(i__).at(j__).at(k__).getD(0)==1) CHECK *=0;
					}

					if(CHECK!=0){

						bool CHECK2=0;

						for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
						{
							if(fvv.at(__j)==nthC.at(i__).at(j__).at(k__))
							{

								long prefact = coeffC.at(i__).at(j__) * nthC.at(i__).at(j__).at(k__).expo();
								coeffC.at(i__).at(j__)=prefact;

								nthC.at(i__).at(j__).at(k__).setE(nthC.at(i__).at(j__).at(k__).expo()-1); //expo()-1

								CHECK2 += 1;
							}
						}

						coeffC.at(i__).at(j__) *= CHECK2;

					}


				}


				vector<FactVec >temp_Kdelt;
				vector<long >temp_Cdelt;
				for(size_t j__=0; j__<nthC.at(i__).size(); j__++)
				{
					bool CHECK=1;

					for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
					{
						if(nthC.at(i__).at(j__).at(k__).getD(0)==1) CHECK *=0;
						//cout<<" ";
						//nthC.at(i__).at(j__).at(k__).printD();
						//cout<<" ";
					}

					if(CHECK!=0 && coeffC.at(i__).at(j__)!=0 ){

						long prefact = coeffC.at(i__).at(j__);

				//		if(prefact<0){
				//			cout<<" - "<<fabs(prefact);
				//		}else{
				//			cout<<" + "<<fabs(prefact);
				//		}			

						vector<int> vkdelt;
						for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
						{
							if(nthC.at(i__).at(j__).at(k__).expo()>0){
						//		cout<<" Ç_";
						//		nthC.at(i__).at(j__).at(k__).printD();
						//		if(nthC.at(i__).at(j__).at(k__).expo()>1) cout<<"^"<<nthC.at(i__).at(j__).at(k__).expo();
								vkdelt.push_back(nthC.at(i__).at(j__).at(k__).getD(0));
								vkdelt.push_back(nthC.at(i__).at(j__).at(k__).expo());
							}
						}

						FactVec temp_temp_Kdelt;
						temp_temp_Kdelt.setV(vkdelt);

						temp_Kdelt.push_back(temp_temp_Kdelt);
						temp_Cdelt.push_back(prefact);

						vector<int>().swap(vkdelt);

					}


				}
			//	cout<<endl;

				Kdelt.push_back(temp_Kdelt);
				Cdelt.push_back(temp_Cdelt);


				nthC.clear();
				coeffC.clear();
				this->Cumulant2Moment(mdim__,"NOEXPS");
			}



//			cout<<endl;
			cout<<"Var(K_"<<i__+1<<")=";
			for(size_t j_=0; j_<fvv.size(); j_++){
				for(size_t k_=0; k_<fvv.size(); k_++){

					if(fvv.at(k_).sumD()==0) continue;

					if(fvv.at(j_) == fvv.at(k_)) 
					{ 
						if(j_!=0){	
							cout<<"(";
							for(size_t __i=0; __i<Kdelt.at(j_).size(); __i++)
							{
								long prefact = Cdelt.at(j_).at(__i);

								if(prefact==1){ if(__i!=0) cout<<" + ";
								}else if(prefact==-1){ cout<<" - ";
								}else if(prefact<0){
									cout<<" - "<<fabs(prefact);
								}else{
									cout<<" + "<<fabs(prefact);
								}

								for(int __k=0; __k<Kdelt.at(j_).at(__i).size()-1; __k += 2) 
								{
									cout<<" Ç_";
									cout<<Kdelt.at(j_).at(__i).getD(__k);

									if(Kdelt.at(j_).at(__i).getD(__k+1)>1)
									{
										cout<<"^";
										cout<<Kdelt.at(j_).at(__i).getD(__k+1);
									}
								}
							}
							cout<<" )^2 ";
						}

						vector<int> _mth;

						_mth.push_back(fvv.at(j_).getD(0)); 
						_mth.push_back(fvv.at(k_).getD(0)); 

						// call central moment error estimation class
						CentErrorN cerr(_mth);

						//calculate the expression of error 
//come back
						cerr.CalcCentError("");
						cout<<" (";
						cerr.printE();
						cout<<") ";

						//cout<<" ∑(";
						//fvv.at(j_).printD();
						//cout<<",";
						//fvv.at(k_).printD();
						//cout<<")";

					}else{
						if(fvv.at(j_).sumD()>fvv.at(k_).sumD()){

							cout<<"2 ";

							if(j_!=0){
								cout<<"(";
								for(size_t __i=0; __i<Kdelt.at(j_).size(); __i++)
								{
									long prefact = Cdelt.at(j_).at(__i);

									if(prefact==1){ if(__i!=0) cout<<" + ";
									}else if(prefact==-1){ cout<<" - ";
									}else if(prefact<0){
										cout<<" - "<<fabs(prefact);
									}else{
										cout<<" + "<<fabs(prefact);
									}

									for(int __k=0; __k<Kdelt.at(j_).at(__i).size()-1; __k += 2) 
									{
										cout<<" Ç_";
										cout<<Kdelt.at(j_).at(__i).getD(__k);

										if(Kdelt.at(j_).at(__i).getD(__k+1)>1)
										{
											cout<<"^";
											cout<<Kdelt.at(j_).at(__i).getD(__k+1);
										}
									}
								}
								cout<<" )";
							}

							cout<<"(";
							for(size_t __i=0; __i<Kdelt.at(k_).size(); __i++)
							{
								long prefact = Cdelt.at(k_).at(__i);

								if(prefact==1){ if(__i!=0) cout<<" + ";
								}else if(prefact==-1){ cout<<" - ";
								}else if(prefact<0){
									cout<<" - "<<fabs(prefact);
								}else{
									cout<<" + "<<fabs(prefact);
								}

								for(int __k=0; __k<Kdelt.at(k_).at(__i).size()-1; __k += 2) 
								{
									cout<<" Ç_";
									cout<<Kdelt.at(k_).at(__i).getD(__k);

									if(Kdelt.at(k_).at(__i).getD(__k+1)>1)
									{
										cout<<"^";
										cout<<Kdelt.at(k_).at(__i).getD(__k+1);
									}
								}
							}
							cout<<" )";


							vector<int> _mth;

							_mth.push_back(fvv.at(j_).getD(0)); 
							_mth.push_back(fvv.at(k_).getD(0)); 

							// call central moment error estimation class
							CentErrorN cerr(_mth);

							//calculate the expression of error 
							cerr.CalcCentError("");
							cout<<" (";
							cerr.printE();
							cout<<") ";



							//cout<<" ∑(";
							//fvv.at(j_).printD();
							//cout<<",";
							//fvv.at(k_).printD();
							//cout<<")";

						}else{
							cout<<"2 ";


							cout<<"(";
							for(size_t __i=0; __i<Kdelt.at(k_).size(); __i++)
							{
								long prefact = Cdelt.at(k_).at(__i);

								if(prefact==1){ if(__i!=0) cout<<" + ";
								}else if(prefact==-1){ cout<<" - ";
								}else if(prefact<0){
									cout<<" - "<<fabs(prefact);
								}else{
									cout<<" + "<<fabs(prefact);
								}

								for(int __k=0; __k<Kdelt.at(k_).at(__i).size()-1; __k += 2) 
								{
									cout<<" Ç_";
									cout<<Kdelt.at(k_).at(__i).getD(__k);

									if(Kdelt.at(k_).at(__i).getD(__k+1)>1)
									{
										cout<<"^";
										cout<<Kdelt.at(k_).at(__i).getD(__k+1);
									}
								}
							}
							cout<<" )";

							if(j_!=0){
								cout<<"(";
								for(size_t __i=0; __i<Kdelt.at(j_).size(); __i++)
								{
									long prefact = Cdelt.at(j_).at(__i);

									if(prefact==1){ if(__i!=0) cout<<" + ";
									}else if(prefact==-1){ cout<<" - ";
									}else if(prefact<0){
										cout<<" - "<<fabs(prefact);
									}else{
										cout<<" + "<<fabs(prefact);
									}

									for(int __k=0; __k<Kdelt.at(j_).at(__i).size()-1; __k += 2) 
									{
										cout<<" Ç_";
										cout<<Kdelt.at(j_).at(__i).getD(__k);

										if(Kdelt.at(j_).at(__i).getD(__k+1)>1)
										{
											cout<<"^";
											cout<<Kdelt.at(j_).at(__i).getD(__k+1);
										}
									}
								}
								cout<<" )";
							}

							vector<int> _mth;

							_mth.push_back(fvv.at(k_).getD(0)); 
							_mth.push_back(fvv.at(j_).getD(0)); 

							// call central moment error estimation class
							CentErrorN cerr(_mth);

							//calculate the expression of error 
							cerr.CalcCentError("");
							cout<<" (";
							cerr.printE();
							cout<<") ";

							//cout<<" ∑(";
							//fvv.at(k_).printD();
							//cout<<",";
							//fvv.at(j_).printD();
							//cout<<")";
						}
					}


					if(j_==fvv.size()-1 && k_==fvv.size()-1){
						cout<<"";
					}else{
						cout<<" + ";
					}
				}
				fvv.at(j_) *= 0;
			}
			cout<<endl;

//			Kdelt__ = Kdelt;
//			Cdelt__ = Cdelt;
//			fvv__ = fvv;
//
			vector<vector<FactVec > >().swap(Kdelt);
			vector<vector<long > >().swap(Cdelt);
			vector<FactVec >().swap(fvv);
		}


		void CalcVariance(const int mdim__, vector<vector<FactVec > >& Kdelt__, vector<vector<long > >& Cdelt__, vector<FactVec>& fvv__)
		{
			vector<FactVec> fvv;

			if(CALCHECK!=1) this->Cumulant2Moment(mdim__,"NOEXPS");

			size_t i__=nthC.size()-1;


			for(size_t j__=0; j__<nthC.at(i__).size(); j__++)
			{

				bool CHECK=1;

				for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
				{
					if(nthC.at(i__).at(j__).at(k__).getD(0)==1) CHECK *=0;
					//cout<<" ";
					//nthC.at(i__).at(j__).at(k__).printD();
					//cout<<" ";

				}

				if(CHECK!=0){

					for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
					{
						fvv.push_back(nthC.at(i__).at(j__).at(k__));
					}
				}
			}

			for(size_t j_=0; j_<fvv.size(); j_++){
				for(size_t k_=j_+1; k_<fvv.size(); k_++) {
					if(fvv.at(j_) == fvv.at(k_)) { fvv.at(k_) *= 0;}
				}
			}

			vector<FactVec> temp_fvv;

			for(size_t j_=0; j_<fvv.size(); j_++) {if(fvv.at(j_).sumD() != 0) temp_fvv.push_back(fvv.at(j_));}

			fvv.swap(temp_fvv);

		//	cout<<"Cumulants to be calculated are :"<<endl;

		//	for(size_t j_=0; j_<fvv.size(); j_++) {cout<<"Ç_"; fvv.at(j_).printD(); cout<<endl;}


			vector<vector<FactVec > >Kdelt;
			vector<vector<long > >Cdelt;

			for(size_t __j=0; __j<fvv.size(); __j++)
			{

				//cout<<"∂K_"<<i__+1<<"/∂Ç_";
				//fvv.at(__j).printD();
				//cout<<" = ";
				for(size_t j__=0; j__<nthC.at(i__).size(); j__++)
				{

					bool CHECK=1;

					for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
					{
						if(nthC.at(i__).at(j__).at(k__).getD(0)==1) CHECK *=0;
					}

					if(CHECK!=0){

						bool CHECK2=0;

						for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
						{
							if(fvv.at(__j)==nthC.at(i__).at(j__).at(k__))
							{

								long prefact = coeffC.at(i__).at(j__) * nthC.at(i__).at(j__).at(k__).expo();
								coeffC.at(i__).at(j__)=prefact;

								nthC.at(i__).at(j__).at(k__).setE(nthC.at(i__).at(j__).at(k__).expo()-1); //expo()-1

								CHECK2 += 1;
							}
						}

						coeffC.at(i__).at(j__) *= CHECK2;

					}


				}


				vector<FactVec >temp_Kdelt;
				vector<long >temp_Cdelt;
				for(size_t j__=0; j__<nthC.at(i__).size(); j__++)
				{
					bool CHECK=1;

					for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
					{
						if(nthC.at(i__).at(j__).at(k__).getD(0)==1) CHECK *=0;
						//cout<<" ";
						//nthC.at(i__).at(j__).at(k__).printD();
						//cout<<" ";
					}

					if(CHECK!=0 && coeffC.at(i__).at(j__)!=0 ){

						long prefact = coeffC.at(i__).at(j__);

						//if(prefact<0){
						//	cout<<" - "<<fabs(prefact);
						//}else{
						//	cout<<" + "<<fabs(prefact);
						//}			

						vector<int> vkdelt;
						for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
						{
							if(nthC.at(i__).at(j__).at(k__).expo()>0){
								//cout<<" Ç_";
								//nthC.at(i__).at(j__).at(k__).printD();
								//if(nthC.at(i__).at(j__).at(k__).expo()>1) cout<<"^"<<nthC.at(i__).at(j__).at(k__).expo();
								vkdelt.push_back(nthC.at(i__).at(j__).at(k__).getD(0));
								vkdelt.push_back(nthC.at(i__).at(j__).at(k__).expo());
							}
						}

						FactVec temp_temp_Kdelt;
						temp_temp_Kdelt.setV(vkdelt);

						temp_Kdelt.push_back(temp_temp_Kdelt);
						temp_Cdelt.push_back(prefact);

						vector<int>().swap(vkdelt);

					}


				}
				cout<<endl;

				Kdelt.push_back(temp_Kdelt);
				Cdelt.push_back(temp_Cdelt);


				nthC.clear();
				coeffC.clear();
				this->Cumulant2Moment(mdim__,"NOEXPS");
			}


/*
			cout<<endl;
			cout<<"Var(K_"<<i__+1<<")=";
			for(size_t j_=0; j_<fvv.size(); j_++){
				for(size_t k_=0; k_<fvv.size(); k_++){

					if(fvv.at(k_).sumD()==0) continue;

					if(fvv.at(j_) == fvv.at(k_)) 
					{ 
						if(j_!=0){	
							cout<<"(";
							for(size_t __i=0; __i<Kdelt.at(j_).size(); __i++)
							{
								long prefact = Cdelt.at(j_).at(__i);

								if(prefact==1){ if(__i!=0) cout<<" + ";
								}else if(prefact==-1){ cout<<" - ";
								}else if(prefact<0){
									cout<<" - "<<fabs(prefact);
								}else{
									cout<<" + "<<fabs(prefact);
								}

								for(int __k=0; __k<Kdelt.at(j_).at(__i).size()-1; __k += 2) 
								{
									cout<<" Ç_";
									cout<<Kdelt.at(j_).at(__i).getD(__k);

									if(Kdelt.at(j_).at(__i).getD(__k+1)>1)
									{
										cout<<"^";
										cout<<Kdelt.at(j_).at(__i).getD(__k+1);
									}
								}
							}
							cout<<" )^2 ";
						}

						vector<int> _mth;

						_mth.push_back(fvv.at(j_).getD(0)); 
						_mth.push_back(fvv.at(k_).getD(0)); 

						// call central moment error estimation class
						CentErrorN cerr(_mth);

						//calculate the expression of error 
//come back
						cerr.CalcCentError("");
						cout<<" (";
						cerr.printE();
						cout<<") ";

						//cout<<" ∑(";
						//fvv.at(j_).printD();
						//cout<<",";
						//fvv.at(k_).printD();
						//cout<<")";

					}else{
						if(fvv.at(j_).sumD()>fvv.at(k_).sumD()){

							cout<<"2 ";

							if(j_!=0){
								cout<<"(";
								for(size_t __i=0; __i<Kdelt.at(j_).size(); __i++)
								{
									long prefact = Cdelt.at(j_).at(__i);

									if(prefact==1){ if(__i!=0) cout<<" + ";
									}else if(prefact==-1){ cout<<" - ";
									}else if(prefact<0){
										cout<<" - "<<fabs(prefact);
									}else{
										cout<<" + "<<fabs(prefact);
									}

									for(int __k=0; __k<Kdelt.at(j_).at(__i).size()-1; __k += 2) 
									{
										cout<<" Ç_";
										cout<<Kdelt.at(j_).at(__i).getD(__k);

										if(Kdelt.at(j_).at(__i).getD(__k+1)>1)
										{
											cout<<"^";
											cout<<Kdelt.at(j_).at(__i).getD(__k+1);
										}
									}
								}
								cout<<" )";
							}

							cout<<"(";
							for(size_t __i=0; __i<Kdelt.at(k_).size(); __i++)
							{
								long prefact = Cdelt.at(k_).at(__i);

								if(prefact==1){ if(__i!=0) cout<<" + ";
								}else if(prefact==-1){ cout<<" - ";
								}else if(prefact<0){
									cout<<" - "<<fabs(prefact);
								}else{
									cout<<" + "<<fabs(prefact);
								}

								for(int __k=0; __k<Kdelt.at(k_).at(__i).size()-1; __k += 2) 
								{
									cout<<" Ç_";
									cout<<Kdelt.at(k_).at(__i).getD(__k);

									if(Kdelt.at(k_).at(__i).getD(__k+1)>1)
									{
										cout<<"^";
										cout<<Kdelt.at(k_).at(__i).getD(__k+1);
									}
								}
							}
							cout<<" )";


							vector<int> _mth;

							_mth.push_back(fvv.at(j_).getD(0)); 
							_mth.push_back(fvv.at(k_).getD(0)); 

							// call central moment error estimation class
							CentErrorN cerr(_mth);

							//calculate the expression of error 
							cerr.CalcCentError("");
							cout<<" (";
							cerr.printE();
							cout<<") ";



							//cout<<" ∑(";
							//fvv.at(j_).printD();
							//cout<<",";
							//fvv.at(k_).printD();
							//cout<<")";

						}else{
							cout<<"2 ";


							cout<<"(";
							for(size_t __i=0; __i<Kdelt.at(k_).size(); __i++)
							{
								long prefact = Cdelt.at(k_).at(__i);

								if(prefact==1){ if(__i!=0) cout<<" + ";
								}else if(prefact==-1){ cout<<" - ";
								}else if(prefact<0){
									cout<<" - "<<fabs(prefact);
								}else{
									cout<<" + "<<fabs(prefact);
								}

								for(int __k=0; __k<Kdelt.at(k_).at(__i).size()-1; __k += 2) 
								{
									cout<<" Ç_";
									cout<<Kdelt.at(k_).at(__i).getD(__k);

									if(Kdelt.at(k_).at(__i).getD(__k+1)>1)
									{
										cout<<"^";
										cout<<Kdelt.at(k_).at(__i).getD(__k+1);
									}
								}
							}
							cout<<" )";

							if(j_!=0){
								cout<<"(";
								for(size_t __i=0; __i<Kdelt.at(j_).size(); __i++)
								{
									long prefact = Cdelt.at(j_).at(__i);

									if(prefact==1){ if(__i!=0) cout<<" + ";
									}else if(prefact==-1){ cout<<" - ";
									}else if(prefact<0){
										cout<<" - "<<fabs(prefact);
									}else{
										cout<<" + "<<fabs(prefact);
									}

									for(int __k=0; __k<Kdelt.at(j_).at(__i).size()-1; __k += 2) 
									{
										cout<<" Ç_";
										cout<<Kdelt.at(j_).at(__i).getD(__k);

										if(Kdelt.at(j_).at(__i).getD(__k+1)>1)
										{
											cout<<"^";
											cout<<Kdelt.at(j_).at(__i).getD(__k+1);
										}
									}
								}
								cout<<" )";
							}

							vector<int> _mth;

							_mth.push_back(fvv.at(k_).getD(0)); 
							_mth.push_back(fvv.at(j_).getD(0)); 

							// call central moment error estimation class
							CentErrorN cerr(_mth);

							//calculate the expression of error 
							cerr.CalcCentError("");
							cout<<" (";
							cerr.printE();
							cout<<") ";

							//cout<<" ∑(";
							//fvv.at(k_).printD();
							//cout<<",";
							//fvv.at(j_).printD();
							//cout<<")";
						}
					}


					if(j_==fvv.size()-1 && k_==fvv.size()-1){
						cout<<"";
					}else{
						cout<<" + ";
					}
				}
				fvv.at(j_) *= 0;
			}
			cout<<endl;

*/
			Kdelt__ = Kdelt;
			Cdelt__ = Cdelt;
			fvv__ = fvv;

			vector<vector<FactVec > >().swap(Kdelt);
			vector<vector<long > >().swap(Cdelt);
			vector<FactVec >().swap(fvv);
		}



		void CalcDeriv(const int mdim__, FactVec Cnn, vector<FactVec>& Kdelt_, vector<long>& Cdelt_, string FORMAT_="")
		{
			//if(CALCHECK!=1) this->Cumulant2Moment(mdim__,"NOEXPS");
			nthC.clear();
			coeffC.clear();
			this->Cumulant2Moment(mdim__,"NOEXPS");

			vector<FactVec> Kdelt;
			vector<long> Cdelt;

			size_t i__=nthC.size()-1;

			//for(size_t __j=0; __j<fvv.size(); __j++)
			//{

			if(FORMAT_!="NOEXPS" && FORMAT_!="NO"){
				cout<<"∂K_"<<i__+1<<"/∂Ç_";
				Cnn.printD();
				cout<<" = ";
			}
			for(size_t j__=0; j__<nthC.at(i__).size(); j__++)
			{

				bool CHECK=1;

				for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
				{
					if(nthC.at(i__).at(j__).at(k__).getD(0)==1) CHECK *=0;
				}

				if(CHECK!=0){

					bool CHECK2=0;

					for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
					{
						if(Cnn==nthC.at(i__).at(j__).at(k__))
						{

							long prefact = coeffC.at(i__).at(j__) * nthC.at(i__).at(j__).at(k__).expo();
							coeffC.at(i__).at(j__)=prefact;

							nthC.at(i__).at(j__).at(k__).setE(nthC.at(i__).at(j__).at(k__).expo()-1); //expo()-1

							CHECK2 += 1;
						}
					}

					coeffC.at(i__).at(j__) *= CHECK2;

				}


			}


			vector<FactVec >temp_Kdelt;
			vector<long >temp_Cdelt;
			for(size_t j__=0; j__<nthC.at(i__).size(); j__++)
			{
				bool CHECK=1;

				for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
				{
					if(nthC.at(i__).at(j__).at(k__).getD(0)==1) CHECK *=0;
				}

				if(CHECK!=0 /*&& coeffC.at(i__).at(j__)!=0*/){

					long prefact = coeffC.at(i__).at(j__);

					//cout<<"<";
					//if(prefact<0){
					//	cout<<" - "<<fabs(prefact);
					//}else{
					//	cout<<" + "<<fabs(prefact);
					//}			
					//cout<<">";
					vector<int> vkdelt;
					for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
					{
						//	if(nthC.at(i__).at(j__).at(k__).expo()>0 /*&& coeffC.at(i__).at(j__)!=0*/){
						//cout<<" Ç_";
						//nthC.at(i__).at(j__).at(k__).printD();
						//if(nthC.at(i__).at(j__).at(k__).expo()>1) cout<<"^"<<nthC.at(i__).at(j__).at(k__).expo();
						vkdelt.push_back(nthC.at(i__).at(j__).at(k__).getD(0));
						vkdelt.push_back(nthC.at(i__).at(j__).at(k__).expo());
						//	}
					}

					FactVec temp_temp_Kdelt;
					temp_temp_Kdelt.setV(vkdelt);

					temp_Kdelt.push_back(temp_temp_Kdelt);
					temp_Cdelt.push_back(prefact);

				}


			}

			//Kdelt.push_back(temp_Kdelt);
			//Cdelt.push_back(temp_Cdelt);

			//come back
			for(size_t j_=0; j_<temp_Kdelt.size(); j_++)
			{
				if(temp_Cdelt.at(j_) ==0 ) temp_Kdelt.at(j_) *= 0;
			}

			for(size_t j_=0; j_<temp_Kdelt.size(); j_++){
				for(size_t k_=j_+1; k_<temp_Kdelt.size(); k_++) {
					if(temp_Kdelt.at(j_) == temp_Kdelt.at(k_)) { 
						temp_Kdelt.at(k_) *= 0; //.setE(0);
						temp_Cdelt.at(j_) += temp_Cdelt.at(k_); 
					}
				}
			}



			for(size_t j_=0; j_<temp_Kdelt.size(); j_++) 
			{
				//if(temp_Kdelt.at(j_).sumD() != 0)
				if(temp_Kdelt.at(j_).sumD() != 0)
				{
					Kdelt.push_back(temp_Kdelt.at(j_));
					Cdelt.push_back(temp_Cdelt.at(j_));
				}
			}


			vector<FactVec>().swap(temp_Kdelt);
			vector<long>().swap(temp_Cdelt);


			//cout<<Kdelt.size()<<" "<<temp_Kdelt.size()<<endl; 



			//int j_=0;
			if(FORMAT_!="NO"){
				if(Kdelt.size()==0){ cout<<"";
				}else{cout<<"(";
					for(size_t __i=0; __i<Kdelt.size(); __i++)
					{
						long prefact = Cdelt.at(__i);

						//if(prefact==1){ if(__i!=0) cout<<" + ";
						//}else if(prefact==-1){ cout<<" - ";
						//}else 

						if(prefact<0){
							cout<<" - "<<fabs(prefact);
						}else{
							cout<<" + "<<fabs(prefact);
						}

						for(int __k=0; __k<Kdelt.at(__i).size()-1; __k += 2) 
						{
							if(Kdelt.at(__i).getD(__k+1)>0){
								cout<<" Ç_";
								cout<<Kdelt.at(__i).getD(__k);

								if(Kdelt.at(__i).getD(__k+1)>1)
								{
									cout<<"^";
									cout<<Kdelt.at(__i).getD(__k+1);
								}
							}
						}
					}
					cout<<")";
				}
				if(FORMAT_ !="NOEXPS")cout<<endl;
			}

			nthC.clear();
			coeffC.clear();
			this->Cumulant2Moment(mdim__,"NOEXPS");
			//}
			Kdelt_ = Kdelt;
			Cdelt_ = Cdelt;

			vector<FactVec>().swap(Kdelt);
			vector<long>().swap(Cdelt);
		}





};




//Strongly Intensive Cumulant Class
//Observable references from : 
class SICumulantVec : public CentVecN
{
	private:
		int mdim;
		int nsize;

		vector<vector<FactVec> > Kappa;
		vector<vector<FactVec> > Mumn;
		vector<vector<long> > cmn;

		vector<vector<vector<FactVec> > >nthC;
		vector<vector<long> > coeffC;
		vector<vector<short> > MupowC;

		bool CALCHECK;
	public:
		SICumulantVec(const size_t mdim_=0) : CentVecN()
	{
		mdim=mdim_;
		nsize=0;
		CALCHECK=0;
		Kappa.resize(mdim);
		Mumn.resize(mdim);
		cmn.resize(mdim);
		MOMENT="SICumulant";
	}

		~SICumulantVec()
		{
			vector<vector<FactVec> >().swap(Kappa);
			vector<vector<FactVec> >().swap(Mumn);
			vector<vector<long> >().swap(cmn);

			vector<vector<vector<FactVec> > >().swap(nthC);
			vector<vector<long> >().swap(coeffC);
		}


		vector<vector<vector<FactVec> > > getnthC(){return nthC;}
		vector<vector<long> > getcoeffC() {return coeffC;}
		vector<vector<short> > getMupowC() {return MupowC;}

		bool calcheck() const {return CALCHECK;}

		void SICumulant2SICumulant(const int mdim__)
		{
			//
			//
			//		     1                r-1  
			// K_r = K_{r,0} = ------  ( µ_{r,0} - ∑ µ_{r-i,1} K_i )
			//		   µ_{0,1}	       1
			//		   
			//	 µ_{1,0}	   
			//K_1 = --------
			//       µ_{0,1}

			cout<<"Func to convert cumulant K_m to cumulant K_m"<<endl; 

			for(int k_=1; k_<=mdim__; k_++)
			{
				vector<FactVec> Kappa_;
				vector<FactVec> Mumn_;
				vector<long> cmn_;

				for(int i_=1; i_<=k_-1; i_++)
				{
					FactVec TEMP_Kappa_(2);
					TEMP_Kappa_.setD(0,i_);
					TEMP_Kappa_.setD(1,0);
					FactVec TEMP_Mumn_(2);
					TEMP_Mumn_.setD(0,k_-i_);
					TEMP_Mumn_.setD(1,1);

					cmn_.push_back(-1);//-calccomb(k_-1,i_-1));
					Kappa_.push_back(TEMP_Kappa_);
					Mumn_.push_back(TEMP_Mumn_);
				}


				cout<<"K_"<<k_<<"= (µ_"<<k_<<0;
				//if(cmn_.size()!=0)cout<<" + ";
				for(size_t l__=0; l__<cmn_.size(); l__++)
				{

					if(cmn_.at(l__)==1){ if(l__!=0) cout<<" + ";
					}else if(cmn_.at(l__)==-1){ cout<<" - ";
					}else if(cmn_.at(l__)<0){
						cout<<" - "<<fabs(cmn_.at(l__));
					}else{
						cout<<" + "<<fabs(cmn_.at(l__));
					}

					cout<<" K_";
					Kappa_.at(l__).printD();
					cout<<" µ_";
					Mumn_.at(l__).printD();
					if(l__<cmn_.size()-1)cout<<" ";
				}
				cout<<")/µ_01";


				cout<<endl;

				cmn_.clear();
				Kappa_.clear();
				Mumn_.clear();
			}
		}

		void SICumulant2Moment(const int mdim__, string FORMAT="")
		{
			cout<<"-------------------------------------------------------------------------------"<<endl;


			if(FORMAT!="NOEXPS")cout<<"Func to convert strongly intesive cumulant K_m to moments µ_m "<<endl; 

			for(int k_=1; k_<=mdim__; k_++)
			{

				vector<FactVec> Kappa_;
				vector<FactVec> Mumn_;
				vector<long> cmn_;
				vector<short> Mupow_;

				for(int i_=1; i_<=k_-1; i_++)
				{

					FactVec TEMP_Kappa_(2);
					TEMP_Kappa_.setD(0,i_);
					TEMP_Kappa_.setD(1,0);

					FactVec TEMP_Mumn_(2);
					TEMP_Mumn_.setD(0,k_-i_);
					TEMP_Mumn_.setD(1,1);
					TEMP_Mumn_.setE(1);


					cmn_.push_back(-calccomb(k_-1,i_-1));
					Kappa_.push_back(TEMP_Kappa_);
					Mumn_.push_back(TEMP_Mumn_);
					Mupow_.push_back(1);
				}

				vector<vector<FactVec> > temp_nthC;
				vector<long > temp_coeffC;
				vector<short> temp_MupowC;


				FactVec TEMP_Mumn__(2);
				TEMP_Mumn__.setD(0,k_);
				TEMP_Mumn__.setD(1,0);
				TEMP_Mumn__.setE(1);


				vector<FactVec > TEMP_TEMP_Mumn__;
				TEMP_TEMP_Mumn__.push_back(TEMP_Mumn__);
				temp_nthC.push_back(TEMP_TEMP_Mumn__);

				temp_coeffC.push_back(1.);
				temp_MupowC.push_back(1);


				for(size_t i__=0; i__<Mumn_.size(); i__++)
				{

					for(size_t i=0; i<nthC.at(Kappa_.at(i__).getD(0)-1).size(); i++)
					{
						vector<FactVec> temp_temp_nthC;

						for(size_t ll_=0; ll_<nthC.at(Kappa_.at(i__).getD(0)-1).at(i).size(); ll_++)
						{
							temp_temp_nthC.push_back(nthC.at(Kappa_.at(i__).getD(0)-1).at(i).at(ll_));
						}
						temp_temp_nthC.push_back(Mumn_.at(i__));

						temp_nthC.push_back(temp_temp_nthC);
						temp_coeffC.push_back(cmn_.at(i__)*coeffC.at(Kappa_.at(i__).getD(0)-1).at(i));

						//cout<<Mupow_.at(i__)<<" "<<MupowC.at(Kappa_.at(i__).getD(0)-1).at(i);
						temp_MupowC.push_back(Mupow_.at(i__) + MupowC.at(Kappa_.at(i__).getD(0)-1).at(i));

					}
				}

				nthC.push_back(temp_nthC);
				coeffC.push_back(temp_coeffC);
				MupowC.push_back(temp_MupowC);

				vector<long>().swap(cmn_);
				vector<short>().swap(Mupow_);
				vector<FactVec>().swap(Kappa_);
				vector<FactVec>().swap(Mumn_);

				vector<vector<FactVec> >().swap(temp_nthC);
				vector<long >().swap(temp_coeffC);
				vector<short >().swap(temp_MupowC);
				vector<FactVec >().swap(TEMP_TEMP_Mumn__);
			}


			//squeeze example
			for(size_t l__=0; l__<nthC.size(); l__++)
			{
				for(size_t i__=0; i__<nthC.at(l__).size(); i__++)
				{
					for(size_t j__=0; j__<nthC.at(l__).at(i__).size(); j__++)
					{
						if(nthC.at(l__).at(i__).at(j__).expo()==0) continue;
						for(size_t k__=j__+1; k__<nthC.at(l__).at(i__).size(); k__++)
						{

							if(nthC.at(l__).at(i__).at(j__) == nthC.at(l__).at(i__).at(k__))
							{
								nthC.at(l__).at(i__).at(j__).addE(nthC.at(l__).at(i__).at(k__).expo());
								nthC.at(l__).at(i__).at(k__).setE(0);
							}
						}
					}
				}
			}



			CALCHECK=1;
		}


		void printL()//const int mdim__, string FORMAT="")
		{
			//printing the expression
			for(size_t i__=0; i__<nthC.size(); i__++)
			{
				//cout<<" K_"<<i__+1<<"= ";
				cout<<" $kappa_{"<<i__+1<<"} &=& ";
				for(size_t j__=0; j__<nthC.at(i__).size(); j__++)
				{
					long prefact = coeffC.at(i__).at(j__);

					if(prefact==1){ if(j__!=0) cout<<" + ";
					}else if(prefact==-1){ cout<<" - ";
					}else if(prefact<0){
						cout<<" - "<<fabs(prefact);
					}else{
						cout<<" + "<<fabs(prefact);
					}

					cout<<"$frac{";
					for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
					{
						if(nthC.at(i__).at(j__).at(k__).expo()>0){
							//cout<<" µ_";
							//cout<<" $mu_{";
							//nthC.at(i__).at(j__).at(k__).printD();

							vector<int> _mth;
							for(size_t __i=0;__i<nthC.at(i__).at(j__).at(k__).size();__i++) 
							{
								//cout<<" mth "<<nthC.at(i__).at(j__).at(k__).getD(__i)<<endl;
								_mth.push_back(nthC.at(i__).at(j__).at(k__).getD(__i));
							}

							MomVecN mnt_;

							cout<<"(";
							mnt_.Mom2FactN(_mth);//,"NOEXPS");
							mnt_.printD("latex");
							cout<<")";

							if(nthC.at(i__).at(j__).at(k__).expo()>1){
								cout<<"^";
								cout<<nthC.at(i__).at(j__).at(k__).expo();
							}
						}
					}
					cout<<"}{{$rm f}_{01}";
					if(MupowC.at(i__).at(j__)>1){
						cout<<"^"<<MupowC.at(i__).at(j__)<<"}";
					}else{
						cout<<"}";
					}

					if(j__!=0 && j__%2==0){
						cout<<" $$";
						cout<<endl;
						cout<<" && ";
					}
				}
				cout<<" $$";
				cout<<endl;
			}
			cout<<endl;

		}


		void printD()//const int mdim__, string FORMAT="")
		{

			//printing the expression
			for(size_t i__=0; i__<nthC.size(); i__++)
			{
				cout<<" (SI)K_"<<i__+1<<"= ";
				//cout<<" $kappa_{"<<i__+1<<"} &=& ";
				for(size_t j__=0; j__<nthC.at(i__).size(); j__++)
				{
					long prefact = coeffC.at(i__).at(j__);

					if(prefact==1){ if(j__!=0) cout<<" + ";
					}else if(prefact==-1){ cout<<" - ";
					}else if(prefact<0){
						cout<<" - "<<fabs(prefact);
					}else{
						cout<<" + "<<fabs(prefact);
					}

					//cout<<"$frac{";
					for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
					{
						if(nthC.at(i__).at(j__).at(k__).expo()>0){
							cout<<" µ";
							//cout<<" $mu_{";
							nthC.at(i__).at(j__).at(k__).printD();
							//cout<<"}";
							if(nthC.at(i__).at(j__).at(k__).expo()>1){
								cout<<"^";
								cout<<nthC.at(i__).at(j__).at(k__).expo();
							}
						}
					}
					//cout<<"}{$mu_{01}";
					//cout<<"}{$mu_{01}";
					cout<<" /µ01";
					if(MupowC.at(i__).at(j__)>1){
						cout<<"^"<<MupowC.at(i__).at(j__);//<<"}";
					}else{
						//cout<<"}";
					}

					//if(j__!=0 && j__%4==0){
					//	cout<<" $$";
					//	cout<<endl;
					//	cout<<" && ";
					//}
				}
				//cout<<" $$";
				cout<<endl;
			}
			cout<<endl;




			//printing the expression
			for(size_t i__=0; i__<nthC.size(); i__++)
			{
				cout<<" (SI)K_"<<i__+1<<"= ";
				//cout<<" $kappa_{"<<i__+1<<"} &=& ";
				for(size_t j__=0; j__<nthC.at(i__).size(); j__++)
				{
					long prefact = coeffC.at(i__).at(j__);

					if(prefact==1){ if(j__!=0) cout<<" + ";
					}else if(prefact==-1){ cout<<" - ";
					}else if(prefact<0){
						cout<<" - "<<fabs(prefact);
					}else{
						cout<<" + "<<fabs(prefact);
					}

					//cout<<"$frac{";

					for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
					{
						if(nthC.at(i__).at(j__).at(k__).expo()>0){
							//cout<<" µ_";
							//cout<<" $mu_{";
							//nthC.at(i__).at(j__).at(k__).printD();

							vector<int> _mth;
							for(size_t __i=0;__i<nthC.at(i__).at(j__).at(k__).size();__i++) 
							{
								//cout<<" mth "<<nthC.at(i__).at(j__).at(k__).getD(__i)<<endl;
								_mth.push_back(nthC.at(i__).at(j__).at(k__).getD(__i));
							}

							MomVecN mnt_;

							cout<<"(";
							mnt_.Mom2FactN(_mth);//,"NOEXPS");
							mnt_.printD("NOEXPS");
							cout<<")";

							if(nthC.at(i__).at(j__).at(k__).expo()>1){
								cout<<"^";
								cout<<nthC.at(i__).at(j__).at(k__).expo();
							}
						}
					}

					//cout<<"}{{$rm f}_{01}";
					cout<<"/f01";

					if(MupowC.at(i__).at(j__)>1){
						cout<<"^"<<MupowC.at(i__).at(j__);//<<"}";
					}else{
						//cout<<"}";
					}

					//if(j__!=0 && j__%2==0){
					//	cout<<" $$";
					//	cout<<endl;
					//	cout<<" && ";
					//}
				}
				//	cout<<" $$";
				cout<<endl;
			}
			cout<<endl;

		}


		///void SetFmnabN(const vector<int> mthN_, const vector<int> abinN_)
		void printD(const vector<int> abinN_)//const int mdim__, string FORMAT="")
		{
			//printing the expression
			for(size_t i__=0; i__<nthC.size(); i__++)
			{
				cout<<" (SI)K_"<<i__+1<<"= ";
				//cout<<" $kappa_{"<<i__+1<<"} &=& ";
				for(size_t j__=0; j__<nthC.at(i__).size(); j__++)
				{
					long prefact = coeffC.at(i__).at(j__);

					if(prefact==1){ if(j__!=0) cout<<" + ";
					}else if(prefact==-1){ cout<<" - ";
					}else if(prefact<0){
						cout<<" - "<<fabs(prefact);
					}else{
						cout<<" + "<<fabs(prefact);
					}

					//cout<<"$frac{";

					for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
					{
						if(nthC.at(i__).at(j__).at(k__).expo()>0){
							//cout<<" µ_";
							//cout<<" $mu_{";
							//nthC.at(i__).at(j__).at(k__).printD();

							vector<int> _mth;
							for(size_t __i=0;__i<nthC.at(i__).at(j__).at(k__).size();__i++) 
							{
								//cout<<" mth "<<nthC.at(i__).at(j__).at(k__).getD(__i)<<endl;
								_mth.push_back(nthC.at(i__).at(j__).at(k__).getD(__i));
							}

							MomVecN mnt_;

							cout<<"(";
							//mnt_.Mom2FactN(_mth,"NOEXPS");
							mnt_.Mom2FactN(_mth);//,abinN_,"NOEXPS");
							mnt_.printD(abinN_,"NOEXPS");
							cout<<")";

							if(nthC.at(i__).at(j__).at(k__).expo()>1){
								cout<<"^";
								cout<<nthC.at(i__).at(j__).at(k__).expo();
							}

							vector<int>().swap(_mth);
						}
					}

					cout<<"/";

					vector<int> _mth;
					//		for(size_t __i=0;__i<nthC.at(i__).at(j__).at(k__).size();__i++) 
					//		{
					//cout<<" mth "<<nthC.at(i__).at(j__).at(k__).getD(__i)<<endl;
					_mth.push_back(0);//nthC.at(i__).at(j__).at(k__).getD(__i));
					_mth.push_back(1);//nthC.at(i__).at(j__).at(k__).getD(__i));
					//		}

					MomVecN mnt__;

					cout<<"(";
					//mnt_.Mom2FactN(_mth,"NOEXPS");
					mnt__.Mom2FactN(_mth);//,abinN_,"NOEXPS");
					mnt__.printD(abinN_,"NOEXPS");
					cout<<")";


					//cout<<"}{{$rm f}_{01}";
					//cout<<"/f01";

					if(MupowC.at(i__).at(j__)>1){
						cout<<"^"<<MupowC.at(i__).at(j__);//<<"}";
					}else{
						//cout<<"}";
					}

					//if(j__!=0 && j__%2==0){
					//	cout<<" $$";
					//	cout<<endl;
					//	cout<<" && ";
					//}
				}
				//	cout<<" $$";
				cout<<endl;
				cout<<endl;
				cout<<endl;
			}
			cout<<endl;

		}

		//Overriding the following function from the base class
		vector<FactVec> CalcFactVec(const int mdim__) //, string FORMAT="")
		{
			vector<FactVec> fvv;

			if(CALCHECK!=1) this->SICumulant2Moment(mdim__,"NOEXPS");

			size_t i__=nthC.size()-1;

			//	if(i__==0){
			//		//cout<<"f10 - f01"<<endl;
			FactVec F10(2);
			F10.setD(0,1);
			FactVec F01(2);
			F01.setD(1,1);
			fvv.push_back(F10);
			fvv.push_back(F01);
			//continue;
			//	}

			for(size_t j__=0; j__<nthC.at(i__).size(); j__++)
			{

				vector<FactVec> temp_fvv;

				//bool CHECK=1;

				//for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
				//{
				//	if(nthC.at(i__).at(j__).at(k__).getD(0)==1) CHECK *=0;
				//}

				//if(CHECK!=0){
				//long prefact = coeffC.at(i__).at(j__);

				for(size_t k__=0; k__<nthC.at(i__).at(j__).size(); k__++)
				{
					vector<int> _mth;
					for(size_t __i=0;__i<nthC.at(i__).at(j__).at(k__).size();__i++) 
					{
						//cout<<" mth "<<nthC.at(i__).at(j__).at(k__).getD(__i)<<endl;
						_mth.push_back(nthC.at(i__).at(j__).at(k__).getD(__i));
					}

					MomVecN mnt_;

					//cout<<"(";
					//mnt_.Mom2FactN(_mth,"NOEXPS");
					//cout<<")";

					temp_fvv=mnt_.CalcFactVec(_mth);

					//CentVecN cnt_;
					//vector<int> mthC;
					//mthC.push_back(nthC.at(i__).at(j__).at(k__).getD(0));

					//cnt_.Corr2FactN(mthC,"NOEXPS","");
					//temp_fvv=cnt_.CalcFactVec("");
					vector<int>().swap(_mth);
				}

				//}

				fvv.insert(fvv.end(), temp_fvv.begin(), temp_fvv.end());
				vector<FactVec>().swap(temp_fvv);	
			}

			vector<FactVec> temp_fvv_;

			for(size_t j_=0; j_<fvv.size(); j_++){
				if(fvv.at(j_).sumD()!=0) temp_fvv_.push_back(fvv.at(j_)); 
				for(size_t k_=j_+1; k_<fvv.size(); k_++) {
					if(fvv.at(j_) == fvv.at(k_)) fvv.at(k_) *= 0; //{ fvv.erase(fvv.begin() + k_); }
				}
			}

			fvv.swap(temp_fvv_);

			vector<FactVec>().swap(temp_fvv_);

			cout<<"Factorial Moments to be calculated are :"<<endl;

			for(size_t j_=0; j_<fvv.size(); j_++) {cout<<"F_"; fvv.at(j_).printD(); cout<<endl;}

			return fvv;
		}

};



class CumulantRatio : public CumulantVec, CentVecN 
{
	protected:

		vector<int> order;
		vector<float> expo;

		vector<FactVec> Kappa;

		int maxorder;
	public:
		CumulantRatio(const vector<int> order_, const vector<float> expo_) : CentVecN(0,"CumulantRatio") , CumulantVec()
	{


		for(size_t __l=0; __l<order_.size(); __l++)
		{
			if(order_.at(__l)!=0 && expo_.at(__l)!=0)
			{
				order.push_back(order_.at(__l));
				expo.push_back(expo_.at(__l));
			}
		}

		float intensive =0.;
		maxorder=0;

		for(size_t __l=0; __l<order.size(); __l++)
		{
			intensive += (order.at(__l)*expo.at(__l));
			if(order.at(__l)>maxorder) maxorder=order.at(__l);
		}

//		if(intensive!=0)cout<<"WARNING : Cumulant Ratio is not intensive"<<endl;
		if(intensive!=0)cout << "\033[1;35mWARNING !! Cumulant Ratio is not intensive \033[0m\n"<<" "<<endl;

		Cumulant2Moment(maxorder,"NOEXPS");
	}


		CumulantRatio(const string input) : CentVecN(0,"CumulantRatio") , CumulantVec()
	{
		//std::string input = "K4/K2";  //no space 
		std::stringstream ss;
		ss << input;
		int found;
		std::string temp;

		cout<<"input="<<input<<endl; 

		int counter=0;

		while(std::getline(ss, temp,'/')) {

//			cout<<"Im here "<<ss<<" "<<temp<<endl;
			if(counter==0 && temp.size()>3) 
			{
				std::stringstream uu;
				std::string temp_;
				uu << temp;
				int counter_=0;
				while(std::getline(uu, temp_,'C'))
				{
					int order_ = atoi(temp_.c_str() +0);
					float expo_ =0;
					if(temp_.size()<3){ expo_=1;
					}else{expo_ = atof(temp_.c_str() + 2);
					}
					if(counter_!=0) cout<<" order="<<order_<<" expo="<<expo_<<endl;

					order.push_back(order_);
					expo.push_back(expo_);

					counter_ ++;
				}
			}else{
				int order__ = atoi(temp.c_str() +1);
				float expo__ =0;
				if(temp.size()<3){ expo__=-1;
				}else{expo__ = -atof(temp.c_str() + 3);
				}

				if(counter==0) expo__ *= -1;

				cout<<" order="<<order__<<" expo="<<expo__<<endl;

				order.push_back(order__);
				expo.push_back(expo__);
			}
			counter ++;
		}



		float intensive =0.;
		maxorder=0;

		for(size_t __l=0; __l<order.size(); __l++)
		{
			intensive += (order.at(__l)*expo.at(__l));
			if(order.at(__l)>maxorder) maxorder=order.at(__l);

		}

		//if(intensive!=0)cout<<"WARNING : Cumulant Ratio is not intensive"<<endl;
		if(intensive!=0)cout << "\033[1;35mWARNING !! Cumulant Ratio is not intensive \033[0m\n"<<" "<<endl;

//		cout<<" moment "<<MOMENT<<endl;
		if(maxorder>10){
			cout << "\033[1;35mWARNING !! Cumulant order is too high (Ignore this for Koch Ratio) \033[0m\n"<<" "<<endl;
			//cout<<"Im out "<<maxorder<<endl;
			
		}else{
			Cumulant2Moment(maxorder,"NOEXPS");
		}
	}


		~CumulantRatio()
		{

			vector<int>().swap(order);
			vector<float>().swap(expo);
			vector<FactVec>().swap(Kappa);
		}



		void CumulantRatio2Cumulant()
		{

			Kappa.clear();
			for(size_t __i=0; __i<order.size(); __i++)
			{
				FactVec temp_Kappa(1);
				temp_Kappa.setD(0,order.at(__i));
				temp_Kappa.setE(expo.at(__i));

				Kappa.push_back(temp_Kappa);
			}

			cout<<"   ";

			for(size_t __j=0; __j<order.size(); __j++)
			{
				if(Kappa.at(__j).expo()>0 && Kappa.at(__j).getD(0)!=0)
				{
					cout<<" K_";
					Kappa.at(__j).printD();
					if(Kappa.at(__j).expo()>1) cout<<"^"<<Kappa.at(__j).expo();
				}
			}

			cout<<endl;
			cout<<"k =";
			for(size_t __i=0; __i<order.size(); __i++) cout<<"---";
			cout<<endl;

			cout<<"   ";
			for(size_t __j=0; __j<order.size(); __j++)
			{
				if(Kappa.at(__j).expo()<0 && Kappa.at(__j).getD(0)!=0)
				{
					cout<<" K_";
					Kappa.at(__j).printD();
					if(fabs(Kappa.at(__j).expo())>1) cout<<"^"<<fabs(Kappa.at(__j).expo());
				}

			}
			cout<<endl;

		}



		void CumulantRatio2CorrMoment()
		{

			Kappa.clear();
			for(size_t __i=0; __i<order.size(); __i++)
			{
				FactVec temp_Kappa(1);
				temp_Kappa.setD(0,order.at(__i));
				temp_Kappa.setE(expo.at(__i));

				Kappa.push_back(temp_Kappa);
			}

			cout<<"   ";

			for(size_t __j=0; __j<order.size(); __j++)
			{
				if(Kappa.at(__j).expo()>0 && Kappa.at(__j).getD(0)!=0)
				{
					cout<<"(";
					Cumulant2CorrMoment(Kappa.at(__j).getD(0),"NOEXPS");
					cout<<")";
					if(Kappa.at(__j).expo()>1) cout<<"^"<<Kappa.at(__j).expo();
				}
			}

			cout<<endl;
			cout<<"k =";
			for(size_t __i=0; __i<getnthC().at(Kappa.at(0).getD(0)-1).size(); __i++) cout<<"----";
			cout<<endl;

			cout<<"   ";
			for(size_t __j=0; __j<order.size(); __j++)
			{
				if(Kappa.at(__j).expo()<0 && Kappa.at(__j).getD(0)!=0)
				{
					cout<<"(";
					Cumulant2CorrMoment(Kappa.at(__j).getD(0),"NOEXPS");
					cout<<")";
					if(fabs(Kappa.at(__j).expo())>1) cout<<"^"<<fabs(Kappa.at(__j).expo());
				}

			}
			cout<<endl;
		}

		void CR2CM(vector<int> order_, vector<float> expo_, float prefact_)
		{

			Kappa.clear();
			for(size_t __i=0; __i<order_.size(); __i++)
			{
				FactVec temp_Kappa(1);
				temp_Kappa.setD(0,order_.at(__i));
				temp_Kappa.setE(expo_.at(__i));

				Kappa.push_back(temp_Kappa);
			}

			cout<<"("<<prefact_<<")";

			for(size_t __j=0; __j<order_.size(); __j++)
			{
				if(Kappa.at(__j).expo()>0 && Kappa.at(__j).getD(0)!=0)
				{
					cout<<"(";
					//Cumulant2Moment(maxorder,"NOEXPS");
					Cumulant2CorrMoment(Kappa.at(__j).getD(0),"NOEXPS");
					cout<<")";
					if(Kappa.at(__j).expo()>1) cout<<"^"<<Kappa.at(__j).expo();
				}
			}

			cout<<"/("; 
			//cout<<endl;
			//cout<<prefact_;
			//for(size_t __i=0; __i<getnthC().at(Kappa.at(0).getD(0)-1).size(); __i++) cout<<"----";
			//cout<<endl;

			for(size_t __j=0; __j<order_.size(); __j++)
			{
				//CALCHECK=0;
				if(Kappa.at(__j).expo()<0 && Kappa.at(__j).getD(0)!=0)
				{
					cout<<"(";
					Cumulant2CorrMoment(Kappa.at(__j).getD(0),"NOEXPS");
					cout<<")";
					if(fabs(Kappa.at(__j).expo())>1) cout<<"^"<<fabs(Kappa.at(__j).expo());
				}

			}
			cout<<") ";
			//	cout<<endl;
		}

		vector<FactVec> CalcCentVec() //, string FORMAT="")
		{
			vector<FactVec> fvv;

			for(size_t __i=0; __i<order.size(); __i++)
			{

				size_t i__=order.at(__i)-1;

				for(size_t j__=0; j__<getnthC().at(i__).size(); j__++)
				{

					bool CHECK=1;

					for(size_t k__=0; k__<getnthC().at(i__).at(j__).size(); k__++)
					{
						if(getnthC().at(i__).at(j__).at(k__).getD(0)==1) CHECK *=0;
					}

					if(CHECK!=0){

						for(size_t k__=0; k__<getnthC().at(i__).at(j__).size(); k__++)
						{
							fvv.push_back(getnthC().at(i__).at(j__).at(k__));
						}
					}
				}

			}

			for(size_t j_=0; j_<fvv.size(); j_++){
				for(size_t k_=j_+1; k_<fvv.size(); k_++) {
					if(fvv.at(j_) == fvv.at(k_)) { fvv.at(k_) *= 0;}
				}
			}

			vector<FactVec> temp_fvv;

			for(size_t j_=0; j_<fvv.size(); j_++) {if(fvv.at(j_).sumD() != 0) temp_fvv.push_back(fvv.at(j_));}

			fvv.swap(temp_fvv);

			//cout<<"Central Moments to be calculated are :"<<endl;

			//for(size_t j_=0; j_<fvv.size(); j_++) {cout<<"Ç_"; fvv.at(j_).printD(); cout<<endl;}
			return fvv;
		}

		void CalcDkDCn(FactVec Cnn_)
		{

			//vector<FactVec> fvv = CalcFactVec();

			vector<int> order_ = order;
			vector<float> expo_;

			//for(size_t j_=0; j_<fvv.size(); j_++)
			//{
			//CumulantVec clt;
			//clt.CalcVariance(order.at(__i));
			//cout<<"∂k/∂C_"<<fvv.at(j_).getD(0)<<"=";
			for(size_t __i=0; __i<order.size(); __i++)
			{

				//FactVec Cnn_ = fvv.at(j_);

				vector<FactVec> Kdelt;
				vector<long> Cdelt;

				//	cout<<expo.at(__i)<<" ";

				CumulantVec clt;
				clt.CalcDeriv(order.at(__i), Cnn_, Kdelt, Cdelt, "NO");

				if(__i!=0 && Kdelt.size()!=0 )cout<<" + ";

				if(Kdelt.size()!=0){
					expo_ = expo;
					expo_.at(__i) -= 1;

					CR2CM(order_, expo_, expo.at(__i));
				}

				clt.CalcDeriv(order.at(__i), Cnn_, Kdelt, Cdelt, "NOEXPS");

			}
			//cout<<endl;
			//cout<<endl;
			//}



		}

		void CalcVariance()
		{
			this->CalcPDeriv();
		}

		void CalcPDeriv()
		{
			vector<FactVec> fvv = CalcCentVec();

			//vector<int> order_ = order;
			//vector<float> expo_;

			//	for(size_t j_=0; j_<fvv.size(); j_++)
			//	{
			//		CalcDkDCn(fvv.at(j_));
			//		cout<<endl;
			//	}

			cout<<endl;
			cout<<"Var(k)=";
			for(size_t j_=0; j_<fvv.size(); j_++){
				for(size_t k_=0; k_<fvv.size(); k_++){

					if(fvv.at(k_).sumD()==0) continue;

					if(fvv.at(j_) == fvv.at(k_)) 
					{ 
						//if(j_!=0){	
						cout<<"(";
						CalcDkDCn(fvv.at(j_));
						cout<<" )^2 ";
						//}

						vector<int> _mth;

						_mth.push_back(fvv.at(j_).getD(0)); 
						_mth.push_back(fvv.at(k_).getD(0)); 

						// call central moment error estimation class
						CentErrorN cerr(_mth);

						//calculate the expression of error 
						cerr.CalcCentError("");
						cout<<" (";
						cerr.printE();
						cout<<") ";

						//cout<<" ∑(";
						//fvv.at(j_).printD();
						//cout<<",";
						//fvv.at(k_).printD();
						//cout<<")";

					}else{
						if(fvv.at(j_).sumD()>fvv.at(k_).sumD()){

							cout<<"2 ";

							//if(j_!=0){
							cout<<"(";
							CalcDkDCn(fvv.at(j_));
							cout<<" )";
							//}

							cout<<"(";
							CalcDkDCn(fvv.at(k_));
							cout<<" )";


							vector<int> _mth;

							_mth.push_back(fvv.at(j_).getD(0)); 
							_mth.push_back(fvv.at(k_).getD(0)); 

							// call central moment error estimation class
							CentErrorN cerr(_mth);

							//calculate the expression of error 
							cerr.CalcCentError("");
							cout<<" (";
							cerr.printE();
							cout<<") ";



							//cout<<" ∑(";
							//fvv.at(j_).printD();
							//cout<<",";
							//fvv.at(k_).printD();
							//cout<<")";

						}else{
							cout<<"2 ";


							cout<<"(";
							CalcDkDCn(fvv.at(k_));
							cout<<" )";

							//if(j_!=0){
							cout<<"(";
							CalcDkDCn(fvv.at(j_));
							cout<<" )";
							//}

							vector<int> _mth;

							_mth.push_back(fvv.at(k_).getD(0)); 
							_mth.push_back(fvv.at(j_).getD(0)); 

							// call central moment error estimation class
							CentErrorN cerr(_mth);

							//calculate the expression of error 
							cerr.CalcCentError("");
							cout<<" (";
							cerr.printE();
							cout<<") ";

							//cout<<" ∑(";
							//fvv.at(k_).printD();
							//cout<<",";
							//fvv.at(j_).printD();
							//cout<<")";
						}
					}


					if(j_==fvv.size()-1 && k_==fvv.size()-1){
						cout<<"";
					}else{
						cout<<" + ";
					}
				}
				fvv.at(j_) *= 0;
			}
			cout<<endl;



		}
};




class KochRatio : public CumulantRatio
{
	protected:
		vector<FactVec> kratio;
		FactVec maxvec;
		vector<FactVec> fvv_ratio;

	public:
		KochRatio(string input_) : CumulantRatio(input_)
		{

			CumulantRatio2Cumulant();

			vector<int> kvec;
			kvec.push_back(0);
			kvec.push_back(0);
			maxvec.setV(kvec);

			for(size_t _i=0;_i<Kappa.size();_i++)
			{

	//			Kappa.at(_i).printD();
	//			cout<<endl;
	//			cout<<"size "<<Kappa.at(_i).size()<<endl;


				FactVec temp_ratio(2);

				(Kappa.at(_i).getD(0)>=10) ? temp_ratio.setD(1,(Kappa.at(_i).getD(0)%10)) : temp_ratio.setD(1,Kappa.at(_i).getD(0));
				if(Kappa.at(_i).getD(0)>=10) temp_ratio.setD(0,int(Kappa.at(_i).getD(0)/10));


				kratio.push_back(temp_ratio);
				maxvec += temp_ratio;

//				temp_ratio.printD();
//				cout<<endl;
			}

  			if(kratio.size()<2){
                                cout << "\033[1;31mERROR !! Your input is not a valid Koch Ratio \033[0m\n"<<" "<<endl;
                                abort();
                        }
//			maxvec.printD();


			vector<int> mth__;     // order 02

			mth__.push_back(maxvec.getD(0));   // order of 1st positive variable
			mth__.push_back(maxvec.getD(1));   // order of 2nd positive variable

			mth__.push_back(maxvec.getD(0));    // order of 1st negative variable
			mth__.push_back(maxvec.getD(1));    // order of 2nd negative variable


			CentErrorN cerr__(mth__);// = new CentErrorN(mth__);  // call central moment error estimation class
			//freopen("formula.out","w",stdout);

			cerr__.CalcCentError("NOEXPS");    //calculate the expression of error

			cout<<"∑=";
			cerr__.printE();    //print to see the error expression
			cerr__.printD();   //print the efficiency corrected error expression


			fvv_ratio=cerr__.CalcFactVec();

		}

		vector<FactVec> CalcFactVec() {return fvv_ratio;}

		vector<FactVec> getkratio(){return kratio;}
};


#endif
