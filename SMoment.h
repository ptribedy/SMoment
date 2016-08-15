/*
  Class to perform e-by-e analysis for multi-variate moments 
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

//v-1.0 modified 14th May, 2015 calss FactmnN and SMomentN added
//v-1.1 modified 31th May, 2015 Nudyn functions added
//v-1.2 modified 3rd July, 2015 weighted sum added, nweight def changed 
//v-1.3 modified 5th July, 2015 error calculation func CalcCorrError added
//v-1.4 modified 18th July,2015 cumulant calculation CalcCumulant added
//v-1.5 modified 4th August, 2015 modified CalcCorrMom added getNevent() & getNweight() included Nihar's func
//v-1.6 modified 15th Sept, 2015 SIC functionality implemented
//v-1.7 modified 30th Dec, 2015 Bug in printing expression inside CalcCorrError fixed  
//v-1.8 modified Aug 11,2016
//v-1.9 current version

#ifndef SMoment_h
#define SMoment_h

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <vector>

#include "FactMom.h"

using namespace std;


class Factmn
{
	private:
		int mdim;
		int ndim;
		int pos(int m__, int n__) const {return (mdim*m__+n__);}

	protected:
		int fsize;
		vector<long double> fmn;

	public:
		Factmn(int mdim_=0, int ndim_=0)
		{
			mdim=mdim_; ndim=ndim_; fsize=0;
			for(int j=0; j<mdim; j++){
				for(int k=0; k<ndim; k++)
					fmn.push_back(0.0L);
				fsize++;
			}
		}

		~Factmn()
		{
			vector<long double>().swap(fmn);
		}

		void setFmn(int i__, long double fmn_){if(i__<fsize) fmn.at(i__) = fmn_;}
		void setFmn(int i__, int j__, long double fmn_){if(i__<mdim && j__<ndim) fmn.at(pos(i__,j__)) = fmn_;}
		void addFmn(int i__, long double fmn_){if(i__<fsize) fmn.at(i__) += fmn_;}
		void addFmn(int i__,int j__, long double fmn_){if(i__<mdim && j__<ndim) fmn.at(pos(i__,j__)) += fmn_;}
		long double getFmn(int i__) const {if(i__<fsize){return fmn.at(i__);}else{return 0;}}
		long double getFmn(int i__, int j__) const {if(i__<mdim && j__<ndim){return fmn.at(pos(i__,j__));}else{return 0;}}
		int size() const {return fsize;}

		Factmn& operator = (Factmn& a)
		{
			if(a.size() != fsize) {
				cout << "\033[1;31mERROR: Setting vectors of differnt dimensions\033[0m\n"<<" "<<endl;
			}else{

				for(int i=0; i<fsize; i++){
					fmn.at(i) = a.getFmn(i);
				}
			}

			return *this;
		}

		Factmn& operator += (Factmn& a)
		{
			if(a.size() != fsize) {
				cout << "\033[1;31mERROR: Setting vectors of differnt dimensions\033[0m\n"<<" "<<endl;
			}else{

				for(int i=0; i<fsize; i++){
					fmn.at(i) += a.getFmn(i);
				}
			}
			return *this;
		}


		Factmn& operator /= (double& c)
		{
			if(c!=0){
				for(int i=0; i<fsize; i++){
					fmn.at(i) /= c;
				}
			}else{
				cout << "\033[1;31mERROR: division by zero, moment unchanged\033[0m\n"<<" "<<endl;
			}
			return *this;
		}

};

class SMoment
{
	protected :
		int cent;
		vector<long double> nweight;
		vector<long long> nevent;
		vector<Factmn> Fmn;

	private:
		int mdim;
		int ndim;

		long double calcFmn(int mval_, int nval_, int CNT_){
			if(mval_>mdim || nval_>ndim || CNT_>cent) {
				cout << "\033[1;31mERROR: Outside dimension\033[0m\n"<<" "<<endl;
				if(CNT_>cent)cout<<"Need to increase no of centrality bins ("<<CNT_+1<<" > "<<cent<<")"<<endl;
				if(mval_>mdim)cout<<"Need to increase no of mdim ("<<mval_<<" > "<<mdim<<")"<<endl;
				if(nval_>ndim)cout<<"Need to increase no of ndim ("<<nval_<<" > "<<ndim<<")"<<endl;
				return 0;
			}

			long double TEMP_F = Fmn.at(CNT_).getFmn(mval_,nval_);
			TEMP_F /= (nweight.at(CNT_));

			//cout<<"Total Evnt# "<<nweight.at(CNT_)<<" F("<<mval_<<","<<nval_<<","<<CNT_<<")= "<<TEMP_F<<endl;
			return TEMP_F;
		}

	public:

		SMoment()
		{

		}

		SMoment(const int cent_, const int mdim_, const int ndim_)
		{
			cent = cent_; mdim = mdim_; ndim = ndim_;
			for(int i=0; i<cent; i++) 
			{
				Factmn TEMP_Fmn(mdim,ndim);
				Fmn.push_back(TEMP_Fmn);
				nweight.push_back(0.0L);
				nevent.push_back(0L);
			}
		}			

		virtual ~SMoment()
		{

			for(int i=0; i<cent; i++) 
				vector<Factmn>().swap(Fmn);
			for(int i=0; i<cent; i++) 
				vector<long double>().swap(nweight);
			for(int i=0; i<cent; i++) 
				vector<long long>().swap(nevent);
		}

		virtual long long getNevent(int CNT_) const {return nevent.at(CNT_);}

		virtual long double getNweight(int CNT_) const {return nweight.at(CNT_);}

		virtual long double calc_fact(long double Mval_, int mval_, long double Nval_, int nval_)
		{
			long double TEMP_FACT=1.0L;
			for(int i=0; i<mval_; i++) TEMP_FACT *= (Mval_-i);
			for(int j=0; j<nval_; j++) TEMP_FACT *= (Nval_-j);
			return TEMP_FACT;
		}

		void Fill(long double Mval_, long double Nval_, int CNT_, double weight_=1.)
		{
			if(CNT_<cent){
				for(int j=0; j<mdim; j++)
				{
					for(int k=0; k<ndim; k++)
					{
						//cout<<"Going to add "<<calc_fact(Mval_,j,Nval_,k)<<" at "<<j<<" "<<k<<endl;
						Fmn.at(CNT_).addFmn(j,k,calc_fact(Mval_,j,Nval_,k));

						//cout<<"Newvalue "<<Fmn.at(CNT_).getFmn(j,k)<<" at "<<j<<" "<<k<<endl;
					}
				}
				nweight.at(CNT_) += weight_;
				nevent.at(CNT_)++;
			}else{

				cout << "\033[1;31mERROR: Outside dimension\033[0m\n"<<" "<<endl;
				cout<<"Need to increase no of centrality bins ("<<CNT_+1<<" > "<<cent<<")"<<endl;
			}
		}


		long double getFmn(int mval_, int nval_, int CNT_) {return calcFmn(mval_,nval_,CNT_);}

};



class FactmnN : public Factmn{

	private :
		int nsize;
		int isize;
		vector<FactVec > fvv; 
		vector<FactMomN > fact;
		vector<int> abinN;

	public :


		FactmnN (const vector<FactVec > fvv_, const vector<int> abinN_) : Factmn (0,0){

			if(fvv_.size()<1){

				cout << "\033[1;31mERROR: The FactVector is empty CalcFactVec \033[0m\n"<<" "<<endl; abort();}

			fvv = fvv_;
			nsize = fvv.at(0).size();

			if(abinN_.size() != static_cast<size_t>(nsize)){
				cout << "\033[1;31mERROR: Mismatch in input dimension of moment orders (≠ abinN.size())\033[0m\n"<<" "<<endl;abort();}

			for(size_t l_=0; l_<fvv.size(); l_++) fmn.push_back(0.0L); //for IncMom
			for(size_t l_=0; l_<fvv.size(); l_++) fmn.push_back(0.0L); //for ObsMom
			
			fsize= fvv.size();
			
			fact.resize(fsize);

			abinN = abinN_;

			isize=0;
			for(size_t k_=0; k_<abinN.size(); k_++) {isize += abinN.at(k_);}


			cout<<"nsize= "<<nsize<<" abinsize= "<<abinN.size()<<" fsize= "<<fsize<<" fmnsize= "<<fmn.size()<<" isize= "<<isize<<endl;

			for(size_t l_=0; l_<fvv.size(); l_++)
			{
			//	fmn.push_back(0.0L);
				vector<int> mthN_;
				for(int i__=0; i__<fvv.at(l_).size(); i__++) mthN_.push_back(fvv.at(l_).getD(i__));
				fact.at(l_).SetFmnabN(mthN_,abinN);
				//cout<<" (";
				//fact.at(l_).printD();
				//cout<<") ";
			}
			//cout<<endl;
		}

		~FactmnN()
		{
			vector<FactVec >().swap(fvv); 
			vector<FactMomN >().swap(fact);
			vector<int >().swap(abinN);
		}

		void printD() {cout<<"< "; for(size_t l_=0; l_<fvv.size(); l_++){ fvv.at(l_).printD(); cout<<" ";} cout<<">"<<endl;}


		FactVec getfvv(int i__) const {if(static_cast<size_t>(i__)<fvv.size()){ return fvv.at(i__);}else{ cout<<"ERROR: outside dimensions"<<endl; return 0;}}


		void addFmnN(const vector<int> Np_, const vector<double> eff_, const double weight_)
		{
			for(int i__=0; i__<fsize; i__++)
			{
				long double IncMom = 0.0L;
				long long ObsMom = 0LL;
				fact.at(i__).GetMomentN(Np_,eff_,ObsMom,IncMom);  
				fmn.at(i__) += (ObsMom * weight_);
				fmn.at(fsize+i__) += (IncMom * weight_);
				//cout<<"} "<<"Obs= "<<ObsMom<<" Inc= "<<IncMom<<endl;
			}
		}


		long double getFmnN(int i__) const {if(i__<2*fsize){return fmn.at(i__);}else{return 0;}}

		size_t sizeI() const {return isize;}
		int sizeN() const {return nsize;}

};

class SMomentN : public SMoment{

	private :
		int norder;
		vector<FactmnN > FmnN;

	public :
		SMomentN (const int cent_, vector<FactVec> fmn_, const vector<int> abinN_) : SMoment()
		{
			cent = cent_;
			//FmnN.resize(cent_); 

			for(int i=0; i<cent; i++) 
			{
				FactmnN TEMP_FmnN(fmn_, abinN_); 
				FmnN.push_back(TEMP_FmnN);
				nweight.push_back(0.0L);
				nevent.push_back(0L);
			}

			for(size_t l_=0; l_<FmnN.size(); l_++) {FmnN.at(l_).printD(); cout<<endl;}
			cout<<"FmnN.size()= "<<FmnN.size()<<endl;
		}

		~SMomentN ()
		{
			vector<FactmnN>().swap(FmnN);
		}
		
		int centN() const {return cent;}

		void Fill(vector<int > Mval_, vector<double > eff_, int CNT_, double weight_=1.)
		{
			if(CNT_<cent){
				if(Mval_.size()!=FmnN.at(CNT_).sizeI()){
					cout << "\033[1;31mERROR: Adding multiplicity of different dimension\033[0m\n"<<" "<<endl;
					cout<<"Fmn= "<<FmnN.at(CNT_).sizeI()<<" Mval= "<<Mval_.size()<<endl;
				}else if(eff_.size()!=FmnN.at(CNT_).sizeI()){
					cout << "\033[1;31mERROR: Adding efficiency of different dimension\033[0m\n"<<" "<<endl;
				}else{
					FmnN.at(CNT_).addFmnN(Mval_,eff_,weight_);
					nweight.at(CNT_) += weight_;
					nevent.at(CNT_)++;
				}
			}else{
				//cout<<"ERROR: Outside dimension"<<endl;
				cout << "\033[1;31mERROR: Outside dimension\033[0m\n"<<" "<<endl;
				cout<<"Need to increase no of centrality bins ("<<CNT_+1<<" > "<<cent<<")"<<endl;
			}
		}

		void calcFmnN(int mval__, int CNT_){
			if(CNT_>cent) {
				//cout<<"ERROR: Outside dimension"<<endl; //return 0;
				cout << "\033[1;31mERROR: Outside dimension\033[0m\n"<<" "<<endl;
				cout<<"Need to increase no of centrality bins ("<<CNT_+1<<" > "<<cent<<")"<<endl;
			}

			if(mval__> 2*FmnN.at(CNT_).size()) {
				//cout<<"ERROR: Outside dimension"<<endl; //return 0;
				cout << "\033[1;31mERROR: Outside dimension\033[0m\n"<<" "<<endl;
				cout<<"mval__> 2*FmnN.at(CNT_).size()"<<endl;
			}

			long double TEMP_F = FmnN.at(CNT_).getFmnN(mval__);
			TEMP_F /= (nweight.at(CNT_));

			//cout<<"Total Weight/Evnt# "<<nweight.at(CNT_)<<" F(";//<<mval__<<","<<CNT_<<")= "<<TEMP_F<<endl;
			cout<<"Total Evnt# "<<nevent.at(CNT_)<<" F(";//<<mval__<<","<<CNT_<<")= "<<TEMP_F<<endl;


			if(mval__>=FmnN.at(CNT_).size()) FmnN.at(CNT_).getfvv(mval__-FmnN.at(CNT_).size()).printD();
			if(mval__<FmnN.at(CNT_).size()) FmnN.at(CNT_).getfvv(mval__).printD();

			cout<<","<<CNT_<<")= "<<TEMP_F<<endl;

		//	return TEMP_F;
		}


		void getFmnN( FactVec fvv__, int CNT_, long double& ObsMom_, long double& IncMom_)
		{
			int TEMP_pos =0;
			bool CHECK=0;
			
		//	cout<<endl;
		//	cout<<" We are looking for ";
		//	fvv__.printD();
		//	cout<<endl;

			for(int i__=0; i__<FmnN.at(CNT_).size(); i__++)
			{ 
				FactVec TEMP_fvv =FmnN.at(CNT_).getfvv(i__);
				if(TEMP_fvv==fvv__)
				{	
					TEMP_pos =i__;
					//cout<<" i= "<<i__<<" ";
					//TEMP_fvv.printD();
					//cout<<" = ";
					//fvv__.printD();
					//cout<<endl;
					CHECK += 1;
				}
			}

			if(CHECK==0) abort();

			long double TEMP_Obs = FmnN.at(CNT_).getFmnN(TEMP_pos);
			TEMP_Obs /= (nweight.at(CNT_));
			ObsMom_ = TEMP_Obs; //FmnN.at(CNT_).getFmnN(TEMP_pos);


			long double TEMP_Inc = FmnN.at(CNT_).getFmnN(TEMP_pos +  FmnN.at(CNT_).size());
			TEMP_Inc /= (nweight.at(CNT_));
			IncMom_ = TEMP_Inc; //FmnN.at(CNT_).getFmnN(TEMP_pos);

			//IncMom_ = FmnN.at(CNT_).getFmnN(TEMP_pos + FmnN.at(CNT_).size());

//			cout<<"ObsMom_= "<<ObsMom_<<" IncMom_= "<<IncMom_<<endl;
		}


		int size(const int CNT_=0) const {return FmnN.at(CNT_).size();}


//				void printD(string FORMAT="")
//				{
		void CalcMoment(const int Centrality_, const vector<int> norder_, long double& IncMom_, long double& ObsMom_,const string VERBOSE_="")
		{
			long double _IncMom__ = 0.0L;
			long double _ObsMom__ = 0.0L;

			MomVecN mnt;
			mnt.Mom2FactN(norder_);

			//if(VERBOSE_=="LOUD") mnt.printD();

			if(VERBOSE_=="LOUD"){
				cout<<"µ_";
				for(size_t __i=0; __i<norder_.size(); __i++) cout<<norder_.at(__i);
				cout<<"=";
			}

			for(size_t k_=0;k_<mnt.getvvmn().size();k_++)
			{

			if(VERBOSE_=="LOUD"){
				if(mnt.getvvmn().at(k_)==1){ if(k_!=0) cout<<" + ";
				}else if(mnt.getvvmn().at(k_)==-1){ cout<<" - ";
				}else if(mnt.getvvmn().at(k_)<0){
					cout<<" - "<<fabs(mnt.getvvmn().at(k_));
				}else{
					cout<<" + "<<fabs(mnt.getvvmn().at(k_));
				}
			}
				long double IncMom__ = 1.0L;
				long double ObsMom__ = 1.0L;
				FactVec temp_ = mnt.getCCmn().at(k_);


				//if(FORMAT!="NOEXPS"){
			if(VERBOSE_=="LOUD"){
				cout<<" f";
				mnt.getCCmn().at(k_).printD();
			}
				//}else if(FORMAT=="latex"){
				//	cout<<"{$rm f}_{";
				//	CCmn.at(k_).printD();
				//	cout<<"}";
				//}else{
				//	cout<<"f";
				//	CCmn.at(k_).printD();
				//}

				this->getFmnN(temp_, Centrality_, ObsMom__, IncMom__);
				_IncMom__ += (IncMom__*mnt.getvvmn().at(k_));	
				_ObsMom__ += (ObsMom__*mnt.getvvmn().at(k_));
			}

			ObsMom_ = _ObsMom__; IncMom_ = _IncMom__;

			if(VERBOSE_=="LOUD") cout<<" =("<<ObsMom_<<","<<IncMom_<<")"<<endl;

			//if(FORMAT!="NOEXPS")cout<<endl;
		}

		void CalcCentMom(const int Centrality_, const vector<int> norder_, long double& IncMom_, long double& ObsMom_,const string VERBOSE_="")
		{
			long double _IncMom__ = 0.0L;
			long double _ObsMom__ = 0LL;

			CentVecN * cnt =new CentVecN();
			cnt->Cent2FactN(norder_);

			if(VERBOSE_=="LOUD") cnt->printD();

			for(int i=0; i<cnt->size(); i++)
			{
				long double IncMom__ = 1.0L;
				long double ObsMom__ = 1LL;

				bool TEMP=1;

				for(size_t j=0; j<cnt->norder(); j++) 
				{

					int ADDPOW=0;
					FactVec temp(cnt->norder());
					temp.setD(j,1);

					FactVec temp2(cnt->norder());
					for(size_t k=0; k<cnt->norder(); k++) temp2.setD(k,cnt->getcfmn(i).getD(k));

					ADDPOW += cnt->getfpow(j,i);

					if(temp == temp2) {ADDPOW += 1; TEMP=0;}

					if(ADDPOW!=0){
						this->getFmnN(temp, Centrality_, ObsMom_, IncMom_);
						IncMom__ *= pow(IncMom_,ADDPOW);	
						ObsMom__ *= pow(ObsMom_,ADDPOW);
					}else{
						///cout<<endl;
						///cout<<" Skipping the calculation of f_";
						///temp.printD();
						///cout<<endl;
					}

					//if(j<cnt->norder()-1)	cout<<" {"<<ObsMom__<<","<<IncMom__<<"} ";//<<endl;
				}


				if(i<cnt->size()-1 && TEMP==1)
				{
					this->getFmnN(cnt->getcfmn(i), Centrality_, ObsMom_, IncMom_);
					//	cout<<" ["<<ObsMom_<<","<<IncMom_<<"] ";
					ObsMom__ *= ObsMom_; 
					IncMom__ *= IncMom_;
				} 

				IncMom__ *= (cnt->getcmn(i));
				ObsMom__ *= (cnt->getcmn(i));

				_IncMom__ += IncMom__;
				_ObsMom__ += ObsMom__;
			}


			ObsMom_ = _ObsMom__; IncMom_ = _IncMom__;

			if(VERBOSE_=="LOUD") cout<<" =("<<ObsMom_<<","<<IncMom_<<")"<<endl;

		}

		void CalcNudynMom(const int Centrality_, const int mth_, const NudynVecN * ndyn__, long double& IncMom__, long double& ObsMom__)
		{
			long double _IncMom__ = 0.0L;
			long double _ObsMom__ = 0LL;

			FactVec F10(2);
			F10.setD(0,1);
			FactVec F01(2);
			F01.setD(1,1);


			long double IncF10_ = 0.0L;
			long double ObsF10_ = 0LL;

			long double IncF01_ = 0.0L;
			long double ObsF01_ = 0LL;


			this->getFmnN(F10, Centrality_, ObsF10_, IncF10_);
			this->getFmnN(F01, Centrality_, ObsF01_, IncF01_);

			for(size_t l__=0; l__<ndyn__->getCCmn().size(); l__++){ 

				long double IncMom_ = 0.0L;
				long double ObsMom_ = 0LL;

				FactVec temp2 = ndyn__->getCCmn().at(l__);

				this->getFmnN(temp2, Centrality_, ObsMom_, IncMom_);

				IncMom_ *= ndyn__->getvvmn().at(l__); 
				ObsMom_ *= ndyn__->getvvmn().at(l__); 


				if(ndyn__->getvvmn().at(l__)==1){ if(l__!=0) cout<<" + ";
				}else if(ndyn__->getvvmn().at(l__)==-1){ cout<<" - ";
				}else if(ndyn__->getvvmn().at(l__)<0){
					cout<<" - "<<fabs(ndyn__->getvvmn().at(l__));
				}else{
					cout<<" + "<<fabs(ndyn__->getvvmn().at(l__));
				}

				cout<<" F_"; 
				ndyn__->getCCmn().at(l__).printD();

				cout<<"/";
				cout<<"( ";

				if(ndyn__->getCCmn().at(l__).getD(0) !=0){
					cout<<"F_10^"<<ndyn__->getCCmn().at(l__).getD(0);
					IncMom_ /= pow(IncF10_,ndyn__->getCCmn().at(l__).getD(0));
					ObsMom_ /= pow(ObsF10_,ndyn__->getCCmn().at(l__).getD(0));
				}

				cout<<" ";
				if(ndyn__->getCCmn().at(l__).getD(1) !=0)
				{
					cout<<"F_01^"<<ndyn__->getCCmn().at(l__).getD(1);
					IncMom_ /= pow(IncF01_,ndyn__->getCCmn().at(l__).getD(1));
					ObsMom_ /= pow(ObsF01_,ndyn__->getCCmn().at(l__).getD(1));
				}

				cout<<")";

				cout<<" == {"<<ObsMom_<<" "<<IncMom_<<"}"<<endl;

				_IncMom__ += IncMom_ ;
				_ObsMom__ += ObsMom_ ;

			}

			cout<<" Nudyn_"<<mth_;
			cout<<"("<<Centrality_<<")";
			cout<<" =["<<_IncMom__<<" "<<_ObsMom__<<"] ";

			cout<<endl;

			IncMom__ = _IncMom__; ObsMom__ = _ObsMom__;
		}


		void CalcCorr(const int Centrality_, const vector<int> mth_, CentVecN * cnt__, long double& IncMom__, long double& ObsMom__, const string VERBOSE_="")
		{
			for(size_t i__=0;i__<mth_.size();i__++){ 
				for(int j_=0; j_<=mth_.at(i__);j_++)
				{
					vector<int> mth__;
					mth__.push_back(j_);
					//cout<<"i__= "<<i__<<" "<<j_<<endl;
					cnt__->Corr2FactN(mth__,"NOEXPS","");  
					this->CalcCorrMom(Centrality_,mth__,cnt__,IncMom__,ObsMom__,VERBOSE_);
				}
			}
		}

		void CalcCorrMom(const int Centrality_, const vector<int> mth_, const CentVecN * cnt__, long double& IncMom__, long double& ObsMom__, const string VERBOSE_="")
		{
			long double _IncMom__ = 0.0L;
			long double _ObsMom__ = 0.0L;

			for(size_t l__=0; l__<cnt__->getCCmn().size(); l__++){ 

				long double IncMom_ = 0.0L;
				long double ObsMom_ = 0.0L;

				vector<int> norder__;
				for(int i__=0; i__<cnt__->getCCmn().at(l__).size(); i__++){
					norder__.push_back(cnt__->getCCmn().at(l__).getD(i__));
				}

//
//							cout<<" The central moment is = C_";
//							cnt__->getCCmn().at(l__).printD();
//							cout<<endl;

				this->CalcCentMom(Centrality_, norder__, IncMom_, ObsMom_,VERBOSE_);

				IncMom_ *= cnt__->getvvmn().at(l__); 
				ObsMom_ *= cnt__->getvvmn().at(l__); 

//				if(VERBOSE_=="LOUD") cout<<" == {"<<ObsMom_<<" "<<IncMom_<<"}"<<endl;

				_IncMom__ += IncMom_ ;
				_ObsMom__ += ObsMom_ ;

			}

			cout<<endl;
			cout<<endl;
			cout<<endl;
			
			if(VERBOSE_=="LOUD") {
				cout<<" Ç_";
				for(size_t j_=0; j_<mth_.size(); j_++) cout<<mth_.at(j_);	
				cout<<"("<<Centrality_<<")";
				cout<<" =["<<_IncMom__<<" "<<_ObsMom__<<"] ";

				cout<<endl;

				cout<<" √Ç_";
				for(size_t j_=0; j_<mth_.size(); j_++) cout<<mth_.at(j_);	
				cout<<"("<<Centrality_<<")";
				cout<<" =["<<sqrt(_IncMom__)<<" "<<sqrt(_ObsMom__)<<"] ";

				cout<<endl;
	

			}

			IncMom__ = _IncMom__; ObsMom__ = _ObsMom__;
		}

			//void CalcCorrMom(const int Centrality_, const vector<int> mth_, const CentVecN * cnt__, long double& IncMom__, long double& ObsMom__, const string VERBOSE_="")
			//if(FORMAT!="NOEXPS")cout<<"Func to convert cumulant K_m to correlation variables Ç_m = (∆N-<∆N>)^m or central moments w.r.to ∆N"<<endl; 
			//if(cnt__->calcheck()!=1) cnt__->Cumulant2Moment(mdim__,"NOEXPS");
			//

		

		void CalcCorrError(const int Centrality_, const vector<int> mth_, CentErrorN * cerr__, long double& IncMom__, long double& ObsMom__, const string VERBOSE_="")
		{
			cout<<endl;
			cout<<endl;

			long double _IncMom__ = 0.0L;
			long double _ObsMom__ = 0.0L;


			if(VERBOSE_=="CHECK") cout<<"mth_.size()= "<<mth_.size()<<" "<<cerr__->getnthE().size()<<endl;
		//	for(size_t i=0; i<cerr__->getnthE().size(); i++) cout<<cerr__->getnthE().at(i).size()<<endl;

			cout<<"∑="; cerr__->printE();
			cout<<endl;

                        for(size_t j__=0; j__<cerr__->getnthE().size(); j__++)
                        {
                                bool check=1;

                                for(size_t k_=0; k_<cerr__->getnthE().at(j__).size(); k_++)
                                {
                                        if(cerr__->getnthE().at(j__).at(k_).sumD() <=1) check*=0;
                                }

                                if(check){
                                        //if(cerr__->getcoeff().at(j__)==1){ if(j__!=0) cout<<" + ";
                                        //}else if(cerr__->getcoeff().at(j__)==-1){ cout<<" - ";
                                        //}else if(cerr__->getcoeff().at(j__)<0){
                                        //        cout<<" - "<<fabs(cerr__->getcoeff().at(j__));
                                        //}else{
                                        //        cout<<" + "<<fabs(cerr__->getcoeff().at(j__));
                                        //}

					long double __IncMom_ = 1.0L;
					long double __ObsMom_ = 1.0L;
                                        for(size_t k_=0; k_<cerr__->getnthE().at(j__).size(); k_++)
                                        {
                                                //cout<<" Ç_";
                                                //cerr__->getnthE().at(j__).at(k_).printD();
						vector<int> mthE_;
						for(int l__=0; l__<cerr__->getnthE().at(j__).at(k_).size(); l__++)
							mthE_.push_back(cerr__->getnthE().at(j__).at(k_).getD(l__));

						CentVecN * cnt__ = new CentVecN();
						cnt__->Corr2CentN(mthE_,"NOEXPS");
						//cnt__->Corr2CentN(mthE_);

						long double IncMom_ = 0.0L;
						long double ObsMom_ = 0.0L;

						this->CalcCorrMom(Centrality_, mthE_, cnt__, IncMom_, ObsMom_,"");

						//cout<<"__IncMom_="<<__IncMom_<<" IncMom_="<<IncMom_<<" _IncMom__="<<_IncMom__<<endl;

						__IncMom_ *= pow(IncMom_,cerr__->getnthE().at(j__).at(k_).getE());
						__ObsMom_ *= pow(ObsMom_,cerr__->getnthE().at(j__).at(k_).getE());

                                                //if(cerr__->getnthE().at(j__).at(k_).expo()>1) cout<<"^"<<cerr__->getnthE().at(j__).at(k_).getE();
                                                //cout<<" ";
                                        }

					_IncMom__ += (cerr__->getcoeff().at(j__) * __IncMom_);
					_ObsMom__ += (cerr__->getcoeff().at(j__) * __ObsMom_);
                                }
                        }


/*
			for(size_t j__=0; j__<cerr__->getnthE().size(); j__++)
			{
				bool check=1;

				for(size_t k_=0; k_<cerr__->getnthE().at(j__).size(); k_++)
				{
					if(cerr__->getnthE().at(j__).at(k_).sumD() <=1) check*=0;
				}

				if(check){

					//if(cerr__->getcoeff().at(j__)==1){ if(j__!=0) cout<<" + ";
					//}else if(cerr__->getcoeff().at(j__)==-1){ cout<<" - ";
					//}else if(cerr__->getcoeff().at(j__)<0){
					//	cout<<" - "<<fabs(cerr__->getcoeff().at(j__));
					//}else{
					//	cout<<" + "<<fabs(cerr__->getcoeff().at(j__));
					//}

					//cout<<"[";

					long double __IncMom_ = 1.0L;
					long double __ObsMom_ = 1.0L;

					for(size_t k_=0; k_<cerr__->getnthE().at(j__).size(); k_++)
					{
						vector<int> mthE_;
						for(int l__=0; l__<cerr__->getnthE().at(j__).at(k_).size(); l__++)
							mthE_.push_back(cerr__->getnthE().at(j__).at(k_).getD(l__));

						CentVecN * cnt__ = new CentVecN();
						//cnt__->Corr2CentN(mthE_,"NOEXPS");
						cnt__->Corr2CentN(mthE_);

						long double IncMom_ = 0.0L;
						long double ObsMom_ = 0.0L;

						this->CalcCorrMom(Centrality_, mthE_, cnt__, IncMom_, ObsMom_,"");

						cout<<"__IncMom_="<<__IncMom_<<" IncMom_="<<IncMom_<<" _IncMom__="<<_IncMom__<<endl;

						__IncMom_ *= IncMom_;
						__ObsMom_ *= ObsMom_;


//						cout<<"{";
//						cnt__.Corr2FactN(mthE_,"NOEXPS");  
//						cout<<"}";
					}
					//cout<<"]";

					_IncMom__ += (cerr__->getcoeff().at(j__) * __IncMom_);
					_ObsMom__ += (cerr__->getcoeff().at(j__) * __ObsMom_);


				}
			}

			cout<<endl;
			cout<<endl;
*/

			vector<int> mth1__;
			vector<int> mth2__;

			//if(VERBOSE_=="CHECK") cout<<"cerr__->size()= "<<cerr__->size()<<endl;

			for(size_t i_=0; i_<mth_.size()/2; i_++)
			{
				//				cout<<"i_= "<<i_<<" mth_at(i_)= "<<mth_.at(i_)<<endl;
				mth1__.push_back(mth_.at(i_));
				mth2__.push_back(mth_.at(i_+mth_.size()/2));
			}


			if(VERBOSE_!="NOEXPS"){

				cout<<"∑_{";
				for(size_t l__=0; l__<mth1__.size(); l__++) cout<<mth1__.at(l__);
				//cerr__->getnthE().at(1).at(0).printD();
				cout<<",";
				for(size_t l__=0; l__<mth2__.size(); l__++) cout<<mth2__.at(l__);
				//cerr__->getnthE().at(1).at(1).printD();
				cout<<"}";
				cout<<" ("<<Centrality_<<")";
				cout<<" =["<<_IncMom__<<" "<<_ObsMom__<<"] ";

				cout<<endl;
				cout<<endl;

				cout<<"Error{";
				for(size_t l__=0; l__<mth1__.size(); l__++) cout<<mth1__.at(l__);
				//cerr__->getnthE().at(1).at(0).printD();
				cout<<",";
				for(size_t l__=0; l__<mth2__.size(); l__++) cout<<mth2__.at(l__);
				//cerr__->getnthE().at(1).at(1).printD();
				cout<<"}";
				cout<<" ("<<Centrality_<<")";
				cout<<" =["<<sqrt(_IncMom__/nevent.at(Centrality_))<<" "<<sqrt(_ObsMom__/nevent.at(Centrality_))<<"] ";


				cout<<endl;
				cout<<endl;

			}

			IncMom__ = _IncMom__;
			ObsMom__ = _ObsMom__;
		}

		void CalcCorrMom(const int Centrality_, const vector<int> mth_, CentErrorN * cerr__, long double& IncMom__, long double& ObsMom__, const string VERBOSE_="", bool CrossTerm=0)
		{

			vector<int> mth1__;
			vector<int> mth2__;

			if(VERBOSE_=="CHECK") cout<<"cerr__->size()= "<<cerr__->size()<<endl;

			for(size_t i_=0; i_<mth_.size()/2; i_++)
			{
				//				cout<<"i_= "<<i_<<" mth_at(i_)= "<<mth_.at(i_)<<endl;
				mth1__.push_back(mth_.at(i_));
				mth2__.push_back(mth_.at(i_+mth_.size()/2));
			}

			bool CHECK =1;
			for(size_t i_=0; i_<mth1__.size(); i_++)
			{
				//				cout<<mth1__.at(i_)<<" "<<mth2__.at(i_)<<endl;
				if(mth1__.at(i_)!=mth2__.at(i_)) 
				{
					cout<<"Warning : CentErrorN input is for a cross term"<<endl;
					cout<<"Warning : Multiple moments have to be calculated"<<endl;
					CHECK *=0;
				}
			}

			if(CHECK==1) {
				CentVecN * cnt__ = new CentVecN();
				cnt__->Corr2CentN(mth1__);

				this->CalcCorrMom(Centrality_, mth1__, cnt__, IncMom__, ObsMom__,"LOUD");
			}else{
				if(CrossTerm==0){
					CentVecN * cnt1__ = new CentVecN();
					cnt1__->Corr2CentN(mth1__);
					this->CalcCorrMom(Centrality_, mth1__, cnt1__, IncMom__, ObsMom__,"");
				}else{
					CentVecN * cnt2__ = new CentVecN();
					cnt2__->Corr2CentN(mth2__);
					this->CalcCorrMom(Centrality_, mth2__, cnt2__, IncMom__, ObsMom__,"");
				}
			}

		}

		void CalcCumulant(const int Centrality_, const int mdim__, CumulantVec * clt__, long double& IncMom__, long double& ObsMom__, const string VERBOSE_="")
		{

			long double _IncMom__ = 0.0L;
			long double _ObsMom__ = 0.0L;

			clt__->Cumulant2CorrMoment(mdim__);


			size_t i__=mdim__-1; //clt__->getnthC().size()-1; 

			for(size_t j__=0; j__<clt__->getnthC().at(i__).size(); j__++)
			{
				bool CHECK=1;

				for(size_t k__=0; k__<clt__->getnthC().at(i__).at(j__).size(); k__++)
				{
					if(clt__->getnthC().at(i__).at(j__).at(k__).getD(0)==1) CHECK *=0;
				}

				long double IncMom0_ = 0.0L;
				long double ObsMom0_ = 0.0L;

				if(CHECK!=0){

					long prefact = clt__->getcoeffC().at(i__).at(j__);

					if(VERBOSE_=="CHECK"){
						if(prefact==1){ if(j__!=0) cout<<" + ";
						}else if(prefact==-1){ cout<<" - ";
						}else if(prefact<0){
							cout<<" - "<<fabs(prefact);
						}else{
							cout<<" + "<<fabs(prefact);
						}
					}

					IncMom0_ = prefact;
					ObsMom0_ = prefact;
  
					for(size_t k__=0; k__<clt__->getnthC().at(i__).at(j__).size(); k__++)
					{
						long double IncMom_ = 0.0L;
						long double ObsMom_ = 0.0L;


						if(VERBOSE_=="CHECK"){
							cout<<"Ç_";
							clt__->getnthC().at(i__).at(j__).at(k__).printD();
							if(clt__->getnthC().at(i__).at(j__).at(k__).expo()>1) cout<<"^"<<clt__->getnthC().at(i__).at(j__).at(k__).expo();
						}


						vector<int> _mth;
						_mth.push_back(clt__->getnthC().at(i__).at(j__).at(k__).getD(0));

						CentVecN * ___cnt = new CentVecN();
						___cnt->Corr2CentN(_mth,"NOEXPS","");
						//__cnt->Corr2FactN(_mth,"NOEXPS","");
						this->CalcCorrMom(Centrality_, _mth, ___cnt, IncMom_, ObsMom_);//,"LOUD");

						//cout<<"(IncMom_="<<IncMom_<<" "<<"ObsMom_="<<ObsMom_<<") "<<endl;

						 IncMom0_ *= pow(IncMom_,clt__->getnthC().at(i__).at(j__).at(k__).expo());
						 ObsMom0_ *= pow(ObsMom_,clt__->getnthC().at(i__).at(j__).at(k__).expo());
					}

				}
				_IncMom__ += IncMom0_;
				_ObsMom__ += ObsMom0_;

			}



			if(mdim__==1){
				FactVec F10(2);
				F10.setD(0,1);
				FactVec F01(2);
				F01.setD(1,1);

				long double IncMom_ = 0.0L;
				long double ObsMom_ = 0.0L;

				this->getFmnN(F10, Centrality_, ObsMom_, IncMom_);
				_IncMom__ += IncMom_;
				_ObsMom__ += ObsMom_;

				IncMom_ = 0.0L;
				ObsMom_ = 0.0L;

				this->getFmnN(F01, Centrality_, ObsMom_, IncMom_);
				_IncMom__ -= IncMom_;
				_ObsMom__ -= ObsMom_;
			}


			if(VERBOSE_=="LOUD") {
				cout<<" K_"<<mdim__;
				cout<<"("<<Centrality_<<")";
				cout<<" =["<<_IncMom__<<" "<<_ObsMom__<<"] ";

				cout<<endl;

				//cout<<" √K_"<<mdim__;
				//cout<<"("<<Centrality_<<")";
				//cout<<" =["<<sqrt(_IncMom__)<<" "<<sqrt(_ObsMom__)<<"] ";

				//cout<<endl;
			}

			IncMom__ = _IncMom__;
			ObsMom__ = _ObsMom__;
		}



		void CalcCumulantError(const int Centrality_, const int mdim__, CumulantVec * clt__, long double& IncMom__, long double& ObsMom__, const string VERBOSE_="")
		{

			long double _IncMom__ = 0.0L;
			long double _ObsMom__ = 0.0L;

			if(mdim__<=1){ 
				if(VERBOSE_=="CHECK"){		
					cout<<endl;
					cout<<"Var(K_1)=Ç_2";
					cout<<endl;
				}		
				this->CalcCumulant(Centrality_,2,clt__,_IncMom__,_ObsMom__,"LOUD");
			}else{


				vector<vector<FactVec > > Kdelt; 
				vector<vector<long > > Cdelt;
				vector<FactVec > fvv; 
				CumulantVec _clt;
				_clt.CalcVariance(mdim__,Kdelt,Cdelt,fvv);

				size_t i__=mdim__-1;


				if(VERBOSE_=="CHECK")	cout<<endl;
				if(VERBOSE_=="CHECK")	cout<<"Var(K_"<<i__+1<<")=";
				for(size_t j_=0; j_<fvv.size(); j_++){
					for(size_t k_=0; k_<fvv.size(); k_++){

						if(fvv.at(k_).sumD()==0) continue;

						if(fvv.at(j_) == fvv.at(k_)) 
						{ 
							long double IncMom1_ = 0.0L;
							long double ObsMom1_ = 0.0L;
							if(j_!=0){	


								if(VERBOSE_=="CHECK")cout<<"(";
								for(size_t __i=0; __i<Kdelt.at(j_).size(); __i++)
								{
									long prefact = Cdelt.at(j_).at(__i);

									if(VERBOSE_=="CHECK"){
										if(prefact==1){ if(__i!=0) cout<<" + ";
										}else if(prefact==-1){ cout<<" - ";
										}else if(prefact<0){
											cout<<" - "<<fabs(prefact);
										}else{
											cout<<" + "<<fabs(prefact);
										}
									}

									for(int __k=0; __k<Kdelt.at(j_).at(__i).size()-1; __k += 2) 
									{
										if(VERBOSE_=="CHECK"){

											cout<<" Ç_";
											cout<<Kdelt.at(j_).at(__i).getD(__k);
										}

										long double IncMom0_ = 0.0L;
										long double ObsMom0_ = 0.0L;

										//this->CalcCumulant(Centrality_,Kdelt.at(j_).at(__i).getD(__k),clt__,IncMom0_,ObsMom0_,"LOUD");

										vector<int> __mdim;
										__mdim.push_back(Kdelt.at(j_).at(__i).getD(__k));
										CentVecN * cent_= new CentVecN();  //(Kdelt.at(j_).at(__i).getD(__k));
										cent_->Corr2CentN(__mdim);
										this->CalcCorrMom(Centrality_,__mdim,cent_,IncMom0_,ObsMom0_,"");//LOUD");

										//this->CalcCorrMom(Centrality_,Kdelt.at(j_).at(__i).getD(__k),cent_,IncMom0_,ObsMom0_,"LOUD");

										if(Kdelt.at(j_).at(__i).getD(__k+1)>1)
										{

											IncMom1_ += (prefact*pow(IncMom0_,Kdelt.at(j_).at(__i).getD(__k+1)));
											ObsMom1_ += (prefact*pow(ObsMom0_,Kdelt.at(j_).at(__i).getD(__k+1)));
											if(VERBOSE_=="CHECK"){
												cout<<"^";
												cout<<Kdelt.at(j_).at(__i).getD(__k+1);
											}		

										}
									}
								}


								if(VERBOSE_=="CHECK")			cout<<" )^2 ";
							}

							vector<int> _mth;

							_mth.push_back(fvv.at(j_).getD(0)); 
							_mth.push_back(fvv.at(k_).getD(0)); 

							// call central moment error estimation class
							CentErrorN cerr(_mth);

							//calculate the expression of error 
							cerr.CalcCentError("");
							if(VERBOSE_=="CHECK"){
								cout<<" (";
								cerr.printE();
								cout<<") ";
							}

							long double IncMom2_ = 0.0L;
							long double ObsMom2_ = 0.0L;

							vector<int> _mth_;
							//for(size_t l__=0;l__<_mth.size();l__++);

							this->CalcCorrError(Centrality_,_mth,&cerr,IncMom2_,ObsMom2_,"");

							if(IncMom1_!=0){
								_IncMom__ += IncMom1_*IncMom1_*IncMom2_;
							}else{
								_IncMom__ += IncMom2_;
							}	
							if(ObsMom1_!=0){
								_ObsMom__ += ObsMom1_*ObsMom1_*ObsMom2_;
							}else{
								_ObsMom__ += ObsMom2_;
							}

							//cout<<" ∑(";
							//fvv.at(j_).printD();
							//cout<<",";
							//fvv.at(k_).printD();
							//cout<<")";


						}else{
							if(fvv.at(j_).sumD()>fvv.at(k_).sumD()){

								if(VERBOSE_=="CHECK")		cout<<"2 ";

								long double IncMom1_ = 0.0L;
								long double ObsMom1_ = 0.0L;
								if(j_!=0){
									if(VERBOSE_=="CHECK")	cout<<"(";
									for(size_t __i=0; __i<Kdelt.at(j_).size(); __i++)
									{
										long prefact = Cdelt.at(j_).at(__i);

										if(VERBOSE_=="CHECK"){
											if(prefact==1){ if(__i!=0) cout<<" + ";
											}else if(prefact==-1){ cout<<" - ";
											}else if(prefact<0){
												cout<<" - "<<fabs(prefact);
											}else{
												cout<<" + "<<fabs(prefact);
											}

										}

										for(int __k=0; __k<Kdelt.at(j_).at(__i).size()-1; __k += 2)
										{
											if(VERBOSE_=="CHECK"){
												cout<<" Ç_";
												cout<<Kdelt.at(j_).at(__i).getD(__k);
											}

											long double IncMom0_ = 0.0L;
											long double ObsMom0_ = 0.0L;

											//this->CalcCumulant(Centrality_,Kdelt.at(j_).at(__i).getD(__k),clt__,IncMom0_,ObsMom0_,"LOUD");
											//CentVecN cent_(Kdelt.at(j_).at(__i).getD(__k));
											//this->CalcCorrMom(Centrality_,Kdelt.at(j_).at(__i).getD(__k),cent_,IncMom0_,ObsMom0_,"LOUD");

											vector<int> __mdim;
											__mdim.push_back(Kdelt.at(j_).at(__i).getD(__k));
											CentVecN * cent__= new CentVecN();  //(Kdelt.at(j_).at(__i).getD(__k));
											cent__->Corr2CentN(__mdim);
											this->CalcCorrMom(Centrality_,__mdim,cent__,IncMom0_,ObsMom0_,"");//LOUD");






											if(Kdelt.at(j_).at(__i).getD(__k+1)>1)
											{
												IncMom1_ += (prefact*pow(IncMom0_,Kdelt.at(j_).at(__i).getD(__k+1)));
												ObsMom1_ += (prefact*pow(ObsMom0_,Kdelt.at(j_).at(__i).getD(__k+1)));
												if(VERBOSE_=="CHECK"){	
													cout<<"^";
													cout<<Kdelt.at(j_).at(__i).getD(__k+1);
												}
											}
										}
									}
									if(VERBOSE_=="CHECK")	cout<<" )";
								}

								long double IncMom3_ = 0.0L;
								long double ObsMom3_ = 0.0L;
								if(VERBOSE_=="CHECK")	cout<<"(";
								for(size_t __i=0; __i<Kdelt.at(k_).size(); __i++)
								{
									long prefact = Cdelt.at(k_).at(__i);

									if(VERBOSE_=="CHECK"){
										if(prefact==1){ if(__i!=0) cout<<" + ";
										}else if(prefact==-1){ cout<<" - ";
										}else if(prefact<0){
											cout<<" - "<<fabs(prefact);
										}else{
											cout<<" + "<<fabs(prefact);
										}
									}

									for(int __k=0; __k<Kdelt.at(k_).at(__i).size()-1; __k += 2)
									{
										if(VERBOSE_=="CHECK")	cout<<" Ç_";
										if(VERBOSE_=="CHECK")	cout<<Kdelt.at(k_).at(__i).getD(__k);

										long double IncMom0_ = 0.0L;
										long double ObsMom0_ = 0.0L;

										//this->CalcCumulant(Centrality_,Kdelt.at(k_).at(__i).getD(__k),clt__,IncMom0_,ObsMom0_,"LOUD");
										//CentVectN cent_(Kdelt.at(k_).at(__i).getD(__k));
										//this->CalcCorrMom(Centrality_,Kdelt.at(k_).at(__i).getD(__k),cent_,IncMom0_,ObsMom0_,"LOUD");

										vector<int> __mdim;
										__mdim.push_back(Kdelt.at(k_).at(__i).getD(__k));
										CentVecN * __cent__= new CentVecN();  //(Kdelt.at(j_).at(__i).getD(__k));
										__cent__->Corr2CentN(__mdim);
										this->CalcCorrMom(Centrality_,__mdim,__cent__,IncMom0_,ObsMom0_,"");//LOUD");





										if(Kdelt.at(k_).at(__i).getD(__k+1)>1)
										{
											IncMom3_ += (prefact*pow(IncMom0_,Kdelt.at(k_).at(__i).getD(__k+1)));
											ObsMom3_ += (prefact*pow(ObsMom0_,Kdelt.at(k_).at(__i).getD(__k+1)));
											if(VERBOSE_=="CHECK")		cout<<"^";
											if(VERBOSE_=="CHECK")		cout<<Kdelt.at(k_).at(__i).getD(__k+1);
										}
									}
								}
								if(VERBOSE_=="CHECK")	cout<<" )";


								vector<int> _mth;

								_mth.push_back(fvv.at(j_).getD(0));
								_mth.push_back(fvv.at(k_).getD(0));

								// call central moment error estimation class
								CentErrorN cerr(_mth);

								//calculate the expression of error 
								cerr.CalcCentError("");

								if(VERBOSE_=="CHECK"){
									cout<<" (";
									cerr.printE();
									cout<<") ";
								}

								long double IncMom2_ = 0.0L;
								long double ObsMom2_ = 0.0L;

								//vector<int> _mth_;
								//for(size_t l__=0;l__<_mth.size();l__++);

								this->CalcCorrError(Centrality_,_mth,&cerr,IncMom2_,ObsMom2_,"");

								//IncMom2_ can't be zero, but IncMom1_ & IncMom3_ which are the derivative terms can be zero
								if(IncMom1_!=0 && IncMom3_!=0){
									_IncMom__ +=  2.*IncMom3_*IncMom2_*IncMom1_;
								}else if(IncMom1_!=0 && IncMom3_==0){
									_IncMom__ +=  2.*IncMom1_*IncMom2_;
								}else if(IncMom1_==0 && IncMom3_!=0){
									_IncMom__ +=  2.*IncMom3_*IncMom2_;
								}else{
									_IncMom__ +=  2.*IncMom2_;
								}	

								if(ObsMom1_!=0 && ObsMom3_!=0){
									_ObsMom__ +=  2.*ObsMom3_*ObsMom2_*ObsMom1_;
								}else if(ObsMom1_!=0 && ObsMom3_==0){
									_ObsMom__ +=  2.*ObsMom1_*ObsMom2_;
								}else if(ObsMom1_==0 && ObsMom3_!=0){
									_ObsMom__ +=  2.*ObsMom3_*ObsMom2_;
								}else{
									_ObsMom__ +=  2.*ObsMom2_;
								}	




								//cout<<" ∑(";
								//fvv.at(j_).printD();
								//cout<<",";
								//fvv.at(k_).printD();
								//cout<<")";

							}else{
								if(VERBOSE_=="CHECK")	cout<<"2 ";

								long double IncMom1_ = 0.0L;
								long double ObsMom1_ = 0.0L;

								if(VERBOSE_=="CHECK")	cout<<"(";
								for(size_t __i=0; __i<Kdelt.at(k_).size(); __i++)
								{
									long prefact = Cdelt.at(k_).at(__i);

									if(VERBOSE_=="CHECK"){
										if(prefact==1){ if(__i!=0) cout<<" + ";
										}else if(prefact==-1){ cout<<" - ";
										}else if(prefact<0){
											cout<<" - "<<fabs(prefact);
										}else{
											cout<<" + "<<fabs(prefact);
										}

									}
									for(int __k=0; __k<Kdelt.at(k_).at(__i).size()-1; __k += 2)
									{
										if(VERBOSE_=="CHECK")	cout<<" Ç_";
										if(VERBOSE_=="CHECK")	cout<<Kdelt.at(k_).at(__i).getD(__k);

										long double IncMom0_ = 0.0L;
										long double ObsMom0_ = 0.0L;

										//									this->CalcCumulant(Centrality_,Kdelt.at(k_).at(__i).getD(__k),clt__,IncMom0_,ObsMom0_,"LOUD");
										//CentVectN cent_(Kdelt.at(k_).at(__i).getD(__k));
										//this->CalcCorrMom(Centrality_,Kdelt.at(k_).at(__i).getD(__k),cent_,IncMom0_,ObsMom0_,"LOUD");

										vector<int> __mdim;
										__mdim.push_back(Kdelt.at(k_).at(__i).getD(__k));
										CentVecN * _cent___= new CentVecN();  //(Kdelt.at(j_).at(__i).getD(__k));
										_cent___->Corr2CentN(__mdim);
										this->CalcCorrMom(Centrality_,__mdim,_cent___,IncMom0_,ObsMom0_,"");//LOUD");





										if(Kdelt.at(k_).at(__i).getD(__k+1)>1)
										{
											IncMom1_ += (prefact*pow(IncMom0_,Kdelt.at(k_).at(__i).getD(__k+1)));
											ObsMom1_ += (prefact*pow(ObsMom0_,Kdelt.at(k_).at(__i).getD(__k+1)));
											if(VERBOSE_=="CHECK")	cout<<"^";
											if(VERBOSE_=="CHECK")	cout<<Kdelt.at(k_).at(__i).getD(__k+1);
										}
									}
								}
								if(VERBOSE_=="CHECK")	cout<<" )";

								long double IncMom3_ = 0.0L;
								long double ObsMom3_ = 0.0L;
								if(j_!=0){
									if(VERBOSE_=="CHECK")		cout<<"(";
									for(size_t __i=0; __i<Kdelt.at(j_).size(); __i++)
									{
										long prefact = Cdelt.at(j_).at(__i);
										if(VERBOSE_=="CHECK"){
											if(prefact==1){ if(__i!=0) cout<<" + ";
											}else if(prefact==-1){ cout<<" - ";
											}else if(prefact<0){
												cout<<" - "<<fabs(prefact);
											}else{
												cout<<" + "<<fabs(prefact);
											}

										}
										for(int __k=0; __k<Kdelt.at(j_).at(__i).size()-1; __k += 2)
										{
											if(VERBOSE_=="CHECK")	cout<<" Ç_";
											if(VERBOSE_=="CHECK")	cout<<Kdelt.at(j_).at(__i).getD(__k);

											long double IncMom0_ = 0.0L;
											long double ObsMom0_ = 0.0L;

											//this->CalcCumulant(Centrality_,Kdelt.at(j_).at(__i).getD(__k),clt__,IncMom0_,ObsMom0_,"LOUD");
											//CentVectN cent_(Kdelt.at(j_).at(__i).getD(__k));
											//this->CalcCorrMom(Centrality_,Kdelt.at(j_).at(__i).getD(__k),cent_,IncMom0_,ObsMom0_,"LOUD");

											vector<int> __mdim;
											__mdim.push_back(Kdelt.at(j_).at(__i).getD(__k));
											CentVecN * _cent__= new CentVecN();  //(Kdelt.at(j_).at(__i).getD(__k));
											_cent__->Corr2CentN(__mdim);
											this->CalcCorrMom(Centrality_,__mdim,_cent__,IncMom0_,ObsMom0_,"");//LOUD");





											if(Kdelt.at(j_).at(__i).getD(__k+1)>1)
											{
												IncMom3_ += (prefact*pow(IncMom0_,Kdelt.at(j_).at(__i).getD(__k+1)));
												ObsMom3_ += (prefact*pow(ObsMom0_,Kdelt.at(j_).at(__i).getD(__k+1)));
												if(VERBOSE_=="CHECK")	cout<<"^";
												if(VERBOSE_=="CHECK")	cout<<Kdelt.at(j_).at(__i).getD(__k+1);
											}
										}
									}
									if(VERBOSE_=="CHECK")	cout<<" )";
								}

								vector<int> _mth;

								_mth.push_back(fvv.at(k_).getD(0));
								_mth.push_back(fvv.at(j_).getD(0));

								// call central moment error estimation class
								CentErrorN cerr(_mth);

								//calculate the expression of error 
								cerr.CalcCentError("");
								if(VERBOSE_=="CHECK"){
									cout<<" (";
									cerr.printE();
									cout<<") ";

								}
								long double IncMom2_ = 0.0L;
								long double ObsMom2_ = 0.0L;

								//vector<int> _mth_;
								//for(size_t l__=0;l__<_mth.size();l__++);

								this->CalcCorrError(Centrality_,_mth,&cerr,IncMom2_,ObsMom2_,"");

								//IncMom2_ can't be zero, but IncMom1_ & IncMom3_ which are the derivative terms can be zero
								if(IncMom1_!=0 && IncMom3_!=0){
									_IncMom__ +=  2.*IncMom3_*IncMom2_*IncMom1_;
								}else if(IncMom1_!=0 && IncMom3_==0){
									_IncMom__ +=  2.*IncMom1_*IncMom2_;
								}else if(IncMom1_==0 && IncMom3_!=0){
									_IncMom__ +=  2.*IncMom3_*IncMom2_;
								}else{
									_IncMom__ +=  2.*IncMom2_;
								}	

								if(ObsMom1_!=0 && ObsMom3_!=0){
									_ObsMom__ +=  2.*ObsMom3_*ObsMom2_*ObsMom1_;
								}else if(ObsMom1_!=0 && ObsMom3_==0){
									_ObsMom__ +=  2.*ObsMom1_*ObsMom2_;
								}else if(ObsMom1_==0 && ObsMom3_!=0){
									_ObsMom__ +=  2.*ObsMom3_*ObsMom2_;
								}else{
									_ObsMom__ +=  2.*ObsMom2_;
								}	



								//cout<<" ∑(";
								//fvv.at(k_).printD();
								//cout<<",";
								//fvv.at(j_).printD();
								//cout<<")";
							}
						}


						if(j_==fvv.size()-1 && k_==fvv.size()-1){
							if(VERBOSE_=="CHECK")	cout<<"";
						}else{
							if(VERBOSE_=="CHECK")	cout<<" + ";
						}
					}
					fvv.at(j_) *= 0;
				}
				if(VERBOSE_=="CHECK")		cout<<endl;

			}


			IncMom__ = (nevent.at(Centrality_)>0 && _IncMom__>0) ? sqrt(_IncMom__/nevent.at(Centrality_)) : 0;
			ObsMom__ = (nevent.at(Centrality_)>0 && _ObsMom__>0) ? sqrt(_ObsMom__/nevent.at(Centrality_)) : 0;

			if(VERBOSE_!="NOEXPS"){

				cout<<endl;
				cout<<endl;
				cout<<"∑_{";
				//for(size_t l__=0; l__<mth1__.size(); l__++) cout<<mth1__.at(l__);
				//cerr__->getnthE().at(1).at(0).printD();
				//cout<<",";
				//for(size_t l__=0; l__<mth2__.size(); l__++) cout<<mth2__.at(l__);
				//cerr__->getnthE().at(1).at(1).printD();
				cout<<"K_"<<mdim__;				
				cout<<"}";
				cout<<" ("<<Centrality_<<")";
				cout<<" =["<<_IncMom__<<" "<<_ObsMom__<<"] ";

				cout<<endl;
				cout<<endl;

				cout<<"Error{";
				//for(size_t l__=0; l__<mth1__.size(); l__++) cout<<mth1__.at(l__);
				////cerr__->getnthE().at(1).at(0).printD();
				//cout<<",";
				//for(size_t l__=0; l__<mth2__.size(); l__++) cout<<mth2__.at(l__);
				////cerr__->getnthE().at(1).at(1).printD();
				cout<<"K_"<<mdim__;				
				cout<<"}";
				cout<<" ("<<Centrality_<<")";
				//cout<<"}";
				//cout<<" ("<<Centrality_<<")";
				//cout<<" =["<<sqrt(_IncMom__/nevent.at(Centrality_))<<" "<<sqrt(_ObsMom__/nevent.at(Centrality_))<<"] ";
				cout<<" =["<<IncMom__<<" "<<ObsMom__<<"] ";

				cout<<endl;
				cout<<endl;

			}

		}


		void CalcSICumulant(const int Centrality_, const int mdim__, SICumulantVec * sic__, long double& IncMom__, long double& ObsMom__, const string VERBOSE_="")
		{

			long double _IncMom__ = 0.0L;
			long double _ObsMom__ = 0.0L;

			if(mdim__>0){
				size_t i__=mdim__-1; //sic__->getnthC().size()-1; 


				/////////////////////////////////////////////////
				//cout<<" (SI)K_"<<i__+1<<"= ";
				for(size_t j__=0; j__<sic__->getnthC().at(i__).size(); j__++)
				{


					long double IncMom0_ = 0.0L;
					long double ObsMom0_ = 0.0L;


					long prefact = sic__->getcoeffC().at(i__).at(j__);

					//if(prefact==1){ if(j__!=0) cout<<" + ";
					//}else if(prefact==-1){ cout<<" - ";
					//}else if(prefact<0){
					//	cout<<" - "<<fabs(prefact);
					//}else{
					//	cout<<" + "<<fabs(prefact);
					//}


					IncMom0_ = prefact;
					ObsMom0_ = prefact;


					for(size_t k__=0; k__<sic__->getnthC().at(i__).at(j__).size(); k__++)
					{
						if(sic__->getnthC().at(i__).at(j__).at(k__).expo()>0){

							vector<int> _mth;
							for(size_t __i=0;__i<sic__->getnthC().at(i__).at(j__).at(k__).size();__i++) 
							{
								_mth.push_back(sic__->getnthC().at(i__).at(j__).at(k__).getD(__i));
							}

							MomVecN mnt_;

							mnt_.Mom2FactN(_mth);//,abinN_,"NOEXPS");
							//cout<<"(";
							//cout<<"µ_";
							//sic__->getnthC().at(i__).at(j__).at(k__).printD();
							////mnt_.printD("NOEXPS");//abinN_,"NOEXPS");
							//cout<<")";

							//if(sic__->getnthC().at(i__).at(j__).at(k__).expo()>1){
							//	cout<<"^";
							//	cout<<sic__->getnthC().at(i__).at(j__).at(k__).expo();
							//}

							//vector<int> norder_;
							//norder_.push_back(1);
							//norder_.push_back(0);

							long double IncMom_ = 0.0L;
							long double ObsMom_ = 0.0L;
							this->CalcMoment(Centrality_,_mth,IncMom_,ObsMom_,"");


							IncMom0_ *= pow(IncMom_,sic__->getnthC().at(i__).at(j__).at(k__).expo());
							ObsMom0_ *= pow(ObsMom_,sic__->getnthC().at(i__).at(j__).at(k__).expo());


							vector<int>().swap(_mth);
						}
					}


					//cout<<"/";

					vector<int> _mth;
					_mth.push_back(0);
					_mth.push_back(1);

					MomVecN mnt__;
					mnt__.Mom2FactN(_mth);//,abinN_,"NOEXPS");

					//cout<<"(";
					////mnt__.printD("NOEXPS");//abinN_,"NOEXPS");
					//cout<<"µ_01";
					//cout<<")";


					long double IncMom___ = 0.0L;
					long double ObsMom___ = 0.0L;
					this->CalcMoment(Centrality_,_mth,IncMom___,ObsMom___,"");


					vector<int>().swap(_mth);

					//if(sic__->getMupowC().at(i__).at(j__)>1){
					//	cout<<"^"<<sic__->getMupowC().at(i__).at(j__);//<<"}";
					//}else{
					//	//cout<<"}";
					//}


					IncMom0_ /= pow(IncMom___,sic__->getMupowC().at(i__).at(j__));
					ObsMom0_ /= pow(ObsMom___,sic__->getMupowC().at(i__).at(j__));


					_IncMom__ += IncMom0_;
					_ObsMom__ += ObsMom0_;

				}
				cout<<endl;
				cout<<endl;
			}

			if(VERBOSE_=="LOUD") {
				cout<<" (SI)K_"<<mdim__;
				cout<<"("<<Centrality_<<")";
				cout<<" =["<<_IncMom__<<" "<<_ObsMom__<<"] ";

				cout<<endl;
			}

			IncMom__ = _IncMom__;
			ObsMom__ = _ObsMom__;
		}

};

#endif
