/////
//built-in types
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

//ROOT-header
#include "TRandom.h"

//PT-header
//#include "FactMom.h"
#include "SMoment.h"

//PT-header that needs ROOT
#include "ToyModel.h"; 

//ROOT header
#include "TH1D.h"
#include "TRandom.h"


//Example of initializing Koch Ratio
int main_kochratio(){
	std::string input = "C11/C20";  //no space 
	KochRatio * kh = new KochRatio(input);

	vector<FactVec> fvv_;
	fvv_=kh->CalcFactVec();

	return 0;
}

//Example of calculating Koch Ratio for a class of events
int main()
{
	std::string input = "C11/C20";  //no space 
	KochRatio * kh = new KochRatio(input);

	vector<FactVec> fvv_;
	fvv_=kh->CalcFactVec();

	//Set the bins of efficiency corresponding to each variate (N1, \bar{N1}, N2, \bar{N2},...)
	vector<int> abinN;
	abinN.push_back(1); //(e.g. No. of protons)
	abinN.push_back(1); //(e.g. No. of anti-protons)

	abinN.push_back(1); //(e.g. No. of protons)
	abinN.push_back(1); //(e.g. No. of anti-protons)

	vector<double> eff_;
	eff_.push_back(1); //(e.g. No. of protons)
	eff_.push_back(1); //(e.g. No. of anti-protons)

	eff_.push_back(1);
	eff_.push_back(1);

	//Set the SMoment class for e-by-e analysis
	SMomentN * smt =  new SMomentN(1,fvv_,abinN);  //(centrality-bins, moment-vector, efficiency-vector)

	//Event loop
	for(int i=0; i<5e5; i++){
		vector<int> Np_;
		Np_.push_back(gRandom->Binomial(20,0.5));
		Np_.push_back(gRandom->Binomial(20,0.5));

		Np_.push_back(gRandom->Binomial(10,0.5));
		Np_.push_back(gRandom->Binomial(10,0.5));

		smt->Fill(Np_,eff_,0,1.);

		vector<int>().swap(Np_); 
		if(i%10000==0) cout<<" Event # "<<i<<endl;
	}

	long double IncMom_ = 0.0L; 
	long double IncErr_ = 0.0L;
	long double ObsMom_ = 0.0L;
	long double ObsErr_ = 0.0L;

	int Centrality =0;

	smt->CalcKochRatio(Centrality, vector<int>(), kh, IncMom_, ObsMom_,"LOUD");
	smt->CalcKochRatioError(Centrality, vector<int>(), kh, IncErr_, ObsErr_,"LOUD");

	cout<<" KR("<<Centrality<<") : Corrected Mom,Err= "<<IncMom_<<"±"<<IncErr_<<" Un-Corrected Mom,Err= "<<ObsMom_<<"±"<<ObsErr_<<endl;
}


int new_main()
{
	const int mth =3;

	CumulantVec * clt__ =new CumulantVec();//mth);
	vector<FactVec> fvv=clt__->CalcFactVec(mth);

	//Set the bins of efficiency corresponding to each variate (N1, \bar{N1}, N2, \bar{N2},...)
	vector<int> abinN;
	abinN.push_back(2); //(e.g. No. of protons)
	abinN.push_back(2); //(e.g. No. of anti-protons)

	vector<double> eff_;
	eff_.push_back(0.5); eff_.push_back(0.5);//(e.g. No. of protons)
	eff_.push_back(0.5); eff_.push_back(0.5);//(e.g. No. of anti-protons)

	//Set the SMoment class for e-by-e analysis
	SMomentN * smt =  new SMomentN(1,fvv,abinN);  //(centrality-bins, moment-vector, efficiency-vector)

	//Event loop
	for(int i=0; i<1e5; i++){
		vector<int> Np_;
		Np_.push_back(gRandom->Binomial(20,0.5));
		Np_.push_back(gRandom->Binomial(20,0.5));

		Np_.push_back(gRandom->Binomial(10,0.5));
		Np_.push_back(gRandom->Binomial(10,0.5));

		smt->Fill(Np_,eff_,0,1.);

		vector<int>().swap(Np_); 
		if(i%10000==0) cout<<" Event # "<<i<<endl;
	}

	long double IncMom_ = 0.0L; 
	long double IncErr_ = 0.0L;
	long double ObsMom_ = 0.0L;
	long double ObsErr_ = 0.0L;

	int Centrality =0;

	for(int i=0; i<3; i++){
		smt->CalcCumulant(Centrality,1,clt__,IncMom_,ObsMom_,"LOUD");
		smt->CalcCumulantError(Centrality,1,clt__,IncErr_,ObsErr_,"LOUD");
		cout<<" K("<<i<<") : Corrected Mom,Err= "<<IncMom_<<"±"<<IncErr_<<" Un-Corrected Mom,Err= "<<ObsMom_<<"±"<<ObsErr_<<endl;
	}
}

int main_example5()
{
	CumulantRatio cr_("K6/K2"); //Class to handle bi-variate cumulant ratio

	cr_.CumulantRatio2CorrMoment(); //convert to central moments of ∆N

	cr_.CalcVariance(); //Find variance in terms of central moments of ∆N
}

int main_example4()
{
	const int mth =6; //Order of cumulant of ∆N

	CumulantVec  clt__; //Class to handle bi-variate cumulants

	clt__.CalcVariance(mth); ////Find variance in terms of central moments
}



int main_example3()
{
	const int mth =6; //Order of cumulant of ∆N

	CumulantVec  clt__; //Class to handle bi-variate cumulants

	clt__.CalcVariance(mth); ////Find variance in terms of central moments
}

int main_example1()
{

	const int mth =3; //Order of cumulant of ∆N

	CumulantVec  clt__; //Class to handle bi-variate cumulants

	vector<int> abinN;  //vector for no. of efficiency bins
	abinN.push_back(2); //N       (e.g. no. of protons)
	abinN.push_back(2); //\bar{N} (e.g. no. of anti-protons)

	clt__.Cumulant2FactMoment(mth,abinN); //cumulants --> factorial moments with efficiency
}

int main_example2(){
	MomVecN mnt_; //Class to handle moments

	vector<int> mth;   //vector for multi-variate moment 
	mth.push_back(4);  //order of first variate 
	mth.push_back(2);  //order of first variate 

	mnt_.Mom2FactN(mth); //convert to factorial moment

	//mnt_.Error(mth);
}

int toymodel_main()
	//int main()
{
	vector<double> eff;  //Efficiency  array for Np (variate)
	//Set the values of efficiencies in each bins
	eff.push_back(0.5);
	//	eff.push_back(0.9);
	//	eff.push_back(0.6);

	eff.push_back(0.5);
	//	eff.push_back(0.8);
	//	eff.push_back(0.4);

	std::ofstream mfile ("sampleinput.dat", std::ofstream::out);

	//Generating events from a Toy-Model 

	vector<int> mean;

	mean.push_back(30); 
	//	mean.push_back(3);
	//	mean.push_back(5);

	mean.push_back(80); 
	//	mean.push_back(2);
	//	mean.push_back(7);


	vector<int> sigma;

	sigma.push_back(5); 
	sigma.push_back(5);



	GaussianDist dt(mean,sigma);
	//	DeltaDist dt(mean);//,sigma);

	//	PoissonDist dt(mean);

	int nbin=150; double xmin=0; double xmax=150;
	TH1D * phin1  = new TH1D ("phin1","",nbin,xmin,xmax);
	TH1D * phin2  = new TH1D ("phin2","",nbin,xmin,xmax);
	TH1D * phin3  = new TH1D ("phin3","",nbin,xmin,xmax);
	TH1D * phin4  = new TH1D ("phin4","",nbin,xmin,xmax);




	for(int i_=0;i_<1e6; i_++)
	{
		if(i_%50000==0) cout<<"No. of event = "<<i_<<endl;
		vector<int> Np_; vector<int> Npinc_;
		dt.GetEvnt(2,eff,Np_,Npinc_); //no. of variate, efficiency vector, final, initial

		phin1->Fill(Npinc_.at(0));
		phin2->Fill(Npinc_.at(1));
		phin3->Fill(Np_.at(0));
		phin4->Fill(Np_.at(1));

		for(size_t k_=0; k_<Np_.size(); k_++) mfile<<Np_.at(k_)<<" ";
		for(size_t k_=0; k_<Npinc_.size(); k_++) mfile<<Npinc_.at(k_)<<" ";
		mfile<<endl;
	}
	mfile.close();


	FILE* foutputfile = fopen("outputfile.txt","w");


	phin1->Scale(1./phin1->Integral()/((xmax-xmin)/nbin));
	phin2->Scale(1./phin2->Integral()/((xmax-xmin)/nbin));
	phin3->Scale(1./phin3->Integral()/((xmax-xmin)/nbin));
	phin4->Scale(1./phin4->Integral()/((xmax-xmin)/nbin));



	for(int i=1;i<phin1->GetNbinsX();i++){//if(phin424->GetBinContent(i)!=0)
		fprintf(foutputfile,"%g %g %g %g %g \n",phin1->GetBinCenter(i),phin1->GetBinContent(i),phin2->GetBinContent(i),phin3->GetBinContent(i),phin4->GetBinContent(i));

	}
	fclose(foutputfile);
}


int momErr_main()
	//int main()
{
	//	const int mth=6;

	vector<int> mth;
	mth.push_back(6);
	mth.push_back(0);
	//	mth.push_back(0);
	//	mth.push_back(0);
	//	mth.push_back(0);

	MomVecN mnt_;

	mnt_.Mom2FactN(mth);

	mnt_.printD("latex");

	mnt_.Error(mth);


	/*
	   vector<int> mth_;
	   mth_.push_back(2);
	//	     mth_.push_back(4);
	//	     mth_.push_back(2);

	mth_.push_back(1);
	//	     mth_.push_back(4);
	//	     mth_.push_back(2);


	MomErrorN emnt_(mth_);
	emnt_.CalcMomError();

	//	mth.push_back(2);
	//	mth.push_back(2);

	//	CentVecN cnt;
	//	cnt.Cent2FactN(mth);

	//	cnt.printD();

	//	SICumulantVec sic_;
	//
	//	sic_.SICumulant2Cumulant(mth);

	 */

}





int cumulant__main()
	//int main()
{
	//
	const int mth =3;

	CumulantVec * clt__ =new CumulantVec();//mth);

	//	clt__->Cumulant2Cumulant(mth);
	//
	//	cout<<"---------------------------------------------"<<endl;
	//        	clt__->Cumulant2Moment(mth);
	//        //
	//        //	cout<<"---------------------------------------------"<<endl;
	//        	clt__->Cumulant2CorrMoment(mth);
	//        
	//        //	cout<<"---------------------------------------------"<<endl;
	//        	clt__->Cumulant2CentralMoment(mth);
	//        //
	//        //	cout<<"---------------------------------------------"<<endl;
	// 	clt__->Cumulant2FactMoment(mth);

	//	clt__->CalcVariance(mth);


	//Set the bins of efficiency corresponding to each variate (N1, \bar{N1}, N2, \bar{N2},...)
	//THE SEQUENCE IS IMPORTANT -- Input the order like : N, \bar{N}....
	vector<int> abinN;
	//must be in multiples of two for central moments 
	abinN.push_back(2); //N1       (for example No. of protons)
	abinN.push_back(2); //\bar{N1} (for example No. of anti-protons)

	clt__->Cumulant2FactMoment(mth,abinN);


	clt__->printD("latex");

	return 0;

	vector<FactVec> fvv=clt__->CalcFactVec(mth);

	//Set the Moment class which does event-by-event averaging of the factorial moments
	SMomentN * smt =  new SMomentN(2,fvv,abinN);  //(No of centrality bins, Factorial moments vector, efficiency bin vector)


	int Centrality=0;
	/*

	//Generating events from a Toy-Model 
	vector<double> eff;  //Efficiency  array for Np (variate)
	//Set the values of efficiencies in each bins
	eff.push_back(0.80);
	eff.push_back(0.60);
	eff.push_back(0.40);

	eff.push_back(0.75);
	eff.push_back(0.55);
	eff.push_back(0.45);

	cout<<"eff=";
	for(size_t k_=0; k_<eff.size(); k_++) cout<<eff.at(k_)<<" ";
	cout<<endl;

	vector<double> mean;
	mean.push_back(10); 
	mean.push_back(5);
	mean.push_back(4);
	mean.push_back(6); 
	mean.push_back(4);
	mean.push_back(2);

	 */

	vector<int> Np;      //Multiplicity value arrary for Np

	//This part to read the event-by-event values of the variates from a file
	std::ifstream nfile ("sampleinput.dat", std::ifstream::in);
	//        	std::ifstream nfile ("toshi_3+3.dat", std::ifstream::in);
	std::string nline;

	//       	std::ofstream mfile ("sampleoutput.dat", std::ofstream::out);

	//Event loop
	long nEvnt =0L;

	//     const Double_t eff1[A1] = { 0.4, 0.9, 0.6 }; // positive
	//     const Double_t eff2[A2] = { 0.3, 0.8, 0.4 }; // negative

	vector<double> eff_;//(eff);
	eff_.push_back(0.5);
	//		eff_.push_back(0.9);
	//		eff_.push_back(0.6);

	eff_.push_back(0.5);
	//		eff_.push_back(0.8);
	//		eff_.push_back(0.4);
	//		eff_.push_back(0.40);
	//
	//		eff_.push_back(0.75);
	//		eff_.push_back(0.55);
	//		eff_.push_back(0.45);



	//eff_.push_back(0.715); eff_.push_back(0.690);
	//eff_.push_back(0.7025); eff_.push_back(0.7025);
	//eff_.swap(eff);


	vector<double> unity;
	unity.push_back(1.); 
	//		unity.push_back(1.);
	//		unity.push_back(1.);

	unity.push_back(1.);
	//		unity.push_back(1.);
	//		unity.push_back(1.);

	/*
	   PoissonDist dt(mean);
	   for(int i_=0;i_<2e6; i_++)
	   {
	   Np.clear();

	   vector<int> Np_; vector<int> Npinc_;
	   dt.GetEvnt(4,eff,Np_,Npinc_); //no. of variate, efficiency vector, final, initial

	//Np.push_back(Np_.at(0)+Np_.at(1));
	//Np.push_back(Np_.at(2)+Np_.at(3));

	smt->Fill(Np_,eff_,0);//bool(floor(gRandom->Uniform(0,2))));
	Np.clear();

	//Np.push_back(Npinc_.at(0)+Npinc_.at(1));
	//Np.push_back(Npinc_.at(2)+Npinc_.at(3));
	//			Np.swap(Np_);
	smt->Fill(Npinc_,unity,1); //refence estimation


	for(size_t k_=0; k_<Np_.size(); k_++) mfile<<Np_.at(k_)<<" ";
	mfile<<endl;

	if(nEvnt%50000==0) 
	{
	cout<<"nEvnt # "<<nEvnt;
	cout<<endl;
	cout<<"Np=";
	for(size_t k_=0; k_<Np.size(); k_++) cout<<Np.at(k_)<<" ";
	cout<<endl;
	}

	nEvnt++;

	vector<int>().swap(Np_); 
	vector<int>().swap(Npinc_); 

	}
	mfile.close();
	 */

	//Example of reading from a file
	while (std::getline(nfile, nline))
	{
		Np.clear();
		std::istringstream iss(nline);
		double nf;
		while (iss >> nf)
		{
			Np.push_back(nf);///10.);
			//Np.push_back(gRandom->Gaus(10.,5.));//Uniform(0,nf/10.));		
		}

		//		Np.resize(2);

		//cout<<Np.at(0)<<" "<<Np.at(1)<<endl;

		//for(size_t k_=0; k_<Np.size(); k_++) cout<<Np.at(k_)<<" ";
		//cout<<endl;


		vector<int> Np_;vector<int> Npinc_;
		for(size_t k_=0; k_<2; k_++) Np_.push_back(Np.at(k_));
		for(size_t k_=2; k_<Np.size(); k_++) Npinc_.push_back(Np.at(k_));
		/*
		   Np_.push_back(Np.at(0)+Np.at(1)+Np.at(2));
		   Np_.push_back(Np.at(3)+Np.at(4)+Np.at(5));
		//Np_.push_back(Np.at(0)+Np.at(1));
		//Np_.push_back(Np.at(2)+Np.at(3));

		Npinc_.push_back(Np.at(6)+Np.at(7)+Np.at(8));
		Npinc_.push_back(Np.at(9)+Np.at(10)+Np.at(11));
		//Npinc_.push_back(Np.at(4)+Np.at(5));
		//Npinc_.push_back(Np.at(6)+Np.at(7));
		 */
		//for(size_t k_=0; k_<Np_.size(); k_++) cout<<Np_.at(k_)<<" ";
		//for(size_t k_=0; k_<Npinc_.size(); k_++) cout<<Npinc_.at(k_)<<" ";
		//cout<<endl;

		smt->Fill(Np_,eff_,0);//bool(floor(gRandom->Uniform(0,2))));
		smt->Fill(Npinc_,unity,1); //refence estimation

		//		smt->Fill(Np,eff,bool(floor(gRandom->Uniform(0,2))));

		//	smt->Fill(Np,eff,Centrality);
		//	smt->Fill(Np,eff,Centrality1);

		if(nEvnt%50000==0) 
		{
			cout<<"nEvnt # "<<nEvnt;
			cout<<endl;
			cout<<"Np=";
			for(size_t k_=0; k_<Np.size(); k_++) cout<<Np.at(k_)<<" ";
			cout<<endl;
		}

		nEvnt++;

		if(nEvnt>=1e6)break;

		vector<int>().swap(Np_);vector<int>().swap(Npinc_);
		//cnt->Corr2FactN(mth);
		//cnt->Corr2FactN(mth,abinN_);
	}

	//To check what are the values of individual factorial moments enable these two lines (optional)
	for(int i__=0; i__<2*smt->size(); i__++)	smt->calcFmnN(i__,Centrality);

	long double IncMom_ = 0.0L;
	long double ObsMom_ = 0.0L;

	Centrality =0;
	smt->CalcCumulant(Centrality,1,clt__,IncMom_,ObsMom_,"LOUD");
	smt->CalcCumulant(Centrality,2,clt__,IncMom_,ObsMom_,"LOUD");
	//		smt->CalcCumulant(Centrality,3,clt__,IncMom_,ObsMom_,"LOUD");
	//		smt->CalcCumulant(Centrality,4,clt__,IncMom_,ObsMom_,"LOUD");
	//	smt->CalcCumulant(Centrality,5,clt__,IncMom_,ObsMom_,"LOUD");
	//	smt->CalcCumulant(Centrality,6,clt__,IncMom_,ObsMom_,"LOUD");

	smt->CalcCumulantError(Centrality,1,clt__,IncMom_,ObsMom_,"LOUD");
	smt->CalcCumulantError(Centrality,2,clt__,IncMom_,ObsMom_,"LOUD");
	//smt->CalcCumulantError(Centrality,3,clt__,IncMom_,ObsMom_,"LOUD");
	//smt->CalcCumulantError(Centrality,4,clt__,IncMom_,ObsMom_,"LOUD");
	//	smt->CalcCumulantError(Centrality,5,clt__,IncMom_,ObsMom_,"LOUD");
	//	smt->CalcCumulantError(Centrality,6,clt__,IncMom_,ObsMom_,"LOUD");


	Centrality =1;

	smt->CalcCumulant(Centrality,1,clt__,IncMom_,ObsMom_,"LOUD");
	smt->CalcCumulant(Centrality,2,clt__,IncMom_,ObsMom_,"LOUD");
	//		smt->CalcCumulant(Centrality,3,clt__,IncMom_,ObsMom_,"LOUD");
	//		smt->CalcCumulant(Centrality,4,clt__,IncMom_,ObsMom_,"LOUD");
	//smt->CalcCumulant(Centrality,5,clt__,IncMom_,ObsMom_,"LOUD");
	//smt->CalcCumulant(Centrality,6,clt__,IncMom_,ObsMom_,"LOUD");

	smt->CalcCumulantError(Centrality,1,clt__,IncMom_,ObsMom_,"LOUD");
	smt->CalcCumulantError(Centrality,2,clt__,IncMom_,ObsMom_,"LOUD");
	//smt->CalcCumulantError(Centrality,3,clt__,IncMom_,ObsMom_,"LOUD");
	//smt->CalcCumulantError(Centrality,4,clt__,IncMom_,ObsMom_,"LOUD");
	//smt->CalcCumulantError(Centrality,5,clt__,IncMom_,ObsMom_,"LOUD");
	//smt->CalcCumulantError(Centrality,6,clt__,IncMom_,ObsMom_,"LOUD");
}






int main_DerivCumulant()
{
	//
	const int mth =20;

	CumulantVec clt__;// =new CumulantVec(mth);

	//	clt__.Cumulant2Cumulant(mth);
	//
	//	cout<<"---------------------------------------------"<<endl;
	//	clt__.Cumulant2Moment(mth);
	//
	//	cout<<"---------------------------------------------"<<endl;
	clt__.Cumulant2CorrMoment(mth);

	//	cout<<"---------------------------------------------"<<endl;
	//	clt__.Cumulant2CentralMoment(mth);
	//
	//	cout<<"---------------------------------------------"<<endl;
	//	clt__.Cumulant2FactMoment(mth);

	vector<int> abinN;
	//must be in multiples of two for central moments 
	abinN.push_back(1); //N1       (for example No. of protons)
	abinN.push_back(1); //\bar{N1} (for example No. of anti-protons)

	//	clt__.Cumulant2FactMoment(mth,abinN);

	clt__.CalcVariance(mth); //, string FORMAT="")
}

//int main()
int main_CumulantR()
{
	CumulantRatio cr_("K8/K4");

	//vector<int> norder;
	//vector<float> nexpo;

	//norder.push_back(8);
	//norder.push_back(2);

	//nexpo.push_back(1);
	//nexpo.push_back(-1);

	//CumulantRatio cr_(norder,nexpo);


	cr_.CumulantRatio2Cumulant();
	cr_.CumulantRatio2CorrMoment();

	cr_.CalcPDeriv();
}


int CentErr_main()
	//int main()
{
	//Delta-Theorem Error estimation of central moments of any dimension & any order & any no of efficiency bins
	//
	// Define the notations of central moment
	// µ_{m,n,...} = <(∆N1-<∆N1>)^m (∆N2-<∆N2>)^n....>
	// Var(µ_{m,n,...}) = ∑_{m,m,n,n,....}
	// Error(µ_{m,n,...}) = √(Var(µ_{m,n,...})/NEvents)
	//
	//
	//set the powers of the variates of the central moment <(∆N1-<∆N1>)^{?}>
	vector<int> mth;

	//mth.push_back(0); 
	//mth.push_back(2); 

	mth.push_back(1); 
	mth.push_back(1); 

	mth.push_back(2); 
	mth.push_back(0); 

	// call central moment error estimation class
	CentErrorN * cerr = new CentErrorN(mth);

	//calculate the expression of error 
	cerr->CalcCentError("NOEXPS");

	//print to see the error expression
	cout<<"∑=";
	cerr->printE();
	//cout<<endl;
	//print the efficiency corrected error expression
	cerr->printD();//abinN);
	//
	//	return 0;
	//Find which factorial moments to be calculated
	vector<FactVec> fvv;
	fvv=cerr->CalcFactVec();//cnt__);

	//Set the bins of efficiency corresponding to each variate (N1, \bar{N1})
	//THE SEQUENCE IS IMPORTANT -- Input the order like : N, \bar{N}....
	vector<int> abinN;
	//abinN.push_back(1); //N1       (for example No. of protons)
	//abinN.push_back(1); //\bar{N1} (for example No. of anti-protons)

	//abinN.push_back(1); //N1       (for example No. of protons)
	//abinN.push_back(1); //\bar{N1} (for example No. of anti-protons)

	abinN.push_back(1); //N1       (for example No. of protons)
	abinN.push_back(1); //\bar{N1} (for example No. of anti-protons)

	cerr->printD(abinN);
	return 0;

	//		abinN.push_back(1); //N1       (for example No. of protons)
	//		abinN.push_back(1); //\bar{N1} (for example No. of anti-protons)
	//
	//
	//Set the Moment class which does event-by-event averaging of the factorial moments
	SMomentN * smt =  new SMomentN(2,fvv,abinN);  //(No of centrality bins, Factorial moments vector, efficiency bin vector)
	//Specify if the analysis is to be done for a specific centrality (for multiple centrality not needed)
	//
	const int Centrality=0;

	vector<double> eff;  //Efficiency  array for Np (variate)
	//Set the values of efficiencies in each bins
	//eff.push_back(0.5);
	//eff.push_back(0.5);

	//eff.push_back(1.);
	//eff.push_back(1.);

	eff.push_back(0.998);
	eff.push_back(0.999);

	cout<<"eff=";
	for(size_t k_=0; k_<eff.size(); k_++) cout<<eff.at(k_)<<" ";
	cout<<endl;

	vector<int> Np;      //Multiplicity value arrary for Np

	//This part to read the event-by-event values of the variates from a file
	std::ifstream nfile ("sampleinput.dat", std::ifstream::in);
	std::string nline;

	//Event loop
	long nEvnt =0L;
	while (std::getline(nfile, nline))
	{
		Np.clear();
		std::istringstream iss(nline);
		double nf;
		while (iss >> nf)
		{
			//		Np.push_back(nf/10.);
			//		Np.push_back(nf/10.);
			//		Np.push_back(nf/10.);
			//Np.push_back(gRandom->Uniform(0,nf/10.));		
			//Np.push_back(gRandom->Uniform(0,nf/10.));		
			//Np.push_back(gRandom->Uniform(0,nf/10.));		
			//Np.push_back(gRandom->Uniform(0,nf/10.));		
			//Np.push_back(gRandom->Uniform(0,nf/10.));		
			//Np.push_back(gRandom->Uniform(0,nf/10.));		
			Np.push_back(gRandom->Gaus(10.,5.));//Uniform(0,nf/10.));               

		}

		Np.resize(2);


		smt->Fill(Np,eff,bool(floor(gRandom->Uniform(0,2))));

		//	smt->Fill(Np,eff,Centrality);
		//	smt->Fill(Np,eff,Centrality1);

		if(nEvnt%10000==0) 
		{
			cout<<"nEvnt # "<<nEvnt;
			cout<<endl;
			cout<<"Np=";
			for(size_t k_=0; k_<Np.size(); k_++) cout<<Np.at(k_)<<" ";
			cout<<endl;
		}

		nEvnt++;

		//cnt->Corr2FactN(mth);
		//cnt->Corr2FactN(mth,abinN_);
	}

	//To check what are the values of individual factorial moments enable these two lines (optional)
	for(int i__=0; i__<2*smt->size(); i__++)	smt->calcFmnN(i__,Centrality);


	long double IncMom_ = 0.0L;
	long double ObsMom_ = 0LL;

	smt->CalcCorrError(0, mth, cerr, IncMom_, ObsMom_,"CHECK");

	for(int i=1; i<=6; i++){
		mth.clear();
		mth.push_back(i);
		mth.push_back(i);
		smt->CalcCorrMom(0, mth, cerr, IncMom_, ObsMom_,"");
	}

	//		smt->CalcCorrMom(1, mth, cerr, IncMom_, ObsMom_,"");
	//		smt->CalcCorrError(1, mth, cerr, IncMom_, ObsMom_,"");

}



//SIC simulation with toy-model event generator
int SIC_toymodel_main()
{
	int norder=8;

	SICumulantVec sic_;

	sic_.SICumulant2Moment(norder);

	vector<int> abinN;
	abinN.push_back(2); //N1       (for example No. of protons)
	abinN.push_back(2); //\bar{N1} (for example No. of anti-protons)

	//sic_.printD();
	//sic_.printD(abinN);

	vector<FactVec> fvv= sic_.CalcFactVec(norder); //, string FORMAT="")

	SMomentN * smt =  new SMomentN(2,fvv,abinN);  //(No of centrality bins, Factorial moments vector, efficiency bin vector)

	vector<int> MeanP;
	MeanP.push_back(10);
	MeanP.push_back(10);
	MeanP.push_back(10);
	MeanP.push_back(10);

	vector<double> eff;
	vector<double> unit;
	vector<int> Np;
	vector<int> Npinc;

	eff.push_back(0.4);
	eff.push_back(0.4);
	eff.push_back(0.4);
	eff.push_back(0.4);

	unit.push_back(1.);
	unit.push_back(1.);
	unit.push_back(1.);
	unit.push_back(1.);

	DeltaDist dt(MeanP);

	for(int l__=0;l__<1e5; l__++){
		dt.GetEvnt(4,eff,Np,Npinc);
		smt->Fill(Np,eff,0);
		smt->Fill(Npinc,unit,1);
	}

	//To check what are the values of individual factorial moments enable these two lines (optional)
	for(int i__=0; i__<2*smt->size(); i__++)	smt->calcFmnN(i__,0);
	for(int i__=0; i__<2*smt->size(); i__++)	smt->calcFmnN(i__,1);

	long double IncMom_ = 0.0L;
	long double ObsMom_ = 0.0L;

	//sic_.printD();

	for(int i=1; i<=norder; i++)smt->CalcSICumulant(0,i,&sic_,IncMom_,ObsMom_,"LOUD");
	for(int i=1; i<=norder; i++)smt->CalcSICumulant(1,i,&sic_,IncMom_,ObsMom_,"LOUD");
}




int SIC_data_main()
{

	int norder=5;

	SICumulantVec sic_;

	//sic_.SICumulant2SICumulant(3);

	sic_.SICumulant2Moment(norder);

	vector<int> abinN;
	abinN.push_back(1); //N1       (for example No. of protons)
	abinN.push_back(1); //\bar{N1} (for example No. of anti-protons)

	//sic_.printD();
	//sic_.printD(abinN);

	vector<FactVec> fvv= sic_.CalcFactVec(norder); //, string FORMAT="")

	SMomentN * smt =  new SMomentN(2,fvv,abinN);  //(No of centrality bins, Factorial moments vector, efficiency bin vector)

	vector<double> eff;  //Efficiency  array for Np (variate)
	//Set the values of efficiencies in each bins
	eff.push_back(0.8);
	eff.push_back(0.8);

	//		eff.push_back(1.);
	//		eff.push_back(1.);

	cout<<"eff=";
	for(size_t k_=0; k_<eff.size(); k_++) cout<<eff.at(k_)<<" ";
	cout<<endl;

	vector<int> Np;      //Multiplicity value arrary for Np

	//This part to read the event-by-event values of the variates from a file
	std::ifstream nfile ("sampleinput.dat", std::ifstream::in);
	std::string nline;

	//Event loop
	long nEvnt =0L;
	while (std::getline(nfile, nline))
	{
		Np.clear();
		std::istringstream iss(nline);
		double nf;
		while (iss >> nf)
		{
			//		Np.push_back(nf/10.);
			//		Np.push_back(nf/10.);
			//		Np.push_back(nf/10.);
			Np.push_back(gRandom->Uniform(0,nf/10.));		
			Np.push_back(gRandom->Uniform(0,nf/10.));		
			Np.push_back(gRandom->Uniform(0,nf/10.));		
		}

		Np.resize(2);


		smt->Fill(Np,eff,bool(floor(gRandom->Uniform(0,2))));

		if(nEvnt%10000==0) 
		{
			cout<<"nEvnt # "<<nEvnt;
			cout<<endl;
			cout<<"Np=";
			for(size_t k_=0; k_<Np.size(); k_++) cout<<Np.at(k_)<<" ";
			cout<<endl;
		}

		nEvnt++;
	}

	//To check what are the values of individual factorial moments enable these two lines (optional)
	for(int i__=0; i__<2*smt->size(); i__++)	smt->calcFmnN(i__,0);
	for(int i__=0; i__<2*smt->size(); i__++)	smt->calcFmnN(i__,1);

	long double IncMom_ = 0.0L;
	long double ObsMom_ = 0.0L;


	//smt->CalcMoment(Centrality_, const vector<int> norder_, long double& IncMom_, long double& ObsMom_,const string VERBOSE_="")
	//vector<int> norder_;
	//norder_.push_back(1);
	//norder_.push_back(0);

	//smt->CalcMoment(0,norder_,IncMom_,ObsMom_,"LOUD");

	sic_.printD();

	for(int i=1; i<=norder; i++)smt->CalcSICumulant(0,i,&sic_,IncMom_,ObsMom_,"LOUD");

	//smt->CalcCumulant(Centrality,1,clt__,IncMom_,ObsMom_,"LOUD");
	//smt->CalcCumulant(Centrality,2,clt__,IncMom_,ObsMom_,"LOUD");
	//smt->CalcCumulant(Centrality,3,clt__,IncMom_,ObsMom_,"LOUD");
	//smt->CalcCumulant(Centrality,4,clt__,IncMom_,ObsMom_,"LOUD");
	//smt->CalcCumulant(Centrality,5,clt__,IncMom_,ObsMom_,"LOUD");
	//smt->CalcCumulant(Centrality,6,clt__,IncMom_,ObsMom_,"LOUD");



}

int SIC_main()
{
	const int mth =3;
	SICumulantVec sic_;

	//sic_.SICumulant2SICumulant(3);

	sic_.SICumulant2Moment(mth);

	vector<int> abinN;
	abinN.push_back(6); //N1       (for example No. of protons)
	abinN.push_back(2); //\bar{N1} (for example No. of anti-protons)

	//sic_.printD();
	sic_.printD(abinN);

	vector<FactVec> fvv= sic_.CalcFactVec(3); //, string FORMAT="")

}


























//main program for ratio fluctuation analysis 
//
int nudyn_main()
	//int main()
{
	//Calculation of Nudyn \nu^{dyn}_{n}  //n-th order dynmical scaled fluctuation by
	//Pruneau, S. Gavin, and S. Voloshin & Christiansen, E. Haslum, and E. Stenlund
	//More information on definition : Phys Rev. C 80, 034903 (2009), Phys. Rev. C 66, 044904 (2002). 
	//
	const int mth =3;

	//Nudynamic moments vector-space class
	NudynVecN * ndyn__ =new NudynVecN();
	//
	//
	//Convert Nudynamic into Factorial moments F_{m,n,..} = <N1(N1-1)..(N1-m+1) N2(N2-1)..(N2-n+1)...>
	ndyn__->Nudyn2FactN(mth);
	ndyn__->printD();  // print to see the expression (optional)

	cout<<"No of factorial moments : ";
	cout<<ndyn__->getCCmn().size()<<endl;

	//Find which factorial moments to be calculated
	vector<FactVec> fvv;
	fvv=ndyn__->CalcFactVec();//cnt__);

	//Set the bins of efficiency corresponding to each variate (N1, \bar{N1})
	//THE SEQUENCE IS IMPORTANT -- Input the order like : N, \bar{N}....
	vector<int> abinN;
	abinN.push_back(2); //N1       (for example No. of protons)
	abinN.push_back(2); //\bar{N1} (for example No. of anti-protons)
	//
	return 0;
	//
	//Set the Moment class which does event-by-event averaging of the factorial moments
	SMomentN * smt =  new SMomentN(2,fvv,abinN);  //(No of centrality bins, Factorial moments vector, efficiency bin vector)


	//Specify if the analysis is to be done for a specific centrality (for multiple centrality not needed)
	const int Centrality=0;
	const int Centrality1=1;

	vector<double> eff;  //Efficiency  array for Np (variate)
	//Set the values of efficiencies in each bins
	eff.push_back(0.8);
	eff.push_back(0.8);
	eff.push_back(0.8);
	eff.push_back(0.8);
	//eff.push_back(0.8);
	//eff.push_back(0.8);
	//eff.push_back(0.8);
	//eff.push_back(0.8);
	//eff.push_back(0.8);
	//eff.push_back(0.8);
	//eff.push_back(0.8);
	//eff.push_back(0.8);

	cout<<"eff=";
	for(size_t k_=0; k_<eff.size(); k_++) cout<<eff.at(k_)<<" ";
	cout<<endl;

	vector<int> Np;      //Multiplicity value arrary for Np

	//This part to read the event-by-event values of the variates from a file
	std::ifstream nfile ("numbers.dat", std::ifstream::in);
	std::string nline;

	//Event loop
	while (std::getline(nfile, nline))
	{
		Np.clear();
		std::istringstream iss(nline);
		double nf;
		while (iss >> nf)
		{
			Np.push_back(nf/10.);
			Np.push_back(nf/10.);
			Np.push_back(nf/10.);
		}

		Np.resize(4);


		//using two centrality bins bool 0->centrality , 1->centrality1
		smt->Fill(Np,eff,bool(floor(gRandom->Uniform(0,2))));

		//	smt->Fill(Np,eff,Centrality);
		//	smt->Fill(Np,eff,Centrality1);

		cout<<endl;
		cout<<"Np=";
		for(size_t k_=0; k_<Np.size(); k_++) cout<<Np.at(k_)<<" ";
		cout<<endl;

		cout<<endl;

		//cnt->Corr2FactN(mth);
		//cnt->Corr2FactN(mth,abinN_);
	}

	//To check what are the values of individual factorial moments enable these two lines (optional)
	for(int i__=0; i__<2*smt->size(); i__++)	smt->calcFmnN(i__,Centrality);
	for(int i__=0; i__<2*smt->size(); i__++)	smt->calcFmnN(i__,Centrality1);

	long double IncMom_ = 0.0L;
	long double ObsMom_ = 0LL;

	smt->CalcNudynMom(Centrality,mth, ndyn__, IncMom_, ObsMom_);
	smt->CalcNudynMom(Centrality1,mth, ndyn__, IncMom_, ObsMom_);

}


//main program for higher moment analysis
//
int corr_main()
//int main()
{
	//Calculation of Correlation variables Ç_{i,j,...} = <(∆N1-<∆N1>)^{i} (∆N2-<∆N2>)^{j}...>
	//with ∆N1 = N1 - \bar{N1}, ∆N2 = N2 - \bar{N2}, .......
	//
	//set the powers of the variates of the central moment <(∆N1-<∆N1>)^{?} (∆N2-<∆N2>)^{?}...>
	vector<int> mth;
	mth.push_back(2); 
	mth.push_back(2); 
//	mth.push_back(2); 
//	mth.push_back(2); 
	//mth.push_back(2); 
	//mth.push_back(2); 

	//Central moment vector-space class
	CentVecN * cnt__ =new CentVecN();
	//
	//Convert the Correlation variable into Central moments C_{m,n,..} = <(N1-<N1>)^m (N2-<N2>)^n...>, 
	cnt__->Corr2CentN(mth);
	//	cnt__->printD();  // print to see the expression (optional)
	//
	//Convert the Correlation variable into Factorial moments F_{m,n,..} = <N1(N1-1)..(N1-m+1) N2(N2-1)..(N2-n+1)...>
	//
	//

	cnt__->Corr2FactN(mth);  

//	cnt__->Squeeze();
	//	cnt__->printD();  // print to see the expression (optional)


	cout<<"No of central moments : ";
	cout<<cnt__->getCCmn().size()<<endl;

	//Find which factorial moments to be calculated
	vector<FactVec> fvv;
	fvv=cnt__->CalcFactVec();//cnt__);

	//Set the bins of efficiency corresponding to each variate (N1, \bar{N1}, N2, \bar{N2},...)
	//THE SEQUENCE IS IMPORTANT -- Input the order like : N, \bar{N}....
	vector<int> abinN;
	//must be in multiples of two for central moments 
	abinN.push_back(1); //N1       (for example No. of protons)
	abinN.push_back(1); //\bar{N1} (for example No. of anti-protons)
	//
	abinN.push_back(1); //N2        (for example No. of kaons)
	abinN.push_back(1); //\bar{N2}  (for example No. of anti-kaons)
//	//	//
//	abinN.push_back(1); //N3        (for example No. of pions)
//	abinN.push_back(1); //\bar{N3}  (for example No. of anti-pions)
	//
	//repeat the order

	cnt__->Corr2FactN(mth,abinN);  

	return 0;
	//
	//Set the Moment class which does event-by-event averaging of the factorial moments
	SMomentN * smt =  new SMomentN(2,fvv,abinN);  //(No of centrality bins, Factorial moments vector, efficiency bin vector)

	//Specify if the analysis is to be done for a specific centrality (for multiple centrality not needed)
	const int Centrality=0;
	const int Centrality1=1;

	vector<double> eff;  //Efficiency  array for Np (variate)
	//Set the values of efficiencies in each bins
	eff.push_back(1.);
	eff.push_back(1.);
	//eff.push_back(0.8);
	//eff.push_back(0.8);
	//eff.push_back(0.8);
	//eff.push_back(0.8);
	//eff.push_back(0.8);
	//eff.push_back(0.8);
	//eff.push_back(0.8);
	//eff.push_back(0.8);
	//eff.push_back(0.8);
	//eff.push_back(0.8);

	cout<<"eff=";
	for(size_t k_=0; k_<eff.size(); k_++) cout<<eff.at(k_)<<" ";
	cout<<endl;

	vector<int> Np;      //Multiplicity value arrary for Np

	//This part to read the event-by-event values of the variates from a file
	std::ifstream nfile ("sampleinput.dat", std::ifstream::in);
	std::string nline;

	long nEvnt=0;

	//Event loop
	while (std::getline(nfile, nline))
	{
		Np.clear();
		std::istringstream iss(nline);
		double nf;
		while (iss >> nf)
		{
			//Np.push_back(nf/10.);
			//Np.push_back(nf/10.);
			//Np.push_back(nf/10.);

        			Np.push_back(gRandom->Uniform(0,nf/10.));		
        			Np.push_back(gRandom->Uniform(0,nf/10.));		
        			Np.push_back(gRandom->Uniform(0,nf/10.));		
		}

		Np.resize(2);


		smt->Fill(Np,eff,bool(floor(gRandom->Uniform(0,2))));

		//	smt->Fill(Np,eff,Centrality);
		//	smt->Fill(Np,eff,Centrality1);

		if(nEvnt%10000==0) 
		{
			cout<<"nEvnt # "<<nEvnt;
			cout<<endl;
			cout<<"Np=";
			for(size_t k_=0; k_<Np.size(); k_++) cout<<Np.at(k_)<<" ";
			cout<<endl;
		}

		nEvnt++;
		//cnt->Corr2FactN(mth);
		//cnt->Corr2FactN(mth,abinN_);
	}

	//To check what are the values of individual factorial moments enable these two lines (optional)
	for(int i__=0; i__<2*smt->size(); i__++)	smt->calcFmnN(i__,Centrality);
	for(int i__=0; i__<2*smt->size(); i__++)	smt->calcFmnN(i__,Centrality1);

	long double IncMom_ = 0.0L;
	long double ObsMom_ = 0.0L;

	//CalcCorrMom(Centrality, mth, cnt__, smt, IncMom_, ObsMom_);
	//CalcCorrMom(Centrality1, mth, cnt__, smt, IncMom_, ObsMom_);

	//for(int i__=0; i__<smt->centN(); i__++)	smt->CalcCorrMom(Centrality, mth, cnt__, IncMom_, ObsMom_);
//	for(int i__=0; i__<smt->centN(); i__++)	
	smt->CalcCorrMom(Centrality, mth, cnt__, IncMom_, ObsMom_,"LOUD");
}




int main_readstring(){
	std::string input = "K4/K2";  //no space 
	std::stringstream ss;
	ss << input;
	int found;
	std::string temp;
	
	cout<<"input="<<input<<endl; 

	int counter=0;

	while(std::getline(ss, temp,'/')) {

		if(counter==0/* && temp.size()>2*/) 
		{
			std::stringstream uu;
			std::string temp_;
			uu << temp;
			int counter_=0;
			while(std::getline(uu, temp_,'C'))
			{
				//cout<<uu.str()<<" "<<temp_<<"--"<<endl;

				int order_ = atoi(temp_.c_str() +0);
				float expo_ =0;
				if(temp_.size()<3){ expo_=1;
				}else{expo_ = atof(temp_.c_str() + 2);
				}
				//if(counter_!=0) cout<<"temp_="<<temp_<<" order_="<<order_<<" expo_="<<expo_<<" temp_size="<<temp_.size()<<endl;
				if(counter_!=0) cout<<" order="<<order_<<" expo="<<expo_<<endl;

				counter_ ++;
			}
		}else{

			int order = atoi(temp.c_str() +1);
			float expo =0;
			if(temp.size()<3){ expo=-1;
			}else{expo = -atof(temp.c_str() + 3);
			}
			//cout<<temp<<" "<<order<<" "<<expo<<" "<<temp.size()<<endl;
			//cout<<"temp="<<temp<<" order="<<order<<" expo="<<expo<<" temp_size="<<temp.size()<<endl;
			cout<<" order="<<order<<" expo="<<expo<<endl;
		}
		counter ++;
	}
	return 0;
}
