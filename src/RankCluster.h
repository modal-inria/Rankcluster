#ifndef RANKCLUSTER_H_
#define RANKCLUSTER_H_

/**@file RankCluster.h
 * @brief
 */
//#include <RcppEigen.h>
#include <vector>
#include <set>
#include "Eigen/Dense"

#include "functions.h"


// create your own data structures
struct PartialRank
{
    std::vector<int> rank;
    std::vector<int> y;
    bool isPartial;
    std::set<int> missingData;
    std::vector<int> missingIndex;

};

struct SEMparameters
{
	std::vector<int> nGibbsSE;
	std::vector<int> nGibbsM;
	int maxIt;
	int burnAlgo;
	int nGibbsL;
	int burnL;
	int maxTry;
	bool detail;
};

struct OutParameters
{
	double L;
	double bic;
	double icl;
	Eigen::ArrayXXd tik;
	Eigen::ArrayXd entropy;
	Eigen::ArrayXXd probabilities;

	//algorithm initialization
	std::vector<std::vector<std::vector<int> > > initialPartialRank;
	std::vector<std::vector<double> > initialP;
	std::vector<int> initialZ;
	std::vector<double> initialProportion;
	std::vector<std::vector<std::vector<int> > > initialMu;

	//distance between parameters
	std::vector<std::vector<double> > distProp;
	std::vector<std::vector<std::vector<double> > > distP;
	std::vector<std::vector<std::vector<int> > > distMu;
	std::vector<double> distZ;
	std::vector<std::vector<std::vector<int> > > distPartialRank;

};

class RankCluster
{
	public:
		RankCluster();
		RankCluster(std::vector<std::vector<int> > const& X,int g, std::vector<int> const& m, SEMparameters const& param);
		RankCluster(std::vector<std::vector<int> > const& X, std::vector<int> const& m, SEMparameters const& param,
                    std::vector<double> const& proportion, std::vector<std::vector<double> > const& p,
                    std::vector<std::vector<std::vector<int> > > const& mu);
		virtual ~RankCluster();
		void run();

		//getters
		inline std::vector<int> z() const {return z_;}
		inline std::vector<std::vector<double> > p() const {return p_;}
		inline std::vector<std::vector<std::vector<int> > >  mu() const {return mu_;}
		inline std::vector<double> proportion()  const  {return proportion_;}
		inline Eigen::ArrayXXd tik() const {return output_.tik;}
		inline Eigen::ArrayXd entropy() const {return output_.entropy;}
		inline Eigen::ArrayXXd probabilities() const {return output_.probabilities;}
		Eigen::ArrayXd probability() const;
		inline double bic() const {return output_.bic;}
		inline double icl() const {return output_.icl;}
		inline double L() const {return output_.L;}
		inline bool convergence() const {return convergence_;}
		inline bool partial() const {return partial_;}
		inline std::vector<std::vector<std::vector<int> > > initialPartialRank() const {return output_.initialPartialRank;}
		inline std::vector<std::vector<double> > initialP() const {return output_.initialP;}
		inline std::vector<int> initialZ() const {return output_.initialZ;}
		inline std::vector<std::vector<std::vector<int> > > initialMu() const {return output_.initialMu;}
		inline std::vector<double> initialProportion() const {return output_.initialProportion;}
		inline std::vector<std::vector<double> > distProp() const {return output_.distProp;}
		inline std::vector<std::vector<std::vector<double> > > distP() const {return output_.distP;}
		inline std::vector<std::vector<std::vector<int> > > distMu() const {return output_.distMu;}
		inline std::vector<double> distZ() const {return output_.distZ;}
		inline std::vector<std::vector<std::vector<int> > > distPartialRank() const {return output_.distPartialRank;}
		inline std::vector<std::vector<int> > indexPartialData() const {return indexPartialData_;}
		inline std::vector<int> rank(int dim, int index) const {return data_[dim][index].rank;}
        void estimateCriterion(double &L,double &bic,double &icl);

	protected: //or private
		void conversion2data(std::vector<std::vector<int> > const& X);
		void initialization();
		void SEstep();
		void gibbsY(int indexDim);
		void zSimulation();
		void gibbsX(int indexDim);
		void Mstep();
		void simuM(int indexDim,int indCl);
		void likelihood(std::vector<std::vector<std::vector<std::vector<int> > > > &listeMu,std::vector<std::vector<std::vector<double> > > &resP,
						std::vector<std::vector<double> > &resProp);
		double computeLikelihood(std::vector<std::vector<std::vector<int> > > const& mu,std::vector<std::vector<double> > const& p,
				std::vector<double> const& proportion,Eigen::ArrayXXd &tik,std::vector<std::vector<std::vector<int> > > &Y,
				std::vector<std::vector<std::vector<int> > > &xTemp, Eigen::ArrayXXd &probabilities);
		void computePartition();
		void computeDistance(std::vector<std::vector<double> > const& resProp,std::vector<std::vector<std::vector<double> > > const& resP,
				std::vector<std::vector<std::vector<std::vector<int> > > > const& resMu,std::vector<std::vector<int> > const& resZ,
				std::vector<std::vector<std::vector<std::vector<int> > > > const& resDonneesPartiel);


	private:
		std::vector<int> m_;//contains the size of rank for each dim
		int n_;//number of individuals
		int d_;//number of dimension
		int g_;//number of cluster
		std::vector<std::vector<PartialRank> > data_;
		std::vector<int> z_;
		std::vector<std::vector<std::vector<int> > > mu_;// mu_[dimension][cluster][indice]
		std::vector<std::vector<double> > p_;// p_[dimension][cluster]
		std::vector<double> proportion_;
		SEMparameters parameter_;
		OutParameters output_;
		bool partial_;//true if there is partial rank in the data
		std::vector<std::vector<int> > indexPartialData_;//index of partial data
		bool convergence_;
};

#endif /* RANKCLUSTER_H_ */
