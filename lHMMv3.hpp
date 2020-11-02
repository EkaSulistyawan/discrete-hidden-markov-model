
/**
 * Hidden Markov Model - Baum Welch Re-estimation (Expectation-Maximization) Algorithm
 * 3rd
 * using variable type "Double"
 * the equation has been adapted through logarithmic operation
 * */


#include <iostream>
#include <vector>
#include <math.h>
#include <stdexcept>

#ifndef LHMM_H
#define LHMM_H
#endif

/**
 * lHMM.h
 * For the sake of simplicity,
 * I put the header (h) and cpp in a same file
 * */

class lHMM
{
private:
    const double LOGZERO = -INFINITY;
    //  initial parameter
    int nData = 4;  //  the number of state
    int nState;
    std::vector<std::vector<double> > T;//   T[i][j] transition state from i to 
    std::vector<std::vector<double> > E;
    std::vector<double> phi;
    double gprob;
    //  parameter update
    std::vector<std::vector<double> > gAlpha;
    std::vector<std::vector<double> > gBeta;
    double gamma(std::vector<char>& _dna,int i, int t);
    double ksi(std::vector<char>& _dna,int i, int j, int t);
    int nt(char t);
    //  pratical HMM implementation
    //  based on
    //  [1] Numerically Stable Hidden Markov Model Implementation by Tobias P. Mann
    //  [2] Durbin R., Eddy., Krogh A., Mitchison G. Biological Sequence Analysis. Cambridge University Press, 1998.
    double eexp(double x);
    double eln(double x);
    double elnsum(double x, double y);      //  eln(x) + eln(y)
    double elnproduct(double x,double y);   //  eln(x) * eln(y)

public:
    //  constructor
    lHMM(int _nstate,std::vector<std::vector<double> > _T,std::vector<std::vector<double> > _E,std::vector<double> _phi);

    //  setting parameter
    void setParameter(int _nstate,std::vector<std::vector <double> > _T,std::vector<std::vector <double> > _E,std::vector<double> _phi);
    
    /*  
        The Forward Algorithm used to solve the first kind of HMM problem, which also used in the third kind.
    */

    void forward(std::vector<char>& _dna,int panjangdna);
    double probObsrvGivenLmd(int panjangdna);
    void backward(std::vector<char>& _dna,int panjangdna);    
    void update(std::vector<char>& _dna,int panjangdna);
    void learn(std::vector<char>& _dna,int panjangdna, double terminProb, int terminIter);
    void fstUpdate(std::vector<char>& _dna,int panjangdna);
    void getTEphi();
    double getProb();
    std::vector<std::vector<double> > getT();   // T[i][j] State dari i ke j
    std::vector<std::vector<double> > getE();
    std::vector<double> getphi();
};

lHMM::lHMM(int _nstate,std::vector<std::vector<double> > _T,std::vector<std::vector<double> > _E,std::vector<double> _phi){
    setParameter(_nstate,_T,_E,_phi);
}

void lHMM::setParameter(int _nstate,std::vector<std::vector <double> > _T,std::vector<std::vector <double> > _E,std::vector<double> _phi){
    nState = _nstate;
    T = _T;
    E = _E;
    phi = _phi;
    gprob = LOGZERO;
}



/**
 * Simple Implementation,
 * separated into Forward and Backward Algo.
 * */

/**
 * 
 * 
 * 
 * The Expectation
 * 
 * 
 * 
 * */

void lHMM::forward(std::vector<char>& _dna, int panjangdna){
    //  reset alpha
    std::vector<double> row(nState,0);
    std::vector<std::vector<double> > alpha(panjangdna,row);

    //  init data
    //  computing initial probability
    for (int i=0;i<nState;i++){
        alpha[0][i] = elnproduct(eln(phi[i]),eln(E[nt(_dna[0])][i]));//->dalam bentuk log
        //printf("%c %.10f ",_dna[0],eexp(alpha[0][i]));
    }
    //printf("\n");

    //  forward
    
    for(int i=1;i<panjangdna;i++){
        for (int j=0;j<nState;j++){
            double sum = LOGZERO;
            for (int k=0;k<nState;k++){
                //  using non-zero
                sum = elnsum(sum,elnproduct(eln(T[k][j]),alpha[i-1][k]));
            }
            alpha[i][j] = elnproduct(sum,eln(E[nt(_dna[i])][j]));
            if((alpha[i][j]==INFINITY)||(alpha[i][j]==-INFINITY)){
                //printf("Infinite Encountered! \n");
                break;break;break;
            }
        }
        //printf("\n");
    }
    //  print last value
    //std::cout<<"Forward Done!\n";
    gAlpha = alpha;
}

/*


    BACKWARD ALGORITHM


*/
void lHMM::backward(std::vector<char>& _dna,int panjangdna){
    //  reset beta
    std::vector<double> row(nState,0);
    std::vector<std::vector<double> > beta(panjangdna,row);

    //  initialize beta
    for (int i=0;i<nState;i++){
        beta[panjangdna-1][i] = eln(1);
        //std::cout<<_dna[panjangdna-1]<<" "<<beta[panjangdna-1][i]<<" ";
    }
    //std::cout<<"\n";

    for (int i = panjangdna-1;i>0;i--){
        for (int j = 0;j<nState;j++){
            double sum=LOGZERO;
            for (int k=0;k<nState;k++){
                sum = elnsum(sum,elnproduct(elnproduct(beta[i][k],eln(E[nt(_dna[i])][k])),eln(T[j][k])));
                //std::cout<<sum<<"\n";
            }
            //std::cout<<"\n";
            beta[i-1][j] = sum;
            //printf("%c %.10f\n",_dna[i],eexp(beta[i][j]));
            //std::cout<<"Sum "<<sum;
            if((beta[i-1][j]==INFINITY)||(beta[i-1][j]==-INFINITY)){
                //printf("Infinite Encountered! \n");
                break;break;break;
            }
        }
        //std::cout<<"\n";
    }

    //std::cout<<"Backward Done!\n";
    gBeta = beta;
}

std::vector<std::vector<double> > lHMM::getT(){
    return T;
}
std::vector<std::vector<double> > lHMM::getE(){
    return E;
}
std::vector<double> lHMM::getphi(){
    return phi;
}

/*


    The Maximization


*/

double lHMM::gamma(std::vector<char>& _dna,int i, int t){//---> in a form of logarithmic
    /*
        double denum = LOGZERO;
        for (int k=0;k<nState;k++){//   k = state
            denum = elnsum(denum,elnproduct(gAlpha[t][k],gBeta[t][k]));
        }
        double num = elnproduct(gAlpha[t][i],gBeta[t][i]);
    */
   double res = LOGZERO;
   for (int j =0;j<nState;j++){
       res = elnsum(res,ksi(_dna,i,j,t));
   }
    return res;
}

double lHMM::ksi(std::vector<char>& _dna,int i, int j, int t){
    double denum = LOGZERO;
    for (int k=0;k<nState;k++){
        for (int l=0;l<nState;l++){
            double part1 = elnproduct(gAlpha[t][k],eln(T[k][l]));
            double part2 = elnproduct(gBeta[t+1][l],eln(E[nt(_dna[t+1])][l]));
            denum = elnsum(denum,elnproduct(part1,part2));
        }
    }
    double part1 = elnproduct(gAlpha[t][i],eln(T[i][j]));
    double part2 = elnproduct(gBeta[t+1][j],eln(E[nt(_dna[t+1])][j]));
    double num = elnproduct(part1,part2);

    return elnproduct(num,-denum);
}

void lHMM::update(std::vector<char>& _dna,int panjangdna){
    //  get new phi
    //std::cout<<"New PHI\n";
    std::vector<double> NEWphi(6,0);
    for (int i =0;i<nState;i++){
        NEWphi[i] = eexp(gamma(_dna,i,0));
        //std::cout<<NEWphi[i]<<" ";
    }
    //std::cout<<"\n";
    //  get new transition
    std::vector<double> row(nState,0);
    std::vector<std::vector<double> > NEWT(nState,row);
    //std::cout<<"New Transition\n";
    for (int i=0;i<nState;i++){
        for (int j=0;j<nState;j++){
            //  init
            double num = LOGZERO;
            double denum = LOGZERO;
            //  sigma num&denum
            for (int k=0;k<panjangdna-1;k++){
                denum = elnsum(denum,gamma(_dna,i,k));
                num = elnsum(num,ksi(_dna,i,j,k));
            }
            NEWT[i][j] = eexp(elnproduct(num,-denum));
            //std::cout<<NEWT[i][j]<<" ";
        }
        //std::cout<<"\n";
    }
    //  get new emission
    //std::cout<<"New Emission\n";
    std::vector<std::vector<double> > NEWE(nData,row);
    double num = LOGZERO;
    double denum = LOGZERO;
    for (int i =0;i<nData;i++){
        for (int j=0;j<nState;j++){//--> vk, possible values 
            double num = LOGZERO;
            double denum = LOGZERO;

            for (int t = 0;t<panjangdna-1;t++){
                if(nt(_dna[t])==i){
                    num = elnsum(num,gamma(_dna,j,t));
                }
                denum = elnsum(denum,gamma(_dna,j,t));
            }

            NEWE[i][j] = eexp(elnproduct(num,-denum));
            //std::cout<<NEWE[i][j]<<" ";
        }
        //std::cout<<"\n";
    }

    //  renew
    T = NEWT;
    E = NEWE;
    phi = NEWphi;
}

/**
 * 
 * Return the P(O|lambda) for a single Forward Process
 * 
 * */

double lHMM::probObsrvGivenLmd(int panjangdna){
    double res = LOGZERO;
    for (int i=0;i<nState;i++){
        res = elnsum(res,gAlpha[panjangdna-1][i]);
        //std::cout<<gAlpha[panjangdna-1][i]<<" ";
    }
    //  res = eexp(res);
    //  printf("\nCurrent Prob(ln) : %f\n",res);
    return res;
}


/**
 * 
 * learn() -> semi-advance
 * include the termination process
 * 
 * */

void lHMM::learn(std::vector<char>& _dna,int panjangdna, double terminProb, int terminIter){
    double prob = LOGZERO;
    int iter = 0;

    while((prob < eln(terminProb))&&(iter < terminIter)){
        forward(_dna,panjangdna); // -> single procedure statement for Forward
        gprob = probObsrvGivenLmd(panjangdna);
        if(!(prob < eln(terminProb))){break;}
        backward(_dna,panjangdna); // -> single procedure statement for Backward
        update(_dna,panjangdna); // -> single procedure statement for Maximization
        iter++;
    }
}

/**
 * fstUpdate -> Fast Update
 * High memory, reusable variable for the next process
 * such as, storing Forward value instead calculating it through recursive algo.
 * */

void lHMM::fstUpdate(std::vector<char>& _dna, int panjangdna){
    //  backward forward in one run
    //  var declaration 
    
    std::vector<double> row(nState,0.0);
    std::vector<std::vector<double> > alpha(panjangdna,row);
    std::vector<std::vector<double> > beta(panjangdna,row);
    int t=0;
    int revt = panjangdna-2; //-> to reverse time

    //  initialization
    for (int i=0;i<nState;i++){
        alpha[0][i] = elnproduct(eln(phi[i]),eln(E[nt(_dna[0])][i])); 
        beta[panjangdna-1][i] = eln(1);
    }  

    //  induction
    for(int i=1;(i<panjangdna)&&(gprob<0);i++){
        t = i;
        revt = (panjangdna-1)-i;
        for (int j=0;((j<nState)&&(gprob<0));j++){
            double sumAlpha = LOGZERO;
            double sumBeta = LOGZERO;
            for (int k=0;((k<nState)&&(gprob<0));k++){
                sumAlpha    = elnsum(sumAlpha,elnproduct(eln(T[k][j]),alpha[t-1][k]));
                sumBeta     = elnsum(sumBeta,elnproduct(beta[revt+1][k],elnproduct(eln(E[nt(_dna[revt+1])][k]),eln(T[j][k])))); 
            }
            alpha[t][j] = elnproduct(sumAlpha,eln(E[nt(_dna[t])][j]));
            beta[revt][j] = sumBeta;

            //  cout
            //printf("%c %.10f\n",_dna[revt],eexp(beta[revt][j]));
            if((alpha[t][j]==INFINITY)/*||(alpha[t][j]==-INFINITY)*/){
                printf("Alpha Infinite Encountered! %f\n",alpha[t][j]);
                gprob = 2;
                throw std::overflow_error("Infinite Alpha While Updating");
            }else if((beta[revt][j]==INFINITY)/*||(beta[revt][j]==-INFINITY)*/){
                printf("Beta Infinite Encountered!\n");
                gprob = 2;
                throw std::overflow_error("Infinite Beta While Updating");
            }
        }
    }
    //  count prob
    double res = LOGZERO;
    for (int i=0;i<nState;i++){
        res = elnsum(res,alpha[panjangdna-1][i]);
        //std::cout<<gAlpha[panjangdna-1][i]<<" ";
    }
    if(gprob!=2){gprob = res;}


    //  ksi -> gamma -> update
    std::vector<double> row2(nState,LOGZERO);
    std::vector<double> row3(nState,LOGZERO);
    std::vector<std::vector<double> > SumNumKsi(nState,row2);
    std::vector<std::vector<double> > SumNumObj(nData,row3);
    std::vector<double> SumNumGamma(nState,LOGZERO);
    for (int i = 0;i<panjangdna-1;i++){//-> besar T
        //  count for all k and l
        t=i;
        std::vector<std::vector<double> >numKsi(nState,row);
        std::vector<double> numGamma(nState,0.0);
        double denum = LOGZERO;
        for (int k=0;k<nState;k++){
            double sumGamma = LOGZERO;
            for (int l=0;l<nState;l++){
                //  count denum
                double part1 = elnproduct(alpha[t][k],beta[t+1][l]);
                double part2 = elnproduct(eln(T[k][l]),eln(E[nt(_dna[t+1])][l]));
                denum = elnsum(denum,elnproduct(part1,part2));
                //  count numerator KSI
                numKsi[k][l] = elnproduct(part1,part2);
                //  count sumGamma
                sumGamma = elnsum(sumGamma,elnproduct(part1,part2));
            }
            numGamma[k] = sumGamma;
        }

        if (i==0){
            for (int k=0;k<nState;k++){
                //  discomment to update PHI
                phi[k] = eexp(elnproduct(numGamma[k],-denum));
            }
        }

        //  sum all over
        for (int k=0;k<nState;k++){
            //  sum all KSI
            for (int l=0;l<nState;l++){
                SumNumKsi[k][l] = elnsum(SumNumKsi[k][l],elnproduct(numKsi[k][l],-denum));
                //  Overflow Flag
                if(SumNumKsi[k][l]==INFINITY){
                    throw std::overflow_error("SUM Numerator Ksi reach INFINITE!");
                }
            }
            //  sum selective 
            SumNumObj[nt(_dna[t])][k] = elnsum(SumNumObj[nt(_dna[t])][k],elnproduct(numGamma[k],-denum));
            if(SumNumObj[nt(_dna[t])][k]==INFINITY){
                throw std::overflow_error("SUM Numerator Obj reach INFINITE!");
            }
            //  Sum all gamma
            SumNumGamma[k] = elnsum(SumNumGamma[k],elnproduct(numGamma[k],-denum));
            if(SumNumGamma[k]==INFINITY){
                throw std::overflow_error("SUM Numerator Gamma reach INFINITE!");
            }
        }
    }

    //  panjangdna1

    for (int k=0;k<nState;k++){
        for (int l=0;l<nState;l++){
            //  discomment to update Transition
            T[k][l] = eexp(elnproduct(SumNumKsi[k][l],-SumNumGamma[k]));
        }
        for (int l=0;l<nData;l++){
            E[l][k] = eexp(elnproduct(SumNumObj[l][k],-SumNumGamma[k]));
        }
    }

}

void lHMM::getTEphi(){
    //  show PHI
    printf("PHI\n");
    for (int i=0;i<nState;i++){
        printf("%.5f ",phi[i]);
    }
    printf("\n");

    printf("Transition\n");
    for (int i=0;i<nState;i++){
        for (int j=0;j<nState;j++){
            printf("%.5f ",T[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    printf("Emmision\n");
    for (int i=0;i<nData;i++){
        for (int j=0;j<nState;j++){
            printf("%.5f ",E[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

double lHMM::getProb(){
    return gprob;
}

/**
 * You may modify the nt() for different kind of use
 * */

int lHMM::nt(char t){
	if((t == 'A')||(t == 'a')) return 0;
	else if ((t == 'C')||(t == 'c')) return 1;
	else if ((t == 'G')||(t == 'g')) return 2;
	else if ((t == 'T')||(t == 't')) return 3;
    //  else for other than ACGT
	else {return -1;};
};

/*




    Script below are the logarithmic functions





*/
double lHMM::eexp(double x){
    if (x == LOGZERO){
        return 0;
    }else{
        return exp(x);
    }
}

double lHMM::eln(double x){
    if(x==0){
        return LOGZERO;
    }else if(x>0){
        return log(x);
    }else{
        printf("Negative input error!\n");
    }
}

double lHMM::elnsum(double x, double y){//  x = eln(input 1), y= eln(input 2)
    if((x==LOGZERO)||(y==LOGZERO)){
        if(x==LOGZERO){
            return y;
        }else{
            return x;
        }
    }else{
        if(x>y){
            return x + eln(1 + exp(y-x));
        }else{
            return y + eln(1 + exp(x-y));
        }
    }
}

double lHMM::elnproduct(double x,double y){//   x = eln(input 1), y= eln(input 2)
    if ((x==LOGZERO)||(y==LOGZERO)){
        return LOGZERO;
    }else{
        return x + y;
    }
}