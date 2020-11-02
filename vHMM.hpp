#include <iostream>
#include <vector>
#include <math.h>

#ifndef VHMM_H
#define VHMM_H
#endif

class vHMM{
private:
	const double LOGZERO = -INFINITY;
    int nState;
    std::vector<std::vector<double> > T;//   T[i][j] transition state from i to j
    std::vector<std::vector<double> > E;
    std::vector<double> phi;
    //  special
    std::vector<int> forward(std::vector<char>& _dna,int sizedna);
    std::vector<int> backward(std::vector<std::vector<double> > data, int _max, int pj);
    int argmax(std::vector<double> inp);
    int nt(char t);
    /*char revnt(int t);*/
	//	implementation
	double eexp(double x);
    double eln(double x);
    double elnsum(double x, double y);  //eln(x) + eln(y)
    double elnproduct(double x,double y);//eln(x) * eln(y)

public:
    //  constructor
    vHMM(int _nstate,std::vector<std::vector<double> > _T,std::vector<std::vector<double> > _E,std::vector<double> _phi);

    //  set parameter dan data
    void setParameter(int _nstate,std::vector<std::vector<double> > _T,std::vector<std::vector<double> > _E,std::vector<double> _phi);
    std::vector<int> viterbi(std::vector<char>& _dna, int panjangdna);
};

//  every coding goes here
vHMM::vHMM(int _nstate,std::vector<std::vector<double> > _T,std::vector<std::vector<double> > _E,std::vector<double> _phi){
    setParameter(_nstate,_T,_E,_phi);
}

void vHMM::setParameter(int _nstate,std::vector<std::vector<double> > _T,std::vector<std::vector<double> > _E,std::vector<double> _phi){
    nState = _nstate;
    T = _T;
    E = _E;
    phi = _phi;
}

std::vector<int> vHMM::viterbi(std::vector<char>& _dna, int panjangdna){
	return forward(_dna,panjangdna);
}

/*
    PRIVATE FUNCTION
*/

std::vector<int> vHMM::forward(std::vector<char>& _dna,int sizedna){
	//2d _argmax
	//source: http://www-h.eng.cam.ac.uk/help/tpl/languages/C++/vectormemory.html
	std::vector<double> temp(nState,0);
	std::vector<std::vector<double> > _argmax(sizedna,temp); 
	
	//	initial probability
	std::vector<double> curr_prob = phi;
	
	//mult
	for (int i=0; i<nState;i++){
		curr_prob[i] = elnproduct(eln(phi[i]),eln(E[nt( _dna[0])][i]));
        //std::cout<<curr_prob[i]<<" ";
	}

    //std::cout<<std::endl;
    
	//forwarding
	for(int i = 1;i<sizedna;i++){
		for(int j =0;j<nState;j++){
			double mx = -INFINITY;
			double t;
			for(int k= 0;k<nState;k++){
				t = elnproduct(curr_prob[k],elnproduct(eln(T[k][j]),eln(E[nt(_dna[i])][j])));
				if (t>mx) {
					_argmax[i][j]= k;
					mx = t;
				}
			}
			curr_prob[j] = mx;
		}
	}

	//	ubah sesuai nState
	return backward(_argmax,argmax(curr_prob),sizedna-1);	
}

std::vector<int> vHMM::backward(std::vector<std::vector<double> > data, int _max, int pj){
    std::vector<int> _s(pj+1,0);
	for (int i=pj;i>1;i--){
        _s[i] = data[i][_max];
		_max = (int)(data[i][_max]);
    }
	return _s;
}

int vHMM::nt(char t){
	if((t == 'A')||(t == 'a')) return 0;
	else if ((t == 'C')||(t == 'c')) return 1;
	else if ((t == 'G')||(t == 'g')) return 2;
	else if ((t == 'T')||(t == 't')) return 3;
	else return -1;
};

/*
char vHMM::revnt(int t){
	if (t==0) return 'A';
	else if(t==1) return 'C';
	else if(t==2) return 'G';
	else if(t==3) return 'T';
	else return 'N';
}
*/

int vHMM::argmax(std::vector<double> inp){
	int arrSize = inp.size();
	double max = -3.402823e+38;
	int pos;
	for (int i=0;i<arrSize;i++){
		if (inp[i]>max){
			max = inp[i];
			pos = i;
		}
	}
	return pos;
};

/*




    Logarithmic transfer function





*/
double vHMM::eexp(double x){
    if (x == LOGZERO){
        return 0;
    }else{
        return exp(x);
    }
}

double vHMM::eln(double x){
    if(x==0){
        return LOGZERO;
    }else if(x>0){
        return log(x);
    }else{
        printf("Negative input error!\n");
    }
}

double vHMM::elnsum(double x, double y){//  x = eln(input 1), y= eln(input 2)
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

double vHMM::elnproduct(double x,double y){//   x = eln(input 1), y= eln(input 2)
    if ((x==LOGZERO)||(y==LOGZERO)){
        return LOGZERO;
    }else{
        return x + y;
    }
}
