#ifndef __INTEGRAL_H__
#define __INTEGRAL_H__

#include <iostream>
#include "FunzioneBase.h"
#include "random.h"


using namespace std;

// Base class : generic integrator

class Integral {

 public:
  
	Integral(){;};
  Integral (double a, double b){
    checkInterval (a,b);
    m_nstep = 0;
    m_h = 0; 
    m_sum = 0;
    m_integral =0;
  }

  virtual double Integra(unsigned int nstep, const FunzioneBase& f) = 0;

 protected:

  void checkInterval (double a, double b){
		m_a = min(a,b);
		m_b = max(a,b);
		if(a>b){
			m_sign = -1;
		}else if(a<b){
			m_sign = 1;
		}else if(a == b){
			m_sign = 0;
		}
  }

  unsigned int m_nstep;
	double m_prec;
  double m_a, m_b;
  double m_sum, m_integral, m_h, m_error;
  int m_sign;
};

// derived classes

class MeanMC : public Integral {

	public:

	MeanMC(double a, double b, string filename, string fileinput) : Integral(a,b){
        m_gen.SetSeed(filename, fileinput);
		m_error = 0.;
		m_sum = 0.;
		};
	~MeanMC(){;};

	double Integra(unsigned int punti, const FunzioneBase& f) override {
		m_sum = 0.;
		for(unsigned int i = 0; i < punti; i++){
			m_sum += f.Eval(m_gen.Rannyu(m_a, m_b));
		}
		m_integral = m_sum * (m_b - m_a) / punti;
		return m_integral;
	};
    double IntegraDistribution(unsigned int punti, const FunzioneBase& f) {
		m_sum = 0.;
		for(unsigned int i = 0; i < punti; i++){
			m_sum += f.Eval(m_gen.Cumulative());
		}
		m_integral = m_sum/ punti;
		return m_integral;
	};
	double GetError(){return m_error;};

	private:
	Random m_gen;
};

#endif // __INTEGRAL_H__
