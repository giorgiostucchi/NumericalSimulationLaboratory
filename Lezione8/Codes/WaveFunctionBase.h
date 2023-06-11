#ifndef __WaveFunctionBase__
#define __WaveFunctionBase__

class WaveFunctionBase { //classe astratta funzione base 

  public:

  virtual double norm (double) const = 0;//metodo virtuale di valutazione della f in un punto
  virtual double hamiltonian (double) const = 0;//metodo virtuale di valutazione della f in un punto
  virtual ~WaveFunctionBase() {;};//distruttore di tipo virtual, i distruttori delle figlie andranno specializzati

};

class WaveFunction : public WaveFunctionBase {

public:

  WaveFunction() {mu=0.; sigma=0.;}
  WaveFunction (double a, double b) {mu = a; sigma=b;}
  ~WaveFunction() {;}
  double norm (double x) const override {return pow(exp(-(x-mu)*(x-mu)/(2.*sigma*sigma)) + exp(-(x+mu)*(x+mu)/(2.*sigma*sigma)),2);}//squared norm of the wavefunction
  double hamiltonian(double x) const override { 
    double sum1 = ((x*x + mu*mu)/(sigma*sigma)-1.)/(2.*sigma*sigma);
    double sum2 = (x*mu)*(exp(-(x+mu)*(x+mu)/(2.*sigma*sigma))-exp(-(x-mu)*(x-mu)/(2.*sigma*sigma)));
    double sum3 = (sigma*sigma*sigma*sigma*(exp(-(x+mu)*(x+mu)/(2.*sigma*sigma))+exp(-(x-mu)*(x-mu)/(2.*sigma*sigma))));
    return -sum1 - sum2/sum3 + x*x*x*x - (5.)/(2.)*x*x; 
  }
  void Setmu (double a) {mu=a;}
  double Getmu () const {return mu;}
  void Setsigma (double b) {sigma = b;}
  double Getsigma () const {return sigma;}

private:

  double mu, sigma;

};

#endif 