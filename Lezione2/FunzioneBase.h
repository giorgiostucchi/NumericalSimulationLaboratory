#ifndef __FunzioneBase__
#define __FunzioneBase__

#include <cmath>

class FunzioneBase {

  public:

  virtual double Eval (double) const = 0;
  virtual ~FunzioneBase() {;};

};
class Integranda : public FunzioneBase {

public:

  Integranda () {;}
  double Eval (double x) const override {return M_PI_2 * cos(M_PI_2 * x );}

};

class Integranda2 : public FunzioneBase {

public:

  Integranda2 () {;}
  double Eval (double x) const override {return M_PI_2 * cos(M_PI_2 * x ) / (-2. * (x - 1));}

};

#endif // __FunzioneBase__