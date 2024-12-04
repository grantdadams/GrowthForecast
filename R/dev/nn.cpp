
// http://www.nielsensweb.org/fims2022/nn.cpp

#include <TMB.hpp>

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
struct NN{
  matrix<Type> W1;
  matrix<Type> W2;
  NN(matrix<Type> w1, matrix<Type> w2){W1=w1; W2=w2;};

  template<class T>
  T phi(T x){
    return T(2)/M_PI*atan(x);
  }

  template<class T>
  matrix<T> phi(matrix<T> x){
    matrix<T> ret(x.rows(),x.cols());
    for(int i=0; i<x.rows(); ++i){
      for(int j=0; j<x.cols(); ++j){
        ret(i,j)=phi(x(i,j));
      }
    }
    return ret;
  }

  template<class T>
  T operator()(vector<T> x){
    matrix<T> X(1,x.size());
    matrix<T> ww1=W1.template cast<T>();
    matrix<T> ww2=W2.template cast<T>();
    X.row(0)=x;
    matrix<T> H=X*ww1;
    vector<T> pred=vector<T>(phi(H)*ww2);
    return pred(0);
  }
};

template<class Type>
Type objective_function<Type>::operator()()
{
  DATA_INTEGER(noNodes);
  DATA_MATRIX(X);
  DATA_VECTOR(Y);

  int noInput=X.cols();
  int nobs = Y.size();

  PARAMETER_MATRIX(W1);
  PARAMETER_MATRIX(W2);
  vector<Type> pred(nobs),dx(nobs),dy(nobs);
  pred.setZero();
  Type tmp, delta;
  Type tiny = 1.0e-6;
  Type nll = 0;
  NN<Type> nn(W1,W2);
  for(int i=0; i<nobs; ++i){
    pred(i)=nn(vector<Type>(X.row(i)));
    if(!isNA(Y(i))){
      delta=Y(i)-pred(i);
      nll += delta*delta;
    }
  }

  nll += tiny*(W1.array()*W1.array()).sum();
  nll += tiny*(W2.array()*W2.array()).sum();

  REPORT(pred);
  return nll;
}
