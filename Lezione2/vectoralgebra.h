#pragma once

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <algorithm>  // se se vogliono usare algoritmi STL

// ===============================================================================
// Nel namespace myVectAlgebra potremmo ridefinire gli operatori matematici che coinvolgono
// vettori e scalari
// ===============================================================================

// somma di due vettori : somma componente per componente

template <typename T> std::vector<T> operator+(const std::vector<T> &a, const std::vector<T> &b) {
  
  assert(a.size() == b.size());  
  std::vector<T> result(a.size());
  
  for (int i = 0; i < static_cast<int>(a.size()); i++) result[i] = a[i] + b[i];    
  
  // Alternativamente si puo' usare l'algoritmo transform della STL
  //
  //    std::transform(a.begin(), a.end(), b.begin(), result.begin() , std::plus<T>()); 
  
  return result;
  
}

// differenza di due vettori componente per componente
// [ preferisco re-implementarlo perche' si fanno meno operazioni rispetto result = a + (-1.*b) ]

template <typename T> std::vector<T> operator-(const std::vector<T> &a, const std::vector<T> &b) {
  
  assert(a.size() == b.size());  
  std::vector<T> result(a.size());
  
  for (int i = 0; i < static_cast<int>(a.size()); i++) result[i] = a[i] - b[i];    
  
  // Alternativamente si puo' usare l'algoritmo transform della STL
  //
  //    std::transform(a.begin(), a.end(), b.begin(), result.begin() , std::minus<T>()); 
  
  return result;
  
}

// ===============================================================================  
// prodotto scalare tra due vettori 
// ===============================================================================

template <typename T> T operator*(const std::vector<T> &a, const std::vector<T> &b) {
  
  assert(a.size() == b.size());  
  T sum = 0 ;
  for (int i = 0; i < static_cast<int>(a.size()); i++) sum += a[i] * b[i];  
  
  // Alternativamente si puo' usare l'algoritmo inner_product della STL
  //
  // sum = std::inner_product(std::begin(a), std::end(a), std::begin(b), 0.);
  
  return sum;
  
}

// ===============================================================================
// prodotto tra uno scalare e un vettore
// ===============================================================================

template <typename T> std::vector<T> operator*( T c , const std::vector<T> &a) {
  
  std::vector<T> result(a.size());
  
  for (int i = 0; i < static_cast<int>(a.size()); i++) result[i] = c * a[i];  
  
  // Alternativamente si puo' usare l'algoritmo inner product
  //
  //     std::transform(std::begin(a), std::end(a), std::begin(result), [&c](T x){ return x * c; } );
  
  return result;
  
}

// ===============================================================================
// prodotto tra un vettore e uno scalare 
// ===============================================================================

template <typename T> std::vector<T> operator*( const std::vector<T> &a , T c) {
  
  std::vector<T> result(a.size());
  for (int i = 0; i < static_cast<int>(a.size()); i++) result[i] = c * a[i];  
  
  // oppure il ciclo for puo' essere sostituito da ( ~ stesso numero di operazioni con il
  // move constructor del vector altrimenti sarebbe meno efficiente )
  // 
  // result = c * a ;
  
  // Alternativamente si puo' usare l'algoritmo transform della STL con una lambda function
  //    
  //    std::transform(std::begin(a), std::end(a), std::begin(result), [&c](T x){ return x * c; } );
  
  return result;
  
}

// ===============================================================================
// divisione tra un vettore e uno scalare 
// ===============================================================================

template <typename T> std::vector<T> operator/( const std::vector<T> &a , T c) {
  
  std::vector<T> result(a.size());
  for (int i = 0; i < static_cast<int>(a.size()); i++) result[i] = a[i] / c ;  
  
  // oppure il ciclo for puo' essere sostituito da
  
  // double fact = 1./c
  // result = a * fact ;
  
  // Alternativamente si puo' usare l'algoritmo transform della STL con una lambda function
  //    
  //    std::transform(std::begin(a), std::end(a), std::begin(result), [&c](T x){ return x / c; } );
  
  return result;
  
}

// ===============================================================================
// elevazione delle componenti di un vettore ad un esponente
// ===============================================================================

template <typename T> std::vector<T> operator^( const std::vector<T> &a , T c) {
  
  std::vector<T> result(a.size());
  for (int i = 0; i < static_cast<int>(a.size()); i++) result[i] = pow(a[i], c) ;  
  
  // oppure il ciclo for puo' essere sostituito da
  
  // double fact = 1./c
  // result = a * fact ;
  
  // Alternativamente si puo' usare l'algoritmo transform della STL con una lambda function
  //    
  //    std::transform(std::begin(a), std::end(a), std::begin(result), [&c](T x){ return x / c; } );
  
  return result;
  
}

// ===============================================================================
// somma ad a un vettore b e il risultato viene messo in a
// ===============================================================================

template <typename T> std::vector<T>& operator+=(std::vector<T>& a, const std::vector<T>& b) {
  
  assert(a.size() == b.size());  
  
  for (int i = 0; i < static_cast<int>(a.size()); i++) a[i] += b[i];    
  
  // Alternativamente si puo' usare l'algoritmo transform della STL
  //
  //    std::transform(a.begin(), a.end(), b.begin(), a.begin() , std::plus<T>()); 
  
  return a;
}


// ===============================================================================
// sottrai ad a un vettore b e il risultato viene messo in a
// ===============================================================================

template <typename T> std::vector<T>& operator-=(std::vector<T>& a, const std::vector<T>& b) {
  
  assert(a.size() == b.size());  
  
  for (int i = 0; i < static_cast<int>(a.size()); i++) a[i] -= b[i];    
  
  // Alternativamente si puo' usare l'algoritmo transform della STL
  //
  //    std::transform(a.begin(), a.end(), b.begin(), a.begin() , std::minus<T>()); 
  
  return a;
}

// ===============================================================================
// Possiamo usare il namespace myVectUtils per funzioni generiche che agiscono sui vettori
// ===============================================================================

namespace myVectUtils {
  
  // metodo comodo per stampare il vettore
  
  template <typename T> void Print( const std::vector<T> &v ) {
    
    std::cout << "Printing vector" << std::endl;
    for ( auto it = v.begin() ; it != v.end() ; it++ ) std::cout << *it << " " << std::endl;
    //std::cout << std::endl;
    std::cout << "End of printing vector" << std::endl;
    
  };

}

template <typename T> double mod(const std::vector<T> &a) {
  
  double sum = 0;
  for (int i = 0; i < static_cast<int>(a.size()); i++) sum+= pow(a[i], 2);    
  return sqrt(sum);
  
}

template <typename T> T Media( const std::vector<T> & v){
  T media = 0;
  for (unsigned int i=0; i < v.size(); i++){
    media += v[i];
  }
  media = media / double (v.size());
  return media;
};

/*
template <typename T> class linearVector : public std::vector<T> {
  
public :

  linearVector<T> ( ) : std::vector<T>(){} ;
  linearVector<T> ( int size ) : std::vector<T>(size){} ;  

  linearVector<T> operator+( const linearVector<T> &a ) {
    
    assert(a.size() == this->size());
    
    linearVector<T> result(a.size()) ;
    for (int i = 0; i < static_cast<int>(a.size()); i++) result[i] = a[i] + (*this)[i];
    
    return result;
    
  }
  
};
*/
