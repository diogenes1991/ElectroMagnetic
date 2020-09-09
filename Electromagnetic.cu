#include <iostream>
#include <sched.h>
#include <pthread.h>
#include <thread>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <random> 

#include"cloak.h"

#define GRANULARITY 13
#define NDIMENSIONS 8

/// Binary Preprocessor Calculator

#define XOR_0_0 0
#define XOR_0_1 1
#define XOR_1_0 1
#define XOR_1_1 0
#define XOR(X,Y) XOR_##X##_##Y

#define AND_0_0 0
#define AND_0_1 0
#define AND_1_0 0
#define AND_1_1 1
#define AND(X,Y) AND_##X##_##Y

#define REMOVE_PARS(...) EXPAND __VA_ARGS__
#define FIRST(X,...) X
#define TRIM_FIRST(X,...) __VA_ARGS__    

#define BIT_ADD(X,Y,...) (XOR(X,Y)IF(AND(X,Y))(COMMA()AND(X,Y),))

#define MUL(X,Y)*Y
#define POW(X,N) 1 REPEAT(N,MUL,X)

#define BUNCH_SIZE 1<<20
#define CREATE_POTENTIAL(N,J) T A ##J## _ ##N [int(POW(GRANULARITY,NDIMENSIONS-1))];
#define CREATE_COPIES_INDIRECT(J,M) REPEAT(M,CREATE_POTENTIAL,J)
#define CREATE_COPIES(J,M) EVAL(REPEAT(J,CREATE_COPIES_INDIRECT,M))

template <class T>
class ElectroMagneticField{
    public:
        
    T Potential[NDIMENSIONS][GRANULARITY][GRANULARITY];
    T DPotential[NDIMENSIONS][GRANULARITY][GRANULARITY];
            
    ElectroMagneticField();
    
    // The class has several types of grid topologies 
    // the default for now is that the (N-1)-Space part 
    // is a (N-1)-Genus surface. Essentially the grid 
    // has periodic boundary conditions.
    
    void Update_Potential(int topology = 0);
    void Show_Potential(int Component);
    void Show_DPotential(int Component);
    void Show_Field(int mu, int nu);
    T Compute_Field(int mu, int nu, int x, int y);
    T DmuAnu(int mu, int nu, int x, int y);
    T Energy();
    
    
};

template <class T>
ElectroMagneticField<T>::ElectroMagneticField(){
    for (int k=0;k<NDIMENSIONS;k++){
            for(int i=0;i<GRANULARITY;i++){
                for(int j=0;j<GRANULARITY;j++){
                    Potential[k][i][j] = 0.;
                    DPotential[k][i][j] = 0.;
                }
            }
    }
}

template <class T>
void ElectroMagneticField<T>::Update_Potential(int Top){
    
    
    /// Periodic Boundary
    if(Top==0){
        for (int k=0;k<NDIMENSIONS;k++){
            for(int i=0;i<GRANULARITY;i++){
                for(int j=0;j<GRANULARITY;j++){
                
                    /// First we update the Value of the Potentials
                    Potential[k][i][j] += DPotential[k][i][j];
                
                    T val = 0;
                    /// First Component
                    val += (Potential[k][(i+2)%GRANULARITY][j] + Potential[k][(GRANULARITY+i-2)%GRANULARITY][j] - 2*Potential[k][i][j])/4;
                    val += (Potential[k][i][(j+2)%GRANULARITY] + Potential[k][i][(GRANULARITY+j-2)%GRANULARITY] - 2*Potential[k][i][j])/4;
                    DPotential[k][i][j] = val;
                }
            }
        
        }
    }
    
}

template <class T>
void ElectroMagneticField<T>::Show_Potential(int Component){
    
    std::cout << Component << "-component of the potential" <<std::endl;
    for(int i=0;i<GRANULARITY;i++){
            std::cout << "| ";
            for(int j=0;j<GRANULARITY;j++){
                    std::cout << Potential[Component][i][j] << " ";
            }
            std::cout << "|" <<std::endl;
    }
    
    
}

template <class T>
void ElectroMagneticField<T>::Show_DPotential(int Component){
    
    std::cout << Component << " component of the derivative of the potential" <<std::endl;
    for(int i=0;i<GRANULARITY;i++){
            std::cout << "| ";
            for(int j=0;j<GRANULARITY;j++){
                    std::cout << DPotential[Component][i][j] << " ";
            }
            std::cout << "|" <<std::endl;
    }
    
    
}

template <class T>
void ElectroMagneticField<T>::Show_Field(int mu, int nu){
    std::cout << "("<<mu<<","<<nu<<")-component of the field" <<std::endl;
    for(int i=0;i<GRANULARITY;i++){
            std::cout << "| ";
            for(int j=0;j<GRANULARITY;j++){
                std::cout << Compute_Field(mu,nu,i,j) << " ";
           }
            std::cout << "|" <<std::endl;
    }
}

template <class T>
T ElectroMagneticField<T>::DmuAnu(int mu, int nu, int x, int y){
    T rval = 0;
    switch (mu){
        case 0:
            rval = DPotential[nu][x][y];
            break;
        case 1:
            rval = (Potential[nu][(x+1)%GRANULARITY][y]+Potential[nu][(GRANULARITY+x-1)%GRANULARITY][y])/2;
            break;
        case 2:
            rval = (Potential[nu][x][(y+1)%GRANULARITY]+Potential[nu][x][(GRANULARITY+y-1)%GRANULARITY])/2;
            break;
    }
    return rval;
}

template <class T>
T ElectroMagneticField<T>::Compute_Field(int mu, int nu, int x, int y){
    T rval = DmuAnu(mu,nu,x,y)-DmuAnu(nu,mu,x,y);
    return rval;
}

template <class T>
T ElectroMagneticField<T>::Energy(){
    T U = 0;
    for(int i=0;i<GRANULARITY;i++){
            for(int j=0;j<GRANULARITY;j++){
                for(int k = 0;k<NDIMENSIONS;k++){
                    for(int m = k+1;m<NDIMENSIONS;m++){
                        T aux = Compute_Field(k,m,i,j);
                        U += aux*aux/2;
                        }
                    }
           }
    }
    return U;
}

void VECT_TO_FLAT(int* vect_x, int* flat_x){
        *flat_x = 0;
        for(int i=0;i<NDIMENSIONS-1;i++){
            *flat_x += vect_x[i]*pow(GRANULARITY,i);
        }
}

void FLAT_TO_VECT(int* vect_x, int* flat_x){
    int aux = *flat_x;
    for(int i=0;i<NDIMENSIONS-1;i++){
            vect_x[i] = aux%GRANULARITY;
            aux -= aux%GRANULARITY;
            aux /= GRANULARITY;
    }
}

int main(void){
    
    typedef double T;
    CREATE_COPIES(4,1)
     
    for (int i=0;i<pow(GRANULARITY,NDIMENSIONS-1);i++){
        A0_0[i] = i;
        
    }
//     
    std::cout<<"| ";
    for (int i=0;i<pow(GRANULARITY,NDIMENSIONS-1);i++){
        std::cout << A0_0[i] << " ";
        if(i%GRANULARITY == GRANULARITY-1)
            if(i==pow(GRANULARITY,NDIMENSIONS-1)-1)std::cout << "|\n";
            else std::cout << "|\n| ";
        
    }
//     
//     dA0_0[x] ~ SUM(2*NDIM Terms)/(2*NDIM) - A0_0[x] 
//     
//     A0_0[x][y][z] ~ A0_0[x+2][y][z] + A0_0[x-2][y][z] + ...
//     
//     100 ~ x[2] = {15,4}
//     {15,4} ~ 100 
//     
//     INV_MAP(100) = {15,4}
//     A0_0(MAP({15+2,4}) + A0_0({15-2,4})
//     
    
//     ElectroMagneticField<double> E;
//     E.Potential[0][10][10] = 100;
//     for(int i=0;i<1000;i++){
//         E.Update_Potential();
//         if(i%10==0){
//             std::cout<<"U = ("<<i/10<<"%) = "<< E.Energy()<<std::endl;
//         }
//     }
//     
//     E.Show_Potential(0);
//     E.Show_Field(0,1);
//     E.Show_Field(0,2);
//     E.Show_Field(1,2);
    
    int h[NDIMENSIONS-1];
    int g = 170;
    #if 3>110
    CREATE_COPIES(3,5)
    #endif
    
    
    FLAT_TO_VECT(h,&g);
    VECT_TO_FLAT(h,&g);
    std::cout<<"("<<h[0]<<","<<h[1]<<","<<h[2]<<")"<<std::endl;
    std::cout<<"("<<g<<")"<<std::endl;
    
}
