#include <iostream>
#include <sched.h>
#include <pthread.h>
#include <thread>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <random> 

#define GRANULARITY 10
#define NDIMENSIONS 3

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
    std::cout << "("<<mu<<","<<nu<<")-component of the derivative of the potential" <<std::endl;
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
            rval = (Potential[nu][x+1][y]+Potential[nu][x-1][y])/2;
            break;
        case 2:
            rval = (Potential[nu][x][y+1]+Potential[nu][x][y-1])/2;
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

int main(void){

    std::cout.precision(4);
    ElectroMagneticField<double> E;
    E.Potential[0][5][0] = 10;
//     E.Potential[0][6][1] = 13;
//     E.Potential[0][5][2] = 11;
//     E.Potential[0][6][3] = 9;
    E.Show_Potential(0);
    std::cout << "Energy = "<<E.Energy()<<std::endl;
    for(int i=0;i<50;i++){
    E.Update_Potential();
//     E.Show_Potential(0);
    }
    E.Show_Potential(0);
    E.Show_Potential(1);
    E.Show_Potential(2);
    std::cout << "Energy = "<<E.Energy()<<std::endl;
    
    
}
