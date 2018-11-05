#include <Rcpp.h>
#include <random>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix MakeBootstrap(NumericMatrix BootstrapData, int BootstrapNumber, 
                            int BootstrapGroupSize, double BootstrapSeed) {
    int k = 0;
    int StartRow = 0;
    double SumExperiment=0, SumControl=0, pExperiment = 0, pControl = 0;
    double HControl=0, HExperiment=0, h=0, n=0;
    NumericMatrix out(BootstrapGroupSize,2);
    NumericMatrix result(BootstrapNumber,2);
    
    srand(BootstrapSeed);
    
    for(int i=0; i<BootstrapNumber; i++){
        SumControl = 0;
        SumExperiment = 0;
        pControl = 0;
        pExperiment = 0;
        
        for(int j=0; j<BootstrapGroupSize; j++){
            StartRow = rand() % BootstrapGroupSize;
            out(j,0) = BootstrapData(StartRow,0);
            out(j,1) = BootstrapData(StartRow,1);
            k+=1;
            SumControl += out(j,0);
            SumExperiment += out(j,1);
        }
        
        for(int j=0; j<BootstrapGroupSize; j++){
            if(out(j,0)!=0){
                pControl += (out(j,0)/SumControl)*log((out(j,0)/SumControl));
            }
            
            if(out(j,1)!=0){
                pExperiment += (out(j,1)/SumExperiment)*
                                log((out(j,1)/SumExperiment));
            }
        }
        
        HControl= -pControl/log(BootstrapGroupSize);
        HExperiment = -pExperiment/log(BootstrapGroupSize);
        
        if((HExperiment + HControl)>0){
            h = HExperiment/(HExperiment + HControl);
        }else{
            h = 0;
        }
        
        if((SumExperiment + SumControl)>0){
            n = SumExperiment/(SumExperiment + SumControl);
        }else{
            n = 0;
        }
        
        result(i,0) = h;
        result(i,1) = n;
    }
    return result;
}