#ifndef DIFFUSION1DP_H
#define DIFFUSION1DP_H

#include "imc_dfm.h"
#include "tdma.h"

#include "functor1d.h"
#include "boundary.h"


// Classe para coleta de dados dos problemas
class DiffusionData
{
    public:
        Functor1D *k;
        tFloat H;
        tFloat Tinf;
        tFloat Tb;
        tFloat Ab;
        tFloat P;
        tFloat mi;
        tFloat C;
        Functor1D *heatSource;
        tInteger nmax;

        DiffusionData(){
            k = new PolynomialConstant(0.0);
            H = 0.0;
            Tinf = 0.0;
            Tb = 0.0;
            Ab = 0.0;
            P = 0.0;
            mi = 0.0;
            C = 0.0;
            heatSource = new PolynomialConstant(0.0);
            nmax = 50;
        }
};


// Tipos de problemas
enum DiffusionProblem{
    ParedePlana,
    Aleta,
    QML
};

// Classe para resolução dos problemas de Difusão 1D Permanente
class Diffusion1Dp
{
    private:
        Boundary1D boundaryLeft, boundaryRight;
        DiffusionData data;
        DiffusionProblem problemType;
        tFloat l, h;
        tInteger n, nmax;
        tFloat *L;

    public:
        TDMA equationsSystem;
        tFloat  averageValue, maxValue, heatFlowLeft, heatFlowRight;

        Diffusion1Dp(){}

        Diffusion1Dp(DiffusionProblem type, tInteger nodes, Boundary1D left, Boundary1D right, DiffusionData data);

        void plotSolution(Functor1D &analyticalSolution);
        void plotInteractionLog(void);
        void printSolution(Functor1D &analyticalSolution);
        void printSecondaryResults(tFloat AS_average = 0.0q, tFloat AS_left = 0.0q, tFloat AS_right = 0.0q,
                                   tFloat AS_left_1 = 0.0q, tFloat AS_right_1 = 0.0q);

        tFloat get_h(){return h;}

        void solver(void);

};



#endif // DIFFUSION1DP_H
