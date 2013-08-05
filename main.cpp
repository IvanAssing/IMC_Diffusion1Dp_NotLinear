#include <iostream>

#include "tdma.h"
#include "functor1d.h"
#include "diffusion1dp.h"

int main()
{
    // >[QUESTÃO 3.2] - DIFUSÃO DE CALOR NÃO LINEAR EM PAREDE PLANA COM TEMPERATURA PRESCRITA
    tInteger n = 5;
    tFloat s0 = -0.5q;
    tFloat s1 = -1.5q;
    tFloat s2 = -1.0q;
    DiffusionData data;
    data.k = new Exponential(1.0q, 1.0q); // 1.0*exp(1.0*T)
    data.heatSource = new PolynomialConstant(0.0q);

    SpecialFunctorTA2 AS_Temperature; // solução analítica

    data.nmax = 50;

    tFloat AS_tm = 0.5q - s0/12.0q - s1/24.0q - s2/40.0q; // temperatura média
    tFloat AS_q0 = -1.0q + s0/2.0q + s1/6.0q + s2/12.0q; // fluxo de calor em x=0
    tFloat AS_ql = -(1.0q + s0/2.0q + s1/3.0q + s2/4.0q); // fluxo de calor em x=1

    Boundary1D left(0.0q, Dirichlet, 0.0q); // x=0, Condição de Contorno de Dirichlet, T(0) = 0
    Boundary1D right(1.0q, Dirichlet, 1.0q); // x=1, Condição de Contorno de Dirichlet, T(1) = 1

    Diffusion1Dp mesh(ParedePlana, n, left, right, data);

    mesh.solver();

    mesh.printSolution(AS_Temperature);
    mesh.plotInteractionLog();

    //mesh.printSecondaryResults(AS_tm, AS_q0, AS_ql, 0.0q, 1.0q);

    mesh.plotSolution(AS_Temperature);


}

