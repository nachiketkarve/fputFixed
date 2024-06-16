#include <iostream>
#include <iomanip>
#include <cmath>
#include "eigen/Eigen/Dense"
#include "eigen/unsupported/Eigen/CXX11/Tensor"
#include "eigen/unsupported/Eigen/Polynomials"
#include "eigen/Eigen/LU"
#include "fputAlphaFixed.h"
#include "eigen/Eigen/Eigenvalues"
#include <fstream>
#include <string>
#include <complex>

using namespace fputAlphaFixed;

int main()
{
    int N = 8;
    std::cout << "Enter System Size ";
    std::cin >> N;
    double E0 = 0.077;
    std::cout << "Enter Total Energy ";
    std::cin >> E0;
    double k0 = 1;
    std::cout << "Enter Seed Mode ";
    std::cin >> k0;
    double alpha = 0.01;
    std::cout << "Enter Initial alpha ";
    std::cin >> alpha;
    double deltaAlpha = 0.01;
    std::cout<<"Enter delta alpha ";
    std::cin>>deltaAlpha;
    double alphaMax = 1.0;
    std::cout << "Enter Max alpha ";
    std::cin >> alphaMax;

    std::string  FileName;

    double tmax = 1000.0;
    double deltaT = 0.01;
    std::cout<<"Enter time step ";
    std::cin >> deltaT;


    Matrix FourierComponents(N + 2, N + 2);

    initialize(N, FourierComponents);

    Vector q0(N + 2);
    Vector p0(N + 2);
    Vector Q0(N + 2);
    Vector P0(N + 2);
    Vector Z0(N - 1);

    initialState(q0, p0, Q0, P0, N, E0, k0);

    Vector Qb(N+2);
    Qb = Q0;
    double E0b = E0;
    int count = 0;

    while (alpha <= alphaMax)
    {
        Pop(Z0, Q0, k0);
        Vector Z(N - 1);
        Z = Z0;

        double period = 0;
        period = PoincareSection_rk(Z, k0, E0, alpha, tmax, deltaT);
        std::cout << period << "\n";

        Matrix QeM(N + 2, N + 2);
        Matrix PeM(N + 2, N + 2);
        PerturbationMatrix(QeM, PeM, Q0, k0, alpha);
        PoincareSectione_rk(QeM, PeM, Q0, k0, E0, alpha, tmax, deltaT);

        Matrix DF(Z0.size(), Z0.size());
        DF.setZero();

        int k = 0;

        for (int i = 0; i < N + 2; i++)
        {
            if (i != 0 & i != N + 1 & i != int(k0))
            {
                Vector Qe(N + 2);
                Vector Ze(N - 1);
                Qe.setZero();
                for (int j = 0; j < N + 2; j++)
                {
                    Qe(j) = QeM(j, i);
                }
                Pop(Ze, Qe, k0);
                for (int j = 0; j < N - 1; j++)
                {
                    DF(j, k) = Ze(j);
                }
                k = k + 1;
            }
        }
       
        Matrix Mat(Z0.size(), Z0.size());
        Mat.setZero();
        Matrix Id(Z0.size(), Z0.size());
        Id.setZero();
        for (int i = 0; i < Z0.size(); i++)
        {
            Id(i, i) = 1.0;
        }

        Mat = DF - Id;
        Vector Zn(Z0.size());
        Vector Qn(N + 2);
        Zn.setZero();

        Zn = Z0 - (Mat.inverse()) * (Z - Z0);
        
        Push(Qn, Zn, k0, E0, alpha);

        Vector Znn(Z.size());
        Znn = Zn;
        PoincareSection_rk(Znn, k0, E0, alpha, tmax, deltaT);
        Vector Qnn(N + 2);
        Vector qnn(N + 2);
        Vector pnn(N + 2);
        Push(Qnn, Znn, k0, E0, alpha);
        qnn = FourierComponents * Qnn;
        pnn.setZero();
        double lmd = Mod(Qnn - Qn) / Mod(Qn);
        std::cout << "Energy = " << TotalEnergy(qnn, pnn, alpha) << " lmd = " << lmd << "\n";
        Q0 = Qn;
        
        if(lmd > pow(10,-1) || count > 15 || std::isnan(lmd) || std::isnan(TotalEnergy(qnn, pnn, alpha))){
            count = 0;
            Q0 = Qb;
            E0 = E0b;
            alpha = alpha + deltaAlpha;
            std::cout << "alpha = " << alpha << "\n";
        }

        if (lmd < pow(10, -12) && lmd > 0)
        {
            count = 0;
            q0 = FourierComponents*Q0;
            P0.setZero();
            p0.setZero();

            Matrix Phi(2*N, 2*N);

            Matrix QeFull(N+2,2*N);
            Matrix PeFull(N+2,2*N);
            QeFull.setZero();
            PeFull.setZero();
            Phi.setZero();
            for(int i = 0; i < 2*N; i++){
                Phi(i,i) = 1;
            }
            
           for(int i = 1; i < N+1; i++){
                for(int j = 0; j < 2*N; j++){
                    QeFull(i,j) = Phi(i-1,j);
                    PeFull(i,j) = Phi(i+N-1,j);
                }
           }

           period = PoincareSectione_rk(QeFull, PeFull, Q0, k0, E0, alpha, tmax, deltaT);

           for(int i = 1; i < N+1; i++){
                for(int j = 0; j < 2*N; j++){
                    Phi(i-1,j) = QeFull(i,j);
                    Phi(i+N-1,j) = PeFull(i,j);
                }
           }

            std::cout<<period<<"\n";

            Eigen::EigenSolver<Matrix> eigensolver(Phi);
            Eigen::Vector<std::complex<double>,Eigen::Dynamic> eigVals = eigensolver.eigenvalues();
            /*
            Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> eigVects = eigensolver.eigenvectors();

            double maxEigVal = std::abs(eigVals[0]);
            int maxEigValPos = 0;
            for(int i = 0; i < 2*N; i++){
                if(std::abs(eigVals[i]) > maxEigVal){
                    maxEigVal = std::abs(eigVals[i]);
                    maxEigValPos = i;
                }
            }
            std::cout<<maxEigValPos<<"\n\n";
            Eigen::Vector<std::complex<double>,Eigen::Dynamic> maxEigVect = eigVects.col(maxEigValPos);
            std::cout<<maxEigVect<<"\n\n";
            */
            std::ofstream file;
            FileName = "Fixed-N" + std::to_string(N) + "-K" + std::to_string(int(k0)) + "-E" + std::to_string(E0) + "-A" + std::to_string(alpha) + ".csv";
            file.open(FileName);
            file<<std::setprecision(15);
            file<<"Q"<<","<<"EigValReal"<<","<<"EigValImag"<<","<<"N"<<","<<"alpha"<<","<<"k"<<","<<"E0"<<","<<"lmd"<<","<<"period"<<"\n";
            file<<Qn[0]<<","<<std::real(eigVals(0))<<","<<std::imag(eigVals(0))<<","<<N<<","<<alpha<<","<<k0<<","<<E0<<","<<lmd<<","<<period<<"\n";
            for(int i = 1; i < eigVals.size(); i++){
                if(i < Qn.size()){
                    file<<Qn(i);
                } else {
                    file<<"";
                }

                if(i < eigVals.size()){
                    file<<","<<std::real(eigVals(i))<<","<<std::imag(eigVals(i))<<"\n";
                } else {
                    file<<"\n";
                }
            }

            alpha = alpha + deltaAlpha;
            std::cout << "alpha = " << alpha << "\n";
            Qb = Qn;
            E0b = E0;

        }
        
    
        count = count + 1;
    }

    system("pause");
}