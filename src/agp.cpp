#include "agp.hpp"
#include "standards.hpp"

// Calculate the AGP of a periodic state
double AGPperiodic(const State &xPeriodic, int seedMode, double deltaT, int iterations, double tmax, std::string integrator)
{
    State x = xPeriodic;
    int N = x.N;
    double period = 0.0;
    int count = 0;
    Matrix FourierComponentsSeedMode(1, N + 2);
    for (int i = 0; i < N + 2; i++)
    {
        FourierComponentsSeedMode(0, i) = sqrt(2.0 / double(N + 1)) * sin(double(i) * double(seedMode) * pi / double(N + 1));
    }

    double QSeedModeInit = (FourierComponentsSeedMode * x.q)(0, 0);
    double PSeedModeInit = (FourierComponentsSeedMode * x.p)(0, 0);

    double QSeedMode = QSeedModeInit;
    double PSeedMode = PSeedModeInit;
    double PSeedModePrev = PSeedModeInit;

    double A1 = 0;
    double A2 = 0;
    double dH = x.interactionPotential();
    double dHPrev = x.interactionPotential();

    for (int i = 0; i < int(tmax / deltaT); i++)
    {
        PSeedModePrev = PSeedMode;
        if (integrator == "rk6")
        {
            x.rkn6(deltaT);
        } else
        {
            x.rkn4(deltaT);
        }
        
        period = period + deltaT;
        dH = x.interactionPotential();
        A1 = A1 + 0.5*deltaT*(dH + dHPrev);
        A2 = A2 + 0.5*deltaT*(dH*period + dHPrev*(period-deltaT));
        dHPrev = dH;

        QSeedMode = (FourierComponentsSeedMode * x.q)(0, 0);
        PSeedMode = (FourierComponentsSeedMode * x.p)(0, 0);

        if (i > 1 && sgn(QSeedMode) == sgn(QSeedModeInit) && sgn(PSeedMode - PSeedModeInit) != sgn(PSeedModePrev - PSeedModeInit))
        {
            count = count + 1;
        }
        if (count >= iterations)
        {
            break;
        }
    }

    return A1/2.0 - A2/period;
}

// Calculate the AGP gradient of a periodic orbit starting at the time reversal symmetric point
void AGPGradientPeriodic(Vector &AGPGrad, const State &periodic, double period, int seedMode, Eigen::Vector<std::complex<double>, Eigen::Dynamic> &eigVals, Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> &eigVecs, double deltaT, double tmax, int iterations, std::string integrator)
{
    int N = periodic.N;
    State x = periodic;
    Matrix qTangent(N + 2, 2 * N);
    Matrix pTangent(N + 2, 2 * N);
    qTangent.setZero();
    pTangent.setZero();
    Matrix Phi(2*N,2*N);
    Phi.setZero();
    for (int i = 0; i < 2 * N; i++)
    {
        Phi(i, i) = 1;
    }
    for (int i = 1; i < N + 1; i++)
    {
        for (int j = 0; j < 2 * N; j++)
        {
            qTangent(i, j) = Phi(i - 1, j);
            pTangent(i, j) = Phi(i + N - 1, j);
        }
    }

    Vector angles(2*N);
    angles.setZero();
    for (int i = 0; i < 2*N; i++)
    {
        angles(i) = std::abs(std::arg(eigVals(i)));
    }

    int skip1 = 0;
    int skip2 = 1;
    if (angles(skip2) < angles(skip1))
    {
        skip1 = 1;
        skip2 = 0;
    }

    for (int i = 0; i < 2*N; i++)
    {
        if (angles(i) <= angles(skip1))
        {
            skip2 = skip1;
            skip1 = i;
        } else if (angles(i) < angles(skip2))
        {
            skip2 = i;
        }
    }

    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> G(2*N, 2*N);
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> eigVecsInv(2*N, 2*N);
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> GReal(2*N, 2*N);
    G.setZero();
    eigVecsInv.setZero();
    GReal.setZero();
    eigVecsInv = eigVecs.inverse();
    int floqSkip = 0;
    for (int i = 0; i < eigVals.size(); i++)
    {
        //if (std::abs(std::imag(eigVals(i))) > pow(10,-3) && std::abs(1.0 - std::real(eigVals(i))) > pow(10,-3))
        if(i != skip1 && i!= skip2)
        {
            G(i,i) = 0.5/(1.0 - eigVals(i));
        } else
        {
            G(i,i) = 0.5;
            floqSkip = floqSkip + 1;
        }
    }

    std::cout << "EigVals Skipped = " << floqSkip << "\n";

    GReal = eigVecs * G * eigVecsInv;
    //std::cout << eigVecs.inverse() * Phi * eigVecs << "\n";
    //std::cout << G.real() << "\n";
    //std::cout << eigVals << "\n";

    int count = 0;
    Matrix FourierComponentsSeedMode(1, N + 2);
    for (int i = 0; i < N + 2; i++)
    {
        FourierComponentsSeedMode(0, i) = sqrt(2.0 / double(N + 1)) * sin(double(i) * double(seedMode) * pi / double(N + 1));
    }

    double QSeedModeInit = (FourierComponentsSeedMode * x.q)(0, 0);
    double PSeedModeInit = (FourierComponentsSeedMode * x.p)(0, 0);

    double QSeedMode = QSeedModeInit;
    double PSeedMode = PSeedModeInit;
    double PSeedModePrev = PSeedModeInit;

    Vector integrand(2*N);
    Vector integrandPrev(2*N);
    Vector potGrad(2*N);
    integrand.setZero();
    integrandPrev.setZero();
    potGrad.setZero();
    AGPGrad.setZero();

    x.potentialGradient(potGrad);
    integrand = potGrad.transpose() * GReal.real() * Phi;

    for (int i = 0; i < int(tmax / deltaT); i++)
    {
        PSeedModePrev = PSeedMode;
        integrandPrev = integrand;
        if (integrator == "rk6")
        {
            x.rknTangent6(qTangent, pTangent, deltaT);
        } else
        {
            x.rknTangent4(qTangent, pTangent, deltaT);
        }
        
        QSeedMode = (FourierComponentsSeedMode * x.q)(0, 0);
        PSeedMode = (FourierComponentsSeedMode * x.p)(0, 0);

        for (int j = 1; j < N + 1; j++)
        {
            for (int k = 0; k < 2 * N; k++)
            {
                Phi(j - 1, k) = qTangent(j, k);
                Phi(j + N - 1, k) = pTangent(j, k);
            }
        }

        G(skip1,skip1) = 0.5 - double(i)*deltaT/period;
        G(skip2,skip2) = 0.5 - double(i)*deltaT/period;

        GReal = eigVecs * G * eigVecsInv;

        x.potentialGradient(potGrad);
        integrand = potGrad.transpose() * GReal.real() * Phi;
        AGPGrad = AGPGrad + 0.5*deltaT*(integrand + integrandPrev);

        if (i > 1 && sgn(QSeedMode) == sgn(QSeedModeInit) && sgn(PSeedMode - PSeedModeInit) != sgn(PSeedModePrev - PSeedModeInit))
        {
            count = count + 1;
        }
        if (count >= iterations)
        {
            break;
        }
    }


}