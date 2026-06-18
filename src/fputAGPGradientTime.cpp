#include "fputFixed.hpp"
#include "nlohmann/json.hpp"
#include <random>

using json = nlohmann::json;

int main(int argc, char *argv[])
{
    int N = 8;
    double E0 = 1;
    int k0 = 1;
    double nonLin = 0.0;
    int seed = 1;

    if (argc != 6)
    {
        std::cout << argc << "\n";
        return -1;
    }

    N = int(std::atof(argv[1]));
    E0 = std::atof(argv[2]);
    k0 = std::atof(argv[3]);
    nonLin = std::atof(argv[4]);
    seed = std::atof(argv[5]);

    std::mt19937                        generator(seed);
    std::uniform_real_distribution<double>  distr(0,2.0*pi);
    
    std::ifstream dataFile("_params.json");
    json initData = json::parse(dataFile);

    std::string model = initData["model"];
    double deltaT = initData["deltaT"];
    double tmax = initData["tmax"];
    int dataPoints = initData["dataPoints"];
    std::string integrator = initData["integrator"];
    std::string saveFolder = initData["saveFolder"];
    int averages = initData["averages"];

    Matrix FourierComponents(N + 2, N + 2);
    for (int i = 0; i < N + 2; i++)
    {
        for (int j = 0; j < N + 2; j++)
        {
            FourierComponents(i, j) = sqrt(2.0 / double(N + 1)) * sin(double(i) * double(j) * pi / double(N + 1));
        }
    }

    Matrix AGPgradAvg(2*N,int(tmax/deltaT));
    AGPgradAvg.setZero();

    for (int average = 0; average < averages; average ++)
    {

        State x(N, model);
        x.nonLin = nonLin;
        
        double theta = distr(generator);

        Vector Q(N+2);
        Vector P(N+2);
        Q.setZero();
        P.setZero();
        Q(k0) = std::sqrt(2.0*E0)*std::cos(theta)/freq(k0,N);
        P(k0) = std::sqrt(2.0*E0)*std::sin(theta);
        if (x.model == "beta" && nonLin*CouplingBeta(k0,k0,k0,k0,N)*Q(k0)*Q(k0)*Q(k0)*Q(k0) != 0)
        {
            Vector coeff(5);
            coeff.setZero();
            coeff(0) = -E0;
            coeff(2) = 0.5*P(k0)*P(k0) + 0.5*freq(k0,N)*freq(k0,N)*Q(k0)*Q(k0);
            coeff(4) = 3.0/4.0*nonLin*CouplingBeta(k0,k0,k0,k0,N)*Q(k0)*Q(k0)*Q(k0)*Q(k0);
            Eigen::PolynomialSolver<double, 4> solver;
            solver.compute(coeff);
            bool realEx;
            double scale = solver.greatestRealRoot(realEx, 0.000001);
            Q(k0) = scale*Q(k0);
            P(k0) = scale*P(k0);
        }
        x.q = FourierComponents * Q;
        x.p = FourierComponents * P;
        

        //x.initializeNormalMode(E0, k0);

        State xRev = x;

        Matrix qTangent(N + 2, 2*N);
        Matrix pTangent(N + 2, 2*N);
        Matrix qTangentRev(N + 2, 2*N);
        Matrix pTangentRev(N + 2, 2*N);

        qTangent.setZero();
        pTangent.setZero();
        qTangentRev.setZero();
        pTangentRev.setZero();

        for (int i = 1; i < N+1; i++)
        {
            qTangent(i,i-1) = 1.0;
            pTangent(i,i+N-1) = 1.0;
            qTangentRev(i,i-1) = 1.0;
            pTangentRev(i,i+N-1) = 1.0;
        }

        Matrix Phi(2*N,2*N);
        Matrix PhiRev(2*N,2*N);

        Phi.setZero();
        PhiRev.setZero();

        Matrix AGPgrad(2*N, int(tmax/deltaT));

        Vector A1(2*N);
        Vector A2(2*N);
        A1.setZero();
        A2.setZero();

        AGPgrad.setZero();

        Vector potGrad(2*N);
        x.potentialGradient(potGrad);

        Vector potGradRev(2*N);
        xRev.potentialGradient(potGradRev);

        double currT = 0;

        for (int i = 1; i < int(tmax/deltaT); i++)
        {
            x.rknTangent4(qTangent, pTangent, deltaT);
            xRev.rknTangent4(qTangentRev, pTangentRev, -deltaT);
            x.potentialGradient(potGrad);
            xRev.potentialGradient(potGradRev);
            currT = currT + deltaT;

            for(int i2 = 0; i2 < 2*N; i2++)
            {
                for (int i1 = 0; i1 < N; i1 ++)
                {
                    Phi(i1,i2) = qTangent(i1+1,i2);
                    Phi(i1+N,i2) = pTangent(i1+1,i2);
                    PhiRev(i1,i2) = qTangentRev(i1+1,i2);
                    PhiRev(i1+N,i2) = pTangentRev(i1+1,i2);
                }
            }

            A1 = A1 + 0.5 * deltaT * (Phi.transpose() * potGrad) - 0.5 * deltaT * (PhiRev.transpose() * potGradRev);
            A2 = A2 + currT * 0.5 * deltaT * (Phi.transpose() * potGrad) - currT * 0.5 * deltaT * (PhiRev.transpose() * potGradRev);

            for (int j = 0; j < 2*N; j++)
            {
                AGPgrad(j,i) = A1(j) - A2(j)/currT;
            }
        }

        AGPgradAvg = AGPgradAvg + AGPgrad/averages;

    }

    int iStep = int(tmax/(deltaT*dataPoints));

    std::string FileName;

    if (model == "alpha")
    {
        FileName = saveFolder + "alphaInit-N" + std::to_string(N) + "-K" + std::to_string(int(k0)) + "-E" + std::to_string(E0) + "-A" + std::to_string(nonLin) + ".csv";
    } else if (model == "toda")
    {
        FileName = saveFolder + "todaInit-N" + std::to_string(N) + "-K" + std::to_string(int(k0)) + "-E" + std::to_string(E0) + "-A" + std::to_string(nonLin) + ".csv";
    } else if (model == "beta")
    {
        FileName = saveFolder + "betaInit-N" + std::to_string(N) + "-K" + std::to_string(int(k0)) + "-E" + std::to_string(E0) + "-B" + std::to_string(nonLin) + ".csv";
    }

    std::ofstream file;
    file.open(FileName);
    file << std::setprecision(15);
    file << "Time" << "," << "AGPGradNorm";

    for (int i = 0; i < int(tmax/deltaT); i = i + iStep)
    {
        file << "\n" << i*deltaT << "," << AGPgradAvg.col(i).norm();
    }

    file.close();

    system("pause");
    return 0;
}