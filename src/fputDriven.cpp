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
    double Tstart = 100.0;
    double Tmax = 100.0;
    double lmd = 0.1;

    if (argc != 10)
    {
        std::cout << argc << "\n";
        return -1;
    }

    N = int(std::atof(argv[1]));
    E0 = std::atof(argv[2]);
    k0 = std::atof(argv[3]);
    double nonLinNum = std::atof(argv[4]);
    double nonLinDen = std::atof(argv[5]);
    Tstart = std::atof(argv[6]);
    Tmax = std::atof(argv[7]);
    lmd = std::atof(argv[8]);
    seed = int(std::atof(argv[9]));
    nonLin = nonLinNum/nonLinDen;

    std::mt19937                        generator(seed);
    std::uniform_real_distribution<double>  distr(0,2.0*pi);
    
    std::ifstream dataFile("_params.json");
    json initData = json::parse(dataFile);

    double deltaT = initData["deltaT"];
    int dataPoints = initData["dataPoints"];
    std::string saveFolder = initData["saveFolder"];
    int averages = initData["averages"];

    std::string fileName = saveFolder + "betaDriven-N" + std::to_string(N) + "-K" + std::to_string(int(k0)) + "-E" + std::to_string(E0) + "-B" + std::to_string(nonLin) + "-T" + std::to_string(Tmax) + "-L" + std::to_string(lmd) + ".csv";

    Matrix FourierComponents(N + 2, N + 2);
    for (int i = 0; i < N + 2; i++)
    {
        for (int j = 0; j < N + 2; j++)
        {
            FourierComponents(i, j) = sqrt(2.0 / double(N + 1)) * sin(double(i) * double(j) * pi / double(N + 1));
        }
    }

    Vector Eis(averages);
    Vector Efs(averages);
    Eis.setZero();
    Efs.setZero();

    double mu = 2.0*pi/Tmax;
    deltaT = Tmax/int(Tmax/deltaT);

    #pragma omp parallel for schedule(dynamic)
    for (int average = 0; average < averages; average++)
    {
        State x(N, "beta");
        x.nonLin = nonLin;

        double theta = distr(generator);

        Vector Q(N+2);
        Vector P(N+2);
        Q.setZero();
        P.setZero();
        Q(k0) = std::sqrt(2.0*E0)*std::cos(theta)/freq(k0,N);
        P(k0) = std::sqrt(2.0*E0)*std::sin(theta);
        if (nonLin*CouplingBeta(k0,k0,k0,k0,N)*Q(k0)*Q(k0)*Q(k0)*Q(k0) != 0)
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

        x.Evolve(Tstart, deltaT, "rk4");
        double Ei = x.totalEnergy();
        x.EvolveDriven(Tmax, deltaT, lmd, mu);
        double Ef = x.totalEnergy();
        Eis(average) = Ei;
        Efs(average) = Ef;
    }

    std::ofstream file;
    file.open(fileName);
    file << std::setprecision(15);
    file << "Ei,Ef,deltaE";
    
    for(int average = 0; average < averages; average++)
    {
        file << "\n" << Eis(average) << "," << Efs(average) << "," << Efs(average) - Eis(average);
    }
    file.close();

}