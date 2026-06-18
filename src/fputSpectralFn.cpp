#include "fputFixed.hpp"
#include "nlohmann/json.hpp"
#include <random>
#include <unsupported/Eigen/FFT>

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

    if (argc != 9)
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
    seed = int(std::atof(argv[8]));
    nonLin = nonLinNum/nonLinDen;

    std::mt19937                        generator(seed);
    std::uniform_real_distribution<double>  distr(0,2.0*pi);
    
    std::ifstream dataFile("_params.json");
    json initData = json::parse(dataFile);

    double deltaT = initData["deltaT"];
    int dataPoints = initData["dataPoints"];
    std::string saveFolder = initData["saveFolder"];
    int averages = initData["averages"];
    std::string model = initData["model"];
    std::string integrator = initData["integrator"];

    Matrix FourierComponents(N + 2, N + 2);
    for (int i = 0; i < N + 2; i++)
    {
        for (int j = 0; j < N + 2; j++)
        {
            FourierComponents(i, j) = sqrt(2.0 / double(N + 1)) * sin(double(i) * double(j) * pi / double(N + 1));
        }
    }

    std::string FileName;

    if (model == "alpha")
    {
        FileName = saveFolder + "alphaSpecFn-N" + std::to_string(N) + "-K" + std::to_string(int(k0)) + "-E" + std::to_string(E0) + "-A" + std::to_string(nonLin) + ".csv";
    } else if (model == "toda")
    {
        FileName = saveFolder + "todaSpecFn-N" + std::to_string(N) + "-K" + std::to_string(int(k0)) + "-E" + std::to_string(E0) + "-A" + std::to_string(nonLin) + ".csv";
    } else if (model == "beta")
    {
        FileName = saveFolder + "betaSpecFn-N" + std::to_string(N) + "-K" + std::to_string(int(k0)) + "-E" + std::to_string(E0) + "-B" + std::to_string(nonLin) + ".csv";
    }

    std::ofstream file;
    file.open(FileName);
    file << std::setprecision(15);
    file << "Frequency" << "," << "Spectral Function";

    Vector dH(int(Tmax/deltaT));
    Eigen::VectorXcd dhFT(int(Tmax/deltaT));
    Vector specFn(int(Tmax/deltaT));
    Vector t(int(Tmax/deltaT));

    dH.setZero();
    dhFT.setZero();
    specFn.setZero();
    t.setZero();

    for (int average = 0; average < averages; average++)
    {
        dhFT.setZero();
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

        x.Evolve(Tstart, deltaT, integrator);

        double currT = 0;

        for (int i = 0; i < int(Tmax/deltaT); i++)
        {
            dH(i) = x.interactionPotential();
            t(i) = currT;
            if (integrator == "leapFrog")
            {
                x.leapFrog(deltaT);
            } else if (integrator == "rk4")
            {
                x.rkn4(deltaT);
            } else if (integrator == "rk6")
            {
                x.rkn6(deltaT);
            }
            currT = currT + deltaT;
            if (x.model == "beta")
            {
                Vector Q(N+2);
                Vector P(N+2);
                Q = FourierComponents * x.q;
                P = FourierComponents * x.p;

                for(int k = k0%2+1; k < N+2; k=k+2)
                {
                    Q(k) = 0;
                    P(k) = 0;
                }

                x.q = FourierComponents * Q;
                x.p = FourierComponents * P;
            }
        }

        dH = dH.array() - dH.mean();

        Eigen::FFT<double> fft;
        fft.fwd(dhFT, dH);

        specFn = specFn + dhFT.cwiseAbs2()/averages;
    }

    Vector freqs(int(Tmax/deltaT));
    for (int i = 0; i < int(Tmax/deltaT); i++)
    {
        freqs(i) = 2.0*pi*double(i)/Tmax;
    }

    for (int i = 0; i < int(Tmax/(2.0*deltaT)); i = i + 1)
    {
        file << "\n" << freqs(i) << "," << specFn(i);
    }

    file.close();

    // return 1;
}