#include <iostream>
#include <iomanip>
#include <cmath>
#include "eigen/Eigen/Dense"
#include "eigen/unsupported/Eigen/CXX11/Tensor"
#include "eigen/unsupported/Eigen/Polynomials"
#include "eigen/Eigen/LU"

namespace fputAlphaFixed
{

    const double pi = 3.14159265358979323846;

    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Vector<double, Eigen::Dynamic> Vector;

    // rkn constants
    const double g1 = 0.0, g2 = 0.205177662, g3 = 0.608198943, g4 = 0.487278067, g5 = 1.0;
    const double b1 = 0.0617588581, b2 = 0.338978026, b3 = 0.614791307, b4 = -0.140548014, b5 = 0.125019823;
    const double a21 = b1 * (g2 - g1);
    const double a31 = b1 * (g3 - g1), a32 = b2 * (g3 - g2);
    const double a41 = b1 * (g4 - g1), a42 = b2 * (g4 - g2), a43 = b3 * (g4 - g3);
    const double a51 = b1 * (g5 - g1), a52 = b2 * (g5 - g2), a53 = b3 * (g5 - g3), a54 = b4 * (g5 - g4);
    const double B1 = b1 * (1.0 - g1), B2 = b2 * (1.0 - g2), B3 = b3 * (1.0 - g3), B4 = b4 * (1.0 - g4), B5 = b5 * (1.0 - g5);

    //--------------------------------------------------------------------------------------------------------------------------

    // returns the sign of the input
    template <typename T>
    int sgn(T val)
    {
        return (T(0) < val) - (val < T(0));
    }

    // returns the mode frequency
    double freq(int mode, int N)
    {
        return 2.0 * sin(pi * double(mode) / (2.0 * double(N + 1)));
    }

    //-------------------------------------------------------------------------------------------------------------------------

    double Coupling(int i, int j, int k, int N)
    {

        double C = 0;
        if (i + j + k == 0)
        {
            C = 1.0;
        }
        if (i + j - k == 0)
        {
            C = 1.0;
        }
        if (i - j + k == 0)
        {
            C = 1.0;
        }
        if (i - j - k == 0)
        {
            C = 1.0;
        }

        if (i + j + k == 2 * (N + 1))
        {
            C = -1.0;
        }
        if (i + j - k == 2 * (N + 1))
        {
            C = -1.0;
        }
        if (i - j + k == 2 * (N + 1))
        {
            C = -1.0;
        }
        if (i - j - k == 2 * (N + 1))
        {
            C = -1.0;
        }
        C = C * freq(i, N) * freq(j, N) * freq(k, N) / sqrt(2.0 * double(N + 1));

        return C;
    }

    void initialize(int N, Matrix &FourierComponents)
    {
        FourierComponents.setZero();

        for (int i = 0; i < N + 2; i++)
        {
            for (int j = 0; j < N + 2; j++)
            {
                FourierComponents(i, j) = sqrt(2.0 / double(N + 1)) * sin(double(i) * double(j) * pi / double(N + 1));
            }
        }
    }

    //----------------------------------------------------------------------------------------------------------------------

    // returns the total energy of the system
    double TotalEnergy(Vector &q, Vector &p, double alpha)
    {
        double Etot = 0;
        int N = q.size() - 2;

        if (p.size() != N + 2)
        {
            std::cout << "Vector length does not match";
            system("pause");
            exit(0);
        }

        for (int i = 0; i < N + 2; i++)
        {
            Etot = Etot + 0.5 * pow(p(i), 2);
            if (i < N + 1)
            {
                Etot = Etot + 0.5 * pow((q(i + 1) - q(i)), 2) + 1.0 / 3.0 * alpha * pow((q(i + 1) - q(i)), 3);
            }
        }

        return Etot;
    }

    // returns the energy of the system in the normal mode space
    double EnergyNormalSpace(const Vector &Q, const Vector &P, double alpha)
    {
        int N = Q.size() - 2;
        double Etot = 0;

        if (P.size() != N + 2)
        {
            std::cout << "Vector length does not match";
            system("pause");
            exit(0);
        }

        for (int i = 0; i < N + 2; i++)
        {
            Etot = Etot + 0.5 * pow(P(i), 2) + 0.5 * pow(freq(i, N) * Q(i), 2);
        }

        for (int j = 0; j < N + 2; j++)
        {
            for (int k = 0; k < N + 2; k++)
            {
                int i = j + k;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + alpha / 3.0 * Coupling(i, j, k, N) * Q(i) * Q(j) * Q(k);
                }
                i = j - k;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + alpha / 3.0 * Coupling(i, j, k, N) * Q(i) * Q(j) * Q(k);
                }
                i = -j + k;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + alpha / 3.0 * Coupling(i, j, k, N) * Q(i) * Q(j) * Q(k);
                }
                i = -j - k;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + alpha / 3.0 * Coupling(i, j, k, N) * Q(i) * Q(j) * Q(k);
                }
                i = 2 * (N + 1) + j + k;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + alpha / 3.0 * Coupling(i, j, k, N) * Q(i) * Q(j) * Q(k);
                }
                i = 2 * (N + 1) + j - k;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + alpha / 3.0 * Coupling(i, j, k, N) * Q(i) * Q(j) * Q(k);
                }
                i = 2 * (N + 1) - j + k;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + alpha / 3.0 * Coupling(i, j, k, N) * Q(i) * Q(j) * Q(k);
                }
                i = 2 * (N + 1) - j - k;
                if (0 <= i && i < N + 2)
                {
                    Etot = Etot + alpha / 3.0 * Coupling(i, j, k, N) * Q(i) * Q(j) * Q(k);
                }
            }
        }
        return Etot;
    }

    //-----------------------------------------------------------------------------------------------------------------------

    // computes the acceleration at every site
    void EOM(Vector &ddq, const Vector &q, double alpha)
    {
        int N = q.size() - 2;

        if (ddq.size() != N + 2)
        {
            std::cout << "Vector length does not match";
            system("pause");
            exit(0);
        }

        ddq(0) = 0;
        ddq(N + 1) = 0;
        for (int i = 1; i <= N; i++)
        {
            ddq(i) = q(i + 1) + q(i - 1) - 2.0 * q(i) + alpha * (pow(q(i + 1) - q(i), 2) - pow(q(i) - q(i - 1), 2));
        }
    }

    void EOMe(Vector &ddq, Matrix &ddqe, const Vector &q, const Matrix &qe, double alpha)
    {
        int N = q.size() - 2;

        if (ddq.size() != N + 2)
        {
            std::cout << "Vector length does not match";
            system("pause");
            exit(0);
        }

        if (ddqe.rows() != N + 2 || qe.rows() != N + 2)
        {
            std::cout << "Matrix size does not match";
            system("pause");
            exit(0);
        }

        if (ddqe.cols() != qe.cols())
        {
            std::cout << "Matrix size does not match";
            system("pause");
            exit(0);
        }

        ddq(0) = 0;
        ddq(N + 1) = 0;
        ddqe.setZero();
        for (int i = 1; i <= N; i++)
        {
            ddq(i) = q(i + 1) + q(i - 1) - 2.0 * q(i) + alpha * (pow(q(i + 1) - q(i), 2) - pow(q(i) - q(i - 1), 2));
            for (int j = 0; j < qe.cols(); j++)
            {
                ddqe(i, j) = qe(i + 1, j) + qe(i - 1, j) - 2.0 * qe(i, j) + 2.0 * alpha * ((q(i + 1) - q(i)) * (qe(i + 1, j) - qe(i, j)) - (q(i) - q(i - 1)) * (qe(i, j) - qe(i - 1, j)));
            }
        }
    }

    //----------------------------------------------------------------------------------------------------------------------

    // evolves the system for one time step
    void rkn(Vector &q, Vector &p, double alpha, double deltaT)
    {
        int N = q.size() - 2;

        if (p.size() != N + 2)
        {
            std::cout << "Vector length does not match";
            system("pause");
            exit(0);
        }

        Vector pPrev(q.size());
        Vector EOM1(q.size());
        Vector EOM2(q.size());
        Vector EOM3(q.size());
        Vector EOM4(q.size());
        Vector EOM5(q.size());

        Vector Q1(q.size());
        Vector Q2(q.size());
        Vector Q3(q.size());
        Vector Q4(q.size());
        Vector Q5(q.size());

        Q1 = q + deltaT * g1 * p;
        EOM(EOM1, Q1, alpha);
        Q2 = q + deltaT * g2 * p + pow(deltaT, 2) * a21 * EOM1;
        EOM(EOM2, Q2, alpha);
        Q3 = q + deltaT * g3 * p + pow(deltaT, 2) * a31 * EOM1 + pow(deltaT, 2) * a32 * EOM2;
        EOM(EOM3, Q3, alpha);
        Q4 = q + deltaT * g4 * p + pow(deltaT, 2) * a41 * EOM1 + pow(deltaT, 2) * a42 * EOM2 + pow(deltaT, 2) * a43 * EOM3;
        EOM(EOM4, Q4, alpha);
        Q5 = q + deltaT * g5 * p + pow(deltaT, 2) * a51 * EOM1 + pow(deltaT, 2) * a52 * EOM2 + pow(deltaT, 2) * a53 * EOM3 + pow(deltaT, 2) * a54 * EOM4;
        EOM(EOM5, Q5, alpha);

        pPrev = p;
        p = p + deltaT * b1 * EOM1 + +deltaT * b2 * EOM2 + deltaT * b3 * EOM3 + deltaT * b4 * EOM4 + deltaT * b5 * EOM5;
        q = q + deltaT * pPrev + pow(deltaT, 2) * B1 * EOM1 + pow(deltaT, 2) * B2 * EOM2 + pow(deltaT, 2) * B3 * EOM3 + pow(deltaT, 2) * B4 * EOM4 + pow(deltaT, 2) * B5 * EOM5;
    }

    // evolves system and perturbation for one time step
    void rkne(Vector &q, Matrix &qe, Vector &p, Matrix &pe, double alpha, double deltaT)
    {
        int N = q.size() - 2;

        if (p.size() != N + 2)
        {
            std::cout << "Vector length does not match";
            system("pause");
            exit(0);
        }

        if (qe.rows() != N + 2 || pe.rows() != N + 2)
        {
            std::cout << "Matrix size does not match";
            system("pause");
            exit(0);
        }

        if (qe.cols() != pe.cols())
        {
            std::cout << "Matrix size does not match";
            system("pause");
            exit(0);
        }

        Vector pPrev(q.size());
        Vector EOM1(q.size());
        Vector EOM2(q.size());
        Vector EOM3(q.size());
        Vector EOM4(q.size());
        Vector EOM5(q.size());

        Vector Q1(q.size());
        Vector Q2(q.size());
        Vector Q3(q.size());
        Vector Q4(q.size());
        Vector Q5(q.size());

        Matrix pPreve(qe.rows(), qe.cols());
        Matrix EOM1e(qe.rows(), qe.cols());
        Matrix EOM2e(qe.rows(), qe.cols());
        Matrix EOM3e(qe.rows(), qe.cols());
        Matrix EOM4e(qe.rows(), qe.cols());
        Matrix EOM5e(qe.rows(), qe.cols());

        Matrix Q1e(qe.rows(), qe.cols());
        Matrix Q2e(qe.rows(), qe.cols());
        Matrix Q3e(qe.rows(), qe.cols());
        Matrix Q4e(qe.rows(), qe.cols());
        Matrix Q5e(qe.rows(), qe.cols());

        Q1 = q + deltaT * g1 * p;
        Q1e = qe + deltaT * g1 * pe;
        EOMe(EOM1, EOM1e, Q1, Q1e, alpha);
        Q2 = q + deltaT * g2 * p + pow(deltaT, 2) * a21 * EOM1;
        Q2e = qe + deltaT * g2 * pe + pow(deltaT, 2) * a21 * EOM1e;
        EOMe(EOM2, EOM2e, Q2, Q2e, alpha);
        Q3 = q + deltaT * g3 * p + pow(deltaT, 2) * a31 * EOM1 + pow(deltaT, 2) * a32 * EOM2;
        Q3e = qe + deltaT * g3 * pe + pow(deltaT, 2) * a31 * EOM1e + pow(deltaT, 2) * a32 * EOM2e;
        EOMe(EOM3, EOM3e, Q3, Q3e, alpha);
        Q4 = q + deltaT * g4 * p + pow(deltaT, 2) * a41 * EOM1 + pow(deltaT, 2) * a42 * EOM2 + pow(deltaT, 2) * a43 * EOM3;
        Q4e = qe + deltaT * g4 * pe + pow(deltaT, 2) * a41 * EOM1e + pow(deltaT, 2) * a42 * EOM2e + pow(deltaT, 2) * a43 * EOM3e;
        EOMe(EOM4, EOM4e, Q4, Q4e, alpha);
        Q5 = q + deltaT * g5 * p + pow(deltaT, 2) * a51 * EOM1 + pow(deltaT, 2) * a52 * EOM2 + pow(deltaT, 2) * a53 * EOM3 + pow(deltaT, 2) * a54 * EOM4;
        Q5e = qe + deltaT * g5 * pe + pow(deltaT, 2) * a51 * EOM1e + pow(deltaT, 2) * a52 * EOM2e + pow(deltaT, 2) * a53 * EOM3e + pow(deltaT, 2) * a54 * EOM4e;
        EOMe(EOM5, EOM5e, Q5, Q5e, alpha);

        pPrev = p;
        pPreve = pe;
        p = p + deltaT * b1 * EOM1 + deltaT * b2 * EOM2 + deltaT * b3 * EOM3 + deltaT * b4 * EOM4 + deltaT * b5 * EOM5;
        pe = pe + deltaT * b1 * EOM1e + deltaT * b2 * EOM2e + deltaT * b3 * EOM3e + deltaT * b4 * EOM4e + deltaT * b5 * EOM5e;
        q = q + deltaT * pPrev + pow(deltaT, 2) * B1 * EOM1 + pow(deltaT, 2) * B2 * EOM2 + pow(deltaT, 2) * B3 * EOM3 + pow(deltaT, 2) * B4 * EOM4 + pow(deltaT, 2) * B5 * EOM5;
        qe = qe + deltaT * pPreve + pow(deltaT, 2) * B1 * EOM1e + pow(deltaT, 2) * B2 * EOM2e + pow(deltaT, 2) * B3 * EOM3e + pow(deltaT, 2) * B4 * EOM4e + pow(deltaT, 2) * B5 * EOM5e;
    }

    //------------------------------------------------------------------------------------------------------------------------

    // evolves the system for tmax amount of time
    void Evolve(Vector &q, Vector &p, double alpha, double tmax, double deltaT)
    {
        for (int i = 0; i < int(tmax / deltaT); i++)
        {
            rkn(q, p, alpha, deltaT);
        }
    }

    //-----------------------------------------------------------------------------------------------------------------------

    void ModeAmplitude(Vector &Q, const Vector &P, double k0, double E0, double alpha)
    {
        int N = Q.size() - 2;
        Vector coeff(3);
        coeff.setZero();
        Q(int(k0)) = 0.0;
        coeff(0) = EnergyNormalSpace(Q, P, alpha) - E0;
        for (int j = 0; j < N + 2; j++)
        {
            int i = j + int(k0);
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + alpha * Coupling(i, j, int(k0), N) * Q(i) * Q(j);
            }
            i = j - int(k0);
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + alpha * Coupling(i, j, int(k0), N) * Q(i) * Q(j);
            }
            i = -j + int(k0);
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + alpha * Coupling(i, j, int(k0), N) * Q(i) * Q(j);
            }
            i = -j - int(k0);
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + alpha * Coupling(i, j, int(k0), N) * Q(i) * Q(j);
            }
            i = 2 * (N + 1) + j + int(k0);
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + alpha * Coupling(i, j, int(k0), N) * Q(i) * Q(j);
            }
            i = 2 * (N + 1) + j - int(k0);
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + alpha * Coupling(i, j, int(k0), N) * Q(i) * Q(j);
            }
            i = 2 * (N + 1) - j + int(k0);
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + alpha * Coupling(i, j, int(k0), N) * Q(i) * Q(j);
            }
            i = 2 * (N + 1) - j - int(k0);
            if (0 <= i && i < N + 2)
            {
                coeff(1) = coeff(1) + alpha * Coupling(i, j, int(k0), N) * Q(i) * Q(j);
            }
        }
        coeff(2) = coeff(2) + 0.5 * pow(freq(int(k0), N), 2);

        int i = int(k0) + int(k0);
        if (0 <= i && i < N + 2)
        {
            coeff(2) = coeff(2) + alpha * Coupling(i, int(k0), int(k0), N) * Q(i);
        }
        i = 2 * (N + 1) - int(k0) - int(k0);
        if (0 <= i && i < N + 2)
        {
            coeff(2) = coeff(2) + alpha * Coupling(i, int(k0), int(k0), N) * Q(i);
        }

        // coeff(3) = coeff(3) + alpha / 3.0 * C(int(k0), int(k0), int(k0));
        Eigen::PolynomialSolver<double, 2> solver;
        solver.compute(coeff);
        bool realEx;
        Q(int(k0)) = solver.greatestRealRoot(realEx, 0.000001);
    }

    void Pop(Vector &Z, const Vector &Q, double k0)
    {
        Z.setZero();
        int N = Q.size() - 2;
        int j = 0;
        for (int i = 0; i < N + 2; i++)
        {
            if (i != 0 && i != N + 1 && i != int(k0))
            {
                Z(j) = Q(i);
                j = j + 1;
            }
        }
    }

    void Push(Vector &Q, const Vector &Z, double k0, double E0, double alpha)
    {
        int N = Q.size() - 2;
        Q.setZero();
        int j = 0;
        for (int i = 0; i < N + 2; i++)
        {
            if (i != 0 && i != N + 1 && i != int(k0))
            {
                Q(i) = Z(j);
                j = j + 1;
            }
        }
        Vector P(N + 2);
        P.setZero();
        ModeAmplitude(Q, P, k0, E0, alpha);
    }

    void Pushe(Vector &Qe, const Vector &Ze, const Vector &Q, double k0, double alpha)
    {
        int N = Q.size() - 2;
        Qe.setZero();
        int j = 0;
        for (int i = 0; i < N + 2; i++)
        {
            if (i != 0 && i != N + 1 && i != int(k0))
            {
                Qe(i) = Ze(j);
                j = j + 1;
            }
        }

        Vector coeff(2);
        coeff.setZero();
        coeff(0) = coeff(0) + pow(freq(k0, N), 2) * Q(int(k0));
        for (int i = 0; i < N + 2; i++)
        {
            for (int j = 0; j < N + 2; j++)
            {
                coeff(0) = coeff(0) + alpha * Coupling(i, j, int(k0), N) * Q(i) * Q(j);
            }
        }

        for (int i = 0; i < N + 2; i++)
        {
            coeff(1) = coeff(1) + pow(freq(i, N), 2) * Q(i) * Qe(i);
        }
        for (int i = 0; i < N + 2; i++)
        {
            for (int j = 0; j < N + 2; j++)
            {
                for (int k = 0; k < N + 2; k++)
                {
                    coeff(1) = coeff(1) + alpha * Coupling(i, j, k, N) * Q(i) * Q(j) * Qe(k);
                }
            }
        }
        Qe(int(k0)) = -coeff(1) / coeff(0);
    }

    //-----------------------------------------------------------------------------------------------------------------------

    double PoincareSection_rk(Vector &Z, double k0, double E0, double alpha, double tStop, double deltaT)
    {
        double period = 0;
        int N = Z.size() + 1;
        int crossing = 0;
        Vector Q(N + 2);
        Vector P(N + 2);
        Matrix FourierComponents(N + 2, N + 2);
        initialize(N, FourierComponents);
        Q.setZero();
        P.setZero();
        Push(Q, Z, k0, E0, alpha);
        double Qk0 = Q(int(k0));
        double Pk0 = P(int(k0));

        Vector q(N + 2);
        Vector p(N + 2);
        q = FourierComponents * Q;
        p = FourierComponents * P;

        for (int i = 0; i < int(tStop / deltaT); i++)
        {
            Pk0 = P(int(k0));
            rkn(q, p, alpha, deltaT);
            Q = FourierComponents * q;
            P = FourierComponents * p;
            period = period + deltaT;

            if (i >= 2 && sgn(Q(int(k0))) == sgn(Qk0) && sgn(P(int(k0))) != sgn(Pk0))
            {
                crossing = crossing + 1;
            }
            if (crossing > 0)
            {
                break;
            }
        }

        Pop(Z, Q, k0);
        return period;
    }

    double PoincareSectione_rk(Matrix &Qe, Matrix &Pe, const Vector Q0, double k0, double E0, double alpha, double tStop, double deltaT)
    {
        int N = Q0.size() - 2;
        int crossing = 0;
        double tCurr = 0;
        Vector Q(N + 2);
        Vector P(N + 2);
        Matrix FourierComponents(N + 2, N + 2);
        initialize(N, FourierComponents);
        Q = Q0;
        P.setZero();

        double Qk0 = Q(int(k0));
        double Pk0 = P(int(k0));

        Vector q(N + 2);
        Vector p(N + 2);
        q = FourierComponents * Q;
        p = FourierComponents * P;

        Matrix qe(Qe.rows(), Qe.cols());
        Matrix pe(Qe.rows(), Qe.cols());
        qe = FourierComponents * Qe;
        pe = FourierComponents * Pe;

        for (int i = 0; i < int(tStop / deltaT); i++)
        {
            Pk0 = P(int(k0));

            rkne(q, qe, p, pe, alpha, deltaT);
            Q = FourierComponents * q;
            P = FourierComponents * p;
            Qe = FourierComponents * qe;
            Pe = FourierComponents * pe;

            tCurr = tCurr + deltaT;

            if (i >= 2 && sgn(Q(int(k0))) == sgn(Qk0) && sgn(P(int(k0))) != sgn(Pk0))
            {
                crossing = crossing + 1;
            }
            if (crossing > 0)
            {
                break;
            }
        }
        return tCurr;
    }

    void PerturbationMatrix(Matrix &QeM, Matrix &PeM, const Vector &Q, double k0, double alpha)
    {
        int N = Q.size() - 2;
        QeM.setZero();
        PeM.setZero();
        Vector Ze(N - 1);
        Vector Qe(N + 2);
        Ze.setZero();
        Qe.setZero();
        int k = 0;
        for (int i = 0; i < N + 2; i++)
        {
            if (i != 0 && i != N + 1 && i != int(k0))
            {
                Ze.setZero();
                Ze(k) = 1.0;
                Pushe(Qe, Ze, Q, k0, alpha);
                for (int j = 0; j < N + 2; j++)
                {
                    QeM(j, i) = Qe(j);
                }
                k = k + 1;
            }
        }
    }

    //-----------------------------------------------------------------------------------------------------------------------

    void initialState(Vector &q, Vector &p, Vector &Q, Vector &P, int N, double E0, double k0)
    {
        double Amp = sqrt(2.0 * E0) / freq(k0, N);
        for (int i = 0; i < N + 2; i++)
        {
            q(i) = sqrt(2.0 / (N + 1)) * Amp * sin(pi * double(i) * double(k0) / double(N + 1));
            p(i) = 0;
            Q(i) = 0;
            P(i) = 0;
        }

        Q[int(k0)] = Amp;
    }

    double Mod(const Vector &X)
    {
        return ((X.array()).abs()).maxCoeff();
    }

    //----------------------------------------------------------------------------------------------------------------------------

    double ent(const Vector &E)
    {
        double Ent = 0;
        int N = E.size() - 2;
        double Etot = 0;
        for (int i = 1; i < N + 1; i++)
        {
            Etot = Etot + E(i);
        }
        for (int i = 1; i < N + 1; i++)
        {
            if (E(i) > pow(10, -15))
            {
                Ent = Ent - (E(i) / Etot) * std::log((E(i) / Etot));
            }
        }
        return Ent;
    }

    void Evolution(Vector &time, Vector &Entropy, Vector &q, Vector &p, double alpha, double tmax, double deltaT)
    {
        int N = q.size() - 2;
        int frames = time.size();
        time.setZero();
        Entropy.setZero();

        Vector Q(N + 2);
        Vector P(N + 2);
        Vector E(N + 2);
        E.setZero();
        Matrix FourierComponents(N + 2, N + 2);
        FourierComponents.setZero();
        initialize(N, FourierComponents);
        Q = FourierComponents * q;
        P = FourierComponents * p;

        Vector freqs(N + 2);
        freqs.setZero();
        for (int i = 1; i < N + 1; i++)
        {
            freqs(i) = freq(i, N);
        }

        E = 0.5 * P.cwiseProduct(P) + 0.5 * (Q.cwiseProduct(freqs)).cwiseProduct(Q.cwiseProduct(freqs));

        int imax = int(tmax / deltaT);
        int interval = (imax / frames);
        int frame = 0;
        double tCurr = 0;

        for (int i = 0; i < int(tmax / deltaT); i++)
        {
            if (i % interval == 0 && frame < frames)
            {
                time(frame) = tCurr;
                Q = FourierComponents * q;
                P = FourierComponents * p;

                E = 0.5 * P.cwiseProduct(P) + 0.5 * (Q.cwiseProduct(freqs)).cwiseProduct(Q.cwiseProduct(freqs));

                Entropy(frame) = ent(E);
                frame = frame + 1;
            }

            rkn(q, p, alpha, deltaT);
            tCurr = tCurr + deltaT;
        }
    }

}

//---------------------------------------------------------------------------------------------------------------------------------
