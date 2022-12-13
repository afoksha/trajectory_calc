#include <iostream>
#include <cmath>
#include <functional>
#include <random>
#include <ctime>

//====================================================================================================
//  n^2 = (x'(t)) ^ 2 + (y'(t)) ^ 12 + (z'(t)) ^ 2
//  rho = exp(-c*z)
//
//
//  m * x''(t) =          w_x - (k_1 + k_2 * n + k_3 * n ^ 2) * x'(t)
//  m * y''(t) =          w_y - (k_1 + k_2 * n + k_3 * n ^ 2) * y'(t)
//  m * z''(t) = -m * g + w_z - (k_1 + k_2 * n + k_3 * n ^ 2) * z'(t)
//
//  Rewrite it as a first order system to apply Runge-Kutta method
//  https://ru.wikipedia.org/wiki/%D0%9C%D0%B5%D1%82%D0%BE%D0%B4_%D0%A0%D1%83%D0%BD%D0%B3%D0%B5_%E2%80%94_%D0%9A%D1%83%D1%82%D1%82%D1%8B
//
//  x'(t) = p
//  y'(t) = q
//  z'(t) = r
//  p'(t) =      (w_x - (k_1 + k_2 * n + k_3 * n ^ 2) * p) / m
//  q'(t) =      (w_y - (k_1 + k_2 * n + k_3 * n ^ 2) * q) / m
//  r'(t) = -g + (w_z - (k_1 + k_2 * n + k_3 * n ^ 2) * r) / m
//
//====================================================================================================

const double pi = 3.14159265358979324;

std::mt19937 mt_rand(std::time(0));

double uniform_srand()
{
    return (1.0 / 2147483648.0) * int(mt_rand());
}

double ball_rand2d(double& x, double& y)
{
    double r;
    do
    {
        x = uniform_srand();
        y = uniform_srand();
        r = x * x + y * y;
    }
    while (r >= 1.0);
    return r;
}

void normal_rand2d(double& x, double& y)
{
    double r = ball_rand2d(x, y);
    double f = std::sqrt(-2.0 * std::log(r) / r);
    x *= f;
    y *= f;
}

struct phasepoint
{
    double x, y, z;
    double p, q, r;

    phasepoint() {}
    phasepoint(double x, double y, double z, double p, double q, double r) : x(x), y(y), z(z), p(p), q(q), r(r) {}

    void initialize(double alpha, double speed)
    {
        x = 0.0;
        y = 0.0;
        z = 0.0;
        p = speed * std::cos(alpha);
        q = 0.0;
        r = speed * std::sin(alpha);
    }
};

phasepoint operator * (double a, const phasepoint& U)
    { return phasepoint(a * U.x, a * U.y, a * U.z, a * U.p, a * U.q, a * U.r); }

phasepoint operator + (const phasepoint& U, const phasepoint& V)
    { return phasepoint(U.x + V.x, U.y + V.y, U.z + V.z, U.p + V.p, U.q + V.q, U.r + V.r); }

phasepoint operator - (const phasepoint& U, const phasepoint& V)
    { return phasepoint(U.x - V.x, U.y - V.y, U.z - V.z, U.p - V.p, U.q - V.q, U.r - V.r); }


struct params
{
    const double g = 9.81;
    const double c = 1.25;

    double w_x, w_y, w_z;
    double m, inv_m;
    double k_1, k_2, k_3;

    void init_mass(double m)
    {
        params::m = m;
        inv_m = 1.0 / m;
    }

    /* choose some constant normally distributed vector in xy plane */
    void init_random_wind(double delta_xy, double delta_z)
    {
        double dummy;
        normal_rand2d(w_x, w_y);
        normal_rand2d(w_z, dummy);
        w_x *= delta_xy;
        w_y *= delta_xy;
        w_z *= delta_z;
    }

    void init_air_resistance(double k_1, double k_2, double k_3)
    {
        params::k_1 = k_1;
        params::k_2 = k_2;
        params::k_3 = k_3;
    }
};

void RK_iterate(
        const params& s,                                                /* calculation parameters */
        const phasepoint& ip,                                           /* initial position and velocity */
        std::function<phasepoint(const params&, const phasepoint&)> F   /* function that defines the differential equation */
    )
{
    double t = 0.0;
    double h6 = 1.0 / 1024.0;   /* h / 6 in Runge-Kutta iteration formulas */
    double h3 = 3.0 * h6;       /* h / 2 in Runge-Kutta iteration formulas */
    double h1 = 6.0 * h6;       /* h, the equation step in Runge-Kutta formula */

    double max_z = 0.0;

    phasepoint p = ip;

    /* header */
    std::cout << "# ================================================================================= " << std::endl;
    std::cout << "# w_x = " << s.w_x << ", w_y = " << s.w_y << ", w_z = " << s.w_z << std::endl;
    std::cout << "# m = " << s.m << std::endl;
    std::cout << "# c = " << s.c << std::endl;
    std::cout << "# k_1 = " << s.k_1 << ", k_2 = " << s.k_2 << "# k_3 = " << s.k_3 << std::endl;
    std::cout << "# ================================================================================= " << std::endl;
    std::cout << "# v0_x = " << ip.p << std::endl;
    std::cout << "# v0_z = " << ip.r << std::endl;
    std::cout << "# ================================================================================= " << std::endl;
    std::cout << std::endl;

    do {
        std::cout << p.x << " " << p.y << " " << p.z << std::endl;

        t += h1;

        phasepoint K1 = F(s, p);
        phasepoint K2 = F(s, p + h3 * K1);
        phasepoint K3 = F(s, p + h3 * K2);
        phasepoint K4 = F(s, p + h1 * K3);

        p = p + h6 * (K1 + 2.0 * K2 + 2.0 * K3 + K4);

        max_z = std::max(max_z, p.z);

    } while (p.z > 0.0);

    double cos_alpha = std::sqrt(p.p * p.p + p.q * p.q) / std::sqrt(p.p * p.p + p.q * p.q + p.r * p.r);
    double distance = std::sqrt(p.x * p.x + p.y * p.y);
    double alpha = 180.0 * std::acos(cos_alpha) / pi;

    std::cout << std::endl;
    std::cout << "# ================================================================================= " << std::endl;
    std::cout << "# max height = " << max_z << std::endl;
    std::cout << "# t = " << t << std::endl;
    std::cout << "# distance = " << distance << std::endl;
    std::cout << "# angle = " << alpha << std::endl;
    std::cout << "# ================================================================================= " << std::endl;
}

int main(int argc, char* argv[])
{
    /* paramaters */
    params s;
    s.init_mass(6.3);                   /* mass = 6.3 kg */
    s.init_random_wind(0.5, 0.001);     /* random wind */
    s.init_air_resistance(0.123, 0.401, 0.00127);

    /* equation */
    auto F = [](const params& S, const phasepoint& V) -> phasepoint
    {
        double n2 = V.p * V.p + V.q * V.q + V.r * V.r;
        double n1 = std::sqrt(n2);
        double rho = std::exp(-S.c * V.z);

        double R = rho * (S.k_1 + S.k_2 * n1 + S.k_3 * n2);

        return phasepoint(
            V.p,
            V.q,
            V.r,
            (S.w_x - R * V.p) * S.inv_m,
            (S.w_y - R * V.q) * S.inv_m,
            (S.w_z - R * V.r) * S.inv_m - S.g
        );
    };

    /* initial conditions */
    phasepoint ip;
    ip.initialize(pi * 0.37, 662);       /* pi * 0.4 = 72 degrees */

    RK_iterate(s, ip, F);

    return 0;
}
