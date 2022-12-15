#include <iostream>
#include <cmath>
#include <functional>
#include <random>
#include <ctime>
#include <iomanip>

const double pi = 3.14159265358979324;

std::mt19937 mt_rand(std::time(0));

/* uniformly distributed random variable in the interval [-1..1] */
double uniform_srand()
{
    return (1.0 / 2147483648.0) * int(mt_rand());
}

/* random variable uniformly distributed in the ball x^2 + y^2 < 1 */
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

/* normally distributed random variable in the plane, with variance = 1 */
void normal_rand2d(double& x, double& y)
{
    double r = ball_rand2d(x, y);
    double f = std::sqrt(-2.0 * std::log(r) / r);
    x *= f;
    y *= f;
}

struct phasepoint_t
{
    double x, y, z;
    double p, q, r;

    phasepoint_t() {}
    phasepoint_t(double x, double y, double z, double p, double q, double r) : x(x), y(y), z(z), p(p), q(q), r(r) {}

    void initialize(double z, double alpha, double speed)
    {
        x = 0.0;
        y = 0.0;
        phasepoint_t::z = z;
        p = speed * std::cos(alpha);
        q = 0.0;
        r = speed * std::sin(alpha);
    }
};

phasepoint_t operator * (double a, const phasepoint_t& U)
    { return phasepoint_t(a * U.x, a * U.y, a * U.z, a * U.p, a * U.q, a * U.r); }

phasepoint_t operator + (const phasepoint_t& U, const phasepoint_t& V)
    { return phasepoint_t(U.x + V.x, U.y + V.y, U.z + V.z, U.p + V.p, U.q + V.q, U.r + V.r); }

phasepoint_t operator - (const phasepoint_t& U, const phasepoint_t& V)
    { return phasepoint_t(U.x - V.x, U.y - V.y, U.z - V.z, U.p - V.p, U.q - V.q, U.r - V.r); }

//=================================================================================================================
//   Air resistance creates a force that is always directed against the direction of motion in the surrounding
// medium and has a magnitude that depends on the absolute speed: F_air = − f(|v|) * v. The speed-dependence of
// the friction force is linear at very low speeds (Stokes drag) and quadratic at larger speeds (Newton drag).
// The transition between these behaviours is determined by the Reynolds number, which depends on speed,
// object size and kinematic viscosity of the medium. For Reynolds numbers below about 1000, the dependence is
// linear, above it becomes quadratic. In air, which has a kinematic viscosity around 0.000015 m^2 / sec, this
// means that the drag force becomes quadratic in v when the product of speed and diameter is more than about
// 0.015 m^2 / s, which is typically the case for projectiles.
//=================================================================================================================

/* this function will define smooth transition between Stokes drag and Newton drag, mixing the two formulas,
   the transition will happen according as described above, near 0.015 m^2 / s" */

double smootherstep(double edge0, double edge1, double x)
{
    x = (x - edge0) / (edge1 - edge0);
    x = std::max(0.0, std::min(1.0, x));
    return x * x * x * (x * (x * 6.0 - 15.0) + 10.0);
}

/* Drag coefficients for typical objects,
   from https://en.wikipedia.org/wiki/Drag_coefficient */

const double c0 = 0.1;          /* Smooth sphere (Re = 10^6) */
const double c1 = 0.47;         /* Rough sphere (Re = 10^6) */
const double c2 = 0.81;         /* Triangular trapeze (45°) */
const double c3 = 1.3;          /* Trapeze with triangular basis (45°) */
const double c4 = 0.295;        /* Bullet (not ogive, at subsonic velocity) */

const double rho_0 = 1.2754;     /* air density near Earth surface at 101.325 kPa and 20°C, in kg / m3 */

const double TP = 0.015;        /* transition point, where air resistance becomes quadratic */
const double TD = 0.002;        /* transition delta, transition will happen in the interval [TP - TD, TP + TD] */
const double Tmin = TP - TD;
const double Tmax = TP + TD;

/* Dynamic viscosity of air for different temperatures
   https://www.engineersedge.com/physics/viscosity_of_air_dynamic_and_kinematic_14483.htm */

const double mu0 = 1.579e-5;                /* -30°C */
const double mu1 = 1.630e-5;                /* -20°C */
const double mu2 = 1.689e-5;                /* -10°C */
const double mu3 = 1.729e-5;                /*   0°C */
const double mu4 = 1.754e-5;                /*   5°C */
const double mu5 = 1.778e-5;                /*  10°C */
const double mu6 = 1.802e-5;                /*  15°C */


struct params_t
{
    /* universal constants */
    const double g = 9.81;                          /* small gravitational constant */
    const double L = -0.0065;                       /* the rate of temperature decay in the troposphere, in °К / m */
    const double R = 8.314462618153;                /* universal gas constant, in J / (K * mol) */
    const double M = 0.0289644;                     /* molar mass of the air, in kg / mol */
    const double T_0 = 288.15;                      /* standard temperature at the sea level, in °К */
    const double P_0 = 101325.0;                    /* standard air pressure at the sea level, in Pa */
    const double K = g * M / (R * L);               /* combination of the constants, used in calculations */

    double h_0;                 /* altitude above sea level */

    double c;                   /* the drag coefficient */
    double mu;                  /* the value of dynamic viscosity that we will use in calculations */

    double D, A;                /* diameter and area of the (circular) cross section of the ballistic shell */

    double w_x, w_y, w_z;       /* wind component */
    double m, inv_m;            /* mass of the body and its inverse */

    void init_altitude(double h_0)
    {
        params_t::h_0 = h_0;
    }

    void init_mass(double m)
    {
        params_t::m = m;
        inv_m = 1.0 / m;
    }

    void init_diameter(double D)
    {
        params_t::D = D;
        A = 0.25 * pi * D * D;
    }

    /* choose some constant normally distributed vector in xy plane with some very small z component */
    void init_random_wind(double delta_xy, double delta_z)
    {
        double dummy;
        normal_rand2d(w_x, w_y);
        normal_rand2d(w_z, dummy);
        w_x *= delta_xy;
        w_y *= delta_xy;
        w_z *= delta_z;
    }

    void init_air_resistance(double c, double mu)
    {
        params_t::c = c;
        params_t::mu = mu;
    }
};

void RK_iterate(bool single_run,
        double step,
        const params_t& params,                                                 /* calculation parameters */
        const phasepoint_t& ip,                                                 /* initial position and velocity */
        std::function<phasepoint_t(const params_t&, const phasepoint_t&)> F     /* function that defines the differential equation */
    )
{
    double t = 0.0;
    double h = step;
    double h2 = 0.5 * h;
    double h6 = h / 6.0;

    double max_z = 0.0;

    phasepoint_t p = ip;

    if (single_run)
    {
        /* header */
        std::cout << "# ================================================================================= " << std::endl;
        std::cout << "# w_x = " << params.w_x << ", w_y = " << params.w_y << ", w_z = " << params.w_z << std::endl;
        std::cout << "# m = " << params.m << std::endl;
        std::cout << "# c = " << params.c << std::endl;
        std::cout << "# mu = " << params.mu << ", c = " << params.c << std::endl;
        std::cout << "# ================================================================================= " << std::endl;
        std::cout << "# v0_x = " << ip.p << std::endl;
        std::cout << "# v0_z = " << ip.r << std::endl;
        std::cout << "# ================================================================================= " << std::endl;
        std::cout << std::endl;
    }

    do {
        /* main data output */
        if (single_run)
            std::cout << p.x << " " << p.y << " " << p.z << std::endl;

        t += h;

        phasepoint_t K1 = F(params, p);
        phasepoint_t K2 = F(params, p + h2 * K1);
        phasepoint_t K3 = F(params, p + h2 * K2);
        phasepoint_t K4 = F(params, p + h * K3);

        p = p + h6 * (K1 + 2.0 * K2 + 2.0 * K3 + K4);

        max_z = std::max(max_z, p.z);

    } while (p.z > 0.0);

    double cos_alpha = std::sqrt(p.p * p.p + p.q * p.q) / std::sqrt(p.p * p.p + p.q * p.q + p.r * p.r);
    double distance = std::sqrt(p.x * p.x + p.y * p.y);
    double alpha = 180.0 * std::acos(cos_alpha) / pi;

    if (single_run)
    {
        /* footer */
        std::cout << std::endl;
        std::cout << "# ================================================================================= " << std::endl;
        std::cout << "# max height = " << max_z << std::endl;
        std::cout << "# t = " << t << std::endl;
        std::cout << "# distance = " << distance << std::endl;
        std::cout << "# angle = " << alpha << std::endl;
        std::cout << "# ================================================================================= " << std::endl;
    }
    else
        std::cout << ", max height: " << max_z << ", time: " << t << ", distance: " << distance << ", arrival angle: " << alpha << std::endl;
}

int main(int argc, char* argv[])
{
    /* if you launched the program without params that means you want a single simulation,
       otherwise we will create a readable table for you */
    bool single_run = argc == 1;

    /* paramaters */
    params_t params;
    params.init_altitude(169.0);            /* Makeevka is 169 meters above sea level */
    params.init_mass(6.3);                  /* mass = 6.3 kg */
    params.init_diameter(0.076);            /* 76 mm = 0.076 m */
    params.init_random_wind(0.5, 0.001);    /* random wind */
    params.init_air_resistance(c1, mu3);    /* some predefined values: c1 - rough sphere, mu3 - viscosity at 0°C */

    /* equation */
    auto F = [](const params_t& params, const phasepoint_t& V) -> phasepoint_t
    {
        double s = std::sqrt(V.p * V.p + V.q * V.q + V.r * V.r);                    /* speed */
        double T = std::max(params.T_0 + params.L * (V.z + params.h_0), 0.01);      /* make sure temperature is never zero */
        double P = params.P_0 * std::exp(-params.K * std::log(T / params.T_0));     /* air pressure at the current height */
        double rho = P * params.M / (params.R * T);                                 /* air density at the current height */

        //=========================================================================================================
        // Stokes drag, small speeds
        //=========================================================================================================
        //
        //    The force of viscosity on a small sphere moving through a viscous fluid is given by:
        //        F = 6πμRv = 3πμDv
        //    where:
        //
        //        F is the frictional force – known as Stokes' drag – acting on the interface between
        //            the fluid and the particle
        //        μ is the dynamic viscosity (some authors use the symbol η)
        //        D is the diameter of the spherical object
        //        v is the flow velocity relative to the object
        //
        //=========================================================================================================
        double SD = 3.0 * pi * params.mu * params.D;

        //=========================================================================================================
        // Newton drag, large speeds
        //=========================================================================================================
        //
        //      Air resistance is given by the following drag formula:
        //        F = − 1/2cρAsv
        //    where
        //        F is the drag force
        //        c is the drag coefficient
        //        ρ is the air density
        //        A is the cross sectional area of the projectile
        //        v is the velocity (vector)
        //        s is the speed (the number), e.g. the norm of v
        //
        //=========================================================================================================
        double ND = 0.5 * params.c * rho * params.A * s;

        double RN = params.R * s;                   /* not a Reynolds number, but something closely related to it */
        double q = smootherstep(Tmin, Tmax, RN);    /* q is the proportion of the Newton drag in the total drag */
        double DF = (1.0 - q) * SD + q * ND;        /* combined force, yet to be multiplied by the velocity vector */

        //=========================================================================================================
        //
        //      The system of differential equations that we will solve is
        //        m * x''(t) =          w_x - DF * x'(t)
        //        m * y''(t) =          w_y - DF * y'(t)
        //        m * z''(t) = -m * g + w_z - DF * z'(t)
        //    where
        //        DF is a combination of Stokes drag force and Newton drag force
        //      Newton drag force contains speed factor and therefore, in total, gives quadratic contribution into
        //    the differential equation
        //
        //      Rewrite it as a first order system to apply Runge-Kutta method
        //  https://ru.wikipedia.org/wiki/%D0%9C%D0%B5%D1%82%D0%BE%D0%B4_%D0%A0%D1%83%D0%BD%D0%B3%D0%B5_%E2%80%94_%D0%9A%D1%83%D1%82%D1%82%D1%8B
        //
        //        x'(t) = p
        //        y'(t) = q
        //        z'(t) = r
        //        p'(t) =      (w_x - DF * p) / m
        //        q'(t) =      (w_y - DF * q) / m
        //        r'(t) = -g + (w_z - DF * r) / m
        //
        //=================================================================================================================

        return phasepoint_t(
            V.p,
            V.q,
            V.r,
            (params.w_x - DF * V.p) * params.inv_m,
            (params.w_y - DF * V.q) * params.inv_m,
            (params.w_z - DF * V.r) * params.inv_m - params.g
        );
    };

    /* initial conditions */
    phasepoint_t ip;
    const double step = 1.0 / 2048.0;
    const double z0 = 0.5;              /* fire from 0.5 meters above the ground */
    const double speed = 662.0;         /* with 662 m / s speed */

    if (!single_run)   /* you want a table */
    {
        std::cout << std::setprecision(2) << std::fixed;

        for (double angle = 15.0; angle < 88.0; angle += 0.5)
        {
            std::cout << "departure angle: " << angle;
            ip.initialize(z0, pi * angle / 180.0, speed);

            /* Runge-Kutta iterative solution */
            RK_iterate(single_run, step, params, ip, F);
        }
    }
    else                /* single simulation */
    {
        const double angle = 67.0;
        ip.initialize(z0, pi * angle / 180.0, speed);          /* fire from 0.5 meters above the ground */
        RK_iterate(single_run, step, params, ip, F);
    }

    return 0;
}
