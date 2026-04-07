//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "smooth-power-law-syn.h"

#include <algorithm>

namespace {
inline Real clamp_smoothing(Real s) noexcept {
    return std::max(s, 0.05);
}

inline Real log2_sum(Real a, Real b) noexcept {
    const Real hi = std::max(a, b);
    const Real lo = std::min(a, b);
    return hi + log2_softplus(lo - hi);
}

inline Real gs02_s12_slow(Real k, Real p) noexcept {
    return clamp_smoothing(1.84 - 0.040 * k - (0.40 - 0.010 * k) * p);
}

inline Real gs02_s23_slow(Real k, Real p) noexcept {
    return clamp_smoothing(1.15 - 0.125 * k - (0.06 - 0.015 * k) * p);
}

inline Real gs02_s12_fast(Real /*k*/, Real /*p*/) noexcept {
    return 0.597;
}

inline Real gs02_s23_fast(Real k, Real p) noexcept {
    return clamp_smoothing(3.34 + 0.17 * k - (0.82 + 0.035 * k) * p);
}
} // namespace

//========================================================================================================
//                                  SmoothPowerLawSyn Class Methods
//========================================================================================================
inline Real log2_smooth_one(Real log2_a, Real log2_b) noexcept {
    return log2_a - log2_softplus(log2_a - log2_b);
}

Real SmoothPowerLawSyn::log2_optical_thin(Real log2_nu) const noexcept {
    const Real log2_x12 = log2_nu - log2_nu12_;
    const Real log2_x23 = log2_nu - log2_nu23_;

    const Real term1 = (s23_ / s12_) * log2_softplus(-s12_ * (b1_ - b2a_) * log2_x12) - s23_ * b2a_ * log2_x12;
    const Real term2 = -s23_ * b2b_ * (log2_nu23_ - log2_nu12_) - s23_ * b3_ * log2_x23;

    return -log2_sum(term1, term2) / s23_;
}

Real SmoothPowerLawSyn::log2_optical_thick(Real log2_nu) const noexcept {
    const Real log2_x = log2_nu - log2_nu_m;
    const Real s = -smooth_a_ * fast_exp2(2. / 3 * log2_x);

    return 2.5 * log2_x + log2_softplus(-0.5 * log2_x + s);
}

Real SmoothPowerLawSyn::log2_optical_thick_sharp(Real log2_nu) const noexcept {
    if (log2_nu < log2_nu_m) {
        return 2. * (log2_nu - log2_nu_m);
    } else {
        return 2.5 * (log2_nu - log2_nu_m);
    }
}

Real SmoothPowerLawSyn::log2_optical_thin_sharp(Real log2_nu) const noexcept {
    if (log2_nu_m < log2_nu_c) {
        if (log2_nu < log2_nu_m) {
            return (log2_nu - log2_nu_m) / 3.0;
        } else if (log2_nu < log2_nu_c) {
            return 0.5 * (1.0 - p) * (log2_nu - log2_nu_m);
        } else {
            return 0.5 * (1.0 - p) * (log2_nu_c - log2_nu_m) - 0.5 * p * (log2_nu - log2_nu_c);
        }
    } else {
        if (log2_nu < log2_nu_c) {
            return (log2_nu - log2_nu_c) / 3.0;
        } else if (log2_nu < log2_nu_m) {
            return -0.5 * (log2_nu - log2_nu_c);
        } else {
            return -0.5 * (log2_nu_m - log2_nu_c) - 0.5 * p * (log2_nu - log2_nu_m);
        }
    }
}

Real SmoothPowerLawSyn::compute_spectrum(Real nu) const noexcept {
    return fast_exp2(compute_log2_spectrum(fast_log2(nu)));
}

Real SmoothPowerLawSyn::compute_log2_spectrum(Real log2_nu) const noexcept {
    Real log2_f_thin = log2_optical_thin(log2_nu);
    Real log2_f_thick = log2_optical_thick(log2_nu);
    return log2_norm_ + log2_smooth_one(log2_f_thin, log2_f_thick + log2_thick_norm_);
}

void SmoothPowerLawSyn::build() noexcept {
    log2_I_nu_max = fast_log2(I_nu_max);
    log2_nu_m = fast_log2(nu_m);
    log2_nu_c = fast_log2(nu_c);
    log2_nu_a = fast_log2(nu_a);
    log2_nu_M = fast_log2(nu_M);

    constexpr Real ln2 = std::numbers::ln2;
    smooth_a_ = (3.5 * p - 1.5) / ln2;
    inv_nu_M_ = 1.0 / nu_M;
    log2_norm_ = 0;

    log2_nu12_ = std::min(log2_nu_m, log2_nu_c);
    log2_nu23_ = std::max(log2_nu_m, log2_nu_c);

    const Real k = std::isfinite(k_eff) ? k_eff : 2.0;
    const Real s12_slow = gs02_s12_slow(k, p);
    const Real s23_slow = gs02_s23_slow(k, p);
    const Real s12_fast = gs02_s12_fast(k, p);
    const Real s23_fast = gs02_s23_fast(k, p);

    b1_ = 1.0 / 3.0;
    b3_ = -0.5 * p;

    const bool fast = nu_m >= nu_c;
    const Real s12_initial = fast ? s12_fast : s12_slow;
    const Real s23_initial = fast ? s23_fast : s23_slow;
    const Real nu_ratio = std::max(nu_m / std::max(nu_c, 1e-300), 1e-300);
    const Real delta_b = 0.5 * p + 1.0 / 3.0;
    const Real q12 = clamp_smoothing(s12_initial * delta_b);
    const Real q23 = clamp_smoothing(s23_initial * delta_b);

    s12_ = s12_fast + (s12_slow - s12_fast) / (1.0 + fast_pow(nu_ratio, q12));
    s23_ = s23_fast + (s23_slow - s23_fast) / (1.0 + fast_pow(nu_ratio, q23));
    b2a_ = -0.5 + (((1.0 - p) / 2.0) + 0.5) / (1.0 + fast_pow(nu_ratio, q12));
    b2b_ = -0.5 + (((1.0 - p) / 2.0) + 0.5) / (1.0 + fast_pow(nu_ratio, q23));

    log2_thick_norm_ = log2_optical_thin(log2_nu_a) - log2_optical_thick(log2_nu_a);
}

Real SmoothPowerLawSyn::compute_I_nu(Real nu) const noexcept {
    if (nu <= nu_c) { // Below cooling frequency, simple scaling
        return fast_exp(-nu * inv_nu_M_) * I_nu_max * compute_spectrum(nu);
    } else {
        return fast_exp(-nu * inv_nu_M_) * I_nu_max * compute_spectrum(nu) * inverse_compton_correction(*this, nu);
    }
}

Real SmoothPowerLawSyn::compute_log2_I_nu(Real log2_nu) const noexcept {
    constexpr Real log2e = std::numbers::log2e;
    if (log2_nu <= log2_nu_c) { // Below cooling frequency, simple scaling
        return log2_I_nu_max + compute_log2_spectrum(log2_nu) - log2e * inv_nu_M_ * fast_exp2(log2_nu);
    } else {
        const Real nu = fast_exp2(log2_nu);
        return log2_I_nu_max + compute_log2_spectrum(log2_nu) - log2e * inv_nu_M_ * nu +
               fast_log2(inverse_compton_correction(*this, nu));
    }
}
