//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "physics.h"

#include <algorithm>

namespace physics::scales { 

    Real dec_radius(Real E_iso, Real n0, Real Gamma0, Real engine_dura, Real k, Real r0) {
        return std::max(thin_shell_dec_radius(E_iso, n0, Gamma0, k, r0),
                        thick_shell_dec_radius(E_iso, n0, engine_dura, k, r0));
    }

    // R_dec = [(3-k) E_iso / (4π n0 mp r0^k c² Γ₀²)]^{1/(3-k)}
    // At k=0, r0=1: reduces to cbrt(3 E / (4π n mp c² Γ₀²))
    Real thin_shell_dec_radius(Real E_iso, Real n0, Real Gamma0, Real k, Real r0) {
        const Real r0k   = std::pow(r0, k);
        const Real inner = (3 - k) * E_iso / (4 * con::pi * n0 * con::mp * r0k * con::c2 * Gamma0 * Gamma0);
        return std::pow(inner, 1.0 / (3 - k));
    }

    // R_dec = [(3-k) E_iso c T / (4π n0 mp r0^k c²)]^{1/(4-k)}
    // At k=0, r0=1: reduces to (3 E c T / (4π n mp c²))^{1/4}
    Real thick_shell_dec_radius(Real E_iso, Real n0, Real engine_dura, Real k, Real r0) {
        const Real r0k   = std::pow(r0, k);
        const Real inner = (3 - k) * E_iso * con::c * engine_dura / (4 * con::pi * n0 * con::mp * r0k * con::c2);
        return std::pow(inner, 1.0 / (4 - k));
    }

    Real shell_spreading_radius(Real Gamma0, Real engine_dura) {
        return Gamma0 * Gamma0 * con::c * engine_dura;
    }

    Real RS_transition_radius(Real E_iso, Real n_ism, Real Gamma0, Real engine_dura) {
        return std::pow(sedov_length(E_iso, n_ism), 1.5) / std::sqrt(con::c * engine_dura) / Gamma0 / Gamma0;
    }

    Real shell_thickness_param(Real E_iso, Real n_ism, Real Gamma0, Real engine_dura) {
        const Real Sedov_l = sedov_length(E_iso, n_ism);
        const Real shell_width = con::c * engine_dura;
        return std::sqrt(Sedov_l / shell_width) * std::pow(Gamma0, -4. / 3);
    }

    Real calc_engine_duration(Real E_iso, Real n_ism, Real Gamma0, Real xi) {
        const Real Sedov_l = sedov_length(E_iso, n_ism);
        return Sedov_l / (xi * xi * std::pow(Gamma0, 8. / 3) * con::c);
    }

} // namespace physics::scales
