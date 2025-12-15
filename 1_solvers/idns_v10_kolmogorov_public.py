#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
idns_v10_kolmogorov.py — iDNS with Corrected Forcing Sign
-------------------------------------------------------------------------
Intelligent Direct Numerical Simulation (iDNS) for 2D Kolmogorov Flow

Version: 10.0.0
Date: 14 December 2025
Author: Jeffrey Camlin
ORCID: 0000-0002-5740-4204

Code and data: doi:10.5281/zenodo.17730873
GitHub: RED-DAWN-AI/iDNS

Copyright (c) 2025 Red Dawn Academic Press & AI Lab, Milwaukee WI
License: CC BY-NC 4.0 (Attribution-NonCommercial)
https://creativecommons.org/licenses/by-nc/4.0/

-------------------------------------------------------------------------
v10 FIXES (2024-12-13):
  - CRITICAL: Forcing sign corrected to match C&K eq (2.6)
  - curl(F·sin(kf·y)·x̂) = -F·kf·cos(kf·y)  (was +F·kf·cos(kf·y) in v9)
  - This fix brings D/D_lam from ~0.05 to ~0.4-0.5 matching C&K benchmark

v9 FIXES (2024-12-12):
  - Stream function sign corrected: ψ̂ = +ω̂/k² (was -ω̂/k²)
  - Spectral derivatives replace np.gradient for consistency

Features:
  - RASI temporal lifting (φ′ adaptive controller)
  - Cascade-Completion Clamping (CCC) at spectral cutoff
  - Two CCC modes: periodic (--ccc periodic) or adaptive (--ccc adaptive)
  - Unique run IDs to prevent file overwrites
"""

import argparse, os, csv, json, time, sys
import numpy as np
from numpy.fft import fft2, ifft2, fftfreq
import matplotlib.pyplot as plt
from datetime import datetime
import random

__version__ = "10.0.0"

# ============================================================
# CLI SETUP
# ============================================================
parser = argparse.ArgumentParser(description="iDNS v10 — Kolmogorov Flow (Corrected Forcing)")
parser.add_argument("--mode", choices=["conservative", "idns"], default="idns")
parser.add_argument("--Re", type=float, default=20000.0)
parser.add_argument("--t_end", type=float, default=100.0)
parser.add_argument("--seed", type=int, default=48)
parser.add_argument("--outdir", type=str, default="results")
parser.add_argument("--N", type=int, default=256)
# CCC Options
parser.add_argument("--ccc", choices=["off", "periodic", "adaptive"], default="off",
                    help="Cascade-Completion Clamping mode: off, periodic, or adaptive")
parser.add_argument("--ccc_interval", type=int, default=200,
                    help="Steps between CCC applications (periodic mode)")
parser.add_argument("--ccc_threshold", type=float, default=2.0,
                    help="Pile-up factor threshold for CCC (adaptive mode)")
parser.add_argument("--ccc_fraction", type=float, default=0.5,
                    help="Only clamp modes above this fraction of k_max (default 0.5)")
parser.add_argument("--ccc_exponent", type=float, default=-3.0,
                    help="Spectral envelope exponent (default -3 for 2D enstrophy cascade)")
args = parser.parse_args()

# ============================================================
# UNIQUE RUN ID — prevents file overwrites
# ============================================================
run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
print(f"\n🔑 Run ID: {run_id}")

# ============================================================
# CONFIGURATION
# ============================================================
N = args.N
L = 2 * np.pi
F = 1  # idns standard forcing (use 1.0 for C&K benchmark, 2.5 normally)
kf = 4   # idns standard wavenumber (use 4 for C&K benchmark, 8 normally)
mode = args.mode
Re = args.Re
seed = args.seed
outdir = os.path.join(args.outdir, f"v10_N{N}_Re{int(Re)}")
os.makedirs(outdir, exist_ok=True)

dt = 0.001 if mode == "conservative" else 0.006
nu = L / (kf * Re)
phi_prime = 1.0
print_interval = 100

# CCC configuration
ccc_mode = args.ccc
ccc_interval = args.ccc_interval
ccc_threshold = args.ccc_threshold
ccc_fraction = args.ccc_fraction
ccc_exponent = args.ccc_exponent

# ============================================================
# CONSOLE HEADER
# ============================================================
print("\n" + "=" * 80)
print(f" iDNS v{__version__} Kolmogorov | {mode.upper()} MODE | N={N} Re={Re:.1e}")
print(f" F={F}  kf={kf}  dt={dt:.3f}  t_end={args.t_end}  ν={nu:.6e}")
print(f" Temporal Lifting  : {'ON (RASI)' if mode=='idns' else 'OFF'}")
print(f" CCC Mode          : {ccc_mode.upper()}", end="")
if ccc_mode == "periodic":
    print(f" (every {ccc_interval} steps)")
elif ccc_mode == "adaptive":
    print(f" (threshold={ccc_threshold})")
else:
    print()
print(f" Run ID            : {run_id}")
print("=" * 80)

# ============================================================
# SPECTRAL SETUP
# ============================================================
# Wavenumbers: k = 2π * f where f = fftfreq(N, d=dx)
# With dx = L/N and L = 2π, this gives integer wavenumbers [-N/2, ..., N/2-1]
k1 = fftfreq(N, d=L / N) * 2 * np.pi
kx, ky = np.meshgrid(k1, k1, indexing="xy")  # xy indexing matches physical (x,y) layout
k2 = kx ** 2 + ky ** 2
k_mag = np.sqrt(k2)
k_max = np.max(k_mag)

# 2/3 dealiasing mask
k_dealias = (2.0 / 3.0) * (N / 2) * (2 * np.pi / L)
dealias_mask = k_mag <= k_dealias

# ============================================================
# CASCADE-COMPLETION CLAMPING (CCC) FUNCTIONS
# ============================================================
def cascade_completion_clamp(w_hat, kx, ky, kf, clamp_fraction=0.5, exponent=-3.0):
    """
    Cascade-Completion Clamping: enforce theoretical spectral envelope at cutoff.
    
    "We're just doing what physics would do with more modes."
    
    Clamps modes above clamp_fraction * k_max to theoretical k^exponent envelope.
    Uses inertial range (2*kf to 3*kf) as reference to avoid forcing contamination.
    """
    k_mag_local = np.sqrt(kx**2 + ky**2)
    k_max_local = np.max(k_mag_local)
    k_clamp_start = clamp_fraction * k_max_local
    
    # Reference energy from INERTIAL RANGE (2*kf to 3*kf), not forcing scale
    k_ref = 2.0 * kf
    inertial_mask = (k_mag_local >= 1.5 * kf) & (k_mag_local <= 3.0 * kf)
    
    if np.sum(inertial_mask) > 0:
        E_inertial = np.mean(np.abs(w_hat[inertial_mask])**2)
    else:
        forcing_mask = (k_mag_local >= kf * 0.8) & (k_mag_local <= kf * 1.2)
        E_inertial = np.mean(np.abs(w_hat[forcing_mask])**2) if np.sum(forcing_mask) > 0 else np.mean(np.abs(w_hat)**2)
        k_ref = kf
    
    # Theoretical envelope: E(k) = E_ref * (k/k_ref)^exponent
    with np.errstate(divide='ignore', invalid='ignore'):
        E_envelope = E_inertial * np.power(k_mag_local / k_ref, exponent)
        E_envelope = np.where(k_mag_local > 0, E_envelope, E_inertial)
    
    E_current = np.abs(w_hat)**2
    high_k_mask = k_mag_local >= k_clamp_start
    
    with np.errstate(divide='ignore', invalid='ignore'):
        correction = np.where(
            (E_current > 0) & high_k_mask,
            np.minimum(1.0, E_envelope / E_current),
            1.0
        )
    
    w_hat_clamped = w_hat * np.sqrt(correction)
    
    E_before = np.sum(np.abs(w_hat)**2)
    E_after = np.sum(np.abs(w_hat_clamped)**2)
    energy_removed = E_before - E_after
    pct_removed = 100 * energy_removed / (E_before + 1e-20)
    
    return w_hat_clamped, energy_removed, pct_removed


def ccc_diagnostic(w_hat, kx, ky, kf):
    """
    Diagnose pile-up at spectral cutoff.
    Returns pile_up_factor and whether clamping is needed.
    """
    k_mag_local = np.sqrt(kx**2 + ky**2)
    k_max_local = np.max(k_mag_local)
    
    k_ref = 2.0 * kf
    inertial_mask = (k_mag_local >= 1.5 * kf) & (k_mag_local <= 3.0 * kf)
    if np.sum(inertial_mask) > 0:
        E_inertial = np.mean(np.abs(w_hat[inertial_mask])**2)
    else:
        return 1.0, False
    
    high_k_mask = (k_mag_local >= 0.5 * k_max_local)
    if np.sum(high_k_mask) > 0:
        E_highk = np.mean(np.abs(w_hat[high_k_mask])**2)
    else:
        return 1.0, False
    
    k_high_avg = 0.75 * k_max_local
    theoretical_ratio = (k_high_avg / k_ref) ** (-3)
    measured_ratio = E_highk / (E_inertial + 1e-20)
    pile_up_factor = measured_ratio / (theoretical_ratio + 1e-20)
    
    return pile_up_factor, pile_up_factor > ccc_threshold


def check_resolution(Re, N, kf=8):
    """Check if grid resolves the enstrophy dissipation scale."""
    k_d = kf * np.sqrt(Re / 100)
    k_max = N / 2
    r = k_d / k_max
    
    if r <= 1.0:
        status = "RESOLVED"
    elif r <= 2.0:
        status = "MARGINAL"
    elif r <= 3.0:
        status = "UNDER-RESOLVED"
    else:
        status = "SEVERELY UNDER-RESOLVED"
    
    return r, status, k_d


# Check resolution at startup
r_ratio, res_status, k_d = check_resolution(Re, N, kf)
print(f"\n📐 Resolution check: k_d/k_max = {r_ratio:.2f} ({res_status})")
print(f"   k_d = {k_d:.1f}, k_max = {N//2}")
if ccc_mode == "off" and r_ratio > 1.5:
    print(f"   ⚠️  Consider enabling CCC (--ccc periodic or --ccc adaptive)")
print()


# ============================================================
# CORE SPECTRAL FUNCTIONS (v9 CORRECTED)
# ============================================================

def project_uv_from_w(w_hat, kx, ky, k2):
    """
    Recover velocity (u, v) from vorticity ω via stream function ψ.
    
    CORRECTED v9 FORMULATION:
        ω = -Δψ  →  in Fourier: -k²ψ̂ = -ω̂  →  ψ̂ = +ω̂/k²
        u = ∂ψ/∂y  →  û = ik_y ψ̂
        v = -∂ψ/∂x →  v̂ = -ik_x ψ̂
    
    v8 BUG: Had ψ̂ = -ω̂/k² which negated the entire velocity field.
    """
    psi_hat = np.zeros_like(w_hat, dtype=complex)
    
    # CORRECTED: positive sign (was negative in v8)
    psi_hat[k2 != 0] = w_hat[k2 != 0] / k2[k2 != 0]
    
    # Velocity from stream function
    u_hat = 1j * ky * psi_hat
    v_hat = -1j * kx * psi_hat
    
    return np.real(ifft2(u_hat)), np.real(ifft2(v_hat))


def spectral_gradient(f_hat, kx, ky):
    """
    Compute spectral derivatives ∂f/∂x and ∂f/∂y.
    
    v9: Replaces np.gradient (finite difference) with proper spectral derivatives.
    """
    df_dx = np.real(ifft2(1j * kx * f_hat))
    df_dy = np.real(ifft2(1j * ky * f_hat))
    return df_dx, df_dy


def rhs_hat(w_hat, kx, ky, k2, nu, curlf_hat, dealias_mask):
    """
    Compute RHS of vorticity equation: ∂ω/∂t = -(u·∇)ω + ν∆ω + curl(f)
    
    v9 CORRECTIONS:
        1. Uses corrected project_uv_from_w (positive stream function)
        2. Uses spectral derivatives instead of np.gradient
    """
    # Get velocity from vorticity (CORRECTED)
    u, v = project_uv_from_w(w_hat, kx, ky, k2)
    
    # Spectral derivatives of vorticity (CORRECTED - was using np.gradient)
    w_x, w_y = spectral_gradient(w_hat, kx, ky)
    
    # Advection term: (u·∇)ω = u * ∂ω/∂x + v * ∂ω/∂y
    advection = u * w_x + v * w_y
    advection_hat = fft2(advection)
    
    # Viscous term: ν∆ω → -νk²ω̂
    visc_hat = -nu * k2 * w_hat
    
    # Full RHS: -advection + viscosity + forcing
    rhs = -advection_hat + visc_hat + curlf_hat
    
    # Apply dealiasing
    rhs *= dealias_mask
    
    return rhs


def curl_forcing_hat(F, kf, N, L):
    """
    Kolmogorov forcing: f = F·sin(kf·y)·x̂
    
    v10 CORRECTED:
        curl(f) = ∂f_x/∂y - ∂f_y/∂x = F·kf·cos(kf·y) - 0
        
        BUT for vorticity equation ∂ω/∂t = ... + curl(f), we need the 
        z-component of curl in 2D, which comes from:
        
        ω = ∂v/∂x - ∂u/∂y
        
        For forcing on u only: f_u = F·sin(kf·y), f_v = 0
        The forcing contribution to ∂ω/∂t is:
        ∂/∂t(∂v/∂x - ∂u/∂y) includes -∂f_u/∂y = -F·kf·cos(kf·y)
        
        C&K eq (2.6) confirms: vorticity forcing = -n·cos(ny)
    """
    x = np.linspace(0, L, N, endpoint=False)
    y = np.linspace(0, L, N, endpoint=False)
    X, Y = np.meshgrid(x, y, indexing='xy')
    
    # CORRECTED v10: negative sign to match C&K eq (2.6)
    # curl(f) for vorticity forcing = -F·kf·cos(kf·y)
    curlf = -F * kf * np.cos(kf * Y)
    
    return fft2(curlf)


def spectral_restart(w_hat, kx, ky, k2, a=5.0, p=3.0, cutoff=0.95, t=None, log_events=None):
    """Emergency spectral filtering for NaN recovery."""
    E_before = np.sum(np.abs(w_hat)**2)
    k_mag_local = np.sqrt(k2)
    k_max_local = np.max(k_mag_local)
    zeta = 1.0 / (1.0 + np.exp(a * ((k_mag_local / k_max_local - cutoff) ** p)))
    w_hat_new = np.nan_to_num(w_hat * zeta, nan=0.0, posinf=0.0, neginf=0.0)
    E_after = np.sum(np.abs(w_hat_new)**2)
    if log_events is not None:
        log_events.append({
            "t": float(t) if t is not None else 0.0,
            "event": "spectral_restart",
            "energy_drop": float(E_before - E_after)
        })
    print(f"⚠️  Spectral restart at t={t if t else 0.0:.3f} | ΔE={E_before-E_after:.3e}")
    return w_hat_new


def shell_energy(w_hat, k2, N):
    """Compute energy in spectral shells."""
    k_mag_local = np.sqrt(k2)
    
    if N >= 256:
        k_bins = [0, 4, 8, 16, 32, 64, 128]
    elif N >= 128:
        k_bins = [0, 4, 8, 16, 32, 64]
    else:
        k_bins = [0, 4, 8, 16, 32]
    
    E_shells = []
    for i in range(len(k_bins) - 1):
        mask_shell = (k_mag_local >= k_bins[i]) & (k_mag_local < k_bins[i + 1])
        E_shells.append(np.sum(np.abs(w_hat[mask_shell])**2))
    
    return E_shells, k_bins


# ============================================================
# INITIAL CONDITIONS
# ============================================================
curlf_hat = curl_forcing_hat(F, kf, N, L)

# Laminar Kolmogorov profile + perturbation
x = np.linspace(0, L, N, endpoint=False)
y = np.linspace(0, L, N, endpoint=False)
X, Y = np.meshgrid(x, y, indexing='xy')  # xy indexing: X varies along columns, Y along rows

# Initial vorticity: laminar solution ω = -kf * cos(kf * y) plus noise
w0 = -kf * np.cos(kf * Y)
rng = np.random.default_rng(seed)
perturb_amp = 0.02 * kf
w0 += perturb_amp * rng.standard_normal((N, N))
print(f"✨ Initialization complete: perturbation amplitude={perturb_amp:.4f}, seed={seed}")

w_hat = fft2(w0)
w_hat *= dealias_mask

# ============================================================
# MAIN LOOP
# ============================================================
start_time = time.time()
t = tau = 0.0
step = 0
rows, spectral_snapshots, restart_log = [], [], []
negotiation_log = []
ccc_log = []
shell_energy_log = []
cz_count = 0
ccc_count = 0
max_restarts = 5
consecutive_restarts = 0
max_omega_achieved = 0.0
total_energy_clamped = 0.0

print("\n🚀 Simulation started...\n")

while t < args.t_end:
    step += 1
    w_hat_prev = w_hat.copy()

    # ========================================================
    # TEMPORAL LIFTING CONTROLLER (RASI)
    # ========================================================
    if mode == "idns":
        w = np.real(ifft2(w_hat))
        w_inf = float(np.max(np.abs(w)))
        spec_amp = float(np.sqrt(np.sum(k2 * np.abs(w_hat)**2)))
        
        score = 0.6 * min(1.0, w_inf / 20.0) + 0.4 * min(1.0, spec_amp / 10.0)
        alpha, center, sigma_min = 4.5, 0.55, 0.03
        sigma = sigma_min + (1.0 - sigma_min) / (1.0 + np.exp(alpha * (score - center)))
        phi_prime = 1.0 / sigma
        
        # Emergency override for ultra-high vorticity
        if w_inf > 5000:
            phi_prime = min(30.0, 1.5 * phi_prime)
    else:
        phi_prime = 1.0
        w = np.real(ifft2(w_hat))
        w_inf = float(np.max(np.abs(w)))

    # ========================================================
    # RK4 Integration (lifted time formulation)
    # ========================================================
    fac = 1.0 / max(1.0, float(phi_prime))
    
    k1 = fac * rhs_hat(w_hat, kx, ky, k2, nu, curlf_hat, dealias_mask)
    k2_ = fac * rhs_hat(w_hat + 0.5 * dt * k1, kx, ky, k2, nu, curlf_hat, dealias_mask)
    k3 = fac * rhs_hat(w_hat + 0.5 * dt * k2_, kx, ky, k2, nu, curlf_hat, dealias_mask)
    k4 = fac * rhs_hat(w_hat + dt * k3, kx, ky, k2, nu, curlf_hat, dealias_mask)
    
    w_hat += (dt / 6.0) * (k1 + 2 * k2_ + 2 * k3 + k4)
    w_hat *= dealias_mask

    # ========================================================
    # CASCADE-COMPLETION CLAMPING (CCC)
    # ========================================================
    ccc_applied = False
    
    if ccc_mode == "periodic" and step % ccc_interval == 0:
        w_hat, E_removed, pct_removed = cascade_completion_clamp(
            w_hat, kx, ky, kf,
            clamp_fraction=ccc_fraction,
            exponent=ccc_exponent
        )
        if pct_removed > 0.01:
            ccc_count += 1
            total_energy_clamped += E_removed
            ccc_applied = True
            ccc_log.append({
                "t": float(t),
                "step": int(step),
                "mode": "periodic",
                "energy_removed": float(E_removed),
                "pct_removed": float(pct_removed)
            })
            if step % (print_interval * 5) == 0:
                print(f"   🔧 CCC (periodic): removed {pct_removed:.3f}% energy")
    
    elif ccc_mode == "adaptive":
        pile_up, needs_clamp = ccc_diagnostic(w_hat, kx, ky, kf)
        if needs_clamp:
            w_hat, E_removed, pct_removed = cascade_completion_clamp(
                w_hat, kx, ky, kf,
                clamp_fraction=ccc_fraction,
                exponent=ccc_exponent
            )
            ccc_count += 1
            total_energy_clamped += E_removed
            ccc_applied = True
            ccc_log.append({
                "t": float(t),
                "step": int(step),
                "mode": "adaptive",
                "pile_up_factor": float(pile_up),
                "energy_removed": float(E_removed),
                "pct_removed": float(pct_removed)
            })
            if step % print_interval == 0:
                print(f"   🔧 CCC (adaptive): pile-up={pile_up:.1f}x, removed {pct_removed:.3f}%")

    # ========================================================
    # SINGULARITY HANDLING
    # ========================================================
    if not np.isfinite(np.max(np.abs(w_hat))):
        consecutive_restarts += 1
        print(f"\n🚨 SINGULARITY at t={t:.6f}, step={step}")
        
        if consecutive_restarts > max_restarts:
            print(f"\n❌ Cannot recover after {max_restarts} attempts")
            break
        
        # Progressive spectral negotiation
        filter_strengths = np.logspace(-4, -0.1, 20)
        negotiated = False
        
        for i, strength in enumerate(filter_strengths):
            k_cutoff = k_max * (1.0 - strength)
            zeta = np.exp(-(k_mag / k_cutoff)**6)
            zeta[0, 0] = 1.0
            w_hat_filtered = w_hat_prev * zeta
            
            try:
                test_rhs = rhs_hat(w_hat_filtered, kx, ky, k2, nu, curlf_hat, dealias_mask)
                if np.isfinite(np.max(np.abs(test_rhs))):
                    w_hat = w_hat_filtered
                    print(f"   ✅ Negotiated with strength={strength:.5f}")
                    negotiation_log.append({"t": float(t), "strength": float(strength)})
                    cz_count += 1
                    negotiated = True
                    break
            except:
                continue
        
        if not negotiated:
            w_hat = spectral_restart(w_hat_prev, kx, ky, k2, t=t, log_events=restart_log)
            step -= 1
            t -= dt
            tau -= dt
            continue
    else:
        consecutive_restarts = 0

    # ========================================================
    # Diagnostics
    # ========================================================
    omega = np.real(ifft2(w_hat))
    omega_max = float(np.max(np.abs(omega)))
    max_omega_achieved = max(max_omega_achieved, omega_max)
    
    u, v = project_uv_from_w(w_hat, kx, ky, k2)
    E_phys = 0.5 * np.mean(u**2 + v**2)
    Z_phys = 0.5 * np.mean(omega**2)
    dissipation = 2 * nu * Z_phys
    power_in = F * np.mean(u * np.sin(kf * Y))

    rows.append([step, t, tau, omega_max, E_phys, Z_phys, power_in, dissipation, phi_prime,
                 1 if ccc_applied else 0])

    # Track shell energies periodically
    if step % 100 == 0:
        E_shells, k_bins = shell_energy(w_hat, k2, N)
        shell_energy_log.append({"t": float(t), "step": int(step), "shells": E_shells})

    # ========================================================
    # Console output
    # ========================================================
    if step % print_interval == 0 or step == 1:
        # Compute R_epsilon for energy budget validation
        R_eps = dissipation / (power_in + 1e-20) if power_in > 0 else float('inf')
        ccc_status = f"CCC={ccc_count}" if ccc_mode != "off" else ""
        print(f"{step:7d} | "
              f"t={t:8.3f} | "
              f"E={E_phys:10.6f} | "
              f"Z={Z_phys:10.6f} | "
              f"P_in={power_in:+.4f} | "
              f"R_ε={R_eps:5.3f} | "
              f"‖ω‖∞={omega_max:8.3f} | "
              f"φ′={phi_prime:5.2f} | "
              f"{ccc_status}")

    # ========================================================
    # Advance clocks
    # ========================================================
    if mode == "idns":
        t += phi_prime * dt
        tau += dt
    else:
        t += dt
        tau += dt

# ============================================================
# WRAP-UP
# ============================================================
runtime = time.time() - start_time
print("\n\n✅ Simulation complete.\n")
print("-" * 70)
print(f" Version          : iDNS v{__version__}")
print(f" Mode             : {mode}")
print(f" Grid N           : {N}")
print(f" Re               : {Re:.1e}")
print(f" CCC Mode         : {ccc_mode}")
print(f" Steps completed  : {step}")
print(f" Final time t     : {t:.6f}  (target {args.t_end})")
print(f" Wall-clock (s)   : {runtime:.2f}")
print(f" CCC applications : {ccc_count}")
print(f" Total E clamped  : {total_energy_clamped:.3e}")
print(f" Max ω achieved   : {max_omega_achieved:.2e}")
print(f" Status           : {'PASS' if t > 0.9 * args.t_end else 'INCOMPLETE'}")
print("-" * 70)

# Shell energy final state
E_shells, k_bins = shell_energy(w_hat, k2, N)
print("\n📊 SHELL ENERGY (Final):")
for i in range(len(E_shells)):
    print(f"   [{k_bins[i]:3d},{k_bins[i+1]:3d}) : {E_shells[i]:.6e}")

# ============================================================
# SAVE OUTPUTS — with unique run_id
# ============================================================
base_name = f"idns_v10_{mode}_N{N}_Re{int(Re)}_ccc{ccc_mode}_{run_id}"

# CSV with all data
suffix = run_id[-4:] + ''.join(random.choices('abcdef0123456789', k=2))
csv_file = f"{outdir}/{base_name}_{suffix}.csv"
with open(csv_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["step", "t", "tau", "w_inf", "E", "Z", "P_in", "epsilon", "phi_prime", "ccc_applied"])
    writer.writerows(rows)
print(f"\n📄 Data saved → {csv_file}")

# JSON telemetry
telemetry = {
    "run_id": run_id,
    "version": __version__,
    "configuration": {
        "Re": float(Re), "N": int(N), "F": float(F), "kf": int(kf),
        "nu": float(nu), "dt": float(dt), "mode": mode,
        "ccc_mode": ccc_mode, "ccc_interval": ccc_interval,
        "ccc_threshold": ccc_threshold, "ccc_fraction": ccc_fraction,
        "ccc_exponent": ccc_exponent
    },
    "resolution": {
        "k_d_over_k_max": float(r_ratio),
        "status": res_status,
        "k_d": float(k_d),
        "k_max": int(N // 2)
    },
    "ccc_events": ccc_log,
    "negotiation_events": negotiation_log,
    "shell_energy_history": shell_energy_log,
    "final_state": {
        "t_final": float(t),
        "steps_total": int(step),
        "ccc_total": ccc_count,
        "total_energy_clamped": float(total_energy_clamped),
        "max_omega": float(max_omega_achieved),
        "runtime_seconds": float(runtime),
        "status": "PASS" if t > 0.9 * args.t_end else "INCOMPLETE",
        "final_shell_energies": E_shells,
        "shell_bins": k_bins
    }
}

json_file = f"{outdir}/{base_name}.json"
with open(json_file, "w") as f:
    json.dump(telemetry, f, indent=2)
print(f"📡 Telemetry saved → {json_file}")

print(f"\niDNS v{__version__} COMPLETE | Run ID: {run_id}\n")
