#!/usr/bin/env python3
"""
extract_spectrum.py — Extract energy spectrum E(k) from iDNS velocity snapshots
Usage: python extract_spectrum.py field_t4.5.npz --output spectrum_t4.5.csv

Outputs a small CSV with columns: k, E_k
"""

import numpy as np
import argparse

def compute_spectrum(ux, uy, uz):
    """Compute 1D energy spectrum E(k) from 3D velocity field."""
    N = ux.shape[0]
    
    # 3D FFT of velocity components
    ux_hat = np.fft.fftn(ux) / N**3
    uy_hat = np.fft.fftn(uy) / N**3
    uz_hat = np.fft.fftn(uz) / N**3
    
    # Energy in Fourier space: 0.5 * |u_hat|^2
    E_hat = 0.5 * (np.abs(ux_hat)**2 + np.abs(uy_hat)**2 + np.abs(uz_hat)**2)
    
    # Wavenumber grid
    kx = np.fft.fftfreq(N, d=1/N)
    ky = np.fft.fftfreq(N, d=1/N)
    kz = np.fft.fftfreq(N, d=1/N)
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
    K = np.sqrt(KX**2 + KY**2 + KZ**2)
    
    # Bin by wavenumber magnitude
    k_max = int(N // 2)
    k_bins = np.arange(0.5, k_max + 0.5, 1)
    E_k = np.zeros(k_max)
    
    for i in range(k_max):
        mask = (K >= k_bins[i] - 0.5) & (K < k_bins[i] + 0.5)
        E_k[i] = np.sum(E_hat[mask])
    
    k = np.arange(1, k_max + 1)
    return k, E_k

def main():
    parser = argparse.ArgumentParser(description='Extract energy spectrum from velocity snapshot')
    parser.add_argument('npz_file', help='Input .npz file with ux, uy, uz fields')
    parser.add_argument('--output', '-o', default='spectrum.csv', help='Output CSV file')
    args = parser.parse_args()
    
    print(f"Loading {args.npz_file}...")
    data = np.load(args.npz_file)
    
    ux = data['ux']
    uy = data['uy']
    uz = data['uz']
    t = float(data['t']) if 't' in data else 0.0
    
    print(f"  Grid: {ux.shape[0]}^3, t = {t:.2f}")
    print("Computing spectrum...")
    
    k, E_k = compute_spectrum(ux, uy, uz)
    
    # Save to CSV
    with open(args.output, 'w') as f:
        f.write('k,E_k\n')
        for ki, Ei in zip(k, E_k):
            f.write(f'{ki},{Ei}\n')
    
    print(f"Saved: {args.output} ({len(k)} wavenumbers)")
    print(f"  k range: 1 to {k[-1]}")
    print(f"  E(k) range: {E_k[E_k>0].min():.2e} to {E_k.max():.2e}")

if __name__ == '__main__':
    main()
