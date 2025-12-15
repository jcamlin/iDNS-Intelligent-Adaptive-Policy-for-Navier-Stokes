#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
idns_visualizer.py — Taylor-Green Vortex Visualization Utilities
-------------------------------------------------------------------------
Creates 2D slices and 3D isosurface renders from iDNS snapshots.

Version: 10.0.0
Date: 14 December 2025
Author: Jeffrey Camlin
ORCID: 0000-0002-5740-4204

Code and data: doi:10.5281/zenodo.17730873
GitHub: RED-DAWN-AI/iDNS

Copyright (c) 2025 Red Dawn Academic Press & AI Lab, Milwaukee WI
License: CC BY-NC 4.0 (Attribution-NonCommercial)
https://creativecommons.org/licenses/by-nc/4.0/
"""

import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from skimage.measure import marching_cubes
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

__version__ = "10.0.0"


def load_and_compute_vorticity(npz_path):
    """Load snapshot and compute vorticity field."""
    print(f"Loading: {npz_path}")
    data = np.load(npz_path)
    
    ux = data['ux']
    uy = data['uy']
    uz = data['uz']
    t = float(data['t']) if data['t'].shape == () else float(data['t'][0])
    
    N = ux.shape[0]
    L = 2 * np.pi
    print(f"  N = {N}, t = {t:.4f}")
    
    # Compute vorticity via spectral derivatives
    dx = L / N
    k = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    KX, KY, KZ = np.meshgrid(k, k, k, indexing='ij')
    
    uh = np.fft.fftn(ux)
    vh = np.fft.fftn(uy)
    wh = np.fft.fftn(uz)
    
    # omega = curl(u)
    omega_x = np.fft.ifftn(1j*KY*wh - 1j*KZ*vh).real
    omega_y = np.fft.ifftn(1j*KZ*uh - 1j*KX*wh).real
    omega_z = np.fft.ifftn(1j*KX*vh - 1j*KY*uh).real
    
    omega_mag = np.sqrt(omega_x**2 + omega_y**2 + omega_z**2)
    
    print(f"  ||omega||_max = {omega_mag.max():.3f}")
    
    return ux, uy, uz, omega_mag, t, N


def plot_composite_4panel(ux, uy, uz, omega_mag, t, N, Re="10^5", output_path="1-TG_composite_4panel.png"):
    """XZ slice four ways: velocity, vorticity, and composites with contours."""
    
    vel_mag = np.sqrt(ux**2 + uy**2 + uz**2)
    mid = N // 2
    extent = [0, 2*np.pi, 0, 2*np.pi]
    
    omega_slice = omega_mag[:, mid, :]
    vel_slice = vel_mag[:, mid, :]
    
    vmax_vel = vel_mag.max()
    vmax_omega = omega_mag.max() * 0.8
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    im0 = axes[0,0].imshow(vel_slice, origin='lower', cmap='viridis',
                           extent=extent, vmin=0, vmax=vmax_vel)
    axes[0,0].set_title('Velocity |u| - wave structure', fontsize=12)
    axes[0,0].set_ylabel('z')
    axes[0,0].set_aspect('equal')
    plt.colorbar(im0, ax=axes[0,0], shrink=0.8, label='|u|')
    
    im1 = axes[0,1].imshow(omega_slice, origin='lower', cmap='hot',
                           extent=extent, vmin=0, vmax=vmax_omega)
    axes[0,1].set_title('Vorticity |omega| - funnel cores', fontsize=12)
    axes[0,1].set_aspect('equal')
    plt.colorbar(im1, ax=axes[0,1], shrink=0.8, label='|omega|')
    
    im2 = axes[1,0].imshow(vel_slice, origin='lower', cmap='viridis',
                           extent=extent, vmin=0, vmax=vmax_vel)
    levels = np.linspace(omega_mag.max()*0.15, omega_mag.max()*0.8, 8)
    axes[1,0].contour(omega_slice, levels=levels, colors='white', 
                      linewidths=1.0, extent=extent)
    axes[1,0].set_title('Velocity + omega contours', fontsize=12)
    axes[1,0].set_xlabel('x')
    axes[1,0].set_ylabel('z')
    axes[1,0].set_aspect('equal')
    plt.colorbar(im2, ax=axes[1,0], shrink=0.8, label='|u|')
    
    im3 = axes[1,1].imshow(vel_slice, origin='lower', cmap='inferno',
                           extent=extent, vmin=0, vmax=vmax_vel)
    axes[1,1].contour(omega_slice, levels=levels, colors='cyan', 
                      linewidths=1.0, extent=extent)
    axes[1,1].set_title('Velocity (hot) + omega contours', fontsize=12)
    axes[1,1].set_xlabel('x')
    axes[1,1].set_aspect('equal')
    plt.colorbar(im3, ax=axes[1,1], shrink=0.8, label='|u|')
    
    fig.suptitle(f'Taylor-Green XZ Slice (y = pi) - t = {t:.4f}, N = {N}^3, Re = {Re}',
                 fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()


def plot_funnel_detail(omega_mag, t, N, Re="10^5", output_path="2-TG_funnel_detail.png"):
    """Side-by-side XZ slice: vorticity with and without contours."""
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    mid = N // 2
    extent = [0, 2*np.pi, 0, 2*np.pi]
    vmax = omega_mag.max() * 0.8
    
    omega_slice = omega_mag[:, mid, :]
    
    # Vorticity only (left)
    im0 = axes[0].imshow(omega_slice, origin='lower', cmap='hot',
                          extent=extent, vmin=0, vmax=vmax)
    axes[0].set_title('Vorticity |omega|', fontsize=12)
    axes[0].set_xlabel('x')
    axes[0].set_ylabel('z')
    axes[0].set_aspect('equal')
    plt.colorbar(im0, ax=axes[0], shrink=0.8, label='|omega|')
    
    # Vorticity with contours (right)
    im1 = axes[1].imshow(omega_slice, origin='lower', cmap='inferno',
                          extent=extent, vmin=0, vmax=vmax)
    levels = np.linspace(omega_mag.max()*0.2, omega_mag.max()*0.8, 6)
    axes[1].contour(omega_slice, levels=levels, colors='cyan', linewidths=1.0,
                    extent=extent, alpha=0.8)
    axes[1].set_title('Vorticity |omega| + Contours', fontsize=12)
    axes[1].set_xlabel('x')
    axes[1].set_ylabel('z')
    axes[1].set_aspect('equal')
    plt.colorbar(im1, ax=axes[1], shrink=0.8, label='|omega|')
    
    fig.suptitle(f'Taylor-Green Funnel Structure - XZ plane (y = pi)\n'
                 f't = {t:.4f}, N = {N}^3, Re = {Re}', 
                 fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()


def plot_composite(ux, uy, uz, omega_mag, t, N, Re="10^5", output_path="3-TG_composite.png"):
    """Velocity magnitude (fill) + Vorticity contours (lines) - 3 panels."""
    
    vel_mag = np.sqrt(ux**2 + uy**2 + uz**2)
    
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    mid = N // 2
    extent = [0, 2*np.pi, 0, 2*np.pi]
    vmax = vel_mag.max()
    
    slices = [
        (omega_mag[:, :, mid].T, vel_mag[:, :, mid].T, 'XY plane (z = pi)', 'x', 'y'),
        (omega_mag[:, mid, :], vel_mag[:, mid, :], 'XZ plane (y = pi)', 'x', 'z'),
        (omega_mag[mid, :, :].T, vel_mag[mid, :, :].T, 'YZ plane (x = pi)', 'y', 'z'),
    ]
    
    for ax, (omega_slice, vel_slice, title, xlab, ylab) in zip(axes, slices):
        im = ax.imshow(vel_slice, origin='lower', cmap='inferno',
                       extent=extent, vmin=0, vmax=vmax, alpha=0.9)
        
        levels = np.linspace(omega_mag.max()*0.1, omega_mag.max()*0.8, 8)
        ax.contour(omega_slice, levels=levels, colors='white', linewidths=0.8,
                   extent=extent, alpha=0.7)
        
        ax.set_title(title, fontsize=12)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_aspect('equal')
        plt.colorbar(im, ax=ax, shrink=0.8, label='|u|')
    
    fig.suptitle(f'Taylor-Green Vortex - t = {t:.4f}, N = {N}^3, Re = {Re}\n'
                 f'Velocity |u| (fill) + Vorticity |omega| (contours)', 
                 fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()


def plot_velocity_slices(ux, uy, uz, t, N, Re="10^5", output_path="4-TG_velocity_slices.png"):
    """Plot velocity magnitude slices, each with own colorbar."""
    
    vel_mag = np.sqrt(ux**2 + uy**2 + uz**2)
    
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    mid = N // 2
    vmax = vel_mag.max()
    extent = [0, 2*np.pi, 0, 2*np.pi]
    
    im0 = axes[0].imshow(vel_mag[:, :, mid].T, origin='lower', cmap='viridis',
                          extent=extent, vmin=0, vmax=vmax)
    axes[0].set_title('XY plane (z = pi)')
    axes[0].set_xlabel('x')
    axes[0].set_ylabel('y')
    plt.colorbar(im0, ax=axes[0], shrink=0.8, label='|u|')
    
    im1 = axes[1].imshow(vel_mag[:, mid, :], origin='lower', cmap='viridis',
                          extent=extent, vmin=0, vmax=vmax)
    axes[1].set_title('XZ plane (y = pi)')
    axes[1].set_xlabel('x')
    axes[1].set_ylabel('z')
    plt.colorbar(im1, ax=axes[1], shrink=0.8, label='|u|')
    
    im2 = axes[2].imshow(vel_mag[mid, :, :].T, origin='lower', cmap='viridis',
                          extent=extent, vmin=0, vmax=vmax)
    axes[2].set_title('YZ plane (x = pi)')
    axes[2].set_xlabel('y')
    axes[2].set_ylabel('z')
    plt.colorbar(im2, ax=axes[2], shrink=0.8, label='|u|')
    
    fig.suptitle(f'Taylor-Green Velocity |u| - t = {t:.4f}, N = {N}^3, Re = {Re}', 
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"Saved: {output_path}")
    plt.close()


def plot_vorticity_slices(omega_mag, t, N, Re="10^5", output_path="5-TG_vorticity_slices.png"):
    """Plot vorticity magnitude at multiple z-slices."""
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    z_indices = [0, N//6, N//4, N//3, N//2, 2*N//3]
    vmax = omega_mag.max() * 0.8
    extent = [0, 2*np.pi, 0, 2*np.pi]
    
    for ax, zi in zip(axes.flat, z_indices):
        slice_data = omega_mag[:, :, zi]
        im = ax.imshow(slice_data.T, origin='lower', cmap='inferno',
                       extent=extent, vmin=0, vmax=vmax)
        ax.set_title(f'z = {zi * 2*np.pi / N:.2f}', fontsize=12)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        plt.colorbar(im, ax=ax, shrink=0.8, label='|omega|')
    
    fig.suptitle(f'Taylor-Green Vorticity |omega| at t = {t:.4f}, N = {N}^3, Re = {Re}', 
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    print(f"Saved: {output_path}")
    plt.close()


def plot_vorticity_3panel(omega_mag, t, N, Re="10^5", output_path="6-TG_vorticity_3panel.png"):
    """Plot three orthogonal slices through center, each with own colorbar."""
    
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    mid = N // 2
    vmax = omega_mag.max() * 0.8
    extent = [0, 2*np.pi, 0, 2*np.pi]
    
    # XY slice (z = mid)
    im0 = axes[0].imshow(omega_mag[:, :, mid].T, origin='lower', cmap='inferno',
                          extent=extent, vmin=0, vmax=vmax)
    axes[0].set_title('XY plane (z = pi)', fontsize=12)
    axes[0].set_xlabel('x')
    axes[0].set_ylabel('y')
    plt.colorbar(im0, ax=axes[0], shrink=0.8, label='|omega|')
    
    # XZ slice (y = mid)
    im1 = axes[1].imshow(omega_mag[:, mid, :], origin='lower', cmap='inferno',
                          extent=extent, vmin=0, vmax=vmax)
    axes[1].set_title('XZ plane (y = pi)', fontsize=12)
    axes[1].set_xlabel('x')
    axes[1].set_ylabel('z')
    plt.colorbar(im1, ax=axes[1], shrink=0.8, label='|omega|')
    
    # YZ slice (x = mid)
    im2 = axes[2].imshow(omega_mag[mid, :, :].T, origin='lower', cmap='inferno',
                          extent=extent, vmin=0, vmax=vmax)
    axes[2].set_title('YZ plane (x = pi)', fontsize=12)
    axes[2].set_xlabel('y')
    axes[2].set_ylabel('z')
    plt.colorbar(im2, ax=axes[2], shrink=0.8, label='|omega|')
    
    fig.suptitle(f'Taylor-Green Vorticity |omega| - t = {t:.4f}, N = {N}^3, Re = {Re}', 
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    print(f"Saved: {output_path}")
    plt.close()


def plot_isosurface_colored(omega_mag, t, N, Re="10^5", level=None, output_path="7-TG_isosurface_colored.png"):
    """3D isosurface colored by local vorticity value."""
    
    if level is None:
        level = omega_mag.max() * 0.25
    
    print(f"Extracting isosurface at |omega| = {level:.2f}...")
    verts, faces, normals, values = marching_cubes(omega_mag, level=level)
    verts = verts * (2 * np.pi / N)
    
    # Get vorticity values at face centers for coloring
    face_centers = verts[faces].mean(axis=1)
    face_indices = (face_centers * N / (2 * np.pi)).astype(int)
    face_indices = np.clip(face_indices, 0, N-1)
    face_values = omega_mag[face_indices[:, 0], face_indices[:, 1], face_indices[:, 2]]
    
    # Normalize for colormap
    face_colors = (face_values - face_values.min()) / (face_values.max() - face_values.min() + 1e-10)
    
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    mesh = Poly3DCollection(verts[faces], alpha=0.8)
    cmap = plt.cm.viridis
    mesh.set_facecolor(cmap(face_colors))
    mesh.set_edgecolor('none')
    ax.add_collection3d(mesh)
    
    ax.set_xlim(0, 2*np.pi)
    ax.set_ylim(0, 2*np.pi)
    ax.set_zlim(0, 2*np.pi)
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('y', fontsize=12)
    ax.set_zlabel('z', fontsize=12)
    ax.set_title(f'Taylor-Green Vortex Tubes (colored by |omega|)\n'
                 f't = {t:.4f}, N = {N}^3, Re = {Re}', fontsize=14, fontweight='bold')
    
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(face_values.min(), face_values.max()))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.1)
    cbar.set_label('|omega|', fontsize=11)
    
    ax.view_init(elev=25, azim=45)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=f"iDNS Visualizer v{__version__}")
    parser.add_argument("npz_path", type=str, help="Path to .npz snapshot file")
    parser.add_argument("--Re", type=str, default="10^5", help="Reynolds number for labels")
    parser.add_argument("--iso", action="store_true", help="Generate 3D isosurface (slower)")
    parser.add_argument("--iso_level", type=float, default=None, help="Isosurface level (default: 25%% of max)")
    args = parser.parse_args()
    
    # Load and compute
    ux, uy, uz, omega_mag, t, N = load_and_compute_vorticity(args.npz_path)
    
    print(f"\nGenerating visualizations...\n")
    
    # Generate in order
    plot_composite_4panel(ux, uy, uz, omega_mag, t, N, args.Re)      # 1
    plot_funnel_detail(omega_mag, t, N, args.Re)                     # 2
    plot_composite(ux, uy, uz, omega_mag, t, N, args.Re)             # 3
    plot_velocity_slices(ux, uy, uz, t, N, args.Re)                  # 4
    plot_vorticity_slices(omega_mag, t, N, args.Re)                  # 5
    plot_vorticity_3panel(omega_mag, t, N, args.Re)                  # 6
    
    if args.iso:
        plot_isosurface_colored(omega_mag, t, N, args.Re, args.iso_level)  # 7
    
    print(f"\niDNS Visualizer v{__version__} complete.")


 #Output order:**

#1-TG_composite_4panel.png
#2-TG_funnel_detail.png
#3-TG_composite.png
#4-TG_velocity_slices.png
#5-TG_vorticity_slices.png
#6-TG_vorticity_3panel.png
#7-TG_isosurface_colored.png  (with --iso)
