#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
npz_inspect.py — Quick Inspection of iDNS Snapshot Files
-------------------------------------------------------------------------
Displays contents, shapes, dtypes, and value ranges of .npz velocity fields.

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

__version__ = "10.0.0"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=f"npz_inspect v{__version__} — Inspect iDNS snapshot contents")
    parser.add_argument("npz_path", type=str, help="Path to .npz file")
    args = parser.parse_args()
    
    data = np.load(args.npz_path)
    print(f"File: {args.npz_path}")
    print(f"Keys: {list(data.keys())}")
    
    for k in data.keys():
        arr = data[k]
        if hasattr(arr, 'shape') and arr.shape != ():
            print(f"  {k}: shape={arr.shape}, dtype={arr.dtype}, min={arr.min():.4g}, max={arr.max():.4g}")
        else:
            print(f"  {k}: {arr}")
    
    print(f"\nnpz_inspect v{__version__} complete.")
