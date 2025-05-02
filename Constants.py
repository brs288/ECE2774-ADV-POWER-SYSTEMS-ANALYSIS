"""
Module for constants throughout power system calculations

Filename: Constants.py
Author: Justin Lipner
Date: 2025-02-04
"""

j = 1j
epsilon = 8.854*10**-12
inch2m = 0.0254
feet2m = 0.3048
mi2m = 1609.34
solar_profile = [0, 0, 0, 0, 0, 0.1, 0.3, 0.6, 0.8, 0.95, 1.0, 1.0,
                 0.95, 0.8, 0.6, 0.3, 0.1, 0, 0, 0, 0, 0, 0, 0]
load_profile_factor = [0.45, 0.42, 0.40, 0.40, 0.42, 0.48, 0.60, 0.75, 0.85, 0.88, 0.86, 0.85,
                       0.84, 0.83, 0.84, 0.85, 0.88, 0.95, 1.00, 0.98, 0.92, 0.80, 0.65, 0.55]