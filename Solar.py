"""
Module to implement PV functionality
Disclaimer: ChatGPT and Gemini 2.5 Pro used for assistance

Filename: Solar.py
Author: Bailey Stout
Date: 2025-05-02
"""

import math


class Solar:
    """
    Class to represent a PV Generator operating in Constant Power (PQ) mode.
    """
    def __init__(self, name: str, bus: str, real_power: float, power_factor: float, pf_mode: str):
        """
        Constructor for a Solar object
        :param name: Name of the PV generator
        :param bus: Bus connection to the generator
        :param real_power: Real power supplied by generator in [MW]
        :param power_factor: Constant power factor between 0 -> 1
        :param pf_mode: String for whether power factor is "leading" or "lagging"
        """
        self.name = name
        self.bus = bus
        self.type = 'PV_PQ'

        # Validation for power factor
        if not (0 < power_factor <= 1.0):
            raise ValueError("Target Power Factor must be between 0 (exclusive) and 1 (inclusive).")
        self.power_factor = power_factor

        if pf_mode.lower() not in ['lagging', 'leading']:
            raise ValueError("pf_mode must be either 'lagging' or 'leading'.")
        self.pf_mode = pf_mode.lower()

        # Validation for real power
        if real_power < 0:
            raise ValueError("Real power supplied by PV generator must be positive")
        self.real_power = real_power * 1e6

        # Set attributes and perform initial calculations
        self.reactive_power = self.calculate_q()
        self.apparent_power = self.real_power + 1j * self.reactive_power

    def calculate_q(self) -> float:
        """
        Calculates reactive power based on real power and power factor.
        :return: Reactive power in VARs
        """
        # Edge cases where power is ~0 or pf is 1.0
        if self.real_power < 1e-9:
            return 0.0
        if self.power_factor == 1.0:
            return 0.0

        # Calculating reactive power
        angle_rad = math.acos(self.power_factor)
        q_magnitude = abs(self.real_power) * math.tan(angle_rad)

        # Determine the sign of Q based on mode
        if self.pf_mode == 'lagging':
            q_vars = q_magnitude
        else:
            q_vars = -q_magnitude

        return q_vars

    def set_real_power(self, real_power: float):
        """
        Set function for real power injected
        :param real_power: Real power (P) injected in [MW]
        :return:
        """
        # Validation for real power
        if real_power < 0:
            raise ValueError("Real power supplied by PV generator must be positive")
        self.real_power = real_power * 1e6

        self.reactive_power = self.calculate_q()
        self.apparent_power = self.real_power + 1j * self.reactive_power

    def set_power_factor(self, power_factor: float, pf_mode: str):
        """
        Set function for power factor
        :param power_factor: Power factor between 0 -> 1
        :param pf_mode: String for whether power factor is "leading" or "lagging"
        :return:
        """
        # Validation for power factor
        if not (0 < power_factor <= 1.0):
            raise ValueError("Target Power Factor must be between 0 (exclusive) and 1 (inclusive).")
        self.power_factor = power_factor

        if pf_mode.lower() not in ['lagging', 'leading']:
            raise ValueError("pf_mode must be either 'lagging' or 'leading'.")
        self.pf_mode = pf_mode.lower()

        self.reactive_power = self.calculate_q()
        self.apparent_power = self.real_power + 1j * self.reactive_power

    def get_s(self):
        """
        Get function for apparent power
        :return: Apparent power [MVA]
        """
        return self.apparent_power / 1e6

    def __str__(self):
        """
        Display function
        :return:
        """
        return (f"PVGenerator(Name: {self.name}, Bus: {self.bus}, "
                f"P: {self.real_power / 1e6:.3f} MW, Q: {self.reactive_power / 1e6:.3f} MVAR, "
                f"Type: {self.type})")


if __name__ == '__main__':
    pv_plant_lag = Solar(name="solar1", bus="bus5", real_power=50.0, power_factor=0.95,
                         pf_mode='lagging')
    print(pv_plant_lag)
    print(f"  Complex Power S = {pv_plant_lag.get_s():.3f} MVA")

    pv_plant_lead = Solar(name="solar2", bus="bus8", real_power=20.0, power_factor=0.98,
                          pf_mode='leading')
    print(pv_plant_lead)
    print(f"  Complex Power S = {pv_plant_lead.get_s():.3f} MVA")

    pv_plant_unity = Solar(name="solar3", bus="bus10", real_power=10.0, power_factor=1.0,
                           pf_mode='lagging')
    print(pv_plant_unity)
    print(f"  Complex Power S = {pv_plant_unity.get_s():.3f} MVA")

    pv_plant_lag.set_real_power(real_power=25.0)
    print("After P change (25 MW):")
    print(pv_plant_lag)
    print(f"  Complex Power S = {pv_plant_lag.get_s():.3f} MVA")

    pv_plant_lag.set_real_power(real_power=50.0)
    print("After P change (50 MW):")
    print(pv_plant_lag)
    print(f"  Complex Power S = {pv_plant_lag.get_s():.3f} MVA")

    pv_plant_lag.set_power_factor(power_factor=0.9, pf_mode="lagging")
    print("After pf change (0.9 lagging):")
    print(pv_plant_lag)
    print(f"  Complex Power S = {pv_plant_lag.get_s():.3f} MVA")

    pv_plant_lag.set_real_power(real_power=0.0)
    print("After P change (0 MW):")
    print(pv_plant_lag)
    print(f"  Complex Power S = {pv_plant_lag.get_s():.3f} MVA")
