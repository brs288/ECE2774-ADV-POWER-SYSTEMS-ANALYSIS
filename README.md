# ECE2774 Power Simulation Software - Project 3 PV Enhancement

This software was created in our ECE2774 Advanced Power Systems Analysis class. It can be used to simulate power systems and conduct fault analysis at a specified bus. This particular enhancement introduces photovoltaic generation capability in the simulation.

## Table of Contents

-   [Purpose](#purpose-and-theoretical-background)
-   [Inputs and Outputs](#inputs-and-outputs)
-   [Reference List](#reference-list)
-   [Testing and Validation](#testing-and-validation)

## Purpose and Theoretical Background

This PV simulation is integrated into the simulation in a modular manner, relying on previously built framework. The focus of this addition was for power flow analysis, but faults can still be studied just as they were in the previous simulation. Furthermore, fundamental architecture such as buses, generators, loads, transmission lines, and transformers all cooperate cleanly with the newly added solar class. The primary goal of this addition is to study how PV generation and load fluctuations impact power systems over time and how this leads to the Duck Curve. For more detail, see the overview in the Project 3 documentation pdf.

## Inputs and Outputs

The primary additional inputs for this enhanced simulator included solar generator name, bus connection, maximum real power, and target power factor. These inputs define the characteristics of a PV generator in the Solar class. To perform a sweep over time, optional inputs to the `solar_sweep` method include a boolean to determine whether the load varies over time as well and a boolean dictating whether the resulting data is written to a csv. For more information on specific additions and usage see the Project 3 documentation pdf.

## Reference List

This fundamental software was created by Justin Lipner and Bailey Stout. The PV enhancement was developed by Bailey Stout.

-   ChatGPT 4o was used to develop conceptual understanding and build framework for code
-   Gemini 2.5 Pro was used similarly
-   NREL was used to build conceptual understanding - <https://www.nrel.gov/pv/>
-   Textbook: J. Glover, M. Sarma, T. Overbye, A. Birchfield, Power Systems Analysis and Design, 7th Edition. Cengage, 2023

## Testing and Validation

To perform a surface level test with this simulator, `main.py` can be run. This will run both test scenarios (these functions reside in `Validations.py`) as well as some edge cases. For further testing, some pre-defined booleans are listed and can be set to `True` to examine failure cases. For users to develop their own simulations, first build a desired circuit object. Refer to documentation on Project 2 if need be. To append a solar generator, simply call the circuit method `add_solar()` and include desired parameters. For example, as is done in the second scenario, calling `circ.add_solar("Solar_Scenario2", "bus7", 20, 1, "unity")` will create a solar generator at `"bus7"` with a maximum of 20 MW real power output at unity power factor. To perform a sweep over time, call the `solar_sweep()` method from the circuit class. For example, the second scenario accomplishes this with `circ.sweep_solar(True, True)`, where the first boolean is for time-varying loads and the second boolean is for csv output. To see further details on usage and scenario results/implementation, see the Project 3 documentation pdf.