Welcome to Simulation of fluid flow analysis with Lattice Boltzmann method using Python

v3 includes Zou He inlet, Bounce Back approach applied to top, bottom and corners, right periodic outlet

Why Bounce Back for corners?
Using pure bounce-back at corners is simpler and stable, and gives practically identical results for most channel flows.
Full Zou-He corner treatment requires more computational power. You will not see a significant difference in the velocity profile or pressure drop unless your Re is high.