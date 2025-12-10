Welcome to Simulation of fluid flow analysis with Lattice Boltzmann method using Python

Current Version v3:
includes Zou He inlet, Bounce Back approach applied to top, bottom and corners, right periodic outlet

Why Bounce Back for corners?
Using pure bounce-back at corners is simpler and stable, and gives practically identical results for most channel flows.
Full Zou-He corner treatment requires more computational power. You will not see a significant difference in the velocity profile or pressure drop unless your Re is high.

Placement of Cylinder:
1. Our domain is 250 x 50. At the moment, it seems optimal for the Pouseille flow. A cylinder of diameter D = 20 seems like a good option to start with, considering the domain size.
Place the cylinder exactly in the center to allow for flow velocity upstream and downstream flow to develop.
2. You can increase it if you would like. However, keep the ratio of nx/ny >4 for optimal condition
3. Keep in mind that Reynolds number is affected if ny is changed.
4. More computational power is required for huge domains.


Fix for improper flow connectivity at outlet and inlet:
Flow is diverging instead of smooth flow from periodic outlet to inlet. It could be because of using Zou He inlet and periodic outlet, which is causing imbalance in mass
and likely diverging flow

Fix : Using Periodic inlet and outlet, driving force using body force gives optimal flow in the channel


The main reason is the low Mach number requirement for incompressible flow simulation:
Mach number = u / c_s
where c_s = 1/√3 ≈ 0.577 is the speed of sound in lattice units
For incompressible flow, we need: Ma < 0.3
So: u_max ≈ 0.3 × 0.577 ≈ 0.173 in lattice units
That's why we typically keep velocities around 0.1-0.15 in lattice units.
