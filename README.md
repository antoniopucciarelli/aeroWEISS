# Weissinger method

This program studies aerodynamics properties of a 3D wing. The program is based on **Weissinger** method for an incompressible irrotational steady airflow.

## Physics

Given incompressibility, irrotationality and steadiness of the flow, the problem can be described with as **potential** flow. This flow can be discritized with horseshoe vortices, that are essentially lines with a constant vortex distribution on them. These lines make the external flow irrotational but also allow us the study of **lift** and **drag** (*induced* drag, so nothing related to shear stress or flow separation) for the 3D wing. These lines are rotational in the **core**; the external flow will be irrotational. These vortex lines, that make the horseshoe vortices, will be distributed on the wing **surface** and on the **wake**.

The **total** lift will be computed as the *sum* of all the vortex lines over the wing; instead the **drag** will be computed as a *flow deviation* due the presence of the **wake** and it will be related to angle of deflection of the flow and to the lift.

## Program steps

* Geometry generation
  * 3D wing paneling
  * Computing panel normal
  * Computing horseshoe vortices geometry
* Flow computation with respect to horseshoe vortices distribution
  * Sum of all the vortices distribution over each panel
  * Imposing non penetration condition over each panel
* Solving linear system
* Computing lift as sum of all circulation from the horseshoe vortices
* Computing drag as sum of all panels lift force corrected by the deviation angle due to wake presence (induced drag)
* Results plotting

## Running different simulation

All the inputs must be written in the ```WEISSINGER.m``` file, then launching the program. The program is divided in different subsections due to the presence of different simulation such as **single** 3D wing, **double** 3D wings and also **ground effect** with **mirroring** method.

## Validation

Results are validated through **XFLR5** and **AVL**.
