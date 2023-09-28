# `PyMoS2`
Python Modules for a Static Star Model
---

`PyMoS2` is a program that aims to offer a comprehensive simulation of a simple model of the interior of a static, non-magnetic star of 1 solar mass $(1$ $M_{\odot})$, from the production of energy through nuclear reactions, the radiative and convective transport of energy and the movement of a convective cell.

## Key Features

- `Stellar_core`: models energy production in the star's core. It implements the reactions involved in the PP chain and from this takes the profile of the relevant parameters depending on the radius.

- `Stellar_structure`: models the structure of the star from the core to the surface and delimits the radiative and convective transport zones in the stellar interior.

- `Stellar_convection`: models the hydrodynamic movement of a convective cell in the stellar interior and generates an animation of the movement.

---

### `Stellar_core`
The graphs show the behavior of mass, density, pressure, temperature, luminosity and energy production as the radius changes. Luminosity, mass and radius are scaled with corresponding values for the Sun. Pressure, density and energy production are shown with logarithmic axes. The dashed vertical line in the luminosity graph marks the beginning of the nucleus $(L < 0.995L_0)$.

![img1](https://github.com/GabrielWendell/PyMoS2/blob/main/img/Stellar_core-img1.png)



### `Stellar_structure`
Cross section of the best star model. The convection and radiation energy transport layers are marked on the graph. The total energy production divided by its maximum value, $\varepsilon/\varepsilon_{max}.$, plotted together with the relative energy production of each of the three branches of the PP chain.

![img2](https://github.com/GabrielWendell/PyMoS2/blob/main/img/Stellar_structure-img2.png)
![img3](https://github.com/GabrielWendell/PyMoS2/blob/main/img/Energy_prod.png)



### `Stellar_convection`
Animation snapshots of temperature (indicated by colors) and velocity (indicated by vectors) at different execution times, corresponding (from top to bottom) to the $1^{\circ}$, $45^{\circ}$, $105^{\circ}$ and $180^{\circ}$ second of the simulation.

![img4](https://github.com/GabrielWendell/PyMoS2/blob/main/img/Frames_convection.png)
