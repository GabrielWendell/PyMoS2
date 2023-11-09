![image](https://github.com/GabrielWendell/PyMoS2/blob/main/Logo/PyMoS2_Color_OldVersion.jpg)

- Art maded by: **Vitor Santos** ([Vitorsodb@gmail.com](mailto:Vitorsodb@gmail.com)).
- You can check his portfolio here: [Viktor Freakstein](https://www.behance.net/gentlesuspenders_).

---

`PyMoS2` is a program that aims to offer a comprehensive simulation of a simple model of the interior of a static, non-magnetic star of $1$ solar mass $(1$ $M_{\odot})$, from the production of energy through nuclear reactions, the radiative and convective transport of energy and the movement of a convective cell.

## Key Features

- `Stellar_core`: models energy production in the star's core. It implements the reactions involved in the PP chain and from this takes the profile of the relevant parameters depending on the radius.

- `Stellar_structure`: models the structure of the star from the core to the surface and delimits the radiative and convective transport zones in the stellar interior.

- `Stellar_convection`: models the hydrodynamic movement of a convective cell in the stellar interior and generates an animation of the movement.

---

## Installation
First you need to check if you have git software installed on your machine. You can check by typing the following in your terminal:
```
$ git --version
```

If no version appears, then you will need to download it. Depending on your operating system, you can install it as follows:
> I. **Mac OSX**
> - Install the package from [https://git-scm.com/download/mac](https://git-scm.com/download/mac).

> II. **Linux**
> - Depending on your distribution:
>   
> **A**. Debian and Ubuntu-based distros:
> ```
> $ sudo apt-get install git
> ```
> **B**. Red-Hat-based distros:
> ```
> $ sudo yum install git
> ```

Once done, from a new terminal create a directory dedicated for this tool with the name PyMoS2, and from it, clone this github files:
```
git clone https://github.com/GabrielWendell/PyMoS2.git
```
This will download the code and instructions =D

---

> Author: Gabriel Wendell Celestino Rocha ([gabrielwendell@fisica.ufrn.br](mailto:gabrielwendell@fisica.ufrn.br)).
