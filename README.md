# Advances in Exoplanet Research — Practical Materials

This repository contains practical materials for the **Advances in Exoplanet Research** course. The notebooks accompany hands-on classes about protoplanetary disk evolution and the behavior of dust in disks.

## Materials

- `LBP_full.ipynb` and `LBP_student.ipynb`: viscous evolution of a one-dimensional gas disk, including a dead zone.
- `Dust_drift.ipynb` and `Dust_drift_student.ipynb`: radial dust drift in an evolving gas disk. `Dust_drift.c` contains the model, and `Dust_drift.yaml` contains its configuration.
- `dust_settling.ipynb` and `dust_settling_student.ipynb`: vertical settling of dust grains in a gas disk.

The full notebooks show complete examples; the student notebooks leave selected code lines for students to complete.

## Running the notebooks

Use Jupyter with Python, NumPy, and Matplotlib. The dust-drift notebook also needs a C compiler and libyaml; run it from the repository root so it can find `Dust_drift.c` and `Dust_drift.yaml`.

## License

These materials are provided under the MIT License. See [LICENSE](./LICENSE).
