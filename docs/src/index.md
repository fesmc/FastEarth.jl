# FastEarth.jl

Global 2D Glacial Isostatic Adjustment (GIA) model, generalising
[FastIsostasy](https://gmd.copernicus.org/articles/17/5263/2024/) from a
regional Cartesian domain to the whole sphere.

FastEarth represents the solid Earth as a thin elastic lithosphere over a
laterally variable viscous mantle, with the 3D radial viscosity profile
collapsed into a 2D wavelength-dependent effective field via layer-stacking.
Self-gravitation is handled kinematically through a gravitationally
self-consistent sea-level equation. The model is designed for interactive
coupling to ice-sheet models at 10–50 km resolution over paleo and
future-projection timescales.

## Documentation

- **[Design](design.md)** — governing equations, numerical strategy,
  architecture, validation hierarchy.

## Status

Pre-implementation. The [design document](design.md) is the canonical
specification.

## Repository

[github.com/fesmc/FastEarth.jl](https://github.com/fesmc/FastEarth.jl)
