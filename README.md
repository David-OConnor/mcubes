# Mcubes: Marching Cubes library

[![Crate](https://img.shields.io/crates/v/mcubes.svg)](https://crates.io/crates/mcubes)
[![Docs](https://docs.rs/mcubes/badge.svg)](https://docs.rs/mcubes)


Uses the Marching Cubes algorithm to create isosurfaces in Rust, from volume data. Designed to be easy to integrate.

Based loosely on [PyMarchingCubes](https://github.com/JustusThies/PyMarchingCubes).

For now, depends on the `graphics` library for `Mesh` and `Vertex`. We will remove this requirement later, for a
lighter dependency tree.

Uses [lin-alg](https://github.com/david-oconnor/lin-alg) for its `Vec3` type.

Used by the [Daedalus molecule viewer](https://github.com/David-OConnor/daedalus) to view experimentally-derived
electron density from protein crystals.


Example use:
```rust
use mcubes::{GridPoint, MarchingCubes};

pub struct ElectronDensity {
    pub coords: Vec3,
    pub density: f64,
}

impl GridPoint for ElectronDensity {
    fn value(&self) -> f64 { self.density }
}

fn create_mesh(hdr: &MapHeader, mol: &Molecule, iso_level: f32) {
    let mc = MarchingCubes::new(
        // Number of grid points along each axis
        (hdr.nx as usize, hdr.ny as usize, hdr.nz as usize),
        // Grid dimensions per unit cell along each axis
        (hdr.cell[0], hdr.cell[1], hdr.cell[2]),
        // Sampling interval along each axis. (Usually the same as grid point number of grid points.)
        (hdr.mx as f32, hdr.my as f32, hdr.mz as f32),
        mol.elec_density,
        // The value to draw the isosurface at.
        iso_level,
    );

    let mesh = mc.generate();
}
```

Why another Marching Cubes library? I couldn't figure out how to use the other ones.


[//]: # (todo: Screenshot[s])
