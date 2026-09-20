# 2-D Modal Active Noise Control (ANC) with k-Wave

A MATLAB reference implementation of a **2-D circular-harmonic (modal) multichannel active noise control**
system, simulated on a rigid barrier-with-slit acoustic domain using the
[k-Wave](http://www.k-wave.org/) acoustics toolbox.

The controller is a three-stage virtual-microphone scheme:

1. **Stage 1 – Modal tuning (`ModuleANC_MV_2D`)** tunes the control filters on a circle of
   *virtual* error microphones by minimising a circular-harmonic regional cost.
2. **Stage 2 – Auxiliary fitting (`AuxiliaryLMS`)** fits auxiliary filters that map the virtual
   residuals back to the *physical* error microphones.
3. **Stage 3 – Online control (`ContrFxLMS`)** runs the converged controller on the real primary
   noise and updates the filters on the physical residuals.

A final k-Wave render shows the resulting sound-pressure-level (SPL) map.

---

## Contents

| File | Role |
|------|------|
| [`kwave_Module_test.m`](#kwave_module_testm) | **Main script.** Builds the scenario, identifies secondary paths, runs the three ANC stages, and renders the final field. |
| [`ModuleANC_MV_2D.m`](#moduleanc_mv_2dm) | Stage 1: 2-D circular-harmonic-domain multichannel FxLMS. |
| [`AuxiliaryLMS.m`](#auxiliarylmsm) | Stage 2: auxiliary (virtual-sensing) LMS filter bank. |
| [`ContrFxLMS.m`](#contrfxlmsm) | Stage 3: online multichannel filtered-x LMS. |
| [`Ref_Sig_Gene.m`](#ref_sig_genem) | Generates filtered-reference signals from the identified secondary paths. |
| [`CreatReferenceSignal.m`](#creatreferencesignalm) | Convenience wrapper that also produces the primary disturbances. |
| [`Noise_canbine_field.m`](#noise_canbine_fieldm) | Runs one k-Wave 2-D simulation and returns the source / sensor signals. |
| [`Noise_canbine_field3_Map.m`](#noise_canbine_field3_mapm) | Runs a k-Wave simulation and renders / saves the SPL map. |

---

## Requirements

- MATLAB (R2018b or newer recommended).
- [k-Wave](http://www.k-wave.org/) added to the MATLAB path (provides `kWaveGrid`,
  `kspaceFirstOrder2D`, `getColorMap`, `scaleSI`).
- DSP System Toolbox (provides `dsp.LMSFilter` for secondary-path identification).

## Usage

Open `kwave_Module_test.m` and run it. The script is self-contained: it sets up the grid,
medium, source and sensor positions, identifies the secondary paths, runs the three adaptive
stages, and produces diagnostic plots plus a final SPL map.

```matlab
kwave_Module_test
```

---

## Function reference

### `kwave_Module_test.m`

**Main script (no function interface).**

Sets up a 256 × 128 grid with air properties (`c0 = 343 m/s`, `rho0 = 1.29 kg/m³`) and a
rigid barrier (scaled impedance ×20) containing a slit. It places:

- a **primary source** at the bottom edge,
- two **secondary loudspeakers** beside the slit,
- two **physical error microphones** on the receiver side,
- ten **virtual microphones** uniformly distributed on a circle of radius `R = 0.9 m`.

It then:

1. Records the primary disturbances at the physical and virtual microphones
   (`Noise_canbine_field`).
2. Identifies the secondary paths from each secondary source to every physical and virtual
   microphone using normalised LMS (`dsp.LMSFilter`, length `M = 4096`).
3. Builds the filtered references (`Ref_Sig_Gene`).
4. Runs Stage 1 (`ModuleANC_MV_2D`), Stage 2 (`AuxiliaryLMS`), and Stage 3
   (`ContrFxLMS`).
5. Drives the two secondary speakers with the converged filters and renders the final
   anti-noise field with `kspaceFirstOrder2D`, reporting the average and error-microphone
   SPL.

Key tunables: `L_fil` (filter length), `u1/u2/u3` (adaptation step sizes), `M_vm`, `R`.

---

### `ModuleANC_MV_2D.m`

```matlab
[W, Er, Info] = ModuleANC_MV_2D(LenFilter, NumSource, NumMic, N, ...
                                 Fx, Distur, StepSize, CHOpts)
```

Stage 1 – **2-D circular-harmonic-domain multichannel FxLMS**.

Transforms the boundary residuals into short-time spectra, projects them onto a real
orthonormal circular-harmonic basis on the microphone circle, forms a regional modal cost
(optionally weighted by the area-integrated acoustic energy inside the target disk), and
updates the FIR control filters with an explicit tap-by-tap spectral gradient.

**Inputs**

| Argument | Meaning |
|----------|---------|
| `LenFilter` | `L`, length of each control filter. |
| `NumSource` | `K`, number of secondary sources. |
| `NumMic` | `M`, number of microphones uniformly distributed on the boundary circle. |
| `N` | Number of simulation samples. |
| `Fx` | Filtered references, size `(N+L-1, M, K[, J])`. |
| `Distur` | Primary disturbances at the `M` boundary mics, size `(N, M)`. |
| `StepSize` | Adaptation step size. |
| `CHOpts` | Optional struct: `Fs`, `c`, `R`, `Angle`, `Weights`, `NDFT`, `Hop`, `Band`, `MaxOrder`, `Window`, `Weighting`, `Epsilon`, `GradScale`, `Norm`. |

**Outputs**

| Argument | Meaning |
|----------|---------|
| `W` | Final control filters, size `(K*L*J, 1)`. |
| `Er` | Online residual pressures at the boundary mics, size `(M, N)`. |
| `Info` | Circular-harmonic geometry, retained orders, and convergence diagnostics. |

The file also contains local helpers: `SetDefault`, `BuildFXMatrix`, `RealCHMatrix`
(orthonormal real circular-harmonic basis), and `RadialOrderWeight2D` (area-energy weights
from Bessel-function radial integrals).

---

### `AuxiliaryLMS.m`

```matlab
[H, Er] = AuxiliaryLMS(LenFilter, NumSource, NumPM, N, Fx_p, ...
                       DisturPhysic, PriNoise, W, StepSize)
```

Stage 2 – **auxiliary (virtual-sensing) LMS**.

Fits an auxiliary filter bank `H` that predicts the residual at each physical error
microphone from the primary noise, so that the virtual-microphone controller `W` can be
monitored on the physical microphones during Stage 3.

**Inputs**

| Argument | Meaning |
|----------|---------|
| `LenFilter` | `L`, filter length. |
| `NumSource` | `K`, number of secondary sources. |
| `NumPM` | `J`, number of physical microphones. |
| `N` | Number of samples. |
| `Fx_p` | Filtered references at the physical mics, `(N+L-1, J, K)`. |
| `DisturPhysic` | Primary disturbances at the physical mics, `(N, J)`. |
| `PriNoise` | Primary reference noise, length `N`. |
| `W` | Modal control filters from Stage 1, `(K*L, 1)`. |
| `StepSize` | Adaptation step size. |

**Outputs**

| Argument | Meaning |
|----------|---------|
| `H` | Auxiliary filters, `(J*L, 1)`. |
| `Er` | Auxiliary error signals, `(J, N)`. |

---

### `ContrFxLMS.m`

```matlab
[WC, ErPhysic, ErVirt] = ContrFxLMS(LenFilter, NumSource, NumVM, NumPM, N, ...
                                    Fx_p, Fx_v, DisturPhysic, DisturVirt, ...
                                    PriNoise, H, StepSize)
```

Stage 3 – **online multichannel filtered-x LMS**.

Runs the converged controller on the actual primary noise. The auxiliary filters `H` are
used to monitor the physical residuals, while the control filters `WC` are adapted on those
physical residuals; the virtual residuals are also reported for diagnosis.

**Inputs**

| Argument | Meaning |
|----------|---------|
| `LenFilter` | `L`, filter length. |
| `NumSource` | `K`, number of secondary sources. |
| `NumVM` | `M`, number of virtual microphones. |
| `NumPM` | `J`, number of physical microphones. |
| `N` | Number of samples. |
| `Fx_p` / `Fx_v` | Filtered references at the physical / virtual mics. |
| `DisturPhysic` / `DisturVirt` | Primary disturbances at the physical / virtual mics. |
| `PriNoise` | Primary reference noise, length `N`. |
| `H` | Auxiliary filters from Stage 2, `(J*L, 1)`. |
| `StepSize` | Adaptation step size. |

**Outputs**

| Argument | Meaning |
|----------|---------|
| `WC` | Online control filters, `(K*L, 1)`. |
| `ErPhysic` | Residuals at the physical microphones, `(J, N)`. |
| `ErVirt` | Residuals at the virtual microphones, `(M, N)`. |

---

### `Ref_Sig_Gene.m`

```matlab
[Fx_v, Fx_p] = Ref_Sig_Gene(Sv, Sp, PriNoise, N, L, K, M, J)
```

Convolves the primary noise with the estimated secondary paths (`Sv` to the virtual
mics, `Sp` to the physical mics) and pads the result with `L-1` leading zeros, producing
the filtered-reference arrays consumed by the adaptive stages.

**Outputs:** `Fx_v` `(N+L-1, M, K)` and `Fx_p` `(N+L-1, J, K)`.

---

### `CreatReferenceSignal.m`

```matlab
[Dv, Dp, Fx_v, Fx_p] = CreatReferenceSignal(Pv, Pp, Sv, Sp, PriNoise, N, L, K, M, J)
```

Extended version of `Ref_Sig_Gene` that also convolves the primary noise with the estimated
**primary paths** `Pv` / `Pp` to directly produce the primary disturbances at the virtual
(`Dv`) and physical (`Dp`) microphones. Useful when primary paths are available from
measurement; the main script instead obtains the disturbances directly from k-Wave.

---

### `Noise_canbine_field.m`

```matlab
[Xin, d] = Noise_canbine_field(kgrid, medium, PML_size, source, sensor, slit_mask)
```

Thin wrapper around `kspaceFirstOrder2D`. Runs one forward simulation with the given
source / sensor, plots the final pressure field and the input / output time waveforms, and
returns the source signal `Xin` and the sensor recordings `d` (both transposed to row
vectors, cast to `double`).

---

### `Noise_canbine_field3_Map.m`

```matlab
Noise_canbine_field3_Map(kgrid, medium, PML_size, source, sensor, slit_mask, ...
                         Nx, Ny, Index, dx, dy, slit_x_pos, FileName)
```

Variant that records the **full-field RMS pressure**, converts it to dB (referenced to
20 µPa), prints the average SPL and the SPL at the error microphones, plots the SPL map
and a contour plot over the receiver region, and saves the SPL map (`map`) and the raw
sensor data (`RMS`) to `FileName.mat`.

---

## Notes

- All adaptive stages are sample-by-sample MATLAB loops prioritising readability and
  algorithmic clarity; they are not optimised for speed.
- Secondary paths are identified online by normalised LMS in the main script before the
  adaptive stages run.
- The virtual-microphone circle and the barrier slit are the key geometric knobs to
  experiment with (`M_vm`, `R`, `slit_*`).

## License

See repository for license details.
