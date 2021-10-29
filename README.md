# Generators

Experiments around the EAMC and space-time committors / coherence.
The state of this repository is pretty messy and is left here for archival use.
The Augmented jump chain code has been rewritten for python and julia in `cmdtools` and `birthdeath`

Notable features: Coherence of meteorology data, PCCA+, Bickley Jet

## Contents

### Source
- Galerkin discretization of the Augmented Jump Chain (called Embedded augmented markov chain (EAMC) during that time): `src/galerkin.jl`
- EAMC reconstruction of the Koopman operator by space-time committor: `src/committor.jl`
- EAMC reconstruction of the Perron-Frobenius operator by jump activity: `src/perronfrobenius.jl`
- Gillespie like simulation of trajectories of the EAMC: `src/gillespie.jl`
- Opimization of space-time committors: `src/main.jl`
- PCCA+ implementation: `src/pccap.jl`
- Some interactive (?) Makie visualization: `src/plot.jl`

### (Nonautonomous) Models
- 4 Cell barrier on-off switching: `barrierswitch.jl`
- The bickley jet: `bickleyjet.jl`
- Fluid flow governed by vectorfield: `continuity.jl`
- Overdamped Langevin via SQRA, also shifting and rising doublewells: `langevin.jl`
- 1d sqra: `mysqra.jl`

### Notebooks
- `example*`: EAMC plots of different systems
- `bickley*`: bickley coherence, plots not saved
- `weather-inspect`: look at the weather data
- `weather-kmeans/pcca`: clustering of the eamc for coherent sets
- `weather-optimize`: clustering via optimization criterion
