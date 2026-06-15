# ciliate-mouth-design

[![DOI](https://img.shields.io/badge/DOI-10.1038%2Fs41467--025--59413--x-blue)](https://doi.org/10.1038/s41467-025-59413-x)
[![Journal](https://img.shields.io/badge/Nature%20Communications-2025-green)](https://www.nature.com/articles/s41467-025-59413-x)

**Flow physics of nutrient transport drives functional design of ciliates**
Jingyi Liu, John H. Costello, Eva Kanso — *Nature Communications* **16**, 4154 (2025)

> Code and models for our study of *why ciliates put their cilia where they do* — how the physics of nutrient transport selects the placement of a feeding "mouth" on a single-celled organism.

📄 [Read the paper](https://www.nature.com/articles/s41467-025-59413-x)

---

## Summary

Phagotrophy — a cell's ability to ingest particles — was a pivotal evolutionary step, but it created a mechanical problem: how does a cell move food particles from the surrounding fluid to its interior? We study this in **ciliates**, single-celled eukaryotes that either swim freely or attach to a surface and beat their cilia to drive feeding currents.

Using mechanistic hydrodynamic models combined with a survey of real ciliate morphologies, we show that cells maximize feeding efficiency by devoting only a **specific region of their surface to a "mouth,"** and that the optimal cilia arrangement depends on life strategy.

## Key findings

- Ciliates optimize feeding efficiency by designating a localized portion of the cell surface as a "mouth" rather than feeding over the whole surface.
- **Optimal cilia coverage depends on lifestyle:** sessile (attached) ciliates feed most efficiently with cilia arranged in *bands around oral structures*, while swimming ciliates tolerate more *diverse* ciliary arrangements that still meet their nutritional needs.
- Past a threshold — roughly a doubling of nutrient uptake — further gains in feeding flux are **not a dominant selective force** in cell design, suggesting other pressures take over.
- Model-predicted optima (including the Pareto fronts trading off uptake against cost) line up closely with the morphologies of surveyed real ciliates.

## What's in this repo

<!-- TODO: adjust to your actual layout -->
```
ciliate-mouth-design/
├── models/          # hydrodynamic feeding-current models (envelope / Stokeslet)
├── optimization/    # feeding-efficiency optimization & Pareto-front computation
├── survey/          # digitized morphology/flow data for surveyed ciliates
├── figures/         # scripts to regenerate paper figures
└── notebooks/       # walkthroughs reproducing the main results
```

## Reproducing the results

<!-- TODO: fill in real commands -->
```bash
git clone https://github.com/jingyiliu1900/ciliate-mouth-design.git
cd ciliate-mouth-design
pip install -r requirements.txt   # or: conda env create -f environment.yml
# e.g. reproduce Figure 4 (optimal mouth coverage vs. lifestyle):
python figures/fig4_mouth_coverage.py
```

## Citation

```bibtex
@article{liu2025flow,
  title   = {Flow physics of nutrient transport drives functional design of ciliates},
  author  = {Liu, Jingyi and Costello, John H. and Kanso, Eva},
  journal = {Nature Communications},
  volume  = {16},
  number  = {1},
  pages   = {4154},
  year    = {2025},
  doi     = {10.1038/s41467-025-59413-x}
}
```

## Part of the ciliate feeding series

[`swim-or-stay`](https://github.com/jingyiliu1900/swim-or-stay) · [`optimal-ciliary-feeding`](https://github.com/jingyiliu1900/optimal-ciliary-feeding) · [`feeding-in-gradients`](https://github.com/jingyiliu1900/feeding-in-gradients)

## Authors

Jingyi Liu · John H. Costello (Providence College / MBL) · Eva Kanso (USC)
*Kanso Lab, University of Southern California.*

## License

<!-- TODO: choose one. MIT is a common, permissive default for research code. -->
Released under the MIT License (see `LICENSE`).
