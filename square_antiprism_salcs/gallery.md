# Illustrated D₄d square-antiprism SALCs

These are the **32 ligand s/p SALCs** derived in [the mathematical analysis](derivation.md). The companion [bonding diagram](bonding_diagram.md) shows how they can mix with a central atom’s s/p/d orbitals. From the project root, install the figure-generation dependencies with `python -m pip install -r square_antiprism_salcs/requirements.txt`, then run `python square_antiprism_salcs/generate_figures.py` to regenerate and validate every image. The drawings use ρ = 1 and h/ρ = 2<sup>−3/4</sup>; changing that ratio leaves the symmetry coefficients unchanged.

Blue and orange mean **positive and negative wavefunction phase**, not positive and negative charge. An s AO is drawn as one sphere; a p AO has two opposite-phase lobes. The lower bar plot shows the **signed coefficient of the local AO** at each site. A missing lobe and zero-height bar mean zero coefficient. Lobe sizes are schematic, not electron-density isosurfaces. Site order is T0, B0, T1, B1, T2, B2, T3, B3, alternating top and bottom around the S₈ orbit. The subscript in each gallery heading is a **Fourier mode**, not a site number. The exact coefficients, including the 1/2 and 1/√8 factors, are in the table at the end of the [analysis](derivation.md#complete-coefficient-table).

For each E pair, the cosine and sine drawings are **two partners of one irrep**. The pair can be orthogonally recombined without changing the represented symmetry or the span of the pair.

## Ligand s orbitals

### 01 · A₁ · s₀

All eight s orbitals have equal positive coefficients. This is the totally symmetric ligand-s combination.

![A1 s SALC, uniform positive coefficients on all eight sites](figures/01_s_A1_m0.png)

### 02 · B₂ · s₄

Top-site coefficients are positive and bottom-site coefficients negative, all with equal magnitude.

![B2 s SALC, opposite phase on top and bottom squares](figures/02_s_B2_m4.png)

### 03 · E₁ · s₁,c

The first cosine partner has its strongest positive coefficient at T0 and strongest negative one at T2; T1 and T3 are nodes.

![E1 s cosine partner with maxima at T0 and T2](figures/03_s_E1_m1c.png)

### 04 · E₁ · s₁,s

The sine partner shifts that pattern around the orbit: its extrema are at T1 and T3, with nodes at T0 and T2.

![E1 s sine partner with maxima at T1 and T3](figures/04_s_E1_m1s.png)

### 05 · E₂ · s₂,c

Only the top four sites contribute, alternating positive and negative around the top square.

![E2 s cosine partner on the top square only](figures/05_s_E2_m2c.png)

### 06 · E₂ · s₂,s

Only the bottom four sites contribute, likewise alternating phase around their square.

![E2 s sine partner on the bottom square only](figures/06_s_E2_m2s.png)

### 07 · E₃ · s₃,c

This higher-winding cosine pattern again has extrema at T0 and T2, but its bottom-site signs differ from those of the E₁ cosine partner.

![E3 s cosine partner with threefold Fourier winding](figures/07_s_E3_m3c.png)

### 08 · E₃ · s₃,s

The companion higher-winding sine pattern has opposite extrema at T1 and T3 and a different sequence of bottom-site phases.

![E3 s sine partner with threefold Fourier winding](figures/08_s_E3_m3s.png)

## Bond-directed pσ orbitals

Here each local +p axis points **outward along the line from the central atom to its ligand**. Opposite colors on the two lobes of one site are the normal sign change of a p orbital; the lower bar gives the overall SALC coefficient for that site.

### 09 · A₁ · σ₀

All eight outward-pointing local pσ functions enter with the same positive coefficient. Their differently oriented lobes together make an A₁ combination.

![A1 bond-directed p SALC with all local axes in phase](figures/09_sigma_A1_m0.png)

### 10 · B₂ · σ₄

The four top bond-directed p functions have positive coefficients; the four bottom ones have negative coefficients.

![B2 bond-directed p SALC with top and bottom coefficients opposed](figures/10_sigma_B2_m4.png)

### 11 · E₁ · σ₁,c

The bond-directed p functions follow the m = 1 cosine weights: largest at T0/T2, absent at T1/T3, and smaller on the bottom sites.

![E1 bond-directed p cosine partner](figures/11_sigma_E1_m1c.png)

### 12 · E₁ · σ₁,s

The m = 1 sine partner shifts the extrema to T1/T3 and the nodes to T0/T2.

![E1 bond-directed p sine partner](figures/12_sigma_E1_m1s.png)

### 13 · E₂ · σ₂,c

Only top-site bond-directed p functions remain, with alternating SALC signs.

![E2 bond-directed p cosine partner on the top square](figures/13_sigma_E2_m2c.png)

### 14 · E₂ · σ₂,s

Only bottom-site bond-directed p functions remain, again alternating in sign.

![E2 bond-directed p sine partner on the bottom square](figures/14_sigma_E2_m2s.png)

### 15 · E₃ · σ₃,c

The m = 3 cosine weights give a more rapidly changing phase sequence than E₁ while retaining T0/T2 extrema.

![E3 bond-directed p cosine partner](figures/15_sigma_E3_m3c.png)

### 16 · E₃ · σ₃,s

The m = 3 sine partner has T1/T3 extrema and the complementary higher-winding phase sequence.

![E3 bond-directed p sine partner](figures/16_sigma_E3_m3s.png)

## Meridional pπ orbitals

Each local +p axis lies in the plane formed by the site and the z-axis, **perpendicular to the metal–ligand bond**. Under S₈ that local axis reverses, so identical-looking coefficient sequences can carry different irrep labels from the s or pσ sequences.

### 17 · B₂ · μ₀

Every local meridional p function has a positive coefficient, yet the whole SALC is B₂ because S₈ reverses these local axes.

![B2 meridional p SALC with uniform local coefficients](figures/17_meridional_B2_m0.png)

### 18 · A₁ · μ₄

Alternating top/bottom coefficients compensate the meridional p-axis reversal, producing the totally symmetric A₁ combination.

![A1 meridional p SALC with alternating local coefficients](figures/18_meridional_A1_m4.png)

### 19 · E₃ · μ₁,c

This m = 1 cosine coefficient pattern is E₃, **not** E₁, because the p-axis sign adds a minus under S₈.

![E3 meridional p cosine partner for m equals one](figures/19_meridional_E3_m1c.png)

### 20 · E₃ · μ₁,s

The corresponding m = 1 sine partner has the same E₃ character and complementary extrema and nodes.

![E3 meridional p sine partner for m equals one](figures/20_meridional_E3_m1s.png)

### 21 · E₂ · μ₂,c

Top meridional p functions alone contribute, alternating signs; the pair's S₈ character is zero.

![E2 meridional p cosine partner on the top square](figures/21_meridional_E2_m2c.png)

### 22 · E₂ · μ₂,s

Bottom meridional p functions alone contribute, alternating signs; this is the second E₂ partner.

![E2 meridional p sine partner on the bottom square](figures/22_meridional_E2_m2s.png)

### 23 · E₁ · μ₃,c

The m = 3 cosine pattern becomes E₁ because the meridional axis reversal changes the sign of the S₈ character.

![E1 meridional p cosine partner for m equals three](figures/23_meridional_E1_m3c.png)

### 24 · E₁ · μ₃,s

The m = 3 sine partner completes that E₁ set, with T1/T3 extrema and T0/T2 nodes.

![E1 meridional p sine partner for m equals three](figures/24_meridional_E1_m3s.png)

## Tangential pπ orbitals

Each local +p axis points around its square in the direction of increasing polar angle. Reflection in the xz plane reverses this local axis, so the uniform tangential combination is A₂ rather than A₁.

### 25 · A₂ · τ₀

All tangential p functions have equal positive local coefficients; reflection reverses the resulting SALC's phase.

![A2 tangential p SALC with uniform local coefficients](figures/25_tangential_A2_m0.png)

### 26 · B₁ · τ₄

The tangential p coefficients alternate between top and bottom. The extra mirror sign makes the irrep B₁.

![B1 tangential p SALC with alternating local coefficients](figures/26_tangential_B1_m4.png)

### 27 · E₁ · τ₁,c

The tangential m = 1 cosine partner has T0/T2 extrema and zero coefficients at T1/T3.

![E1 tangential p cosine partner for m equals one](figures/27_tangential_E1_m1c.png)

### 28 · E₁ · τ₁,s

Its sine partner shifts the extrema to T1/T3 and the top-site nodes to T0/T2.

![E1 tangential p sine partner for m equals one](figures/28_tangential_E1_m1s.png)

### 29 · E₂ · τ₂,c

Only the four top tangential p functions contribute, with alternating signs.

![E2 tangential p cosine partner on the top square](figures/29_tangential_E2_m2c.png)

### 30 · E₂ · τ₂,s

Only the four bottom tangential p functions contribute, with alternating signs.

![E2 tangential p sine partner on the bottom square](figures/30_tangential_E2_m2s.png)

### 31 · E₃ · τ₃,c

The m = 3 cosine partner uses the faster phase-winding pattern of E₃, with extrema at T0/T2.

![E3 tangential p cosine partner for m equals three](figures/31_tangential_E3_m3c.png)

### 32 · E₃ · τ₃,s

The final sine partner has extrema at T1/T3 and completes the 32-function SALC basis.

![E3 tangential p sine partner for m equals three](figures/32_tangential_E3_m3s.png)
