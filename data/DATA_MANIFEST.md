# Data manifest — provenance of every shipped model

Generated 2026-09-14 by `python analysis/make_manifest.py` (hashes computed, never typed); verify with `python analysis/make_manifest.py --check`. The machine-readable twin is `DATA_MANIFEST.json`, which `gpe_analysis.Model.stress_frame()` reads. The shipped models predate the driver's `provenance.txt` stamp; this file is their provenance record.

**Stress frame** `massless` for every model: the stored `sigma_*` fields are the 2nd Piola–Kirchhoff stress of a run without gravity or lithostatic prestress (see README §4 for what the fields are).

## `idealized_beam`

- command: `julia --project=. model/idealized_beam.jl`
- origin: re-solved with the release code (model/<model>.jl) on 2026-09-14; identical to the archive copy (finite_strain/out/<model>, commit 10d0d39) on every shared field, plus the three Cauchy-gradient fields the archive export predated

| model | parameters | file | bytes | sha256 |
|---|---|---|---|---|
| `idealized_beam` | nondimensional=True, L=10.0, h=1.0, nx=200, nz=40, E=60.0, nu=0.25, sigma_Y=1.0, delta_max=2.2, nsteps=44, element_order=2, foundation=none, gravity=none | `gpe_model.vtu` | 825354 | `97b9b1de74f75a80e436f10e05661351297730ab4e3680a21eec307a7a7fb70f` |

## `idealized_beam_elastic`

- command: `julia --project=. model/idealized_beam.jl`
- origin: re-solved with the release code (model/idealized_beam.jl) on 2026-09-14; identical to the archive copy (finite_strain/out/<model>, commit 758e4dc) on every shared field, plus the three Cauchy-gradient fields the archive export predated

| model | parameters | file | bytes | sha256 |
|---|---|---|---|---|
| `idealized_beam_elastic` | nondimensional=True, L=10.0, h=1.0, nx=200, nz=40, E=60.0, nu=0.25, sigma_Y=1000000.0, delta_max=2.2, nsteps=44, element_order=2, foundation=none, gravity=none | `gpe_model.vtu` | 711416 | `0c31ab5cfff5de72a8923e428b70d777227be769706c5695126348a826415460` |

## `suite1_strength`

- command: `julia --project=. model/paper_models.jl suite1`
- origin: committed in the private development archive (ferrite_plate_flexure, HEAD 6e7fb08, 2026-09-11) at finite_strain/out/paper/set1_rheology/<model> (commit 0c6fe38, 2026-07-08)

| model | parameters | file | bytes | sha256 |
|---|---|---|---|---|
| `dd_vm_asym_60km_V4` | E_Pa=70000000000.0, nu=0.25, h_km=60, nx=800, nz=48, nsteps=22, L_km=1600, rho_a=3300, rho_w=1000, g=9.81, V_TN=4.0, strength=DD-VM asymmetric (depth-dependent von Mises), mu=0.15, cohesion_MPa=31.534, hardening_H_Pa=700000000.0 | `gpe_model.vtu` | 3748810 | `7420805d6da43078b842accb89949352cc2b9be29f9f832b7022bb9065ef4e27` |
|  |  | `gpe_topo.csv` | 21227 | `ade789c60e0b7543e35b3e6610afed6cd1fe8f5a41a02dd8c2bdd41e34cda621` |
| `dd_vm_sym_60km_V4` | E_Pa=70000000000.0, nu=0.25, h_km=60, nx=800, nz=48, nsteps=22, L_km=1600, rho_a=3300, rho_w=1000, g=9.81, V_TN=4.0, strength=DD-VM symmetric (depth-dependent von Mises), mu=0.35, cohesion_MPa=36.695, hardening_H_Pa=700000000.0 | `gpe_model.vtu` | 3774627 | `3402bedda8e2dbdad860f46b1b5bb04b377edba4c1c33a61cf38b3d64c472c16` |
|  |  | `gpe_topo.csv` | 21229 | `d04ffaa60f997f5e6fbf7fb689a22b92b91f89794bc098bebf161330d48bcc9a` |
| `elastic_deep_60km_V4` | E_Pa=70000000000.0, nu=0.25, h_km=60, nx=800, nz=48, nsteps=24, L_km=1600, rho_a=3300, rho_w=1000, g=9.81, V_TN=4.0, strength=elastic (sigma_Y -> inf: 1e12 Pa) | `gpe_model.vtu` | 3469414 | `4019714c780f917ab4bd80477cc0b4258ac362bfbeae196b505ae91c20f79bdf` |
|  |  | `gpe_topo.csv` | 21235 | `825ff471dd7ddd2df9a899e3ff3a35e1f8d60aadd7ecf506466d5667fa31efbd` |
| `tresca_deep_150_60km_V4` | E_Pa=70000000000.0, nu=0.25, h_km=60, nx=800, nz=48, nsteps=24, L_km=1600, rho_a=3300, rho_w=1000, g=9.81, V_TN=4.0, strength=Tresca, uniform, sigma_Y_MPa=150 | `gpe_model.vtu` | 3791732 | `7ac0001648233e54a5e69446a9116efe6d732455bf1704b9a2a41e0c16b5c9d2` |
|  |  | `gpe_topo.csv` | 21234 | `c04f1a5d6766a915f360e34f8867645784ab71a8776570c9cebe06719b93ea57` |

## `suite2_load`

- command: `julia --project=. model/paper_models.jl v_sweep`
- origin: on-disk output of `paper_models.jl v_sweep` in the working tree of the private development archive (ferrite_plate_flexure, HEAD 6e7fb08, 2026-09-11) at finite_strain/out/paper/set3_v_sweep/<model>; never committed there

| model | parameters | file | bytes | sha256 |
|---|---|---|---|---|
| `tresca_deep_150_60km_V1` | E_Pa=70000000000.0, nu=0.25, h_km=60, nx=800, nz=48, nsteps=24, L_km=1600, rho_a=3300, rho_w=1000, g=9.81, strength=Tresca, uniform, sigma_Y_MPa=150, V_TN=1.0 | `gpe_model.vtu` | 3467621 | `2f412b21ec9348f35c027e58d3f59d9e35da273eada7ea0a43417dc4bb52b8fb` |
|  |  | `gpe_topo.csv` | 21235 | `bd2b3fd9165e1ad84c927b2cda69389ebe8a4bb24e7e2cfbfff308b53d97dc9f` |
| `tresca_deep_150_60km_V1p5` | E_Pa=70000000000.0, nu=0.25, h_km=60, nx=800, nz=48, nsteps=24, L_km=1600, rho_a=3300, rho_w=1000, g=9.81, strength=Tresca, uniform, sigma_Y_MPa=150, V_TN=1.5 | `gpe_model.vtu` | 3468362 | `089ee380d1bc23b9ee72d99cbd4a4197da49d526f90b7247e40fec5226c71b26` |
|  |  | `gpe_topo.csv` | 21235 | `fb7c97cbb09c2586a24918cf097541761a5ccaa04dc66da151e724e5456a442c` |
| `tresca_deep_150_60km_V2` | E_Pa=70000000000.0, nu=0.25, h_km=60, nx=800, nz=48, nsteps=24, L_km=1600, rho_a=3300, rho_w=1000, g=9.81, strength=Tresca, uniform, sigma_Y_MPa=150, V_TN=2.0 | `gpe_model.vtu` | 3468300 | `96163d790f655d264eae9369c3a950d16b2537d211db7fbea96eddebea5502c4` |
|  |  | `gpe_topo.csv` | 21235 | `be2b49e0de285d968894187ad3a722b343891acdc5af4ea97618b64a83a8a121` |
| `tresca_deep_150_60km_V2p5` | E_Pa=70000000000.0, nu=0.25, h_km=60, nx=800, nz=48, nsteps=24, L_km=1600, rho_a=3300, rho_w=1000, g=9.81, strength=Tresca, uniform, sigma_Y_MPa=150, V_TN=2.5 | `gpe_model.vtu` | 3712105 | `6d2c3468fc04763e56c26e8da54e9418145f402f09e9edb3872504370d9bf926` |
|  |  | `gpe_topo.csv` | 21235 | `93a07a4cae32d83f9340f1811b6d441ccbe887c51ba3209e23bdd2e8f89bebbb` |
| `tresca_deep_150_60km_V3` | E_Pa=70000000000.0, nu=0.25, h_km=60, nx=800, nz=48, nsteps=24, L_km=1600, rho_a=3300, rho_w=1000, g=9.81, strength=Tresca, uniform, sigma_Y_MPa=150, V_TN=3.0 | `gpe_model.vtu` | 3754606 | `a5398ac5470b7fccbf55a3424a9539c7851fa50438a2831ee98569c3d86bf85b` |
|  |  | `gpe_topo.csv` | 21234 | `f8c8bbb2fd86eb5a8b8e6e730a92fbd010fd193918c67d719223506a32642b1f` |
| `tresca_deep_150_60km_V3p5` | E_Pa=70000000000.0, nu=0.25, h_km=60, nx=800, nz=48, nsteps=24, L_km=1600, rho_a=3300, rho_w=1000, g=9.81, strength=Tresca, uniform, sigma_Y_MPa=150, V_TN=3.5 | `gpe_model.vtu` | 3778567 | `5a774bdc52b0995f46bec0c887eadade9d7e3c9afdc278c1ae010063c3daca09` |
|  |  | `gpe_topo.csv` | 21234 | `7f1cff9782afc8ca7186c0da4428c7a4117d8094e562e0749dbffcb1b96af6f1` |
| `tresca_deep_150_60km_V4p5` | E_Pa=70000000000.0, nu=0.25, h_km=60, nx=800, nz=48, nsteps=24, L_km=1600, rho_a=3300, rho_w=1000, g=9.81, strength=Tresca, uniform, sigma_Y_MPa=150, V_TN=4.5 | `gpe_model.vtu` | 3795093 | `88a011b3402da7f7e34cb7a8c4634c119508de73b8e615e84851b3701b8c4c61` |
|  |  | `gpe_topo.csv` | 21233 | `2a35f2527cdb28717bb961b31f00b29238b6138dfb1e8ecd0e7149d998114671` |

## `suite3_background`

- command: `julia --project=. model/paper_models.jl nd_sweep`
- origin: on-disk output of `paper_models.jl nd_sweep` in the working tree of the private development archive (ferrite_plate_flexure, HEAD 6e7fb08, 2026-09-11) at finite_strain/out/paper/set2_nd_sweep/<model>; never committed there

| model | parameters | file | bytes | sha256 |
|---|---|---|---|---|
| `tresca_deep_150_60km_V4_mem-1` | E_Pa=70000000000.0, nu=0.25, h_km=60, nx=800, nz=48, nsteps=24, L_km=1600, rho_a=3300, rho_w=1000, g=9.81, strength=Tresca, uniform, sigma_Y_MPa=150, V_TN=4.0, N_mem_TN=-1.0 | `gpe_model.vtu` | 3714623 | `fdbb77376ec3b0447bf09033e48357ff84899e926bdb2589eede23135f0c00cc` |
|  |  | `gpe_topo.csv` | 21313 | `3ce64920d81f5e3374feb3ca80de462d9ddc061baebe2f7ace2ceb87a9efeb44` |
| `tresca_deep_150_60km_V4_mem-2` | E_Pa=70000000000.0, nu=0.25, h_km=60, nx=800, nz=48, nsteps=24, L_km=1600, rho_a=3300, rho_w=1000, g=9.81, strength=Tresca, uniform, sigma_Y_MPa=150, V_TN=4.0, N_mem_TN=-2.0 | `gpe_model.vtu` | 3692085 | `f4c4f38c7901111ece7d58e479bf28ca9d047e81af790d04c22d0172779dbc7d` |
|  |  | `gpe_topo.csv` | 21272 | `525b34110852c07536261cce86a57c217105005ca5e9985bc841f86dc0997b9c` |
| `tresca_deep_150_60km_V4_mem-3` | E_Pa=70000000000.0, nu=0.25, h_km=60, nx=800, nz=48, nsteps=24, L_km=1600, rho_a=3300, rho_w=1000, g=9.81, strength=Tresca, uniform, sigma_Y_MPa=150, V_TN=4.0, N_mem_TN=-3.0 | `gpe_model.vtu` | 3680722 | `bd454842e76ff7c76869ac07effc5c8eeb7c637d1314be8551b0bc604b801bcf` |
|  |  | `gpe_topo.csv` | 21256 | `e15c40ad5eb88532bda2409bcc4fc803fbceb7da5c3e402b35710424a5f59a92` |
| `tresca_deep_150_60km_V4_mem1` | E_Pa=70000000000.0, nu=0.25, h_km=60, nx=800, nz=48, nsteps=24, L_km=1600, rho_a=3300, rho_w=1000, g=9.81, strength=Tresca, uniform, sigma_Y_MPa=150, V_TN=4.0, N_mem_TN=1.0 | `gpe_model.vtu` | 3714364 | `dbe995296b718690f7e358ae81ebd945a43512469075f82a16dac12b2f0d1625` |
|  |  | `gpe_topo.csv` | 21198 | `921f234c6bb75bc1b26928ff180e02728021250a0e27145aa87a36cf01ef522f` |
| `tresca_deep_150_60km_V4_mem2` | E_Pa=70000000000.0, nu=0.25, h_km=60, nx=800, nz=48, nsteps=24, L_km=1600, rho_a=3300, rho_w=1000, g=9.81, strength=Tresca, uniform, sigma_Y_MPa=150, V_TN=4.0, N_mem_TN=2.0 | `gpe_model.vtu` | 3703141 | `91dd1b45819867344c7e6351bb372d62d70a1fc6b25a52c795b0e1c248580316` |
|  |  | `gpe_topo.csv` | 21199 | `bccd499338dd276b3bf8996fcda745abd2fb16f68b1cd94e992f1b1ebcb0545f` |
| `tresca_deep_150_60km_V4_mem3` | E_Pa=70000000000.0, nu=0.25, h_km=60, nx=800, nz=48, nsteps=24, L_km=1600, rho_a=3300, rho_w=1000, g=9.81, strength=Tresca, uniform, sigma_Y_MPa=150, V_TN=4.0, N_mem_TN=3.0 | `gpe_model.vtu` | 3687056 | `2792c39784fd55133b2c5489f23116cefc19dfa80e2a0e446e82869b42001418` |
|  |  | `gpe_topo.csv` | 21200 | `f8645deb576ed0a5cffd7ca93ed7d204aaf625fdec6c4df32561fdd607efeed2` |

## `suite4_thickness`

- command: `julia --project=. model/paper_models.jl thickness`
- origin: committed in the private development archive (ferrite_plate_flexure, HEAD 6e7fb08, 2026-09-11) at finite_strain/out/paper/set4_thickness/<model> (commit 26e6c9b, 2026-08-11)

| model | parameters | file | bytes | sha256 |
|---|---|---|---|---|
| `tresca_150_30km` | E_Pa=70000000000.0, nu=0.25, h_km=30, nx=800, nz=48, nsteps=24, L_km=1600, rho_a=3300, rho_w=1000, g=9.81, strength=Tresca, uniform, sigma_Y_MPa=150, w_target_m=3232.9, V0_TN=2.05, V_TN_tuned=2.0446 | `gpe_model.vtu` | 3847841 | `dee0b236fedace1065e1f2d851d35fc32241773829d5324961ad026024ff8324` |
|  |  | `gpe_topo.csv` | 21219 | `528d9c6dd15b6e850962321d170bf0099531546fb2b63437e930e2a7e7e4ec62` |
|  |  | `tuned_V.txt` | 91 | `533e6d844391eb8632cc75372fd619da950808ea510328b207d151ab36d11cd3` |
| `tresca_150_40km` | E_Pa=70000000000.0, nu=0.25, h_km=40, nx=800, nz=48, nsteps=24, L_km=1600, rho_a=3300, rho_w=1000, g=9.81, strength=Tresca, uniform, sigma_Y_MPa=150, w_target_m=3232.9, V0_TN=2.7, V_TN_tuned=2.7016 | `gpe_model.vtu` | 3824240 | `aea4ceada5d5e0d35f792a1194d587d5ba81f3e88f34c22b19fbae2c656735f9` |
|  |  | `gpe_topo.csv` | 21272 | `5ef60b79185155973011ca44f5d218df69e40200407ae07563514872835603f3` |
|  |  | `tuned_V.txt` | 91 | `cb5d380c2faac7dea5290260b5ee0ae516f9fccdf125d2da9452185696691368` |
| `tresca_150_50km` | E_Pa=70000000000.0, nu=0.25, h_km=50, nx=800, nz=48, nsteps=24, L_km=1600, rho_a=3300, rho_w=1000, g=9.81, strength=Tresca, uniform, sigma_Y_MPa=150, w_target_m=3232.9, V0_TN=3.35, V_TN_tuned=3.3541 | `gpe_model.vtu` | 3804351 | `130ccd753a0ae1028ad60c6a1a52d1fe2dbdef40cca8580a7d16427ec8e2f71f` |
|  |  | `gpe_topo.csv` | 21178 | `9b70dee8841177a20c75538a8373498856921781deabf625d39d610c3ebe920b` |
|  |  | `tuned_V.txt` | 91 | `5447118381b29a151ca8810294045aaf3aaf19fbf4043a7d18acdce8ea687595` |
