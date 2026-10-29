# Research gap: how long should a fusion node keep retrying?

All numbers in this study come from classical simulation on a laptop (statevector-free: the network is
treated at the level of heralded erasures, lattice sizes up to 10 x 10 x 10). Nothing here is a hardware
measurement and no quantum advantage is claimed.

## State of the art from the provided literature

Photonic quantum computing is reviewed as a leading route to large machines because photons are low-noise and
mobile, with loss management, photon indistinguishability and large-scale integration named as the main
bottlenecks (Slussarenko and Pryde, 2019, doi:10.1063/1.5115814; Nusyura et al., 2026,
doi:10.1007/s42484-026-00402-8; Elshaari et al., 2020, doi:10.1038/s41566-020-0609-x). Noisy intermediate-scale
quantum (NISQ) devices motivate loss-aware thinking (Preskill, 2018, doi:10.22331/q-2018-08-06-79), while
fault-tolerant machines rely on error-correcting codes whose erasure tolerance depends on the code geometry
(Terhal, 2015, doi:10.1103/RevModPhys.87.307). The most recent architecture paper in the brief, by Chen et al.
(2026, doi:10.1038/s41534-026-01233-y), builds quantum low-density parity-check (qLDPC) codes from fusions of
emitter-generated resource states and uses a conditional repeat-until-success (RUS) strategy to reach high
photon-loss tolerance. At the system level, Yousif et al. (2027, doi:10.1007/978-3-032-39708-9_27) show that
facet loss and detection efficiency set the scalability of photonic quantum clouds and that tail latency matters,
and the low-loss integrated platforms of Wang et al. (2018, doi:10.1038/s41586-018-0551-y) and the programmable
circuits surveyed by Wang et al. (2020, doi:10.1038/s41566-019-0532-1) set the loss scale per component.

## Related work

| Paper | Year | Role | What it does | What it does NOT cover |
|---|---|---|---|---|
| Chen et al., doi:10.1038/s41534-026-01233-y | 2026 | recent | Emitter-based fusion architecture for qLDPC codes with RUS fusions and loss tolerance | Does not treat the retry window length as a design variable priced by delay-line loss |
| Yousif et al., doi:10.1007/978-3-032-39708-9_27 | 2027 | recent | Maps device loss and detection efficiency to cloud throughput and latency | No fusion-level or code-level loss budget |
| Di Bartolo et al., doi:10.1038/s41534-026-01236-9 | 2026 | recent | Photonic quantum reservoir: indistinguishable photons help forecasting | Not about fault tolerance or fusion loss |
| Nusyura et al., doi:10.1007/s42484-026-00402-8 | 2026 | recent | Review of photonic quantum information, names loss management as a bottleneck | Qualitative, no budget |
| Burger et al., doi:10.1038/s41467-026-75489-5 | 2026 | recent | Emitters in photonic-crystal waveguides with inhibited decay | Source physics only, no network-level tolerance |
| Preskill, doi:10.22331/q-2018-08-06-79 | 2018 | landmark | Frames noise-limited NISQ devices | No photonic fusion model |
| Terhal, doi:10.1103/RevModPhys.87.307 | 2015 | landmark | Stabiliser codes, thresholds, topological codes | No linear-optics fusion failure or retry cost |
| Wang et al., doi:10.1038/s41566-019-0532-1 | 2020 | landmark | Integrated photonic quantum technologies | No retry-window budget |
| Wang et al., doi:10.1038/s41586-018-0551-y | 2018 | landmark | Modulators below 0.5 dB on-chip loss, enabling fast feed-forward | Component level only |
| Zhong et al., doi:10.1126/science.abe8770 | 2020 | landmark | Photonic Gaussian boson sampling with high loss control | Sampling, not error-corrected computing |
| Slussarenko and Pryde, doi:10.1063/1.5115814 | 2019 | review | Concise review of photonic quantum computing | Open problems stated, not quantified |
| Elshaari et al., doi:10.1038/s41566-020-0609-x | 2020 | review | Hybrid integrated photonic circuits and their loss and connectivity trade-offs | No architecture-level optimisation |
| Cai et al., doi:10.1016/j.fmre.2020.12.006 | 2021 | review | Bosonic error correction reaching break-even | A different (microwave) platform, no photon retry |

## The gap

Recent work shows that repeating a failed fusion buys loss tolerance (Chen et al., 2026), and reviews list loss
as the main bottleneck, but the cost side of retrying is rarely priced: every extra attempt keeps photons in delay
lines, and if the fusion network runs on a synchronous clock, a longer retry window makes every photon wait
longer. Taken together, retries lower the fusion erasure probability while the window raises per-photon loss.
Whether the window should be long, short or unbounded is a design question that nobody can answer from the
qualitative statements in the reviews.

**Novelty in one sentence.** This study treats the length of the repeat-until-success window as an explicit design
variable, prices it with a synchronous-clock delay-loss model, and uses a Monte Carlo erasure-correctability
criterion on a periodic cubic fusion network to turn it into a tolerable-photon-loss budget and an optimal window.

## Falsifiable hypotheses and success metrics

- **H1 (finite optimum).** With a synchronous clock (clock coupling c = 1) and fusion success 0.75, there is a delay
  loss per slot at which the loss-budget-optimal window R* is finite and unbounded retrying cannot reach the
  target, whereas with no clock penalty (c = 0) the budget grows with every extra attempt. Falsified if the
  budget is monotone in the window at c = 1, or not monotone at c = 0.
- **H2 (standard baselines).** Single-attempt fusion (standard linear-optics success 0.5 and boosted success 0.75)
  cannot reach a 5 percent logical failure at lattice size L = 8 even with lossless photons, so retries are not an
  optimisation but a requirement in this proxy. Falsified if the budget at R = 1 is positive.
- **H3 (decoder awareness is not free).** Correlated site outages (all fusions at a node erased together) lower the
  tolerable marginal erasure probability, and hence the loss budget, by at least 0.3 dB at fusion success 0.75 and
  delay loss 0.1 dB per slot, relative to independent erasure. Falsified if the reduction is smaller.

**Success metric.** Tolerable base photon loss (dB per photon) at which the logical failure probability of the
L = 8 network stays at or below 5 percent, reported as mean and standard deviation over 5 seeds, plus the window that
maximises it.
