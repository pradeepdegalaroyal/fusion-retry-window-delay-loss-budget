# Method

## Fusion node model (analytic)

A fusion consumes two photons. A node may make up to R attempts inside a window of R time slots; attempt k
(0-indexed) is a usable outcome with probability

    q_k = p_s * eta_k^2,    eta_k = 10^(-(l0 + d * (k + c * (R - 1))) / 10),

where p_s is the single-attempt linear-optics fusion success probability, l0 is the base photon loss in dB (source,
coupling, detector, no waiting), d is the delay-line loss per slot in dB, k slots are spent waiting for the retry,
and c in [0, 1] is the clock coupling: the fraction of the window that every photon of the network idles because
the layers run on a synchronous clock. The fusion outcome is erased with probability e(R) = product over k of
(1 - q_k). Retries are modelled as emitter re-emissions (fresh photons), so a failed or lost attempt is heralded
and can be retried; loss and fusion failure are treated alike. Pauli errors from imperfect photon
indistinguishability are not modelled.

## Network and decoder (Monte Carlo)

The fusion network is an L x L x L periodic cubic lattice (L = 4, 6, 8, 10) whose edges are fusion outcomes. For a
stabiliser code with this geometry an erasure pattern is correctable exactly when the erased edges contain no
non-contractible cycle. A union-find structure that tracks each node's displacement in the infinite cover detects
such a winding cycle (checked in the tests against an independent breadth-first lift). The logical failure
probability P_L is the fraction of failing patterns. Correlated erasure is added through a site-outage share kappa:
a fraction kappa of the marginal edge erasure comes from whole-node outages (all six edges of a node erased) and the
rest from independent edge erasures; the marginal edge erasure is the same for every kappa.

## Baselines

- **Standard quantum baseline:** single-attempt fusion (R = 1), with p_s = 0.5 (bare linear-optics fusion) and
  p_s = 0.75 (boosted with a Bell-pair ancilla, 4 photons per attempt in the cost model). Photon-count costs are
  assumptions of this cost model.
- **Classical baseline:** a lattice-blind design rule that picks the window minimising the independent-erasure
  probability e(R) at the operating loss, without any lattice or decoder simulation; it is compared with the
  decoder-aware optimum.
- **Naive upper baseline:** unbounded retrying (largest window, R = 30).

## Metrics

Tolerable base photon loss l0* (dB per photon) such that e(R) is at or below the erasure level where P_L(L = 8) is
5 percent (target read from the Monte Carlo curve), also reported at the infinite-lattice crossing
(threshold) estimated from the crossing of the L = 8 and L = 10 curves; optimal window R* = argmax of l0*;
cost = photons per attempt times window slots; logical failure P_L.

## Ablation and robustness sweep

- **Ablation:** clock coupling c in {0, 0.5, 1}, removing the synchronous-clock penalty (c = 0).
- **Noise and robustness sweep:** base photon loss from 0 to 3 dB with delay loss 0.01 to 0.3 dB per slot, fusion
  success 0.5 to 0.95, and correlated site-outage share kappa in {0, 0.5, 1}.

## Seeds and sizes

5 seeds for every Monte Carlo curve (1500 trials per point and seed for the erasure sweeps, 800 for the direct
failure-versus-loss sweeps); lattice sizes L = 4, 6, 8, 10; reference lattice L = 8. Mean and standard deviation are
taken over seeds. The run is reproducible with `python experiments/run_all.py`.
