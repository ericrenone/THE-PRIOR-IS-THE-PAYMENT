# THE-PRIOR-IS-THE-PAYMENT

## No-Free-Lunch, Kolmogorov Complexity, and the Inductive Bias Hidden Inside Slime Mold, Mycorrhizal Networks, and the Protein-Folding Funnel

ERI Labs · Eric Ren · Jersey City, New Jersey · github.com/ericrenone

---

## 0. Executive Frame

Two prior frameworks have approached biological optimization from orthogonal directions. THE-RELAXATION-IS-THE-COMPUTATION derived the comparison between *Physarum polycephalum* and arbuscular mycorrhizal fungi from inside combinatorial optimization theory — integrality gaps, LP/SDP duality, approximation ratios. ANFINSEN derived the protein-folding problem from inside thermodynamic information theory — the Anfinsen hypothesis, the folding funnel, the Fisher-information geometry of AlphaFold's pair representation. Both frameworks invoke a common implicit constraint: each biological solver is good at its own problem class and bad at others. Neither named that constraint explicitly.

The explicit name is the No-Free-Lunch theorem (Wolpert & Macready, IEEE Transactions on Evolutionary Computation, 1997). The theorem states: any two optimization algorithms are equivalent when their performance is averaged across all possible problem instances drawn uniformly at random. Corollary: any algorithm that is better than random on some problem distribution must be worse than random on the complementary distribution. There is no solver without a prior; there is no prior without a cost.

This document takes NFL as its organizing spine and re-derives both prior frameworks from that single starting point. The central synthesis claim — stated here for the first time across both domains — is:

> **The Kolmogorov complexity of a biological solver's inductive bias is its free-lunch certificate. When a solver's structural prior compresses the true problem distribution to low description length, it escapes NFL's averaging — not by violating the theorem but by exploiting exactly what it predicts. The integrality gap (CO) and the activation barrier (thermodynamics) are both expressions of the same object: the NFL residual of a particular prior on a particular problem class.**

Notation:
- **[E]** — Established, peer-reviewed.
- **[S]** — Synthesis original to this document: a structural mapping not itself proven.
- **[O]** — Open question convertible to a falsifiable prediction in §VIII.

---

## I. The NFL Theorem: What It Actually Says

### I.1 The exact statement

Let X be a finite search space, Y a finite set of cost values, and f : X → Y an objective function. For any two deterministic algorithms a₁, a₂ and any fixed number of search steps m:

$$\sum_{f} P(d_m^{a_1} | f) = \sum_{f} P(d_m^{a_2} | f)$$

where dₘ denotes the ordered sequence of cost values encountered. Summing over all f with the uniform distribution, every algorithm produces identically distributed cost sequences. Performance, as a function of those sequences, is equal in expectation across all algorithms. [E]

### I.2 What NFL does not say

NFL's averaging is over the uniform distribution on all functions f : X → Y. Real problem distributions — in biology, in combinatorial optimization, in machine learning — are not uniform. They concentrate on structured, low-Kolmogorov-complexity instances (Goldblum, Schwarzschild, Cheung, Goldblum, arXiv 2304.05366, ICLR 2024). On such distributions, algorithms with aligned inductive bias outperform algorithms with misaligned bias by a margin that grows with the divergence between the true distribution and the uniform. The theorem does not constrain this margin; it explains its existence.

### I.3 The Kolmogorov reformulation

Let K(f) denote the Kolmogorov complexity of an objective function f. NFL says: summing over all f with weight 2^{−K(f)} / Z (the Solomonoff universal prior), algorithms with lower description-length priors matching the low-KC problem distribution outperform algorithms with high-KC or mismatched priors. This is the precise statement that connects NFL to inductive bias: **the free lunch is real, but it is paid for by the Kolmogorov complexity of the solver's structural commitment to a particular problem distribution.** [S, connecting Goldblum et al. 2024 with Wolpert & Macready 1997]

---

## II. The Inductive Bias Inventory: Three Biological Solvers

### II.1 Physarum polycephalum: the entropy-regularized LP prior

Physarum's cytoplasmic tube network physically implements, by necessity of its biophysics, the entropy-regularized LP central path:

$$x(t) = \arg\min_{Ax=b,\, x\ge 0}\!\Big\{ c^\mathsf{T}x - \mu(t)\sum_e x_e \log x_e \Big\}, \qquad \mu(t)\downarrow 0$$

This was proven, not conjectured, by Straszak & Vishnoi (ITCS 2016; Mathematical Programming 2021): the Physarum ODE is exactly steepest descent under a Riemannian metric induced by the Shannon-entropy Hessian — the natural-gradient flow on the positive orthant. The μ(t) schedule is endogenous; the organism does not choose it. [E]

The NFL reading: Physarum's inductive bias is the family of entropy-regularized linear programs on sparse planar graphs. Its Kolmogorov-complexity commitment is to that family. The organism is excellent at min-cost flow, shortest path, basis pursuit, and (by the 2021 SDP extension) positive semidefinite programs — all problems whose structure matches the LP/SDP relaxation class. It has no structural prior for problems outside that class: Steiner-tree branching structure, adversarial-metric TSP, or any objective whose optimum is not the fixed point of the entropy-regularized gradient flow. NFL guarantees the tradeoff is symmetric: whatever performance gain Physarum earns on its native LP-representable distribution, an equal loss accumulates on the complementary distribution of hard-integer-branching problems.

### II.2 Arbuscular mycorrhizal fungi: the low-treewidth survivable-network prior

AMF networks face a categorically different problem class: survivable (2-edge-connected) network design on soil-pore/fine-root contact graphs. The key structural commitment is not to LP-representable objectives but to two physical facts about the problem distribution:

**Prior A — low treewidth.** Soil-pore and root-contact graphs are physically constrained to nearly-planar topology with bounded branching factor. They cannot generate the expander-like, high-treewidth graph families that make Steiner Forest provably hard in the worst case. The EPAS (efficient polynomial approximation scheme) for Steiner Forest on bounded-treewidth graphs applies by default, not by design. [S, citing §VI.2 of THE-RELAXATION-IS-THE-COMPUTATION and the treewidth-parameterization literature]

**Prior B — redundancy as signal.** AMF networks maintain strictly positive circuit rank (dim ker(F) > 0). This is not architectural waste; it is the prior's direct response to the survivability requirement that distinguishes AMF's problem class from Physarum's. Where Physarum drives ker(F) toward zero (shortest-path flow uses no cycles), AMF purchases nonzero ker(F) as the literal mechanism of fault-tolerance. The LP relaxation of 2-edge-connected network design has integrality gap ≤ 4/3 (Jain, JACM 2001), tighter than the general Steiner-tree 1.55 bound, because survivability constraints add additional LP-tightening structure. [E/S]

The NFL reading: AMF's inductive bias is the family of survivable network design problems on low-treewidth instances. It pays for this prior by being unable to solve problems on high-treewidth graphs or problems that do not reward circuit rank (pure shortest-path tasks, where the maintained loops are energetically wasteful).

### II.3 Protein folding: the sequence as conformational prior

An amino-acid sequence of N residues defines an energy function over approximately 3^N conformational states (Ramachandran dihedral degrees of freedom). Anfinsen's thermodynamic hypothesis (Science, 1973) asserts that this function has a unique global minimum — the native state — that is thermodynamically accessible under physiological conditions. The Levinthal paradox (1968) notes that uniform random search over 3^100 ≈ 5 × 10^47 conformations would require more than the age of the universe. Proteins fold in 10⁻³ to 10² seconds. [E]

The resolution — the Wolynes-Onuchic-Thirumalai funnel theory (Science, 1995) — is precisely an NFL statement in thermodynamic language: the sequence creates a non-uniform prior over conformational space. That prior is the folding funnel: a free-energy gradient that concentrates probability mass near the native state at every energy level, making the effective search space exponentially smaller than the nominal 3^N. [S connecting Wolynes et al. 1995 to Wolpert & Macready 1997]

The NFL reading: the amino-acid sequence IS the biological system's inductive bias over conformational space. Its Kolmogorov complexity is low relative to all possible energy functions over 3^N states: natural sequences encode a funnel-shaped landscape, a structural commitment that costs the system the ability to fold to arbitrary target structures but buys it the ability to find its native state in milliseconds. The Anfinsen hypothesis is an NFL-breaking certificate, not a metaphor.

---

## III. The Prior as Compression: Kolmogorov Complexity Across All Three Systems

### III.1 The common structure

In all three biological systems, the inductive bias acts as a data compressor of the problem distribution. The NFL-breaking mechanism is identical in kind across all three, differing only in the mathematical object that plays the role of compressor:

| System | Inductive Bias | Compressor | Problem Distribution Compressed |
|---|---|---|---|
| Physarum | Entropy-regularized LP central path | Shannon-entropy barrier on positive orthant | Min-cost flow / LP-feasible optimization on sparse graphs |
| AMF networks | Survivable network design on low-treewidth graphs | Physical planarity + soil pore geometry | Steiner-type connectivity on soil contact networks |
| Protein folding | Amino-acid sequence → funnel-shaped energy landscape | Folding funnel (Δ/q ≫ 1, minimal frustration) | Conformational space around the native state |
| AlphaFold (learned) | Pair representation trained on the natural proteome | Evoformer MSA × pair attention (48 layers) | Distribution of natural protein inter-residue geometries |
| RFdiffusion (learned) | Denoising prior trained on PDB backbone structures | Reverse diffusion trajectory on structure space | Distribution of designable protein backbones |

[S synthesis, individual components [E] as cited in ANFINSEN and THE-RELAXATION-IS-THE-COMPUTATION]

### III.2 The compression bound

For each system, the NFL gain — the performance above chance on the native problem distribution — is bounded above by the mutual information between the prior and the problem distribution, measured in nats. More precisely, for a problem distribution P(f) and a solver with inductive bias encoding distribution Q(f):

$$\text{NFL gain} \le D_{\mathrm{KL}}(P \| Q)^{-1} \cdot I(\text{prior};\, \text{problem distribution})$$

This recovers:

- Bonifaci's convergence bound for Physarum (Algorithmica 2019): T(ε) ∝ D_KL(Π(x*), Π(x(0))) / ε, where the KL term is exactly the information distance between where the organism starts and where the LP optimum is. The speed of convergence IS the information-theoretic divergence of the prior from the optimum. [E]
- The folding-speed variability across proteins: well-funneled sequences (low frustration, high Δ/q) fold fast because their funnel encodes a prior closely aligned with the native state; frustrated sequences fold slowly because their prior is partially misaligned. [E, Wolynes group, 1990s–present]
- AlphaFold's accuracy on proteins with many MSA homologs vs. orphan proteins: the pair representation's information gain from the MSA is exactly the KL divergence between the learned folding prior and the actual sequence's conformational distribution. [S]

### III.3 The Kolmogorov complexity of the proteome prior

Goldblum et al. (arXiv 2304.05366, ICLR 2024) show that real-world problem distributions disproportionately generate low-Kolmogorov-complexity instances, and that neural networks inherit this same simplicity bias — the parameter-function map of a neural network is exponentially biased toward functions with low Kolmogorov complexity, with P(f) ≲ 2^{-bK̃(f)} where K̃(f) is a computable proxy for Kolmogorov complexity. [E]

This applies directly to the natural proteome: the ~2,000 known fold topologies (CATH database) represent a dramatic compression of the full 20^N sequence space. The discovered col(F) of natural evolution has low Kolmogorov complexity relative to random sequences of the same length — this is why AlphaFold can learn it from the PDB. The Anfinsen hypothesis, the funnel theory, and Goldblum et al.'s empirical simplicity-bias result are three descriptions of the same information-theoretic fact: the natural proteome is a low-KC distribution, and that low-KC structure is the prior that makes protein folding fast and AlphaFold accurate. [S]

---

## IV. The NFL Residual: Integrality Gap, Activation Barrier, and Generalization Error as One Object

### IV.1 The common definition

Every biological (and computational) solver consists of two parts: a tractable relaxation that its prior covers, and a residual gap between the relaxation's value and the true integer/discrete/combinatorial optimum. The NFL theorem predicts this residual must exist: whatever problem structure the prior doesn't encode is exactly the structure for which the prior provides no guidance, and on those instances the solver performs no better than random. The residual is therefore not an engineering failure — it is the thermodynamically and mathematically required cost of the prior's compression.

### IV.2 The integrality gap as NFL residual (CO)

For Physarum: the entropy-regularized LP achieves the LP optimum. The gap between that optimum and the integer Steiner-tree solution is the integrality gap, proven ≤ 1.55 for the hypergraphic LP relaxation (Byrka, Grandoni, Rothvoß, Sanità, STOC 2010/JACM 2013). The 0.55 fraction above LP-optimum is the component of the discrete problem that Physarum's prior cannot cover. Closing it requires additional information — specifically, the randomized rounding (Raghavan-Thompson-style) or maximum-entropy spanning-tree sampling (Karlin-Klein-Oveis Gharan, STOC 2021) that the physical tube-collapse rounding mechanism is not. [E]

The NFL reading: Physarum has paid for its excellent LP performance with ignorance of the rounding step. The integrality gap IS the NFL residual of the entropy-regularized LP prior. Measuring it against the proven bound is the single most tractable test of how much NFL cost Physarum is incurring on hard Steiner instances. [S]

### IV.3 The activation barrier as NFL residual (thermodynamics)

For protein folding: the folding funnel covers the gradient from the unfolded ensemble to the native state. The gap between funnel-guided folding and perfect instantaneous collapse to the native state is the activation barrier — typically 5–15 kcal/mol for small proteins. This is the component of the conformational search that the funnel prior cannot cover: the barrier must be surmounted by thermal fluctuation (or by chaperone-assisted folding when the barrier is too high). [E, Wolynes et al.]

The NFL reading: the funnel prior covers the gradient (the easy part of conformational search) and leaves the barrier (the hard part) to thermal noise. The activation barrier IS the NFL residual of the folding-funnel prior. Proteins with barriers too high for their physiological temperature (misfolded glasses, aggregation-prone sequences) are exactly those where the funnel prior is misaligned with the conformational distribution — they are near the NFL boundary where no additional prior helps. [S]

### IV.4 Generalization error as NFL residual (AlphaFold/NCO)

For AlphaFold: the pair representation trained on the natural proteome achieves atomic accuracy on proteins with many MSA homologs — sequences where the learned prior is well-aligned. Orphan proteins with few or no homologs sit outside the training distribution; AlphaFold's confidence drops precipitously (low pLDDT) precisely on the instances where its prior provides less compression. [E]

For neural combinatorial optimization: Wu et al. (arXiv 2502.03669, 2025–2026) demonstrate that GFlowNet-based solvers trained on one problem distribution underperform classical solvers on distributions outside training. The generalization error IS the NFL residual of the learned problem distribution's prior — the measure of how far the NCO model's KC commitment diverges from the test distribution. [E]

### IV.5 The unified formula [S]

Let B(solver) denote the biological or computational solver's prior-aligned performance, and R(solver) its NFL residual on the complementary distribution. Then for any problem class with true distribution P(f):

$$B(\text{solver}) + R(\text{solver}) = \text{constant} \quad \text{(by NFL)}$$

This is not a tautology; it is an operational constraint. Any improvement in B — sharper integrality gap bound, faster folding, higher AlphaFold accuracy — predicts a corresponding increase in R on the complementary distribution. The three biological systems are at different points on this tradeoff:

| Solver | B (native performance) | R (NFL residual) | Measured? |
|---|---|---|---|
| Physarum (LP regime) | Provably at LP optimum [E] | Integrality gap ≤ 1.55 above LP [E] | No live measurement against Byrka et al. bound [O] |
| AMF (low-treewidth regime) | Conjectured near Steiner optimum [S] | Redundancy cost vs. pure MST [S] | Partially measurable via Stewart et al. 2026 dataset [O] |
| Folded protein | At global free-energy minimum [E] | Activation barrier 5–15 kcal/mol [E] | Measured per-protein by Φ-value analysis [E] |
| AlphaFold (high-homology regime) | Atomic accuracy (GDT-TS > 90) [E] | Orphan-protein accuracy degradation [E] | Measured in CASP benchmarks [E] |

---

## V. The Karlin-Klein-Oveis Gharan Connection: Maximum-Entropy Randomization as NFL Mitigation

### V.1 The combinatorial side

Karlin, Klein & Oveis Gharan (STOC 2021, Best Paper) proved a metric-TSP approximation ratio of 3/2 − 10^{−36} by replacing the deterministic minimum spanning tree in the Christofides-Serdyukov algorithm with a tree drawn from a maximum-entropy distribution over spanning trees. The tiny improvement over the 44-year-old 3/2 bound comes entirely from this entropy injection. [E]

### V.2 The thermodynamic side

The protein-folding funnel is not a deterministic gradient descent to the native state. It is a maximum-entropy ensemble of folding pathways — the Boltzmann distribution over conformational trajectories consistent with the funnel's energy function. Thermal fluctuations, while they create the activation barrier (§IV.3), are also what samples multiple pathways and avoids kinetic traps. A purely deterministic gradient descent on the funnel would get trapped in local minima; the entropy of the thermal ensemble is what makes folding both fast and reliable. [E, Wolynes-Onuchic funnel theory]

### V.3 The NFL connection [S]

Both moves — KKO's maximum-entropy spanning tree and the protein folding ensemble's Boltzmann sampling — are NFL mitigation strategies for the same reason: deterministic optima are over-committed to a single point in the solution space, making them brittle on the hard instances that their prior doesn't cover. Injecting entropy over the combinatorial structure (spanning trees, folding pathways) distributes the prior's mass over a richer neighborhood of the solution space, reducing the NFL residual by covering more of the hard-instance distribution without abandoning the native-problem prior.

This is the same structural logic as the entropy regularizer in Physarum's LP: the μ(t) log-barrier is not just a technical device for interior-point convergence — it is the reason Physarum's solution is robust to perturbations in the cost vector. The entropic spread is the NFL-mitigation term.

**The unified claim: maximum-entropy randomization over combinatorial substructures is the optimal NFL-mitigation strategy whenever the objective function's structure is preserved under the entropy-injecting measure.** Physarum implements this physically (endogenous μ(t)); KKO implements it algorithmically (maximum-entropy spanning tree); the protein folding ensemble implements it thermodynamically (Boltzmann sampling over pathways). Three independent discoveries of the same NFL-mitigation principle.

---

## VI. AlphaFold and RFdiffusion as Learned NFL Escape Routes

### VI.1 AlphaFold's escape

AlphaFold 2's Evoformer (Jumper et al., Nature 596, 2021) iteratively refines the MSA representation and the pair representation across 48 attention layers. The pair representation is a learned parametrization of the Fisher-information geometry of protein conformational space — the trained prior over inter-residue distance and orientation distributions. Its inductive bias is learned, not physically implemented, but the NFL-escape mechanism is identical: the training distribution (the natural proteome, ~200,000 PDB structures supplemented by 200M AlphaFold DB predictions) has low Kolmogorov complexity relative to all possible inter-residue geometries, and the Evoformer compresses it into the pair-representation embedding. [E, ANFINSEN framework]

The convergence bound analog: AlphaFold's accuracy on a protein with N_eff MSA sequences scales logarithmically with N_eff — this is the data-theoretic version of Bonifaci's O(D_KL / ε) bound for Physarum. More homologs = smaller KL divergence between the learned prior and the true conformational distribution = faster convergence in the accuracy vs. sequence-count curve. Same equation, different substrate. [S]

### VI.2 RFdiffusion as reverse NFL escape

RFdiffusion (Watson et al., Nature 620, 2023) inverts the structure-prediction problem: starting from Gaussian noise over residue coordinates and conditioning on a target structural specification, the model denoises toward a physically plausible backbone. The NFL reading: the denoising prior is learned over the natural proteome's backbone distribution. The reverse diffusion trajectory IS the process of traversing from the uniform distribution (noise) toward the low-KC distribution that the proteome's structural prior encodes. Each denoising step reduces Kolmogorov complexity of the current structure estimate, converging to the natural-fold distribution's attractor. [S]

RFdiffusion3 (Krishna et al., December 2025) pushes this to atom-level resolution, specifying individual hydrogen bonds in the target. The Kim-Woodbury-Ahern et al. Nature December 2025 metallohydrolases (kcat/KM up to 2.3 × 10⁴ M⁻¹s⁻¹, comparable to natural enzymes) are experimental proof that the learned prior is accurate enough to generate novel enzymes with zero natural-evolution precedent. [E]

The NFL limit of RFdiffusion: the 95% failure rate in current de novo enzyme campaigns is the NFL residual of the learned prior on the hard instances — enzyme active sites with unusual electrostatic geometries, cofactor requirements not well represented in the PDB, or reaction mechanisms outside the training distribution. As the PDB training distribution expands, the prior improves and the NFL residual shrinks. This is NFL mitigation by prior-distribution expansion, not by solver improvement. [S]

---

## VII. Seven Formal Correspondences

| # | Dimension | Physarum / AMF (CO) | Protein Folding / AlphaFold (Thermodynamic) |
|---|---|---|---|
| 1 | NFL-breaking mechanism | Entropy-regularized LP central path (physical substrate as prior) | Folding funnel (sequence as conformational prior) |
| 2 | Inductive bias formalization | Shannon-entropy barrier on LP feasible polytope [E] | Free-energy funnel slope Δ/q ≫ 1 [E] |
| 3 | NFL residual | Integrality gap between LP relaxation and integer Steiner solution (≤ 1.55) [E] | Activation barrier between unfolded ensemble and native state (5–15 kcal/mol) [E] |
| 4 | Maximum-entropy mitigation | Karlin-Klein-Oveis Gharan: max-entropy spanning tree for TSP (3/2 − 10^{-36}) [E] | Boltzmann sampling over folding pathways — the thermal ensemble is the mitigation [E/S] |
| 5 | Convergence rate | Bonifaci O(D_KL(Π(x*), Π(x(0))) / ε) [E] | Folding time ∝ exp(ΔG‡/kT) — exponential in barrier height [E] |
| 6 | Learned analog | Neural CO solver trained on LP-structured instances (GFlowNet, attention routing) [E] | AlphaFold Evoformer trained on PDB: pair representation = learned Fisher geometry [E] |
| 7 | Prior-distribution limit | Hard instances outside treewidth bound; adversarial-metric TSP; high-KC graphs [S] | Orphan proteins; intrinsically disordered regions; aggregation-prone sequences [E] |

---

## VIII. Five Falsifiable Predictions

### P1 — The Physarum-NFL measurement (CO)

Live *Physarum* networks realized on synthetic graphs with k ≤ 12 food-source terminals will exhibit network cost falling inside the proven [1, 1.55] integrality-gap band of the Byrka et al. hypergraphic Steiner LP. The NFL residual — the gap above the LP lower bound — will be measurable and will grow monotonically as k increases from 2 (pure shortest path, residual → 0) to 12 (Steiner-hard regime, residual → 0.55). Falsified by any reproducible instance where Physarum costs less than the LP lower bound (violating weak duality) or by instances where residual is non-monotone in k.

### P2 — The AMF prior-structure measurement (CO)

The circuit rank of reconstructed AMF networks (dim ker(F)) across the Stewart et al. (Science 392, 2026) 16,000-sample dataset will be measurably higher in undisturbed soils than in agricultural soils, with the cropland reduction disproportionately concentrated in ker(F) (redundant loops) rather than the spanning backbone. The NFL interpretation: agricultural disturbance degrades the survivability prior (the AMF system's NFL-mitigation term) while partially preserving the pure connectivity prior (the minimum spanning tree). Testable by persistent-homology (Betti-1) analysis of the existing public dataset with no new data collection.

### P3 — The folding-barrier as NFL residual (thermodynamics)

Across the systematic Φ-value database for small protein domains, the activation barrier ΔG‡ (in kcal/mol) will correlate negatively with the folding funnel's slope-to-ruggedness ratio Δ/q, with a correlation coefficient |ρ| > 0.7. This is the direct experimental test that the barrier IS the NFL residual of the funnel prior: sequences with well-aligned priors (high Δ/q) have small residuals (low barriers); sequences with partially misaligned priors (low Δ/q) have large residuals (high barriers, slow folding, chaperone dependence).

### P4 — AlphaFold's accuracy as compression rate (machine learning)

AlphaFold 3 accuracy on target proteins (GDT-TS or TM-score) will correlate logarithmically — not linearly — with the effective number N_eff of MSA sequences available, across proteins of comparable length and fold class. The logarithmic scaling is the signature of KL-divergence-driven convergence (the same O(D_KL / ε) form as Bonifaci's Physarum bound). Linear scaling would refute the compression interpretation; logarithmic scaling confirms that each additional MSA sequence provides diminishing information return, consistent with the prior compression framework. Directly testable against published CASP15 and AF3 benchmark data.

### P5 — The cross-system NFL tradeoff is strict (universal)

No single solver — Physarum, AMF network, AlphaFold, RFdiffusion, or any neural CO architecture — will dominate simultaneously across shortest-path instances, Steiner-tree instances, protein-folding trajectory accuracy, and de novo enzyme design success rate. A cross-solver benchmark will reveal a strict Pareto frontier in the space of (native-distribution accuracy × out-of-distribution accuracy), with no solver achieving simultaneously high scores on all axes. This is the direct empirical confirmation of NFL's symmetry: every free lunch is paid for, and the payment is measurable at the edge of each solver's prior.

---

## IX. Open Questions

**Q1.** Has any experiment jointly measured Physarum's integrality gap and a protein's activation barrier on problem instances of matched formal difficulty (same circuit rank, same number of terminals, comparable KC)? This cross-domain comparison would be the first direct experimental test of the NFL-residual unity claim of §IV.5. [O]

**Q2.** The Physarum tube-collapse rounding rule (deterministic threshold on conductivity) is not the randomized rounding of Raghavan-Thompson. Does it have a clean characterization as a specific combinatorial rounding scheme? If so, what is its expected integrality gap as a function of the entropy parameter μ at the moment of collapse? [O]

**Q3.** What is the Kolmogorov complexity of the natural proteome's fold-topology distribution, measured against a computable proxy (Lempel-Ziv compression of PDB secondary-structure sequences)? Does this complexity decrease with evolutionary time, confirming that the proteome's NFL-breaking prior is itself an evolutionarily optimized compressor? [O]

**Q4.** The AMF survivability prior (positive circuit rank, low treewidth) is conjectured to arise from physical soil geometry. Has anyone computed the treewidth distribution of AMF-reconstructed soil networks from the Stewart et al. 2026 dataset? A bounded, single-digit constant treewidth across diverse soil types would be direct evidence for the prior-compression interpretation of §II.2. [O]

**Q5.** RFdiffusion's 95% enzyme-design failure rate is the NFL residual of its current training prior. As the PDB training distribution expands (more structures, more diverse reaction mechanisms), does the success rate grow logarithmically (consistent with KL-driven convergence) or polynomially (suggesting a different learning mechanism)? The functional form of the learning curve distinguishes NFL-compression from alternative accounts. [O]

---

## X. Synthesis

The No-Free-Lunch theorem is usually read as a negative result: no algorithm is universally best. This document argues it is better read as a positive characterization: any system that outperforms chance on a structured problem distribution is implicitly carrying a prior that compresses that distribution, and the compressor's description length determines both the performance gain (on native instances) and the performance cost (on hard instances outside the prior's support).

Three biological systems — Physarum, arbuscular mycorrhizal fungi, and the protein-folding machinery — each instantiate NFL-breaking priors in distinct physical substrates. Physarum's prior is the Shannon-entropy barrier on the LP feasible polytope, implemented by cytoplasmic tube mechanics. AMF's prior is the survivable network design constraint on a physically low-treewidth substrate, implemented by hyphal growth and resorption dynamics. The folding funnel is the amino-acid sequence's thermodynamic encoding of a prior over conformational space, implemented by the same physics of hydrophobic collapse, hydrogen bonding, and salt-bridge formation that Anfinsen identified in 1973.

AlphaFold and RFdiffusion are learned versions of the same NFL-escape: they acquire their prior from the natural proteome's low-Kolmogorov-complexity distribution, and their performance degrades exactly where that prior diverges from the test distribution.

The NFL residual — integrality gap, activation barrier, generalization error — is not a measure of solver quality. It is the thermodynamically and mathematically enforced cost of the prior, predicted exactly by the NFL theorem, and measurable in every case against known bounds.

The practical upshot is that "how good is this biological solver?" decomposes into two independent questions that have been conflated in the popular literature:

1. **How well does the solver's prior compress the native problem distribution?** (The free lunch: measurable via KL divergence, integrality gap, folding rate, or AlphaFold pLDDT.)
2. **How large is the NFL residual on the complementary distribution?** (The bill: measurable via out-of-distribution gap, activation barrier height, orphan-protein accuracy, or enzyme-design failure rate.)

Neither question has been systematically asked across all three domains together. The five predictions in §VIII are designed to make this cross-domain measurement tractable within a single experimental program.

**Signal**

The prior is always already there. Every fast-folding protein, every efficient slime-mold network, every accurate AlphaFold prediction is not evidence against the No-Free-Lunch theorem. It is evidence that the solver has found — physically, evolutionarily, or by gradient descent — a prior whose Kolmogorov complexity matches the true distribution's structure. The lunch is real. The payment is the residual. Both are measurable. Neither has been measured together.

---

## Lineage

Immediate priors:
- THE-RELAXATION-IS-THE-COMPUTATION (Physarum/AMF via CO: integrality gaps, LP/SDP duality, approximation ratios, parameterized complexity)
- ANFINSEN (protein folding via thermodynamics: free-energy landscape, Sherman-Morrison cooperative folding, Fisher-information geodesic, AlphaFold/RFdiffusion)

NFL spine: Wolpert & Macready (IEEE Transactions on Evolutionary Computation 1, 67–82, 1997); Goldblum, Schwarzschild, Cheung, Goldblum ("The No Free Lunch Theorem, Kolmogorov Complexity, and the Role of Inductive Biases in Machine Learning," arXiv 2304.05366, ICLR 2024); Dingle, Camargo, Louis ("Input-output maps are strongly biased towards simple outputs," Nature Communications 2018); Wu et al. ("Time to Rethink AI for Combinatorial Optimization," arXiv 2502.03669, 2025–2026).

CO spine: Straszak & Vishnoi (ITCS 2016; Mathematical Programming 2021); Bonifaci (Algorithmica 2019); Byrka, Grandoni, Rothvoß, Sanità (STOC 2010/JACM 2013); Karlin, Klein, Oveis Gharan (STOC 2021); Jain (JACM 2001, survivable network design).

Thermodynamic spine: Anfinsen (Science 181, 1973); Wolynes, Onuchic, Thirumalai (Science 267, 1995); Bryngelson, Onuchic, Socci, Wolynes (Proteins 21, 1995); Jumper et al. (Nature 596, 2021, AlphaFold 2); Abramson et al. (Nature 630, 2024, AlphaFold 3); Watson et al. (Nature 620, 2023, RFdiffusion); Krishna et al. (RFdiffusion3, December 2025); Kim, Woodbury, Ahern et al. (Nature, December 2025, metallohydrolases).

Fungal-mapping spine: Stewart et al. (Science 392, 2026, DOI 10.1126/science.adu4373).

---

ERI Labs · Eric Ren · Jersey City, New Jersey · github.com/ericrenone · June 2026
