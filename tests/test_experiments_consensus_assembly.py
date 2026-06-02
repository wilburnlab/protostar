"""Smoke test for the consensus-assembly glue (the statistics are covered by
constellation's consensus tests; here we just verify the wiring)."""

from __future__ import annotations

import torch

from protostar.experiments import consensus_assembly


def test_basis_for_peptide():
    basis = consensus_assembly.basis_for("PEPTIDE")
    assert basis.K == 24  # 6 b + 6 y per charge, charges 1-2


def test_assemble_consensus_wires_basis_and_replicates():
    modseq = "PEPTIDE"
    basis = consensus_assembly.basis_for(modseq)
    ch = basis.mz_theoretical[:3].clone()  # first 3 channels' theoretical m/z
    spectra = [
        (ch.clone(), torch.tensor([6.0, 4.0, 2.0], dtype=torch.float64)),
        (ch.clone(), torch.tensor([5.0, 5.0, 2.0], dtype=torch.float64)),
    ]
    cons = consensus_assembly.assemble_consensus(modseq, spectra, aggregate="median")
    assert cons.n_replicates == 2
    assert cons.per_replicate.shape == (2, basis.K)
    # intensities aligned (by exact theoretical m/z) into the first 3 channels
    assert cons.per_replicate[0, 0].item() == 6.0
    assert cons.per_replicate[1, 1].item() == 5.0
    assert cons.per_replicate[:, 3:].sum().item() == 0.0
