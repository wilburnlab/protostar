"""Per-group consensus assembly.

Thin glue tying the extracted replicate spectra (``ms2_extract``) to
Constellation's fragment basis + consensus builder. For a replicate group of
one ``(modified_sequence, charge, mode)``, build the fixed fragment basis and
aggregate the per-replicate spectra into a ``ConsensusSpectrum`` (retaining the
``per_replicate[R, K]`` matrix the multinomial variance-law analysis reads).

No model code: the basis, alignment, and consensus statistics all come from
``constellation.massspec.spectra``.
"""

from __future__ import annotations

from collections.abc import Sequence

import torch
from constellation.core.sequence.proforma import parse_proforma
from constellation.massspec.peptide.ions import IonType
from constellation.massspec.spectra.consensus import (
    ConsensusSpectrum,
    FragmentBasis,
    build_consensus,
    fragment_basis,
)

#: The b/y-only v1 channel basis (charges 1–2, no losses), matching the XIC
#: level-2 extraction defaults.
DEFAULT_ION_TYPES: tuple[IonType, ...] = (IonType.B, IonType.Y)
DEFAULT_MAX_FRAGMENT_CHARGE: int = 2


def basis_for(
    modified_sequence: str,
    *,
    ion_types: Sequence[IonType] = DEFAULT_ION_TYPES,
    max_fragment_charge: int = DEFAULT_MAX_FRAGMENT_CHARGE,
) -> FragmentBasis:
    """Fixed fragment-channel basis for a peptidoform (parsed from its ProForma
    modified sequence)."""
    return fragment_basis(
        parse_proforma(modified_sequence),
        ion_types=ion_types,
        max_fragment_charge=max_fragment_charge,
    )


def assemble_consensus(
    modified_sequence: str,
    spectra: Sequence[tuple[torch.Tensor, torch.Tensor]],
    *,
    aggregate: str = "median",
    ion_types: Sequence[IonType] = DEFAULT_ION_TYPES,
    max_fragment_charge: int = DEFAULT_MAX_FRAGMENT_CHARGE,
) -> ConsensusSpectrum:
    """Aggregate replicate ``(mz_theoretical, intensity)`` spectra of one
    peptide onto its fragment basis.

    The spectra carry **theoretical** m/z (from the XIC trace), so they align to
    the basis exactly. ``aggregate="median"`` is the robust default for the
    variance-law analysis; ``"sum"`` gives the pooled count-space estimate used
    by the MSP-provenance forensic."""
    basis = basis_for(
        modified_sequence, ion_types=ion_types, max_fragment_charge=max_fragment_charge
    )
    return build_consensus(spectra, basis, aggregate=aggregate)


__all__ = [
    "DEFAULT_ION_TYPES",
    "DEFAULT_MAX_FRAGMENT_CHARGE",
    "basis_for",
    "assemble_consensus",
]
