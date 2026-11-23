"""Sweep reducers for the extraction pilots (short-lived prototypes; ledger promotion entry:
the screen reducer). They emit measurements only and do no scoring.

Kept in their own module so spawned sweep workers can import them by name."""

from __future__ import annotations

import numpy as np
import pyarrow as pa

ISOTOPES = (-1, 0, 1, 2, 3)
_ISO_NAMES = ("m1", "0", "1", "2", "3")


class IsotopeScanRows:
    """Per (query, MS1 scan) with ≥ 2 of M0..M+2 matched: raw isotope intensities and ppm errors.

    Pilot prototype; measurements only."""

    name = "isotope_scans"
    gated = False
    zero_fills = False
    schema = pa.schema(
        [
            pa.field("target_id", pa.int64()),
            pa.field("scan", pa.int32()),
            pa.field("rt", pa.float64()),
            pa.field("iit_ms", pa.float32()),
            pa.field("tic", pa.float32()),
            pa.field("floor", pa.float32()),
            pa.field("n_iso", pa.int8()),
            *[pa.field(f"i_{n}", pa.float32()) for n in _ISO_NAMES],
            *[pa.field(f"e_{n}", pa.float32()) for n in _ISO_NAMES],
        ]
    )
    outputs = {
        "isotope_scans": schema,
        "query_channels": pa.schema(
            [
                pa.field("target_id", pa.int64()),
                pa.field("isotope", pa.int8()),
                pa.field("mz", pa.float64()),
                pa.field("abundance", pa.float32()),
            ]
        ),
    }

    def reduce(self, batch):
        plan = batch.plan
        ch = plan.channels_for_level(1)
        tid_col = plan.targets.column("target_id")
        # the batch's queries' channels (theoretical m/z and shares), once per query
        q = batch.queries
        counts = ch.counts()[q]
        flat_q = np.repeat(ch.offsets[q], counts) + (
            np.arange(counts.sum()) - np.repeat(np.cumsum(counts) - counts, counts)
        )
        qch = pa.table(
            {
                "target_id": tid_col.take(pa.array(np.repeat(q, counts))),
                "isotope": pa.array(ch.isotope[flat_q].astype(np.int8)),
                "mz": pa.array(ch.mz[flat_q]),
                "abundance": pa.array(ch.abundance[flat_q].astype(np.float32)),
            },
            schema=self.outputs["query_channels"],
        )
        c = batch.cells(1, fill=False)
        if c.n == 0:
            return {"isotope_scans": self.schema.empty_table(), "query_channels": qch}
        flat = ch.offsets[c.query] + c.channel
        iso = ch.isotope[flat].astype(np.int64)
        mzt = ch.mz[flat]
        ppm = (c.mz_observed - mzt) / mzt * 1e6
        n_rows = int(batch.scans.scan.size)
        key = c.query.astype(np.int64) * n_rows + c.row.astype(np.int64)
        uk, inv = np.unique(key, return_inverse=True)
        n = uk.size
        cols_i = {k: np.zeros(n, np.float32) for k in ISOTOPES}
        cols_e = {k: np.full(n, np.nan, np.float32) for k in ISOTOPES}
        for k in ISOTOPES:
            sel = iso == k
            cols_i[k][inv[sel]] = c.intensity[sel]
            cols_e[k][inv[sel]] = ppm[sel]
        n_iso = sum((cols_i[k] > 0).astype(np.int8) for k in (0, 1, 2))
        keep = n_iso >= 2
        qq = (uk // n_rows)[keep]
        rr = (uk % n_rows)[keep]
        sc = batch.scans

        def arr(x, dtype):
            return pa.array(
                np.zeros(rr.size, dtype) if x is None else np.asarray(x)[rr].astype(dtype)
            )

        table = pa.table(
            {
                "target_id": tid_col.take(pa.array(qq)),
                "scan": arr(sc.scan, np.int32),
                "rt": arr(sc.rt, np.float64),
                "iit_ms": arr(sc.iit, np.float32),
                "tic": arr(sc.tic, np.float32),
                "floor": arr(sc.floor, np.float32),
                "n_iso": pa.array(n_iso[keep]),
                **{
                    f"i_{nm}": pa.array(cols_i[k][keep])
                    for k, nm in zip(ISOTOPES, _ISO_NAMES, strict=True)
                },
                **{
                    f"e_{nm}": pa.array(cols_e[k][keep])
                    for k, nm in zip(ISOTOPES, _ISO_NAMES, strict=True)
                },
            },
            schema=self.schema,
        )
        return {"isotope_scans": table, "query_channels": qch}
