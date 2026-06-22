//! Vendor-grade Bruker timsTOF mobility recalibration (scan -> 1/K0) from the `TimsCalibration`
//! `ModelType = 2` model stored in `analysis.tdf`.
//!
//! `timsrust`'s [`Scan2ImConverter`](timsrust::converters::Scan2ImConverter) approximates 1/K0 by
//! linearly interpolating the *nominal* acquisition range, which is ~0.03 Vs.s/cm^2 too low at the
//! high-mobility edge relative to Bruker's `timsdata` SDK. This evaluates the actual vendor model:
//!
//! ```text
//!   V(scan) = C2 + (C3 - C2) * (scan - C0) / (C1 - C0)     // linear TIMS voltage ramp
//!   1/K0    = (V + delta) / (C7 + C6 * V)                  // vendor rational model
//! ```
//!
//! `delta` is anchored so the lowest-mobility scan (V = C3) reproduces
//! `GlobalMetadata.OneOverK0AcqRangeLower`. The model was reverse-engineered from Bruker's SDK and
//! validated against `tims_scannum_to_oneoverk0` on 68 PRIDE timsTOF datasets (nscans 473..2831):
//! per-dataset max error median 1.4e-3 / worst 2.8e-3 Vs.s/cm^2, vs ~3.0e-2 for the linear
//! approximation. Across 161 distinct PRIDE timsTOF datasets every `TimsCalibration` row was
//! `ModelType = 2`.
//!
//! Coefficient laws confirmed across all 68 datasets (the Taylor series of the rational):
//!   linear term = 1/C7, quadratic = -C6/C7^2, cubic = +C6^2/C7^3.
//!
//! Only `ModelType = 2` is handled (the C-columns mean different things for other model types);
//! [`from_connection`](TimsMobilityCalibration::from_connection) returns `None` otherwise, and the
//! reader keeps timsrust's linear conversion.

use rusqlite::{Connection, OptionalExtension};

/// Bruker timsTOF `ModelType = 2` mobility calibration: mobility scan index -> 1/K0 (Vs.s/cm^2).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TimsMobilityCalibration {
    c0: f64,
    c1: f64,
    c2: f64,
    c3: f64,
    c6: f64,
    c7: f64,
    /// Voltage offset, anchored so V = C3 maps to `OneOverK0AcqRangeLower`.
    delta: f64,
}

impl TimsMobilityCalibration {
    /// Build from the raw `TimsCalibration` coefficients + `GlobalMetadata.OneOverK0AcqRangeLower`.
    pub fn new(c0: f64, c1: f64, c2: f64, c3: f64, c6: f64, c7: f64, one_over_k0_lower: f64) -> Self {
        // Anchor: at the lowest-mobility scan V = C3,
        //   one_over_k0_lower = (C3 + delta) / (C7 + C6*C3)  =>  delta = lower*(C7 + C6*C3) - C3.
        let delta = one_over_k0_lower * (c7 + c6 * c3) - c3;
        Self { c0, c1, c2, c3, c6, c7, delta }
    }

    /// Inverse reduced ion mobility 1/K0 (Vs.s/cm^2) for a 0-based mobility scan index. Fractional
    /// indices interpolate, matching the SDK.
    #[inline]
    pub fn one_over_k0(&self, scan: f64) -> f64 {
        let span = self.c1 - self.c0;
        let v = if span == 0.0 {
            self.c2
        } else {
            self.c2 + (self.c3 - self.c2) * (scan - self.c0) / span
        };
        (v + self.delta) / (self.c7 + self.c6 * v)
    }

    /// Load from an open `analysis.tdf` connection. Returns `None` when there is no `ModelType = 2`
    /// row (or the metadata can't be read) so the caller falls back to the linear conversion.
    pub fn from_connection(conn: &Connection) -> Option<Self> {
        let (c0, c1, c2, c3, c6, c7): (f64, f64, f64, f64, f64, f64) = conn
            .query_row(
                "SELECT C0, C1, C2, C3, C6, C7 FROM TimsCalibration WHERE ModelType = 2 ORDER BY Id LIMIT 1",
                [],
                |r| Ok((r.get(0)?, r.get(1)?, r.get(2)?, r.get(3)?, r.get(4)?, r.get(5)?)),
            )
            .optional()
            .ok()
            .flatten()?;
        let lower: f64 = conn
            .query_row(
                "SELECT CAST(Value AS REAL) FROM GlobalMetadata WHERE Key = 'OneOverK0AcqRangeLower'",
                [],
                |r| r.get(0),
            )
            .ok()?;
        Some(Self::new(c0, c1, c2, c3, c6, c7, lower))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reproduces_vendor_sdk_at_endpoints() {
        // SBA415 TimsCalibration row + acq-range lower; SDK 1/K0 ground truth at the endpoints.
        let cal = TimsMobilityCalibration::new(
            1.0, 909.0, 211.45198604901222, 73.95258004355563, 0.00492817555366883,
            131.11541877221117, 0.600,
        );
        assert!((cal.one_over_k0(909.0) - 0.600).abs() < 1e-9); // anchored at last scan
        assert!((cal.one_over_k0(0.0) - 1.6385).abs() < 1e-3); // SDK 1.6383 (linear gives 1.600)
        assert!(cal.one_over_k0(100.0) > cal.one_over_k0(800.0)); // monotonic
    }
}
