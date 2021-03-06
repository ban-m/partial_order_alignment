//! Here we define some useful configurations for Partial Order Hidden Markov Model.
//! Though it hopefully would be a good guess, I recommend to calculate these parameters for each dataset one by one.

/// Configuration for BADREAD read simulator.
/// There parameters are tuned by hand and Last program.
pub const BADREAD_CONFIG: Config = Config {
    mismatch: 0.0344,
    base_freq: [0.25, 0.25, 0.25, 0.25],
    p_match: 0.88,
    p_ins: 0.0549,
    p_del: 0.0651,
    p_extend_ins: 0.0337,
    p_extend_del: 0.1787,
    p_del_to_ins: 0.0,
};

/// Default configuration. Parameters are tuned just by my guess.
pub const DEFAULT_CONFIG: Config = Config {
    mismatch: 0.04,
    base_freq: [0.25, 0.25, 0.25, 0.25],
    p_match: 0.85,
    p_ins: 0.07,
    p_del: 0.04,
    p_extend_ins: 0.08,
    p_extend_del: 0.07,
    p_del_to_ins: 0.08,
};

/// Strict configuration. In other words, this configuration is for
/// CCS reads.
pub const STRICT_CONFIG: Config = Config {
    mismatch: 0.02,
    base_freq: [0.25, 0.25, 0.25, 0.25],
    p_match: 0.95,
    p_ins: 0.03,
    p_del: 0.02,
    p_extend_ins: 0.03,
    p_extend_del: 0.02,
    p_del_to_ins: 0.02,
};

/// Configuration for PacBio's sequel system.
/// There parameters are tuned by Last program.
pub const PACBIO_CONFIG: Config = Config {
    mismatch: 0.0341,
    base_freq: [0.28, 0.22, 0.22, 0.28],
    p_match: 0.9124,
    p_ins: 0.0606,
    p_del: 0.0270,
    p_extend_ins: 0.3014,
    p_extend_del: 0.0695,
    p_del_to_ins: 0.0086,
};

/// A configure struct.
#[derive(Debug, Clone, Default)]
pub struct Config {
    /// Mismatch probability at given position. # mismatch/(#mism + #match)
    pub mismatch: f64,
    /// The base composition of the read. A,C,G, and T.
    pub base_freq: [f64; 4],
    /// probability of matching at match states.
    pub p_match: f64,
    /// probability of starting insertion at match states.
    pub p_ins: f64,
    /// probability of starting deletion at match states.
    pub p_del: f64,
    /// probability of extending insertion at insertion state.
    pub p_extend_ins: f64,
    /// same for deletion
    pub p_extend_del: f64,
    /// probability of jump deletion state to insertion state.
    pub p_del_to_ins: f64,
}

impl std::fmt::Display for Config {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(f, "Mismatch:{:.3}", self.mismatch)?;
        writeln!(f, "P_match:{:.3}", self.p_match)?;
        writeln!(f, "P_insertion:{:.3}", self.p_ins)?;
        writeln!(f, "P_deletion:{:.3}", self.p_del)?;
        writeln!(f, "P_extending_insertion:{:.3}", self.p_extend_ins)?;
        writeln!(f, "P_extending_deletion:{:.3}", self.p_extend_del)?;
        writeln!(f, "P_deletion_to_insertion:{:.3}", self.p_extend_del)?;
        writeln!(
            f,
            "BaseFrequency:[A,C,G,T]=[{:.3},{:.3},{:.3},{:.3}]",
            self.base_freq[0], self.base_freq[1], self.base_freq[2], self.base_freq[3]
        )
    }
}

impl Config {
    /// Create a new configuration from the given parameters.
    pub fn new(
        mismatch: f64,
        (p_match, p_ins, p_del): (f64, f64, f64),
        (p_extend_ins, p_extend_del, p_del_to_ins): (f64, f64, f64),
        base_freq: [f64; 4],
    ) -> Self {
        Self {
            mismatch,
            base_freq,
            p_match,
            p_ins,
            p_del,
            p_extend_ins,
            p_extend_del,
            p_del_to_ins,
        }
    }
    /// Return the likelihood of the null model.
    /// Here, the null model is the single-state 0-th order Markov model
    /// ,which generate each base (A, T, G, and C) by base_freq[base].
    pub fn likelihood_of_null_model(&self, seq: &[u8]) -> f64 {
        let mut counts = [0; 4];
        for &base in seq {
            counts[super::base_table::BASE_TABLE[base as usize]] += 1;
        }
        counts
            .iter()
            .zip(self.base_freq.iter())
            .map(|(&c, f)| (c as f64) * f.ln())
            .sum::<f64>()
    }
}
