//! Module for multiple sequence alignment.
pub type MSA<'a, F> = MultipseSequenceAlignment<'a, F>;
/// A struct for a multiple sequence alignment.
/// # Example
/// ```
/// use poa_hmm::MSA;
/// let mut msa = MSA::new_default();
/// msa.register(b"AACGT", "seq1");
/// msa.register(b"AACTT", "seq2");
/// msa.register(b"ATCGT", "seq3");
/// for (id, seq) in msa.build(){
///     let seq = String::from_utf8_lossy(&seq);
///     eprintln!("{}\t{}",id, seq);
/// }
/// ```
pub struct MultipseSequenceAlignment<'a, F>
where
    F: Fn(u8, u8) -> i32,
{
    seqs: Vec<&'a [u8]>,
    ids: Vec<&'a str>,
    parameters: (i32, i32, F),
}

impl<'a, F> std::fmt::Debug for MSA<'a, F>
where
    F: Fn(u8, u8) -> i32,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        for (seq, id) in self.seqs.iter().zip(self.ids.iter()) {
            writeln!(f, "{}\t{}", id, String::from_utf8_lossy(seq))?;
        }
        let (ins, del, ref score) = self.parameters;
        let (mat, mism) = (score(0, 0), score(1, 0));
        write!(
            f,
            "Param:Ins:{}\tDel:{}\tMat:{}\tMism:{}",
            ins, del, mat, mism
        )
    }
}

fn default_score(x: u8, y: u8) -> i32 {
    if x == y {
        1
    } else {
        -1
    }
}

impl<'a> MSA<'a, fn(u8, u8) -> i32> {
    pub fn new_default() -> Self {
        Self {
            seqs: vec![],
            ids: vec![],
            parameters: (-1, -1, default_score),
        }
    }
}

impl<'a, F> MSA<'a, F>
where
    F: Fn(u8, u8) -> i32,
{
    pub fn new(ins: i32, del: i32, score: F) -> Self {
        Self {
            seqs: vec![],
            ids: vec![],
            parameters: (ins, del, score),
        }
    }
    pub fn register(&mut self, seq: &'a [u8], id: &'a str) {
        self.seqs.push(seq);
        self.ids.push(id);
    }
    pub fn build(self) -> Vec<(&'a str, Vec<u8>)> {
        let param = (self.parameters.0, self.parameters.1, &self.parameters.2);
        let msa = self
            .seqs
            .iter()
            .fold(super::POA::default(), |x, y| x.add(y, 1., param));
        self.ids
            .into_iter()
            .zip(self.seqs.iter())
            .map(|(id, seq)| {
                let q = Self::recover(&msa, seq, param);
                (id, q)
            })
            .collect()
    }
    // Recover the alignment path. It is currently O(|E|n),
    // where E is the edges of self, and n is the length of seq.
    // However, it's very easy to reduce this to approx O(|E|+n) just by
    // using depth first search. This is TODO.
    fn recover(msa: &super::POA, seq: &[u8], param: (i32, i32, &F)) -> Vec<u8> {
        let (_, tb) = msa.align(seq, param);
        let mut res = vec![b'-'; msa.num_nodes()];
        let mut q_pos = 0;
        use super::EditOp;
        for op in tb {
            if let EditOp::Match(g_pos) = op {
                res[g_pos] = seq[q_pos];
                q_pos += 1;
            }
        }
        res
    }
}
