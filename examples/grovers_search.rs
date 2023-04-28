//! Implements Grover's algorithm for unstructured search,
//! using only X, H and multi-controlled Z gates
//!
//! Adapted from grovers_search.c in QuEST repo.

use std::f64::consts::PI;

use quest_rs::{QReal, QuReg, QuestEnv};

#[derive(Debug)]
enum Error {
    KeyOutOfBound,
}

fn oracle(reg: &mut GroverReg, key: i64) -> QReal {
    let num_qubits = reg.num_qubits;
    for q in 0..num_qubits {
        if ((key >> q) & 1) == 0 {
            reg.qureg.pauli_x(q);
        }
    }

    reg.qureg
        .multi_controlled_phase_flip((0..num_qubits).collect());

    for q in 0..num_qubits {
        if ((key >> q) & 1) == 0 {
            reg.qureg.pauli_x(q);
        }
    }

    reg.qureg.probability_amplitude(key)
}

struct GroverReg<'a> {
    qureg: QuReg<'a>,
    num_qubits: i32,
}

impl<'a> GroverReg<'a> {
    fn new(env: &'a QuestEnv, num_qubits: i32) -> Self {
        GroverReg {
            qureg: QuReg::new(num_qubits, env),
            num_qubits,
        }
    }

    fn diffuse(&mut self) {
        // apply H to transform |+> into |0>
        // apply X to transform |11..1> into |00..0>
        (0..self.num_qubits).fold(&mut self.qureg, |q, i| q.hadamard(i).pauli_x(i));

        // effect |11..1> -> -|11..1>
        self.qureg
            .multi_controlled_phase_flip((0..self.num_qubits).collect());

        // apply X to transform |00..0> into |11..1>
        // apply H to transform |0> into |+>
        (0..self.num_qubits).fold(&mut self.qureg, |q, i| q.pauli_x(i).hadamard(i));
    }
}

trait GroverSearch {
    fn key_bound(&self) -> i64;
    fn find_key<'a>(&'a mut self, key: i64) -> Result<&'a GroverReg<'a>, Error>;
}

impl<'a> GroverSearch for GroverReg<'a> {
    fn key_bound(&self) -> i64 {
        1 << self.num_qubits
    }

    fn find_key(&mut self, key: i64) -> Result<&GroverReg<'a>, Error> {
        if key < 0 || key >= self.key_bound() {
            Err(Error::KeyOutOfBound)
        } else {
            self.qureg.init_plus_state();

            let num_elems = 1 << self.num_qubits;
            let num_reps = f64::ceil(PI / 4.0 * (num_elems as QReal).sqrt()) as usize;

            // apply Grover's algorithm
            for _ in 0..num_reps {
                oracle(self, key);
                self.diffuse();
            }

            Ok(self)
        }
    }
}

fn main() {
    let env = QuestEnv::new();

    // choose the system size
    const NUM_QUBITS: i64 = 15;
    const NUM_ELEMS: i64 = 1 << NUM_QUBITS;
    let num_reps = f64::ceil(PI / 4.0 * (NUM_ELEMS as QReal).sqrt()) as usize;

    println!("num_qubits: {NUM_QUBITS}, num_elems: {NUM_ELEMS}, num_reps: {num_reps}");

    // "randomly" choose the element for which to search
    let sol_elem = std::time::SystemTime::now()
        .duration_since(std::time::SystemTime::UNIX_EPOCH)
        .unwrap()
        .subsec_nanos() as i64
        % NUM_ELEMS;

    let mut grov_reg = GroverReg::new(&env, 15);
    if let Ok(reg) = grov_reg.find_key(sol_elem) {
        println!(
            "prob. = {}",
            reg.qureg.probability_amplitude(sol_elem as i64)
        );
    }
}
