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

    // effect |solElem> -> -|solElem> via a
    // multi-controlled phase flip gate
    fn apply_oracle(&mut self, key: i32) {
        // apply X to transform |111> into |solElem>

        for q in 0..self.num_qubits {
            if ((key >> q) & 1) == 0 {
                self.qureg.pauli_x(q);
            }
        }

        // effect |111> -> -|111>
        self.qureg
            .multi_controlled_phase_flip((0..self.num_qubits).collect());

        // apply X to transform |solElem> into |111>
        for q in 0..self.num_qubits {
            if ((key >> q) & 1) == 0 {
                self.qureg.pauli_x(q);
            }
        }
    }

    fn apply_diffuser(&mut self) {
        // apply H to transform |+> into |0>
        // apply X to transform |11..1> into |00..0>
        (0..self.num_qubits).fold(&mut self.qureg, |q, i| q.hadamard(i).pauli_x(i));

        // effect |11..1> -> -|11..1>
        // effect |111> -> -|111>
        self.qureg
            .multi_controlled_phase_flip((0..self.num_qubits).collect());

        // apply X to transform |00..0> into |11..1>
        // apply H to transform |0> into |+>
        (0..self.num_qubits).fold(&mut self.qureg, |q, i| q.pauli_x(i).hadamard(i));
    }
}

trait GroverSearch {
    fn key_bound(&self) -> i32;
    fn find_key<'a>(&'a mut self, key: i32) -> Result<&'a QuReg<'a>, Error>;
}

impl<'a> GroverSearch for GroverReg<'a> {
    fn key_bound(&self) -> i32 {
        1 << self.num_qubits
    }

    fn find_key(&mut self, key: i32) -> Result<&QuReg<'a>, Error> {
        if key < 0 || key >= self.key_bound() {
            Err(Error::KeyOutOfBound)
        } else {
            self.qureg.init_plus_state();

            let num_elems = 1 << self.num_qubits;
            let num_reps = f64::ceil(PI / 4.0 * (num_elems as QReal).sqrt()) as usize;

            // apply Grover's algorithm
            for _ in 0..num_reps {
                self.apply_oracle(key);
                self.apply_diffuser();
            }

            Ok(&self.qureg)
        }
    }
}

fn main() {
    let env = QuestEnv::new();

    // choose the system size
    const NUM_QUBITS: i32 = 15;
    const NUM_ELEMS: i32 = 1 << NUM_QUBITS;
    let num_reps = f64::ceil(PI / 4.0 * (NUM_ELEMS as QReal).sqrt()) as usize;

    println!("num_qubits: {NUM_QUBITS}, num_elems: {NUM_ELEMS}, num_reps: {num_reps}");

    // "randomly" choose the element for which to search
    let sol_elem = std::time::SystemTime::now()
        .duration_since(std::time::SystemTime::UNIX_EPOCH)
        .unwrap()
        .subsec_nanos() as i32
        % NUM_ELEMS;

    let mut grov_reg = GroverReg::new(&env, 15);
    if let Ok(qubit) = grov_reg.find_key(sol_elem) {
        println!("prob. = {}", qubit.probability_amplitude(sol_elem as i64));
    }
}
