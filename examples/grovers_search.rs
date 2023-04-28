//! Implements Grover's algorithm for unstructured search,
//! using only X, H and multi-controlled Z gates
//!
//! Adapted from grovers_search.c in QuEST repo.

use std::f64::consts::PI;

use quest_rs::{QReal, QuReg, QuestEnv};

// effect |solElem> -> -|solElem> via a
// multi-controlled phase flip gate
fn apply_oracle(qureg: &mut QuReg, num_qubits: i32, sol_elem: i32) {
    // apply X to transform |111> into |solElem>
    //  for (int q=0; q<num_qubits; q++)
    for q in 0..num_qubits {
        if ((sol_elem >> q) & 1) == 0 {
            qureg.pauli_x(q);
        }
    }

    // effect |111> -> -|111>
    let control_qubits = (0..num_qubits).collect::<Vec<_>>();
    qureg.multi_controlled_phase_flip(control_qubits);

    // apply X to transform |solElem> into |111>
    for q in 0..num_qubits {
        if ((sol_elem >> q) & 1) == 0 {
            qureg.pauli_x(q);
        }
    }
}

// apply 2|+><+|-I by transforming into the Hadamard basis
// and effecting 2|0><0|-I. We do this, by observing that
//   c..cZ = diag{1,..,1,-1}
//         = I - 2|1..1><1..1|
// and hence
//   X..X c..cZ X..X = I - 2|0..0><0..0|
// which differs from the desired 2|0><0|-I state only by
// the irrelevant global phase pi
fn apply_diffuser(qureg: &mut QuReg, num_qubits: i32) {
    // apply H to transform |+> into |0>
    for q in 0..num_qubits {
        qureg.hadamard(q);
    }

    // apply X to transform |11..1> into |00..0>
    for q in 0..num_qubits {
        qureg.pauli_x(q);
    }

    // effect |11..1> -> -|11..1>
    // effect |111> -> -|111>
    let control_qubits = (0..num_qubits).collect::<Vec<_>>();
    qureg.multi_controlled_phase_flip(control_qubits);

    // apply X to transform |00..0> into |11..1>
    for q in 0..num_qubits {
        qureg.pauli_x(q);
    }

    // apply H to transform |0> into |+>
    for q in 0..num_qubits {
        qureg.hadamard(q);
    }
}

fn main() {
    let env = QuestEnv::new();

    // choose the system size
    const NUM_QUBITS: i32 = 15;
    let num_elems = 1 << NUM_QUBITS;
    let num_reps = f64::ceil(PI / 4.0 * (num_elems as QReal).sqrt()) as usize;

    println!("num_qubits: {NUM_QUBITS}, num_elems: {num_elems}, num_reps: {num_reps}");

    // "randomly" choose the element for which to search
    let sol_elem = std::time::SystemTime::now()
        .duration_since(std::time::SystemTime::UNIX_EPOCH)
        .unwrap()
        .subsec_nanos() as i32
        % num_elems;

    // prepare |+>
    let mut qureg = QuReg::new(NUM_QUBITS, &env);
    let qureg = qureg.init_plus_state();

    // apply Grover's algorithm
    for _ in 0..num_reps {
        apply_oracle(qureg, NUM_QUBITS, sol_elem);
        apply_diffuser(qureg, NUM_QUBITS);
        println!(
            "prob of solution |{}> = {}",
            sol_elem,
            qureg.probability_amplitude(sol_elem as i64)
        );
    }
}
