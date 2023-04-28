//! Implements Grover's algorithm for unstructured search,
//! using only X, H and multi-controlled Z gates
//!
//! Adapted from grovers_search.c in QuEST repo.

use std::f64::consts::PI;

use quest_rs::{QReal, QuReg, QuestEnv};

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

    // prepare |+>
    let mut qureg = QuReg::new(NUM_QUBITS, &env);
    let qureg = qureg.init_plus_state();

    // effect |solElem> -> -|solElem> via a
    // multi-controlled phase flip gate
    let apply_oracle = |qureg: &mut QuReg| {
        // apply X to transform |111> into |solElem>

        for q in 0..NUM_QUBITS {
            if ((sol_elem >> q) & 1) == 0 {
                qureg.pauli_x(q);
            }
        }

        // effect |111> -> -|111>
        qureg.multi_controlled_phase_flip((0..NUM_QUBITS).collect());

        // apply X to transform |solElem> into |111>
        for q in 0..NUM_QUBITS {
            if ((sol_elem >> q) & 1) == 0 {
                qureg.pauli_x(q);
            }
        }
    };

    let apply_diffuser = |qureg: &mut QuReg| {
        // apply H to transform |+> into |0>
        // apply X to transform |11..1> into |00..0>
        let qureg = (0..NUM_QUBITS).fold(qureg, |q, i| q.hadamard(i).pauli_x(i));

        // effect |11..1> -> -|11..1>
        // effect |111> -> -|111>
        qureg.multi_controlled_phase_flip((0..NUM_QUBITS).collect());

        // apply X to transform |00..0> into |11..1>
        // apply H to transform |0> into |+>
        (0..NUM_QUBITS).fold(qureg, |q, i| q.pauli_x(i).hadamard(i));
    };

    // apply Grover's algorithm
    for _ in 0..num_reps {
        apply_oracle(qureg);
        apply_diffuser(qureg);
        println!(
            "prob of solution |{}> = {}",
            sol_elem,
            qureg.probability_amplitude(sol_elem as i64)
        );
    }
}
