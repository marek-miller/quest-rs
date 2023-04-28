// /** @file
//  * Implements Grover's algorithm for unstructured search,
//  * using only X, H and multi-controlled Z gates
//  *
//  * Compile and run from within the build folder, using:
// cmake .. -DUSER_SOURCE=../examples/grovers_search.c \
//         -DOUTPUT_EXE=grovers
// make
// ./grovers
//  *
//  *
//  * @author Tyson Jones
//  */
//  #include <stdio.h>
//  #include <math.h>
//  #include <time.h>
//  #include <stdlib.h>

//  #include "QuEST.h"

//  /* effect |solElem> -> -|solElem> via a
//   * multi-controlled phase flip gate
//   */
//  void applyOracle(Qureg qureg, int numQubits, int solElem) {

//      // apply X to transform |111> into |solElem>
//      for (int q=0; q<numQubits; q++)
//          if (((solElem >> q) & 1) == 0)
//              pauliX(qureg, q);

//      // effect |111> -> -|111>
//      int ctrls[numQubits];
//      for (int q=0; q<numQubits; q++)
//          ctrls[q] = q;
//      multiControlledPhaseFlip(qureg, ctrls, numQubits);

//      // apply X to transform |solElem> into |111>
//      for (int q=0; q<numQubits; q++)
//          if (((solElem >> q) & 1) == 0)
//              pauliX(qureg, q);
//  }

//  /* apply 2|+><+|-I by transforming into the Hadamard basis
//   * and effecting 2|0><0|-I. We do this, by observing that
//   *   c..cZ = diag{1,..,1,-1}
//   *         = I - 2|1..1><1..1|
//   * and hence
//   *   X..X c..cZ X..X = I - 2|0..0><0..0|
//   * which differs from the desired 2|0><0|-I state only by
//   * the irrelevant global phase pi
//   */
//  void applyDiffuser(Qureg qureg, int numQubits) {

//      // apply H to transform |+> into |0>
//      for (int q=0; q<numQubits; q++)
//          hadamard(qureg, q);

//      // apply X to transform |11..1> into |00..0>
//      for (int q=0; q<numQubits; q++)
//          pauliX(qureg, q);

//      // effect |11..1> -> -|11..1>
//      int ctrls[numQubits];
//      for (int q=0; q<numQubits; q++)
//          ctrls[q] = q;
//      multiControlledPhaseFlip(qureg, ctrls, numQubits);

//      // apply X to transform |00..0> into |11..1>
//      for (int q=0; q<numQubits; q++)
//          pauliX(qureg, q);

//      // apply H to transform |0> into |+>
//      for (int q=0; q<numQubits; q++)
//          hadamard(qureg, q);
//  }

//  int main() {

//      // prepare the hardware-agnostic QuEST environment
//      QuESTEnv env = createQuESTEnv();

//      // choose the system size
//      int numQubits = 15;
//      int numElems = (int) pow(2, numQubits);
//      int numReps = ceil(M_PI/4 * sqrt(numElems));

//      printf("numQubits: %d, numElems: %d, numReps: %d\n",
//          numQubits, numElems, numReps);

//      // randomly choose the element for which to search
//      srand(time(NULL));
//      int solElem = rand() % numElems;

//      // prepare |+>
//      Qureg qureg = createQureg(numQubits, env);
//      initPlusState(qureg);

//      // apply Grover's algorithm
//      for (int r=0; r<numReps; r++) {
//          applyOracle(qureg, numQubits, solElem);
//          applyDiffuser(qureg, numQubits);

//          // monitor the probability of the solution state
//          printf("prob of solution |%d> = " REAL_STRING_FORMAT "\n",
//              solElem, getProbAmp(qureg, solElem));
//      }

//      // free memory
//      destroyQureg(qureg, env);
//      destroyQuESTEnv(env);
//      return 0;
//  }

use core::num;
use std::f64::consts::PI;

use quest_rs::{Complex, ComplexMatrix2, ComplexMatrixN, QReal, QuReg, QuestEnv, Vector};

// effect |solElem> -> -|solElem> via a
// multi-controlled phase flip gate
fn applyOracle(qureg: &mut QuReg, num_qubits: i32, sol_elem: i32) {
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

/* apply 2|+><+|-I by transforming into the Hadamard basis
 * and effecting 2|0><0|-I. We do this, by observing that
 *   c..cZ = diag{1,..,1,-1}
 *         = I - 2|1..1><1..1|
 * and hence
 *   X..X c..cZ X..X = I - 2|0..0><0..0|
 * which differs from the desired 2|0><0|-I state only by
 * the irrelevant global phase pi
 */
fn applyDiffuser(qureg: &mut QuReg, num_qubits: i32) {
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
    // let env = QuestEnv::new();
    // let mut qubits = QuReg::new(3, &env);

    // qubits.init_zero_state();

    // println!("Out environment is:");
    // qubits.report_params();
    // env.report();

    // // Set up the circuitry

    // let unitary_alpha = Complex::new(0.5, 0.5);
    // let unitary_beta = Complex::new(0.5, -0.5);

    // let unitary_matrix = ComplexMatrix2 {
    //     real: [[0.5, 0.5], [0.5, 0.5]],
    //     imag: [[0.5, -0.5], [-0.5, 0.5]],
    // };

    // let mut toffoli_gate = ComplexMatrixN::new(3);
    // for i in 0..6 {
    //     toffoli_gate.set_real(i, i, 1.0);
    // }
    // toffoli_gate.set_real(6, 7, 1.0);
    // toffoli_gate.set_real(7, 6, 1.0);

    // qubits
    //     .hadamard(0)
    //     .controlled_not(0, 1)
    //     .rotate_y(2, 0.1)
    //     .multi_controlled_phase_flip(vec![0, 1, 2])
    //     .unitary(0, unitary_matrix)
    //     .compact_unitary(1, unitary_alpha, unitary_beta)
    //     //.rotate_around_axis(2, (PI / 2) as QReal, Vector::new(1.0, 0.0, 0.0))
    //     .rotate_around_axis(2, (3.14 / 2.0) as QReal, Vector::new(1.0, 0.0, 0.0))
    //     .controlled_compact_unitary(0, 1, unitary_alpha, unitary_beta)
    //     .multi_controlled_unitary(vec![0, 1], 2, unitary_matrix)
    //     .multi_qubit_unitary(vec![0, 1, 2], toffoli_gate);

    // // Study the output
    // println!("Circuit output:");
    // println!("---------------");

    // let prob_amp_state_111 = qubits.probability_amplitude(0b111);
    // println!("Probability amplitude of |111> is: {}", prob_amp_state_111);

    // let prob_qubit_two_in_state_1 = qubits.calculate_probability_of_outcome(2, 1);
    // println!(
    //     "Probability of qubit 2 being in state 1: {}",
    //     prob_qubit_two_in_state_1
    // );

    // println!("Qubit 0 was measured in state: {}", qubits.measure(0));
    // let (outcome, outcome_probability) = qubits.measure_with_stats(2);
    // println!(
    //     "Qubit 2 collapsed to {} with probability {}",
    //     outcome, outcome_probability
    // );


    
     // prepare the hardware-agnostic QuEST environment
    let env = QuestEnv::new();
    
    const num_qubits: i32 = 15;
    let num_elems = 1<<num_qubits;
    let num_reps = f64::ceil(PI / 4.0 * (num_elems as QReal).sqrt()) as usize;
    

    println!("num_qubits: {num_qubits}, num_elems: {num_elems}, num_reps: {num_reps}");
    
    
    //  // choose the system size
    //  int numQubits = 15;
    //  int numElems = (int) pow(2, numQubits);
    //  int numReps = ceil(M_PI/4 * sqrt(numElems));
    
    
    // randomly choose the element for which to search
    let sol_elem = 3253523 % num_elems;
    //  srand(time(NULL));
    //  int solElem = rand() % numElems;
    
    // prepare |+>
    let mut qureg = QuReg::new(num_qubits, &env);
    let qureg = qureg.init_plus_state();
     
     // apply Grover's algorithm
     for r in 0..num_reps {
        applyOracle(qureg, num_qubits, sol_elem);
        applyDiffuser(qureg, num_qubits);
        println!("prob of solution |{}> = {}", sol_elem, qureg.probability_amplitude(sol_elem as i64));
     }
    //  for (int r=0; r<numReps; r++) {
    //      applyOracle(qureg, numQubits, solElem);
    //      applyDiffuser(qureg, numQubits);

    //      // monitor the probability of the solution state
    //      printf("prob of solution |%d> = " REAL_STRING_FORMAT "\n",
    //          solElem, getProbAmp(qureg, solElem));
    //  }

     // free memory
    //  destroyQureg(qureg, env);
    //  destroyQuESTEnv(env);
    //  return 0;
}
