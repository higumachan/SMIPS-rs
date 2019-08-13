use std::cmp::Ordering;
use num::{Num, range};
use std::fmt;
use std::fmt::{Debug, Formatter};

struct MIP<T: num::Float> {
    variable_count: usize,
    equation_count: usize,
    objective_coef: Vec<T>,
    equation_left: Vec<Vec<T>>,
    equation_compare: Vec<Ordering>,
    equation_right: Vec<T>,
    integer_flag: Vec<bool>,
}


struct SimplexModule<T: num::Float> {
    table: Vec<Vec<T>>,
    right_coef: Vec<T>,
    z_index: usize,
    base_variables: Vec<Option<usize>>,
}

enum SolveState {
    Optimal,
    Infeasible,
}

trait OrdFloat: num::Float + Ord {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap_or(Ordering::Equal)
    }
}


impl<T: num::Float> SimplexModule<T> {
    fn new(table: Vec<Vec<T>>, right_coef: Vec<T>, z_index: usize) -> Self {
        let base_variables = (0..table.len()).map(|_| None).collect();
        SimplexModule { table, right_coef, z_index, base_variables }
    }

    fn solve(&mut self) -> SolveState {
        while let Some(col_idx) = self.table[self.z_index].iter().enumerate().filter(|x| x.1 < &T::zero()).next().map(|x|x.0) {
            let row_idx = self.table.iter()
                .map(|x| if x[col_idx].is_zero() { T::infinity() } else { self.right_coef[col_idx] / x[col_idx] })
                .enumerate()
                .filter(|x| x.0 != self.z_index)
                .min_by(|x, y| x.1.partial_cmp(&y.1).unwrap_or(Ordering::Equal)).map(|x| x.0);
            if row_idx.is_none() {
                return SolveState::Infeasible;
            }
            let row_idx = row_idx.unwrap();

            let target_value = self.table[row_idx][col_idx];
            self.table[row_idx] = self.table[row_idx].iter().map(|x| *x / target_value).collect();
            self.right_coef[row_idx] = self.right_coef[row_idx] / target_value;

            for i in 0..self.table.len() {
                if i == row_idx {
                    continue;
                }

                let coef = self.table[i][col_idx];
                self.table[i] = self.table[i].iter().zip(self.table[row_idx].clone()).map(|x| *x.0 - (coef * x.1)).collect();
                self.right_coef[i] = self.right_coef[i] - coef * self.right_coef[row_idx];
            }

            self.base_variables[row_idx] = Some(col_idx);
        }

        SolveState::Optimal
    }
}

impl Debug for SimplexModule<f64> {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        writeln!(f, "objective={} ", self.right_coef[self.z_index]);
        for (base, v) in self.base_variables.iter().zip(self.right_coef.clone()) {
            if let Some(base) = base {
                writeln!(f, "x_{}={}", base, v);
            }
        }
        Ok(())
    }
}


fn main() {

    let mut sm = SimplexModule::new(
        vec![
            vec![1.5, 3.0, 1.0, 0.0],
            vec![3.0, 1.0, 0.0, 1.0],
            vec![-5.0, -4.0, 0.0, 0.0],
        ],
        vec![13.5, 10.0, 0.0],
        2,
    );

    sm.solve();

    println!("{:?}", sm);
}
