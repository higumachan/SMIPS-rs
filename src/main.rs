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

#[derive(Debug)]
struct Solution<T: num::Float> {
    variables: Vec<T>,
    objective: T,
}

impl<T: num::Float> MIP<T> {
    pub fn solve(&self) -> Solution<T>{
        let (mut sm, slack_count, artificial_count) = self.build_simplex_module_step1();
        sm.solve();
        self.transform_simplex_module_step2(&mut sm, slack_count, artificial_count);
        sm.solve();

        self.build_solution(&sm)
    }

    fn build_solution(&self, simplex_module: &SimplexModule<T>) -> Solution<T> {
        let mut variables = vec![T::zero(); self.variable_count];

        for i in 0..self.equation_count {
            if 0 <= simplex_module.base_variables[i] && simplex_module.base_variables[i] < self.variable_count {
                variables[simplex_module.base_variables[i]] = simplex_module.right_coef[i];
            }
        }

        Solution { variables, objective: simplex_module.right_coef[simplex_module.z_index] }
    }

    fn build_simplex_module_step1(&self) -> (SimplexModule<T>, usize, usize) {
        let mut s_count = 0usize;
        let mut a_count = 0usize;


        for i in 0..self.equation_count {
            if self.equation_compare[i] == Ordering::Equal {
                a_count += 1;
            }
            else {
                s_count += 1;
                if (self.equation_right[i] < T::zero() && self.equation_compare[i] == Ordering::Less) ||
                    (self.equation_right[i] >= T::zero() && self.equation_compare[i] == Ordering::Greater) {
                    a_count += 1;
                }
            }
        }


        let mut table: Vec<Vec<T>> = vec![];
        let mut s_idx = 0;
        let mut a_idx = 0;
        let mut base_variables = vec![];


        for i in 0..self.equation_count {
            table.push(self.equation_left[i].clone());
            let row = table.last_mut().unwrap();
            row.append(&mut vec![T::zero(); s_count + a_count]);
            if self.equation_compare[i] == Ordering::Equal {
                let idx = self.variable_count + s_count + a_idx;
                row[idx] = T::one();
                base_variables.push(idx);
                a_idx += 1;
            }
            else {
                row[self.variable_count + s_idx] = if self.equation_compare[i] == Ordering::Less { T::one() } else { T::one().neg() };
                s_idx += 1;

                if (self.equation_right[i] < T::zero() && self.equation_compare[i] == Ordering::Less) ||
                    (self.equation_right[i] >= T::zero() && self.equation_compare[i] == Ordering::Greater) {
                    row[self.variable_count + s_count + a_idx] = T::one();
                    base_variables.push(self.variable_count + s_count + a_idx);
                    a_idx += 1;
                }
                else {
                    base_variables.push(self.variable_count + s_idx - 1);
                }
            }
        }


        let mut right_coef = self.equation_right.clone();

        let mut t = vec![T::zero(); self.variable_count + s_count + a_count];

        let mut right_coef_z = T::zero();
        let col_count = self.variable_count + s_count + a_count;
        for row in 0..(s_count) {
            if table[row][(self.variable_count + s_count)..col_count].iter().any(|x| !x.is_zero()) {
                for col in 0..(self.variable_count + s_count) {
                    t[col] = t[col] - table[row][col];
                }
                right_coef_z = right_coef_z - right_coef[row];
            }
        }
        right_coef.push(right_coef_z);

        table.push(t);

        let z_index = table.len() - 1;

        (SimplexModule::new(table, right_coef, z_index, base_variables), s_count, a_count)
    }

    fn transform_simplex_module_step2(&self, simplex_module: &mut SimplexModule<T>, slack_count: usize, artificial_count: usize) {
        for j in 0..(self.variable_count + slack_count) {
            simplex_module.table[simplex_module.z_index][j] = T::zero();
        }
        simplex_module.right_coef[simplex_module.z_index] = T::zero();
        for j in 0..self.variable_count {
            simplex_module.table[simplex_module.z_index][j] = -self.objective_coef[j];
        }

        for i in 0..self.equation_count {
            if simplex_module.base_variables[i] < self.variable_count {
                simplex_module.right_coef[simplex_module.z_index] = simplex_module.right_coef[simplex_module.z_index] + self.objective_coef[simplex_module.base_variables[i]] * simplex_module.right_coef[i];
            }
            else if simplex_module.base_variables[i] < self.variable_count + slack_count {

            }
            else {
                simplex_module.right_coef[simplex_module.z_index] = simplex_module.right_coef[simplex_module.z_index] + -T::one() * simplex_module.right_coef[i];
            }
        }

        for j in 0..(self.variable_count + slack_count) {
            for i in 0..self.equation_count {
                if simplex_module.base_variables[i] < self.variable_count {
                    simplex_module.table[simplex_module.z_index][j] = simplex_module.table[simplex_module.z_index][j] + self.objective_coef[simplex_module.base_variables[i]] * simplex_module.table[i][j];
                }
                else if simplex_module.base_variables[i] < self.variable_count + slack_count {

                }
                else {
                    simplex_module.table[simplex_module.z_index][j] = simplex_module.table[simplex_module.z_index][j] + -T::one() * simplex_module.table[i][j];
                }
            }
        }
        for i in 0..self.equation_count+1 {
            for l in 0..artificial_count {
                simplex_module.table[i].pop();
            }
        }
    }
}


struct SimplexModule<T: num::Float> {
    table: Vec<Vec<T>>,
    right_coef: Vec<T>,
    z_index: usize,
    base_variables: Vec<usize>,
}

#[derive(Debug)]
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
    fn new(table: Vec<Vec<T>>, right_coef: Vec<T>, z_index: usize, base_variables: Vec<usize>) -> Self {
        SimplexModule { table, right_coef, z_index, base_variables }
    }

    fn solve(&mut self) -> SolveState {
        while let Some(col_idx) = self.table[self.z_index].iter().enumerate().filter(|x| x.1 < &T::zero()).next().map(|x|x.0) {
            let row_idx = self.table.iter()
                .enumerate()
                .map(|x| (x.0, if x.1[col_idx] < T::zero() { T::infinity() } else { self.right_coef[x.0] / x.1[col_idx] }))
                .filter(|x| x.0 != self.z_index)
                .min_by(|x, y| x.1.partial_cmp(&y.1).unwrap_or(Ordering::Equal))
                .map(|x| x.0);
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
                self.table[i] = self.table[i]
                    .iter()
                    .zip(self.table[row_idx].clone())
                    .map(|x| *x.0 - (coef * x.1))
                    //.map(|x| if x.abs() > T::epsilon() { x } else { T::zero() })
                    .collect();
                self.right_coef[i] = self.right_coef[i] - coef * self.right_coef[row_idx];
            }

            self.base_variables[row_idx] = col_idx;
        }

        SolveState::Optimal
    }
}

impl Debug for SimplexModule<f64> {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        writeln!(f, "objective={} ", self.right_coef[self.z_index]);
        for (base, v) in self.base_variables.iter().zip(self.right_coef.clone()) {
            writeln!(f, "x_{}={}", base, v);
        }
        Ok(())
    }
}


fn main() {
    let mut mip = MIP {
        variable_count: 2,
        equation_count: 3,
        objective_coef: vec![5.0, 4.0],
        equation_left: vec![
            vec![1.5, 3.0],
            vec![3.0, 1.0],
            vec![1.0, 2.0],
        ],
        equation_compare: vec![Ordering::Less, Ordering::Less, Ordering::Greater],
        equation_right: vec![ 13.5, 10.0, 7.0 ],
        integer_flag: vec![true, true],
    };
    let solution = mip.solve();
    println!("{:?}", solution);
}
