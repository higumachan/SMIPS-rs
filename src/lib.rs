use std::cmp::Ordering;
use std::fmt;
use std::fmt::{Debug, Formatter};


#[derive(Clone)]
struct Equation<F: num::Float> {
    left_coefs: Vec<F>,
    right: F,
    compare: Ordering,
}

#[derive(Clone)]
struct MIP<F: num::Float> {
    objective_coef: Vec<F>,
    equations: Vec<Equation<F>>,
    integer_flag: Vec<bool>,
}

#[derive(Debug, PartialEq, Clone)]
struct LPSolution<T: num::Float> {
    variables: Vec<T>,
    objective: T,
}

#[derive(PartialEq)]
struct IPSolution<F: num::Float, I: num::Integer> {
    variables: Vec<I>,
    objective: F,
}

enum SolveStateIP<F: num::Float, I: num::Integer> {
    Optimal(IPSolution<F, I>),
    Infeasible,
}

#[derive(PartialEq, Clone, Debug)]
enum SolveStateLP<F :num::Float> {
    Optimal(LPSolution<F>),
    Infeasible,
}

enum Direction {
    Left,
    Right,
}

impl<F: num::Float> PartialOrd for SolveStateLP<F> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        match (self, other) {
            (SolveStateLP::Infeasible, SolveStateLP::Infeasible) => None,
            (SolveStateLP::Optimal(_), SolveStateLP::Infeasible) => Some(Ordering::Greater),
            (SolveStateLP::Infeasible, SolveStateLP::Optimal(_)) => Some(Ordering::Less),
            (SolveStateLP::Optimal(a), SolveStateLP::Optimal(b)) => a.objective.partial_cmp(&b.objective),
        }
    }
}

impl<F: num::Float> SolveStateLP<F> {
    pub fn solution(&self) -> &LPSolution<F> {
        match self {
            SolveStateLP::Optimal(r) => r,
            SolveStateLP::Infeasible => panic!(),
        }
    }

    pub fn is_integer_condition(&self, integer_flags: &Vec<bool>) -> bool {
        match self {
            SolveStateLP::Infeasible => false,
            SolveStateLP::Optimal(r) => r.variables.iter().zip(integer_flags).all(|(x, integer_flag)| !*integer_flag || is_integer(*x)),
        }
    }
}

fn is_integer<F: num::Float>(f: F) -> bool {
    (f - f.round()) < F::epsilon()
}

impl<F: num::Float> MIP<F> {
    pub fn solve(&self) -> SolveStateLP<F> {
        self.branch_bound(&SolveStateLP::Infeasible)
    }

    fn branch_bound(&self, provisional: &SolveStateLP<F>) -> SolveStateLP<F> {
        let answer_lp = self.solve_lp();

        if answer_lp == SolveStateLP::Infeasible {
            return provisional.clone();
        }

        if answer_lp.partial_cmp(provisional) == Some(Ordering::Less) {
            return provisional.clone();
        }

        if answer_lp.is_integer_condition(&self.integer_flag) {
            return answer_lp;
        }

        let (branch_right, branch_left) = self.make_branches(answer_lp.solution()).unwrap();
        let right_provisional = branch_right.branch_bound(provisional);

        branch_left.branch_bound(provisional)
    }

    fn make_branches(&self, solution: &LPSolution<F>) -> Option<(Self, Self)> {
        solution.variables
            .iter()
            .zip(self.integer_flag.iter())
            .enumerate()
            .filter(|(i, (v, is_int))| **is_int && !is_integer(**v))
            .map(|(i, (v, is_integer))| i)
            .next()
            .map(|branch_variable_idx| {
                let mut left_branch = self.clone();
                let mut right_branch = self.clone();

                left_branch.add_equation_with_variable_index(branch_variable_idx, solution.variables[branch_variable_idx].floor(), Ordering::Less);
                right_branch.add_equation_with_variable_index(branch_variable_idx, solution.variables[branch_variable_idx].ceil(), Ordering::Greater);
                (
                    left_branch,
                    right_branch,
                )
            })
    }

    pub fn add_equation(&mut self, equation: Equation<F>) -> &Self {
        self.equations.push(equation);
        self
    }

    pub fn add_equation_with_variable_index(&mut self, variable_index: usize, equation_right: F, equation_compare: Ordering) -> &Self {
        let left_coefs = vec![F::zero(); self.variable_count()];
        self.add_equation(Equation{ left_coefs, right: equation_right, compare: equation_compare })
    }

    pub fn solve_lp(&self) -> SolveStateLP<F> {
        let (mut sm, slack_count, artificial_count) = self.build_simplex_module_step1();
        let st = sm.solve();

        if let SolveState::Infeasible = st {
            return SolveStateLP::Infeasible;
        }

        self.transform_simplex_module_step2(&mut sm, slack_count, artificial_count);
        let st = sm.solve();

        match st {
            SolveState::Infeasible => return SolveStateLP::Infeasible,
            _ => {},
        }


        SolveStateLP::Optimal(self.build_solution(&sm))
    }

    fn build_solution(&self, simplex_module: &SimplexModule<F>) -> LPSolution<F> {
        let mut variables = vec![F::zero(); self.variable_count()];

        for i in 0..self.equation_count() {
            if 0 <= simplex_module.base_variables[i] && simplex_module.base_variables[i] < self.variable_count() {
                variables[simplex_module.base_variables[i]] = simplex_module.right_coef[i];
            }
        }

        LPSolution { variables, objective: simplex_module.right_coef[simplex_module.z_index] }
    }

    fn build_simplex_module_step1(&self) -> (SimplexModule<F>, usize, usize) {
        let mut s_count = 0usize;
        let mut a_count = 0usize;


        for i in 0..self.equation_count() {
            if self.equations[i].compare == Ordering::Equal {
                a_count += 1;
            }
            else {
                s_count += 1;
                if (self.equations[i].right < F::zero() && self.equations[i].compare == Ordering::Less) ||
                    (self.equations[i].right >= F::zero() && self.equations[i].compare == Ordering::Greater) {
                    a_count += 1;
                }
            }
        }


        let mut table: Vec<Vec<F>> = vec![];
        let mut s_idx = 0;
        let mut a_idx = 0;
        let mut base_variables = vec![];


        for i in 0..self.equation_count() {
            table.push(self.equations[i].left_coefs.clone());
            let row = table.last_mut().unwrap();
            row.append(&mut vec![F::zero(); s_count + a_count]);
            if self.equations[i].compare == Ordering::Equal {
                let idx = self.variable_count() + s_count + a_idx;
                row[idx] = F::one();
                base_variables.push(idx);
                a_idx += 1;
            }
            else {
                row[self.variable_count() + s_idx] = if self.equations[i].compare == Ordering::Less { F::one() } else { F::one().neg() };
                s_idx += 1;

                if (self.equations[i].right < F::zero() && self.equations[i].compare == Ordering::Less) ||
                    (self.equations[i].right >= F::zero() && self.equations[i].compare == Ordering::Greater) {
                    row[self.variable_count() + s_count + a_idx] = F::one();
                    base_variables.push(self.variable_count() + s_count + a_idx);
                    a_idx += 1;
                }
                else {
                    base_variables.push(self.variable_count() + s_idx - 1);
                }
            }
        }


        let mut right_coef:Vec<F> = self.equations.iter().map(|x| x.right).collect();

        let mut t = vec![F::zero(); self.variable_count() + s_count + a_count];

        let mut right_coef_z = F::zero();
        let col_count = self.variable_count() + s_count + a_count;
        for row in 0..(s_count) {
            if table[row][(self.variable_count() + s_count)..col_count].iter().any(|x| !x.is_zero()) {
                for col in 0..(self.variable_count() + s_count) {
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

    fn transform_simplex_module_step2(&self, simplex_module: &mut SimplexModule<F>, slack_count: usize, artificial_count: usize) {
        for j in 0..(self.variable_count() + slack_count) {
            simplex_module.table[simplex_module.z_index][j] = F::zero();
        }
        simplex_module.right_coef[simplex_module.z_index] = F::zero();
        for j in 0..self.variable_count() {
            simplex_module.table[simplex_module.z_index][j] = -self.objective_coef[j];
        }

        for i in 0..self.equation_count() {
            if simplex_module.base_variables[i] < self.variable_count() {
                simplex_module.right_coef[simplex_module.z_index] = simplex_module.right_coef[simplex_module.z_index] + self.objective_coef[simplex_module.base_variables[i]] * simplex_module.right_coef[i];
            }
            else if simplex_module.base_variables[i] < self.variable_count() + slack_count {

            }
            else {
                simplex_module.right_coef[simplex_module.z_index] = simplex_module.right_coef[simplex_module.z_index] + -F::one() * simplex_module.right_coef[i];
            }
        }

        for j in 0..(self.variable_count() + slack_count) {
            for i in 0..self.equation_count() {
                if simplex_module.base_variables[i] < self.variable_count() {
                    simplex_module.table[simplex_module.z_index][j] = simplex_module.table[simplex_module.z_index][j] + self.objective_coef[simplex_module.base_variables[i]] * simplex_module.table[i][j];
                }
                else if simplex_module.base_variables[i] < self.variable_count() + slack_count {

                }
                else {
                    simplex_module.table[simplex_module.z_index][j] = simplex_module.table[simplex_module.z_index][j] + -F::one() * simplex_module.table[i][j];
                }
            }
        }
        for i in 0..self.equation_count()+1 {
            for l in 0..artificial_count {
                simplex_module.table[i].pop();
            }
        }
    }

    pub fn equation_count(&self) -> usize {
        self.equations.len()
    }

    pub fn variable_count(&self) -> usize {
        self.objective_coef.len()
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

#[cfg(test)]
mod tests {
    use super::MIP;
    use std::cmp::Ordering;
    use crate::{SolveStateLP, LPSolution, Equation};

    #[test]
    fn solve_lp() {
        let mut mip = MIP {
            objective_coef: vec![5.0, 4.0],
            equations: vec![
                Equation{left_coefs: vec![1.5, 3.0], right: 13.5, compare: Ordering::Less},
                Equation{left_coefs: vec![3.0, 1.0], right: 10.0, compare: Ordering::Less},
                Equation{left_coefs: vec![1.0, 2.0], right: 7.0, compare: Ordering::Greater},
            ],
            integer_flag: vec![false, false],
        };
        let st = mip.solve();

        assert_eq!(SolveStateLP::Optimal(LPSolution{variables: vec![2.2000000000000006, 3.4], objective: 24.600000000000005}), st);
    }

    #[test]
    fn solve_mip1() {
        let mut mip = MIP {
            objective_coef: vec![100.0, 4.0],
            equations: vec![
                Equation{left_coefs: vec![1.0, 1.0], right: 1.0, compare: Ordering::Equal},
                Equation{left_coefs: vec![1.0, 0.0], right: 1.0, compare: Ordering::Less},
                Equation{left_coefs: vec![1.0, 0.0], right: 0.0, compare: Ordering::Greater},
                Equation{left_coefs: vec![0.0, 1.0], right: 1.0, compare: Ordering::Less},
                Equation{left_coefs: vec![0.0, 1.0], right: 0.0, compare: Ordering::Greater},
            ],
            integer_flag: vec![false, false],
        };
        let st = mip.solve();

        assert_eq!(SolveStateLP::Optimal(LPSolution{variables: vec![1.0, 0.0], objective: 100.0}), st)
    }

    #[test]
    fn solve_mip2() {
        let mut mip = MIP {
            objective_coef: vec![100.0, 400.0],
            equations: vec![
                Equation{left_coefs: vec![1.0, 1.0], right: 1.0, compare: Ordering::Equal},
                Equation{left_coefs: vec![1.0, 0.0], right: 1.0, compare: Ordering::Less},
                Equation{left_coefs: vec![1.0, 0.0], right: 0.0, compare: Ordering::Greater},
                Equation{left_coefs: vec![0.0, 1.0], right: 1.0, compare: Ordering::Less},
                Equation{left_coefs: vec![0.0, 1.0], right: 0.0, compare: Ordering::Greater},
            ],
            integer_flag: vec![false, false],
        };
        let st = mip.solve();

        assert_eq!(SolveStateLP::Optimal(LPSolution { variables: vec![0.0, 1.0], objective: 400.0 }), st)
    }
}
