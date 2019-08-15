
#[cfg(test)]
mod tests {
    use super::MIP;
    use std::cmp::Ordering;
    use crate::{SolveStateLP, LPSolution};

    #[test]
    fn solve_lp() {
        let mut mip = MIP {
            objective_coef: vec![5.0, 4.0],
            equation_left: vec![
                vec![1.5, 3.0],
                vec![3.0, 1.0],
                vec![1.0, 2.0],
            ],
            equation_compare: vec![Ordering::Less, Ordering::Less, Ordering::Greater],
            equation_right: vec![ 13.5, 10.0, 7.0 ],
            integer_flag: vec![false, false],
        };
        let st = mip.solve();

        assert_eq!(SolveStateLP::Optimal(LPSolution{variables: vec![2.2000000000000006, 3.4], objective: 24.600000000000005}), st);
    }

    #[test]
    fn solve_mip() {

    }
}