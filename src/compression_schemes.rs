use itertools::iproduct;
use na::*;
use nalgebra as na;

type Row3 = RowVector3<u8>;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SchemeType {
    Degenerate,
    A,
    B,
    Secure,
}

pub struct SingleQueryScheme {
    m: Row3,
    k: Row3,
    x: Row3,
    y: Row3,
}
impl SingleQueryScheme {
    fn is_y_unconstrained(&self) -> bool {
        let matrix = na::Matrix3::from_rows(&[self.m, self.k, self.x]).cast::<f64>();
        let linear_combination = matrix.lu().solve(&self.y.cast().transpose());
        linear_combination.is_none()
    }
    fn is_x_unconstrained(&self) -> bool {
        let matrix = na::Matrix3::from_rows(&[self.m, self.k, self.y]).cast::<f64>();
        let linear_combination = matrix.lu().solve(&self.x.cast().transpose());
        linear_combination.is_none()
    }
    pub fn collision_structure_type(&self) -> SchemeType {
        let y_free = self.is_y_unconstrained();
        let x_free = self.is_x_unconstrained();
        match (y_free, x_free) {
            (true, true) => SchemeType::Degenerate,
            (true, false) => SchemeType::A,
            (false, true) => SchemeType::B,
            (false, false) => SchemeType::Secure,
        }
    }
}

pub fn generate_all_schemes() -> Vec<SingleQueryScheme> {
    let ms = [(0, 0, 1), (0, 1, 1), (1, 0, 1), (1, 1, 1)].map(|(a, b, c)| Row3::new(a, b, c));
    let ks = [(0, 0, 0), (0, 1, 0), (1, 0, 0), (1, 1, 0)].map(|(a, b, c)| Row3::new(a, b, c));
    let xs = [(0, 0, 0), (0, 1, 0), (1, 0, 0), (1, 1, 0)].map(|(a, b, c)| Row3::new(a, b, c));
    let y = Row3::new(0, 0, 1);

    iproduct!(ms, ks, xs)
        .map(|(m, k, x)| SingleQueryScheme { m, k, x, y })
        .collect()
}

// pub fn generate_non_constant_schemes() -> Vec<SingleQueryScheme> {
//     let ms = [(0, 1, 1), (1, 0, 1), (1, 1, 1)].map(|(a, b, c)| Row3::new(a, b, c));
//     let ks = [(0, 1, 0), (1, 0, 0), (1, 1, 0)].map(|(a, b, c)| Row3::new(a, b, c));
//     let xs = [(0, 1, 0), (1, 0, 0), (1, 1, 0)].map(|(a, b, c)| Row3::new(a, b, c));
//     let y = Row3::new(0, 0, 1);

//     iproduct!(ms, ks, xs)
//         .map(|(m, k, x)| SingleQueryScheme { m, k, x, y })
//         .collect()
// }
pub fn scheme_to_lines(scheme: &SingleQueryScheme) -> Vec<String> {
    let m_line = format!("M={}{}{}", scheme.m[0], scheme.m[1], scheme.m[2]);
    let k_line = format!("k={}{}{}", scheme.k[0], scheme.k[1], scheme.k[2]);
    let x_line = format!("x={}{}{}", scheme.x[0], scheme.x[1], scheme.x[2]);
    let y_line = format!("y={}{}{}", scheme.y[0], scheme.y[1], scheme.y[2]);
    let type_line = format!("{:?}", scheme.collision_structure_type());
    vec![m_line, k_line, x_line, y_line, type_line]
}
