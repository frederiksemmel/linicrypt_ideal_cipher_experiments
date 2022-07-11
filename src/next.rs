use na::*;
use nalgebra as na;

// type Query = DMatrix<u8>;
type Query<const B: usize> = Matrix<Dynamic, Const<B>, u8, VecStorage<u8, Dynamic, Const<B>>>;
type Answer = Matrix1xX<u8>;
type Nonce = Vec<u8>;

pub struct ROMConstraint<const B: usize> {
    pub t: Nonce,
    pub q: Query<B>,
    pub a: Answer,
}

pub struct Constraints(Vec<ROMConstraint>);

impl Constraints {
    pub fn new(c: Vec<ROMConstraint>) -> Self {
        let c = Constraints(c);
        if !c.well_defined() {
            panic!("Needs to be well defined!")
        }
        c
    }

    pub fn well_defined(&self) -> bool {
        // TODO
        true
    }
    
    pub fn deterministically_solvable(fixing: DMatrix<)
}

pub struct AlgebraicRepresentation {
    pub i: DMatrix<u8>,
    pub o: DMatrix<u8>,
    pub c: Constraints,
}

impl AlgebraicRepresentation {
    
}
