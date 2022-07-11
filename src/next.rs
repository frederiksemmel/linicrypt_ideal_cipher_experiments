use na::*;
use nalgebra as na;

type Query = DMatrix<u8>;
type Answer = Matrix1xX<u8>;
type Nonce = Vec<u8>;

pub struct ROMConstraint {
    pub t: Nonce,
    pub q: Query,
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

    pub fn deterministically_solvable(fixing: DMatrix<u8>, ordering: Vec<usize>) -> bool {
        todo!()
    }
}

pub struct AlgebraicRepresentation {
    pub i: DMatrix<u8>,
    pub o: DMatrix<u8>,
    pub c: Constraints,
}

impl AlgebraicRepresentation {}
