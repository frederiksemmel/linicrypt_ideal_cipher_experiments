use na::{ArrayStorage, Const, Dynamic, VecStorage};
use nalgebra as na;

type DConstraint<const B: usize> =
    na::Matrix<Dynamic, Const<B>, u8, VecStorage<u8, Dynamic, Const<B>>>;
type SConstraint<const B: usize, const K: usize> =
    na::Matrix<Const<K>, Const<B>, u8, VecStorage<u8, Dynamic, Const<B>>>;
type Matrix<const B: usize> = DConstraint<B>;
type Row<const B: usize> = SConstraint<B, 1>;
type Transformation<const B: usize> = na::Matrix<Const<B>, Const<B>, u8, ArrayStorage<u8, B, B>>;
type Nonce = Vec<u8>;

// Only B = base can be generic, otherwise need to box it for the set of constraints
pub struct ROMConstraint<const B: usize> {
    pub t: Nonce,
    pub q: Matrix<B>,
    pub a: Row<B>,
}

pub struct Constraints<const B: usize>(Vec<ROMConstraint<B>>);

impl<const B: usize> Constraints<B> {
    pub fn well_defined(&self) -> bool {
        todo!()
    }

    pub fn det_solvable(_fixing: DConstraint<B>, _ordering: Vec<usize>) -> bool {
        todo!()
    }
}

pub struct AlgebraicRepresentation<const B: usize> {
    pub i: Matrix<B>,
    pub o: Matrix<B>,
    pub c: Constraints<B>,
}

impl<const B: usize> AlgebraicRepresentation<B> {
    pub fn det_invertible(&self) -> bool {
        todo!()
    }

    pub fn has_collision_structure(&self, _coll: &CollisionStructure<B>) -> bool {
        todo!()
    }

    pub fn find_collision_structure(&self) -> Option<CollisionStructure<B>> {
        todo!()
    }
}

pub struct CollisionStructure<const B: usize> {
    // this contains the indices of the divergent constraints
    pub divergent: Vec<usize>,
    pub q_star: Matrix<B>,
}
