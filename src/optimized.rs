//! Optimized functions using static tables

use crate::Polynomial;

/// Primitive element of GF(8) used to generate the log tables.
/// Here, we choose X+1
const GENERATOR: Polynomial = Polynomial(0b11);

/// Generate LOG256 and ALOG256 tables
fn tables_log256() -> ([u8; 256], [u8; 256]) {
    let (mut log256, mut alog256) = ([0; 256], [0; 256]);
    let mut elt = Polynomial::ONE;

    alog256[255] = 1;
    log256[0] = 255; // special case

    for (i, alog) in &mut alog256[0..255].iter_mut().enumerate() {
        log256[elt.0 as usize] = i as u8;
        *alog = elt.0 as u8;
        elt = elt.mul256(GENERATOR); // use the standard implementation of mul256
    }

    (log256, alog256)
}

use std::sync::LazyLock; // use a "lazy" container to generate the tables on first access

// Generate both tables in a single pass
static TABLES: LazyLock<([u8; 256], [u8; 256])> = LazyLock::new(tables_log256); // takes a function argument for delayed initialization

/// Macro for easy log256 access
macro_rules! log256 {
    ($x:expr) => {
        TABLES.0[($x).0 as usize]
    };
}
/// Macro for easy alog256 access
macro_rules! alog256 {
    ($x:expr) => {
        TABLES.1[($x) as usize]
    };
}

/// Perform multiplication of `a` and `b` in GF(8)
/// using static log tables
pub fn mul256(a: Polynomial, b: Polynomial) -> Polynomial {
    if a == Polynomial::ZERO || b == Polynomial::ZERO {
        return Polynomial::ZERO;
    }
    let ai = log256![a];
    let bi = log256![b];
    let ci = (ai as usize + bi as usize).rem_euclid(255);
    Polynomial::from(alog256![ci])
}

/// Perform division of `a` and `b` in GF(8)
/// using static log tables
pub fn div256(a: Polynomial, b: Polynomial) -> Polynomial {
    assert_ne!(b, Polynomial::ZERO, "Cannot divide by 0");
    if a == Polynomial::ZERO {
        return Polynomial::ZERO;
    }
    let ai = log256![a];
    let bi = log256![b];
    let ci = (ai as isize - bi as isize).rem_euclid(255);
    Polynomial::from(alog256![ci])
}

/// Perform inversion of `a` in GF(8)
/// using static log tables
pub fn inv256(a: Polynomial) -> Polynomial {
    assert_ne!(a, Polynomial::ZERO, "Cannot get the inverse of 0");
    let ai = log256![a];
    Polynomial::from(alog256![255 - ai])
}

/// Unit tests for optimized multiplication, division and inversion
#[cfg(test)]
mod tests {
    use super::*;
    use crate::Polynomial;

    #[test]
    fn test_mul256() {
        let a = Polynomial::from("10001100");
        let b = Polynomial::from("101010");

        assert_eq!(mul256(a, b), Polynomial::from("111111"));
    }

    #[test]
    fn test_div256() {
        let a = Polynomial::from("100101");
        let b = Polynomial::from("11010");

        assert_eq!(div256(a, b), Polynomial::from("10110110"))
    }

    #[test]
    fn test_inv256() {
        let a = Polynomial::from("111");

        assert_eq!(inv256(a), Polynomial::from("11010001"))
    }
}
