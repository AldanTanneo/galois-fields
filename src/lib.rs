#![feature(once_cell)]
pub mod optimized;

use std::fmt::{Debug, Display};

/// A binary polynomial represented using a 64-bits unsigned integer.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct Polynomial(pub u64);

impl Polynomial {
    pub const SYMBOL: &'static str = "X";

    pub const ZERO: Polynomial = Polynomial(0);
    pub const ONE: Polynomial = Polynomial(1);
    pub const X: Polynomial = Polynomial(0b10);

    #[inline]
    /// Construct a Polynomial from a list of exponents
    pub fn new<const N: usize>(pol: [u8; N]) -> Self {
        Polynomial(pol.into_iter().fold(0, |res, n| res ^ (1 << n)))
    }

    #[inline]
    /// Construct the X^n polynomial
    pub const fn x_pow(n: u32) -> Self {
        Polynomial(1 << n)
    }

    #[inline]
    /// Compute the degree of a non zero polynomial
    pub const fn degree(self) -> u32 {
        assert!(self.0 != 0, "Cannot compute the degree of 0");
        // cannot truncate: degree will always be less than 64
        self.0.ilog2()
    }

    /// Get the list of exponents of non-zero terms
    pub fn to_exponents(mut self) -> Vec<u8> {
        let mut vec = Vec::new();
        while self.0 != 0 {
            let pow = self.0.trailing_zeros();
            vec.push(pow as u8);
            self.0 -= 1 << pow;
        }
        vec
    }

    /// Compute the euclidean division of `self` by `other`
    ///
    /// Returns `(q, r)` where `self = q * other + r`, with `r.degree() < other.degree()`
    pub fn div_euclid(self, other: Self) -> (Self, Self) {
        assert_ne!(other, Self::ZERO, "Cannot perform a division by 0");

        let (mut q, mut r) = (Self::ZERO, self);
        while r != Self::ZERO && r.degree() >= other.degree() {
            let delta = r.degree() - other.degree();
            q += Polynomial(1 << delta);
            r += Polynomial(other.0 << delta);
        }
        (q, r)
    }

    /// Compute the gcd of `self` and `other` using the extended Euclid algorithm.
    ///
    /// Returns `(gcd, s, t)`, where `s * self + t * other = gcd`
    pub fn extended_gcd(self, other: Self) -> (Self, Self, Self) {
        let mut a = self;
        let mut b = other;

        let mut s = [Self::ONE, Self::ZERO];
        let mut t = [Self::ZERO, Self::ONE];

        while b != Self::ZERO {
            let (q, r) = a.div_euclid(b);

            (a, b) = (b, r);
            t = [t[1], t[0] - q * t[1]];
            s = [s[1], s[0] - q * s[1]];
        }

        (a, s[0], t[0])
    }

    /// X^8 + X^4 + X^3 + X + 1
    pub const POLY256: Self = Polynomial(0b100011011);

    #[inline]
    /// Perform addition in GF(8)
    pub fn add256(self, other: Self) -> Self {
        self + other // addition in GF(8) is just addition of two polynomials
    }

    #[inline]
    /// Perform multiplication of `self` and `other` in GF(8)
    /// using X^8 + X^4 + X^3 + X + 1 for reduction
    pub fn mul256(self, other: Self) -> Self {
        // can span up to 16 bits
        let full_multiplication = self * other;

        // reduce using the remainder of the euclidean division
        let (_q, r) = full_multiplication.div_euclid(Self::POLY256);

        r
    }

    #[inline]
    /// Compute the inverse of a polynomial in GF(8).
    pub fn inv256(self) -> Self {
        debug_assert!(self.degree() < 8);
        // use the extended euclid algorithm
        // since `X^8 + X^4 + X^3 + X + 1` is irreducible, and `self.degree() < 8`,
        // `_gcd = 1`, and thus `Self::POLY256 * s + self * t = 1`
        let (_gcd, _s, t) = Self::POLY256.extended_gcd(self);
        t
    }

    #[inline]
    /// Divide `self` by `other` in GF(8)
    pub fn div256(self, other: Self) -> Self {
        other.inv256().mul256(self)
    }
}

impl From<&str> for Polynomial {
    /// Implement string parsing for polynomials.
    /// Panics when provided an invalid binary string.
    #[inline]
    fn from(value: &str) -> Self {
        Polynomial(u64::from_str_radix(value, 2).expect("Invalid binary value"))
    }
}

impl From<u8> for Polynomial {
    #[inline]
    fn from(value: u8) -> Self {
        Polynomial(u64::from(value))
    }
}

impl Display for Polynomial {
    /// Implement full form display for polynomials
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.0 == 0 {
            write!(f, "0")
        } else {
            let mut x = self.0;
            while x != 0 {
                if x != self.0 {
                    write!(f, " + ")?;
                }
                let pow = x.ilog2();
                match pow {
                    0 => write!(f, "1"),
                    1 => write!(f, "{}", Self::SYMBOL),
                    _ => write!(f, "{}^{}", Self::SYMBOL, pow),
                }?;
                x ^= 1 << pow;
            }

            Ok(())
        }
    }
}

impl std::fmt::Binary for Polynomial {
    /// Implement binary display for polynomials
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Binary::fmt(&self.0, f)
    }
}

impl Debug for Polynomial {
    /// Implement Debug display for polynomials
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_tuple("Polynomial")
            .field(&format_args!("{:#b}", self.0))
            .finish()
    }
}

impl std::ops::Add for Polynomial {
    type Output = Self;

    /// Overload the `+` operator
    #[inline]
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0 ^ rhs.0)
    }
}

impl std::ops::AddAssign for Polynomial {
    /// Overload the `+=` operator
    #[inline]
    #[allow(clippy::suspicious_op_assign_impl)]
    fn add_assign(&mut self, rhs: Self) {
        self.0 ^= rhs.0
    }
}

impl std::ops::Sub for Polynomial {
    type Output = Self;

    /// Overload the `-` operator
    ///
    /// Substraction is just addition since we are dealing with binary coefficients
    #[inline]
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn sub(self, rhs: Self) -> Self::Output {
        self + rhs
    }
}

impl std::ops::SubAssign for Polynomial {
    /// Overload the `-=` operator
    #[inline]
    #[allow(clippy::suspicious_op_assign_impl)]
    fn sub_assign(&mut self, rhs: Self) {
        *self += rhs
    }
}

impl std::ops::Mul for Polynomial {
    type Output = Self;

    /// Overload the `*` operator
    fn mul(self, rhs: Self) -> Self::Output {
        // On x86/x64: use CPU specific instructions
        //
        // The `pclmulqdq` CPU flag is enabled on the vast majority of modern x86/x64 processors.
        #[cfg(any(target_arch = "x86_64", target_arch = "x86"))]
        {
            #[cfg(target_arch = "x86")]
            use core::arch::x86::{_mm_clmulepi64_si128, _mm_cvtsi128_si64, _mm_cvtsi64_si128};
            #[cfg(target_arch = "x86_64")]
            use core::arch::x86_64::{_mm_clmulepi64_si128, _mm_cvtsi128_si64, _mm_cvtsi64_si128};

            // SAFETY: use of x86_64 specific instructions
            //
            // _mm_clmulepi64_si128 performs a "carry-less multiplication of two 64-bits polynomials"
            // in SIMD registers
            Self(unsafe {
                _mm_cvtsi128_si64(_mm_clmulepi64_si128::<0>(
                    _mm_cvtsi64_si128(self.0 as i64),
                    _mm_cvtsi64_si128(rhs.0 as i64),
                ))
            } as u64)
        }

        // On other architectures: compute the multiplication by hand
        #[cfg(not(any(target_arch = "x86_64", target_arch = "x86")))]
        {
            let (mut a, mut b) = (self.0, rhs.0);
            let mut res = Polynomial(0);
            while a != 0 && b != 0 {
                res.0 ^= (a & 1) * b;
                b <<= 1;
                a >>= 1;
            }
            res
        }
    }
}

impl std::ops::MulAssign for Polynomial {
    /// Overload the  `*=` operator
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl std::ops::Rem for Polynomial {
    type Output = Self;

    /// Overload the `%` operator
    #[inline]
    fn rem(self, rhs: Self) -> Self::Output {
        self.div_euclid(rhs).1
    }
}

impl<T> std::ops::Shl<T> for Polynomial
where
    u64: std::ops::Shl<T, Output = u64>,
{
    type Output = Self;

    fn shl(self, rhs: T) -> Self::Output {
        Polynomial(self.0 << rhs)
    }
}

/// A few unit tests
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_degree() {
        let a = Polynomial::from("1101");
        let b = Polynomial::from("1");
        let c = Polynomial::from("10000");

        assert_eq!(a.degree(), 3);
        assert_eq!(b.degree(), 0);
        assert_eq!(c.degree(), 4);
    }

    #[test]
    fn test_div_euclid() {
        let a = Polynomial::from("10011");
        let b = Polynomial::from("111");

        assert_eq!(
            a.div_euclid(b),
            (Polynomial::from("110"), Polynomial::from("1"))
        )
    }

    #[test]
    fn test_mul256() {
        let a = Polynomial::from("10001100");
        let b = Polynomial::from("101010");

        assert_eq!(a.mul256(b), Polynomial::from("111111"));
    }

    #[test]
    fn test_extended_gcd() {
        let a = Polynomial::from("1010111");
        let b = Polynomial::from("1011");

        let (gcd, s, t) = a.extended_gcd(b);

        assert_eq!(
            (gcd, s, t),
            (
                Polynomial::from("1"),
                Polynomial::from("111"),
                Polynomial::from("111100")
            )
        );

        assert_eq!(gcd, s * a + t * b);
    }

    #[test]
    fn test_inv256() {
        let a = Polynomial::from("111");

        assert_eq!(a.inv256(), Polynomial::from("11010001"))
    }

    #[test]
    fn test_div256() {
        let a = Polynomial::from("100101");
        let b = Polynomial::from("11010");

        assert_eq!(a.div256(b), Polynomial::from("10110110"))
    }
}
