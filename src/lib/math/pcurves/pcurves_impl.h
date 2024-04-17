/*
* (C) 2024 Jack Lloyd
*
* Botan is released under the Simplified BSD License (see license.txt)
*/

#ifndef BOTAN_PCURVES_IMPL_H_
#define BOTAN_PCURVES_IMPL_H_

#include <botan/rng.h>
#include <botan/internal/loadstor.h>
#include <botan/internal/pcurves_util.h>
#include <botan/internal/stl_util.h>
#include <botan/internal/xmd.h>
#include <vector>

namespace Botan {

template <typename Params>
class MontgomeryRep final {
   public:
      typedef MontgomeryRep<Params> Self;

      static const constexpr auto P = Params::P;
      static const constexpr size_t N = Params::N;
      typedef typename Params::W W;

      static_assert(N > 0 && (Params::P[0] & 1) == 1, "Invalid Montgomery modulus");

      static const constinit auto P_dash = monty_inverse(P[0]);

      static const constexpr auto R1 = montygomery_r(P);
      static const constexpr auto R2 = mul_mod(R1, R1, P);
      static const constexpr auto R3 = mul_mod(R1, R2, P);

      constexpr static std::array<W, N> one() { return R1; }

      constexpr static std::array<W, N> redc(const std::array<W, 2 * N>& z) {
         if constexpr(P_dash == 1) {
            return monty_redc_pdash1(z, P);
         } else {
            return monty_redc(z, P, P_dash);
         }
      }

      constexpr static std::array<W, N> to_rep(const std::array<W, N>& x) {
         std::array<W, 2 * N> z;
         comba_mul<N>(z.data(), x.data(), R2.data());
         return Self::redc(z);
      }

      constexpr static std::array<W, N> wide_to_rep(const std::array<W, 2 * N>& x) {
         auto redc_x = Self::redc(x);
         std::array<W, 2 * N> z;
         comba_mul<N>(z.data(), redc_x.data(), R3.data());
         return Self::redc(z);
      }

      constexpr static std::array<W, N> from_rep(const std::array<W, N>& z) {
         std::array<W, 2 * N> ze = {};
         std::copy(z.begin(), z.end(), ze.begin());
         return Self::redc(ze);
      }
};

template <typename Params>
class P521Rep final {
   private:
      template <WordType W, size_t N>
      static consteval bool is_p521(const std::array<W, N>& p) {
         if constexpr(N != (521 + WordInfo<W>::bits - 1) / WordInfo<W>::bits) {
            return false;
         }

         if(p[N - 1] != 0x1FF) {
            return false;
         }
         for(size_t i = 0; i != N - 1; ++i) {
            if(~p[i] != 0) {
               return false;
            }
         }

         return true;
      }

   public:
      static const constexpr auto P = Params::P;
      static const constexpr size_t N = Params::N;
      typedef typename Params::W W;

      static_assert(is_p521(P));

      constexpr static std::array<W, N> one() {
         std::array<W, N> one = {};
         one[0] = 1;
         return one;
      }

      constexpr static std::array<W, N> redc(const std::array<W, 2 * N>& z) {
         constexpr W TOP_MASK = static_cast<W>(0x1FF);

         std::array<W, N> hi = {};
         std::copy(z.begin() + N - 1, z.begin() + 2 * N - 1, hi.begin());
         shift_right<9>(hi);

         std::array<W, N> lo = {};
         std::copy(z.begin(), z.begin() + N, lo.begin());
         lo[N - 1] &= TOP_MASK;

         // s = hi + lo
         std::array<W, N> s = {};
         // Will never carry out
         word carry = bigint_add3_nc(s.data(), lo.data(), N, hi.data(), N);

         // But might be greater than modulus:
         std::array<W, N> r = {};
         bigint_monty_maybe_sub<N>(r.data(), carry, s.data(), P.data());

         return r;
      }

      constexpr static std::array<W, N> to_rep(const std::array<W, N>& x) { return x; }

      constexpr static std::array<W, N> wide_to_rep(const std::array<W, 2 * N>& x) { return redc(x); }

      constexpr static std::array<W, N> from_rep(const std::array<W, N>& z) { return z; }
};

template <typename Rep>
class IntMod final {
   private:
      static const constexpr auto P = Rep::P;
      static const constexpr size_t N = Rep::N;
      typedef typename Rep::W W;

      static const constexpr auto P_MINUS_2 = p_minus<2>(P);
      static const constexpr auto P_PLUS_1_OVER_4 = p_plus_1_over_4(P);
      static const constexpr auto P_MINUS_1_OVER_2 = p_minus_1_over_2(P);

   public:
      static const constexpr size_t BITS = count_bits(P);
      static const constexpr size_t BYTES = (BITS + 7) / 8;

      static const constexpr auto P_MOD_4 = P[0] % 4;

      typedef IntMod<Rep> Self;

      // Default value is zero
      constexpr IntMod() : m_val({}) {}

      IntMod(const Self& other) = default;
      IntMod(Self&& other) = default;
      IntMod& operator=(const Self& other) = default;
      IntMod& operator=(Self&& other) = default;

      // ??
      //~IntMod() { secure_scrub_memory(m_val); }

      static constexpr Self zero() { return Self(std::array<W, N>{0}); }

      static constexpr Self one() { return Self(Rep::one()); }

      static constexpr Self from_word(W x) {
         std::array<W, 1> v{x};
         return Self::from_words(v);
      }

      template <size_t L>
      static constexpr Self from_words(std::array<W, L> w) {
         if constexpr(L == N) {
            return Self(Rep::to_rep(w));
         } else {
            static_assert(L < N);
            std::array<W, N> ew = {};
            for(size_t i = 0; i != L; ++i) {
               ew[i] = w[i];
            }
            return Self(Rep::to_rep(ew));
         }
      }

      constexpr bool is_zero() const { return CT::all_zeros(m_val.data(), m_val.size()).as_bool(); }

      constexpr bool is_nonzero() const { return !is_zero(); }

      constexpr bool is_one() const { return (*this == Self::one()); }

      constexpr bool is_even() const {
         auto v = Rep::from_rep(m_val);
         return (v[0] & 0x01) == 0;
      }

      friend constexpr Self operator+(const Self& a, const Self& b) {
         std::array<W, N> t;
         W carry = bigint_add3_nc(t.data(), a.data(), N, b.data(), N);

         std::array<W, N> r;
         bigint_monty_maybe_sub<N>(r.data(), carry, t.data(), P.data());
         return Self(r);
      }

      constexpr Self& operator+=(const Self& other) {
         std::array<W, N> t;
         W carry = bigint_add3_nc(t.data(), this->data(), N, other.data(), N);
         bigint_monty_maybe_sub<N>(m_val.data(), carry, t.data(), P.data());
         return (*this);
      }

      friend constexpr Self operator-(const Self& a, const Self& b) { return a + b.negate(); }

      friend constexpr Self operator*(uint8_t a, const Self& b) { return b * a; }

      friend constexpr Self operator*(const Self& a, uint8_t b) {
         // We assume b is a small constant and allow variable time
         // computation

         Self z = Self::zero();
         Self x = a;

         while(b > 0) {
            if(b & 1) {
               z = z + x;
            }
            x += x;
            b >>= 1;
         }

         return z;
      }

      friend constexpr Self operator*(const Self& a, const Self& b) {
         std::array<W, 2 * N> z;
         comba_mul<N>(z.data(), a.data(), b.data());
         return Self(Rep::redc(z));
      }

      constexpr Self& operator-=(const Self& other) {
         (*this) += other.negate();
         return (*this);
      }

      constexpr Self& operator*=(const Self& other) {
         std::array<W, 2 * N> z;
         comba_mul<N>(z.data(), data(), other.data());
         m_val = Rep::redc(z);
         return (*this);
      }

      constexpr void conditional_add(bool cond, const Self& other) { conditional_assign(cond, *this + other); }

      constexpr void conditional_mul(bool cond, const Self& other) { conditional_assign(cond, *this * other); }

      constexpr void conditional_sub(bool cond, const Self& other) { conditional_add(cond, other.negate()); }

      // if cond is true, assigns other to *this
      constexpr void conditional_assign(bool cond, const Self& other) {
         CT::conditional_assign_mem(static_cast<W>(cond), m_val.data(), other.data(), N);
      }

      constexpr Self square() const {
         std::array<W, 2 * N> z;
         comba_sqr<N>(z.data(), this->data());
         return Self(Rep::redc(z));
      }

      // Negation modulo p
      constexpr Self negate() const {
         auto x_is_zero = CT::all_zeros(this->data(), N);

         std::array<W, N> r;
         bigint_sub3(r.data(), P.data(), N, this->data(), N);
         x_is_zero.if_set_zero_out(r.data(), N);
         return Self(r);
      }

      constexpr Self pow_vartime(const std::array<W, N>& exp) const {
         auto x = (*this);
         auto y = Self::one();

         for(size_t i = 0; i != Self::BITS; ++i) {
            if(get_bit(i, exp)) {
               y = y * x;
            }
            x = x.square();
         }

         return y;
      }

      /**
      * Returns the modular inverse, or 0 if no modular inverse exists.
      *
      * If the modulus is prime the only value that has no modular inverse is 0.
      *
      * This uses Fermat's little theorem, and so assumes that p is prime
      */
      constexpr Self invert() const { return pow_vartime(Self::P_MINUS_2); }

      /**
      * Return the modular square root, or zero if no root exists
      *
      * Current impl assumes p == 3 (mod 4)
      */
      constexpr Self sqrt() const {
         static_assert(Self::P_MOD_4 == 3);
         auto z = pow_vartime(Self::P_PLUS_1_OVER_4);
         const bool correct = (z * z) == *this;
         z.conditional_assign(!correct, Self::zero());
         return z;
      }

      constexpr bool is_square() const {
         static_assert(Self::P_MOD_4 == 3);
         auto z = pow_vartime(Self::P_MINUS_1_OVER_2);
         const bool is_one = z.is_one();
         const bool is_zero = z.is_zero();
         return (is_one || is_zero);
      }

      constexpr bool operator==(const Self& other) const {
         return CT::is_equal(this->data(), other.data(), N).as_bool();
      }

      constexpr bool operator!=(const Self& other) const {
         return CT::is_not_equal(this->data(), other.data(), N).as_bool();
      }

      constexpr std::array<W, Self::N> to_words() const { return Rep::from_rep(m_val); }

      std::vector<uint8_t> serialize_to_vec() const {
         const auto b = this->serialize();
         return std::vector(b.begin(), b.end());
      }

      constexpr std::array<uint8_t, Self::BYTES> serialize() const {
         auto v = Rep::from_rep(m_val);
         std::reverse(v.begin(), v.end());
         auto bytes = store_be(v);

         if constexpr(Self::BYTES == N * WordInfo<W>::bytes) {
            return bytes;
         } else {
            // Remove leading zero bytes
            constexpr size_t extra = N * WordInfo<W>::bytes - Self::BYTES;
            std::array<uint8_t, Self::BYTES> out;
            copy_mem(out, std::span{bytes}.template subspan<extra, Self::BYTES>());
            return out;
         }
      }

      template <size_t L>
      std::array<W, L> stash_value() const {
         static_assert(L >= N);
         std::array<W, L> stash = {};
         for(size_t i = 0; i != N; ++i) {
            stash[i] = m_val[i];
         }
         return stash;
      }

      template <size_t L>
      static Self from_stash(const std::array<W, L>& stash) {
         static_assert(L >= N);
         std::array<W, N> val = {};
         for(size_t i = 0; i != N; ++i) {
            val[i] = stash[i];
         }
         return Self(val);
      }

      // Returns nullopt if the input is an encoding greater than or equal P
      constexpr static std::optional<Self> deserialize(std::span<const uint8_t> bytes) {
         // We could allow either short inputs or longer zero padded
         // inputs here, however it seems best to avoid non-canonical
         // representations unless required
         if(bytes.size() != Self::BYTES) {
            return {};
         }

         const auto words = bytes_to_words<W, N, BYTES>(&bytes[0]);

         if(!bigint_ct_is_lt(words.data(), N, P.data(), N).as_bool()) {
            return {};
         }

         return Self::from_words(words);
      }

      // Reduces large input modulo the order
      template <size_t L>
      static constexpr Self from_wide_bytes(std::span<const uint8_t, L> bytes) {
         static_assert(8 * L <= 2 * Self::BITS);
         std::array<uint8_t, 2 * BYTES> padded_bytes = {};
         copy_mem(padded_bytes.data() + 2 * BYTES - L, bytes.data(), L);
         return Self(Rep::wide_to_rep(bytes_to_words<W, 2 * N, 2 * BYTES>(&padded_bytes[0])));
      }

      static constexpr Self random(RandomNumberGenerator& rng) {
         const size_t R_bytes = Self::BYTES + 16;
         std::array<uint8_t, R_bytes> buf;
         rng.randomize(buf);
         return Self::from_wide_bytes(std::span<const uint8_t, R_bytes>{buf});
      }

      static consteval Self constant(int8_t x) {
         std::array<W, 1> v;
         v[0] = (x >= 0) ? x : -x;
         auto s = Self::from_words(v);
         return (x >= 0) ? s : s.negate();
      }

   private:
      constexpr const std::array<W, N>& value() const { return m_val; }

      constexpr const W* data() const { return m_val.data(); }

      explicit constexpr IntMod(std::array<W, N> v) : m_val(v) {}

      std::array<W, N> m_val;
};

template <typename FieldElement, typename Params>
class AffineCurvePoint {
   public:
      // We can't pass a FieldElement directly because FieldElement is
      // not "structural" due to having private members, so instead
      // recreate it here from the words.
      static const constexpr FieldElement A = FieldElement::from_words(Params::AW);
      static const constexpr FieldElement B = FieldElement::from_words(Params::BW);

      static const constinit size_t BYTES = 1 + 2 * FieldElement::BYTES;
      static const constinit size_t COMPRESSED_BYTES = 1 + FieldElement::BYTES;

      typedef AffineCurvePoint<FieldElement, Params> Self;

      constexpr AffineCurvePoint(const FieldElement& x, const FieldElement& y) : m_x(x), m_y(y) {}

      constexpr AffineCurvePoint() : m_x(FieldElement::zero()), m_y(FieldElement::zero()) {}

      static constexpr Self identity() { return Self(FieldElement::zero(), FieldElement::zero()); }

      constexpr bool is_identity() const { return x().is_zero() && y().is_zero(); }

      AffineCurvePoint(const Self& other) = default;
      AffineCurvePoint(Self&& other) = default;
      AffineCurvePoint& operator=(const Self& other) = default;
      AffineCurvePoint& operator=(Self&& other) = default;

      constexpr Self negate() const { return Self(x(), y().negate()); }

      std::vector<uint8_t> serialize_to_vec(bool compress) const {
         if(compress) {
            const auto b = this->serialize_compressed();
            return std::vector(b.begin(), b.end());
         } else {
            const auto b = this->serialize();
            return std::vector(b.begin(), b.end());
         }
      }

      constexpr std::array<uint8_t, Self::BYTES> serialize() const {
         std::array<uint8_t, Self::BYTES> r = {};
         BufferStuffer pack(r);
         pack.append(0x04);
         pack.append(x().serialize());
         pack.append(y().serialize());
         return r;
      }

      constexpr std::array<uint8_t, Self::COMPRESSED_BYTES> serialize_compressed() const {
         const bool y_is_even = y().is_even();
         const uint8_t hdr = y_is_even ? 0x02 : 0x03;

         std::array<uint8_t, Self::COMPRESSED_BYTES> r = {};
         BufferStuffer pack(r);
         pack.append(hdr);
         pack.append(x().serialize());
         return r;
      }

      /**
      * If idx is zero then return the identity element. Otherwise return pts[idx - 1]
      *
      * Returns the identity element also if idx is out of range
      */
      static constexpr auto ct_select(std::span<const Self> pts, size_t idx) {
         auto result = Self::identity();

         // Intentionally wrapping; set to maximum size_t if idx == 0
         const size_t idx1 = static_cast<size_t>(idx - 1);
         for(size_t i = 0; i != pts.size(); ++i) {
            const bool found = (idx1 == i);
            result.conditional_assign(found, pts[i]);
         }

         return result;
      }

      static constexpr std::optional<Self> deserialize(std::span<const uint8_t> bytes) {
         auto x3_ax_b = [](const FieldElement& x) { return (x.square() + Self::A) * x + Self::B; };

         auto valid_xy = [x3_ax_b](const FieldElement& x, const FieldElement& y) { return (y.square() == x3_ax_b(x)); };

         if(bytes.size() == Self::BYTES) {
            if(bytes[0] != 0x04) {
               return {};
            }
            auto x = FieldElement::deserialize(bytes.subspan(1, Self::BYTES));
            auto y = FieldElement::deserialize(bytes.subspan(1 + Self::BYTES, Self::BYTES));

            if(x && y) {
               if(valid_xy(*x, *y)) {
                  return Self(*x, *y);
               }
            }

            return {};
         } else if(bytes.size() == Self::COMPRESSED_BYTES) {
            if(bytes[0] != 0x02 && bytes[0] != 0x03) {
               return {};
            }
            const bool y_is_even = (bytes[0] == 0x02);

            if(auto x = FieldElement::deserialize(bytes.subspan(1, Self::BYTES))) {
               const auto y2 = x3_ax_b(*x);
               auto y = y2.sqrt();
               if(y_is_even && !y.is_even()) {
                  y = y.negate();
               }
               return Self(*x, y);
            }

            return {};
         } else {
            return {};
         }
      }

      constexpr const FieldElement& x() const { return m_x; }

      constexpr const FieldElement& y() const { return m_y; }

      constexpr void conditional_assign(bool cond, const Self& pt) {
         m_x.conditional_assign(cond, pt.x());
         m_y.conditional_assign(cond, pt.y());
      }

   private:
      FieldElement m_x;
      FieldElement m_y;
};

template <typename FieldElement, typename Params>
class ProjectiveCurvePoint {
   public:
      // We can't pass a FieldElement directly because FieldElement is
      // not "structural" due to having private members, so instead
      // recreate it here from the words.
      static const constexpr FieldElement A = FieldElement::from_words(Params::AW);

      static const constinit bool A_is_zero = A.is_zero();
      // Should be constinit but this triggers a bug in Clang
      static const constexpr bool A_is_minus_3 = (A == FieldElement::constant(-3));

      typedef ProjectiveCurvePoint<FieldElement, Params> Self;
      typedef AffineCurvePoint<FieldElement, Params> AffinePoint;

      static constexpr Self from_affine(const AffinePoint& pt) {
         if(pt.is_identity()) {
            return Self::identity();
         } else {
            return ProjectiveCurvePoint(pt.x(), pt.y());
         }
      }

      static constexpr Self identity() { return Self(FieldElement::zero(), FieldElement::one(), FieldElement::zero()); }

      constexpr ProjectiveCurvePoint() :
            m_x(FieldElement::zero()), m_y(FieldElement::one()), m_z(FieldElement::zero()) {}

      constexpr ProjectiveCurvePoint(const FieldElement& x, const FieldElement& y) :
            m_x(x), m_y(y), m_z(FieldElement::one()) {}

      constexpr ProjectiveCurvePoint(const FieldElement& x, const FieldElement& y, const FieldElement& z) :
            m_x(x), m_y(y), m_z(z) {}

      ProjectiveCurvePoint(const Self& other) = default;
      ProjectiveCurvePoint(Self&& other) = default;
      ProjectiveCurvePoint& operator=(const Self& other) = default;
      ProjectiveCurvePoint& operator=(Self&& other) = default;

      friend constexpr Self operator+(const Self& a, const Self& b) { return Self::add(a, b); }

      friend constexpr Self operator+(const Self& a, const AffinePoint& b) { return Self::add_mixed(a, b); }

      friend constexpr Self operator+(const AffinePoint& a, const Self& b) { return Self::add_mixed(b, a); }

      constexpr Self& operator+=(const Self& other) {
         (*this) = (*this) + other;
         return (*this);
      }

      constexpr Self& operator+=(const AffinePoint& other) {
         (*this) = (*this) + other;
         return (*this);
      }

      friend constexpr Self operator-(const Self& a, const Self& b) { return a + b.negate(); }

      constexpr bool is_identity() const { return z().is_zero(); }

      template <typename Pt>
      constexpr void conditional_add(bool cond, const Pt& pt) {
         conditional_assign(cond, *this + pt);
      }

      constexpr void conditional_assign(bool cond, const Self& pt) {
         m_x.conditional_assign(cond, pt.x());
         m_y.conditional_assign(cond, pt.y());
         m_z.conditional_assign(cond, pt.z());
      }

      constexpr static Self add_mixed(const Self& a, const AffinePoint& b) {
         // TODO avoid these early returns by masking instead
         if(a.is_identity()) {
            return Self::from_affine(b);
         }

         if(b.is_identity()) {
            return a;
         }

         /*
         https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-3.html#addition-add-1998-cmo-2

         // 12M + 4S + 6add + 1*2
         TODO rename these vars

         TODO reduce vars

         TODO use a complete addition formula??? (YES)
         https://eprint.iacr.org/2015/1060.pdf
         */

         const auto Z1Z1 = a.z().square();
         const auto U2 = b.x() * Z1Z1;
         const auto S2 = b.y() * a.z() * Z1Z1;
         const auto H = U2 - a.x();
         const auto r = S2 - a.y();

         if(H.is_zero()) {
            if(r.is_zero()) {
               return a.dbl();
            } else {
               return Self::identity();
            }
         }

         const auto HH = H.square();
         const auto HHH = H * HH;
         const auto V = a.x() * HH;
         const auto t2 = r.square();
         const auto t3 = V + V;
         const auto t4 = t2 - HHH;
         const auto X3 = t4 - t3;
         const auto t5 = V - X3;
         const auto t6 = a.y() * HHH;
         const auto t7 = r * t5;
         const auto Y3 = t7 - t6;
         const auto Z3 = a.z() * H;

         return Self(X3, Y3, Z3);
      }

      constexpr static Self add(const Self& a, const Self& b) {
         // TODO avoid these early returns by masking instead
         if(a.is_identity()) {
            return b;
         }

         if(b.is_identity()) {
            return a;
         }

         /*
         https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-3.html#addition-add-1998-cmo-2

         TODO rename these vars

         TODO reduce vars

         TODO use a complete addition formula??? (YES)
         https://eprint.iacr.org/2015/1060.pdf
         */

         const auto Z1Z1 = a.z().square();
         const auto Z2Z2 = b.z().square();
         const auto U1 = a.x() * Z2Z2;
         const auto U2 = b.x() * Z1Z1;
         const auto S1 = a.y() * b.z() * Z2Z2;
         const auto S2 = b.y() * a.z() * Z1Z1;
         const auto H = U2 - U1;
         const auto r = S2 - S1;

         if(H.is_zero()) {
            if(r.is_zero()) {
               return a.dbl();
            } else {
               return Self::identity();
            }
         }

         const auto HH = H.square();
         const auto HHH = H * HH;
         const auto V = U1 * HH;
         const auto t2 = r.square();
         const auto t3 = V + V;
         const auto t4 = t2 - HHH;
         const auto X3 = t4 - t3;
         const auto t5 = V - X3;
         const auto t6 = S1 * HHH;
         const auto t7 = r * t5;
         const auto Y3 = t7 - t6;
         const auto t8 = b.z() * H;
         const auto Z3 = a.z() * t8;

         return Self(X3, Y3, Z3);
      }

      constexpr Self dbl_n(size_t n) const {
         Self pt = (*this);
         for(size_t i = 0; i != n; ++i) {
            pt = pt.dbl();
         }
         return pt;
      }

      constexpr Self dbl() const {
         //https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#doubling-dbl-1998-cmo-2

         FieldElement m = FieldElement::zero();

         if constexpr(Self::A_is_minus_3) {
            /*
            if a == -3 then
            3*x^2 + a*z^4 == 3*x^2 - 3*z^4 == 3*(x^2-z^4) == 3*(x-z^2)*(x+z^2)

            Cost: 2M + 2A + 1*3
            */
            const auto z2 = z().square();
            m = 3 * (x() - z2) * (x() + z2);
         } else if constexpr(Self::A_is_zero) {
            // If a == 0 then 3*x^2 + a*z^4 == 3*x^2
            // Cost: 1S + 1*3
            m = 3 * x().square();
         } else {
            // Cost: 1M + 3S + 1A + 1*3
            const auto z2 = z().square();
            m = 3 * x().square() + A * z2.square();
         }

         const auto y2 = y().square();
         const auto s = 4 * x() * y2;
         const auto nx = m.square() - 2 * s;
         const auto ny = m * (s - nx) - 8 * y2.square();
         const auto nz = 2 * y() * z();

         return Self(nx, ny, nz);
      }

      constexpr Self negate() const { return Self(x(), y().negate(), z()); }

      constexpr AffinePoint to_affine() const {
         // Not strictly required right? - default should work as long
         // as (0,0) is identity and invert returns 0 on 0
         if(this->is_identity()) {
            return AffinePoint::identity();
         }

         // Maybe also worth skipping ...
         if(m_z.is_one()) {
            return AffinePoint(m_x, m_y);
         }

         const auto z_inv = m_z.invert();
         const auto z2_inv = z_inv.square();
         const auto z3_inv = z_inv * z2_inv;

         const auto x = m_x * z2_inv;
         const auto y = m_y * z3_inv;
         return AffinePoint(x, y);
      }

      static std::vector<AffinePoint> to_affine_batch(std::span<const Self> projective) {
         const size_t N = projective.size();
         std::vector<AffinePoint> affine(N, AffinePoint::identity());

         bool any_identity = false;
         for(size_t i = 0; i != N; ++i) {
            if(projective[i].is_identity()) {
               any_identity = true;
            }
         }

         if(N <= 2 || any_identity) {
            for(size_t i = 0; i != N; ++i) {
               affine[i] = projective[i].to_affine();
            }
         } else {
            std::vector<FieldElement> c(N);

            /*
            Batch projective->affine using Montgomery's trick

            See Algorithm 2.26 in "Guide to Elliptic Curve Cryptography"
            (Hankerson, Menezes, Vanstone)
            */

            c[0] = projective[0].z();
            for(size_t i = 1; i != N; ++i) {
               c[i] = c[i - 1] * projective[i].z();
            }

            auto s_inv = c[N - 1].invert();

            for(size_t i = N - 1; i > 0; --i) {
               const auto& p = projective[i];

               const auto z_inv = s_inv * c[i - 1];
               const auto z2_inv = z_inv.square();
               const auto z3_inv = z_inv * z2_inv;

               s_inv = s_inv * p.z();

               affine[i] = AffinePoint(p.x() * z2_inv, p.y() * z3_inv);
            }

            const auto z2_inv = s_inv.square();
            const auto z3_inv = s_inv * z2_inv;
            affine[0] = AffinePoint(projective[0].x() * z2_inv, projective[0].y() * z3_inv);
         }

         return affine;
      }

      void randomize_rep(RandomNumberGenerator& rng) {
         auto r = FieldElement::random(rng);

         auto r2 = r.square();
         auto r3 = r2 * r;

         m_x *= r2;
         m_y *= r3;
         m_z *= r;
      }

      constexpr const FieldElement& x() const { return m_x; }

      constexpr const FieldElement& y() const { return m_y; }

      constexpr const FieldElement& z() const { return m_z; }

   private:
      FieldElement m_x;
      FieldElement m_y;
      FieldElement m_z;
};

template <StringLiteral PS,
          StringLiteral AS,
          StringLiteral BS,
          StringLiteral NS,
          StringLiteral GXS,
          StringLiteral GYS,
          int8_t ZI = 0>
class EllipticCurveParameters {
   public:
      typedef word W;

      static const constexpr auto PW = hex_to_words<W>(PS.value);
      static const constexpr auto NW = hex_to_words<W>(NS.value);
      static const constexpr auto AW = hex_to_words<W>(AS.value);
      static const constexpr auto BW = hex_to_words<W>(BS.value);
      static const constexpr auto GXW = hex_to_words<W>(GXS.value);
      static const constexpr auto GYW = hex_to_words<W>(GYS.value);

      static const constexpr int8_t Z = ZI;
};

template <WordType WI, size_t NI, std::array<WI, NI> PI>
struct IntParams {
   public:
      typedef WI W;
      static constexpr size_t N = NI;
      static constexpr auto P = PI;
};

template <typename Params, template <typename FieldParams> typename FieldRep = MontgomeryRep>
class EllipticCurve {
   public:
      typedef word W;

      static const constexpr auto PW = Params::PW;
      static const constexpr auto NW = Params::NW;
      static const constexpr auto AW = Params::AW;

      // Simplifying assumption
      static_assert(PW.size() == NW.size());

      class ScalarParams final : public IntParams<W, NW.size(), NW> {};

      class FieldParams final : public IntParams<W, PW.size(), PW> {};

      typedef IntMod<MontgomeryRep<ScalarParams>> Scalar;
      typedef IntMod<FieldRep<FieldParams>> FieldElement;

      static const constinit size_t OrderBits = Scalar::BITS;
      static const constinit size_t PrimeFieldBits = FieldElement::BITS;

      // Use 1/3 the order, rounded up to the next word for blinding
      static const constinit size_t BlindingBits =
         ((OrderBits / 3 + WordInfo<W>::bits - 1) / WordInfo<W>::bits) * WordInfo<W>::bits;

      static const constexpr FieldElement A = FieldElement::from_words(Params::AW);
      static const constexpr FieldElement B = FieldElement::from_words(Params::BW);
      static const constexpr FieldElement Gx = FieldElement::from_words(Params::GXW);
      static const constexpr FieldElement Gy = FieldElement::from_words(Params::GYW);

      static const constexpr FieldElement SSWU_Z = FieldElement::constant(Params::Z);

      static const constinit bool ValidForSswuHash =
         (SSWU_Z.is_nonzero() && A.is_nonzero() && B.is_nonzero() && FieldElement::P_MOD_4 == 3);

      typedef AffineCurvePoint<FieldElement, Params> AffinePoint;
      typedef ProjectiveCurvePoint<FieldElement, Params> ProjectivePoint;

      static const constexpr AffinePoint G = AffinePoint(Gx, Gy);

      // (-B / A), will be zero if A == 0 or B == 0 or Z == 0
      template <typename = typename std::enable_if<ValidForSswuHash>>
      static const FieldElement& SSWU_C1() {
         // We derive it from C2 to avoid a second inversion
         static const auto C1 = (SSWU_C2() * SSWU_Z).negate();
         return C1;
      }

      // (B / (Z * A)), will be zero if A == 0 or B == 0 or Z == 0
      template <typename = typename std::enable_if<ValidForSswuHash>>
      static const FieldElement& SSWU_C2() {
         // This could use a variable time inversion
         static const auto C2 = (B * (SSWU_Z * A).invert());
         return C2;
      }
};

template <typename C, size_t WindowBits>
class BlindedScalarBits final {
   private:
      typedef typename C::W W;

      static_assert(C::BlindingBits % WordInfo<W>::bits == 0);
      static_assert(C::BlindingBits < C::Scalar::BITS);

      // A bitmask that has the lowest WindowBits bits set
      static const constinit uint8_t WindowMask = static_cast<uint8_t>(1 << WindowBits) - 1;

   public:
      static constexpr size_t Bits = C::Scalar::BITS + C::BlindingBits;

      BlindedScalarBits(const typename C::Scalar& scalar, RandomNumberGenerator& rng) {
         const size_t mask_words = C::BlindingBits / WordInfo<W>::bits;
         const size_t mask_bytes = mask_words * WordInfo<W>::bytes;

         const size_t n_words = C::NW.size();

         uint8_t maskb[mask_bytes] = {0};
         rng.randomize(maskb, mask_bytes);

         W mask[n_words] = {0};
         load_be(mask, maskb, mask_words);

         W mask_n[2 * n_words] = {0};

         const auto sw = scalar.to_words();

         // Compute masked scalar s + k*n
         comba_mul<n_words>(mask_n, mask, C::NW.data());
         bigint_add2_nc(mask_n, 2 * n_words, sw.data(), sw.size());

         std::reverse(mask_n, mask_n + 2 * n_words);
         m_bytes = store_be<std::vector<uint8_t>>(mask_n);
      }

      // Extract a WindowBits sized window out of s, depending on offset.
      size_t get_window(size_t offset) const {
         const auto bit_shift = offset % 8;
         const auto byte_offset = m_bytes.size() - 1 - (offset / 8);

         const bool single_byte_window = bit_shift <= (8 - WindowBits) || byte_offset == 0;

         const auto w0 = m_bytes[byte_offset];

         if(single_byte_window) {
            return (w0 >> bit_shift) & WindowMask;
         } else {
            // Otherwise we must join two bytes and extract the result
            const auto w1 = m_bytes[byte_offset - 1];
            const auto combined = ((w0 >> bit_shift) | (w1 << (8 - bit_shift)));
            return combined & WindowMask;
         }
      }

   private:
      // secure_vector?
      std::vector<uint8_t> m_bytes;
};

template <typename C>
class PrecomputedMulTable final {
   public:
      typedef typename C::Scalar Scalar;
      typedef typename C::AffinePoint AffinePoint;
      typedef typename C::ProjectivePoint ProjectivePoint;

      // TODO allow config?
      static const constinit size_t WindowBits = 4;
      static_assert(WindowBits >= 1 && WindowBits < 8);

      typedef BlindedScalarBits<C, WindowBits> BlindedScalar;

      static const constinit size_t Windows = (BlindedScalar::Bits + WindowBits - 1) / WindowBits;

      // 2^W elements, less the identity element
      static const constinit size_t WindowElements = (1 << WindowBits) - 1;

      static const constinit size_t TableSize = Windows * WindowElements;

      PrecomputedMulTable(const AffinePoint& p) : m_table{} {
         std::vector<ProjectivePoint> table;
         table.reserve(TableSize);

         auto accum = ProjectivePoint::from_affine(p);

         for(size_t i = 0; i != TableSize; i += WindowElements) {
            table.push_back(accum);

            for(size_t j = 1; j != WindowElements; ++j) {
               if(j % 2 == 1) {
                  table.push_back(table[i + j / 2].dbl());
               } else {
                  table.push_back(table[i + j - 1] + table[i]);
               }
            }

            accum = table[i + (WindowElements / 2)].dbl();
         }

         m_table = ProjectivePoint::to_affine_batch(table);
      }

      ProjectivePoint mul(const Scalar& s, RandomNumberGenerator& rng) const {
         const BlindedScalar bits(s, rng);

         auto accum = ProjectivePoint::identity();

         for(size_t i = 0; i != Windows; ++i) {
            const size_t w_i = bits.get_window(WindowBits * i);

            auto tbl_i = std::span{m_table.begin() + WindowElements * i, WindowElements};
            accum += AffinePoint::ct_select(tbl_i, w_i);

            if(i == 0 || i == Windows / 2) {
               accum.randomize_rep(rng);
            }
         }

         return accum;
      }

   private:
      std::vector<AffinePoint> m_table;
};

template <typename C, size_t W = 5>
class WindowedMulTable final {
   public:
      typedef typename C::Scalar Scalar;
      typedef typename C::AffinePoint AffinePoint;
      typedef typename C::ProjectivePoint ProjectivePoint;

      static const constinit size_t WindowBits = W;
      static_assert(WindowBits >= 1 && WindowBits < 6);

      typedef BlindedScalarBits<C, WindowBits> BlindedScalar;

      static const constinit size_t Windows = (BlindedScalar::Bits + WindowBits - 1) / WindowBits;

      // 2^W elements, less the identity element
      static const constinit size_t TableSize = (1 << WindowBits) - 1;

      WindowedMulTable(const AffinePoint& p) : m_table{} {
         std::vector<ProjectivePoint> table;
         table.reserve(TableSize);

         table.push_back(ProjectivePoint::from_affine(p));
         for(size_t i = 1; i != TableSize; ++i) {
            if(i % 2 == 1) {
               table.push_back(table[i / 2].dbl());
            } else {
               table.push_back(table[i - 1] + table[0]);
            }
         }

         m_table = ProjectivePoint::to_affine_batch(table);
      }

      ProjectivePoint mul(const Scalar& s, RandomNumberGenerator& rng) const {
         const BlindedScalar bits(s, rng);

         auto accum = ProjectivePoint::identity();

         for(size_t i = 0; i != Windows; ++i) {
            if(i > 0) {
               accum = accum.dbl_n(WindowBits);
            }
            const size_t w_i = bits.get_window((Windows - i - 1) * WindowBits);

            accum += AffinePoint::ct_select(m_table, w_i);

            if(i == 0 || i == Windows / 2) {
               accum.randomize_rep(rng);
            }
         }

         return accum;
      }

   private:
      std::vector<AffinePoint> m_table;
};

template <typename C>
inline auto map_to_curve_sswu(const typename C::FieldElement& u) -> typename C::AffinePoint {
   const auto z_u2 = C::SSWU_Z * u.square();  // z * u^2
   const auto z2_u4 = z_u2.square();
   const auto tv1 = (z2_u4 + z_u2).invert();
   auto x1 = C::SSWU_C1() * (C::FieldElement::one() + tv1);
   x1.conditional_assign(tv1.is_zero(), C::SSWU_C2());
   const auto gx1 = (x1.square() + C::A) * x1 + C::B;

   const auto x2 = C::SSWU_Z * u.square() * x1;
   const auto gx2 = (x2.square() + C::A) * x2 + C::B;

   const auto gx1_is_square = gx1.is_square();

   auto x = x2;
   auto y = gx2.sqrt();

   x.conditional_assign(gx1_is_square, x1);
   y.conditional_assign(gx1_is_square, gx1.sqrt());

   const bool flip_y = y.is_even() != u.is_even();
   y.conditional_assign(flip_y, y.negate());

   return typename C::AffinePoint(x, y);
}

template <typename C>
inline auto hash_to_curve_sswu(std::string_view hash,
                               bool random_oracle,
                               std::span<const uint8_t> pw,
                               std::span<const uint8_t> dst) -> typename C::ProjectivePoint {
   static_assert(C::ValidForSswuHash);

   const size_t SecurityLevel = (C::OrderBits + 1) / 2;
   const size_t L = (C::PrimeFieldBits + SecurityLevel + 7) / 8;

   const size_t Cnt = (random_oracle ? 2 : 1);

   std::vector<uint8_t> xmd(L * Cnt);
   expand_message_xmd(hash, xmd, pw, dst);

   auto pt = C::ProjectivePoint::identity();

   for(size_t i = 0; i != Cnt; ++i) {
      const auto u = C::FieldElement::from_wide_bytes(std::span<const uint8_t, L>(xmd.data() + i * L, L));
      pt += map_to_curve_sswu<C>(u);
   }

   return pt;
}

}  // namespace Botan

#endif
