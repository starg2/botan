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
         W carry = bigint_add<W, N>(s, lo, hi);

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
         W carry = bigint_add<W, N>(t, a.value(), b.value());

         std::array<W, N> r;
         bigint_monty_maybe_sub<N>(r.data(), carry, t.data(), P.data());
         return Self(r);
      }

      Self dbl() const {
         std::array<W, N> t = value();
         W carry = shift_left<1>(t);

         std::array<W, N> r;
         bigint_monty_maybe_sub<N>(r.data(), carry, t.data(), P.data());
         return Self(r);
      }

      constexpr Self& operator+=(const Self& other) {
         std::array<W, N> t;
         W carry = bigint_add<W, N>(t, this->value(), other.value());
         bigint_monty_maybe_sub<N>(m_val.data(), carry, t.data(), P.data());
         return (*this);
      }

      friend constexpr Self operator-(const Self& a, const Self& b) { return a + b.negate(); }

      friend constexpr Self operator*(uint8_t a, const Self& b) { return b * a; }

      friend constexpr Self operator*(const Self& a, uint8_t b) {
         // We assume b is a small constant and allow variable time
         // computation

         // In practice this function is called for 2, 3, 4, or 8

         if(b == 2) {
            return a.dbl();
         } else if(b == 3) {
            return a.dbl() + a;
         } else if(b == 4) {
            return a.dbl().dbl();
         } else if(b == 8) {
            return a.dbl().dbl().dbl();
         } else {
            Self z = Self::zero();
            Self x = a;

            while(b > 0) {
               if(b & 1) {
                  z = z + x;
               }
               x = x.dbl();
               b >>= 1;
            }

            return z;
         }
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
         constexpr size_t WindowBits = 4;
         constexpr size_t WindowElements = (1 << WindowBits) - 1;

         constexpr size_t Windows = (Self::BITS + WindowBits - 1) / WindowBits;

         std::array<Self, WindowElements> tbl;

         tbl[0] = (*this);

         for(size_t i = 1; i != WindowElements; ++i) {
            if(i % 2 == 1) {
               tbl[i] = tbl[i / 2].square();
            } else {
               tbl[i] = tbl[i - 1] * tbl[0];
            }
         }

         auto r = Self::one();

         const size_t w0 = read_window_bits<WindowBits>(std::span{exp}, (Windows - 1) * WindowBits);

         if(w0 > 0) {
            r = tbl[w0 - 1];
         }

         for(size_t i = 1; i != Windows; ++i) {
            for(size_t j = 0; j != WindowBits; ++j) {
               r = r.square();
            }

            const size_t w = read_window_bits<WindowBits>(std::span{exp}, (Windows - i - 1) * WindowBits);

            if(w > 0) {
               r *= tbl[w - 1];
            }
         }

         return r;
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
         const bool correct = (z.square() == *this);
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

      // ECDSA style hash->scalar conversion
      //
      // This must accept inputs of any length
      static Self from_bits_with_trunc(std::span<const uint8_t> bytes) {
         const size_t bit_length = bytes.size() * 8;

         if(bit_length <= Self::BITS) {
            // No shifting required, but might still need to reduce by modulus
            std::array<uint8_t, 2 * BYTES> padded_bytes = {};
            copy_mem(padded_bytes.data() + 2 * BYTES - bytes.size(), bytes.data(), bytes.size());
            return Self(Rep::wide_to_rep(bytes_to_words<W, 2 * N, 2 * BYTES>(&padded_bytes[0])));
         } else {
            const size_t shift = bit_length - Self::BITS;

            if(shift % 8 == 0) {
               // Easy case just copy different bytes
               const size_t new_length = bytes.size() - (shift / 8);
               return Self::from_bits_with_trunc(bytes.subspan(0, new_length));
            } else {
               // fixme
               throw Not_Implemented("Bit shifting for hash to scalar conversion not implemented");
            }
         }
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
         std::array<uint8_t, Self::BYTES> buf;

         for(;;) {
            rng.randomize(buf);

            // Zero off high bits that if set would certainly cause us
            // to be out of range
            if constexpr(Self::BITS % 8 != 0) {
               constexpr uint8_t mask = 0xFF >> (8 - (Self::BITS % 8));
               buf[0] &= mask;
            }

            if(auto s = Self::deserialize(buf)) {
               if(!s.value().is_zero()) {
                  return s.value();
               }
            }
         }
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

      secure_vector<uint8_t> serialize_x_coordinate() const {
         auto xb = x().serialize();
         return secure_vector<uint8_t>(xb.begin(), xb.end());
      }

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

      static constexpr FieldElement x3_ax_b(const FieldElement& x) { return (x.square() + Self::A) * x + Self::B; }

      static constexpr std::optional<Self> deserialize(std::span<const uint8_t> bytes) {
         if(bytes.size() == Self::BYTES) {
            if(bytes[0] != 0x04) {
               return {};
            }
            auto x = FieldElement::deserialize(bytes.subspan(1, FieldElement::BYTES));
            auto y = FieldElement::deserialize(bytes.subspan(1 + FieldElement::BYTES, FieldElement::BYTES));

            if(x && y) {
               if((*y).square() == Self::x3_ax_b(*x)) {
                  return Self(*x, *y);
               }
            }

            return {};
         } else if(bytes.size() == Self::COMPRESSED_BYTES) {
            if(bytes[0] != 0x02 && bytes[0] != 0x03) {
               return {};
            }
            const bool y_is_even = (bytes[0] == 0x02);

            if(auto x = FieldElement::deserialize(bytes.subspan(1, FieldElement::BYTES))) {
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
         const auto a_is_identity = a.is_identity();
         const auto b_is_identity = b.is_identity();
         if(a_is_identity && b_is_identity) {
            return Self::identity();
         }

         /*
         https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-3.html#addition-add-1998-cmo-2

         12M + 4S + 6add + 1*2
         */

         const auto Z1Z1 = a.z().square();
         const auto U2 = b.x() * Z1Z1;
         const auto S2 = b.y() * a.z() * Z1Z1;
         const auto H = U2 - a.x();
         const auto r = S2 - a.y();

         // If r is zero then we are in the doubling case
         if(r.is_zero()) {
            return a.dbl();
         }

         const auto HH = H.square();
         const auto HHH = H * HH;
         const auto V = a.x() * HH;
         const auto t2 = r.square();
         const auto t3 = V + V;
         const auto t4 = t2 - HHH;
         auto X3 = t4 - t3;
         const auto t5 = V - X3;
         const auto t6 = a.y() * HHH;
         const auto t7 = r * t5;
         auto Y3 = t7 - t6;
         auto Z3 = a.z() * H;

         // TODO these could be combined
         // if a is identity then return b
         X3.conditional_assign(a_is_identity, b.x());
         Y3.conditional_assign(a_is_identity, b.y());
         Z3.conditional_assign(a_is_identity, FieldElement::one());

         // if b is identity then return a
         X3.conditional_assign(b_is_identity, a.x());
         Y3.conditional_assign(b_is_identity, a.y());
         Z3.conditional_assign(b_is_identity, a.z());

         return Self(X3, Y3, Z3);
      }

      constexpr static Self add(const Self& a, const Self& b) {
         const auto a_is_identity = a.is_identity();
         const auto b_is_identity = b.is_identity();
         if(a_is_identity && b_is_identity) {
            return Self::identity();
         }

         /*
         https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-3.html#addition-add-1998-cmo-2
         */

         const auto Z1Z1 = a.z().square();
         const auto Z2Z2 = b.z().square();
         const auto U1 = a.x() * Z2Z2;
         const auto U2 = b.x() * Z1Z1;
         const auto S1 = a.y() * b.z() * Z2Z2;
         const auto S2 = b.y() * a.z() * Z1Z1;
         const auto H = U2 - U1;
         const auto r = S2 - S1;

         if(r.is_zero()) {
            return a.dbl();
         }

         const auto HH = H.square();
         const auto HHH = H * HH;
         const auto V = U1 * HH;
         const auto t2 = r.square();
         const auto t3 = V + V;
         const auto t4 = t2 - HHH;
         auto X3 = t4 - t3;
         const auto t5 = V - X3;
         const auto t6 = S1 * HHH;
         const auto t7 = r * t5;
         auto Y3 = t7 - t6;
         const auto t8 = b.z() * H;
         auto Z3 = a.z() * t8;

         // TODO these could be combined
         // if a is identity then return b
         X3.conditional_assign(a_is_identity, b.x());
         Y3.conditional_assign(a_is_identity, b.y());
         Z3.conditional_assign(a_is_identity, b.z());

         // if b is identity then return a
         X3.conditional_assign(b_is_identity, a.x());
         Y3.conditional_assign(b_is_identity, a.y());
         Z3.conditional_assign(b_is_identity, a.z());

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
      typedef typename Params::W W;

      static const constexpr auto PW = Params::PW;
      static const constexpr auto NW = Params::NW;
      static const constexpr auto AW = Params::AW;

      // Simplifying assumption
      static_assert(PW.size() == NW.size());

      class ScalarParams final : public IntParams<W, NW.size(), NW> {};

      typedef IntMod<MontgomeryRep<ScalarParams>> Scalar;

      class FieldParams final : public IntParams<W, PW.size(), PW> {};

      typedef IntMod<FieldRep<FieldParams>> FieldElement;

      typedef AffineCurvePoint<FieldElement, Params> AffinePoint;
      typedef ProjectiveCurvePoint<FieldElement, Params> ProjectivePoint;

      static const constinit size_t OrderBits = Scalar::BITS;
      static const constinit size_t PrimeFieldBits = FieldElement::BITS;

      static const constexpr FieldElement A = FieldElement::from_words(Params::AW);
      static const constexpr FieldElement B = FieldElement::from_words(Params::BW);

      static const constexpr AffinePoint G =
         AffinePoint(FieldElement::from_words(Params::GXW), FieldElement::from_words(Params::GYW));

      static const constexpr FieldElement SSWU_Z = FieldElement::constant(Params::Z);

      static const constinit bool ValidForSswuHash =
         (SSWU_Z.is_nonzero() && A.is_nonzero() && B.is_nonzero() && FieldElement::P_MOD_4 == 3);

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

      // Use 1/3 the order, rounded up to the next word for blinding
      static const constinit size_t BlindingBits =
         ((C::OrderBits / 3 + WordInfo<W>::bits - 1) / WordInfo<W>::bits) * WordInfo<W>::bits;

      static_assert(BlindingBits % WordInfo<W>::bits == 0);
      static_assert(BlindingBits < C::Scalar::BITS);

   public:
      static constexpr size_t Bits = C::Scalar::BITS + BlindingBits;

      BlindedScalarBits(const typename C::Scalar& scalar, RandomNumberGenerator& rng) {
         constexpr size_t mask_words = BlindingBits / WordInfo<W>::bits;
         constexpr size_t mask_bytes = mask_words * WordInfo<W>::bytes;

         constexpr size_t n_words = C::NW.size();

         uint8_t maskb[mask_bytes] = {0};
         rng.randomize(maskb, mask_bytes);

         W mask[n_words] = {0};
         load_le(mask, maskb, mask_words);
         mask[mask_words-1] |= WordInfo<W>::top_bit;

         W mask_n[2 * n_words] = {0};

         const auto sw = scalar.to_words();

         // Compute masked scalar s + k*n
         comba_mul<n_words>(mask_n, mask, C::NW.data());
         bigint_add2_nc(mask_n, 2 * n_words, sw.data(), sw.size());

         std::reverse(mask_n, mask_n + 2 * n_words);
         m_bytes = store_be<secure_vector<uint8_t>>(mask_n);
      }

      // Extract a WindowBits sized window out of s, depending on offset.
      size_t get_window(size_t offset) const { return read_window_bits<WindowBits>(std::span{m_bytes}, offset); }

   private:
      secure_vector<uint8_t> m_bytes;
};

template <typename C, size_t WindowBits>
class UnblindedScalarBits final {
   public:
      static constexpr size_t Bits = C::Scalar::BITS;

      UnblindedScalarBits(const typename C::Scalar& scalar) { m_bytes = scalar.serialize_to_vec(); }

      // Extract a WindowBits sized window out of s, depending on offset.
      size_t get_window(size_t offset) const { return read_window_bits<WindowBits>(std::span{m_bytes}, offset); }

   private:
      std::vector<uint8_t> m_bytes;
};

template <typename C, size_t W>
class PrecomputedBaseMulTable final {
   public:
      typedef typename C::Scalar Scalar;
      typedef typename C::AffinePoint AffinePoint;
      typedef typename C::ProjectivePoint ProjectivePoint;

      static const constinit size_t WindowBits = W;
      static_assert(WindowBits >= 1 && WindowBits <= 8);

      typedef BlindedScalarBits<C, WindowBits> BlindedScalar;

      static const constinit size_t Windows = (BlindedScalar::Bits + WindowBits - 1) / WindowBits;

      static_assert(Windows > 1);

      // 2^W elements, less the identity element
      static const constinit size_t WindowElements = (1 << WindowBits) - 1;

      static const constinit size_t TableSize = Windows * WindowElements;

      PrecomputedBaseMulTable(const AffinePoint& p) : m_table{} {
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

         auto accum = [&]() {
            const size_t w_i = bits.get_window(0);
            const auto tbl_i = std::span{m_table.begin(), WindowElements};
            auto pt = ProjectivePoint::from_affine(AffinePoint::ct_select(tbl_i, w_i));
            pt.randomize_rep(rng);
            return pt;
         }();

         for(size_t i = 1; i != Windows; ++i) {
            const size_t w_i = bits.get_window(WindowBits * i);
            const auto tbl_i = std::span{m_table.begin() + WindowElements * i, WindowElements};

            /*
            None of these additions can be doublings, because in each iteration, the
            discrete logarithms of the points we're selecting out of the table are
            larger than the largest possible dlog of accum.
            */
            accum += AffinePoint::ct_select(tbl_i, w_i);

            if(i <= 3) {
               accum.randomize_rep(rng);
            }
         }

         return accum;
      }

   private:
      std::vector<AffinePoint> m_table;
};

template <typename C, size_t W>
class WindowedMulTable final {
   public:
      typedef typename C::Scalar Scalar;
      typedef typename C::AffinePoint AffinePoint;
      typedef typename C::ProjectivePoint ProjectivePoint;

      static const constinit size_t WindowBits = W;
      static_assert(WindowBits >= 1 && WindowBits <= 8);

      typedef BlindedScalarBits<C, WindowBits> BlindedScalar;

      static const constinit size_t Windows = (BlindedScalar::Bits + WindowBits - 1) / WindowBits;

      static_assert(Windows > 1);

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

         auto accum = [&]() {
            const size_t w_0 = bits.get_window((Windows - 1) * WindowBits);
            // Guaranteed because we set the high bit of the randomizer
            BOTAN_DEBUG_ASSERT(w_0 != 0);
            auto pt = ProjectivePoint::from_affine(AffinePoint::ct_select(m_table, w_0));
            pt.randomize_rep(rng);
            return pt;
         }();

         for(size_t i = 1; i != Windows; ++i) {
            accum = accum.dbl_n(WindowBits);
            const size_t w_i = bits.get_window((Windows - i - 1) * WindowBits);

            /*
            It's impossible for this addition to ever be a doubling.

            Consider the sequence of points that are operated on, and specifically
            their discrete logarithms. We start out at at the point at infinity
            (dlog 0) and then add the initial window which is precisely P*w_0

            We then perform WindowBits doublings, so accum's dlog at the point
            of the addition in the first iteration of the loop (when i == 1) is
            at least 2^W * w_0.

            Since we know w_0 > 0, then in every iteration of the loop, accums
            dlog will always be greater than the dlog of the table element we
            just looked up (something between 0 and 2^W-1), and thus the
            addition into accum cannot be a doubling.
            */
            accum += AffinePoint::ct_select(m_table, w_i);

            if(i <= 3) {
               accum.randomize_rep(rng);
            }
         }

         return accum;
      }

   private:
      std::vector<AffinePoint> m_table;
};

template <typename C, size_t W>
class WindowedMul2Table final {
   public:
      // We look at 2*W bits of scalar per iteration
      static_assert(W >= 1 && W <= 4);

      typedef typename C::Scalar Scalar;
      typedef typename C::AffinePoint AffinePoint;
      typedef typename C::ProjectivePoint ProjectivePoint;

      static const constinit size_t WindowBits = W;

      static const constinit size_t Windows = (Scalar::BITS + WindowBits - 1) / WindowBits;

      static const constinit size_t WindowSize = (1 << WindowBits);

      // 2^(2*W) elements, less the identity element
      static const constinit size_t TableSize = (1 << (2 * WindowBits)) - 1;

      WindowedMul2Table(const AffinePoint& x, const AffinePoint& y) {
         std::vector<ProjectivePoint> table;
         table.reserve(TableSize);

         for(size_t i = 0; i != TableSize; ++i) {
            const size_t t_i = (i + 1);
            const size_t x_i = t_i % WindowSize;
            const size_t y_i = (t_i >> WindowBits) % WindowSize;

            auto next_tbl_e = [&]() {
               if(x_i % 2 == 0 && y_i % 2 == 0) {
                  return table[(t_i / 2) - 1].dbl();
               } else if(x_i > 0 && y_i > 0) {
                  return table[x_i - 1] + table[(y_i << WindowBits) - 1];
               } else if(x_i > 0 && y_i == 0) {
                  if(x_i == 1) {
                     return ProjectivePoint::from_affine(x);
                  } else {
                     return table[x_i - 1 - 1] + x;
                  }
               } else if(x_i == 0 && y_i > 0) {
                  if(y_i == 1) {
                     return ProjectivePoint::from_affine(y);
                  } else {
                     return table[((y_i - 1) << WindowBits) - 1] + y;
                  }
               } else {
                  BOTAN_ASSERT_UNREACHABLE();
               }
            };

            table.push_back(next_tbl_e());
         }

         m_table = ProjectivePoint::to_affine_batch(table);
      }

      ProjectivePoint mul2_vartime(const Scalar& s1, const Scalar& s2) const {
         const UnblindedScalarBits<C, W> bits1(s1);
         const UnblindedScalarBits<C, W> bits2(s2);

         auto accum = ProjectivePoint::identity();

         for(size_t i = 0; i != Windows; ++i) {
            if(i > 0) {
               accum = accum.dbl_n(WindowBits);
            }

            const size_t w_1 = bits1.get_window((Windows - i - 1) * WindowBits);
            const size_t w_2 = bits2.get_window((Windows - i - 1) * WindowBits);

            const size_t window = w_1 + (w_2 << WindowBits);

            if(window > 0) {
               accum += m_table[window - 1];
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
   const auto gx1 = C::AffinePoint::x3_ax_b(x1);

   const auto x2 = C::SSWU_Z * u.square() * x1;
   const auto gx2 = C::AffinePoint::x3_ax_b(x2);

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

   if(Cnt == 1) {
      const auto u = C::FieldElement::from_wide_bytes(std::span<const uint8_t, L>(xmd));
      return C::ProjectivePoint::from_affine(map_to_curve_sswu<C>(u));
   } else {
      const auto u0 = C::FieldElement::from_wide_bytes(std::span<const uint8_t, L>(xmd.data(), L));
      const auto u1 = C::FieldElement::from_wide_bytes(std::span<const uint8_t, L>(xmd.data() + L, L));

      auto accum = C::ProjectivePoint::from_affine(map_to_curve_sswu<C>(u0));
      accum += map_to_curve_sswu<C>(u1);
      return accum;
   }
}

}  // namespace Botan

#endif
