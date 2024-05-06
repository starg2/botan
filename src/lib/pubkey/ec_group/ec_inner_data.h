/*
* (C) 2024 Jack Lloyd
*
* Botan is released under the Simplified BSD License (see license.txt)
*/

#ifndef BOTAN_EC_INNER_DATA_H_
#define BOTAN_EC_INNER_DATA_H_

#include <botan/ec_group.h>

#include <botan/asn1_obj.h>
#include <botan/bigint.h>
#include <botan/reducer.h>
#include <botan/internal/point_mul.h>

namespace Botan {

class EC_Scalar;
class EC_AffinePoint;

class EC_Scalar_Data {
   public:
      virtual ~EC_Scalar_Data() = default;

      virtual std::unique_ptr<EC_Scalar_Data> clone() const = 0;
};

class EC_Scalar_Data_BigInt final : public EC_Scalar_Data {
   public:
      EC_Scalar_Data_BigInt(BigInt v) : m_v(std::move(v)) {}

      std::unique_ptr<EC_Scalar_Data> clone() const override { return std::make_unique<EC_Scalar_Data_BigInt>(m_v); }

      const BigInt& value() const { return m_v; }

      void set_value(const BigInt& v) { m_v = v; }

   private:
      BigInt m_v;
};

class EC_Point_Data final {
   public:
      EC_Point_Data(EC_Point&& pt) : m_pt(std::move(pt)) {
         m_pt.force_affine();
         m_xy = m_pt.xy_bytes();
      }

      EC_Point_Data(const EC_Point& pt) : m_pt(pt) {
         m_pt.force_affine();
         m_xy = m_pt.xy_bytes();
      }

      const EC_Point& value() const { return m_pt; }

      secure_vector<uint8_t> x_bytes() const {
         const size_t half = m_xy.size() / 2;
         return secure_vector<uint8_t>(m_xy.begin(), m_xy.begin() + half);
      }

      secure_vector<uint8_t> y_bytes() const {
         const size_t half = m_xy.size() / 2;
         return secure_vector<uint8_t>(m_xy.begin() + half, m_xy.end());
      }

      const secure_vector<uint8_t>& xy_bytes() const { return m_xy; }

   private:
      EC_Point m_pt;
      secure_vector<uint8_t> m_xy;
};

class EC_Mul2Table_Data {
   public:
      virtual ~EC_Mul2Table_Data() = default;

      virtual std::optional<EC_AffinePoint> mul2(const EC_Scalar& x, const EC_Scalar& y) const = 0;

      virtual std::optional<EC_Scalar> mul2_x_mod_order(const EC_Scalar& x, const EC_Scalar& y) const = 0;
};

class EC_Mul2Table_Data_BigInt final : public EC_Mul2Table_Data {
   public:
      EC_Mul2Table_Data_BigInt(const EC_AffinePoint& h);

      std::optional<EC_AffinePoint> mul2(const EC_Scalar& x, const EC_Scalar& y) const override;

      std::optional<EC_Scalar> mul2_x_mod_order(const EC_Scalar& x, const EC_Scalar& y) const override;
   private:
      std::shared_ptr<EC_Group_Data> m_group;
      EC_Point_Multi_Point_Precompute m_tbl;
};

class EC_Group_Data final {
   public:
      EC_Group_Data(const BigInt& p,
                    const BigInt& a,
                    const BigInt& b,
                    const BigInt& g_x,
                    const BigInt& g_y,
                    const BigInt& order,
                    const BigInt& cofactor,
                    const OID& oid,
                    EC_Group_Source source) :
            m_curve(p, a, b),
            m_base_point(m_curve, g_x, g_y),
            m_g_x(g_x),
            m_g_y(g_y),
            m_order(order),
            m_cofactor(cofactor),
            m_mod_order(order),
            m_base_mult(m_base_point, m_mod_order),
            m_oid(oid),
            m_p_bits(p.bits()),
            m_order_bits(order.bits()),
            m_order_bytes((m_order_bits + 7) / 8),
            m_a_is_minus_3(a == p - 3),
            m_a_is_zero(a.is_zero()),
            m_source(source) {}

      bool params_match(const BigInt& p,
                        const BigInt& a,
                        const BigInt& b,
                        const BigInt& g_x,
                        const BigInt& g_y,
                        const BigInt& order,
                        const BigInt& cofactor) const {
         return (this->p() == p && this->a() == a && this->b() == b && this->order() == order &&
                 this->cofactor() == cofactor && this->g_x() == g_x && this->g_y() == g_y);
      }

      bool params_match(const EC_Group_Data& other) const {
         return params_match(
            other.p(), other.a(), other.b(), other.g_x(), other.g_y(), other.order(), other.cofactor());
      }

      void set_oid(const OID& oid) {
         BOTAN_STATE_CHECK(m_oid.empty());
         m_oid = oid;
      }

      const OID& oid() const { return m_oid; }

      const BigInt& p() const { return m_curve.get_p(); }

      const BigInt& a() const { return m_curve.get_a(); }

      const BigInt& b() const { return m_curve.get_b(); }

      const BigInt& order() const { return m_order; }

      const BigInt& cofactor() const { return m_cofactor; }

      const BigInt& g_x() const { return m_g_x; }

      const BigInt& g_y() const { return m_g_y; }

      size_t p_bits() const { return m_p_bits; }

      size_t p_bytes() const { return (m_p_bits + 7) / 8; }

      size_t order_bits() const { return m_order_bits; }

      size_t order_bytes() const { return (m_order_bits + 7) / 8; }

      const CurveGFp& curve() const { return m_curve; }

      const EC_Point& base_point() const { return m_base_point; }

      bool a_is_minus_3() const { return m_a_is_minus_3; }

      bool a_is_zero() const { return m_a_is_zero; }

      BigInt mod_order(const BigInt& x) const { return m_mod_order.reduce(x); }

      BigInt square_mod_order(const BigInt& x) const { return m_mod_order.square(x); }

      BigInt multiply_mod_order(const BigInt& x, const BigInt& y) const { return m_mod_order.multiply(x, y); }

      BigInt multiply_mod_order(const BigInt& x, const BigInt& y, const BigInt& z) const {
         return m_mod_order.multiply(m_mod_order.multiply(x, y), z);
      }

      BigInt inverse_mod_order(const BigInt& x) const { return inverse_mod(x, m_order); }

      EC_Point blinded_base_point_multiply(const BigInt& k, RandomNumberGenerator& rng, std::vector<BigInt>& ws) const {
         return m_base_mult.mul(k, rng, m_order, ws);
      }

      EC_Group_Source source() const { return m_source; }

      std::unique_ptr<EC_Mul2Table_Data> make_mul2_table(const EC_AffinePoint& h) const {
         return std::make_unique<EC_Mul2Table_Data_BigInt>(h);
      }

      std::unique_ptr<EC_Scalar_Data> scalar_from_bytes_with_trunc(std::span<const uint8_t> bytes) const {
         auto bn = BigInt::from_bytes_with_max_bits(bytes.data(), bytes.size(), m_order_bits);
         return std::make_unique<EC_Scalar_Data_BigInt>(mod_order(bn));
      }

      std::unique_ptr<EC_Scalar_Data> scalar_from_bytes_mod_order(std::span<const uint8_t> bytes) const {
         return std::make_unique<EC_Scalar_Data_BigInt>(mod_order(BigInt(bytes)));
      }

      std::unique_ptr<EC_Scalar_Data> scalar_random(RandomNumberGenerator& rng) const {
         return std::make_unique<EC_Scalar_Data_BigInt>(BigInt::random_integer(rng, BigInt::one(), m_order));
      }

      bool scalar_is_zero(const EC_Scalar_Data& s) const { return as_bn(s).value().is_zero(); }

      bool scalar_is_eq(const EC_Scalar_Data& x, const EC_Scalar_Data& y) const { return as_bn(x).value() == as_bn(y).value(); }

      std::unique_ptr<EC_Scalar_Data> scalar_zero() const {
         return std::make_unique<EC_Scalar_Data_BigInt>(BigInt::zero());
      }

      std::unique_ptr<EC_Scalar_Data> scalar_one() const {
         return std::make_unique<EC_Scalar_Data_BigInt>(BigInt::one());
      }

      std::unique_ptr<EC_Scalar_Data> scalar_invert(const EC_Scalar_Data& s) const {
         return std::make_unique<EC_Scalar_Data_BigInt>(inverse_mod_order(as_bn(s).value()));
      }

      void scalar_assign(EC_Scalar_Data& x, const EC_Scalar_Data& y) { as_bn(x).set_value(as_bn(y).value()); }

      void scalar_square_self(EC_Scalar_Data& s) { as_bn(s).set_value(square_mod_order(as_bn(s).value())); }

      std::unique_ptr<EC_Scalar_Data> scalar_negate(const EC_Scalar_Data& s) const {
         return std::make_unique<EC_Scalar_Data_BigInt>(mod_order(-as_bn(s).value()));
      }

      std::unique_ptr<EC_Scalar_Data> scalar_add(const EC_Scalar_Data& a, const EC_Scalar_Data& b) const {
         return std::make_unique<EC_Scalar_Data_BigInt>(mod_order(as_bn(a).value() + as_bn(b).value()));
      }

      std::unique_ptr<EC_Scalar_Data> scalar_sub(const EC_Scalar_Data& a, const EC_Scalar_Data& b) const {
         return std::make_unique<EC_Scalar_Data_BigInt>(mod_order(as_bn(a).value() - as_bn(b).value()));
      }

      std::unique_ptr<EC_Scalar_Data> scalar_mul(const EC_Scalar_Data& a, const EC_Scalar_Data& b) const {
         return std::make_unique<EC_Scalar_Data_BigInt>(multiply_mod_order(as_bn(a).value(), as_bn(b).value()));
      }

      std::unique_ptr<EC_Scalar_Data> scalar_from_bigint(const BigInt& bn) const {
         // Assumed to have been already checked as in range
         return std::make_unique<EC_Scalar_Data_BigInt>(bn);
      }

      std::unique_ptr<EC_Scalar_Data> gk_x_mod_order(const EC_Scalar_Data& scalar,
                                                     RandomNumberGenerator& rng,
                                                     std::vector<BigInt>& ws) const {
         const auto pt = m_base_mult.mul(as_bn(scalar).value(), rng, m_order, ws);

         if(pt.is_zero()) {
            return scalar_zero();
         } else {
            return std::make_unique<EC_Scalar_Data_BigInt>(mod_order(pt.get_affine_x()));
         }
      }

      std::unique_ptr<EC_Scalar_Data> scalar_deserialize(std::span<const uint8_t> bytes) {
         if(bytes.size() != m_order_bytes) {
            return nullptr;
         }

         BigInt r(bytes.data(), bytes.size());

         if(r.is_zero() || r >= m_order) {
            return nullptr;
         }

         return std::make_unique<EC_Scalar_Data_BigInt>(std::move(r));
      }

      std::vector<uint8_t> scalar_serialize(const EC_Scalar_Data& s) const {
         std::vector<uint8_t> bytes(m_order_bytes);
         as_bn(s).value().binary_encode(bytes.data(), m_order_bytes);
         return bytes;
      }

      std::vector<uint8_t> scalar_serialize_pair(const EC_Scalar_Data& r, const EC_Scalar_Data& s) const {
         std::vector<uint8_t> bytes(2 * m_order_bytes);
         as_bn(r).value().binary_encode(bytes.data(), m_order_bytes);
         as_bn(s).value().binary_encode(bytes.data() + m_order_bytes, m_order_bytes);
         return bytes;
      }

   private:
      const EC_Scalar_Data_BigInt& as_bn(const EC_Scalar_Data& s) const {
         try {
            return dynamic_cast<const EC_Scalar_Data_BigInt&>(s);
         }
         catch(std::bad_cast&) {
            throw Invalid_State("Bad cast to EC_Scalar_Data_BigInt");
         }
      }

      EC_Scalar_Data_BigInt& as_bn(EC_Scalar_Data& s) const {
         try {
            return dynamic_cast<EC_Scalar_Data_BigInt&>(s);
         }
         catch(std::bad_cast&) {
            throw Invalid_State("Bad cast to EC_Scalar_Data_BigInt");
         }
      }

      CurveGFp m_curve;
      EC_Point m_base_point;

      BigInt m_g_x;
      BigInt m_g_y;
      BigInt m_order;
      BigInt m_cofactor;
      Modular_Reducer m_mod_order;
      EC_Point_Base_Point_Precompute m_base_mult;
      OID m_oid;
      size_t m_p_bits;
      size_t m_order_bits;
      size_t m_order_bytes;
      bool m_a_is_minus_3;
      bool m_a_is_zero;
      EC_Group_Source m_source;
};

}  // namespace Botan

#endif
