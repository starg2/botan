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

class EC_Scalar_Data final {
   public:
      EC_Scalar_Data(BigInt v) : m_v(std::move(v)) {}

      std::unique_ptr<EC_Scalar_Data> clone() const { return std::make_unique<EC_Scalar_Data>(m_v); }

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
                    EC_Group_Source source);

      bool params_match(const BigInt& p,
                        const BigInt& a,
                        const BigInt& b,
                        const BigInt& g_x,
                        const BigInt& g_y,
                        const BigInt& order,
                        const BigInt& cofactor) const;

      bool params_match(const EC_Group_Data& other) const;

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

      std::unique_ptr<EC_Scalar_Data> scalar_from_bytes_with_trunc(std::span<const uint8_t> bytes) const;

      std::unique_ptr<EC_Scalar_Data> scalar_from_bytes_mod_order(std::span<const uint8_t> bytes) const;

      std::unique_ptr<EC_Scalar_Data> scalar_random(RandomNumberGenerator& rng) const;

      bool scalar_is_zero(const EC_Scalar_Data& s) const;

      bool scalar_is_eq(const EC_Scalar_Data& x, const EC_Scalar_Data& y) const;

      std::unique_ptr<EC_Scalar_Data> scalar_zero() const;

      std::unique_ptr<EC_Scalar_Data> scalar_one() const;

      std::unique_ptr<EC_Scalar_Data> scalar_invert(const EC_Scalar_Data& s) const;

      void scalar_assign(EC_Scalar_Data& x, const EC_Scalar_Data& y);

      void scalar_square_self(EC_Scalar_Data& s);

      std::unique_ptr<EC_Scalar_Data> scalar_negate(const EC_Scalar_Data& s) const;

      std::unique_ptr<EC_Scalar_Data> scalar_add(const EC_Scalar_Data& a, const EC_Scalar_Data& b) const;

      std::unique_ptr<EC_Scalar_Data> scalar_sub(const EC_Scalar_Data& a, const EC_Scalar_Data& b) const;

      std::unique_ptr<EC_Scalar_Data> scalar_mul(const EC_Scalar_Data& a, const EC_Scalar_Data& b) const;

      std::unique_ptr<EC_Scalar_Data> scalar_from_bigint(const BigInt& bn) const;

      std::unique_ptr<EC_Scalar_Data> gk_x_mod_order(const EC_Scalar_Data& scalar,
                                                     RandomNumberGenerator& rng,
                                                     std::vector<BigInt>& ws) const;

      std::unique_ptr<EC_Scalar_Data> scalar_deserialize(std::span<const uint8_t> bytes);

      std::vector<uint8_t> scalar_serialize(const EC_Scalar_Data& s) const;

      std::vector<uint8_t> scalar_serialize_pair(const EC_Scalar_Data& r, const EC_Scalar_Data& s) const;

   private:
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
