/*
* (C) 2024 Jack Lloyd
*
* Botan is released under the Simplified BSD License (see license.txt)
*/

#include <botan/internal/ec_inner_data.h>

namespace Botan {

EC_Mul2Table_Data_BigInt::EC_Mul2Table_Data_BigInt(const EC_AffinePoint& h) :
   m_group(h.m_group), m_tbl(h.m_group->base_point(), h.to_legacy_point()) {}

std::optional<EC_AffinePoint> EC_Mul2Table_Data_BigInt::mul2(const EC_Scalar& x, const EC_Scalar& y) const {
   auto pt = m_tbl.multi_exp(x.m_scalar->value(), y.m_scalar->value());

   if(pt.is_zero()) {
      return std::nullopt;
   }
   return EC_AffinePoint(m_group, std::move(pt));
}

std::optional<EC_Scalar> EC_Mul2Table_Data_BigInt::mul2_x_mod_order(const EC_Scalar& x, const EC_Scalar& y) const {
   auto pt = m_tbl.multi_exp(x.m_scalar->value(), y.m_scalar->value());

   if(pt.is_zero()) {
      return std::nullopt;
   }
   return EC_Scalar(m_group, m_group->scalar_from_bigint(m_group->mod_order(pt.get_affine_x())));
}

}
