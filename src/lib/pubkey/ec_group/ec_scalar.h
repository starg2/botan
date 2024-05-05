/*
* (C) 2024 Jack Lloyd
*
* Botan is released under the Simplified BSD License (see license.txt)
*/

#ifndef BOTAN_EC_SCALAR_H_
#define BOTAN_EC_SCALAR_H_

#include <botan/types.h>

#include <optional>
#include <span>
#include <vector>

namespace Botan {

class BigInt;
class RandomNumberGenerator;
class EC_Group;
class EC_Group_Data;
class EC_Scalar_Data;

/// Represents an integer modulo the group order of an elliptic curve
class BOTAN_UNSTABLE_API EC_Scalar final {
   public:
      /// With ECDSA style shift followed by reduction mod order
      static EC_Scalar from_bytes_with_trunc(const EC_Group& group, std::span<const uint8_t> bytes);

      /// Reduces input modulo the group order
      static EC_Scalar from_bytes_mod_order(const EC_Group& group, std::span<const uint8_t> bytes);

      /// Returns nullopt if wrong size or if either scalar out of range
      static std::optional<std::pair<EC_Scalar, EC_Scalar>> deserialize_pair(const EC_Group& group,
                                                                             std::span<const uint8_t> bytes);

      static std::optional<EC_Scalar> deserialize(const EC_Group& group, std::span<const uint8_t> bytes);

      static EC_Scalar random(const EC_Group& group, RandomNumberGenerator& rng);

      static EC_Scalar one(const EC_Group& group);

      static EC_Scalar from_bigint(const EC_Group& group, const BigInt& bn);

      // Workspace argument is transitional
      static EC_Scalar gk_x_mod_order(const EC_Scalar& scalar, RandomNumberGenerator& rng, std::vector<BigInt>& ws);

      static std::vector<uint8_t> serialize_pair(const EC_Scalar& r, const EC_Scalar& s);

      std::vector<uint8_t> serialize() const;

      bool is_zero() const;

      bool is_nonzero() const { return !is_zero(); }

      EC_Scalar invert() const;

      EC_Scalar negate() const;

      EC_Scalar add(const EC_Scalar& x) const;

      EC_Scalar sub(const EC_Scalar& x) const;

      EC_Scalar mul(const EC_Scalar& x) const;

      void assign(const EC_Scalar& x);

      void square_self();

      bool is_eq(const EC_Scalar& x) const;

      friend EC_Scalar operator+(const EC_Scalar& x, const EC_Scalar& y) { return x.add(y); }

      friend EC_Scalar operator-(const EC_Scalar& x, const EC_Scalar& y) { return x.sub(y); }

      friend EC_Scalar operator*(const EC_Scalar& x, const EC_Scalar& y) { return x.mul(y); }

      friend bool operator==(const EC_Scalar& x, const EC_Scalar& y) { return x.is_eq(y); }

      EC_Scalar(const EC_Scalar& other);
      EC_Scalar(EC_Scalar&& other) noexcept;

      EC_Scalar& operator=(const EC_Scalar& x) {
         this->assign(x);
         return (*this);
      }

      ~EC_Scalar();

   private:
      friend class EC_AffinePoint;
      friend class EC_Mul2Table_Data;

      EC_Scalar(std::shared_ptr<EC_Group_Data> group, std::unique_ptr<EC_Scalar_Data> scalar);

      EC_Scalar(const EC_Group& group, std::unique_ptr<EC_Scalar_Data> scalar);

      const EC_Scalar_Data& inner() const { return *m_scalar; }

      const std::shared_ptr<EC_Group_Data>& group() const { return m_group; }

      std::shared_ptr<EC_Group_Data> m_group;
      std::unique_ptr<EC_Scalar_Data> m_scalar;
};

}  // namespace Botan

#endif
