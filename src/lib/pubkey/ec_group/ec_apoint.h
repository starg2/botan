/*
* (C) 2024 Jack Lloyd
*
* Botan is released under the Simplified BSD License (see license.txt)
*/

#ifndef BOTAN_EC_APOINT_H_
#define BOTAN_EC_APOINT_H_

#include <botan/secmem.h>
#include <botan/types.h>
#include <optional>
#include <span>
#include <string_view>
#include <vector>

namespace Botan {

class BigInt;
class RandomNumberGenerator;
class EC_Group;
class EC_Scalar;
class EC_Point;

class EC_Group_Data;
class EC_Point_Data;

class BOTAN_UNSTABLE_API EC_AffinePoint final {
   public:
      /// Point deserialization. Returns nullopt if wrong length or not a valid point
      static std::optional<EC_AffinePoint> deserialize(const EC_Group& group, std::span<const uint8_t> bytes);

      /// Multiply by the group generator returning a complete point
      ///
      /// Workspace argument is transitional
      static EC_AffinePoint g_mul(const EC_Scalar& scalar, RandomNumberGenerator& rng, std::vector<BigInt>& ws);

      /// Hash to curve, random oracle variant
      ///
      /// Only supported for specific groups
      static EC_AffinePoint hash_to_curve_ro(const EC_Group& group,
                                             std::string_view hash_fn,
                                             std::span<const uint8_t> input,
                                             std::span<const uint8_t> domain_sep);

      /// Hash to curve, non uniform variant
      static EC_AffinePoint hash_to_curve_nu(const EC_Group& group,
                                             std::string_view hash_fn,
                                             std::span<const uint8_t> input,
                                             std::span<const uint8_t> domain_sep);

      /// Multiply a point by a scalar returning a complete point
      ///
      /// Workspace argument is transitional
      EC_AffinePoint mul(const EC_Scalar& scalar, RandomNumberGenerator& rng, std::vector<BigInt>& ws) const;

      /// Return the fixed length encoding of affine x coordinate
      secure_vector<uint8_t> x_bytes() const;

      /// Return the fixed length encoding of affine y coordinate
      secure_vector<uint8_t> y_bytes() const;

      /// Return the fixed length encoding of affine x and y coordinates
      secure_vector<uint8_t> xy_bytes() const;

      /// Return the fixed length encoding of SEC1 compressed encoding
      std::vector<uint8_t> serialize_compressed();

      /// Return the fixed length encoding of SEC1 uncompressed encoding
      std::vector<uint8_t> serialize_uncompressed();

      EC_AffinePoint(const EC_AffinePoint& other);
      EC_AffinePoint(EC_AffinePoint&& other) noexcept;

      EC_AffinePoint(const EC_Group& group, const EC_Point& pt);
      EC_AffinePoint(const std::shared_ptr<EC_Group_Data>, EC_Point&& pt);

      EC_Point to_legacy_point() const;

      ~EC_AffinePoint();

   private:
      friend class EC_Mul2Table_Data;

      EC_AffinePoint(std::shared_ptr<EC_Group_Data> group, std::unique_ptr<EC_Point_Data> point);

      std::shared_ptr<EC_Group_Data> m_group;
      std::unique_ptr<EC_Point_Data> m_point;
};

}  // namespace Botan

#endif
