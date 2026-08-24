#include <NTL/GF2X.h>

#include <cstddef>
#include <cstdint>
#include <cstring>

namespace {

inline std::uint32_t field_multiply(
    std::uint32_t left,
    std::uint32_t right,
    std::uint32_t p,
    std::uint32_t ia,
    std::uint32_t ib) {
    const std::uint64_t a0 = left % p;
    const std::uint64_t a1 = left / p;
    const std::uint64_t b0 = right % p;
    const std::uint64_t b1 = right / p;
    const std::uint32_t c0 = static_cast<std::uint32_t>(
        (a0 * b0 + a1 * b1 * ib) % p);
    const std::uint32_t c1 = static_cast<std::uint32_t>(
        (a0 * b1 + a1 * b0 + a1 * b1 * ia) % p);
    return c0 + p * c1;
}

}  // namespace

extern "C" long qml_gf2x_gcd(
    const std::uint8_t* a,
    long a_bytes,
    const std::uint8_t* b,
    long b_bytes,
    std::uint8_t* output,
    long output_capacity) {
    if (!a || !b || !output || a_bytes < 0 || b_bytes < 0 || output_capacity < 1) {
        return -1;
    }
    NTL::GF2X left;
    NTL::GF2X right;
    NTL::GF2X result;
    NTL::GF2XFromBytes(left, a, a_bytes);
    NTL::GF2XFromBytes(right, b, b_bytes);
    NTL::GCD(result, left, right);
    const long degree = NTL::deg(result);
    const long result_bytes = degree < 0 ? 1 : (degree + 8) / 8;
    if (result_bytes > output_capacity) {
        return -result_bytes;
    }
    std::memset(output, 0, static_cast<std::size_t>(output_capacity));
    if (degree >= 0) {
        NTL::BytesFromGF2X(output, result, result_bytes);
    }
    return result_bytes;
}

extern "C" int qml_field_orbits(
    std::uint32_t p,
    std::uint32_t ia,
    std::uint32_t ib,
    std::uint32_t generator,
    std::uint32_t omega,
    std::uint32_t* square,
    std::uint32_t* nonsquare,
    std::uint32_t* logarithm,
    long orbit_length) {
    if (!square || !nonsquare || p < 3 || orbit_length < 1 ||
        orbit_length >= (1L << 31)) {
        return -1;
    }
    const std::uint64_t q = static_cast<std::uint64_t>(p) * p;
    if (logarithm) {
        std::memset(
            logarithm,
            0,
            static_cast<std::size_t>(q) * sizeof(std::uint32_t));
    }
    std::uint32_t point = 1;
    for (long j = 0; j < orbit_length; ++j) {
        const std::uint32_t other = field_multiply(omega, point, p, ia, ib);
        square[j] = point;
        nonsquare[j] = other;
        if (logarithm) {
            logarithm[point] = static_cast<std::uint32_t>(j);
            logarithm[other] = static_cast<std::uint32_t>(j) | 0x80000000U;
        }
        point = field_multiply(point, generator, p, ia, ib);
    }
    return point == 1 ? 0 : -2;
}
