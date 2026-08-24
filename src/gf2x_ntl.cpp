#include <NTL/GF2X.h>

#include <cstddef>
#include <cstdint>
#include <cstring>

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
