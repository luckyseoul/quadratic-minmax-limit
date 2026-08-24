#include <NTL/GF2X.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <utility>
#include <vector>

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

inline std::uint32_t field_power(
    std::uint32_t base,
    std::uint64_t exponent,
    std::uint32_t p,
    std::uint32_t ia,
    std::uint32_t ib) {
    std::uint32_t result = 1;
    while (exponent) {
        if (exponent & 1U) {
            result = field_multiply(result, base, p, ia, ib);
        }
        base = field_multiply(base, base, p, ia, ib);
        exponent >>= 1U;
    }
    return result;
}

inline std::uint32_t prime_field_power(
    std::uint32_t base,
    std::uint32_t exponent,
    std::uint32_t p) {
    std::uint64_t result = 1;
    std::uint64_t current = base;
    while (exponent) {
        if (exponent & 1U) {
            result = result * current % p;
        }
        current = current * current % p;
        exponent >>= 1U;
    }
    return static_cast<std::uint32_t>(result);
}

inline std::uint32_t field_inverse(
    std::uint32_t value,
    std::uint32_t p,
    std::uint32_t ia,
    std::uint32_t ib) {
    const std::uint64_t a = value % p;
    const std::uint64_t b = value / p;
    const std::uint32_t norm = static_cast<std::uint32_t>(
        (a * a + a * b * ia + p - (b * b * ib) % p) % p);
    if (!norm) {
        return 0;
    }
    const std::uint64_t inverse_norm = prime_field_power(norm, p - 2, p);
    const std::uint32_t inverse_a = static_cast<std::uint32_t>(
        ((a + b * ia) % p) * inverse_norm % p);
    const std::uint32_t inverse_b = static_cast<std::uint32_t>(
        ((p - b) % p) * inverse_norm % p);
    return inverse_a + p * inverse_b;
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

extern "C" int qml_selected_line_bins(
    std::uint32_t p,
    std::uint32_t ia,
    std::uint32_t ib,
    std::uint32_t generator,
    std::uint32_t omega,
    std::uint32_t sigma,
    std::uint32_t pole_t,
    const std::uint32_t* levels,
    long n_levels,
    const std::uint32_t* orders,
    const std::uint32_t* order_offsets,
    long n_orders,
    long total_order,
    std::uint8_t* output,
    long output_capacity) {
    if (p < 3 || !levels || n_levels < 1 || !orders || !order_offsets ||
        n_orders < 1 || total_order < 1 || !output || output_capacity < 1) {
        return -1;
    }
    const std::uint64_t q = static_cast<std::uint64_t>(p) * p;
    const std::uint64_t orbit_length = (q - 1) / 2;
    const std::uint64_t needed =
        2ULL * static_cast<std::uint64_t>(n_levels) * total_order;
    if (q >= (1ULL << 32) || needed > static_cast<std::uint64_t>(output_capacity)) {
        return -2;
    }
    std::memset(output, 0, static_cast<std::size_t>(output_capacity));

    const std::uint32_t omega_inverse = field_inverse(omega, p, ia, ib);
    if (!omega_inverse) {
        return -3;
    }
    std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>> tables;
    tables.reserve(static_cast<std::size_t>(n_orders));
    for (long index = 0; index < n_orders; ++index) {
        const std::uint32_t order = orders[index];
        if (order < 1 || orbit_length % order != 0 ||
            order_offsets[index] + order > static_cast<std::uint32_t>(total_order)) {
            return -4;
        }
        const std::uint32_t root = field_power(
            generator, orbit_length / order, p, ia, ib);
        std::vector<std::pair<std::uint32_t, std::uint32_t>> table;
        table.reserve(order);
        std::uint32_t point = 1;
        for (std::uint32_t residue = 0; residue < order; ++residue) {
            table.emplace_back(point, residue);
            point = field_multiply(point, root, p, ia, ib);
        }
        if (point != 1) {
            return -5;
        }
        std::sort(table.begin(), table.end());
        tables.push_back(std::move(table));
    }

    for (long level_index = 0; level_index < n_levels; ++level_index) {
        const std::uint32_t level = levels[level_index];
        if (level >= p) {
            return -6;
        }
        for (std::uint32_t z = 0; z < p; ++z) {
            const std::uint32_t line_coordinate = z + p * level;
            const std::uint32_t image = field_multiply(
                sigma, line_coordinate, p, ia, ib);
            std::uint32_t denominator = field_multiply(
                pole_t, image, p, ia, ib);
            denominator = (denominator % p + p - 1) % p +
                p * (denominator / p);
            const std::uint32_t denominator_inverse = field_inverse(
                denominator, p, ia, ib);
            if (!denominator_inverse) {
                continue;
            }
            const std::uint32_t point = field_multiply(
                image, denominator_inverse, p, ia, ib);
            if (!point) {
                continue;
            }
            const std::uint64_t a = point % p;
            const std::uint64_t b = point / p;
            const std::uint32_t norm = static_cast<std::uint32_t>(
                (a * a + a * b * ia + p - (b * b * ib) % p) % p);
            const bool square = prime_field_power(norm, (p - 1) / 2, p) == 1;
            const std::uint32_t base = square ? point : field_multiply(
                omega_inverse, point, p, ia, ib);
            const std::uint64_t component = square ? 0 : 1;

            for (long order_index = 0; order_index < n_orders; ++order_index) {
                const std::uint32_t order = orders[order_index];
                const std::uint32_t value = field_power(
                    base, orbit_length / order, p, ia, ib);
                const auto& table = tables[order_index];
                const auto found = std::lower_bound(
                    table.begin(), table.end(),
                    std::make_pair(value, static_cast<std::uint32_t>(0)));
                if (found == table.end() || found->first != value) {
                    return -7;
                }
                const std::uint64_t output_index =
                    (component * n_levels + level_index) * total_order +
                    order_offsets[order_index] + found->second;
                output[output_index] ^= 1U;
            }
        }
    }
    return 0;
}
