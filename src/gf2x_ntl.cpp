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

inline std::uint64_t field_multiply_wide(
    std::uint64_t left,
    std::uint64_t right,
    std::uint32_t p,
    std::uint32_t ia,
    std::uint32_t ib) {
    const std::uint64_t a0 = left % p;
    const std::uint64_t a1 = left / p;
    const std::uint64_t b0 = right % p;
    const std::uint64_t b1 = right / p;
    const std::uint64_t c0 = (a0 * b0 + a1 * b1 * ib) % p;
    const std::uint64_t c1 =
        (a0 * b1 + a1 * b0 + a1 * b1 * ia) % p;
    return c0 + static_cast<std::uint64_t>(p) * c1;
}

inline std::uint64_t field_power_wide(
    std::uint64_t base,
    std::uint64_t exponent,
    std::uint32_t p,
    std::uint32_t ia,
    std::uint32_t ib) {
    std::uint64_t result = 1;
    while (exponent) {
        if (exponent & 1U) {
            result = field_multiply_wide(result, base, p, ia, ib);
        }
        base = field_multiply_wide(base, base, p, ia, ib);
        exponent >>= 1U;
    }
    return result;
}

inline std::uint64_t field_inverse_wide(
    std::uint64_t value,
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
    const std::uint64_t inverse_a = static_cast<std::uint64_t>(
        (((a + b * ia) % p) * inverse_norm) % p);
    const std::uint64_t inverse_b = static_cast<std::uint64_t>(
        (((p - b) % p) * inverse_norm) % p);
    return inverse_a + static_cast<std::uint64_t>(p) * inverse_b;
}

inline void prepare_prime_field_tables(
    std::uint32_t p,
    std::vector<std::uint32_t>& inverse,
    std::vector<std::uint8_t>& square) {
    inverse.assign(p, 0U);
    square.assign(p, 0U);
    inverse[1] = 1U;
    for (std::uint32_t value = 2; value < p; ++value) {
        inverse[value] = static_cast<std::uint32_t>(
            p - (static_cast<std::uint64_t>(p / value) *
                 inverse[p % value]) % p);
    }
    for (std::uint32_t value = 1; value <= (p - 1U) / 2U; ++value) {
        square[(static_cast<std::uint64_t>(value) * value) % p] = 1U;
    }
}

inline std::uint32_t field_inverse_from_table(
    std::uint32_t value,
    std::uint32_t p,
    std::uint32_t ia,
    std::uint32_t ib,
    const std::vector<std::uint32_t>& inverse) {
    const std::uint64_t a = value % p;
    const std::uint64_t b = value / p;
    const std::uint32_t norm = static_cast<std::uint32_t>(
        (a * a + a * b * ia + p - (b * b * ib) % p) % p);
    if (!norm) {
        return 0;
    }
    const std::uint64_t inverse_norm = inverse[norm];
    const std::uint32_t inverse_a = static_cast<std::uint32_t>(
        ((a + b * ia) % p) * inverse_norm % p);
    const std::uint32_t inverse_b = static_cast<std::uint32_t>(
        ((p - b) % p) * inverse_norm % p);
    return inverse_a + p * inverse_b;
}

inline std::uint64_t field_inverse_wide_from_table(
    std::uint64_t value,
    std::uint32_t p,
    std::uint32_t ia,
    std::uint32_t ib,
    const std::vector<std::uint32_t>& inverse) {
    const std::uint64_t a = value % p;
    const std::uint64_t b = value / p;
    const std::uint32_t norm = static_cast<std::uint32_t>(
        (a * a + a * b * ia + p - (b * b * ib) % p) % p);
    if (!norm) {
        return 0;
    }
    const std::uint64_t inverse_norm = inverse[norm];
    const std::uint64_t inverse_a =
        ((a + b * ia) % p) * inverse_norm % p;
    const std::uint64_t inverse_b = (p - b) * inverse_norm % p;
    return inverse_a + static_cast<std::uint64_t>(p) * inverse_b;
}

inline std::uint32_t order21_character(
    std::uint32_t base,
    std::uint32_t p,
    std::uint32_t ia,
    std::uint32_t ib,
    const std::vector<std::uint32_t>& inverse) {
    const std::uint64_t a = base % p;
    const std::uint64_t b = base / p;
    const std::uint32_t norm = static_cast<std::uint32_t>(
        (a * a + a * b * ia + p - (b * b * ib) % p) % p);
    const std::uint32_t scalar =
        prime_field_power(norm, (p - 1U) / 14U, p);
    const std::uint64_t scalar2 =
        static_cast<std::uint64_t>(scalar) * scalar % p;
    const std::uint64_t scalar4 = scalar2 * scalar2 % p;
    const std::uint32_t scalar5 = static_cast<std::uint32_t>(
        scalar4 * scalar % p);

    const std::uint32_t inverse_base = field_inverse_from_table(
        base, p, ia, ib, inverse);
    const std::uint32_t conjugate = static_cast<std::uint32_t>(
        (a + b * ia) % p) + p * static_cast<std::uint32_t>(
            b == 0 ? 0 : p - b);
    const std::uint32_t torus = field_multiply(
        conjugate, inverse_base, p, ia, ib);
    const std::uint32_t projective = field_power(
        torus, (p + 1U) / 6U, p, ia, ib);
    return static_cast<std::uint32_t>(
        (static_cast<std::uint64_t>(projective % p) * scalar5) % p) +
        p * static_cast<std::uint32_t>(
            (static_cast<std::uint64_t>(projective / p) * scalar5) % p);
}

inline std::uint64_t order21_character_wide(
    std::uint64_t base,
    std::uint32_t p,
    std::uint32_t ia,
    std::uint32_t ib,
    const std::vector<std::uint32_t>& inverse) {
    const std::uint64_t a = base % p;
    const std::uint64_t b = base / p;
    const std::uint32_t norm = static_cast<std::uint32_t>(
        (a * a + a * b * ia + p - (b * b * ib) % p) % p);
    const std::uint32_t scalar =
        prime_field_power(norm, (p - 1U) / 14U, p);
    const std::uint64_t scalar2 =
        static_cast<std::uint64_t>(scalar) * scalar % p;
    const std::uint64_t scalar4 = scalar2 * scalar2 % p;
    const std::uint64_t scalar5 = scalar4 * scalar % p;

    const std::uint64_t inverse_base = field_inverse_wide_from_table(
        base, p, ia, ib, inverse);
    const std::uint64_t conjugate = (a + b * ia) % p +
        static_cast<std::uint64_t>(p) * (b == 0 ? 0 : p - b);
    const std::uint64_t torus = field_multiply_wide(
        conjugate, inverse_base, p, ia, ib);
    const std::uint64_t projective = field_power_wide(
        torus, (p + 1U) / 6U, p, ia, ib);
    return (projective % p) * scalar5 % p +
        static_cast<std::uint64_t>(p) *
            ((projective / p) * scalar5 % p);
}

}  // namespace

extern "C" std::uint64_t qml_field_primitive(
    std::uint32_t p,
    std::uint32_t ia,
    std::uint32_t ib) {
    if (p < 3 || p > 1000000U) {
        return 0;
    }
    const std::uint64_t order = static_cast<std::uint64_t>(p) * p - 1U;
    std::uint64_t remaining = order;
    std::vector<std::uint64_t> prime_factors;
    for (std::uint64_t divisor = 2; divisor * divisor <= remaining;
         divisor += divisor == 2 ? 1 : 2) {
        if (remaining % divisor != 0) {
            continue;
        }
        prime_factors.push_back(divisor);
        do {
            remaining /= divisor;
        } while (remaining % divisor == 0);
    }
    if (remaining > 1) {
        prime_factors.push_back(remaining);
    }
    const std::uint64_t q = static_cast<std::uint64_t>(p) * p;
    for (std::uint64_t value = 2; value < q; ++value) {
        bool primitive = true;
        for (const std::uint64_t factor : prime_factors) {
            if (field_power_wide(
                    value, order / factor, p, ia, ib) == 1) {
                primitive = false;
                break;
            }
        }
        if (primitive) {
            return value;
        }
    }
    return 0;
}

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

extern "C" long qml_gf2x_cyclic_product(
    const std::uint8_t* a,
    long a_bytes,
    const std::uint8_t* b,
    long b_bytes,
    long order,
    std::uint8_t* output,
    long output_capacity) {
    if (!a || !b || !output || a_bytes < 0 || b_bytes < 0 || order < 1 ||
        output_capacity < (order + 7) / 8) {
        return -1;
    }
    NTL::GF2X left;
    NTL::GF2X right;
    NTL::GF2X product;
    NTL::GF2X modulus;
    NTL::GF2X result;
    NTL::GF2XFromBytes(left, a, a_bytes);
    NTL::GF2XFromBytes(right, b, b_bytes);
    NTL::mul(product, left, right);
    NTL::SetCoeff(modulus, order);
    NTL::SetCoeff(modulus, 0);
    NTL::rem(result, product, modulus);
    const long result_bytes = (order + 7) / 8;
    std::memset(output, 0, static_cast<std::size_t>(output_capacity));
    NTL::BytesFromGF2X(output, result, result_bytes);
    return result_bytes;
}

extern "C" long qml_gf2x_cyclic_star_product(
    const std::uint8_t* a,
    long a_bytes,
    const std::uint8_t* b,
    long b_bytes,
    long order,
    std::uint8_t* output,
    long output_capacity) {
    if (!a || !b || !output || a_bytes < 0 || b_bytes < 0 || order < 1 ||
        output_capacity < (order + 7) / 8) {
        return -1;
    }
    NTL::GF2X left;
    NTL::GF2X right;
    NTL::GF2X right_star;
    NTL::GF2X product;
    NTL::GF2X modulus;
    NTL::GF2X result;
    NTL::GF2XFromBytes(left, a, a_bytes);
    NTL::GF2XFromBytes(right, b, b_bytes);
    for (long exponent = 0; exponent < order; ++exponent) {
        if (NTL::IsOne(NTL::coeff(right, exponent))) {
            NTL::SetCoeff(right_star, exponent == 0 ? 0 : order - exponent);
        }
    }
    NTL::mul(product, left, right_star);
    NTL::SetCoeff(modulus, order);
    NTL::SetCoeff(modulus, 0);
    NTL::rem(result, product, modulus);
    const long result_bytes = (order + 7) / 8;
    std::memset(output, 0, static_cast<std::size_t>(output_capacity));
    NTL::BytesFromGF2X(output, result, result_bytes);
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
    // These buffers persist within each worker process.  One O(p) setup
    // replaces norm inversion and Euler-criterion exponentiations on every
    // point of every selected line.
    static thread_local std::vector<std::uint32_t> prime_inverse;
    static thread_local std::vector<std::uint8_t> quadratic_residue;
    prepare_prime_field_tables(p, prime_inverse, quadratic_residue);
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
            const std::uint32_t denominator_inverse = field_inverse_from_table(
                denominator, p, ia, ib, prime_inverse);
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
            const bool square = quadratic_residue[norm] != 0;
            const std::uint32_t base = square ? point : field_multiply(
                omega_inverse, point, p, ia, ib);
            const std::uint64_t component = square ? 0 : 1;

            for (long order_index = 0; order_index < n_orders; ++order_index) {
                const std::uint32_t order = orders[order_index];
                const std::uint32_t value =
                    order == 21U && p % 84U == 29U
                    ? order21_character(
                          base, p, ia, ib, prime_inverse)
                    : field_power(
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

extern "C" int qml_selected_line_bins_wide(
    std::uint32_t p,
    std::uint32_t ia,
    std::uint32_t ib,
    std::uint64_t generator,
    std::uint64_t omega,
    std::uint64_t sigma,
    std::uint64_t pole_t,
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
    // The wide formulas use products bounded by roughly p^3.  This guard is
    // comfortably below the uint64 overflow threshold and above our scans.
    if (p > 1000000U ||
        needed > static_cast<std::uint64_t>(output_capacity)) {
        return -2;
    }
    std::memset(output, 0, static_cast<std::size_t>(output_capacity));

    const std::uint64_t omega_inverse =
        field_inverse_wide(omega, p, ia, ib);
    if (!omega_inverse) {
        return -3;
    }
    static thread_local std::vector<std::uint32_t> prime_inverse;
    static thread_local std::vector<std::uint8_t> quadratic_residue;
    prepare_prime_field_tables(p, prime_inverse, quadratic_residue);
    std::vector<std::vector<std::pair<std::uint64_t, std::uint32_t>>> tables;
    tables.reserve(static_cast<std::size_t>(n_orders));
    for (long index = 0; index < n_orders; ++index) {
        const std::uint32_t order = orders[index];
        if (order < 1 || orbit_length % order != 0 ||
            order_offsets[index] + order >
                static_cast<std::uint32_t>(total_order)) {
            return -4;
        }
        const std::uint64_t root = field_power_wide(
            generator, orbit_length / order, p, ia, ib);
        std::vector<std::pair<std::uint64_t, std::uint32_t>> table;
        table.reserve(order);
        std::uint64_t point = 1;
        for (std::uint32_t residue = 0; residue < order; ++residue) {
            table.emplace_back(point, residue);
            point = field_multiply_wide(point, root, p, ia, ib);
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
            const std::uint64_t line_coordinate =
                z + static_cast<std::uint64_t>(p) * level;
            const std::uint64_t image = field_multiply_wide(
                sigma, line_coordinate, p, ia, ib);
            std::uint64_t denominator = field_multiply_wide(
                pole_t, image, p, ia, ib);
            denominator = (denominator % p + p - 1) % p +
                static_cast<std::uint64_t>(p) * (denominator / p);
            const std::uint64_t denominator_inverse =
                field_inverse_wide_from_table(
                    denominator, p, ia, ib, prime_inverse);
            if (!denominator_inverse) {
                continue;
            }
            const std::uint64_t point = field_multiply_wide(
                image, denominator_inverse, p, ia, ib);
            if (!point) {
                continue;
            }
            const std::uint64_t a = point % p;
            const std::uint64_t b = point / p;
            const std::uint32_t norm = static_cast<std::uint32_t>(
                (a * a + a * b * ia + p - (b * b * ib) % p) % p);
            const bool square = quadratic_residue[norm] != 0;
            const std::uint64_t base = square ? point : field_multiply_wide(
                omega_inverse, point, p, ia, ib);
            const std::uint64_t component = square ? 0 : 1;

            for (long order_index = 0; order_index < n_orders; ++order_index) {
                const std::uint32_t order = orders[order_index];
                const std::uint64_t value =
                    order == 21U && p % 84U == 29U
                    ? order21_character_wide(
                          base, p, ia, ib, prime_inverse)
                    : field_power_wide(
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
