// Exact finite scout for the p+3 one-profile energy, two-double case.
//
// Affine-normalize the three omitted roots to 0,1,rho.  For each pair of
// repeated roots alpha,beta, form
//
//   P=(x^p-x)(x-alpha)(x-beta)/(x(x-1)(x-rho)).
//
// The profile exists iff P+d is a monic square T^2 for a nonzero square d.
// High coefficients determine T uniquely, so each tuple costs O(p^2) exact
// field operations.  Absence is a finite audit, not an all-prime proof.
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

static int mod(long long value, int p) {
    value %= p;
    if (value < 0) value += p;
    return static_cast<int>(value);
}

static int mod_pow(int base, int exponent, int p) {
    long long out = 1;
    long long x = mod(base, p);
    while (exponent) {
        if (exponent & 1) out = out * x % p;
        x = x * x % p;
        exponent >>= 1;
    }
    return static_cast<int>(out);
}

static int eval_poly(const std::vector<int>& coefficients, int x, int p) {
    long long out = 0;
    for (auto it = coefficients.rbegin(); it != coefficients.rend(); ++it) {
        out = (out * x + *it) % p;
    }
    return static_cast<int>(out);
}

struct Example {
    int rho;
    int alpha;
    int beta;
    int d;
    int root_d;
    std::vector<int> profile;
};

static std::vector<int> quotient_for_rho(int p, int rho) {
    // (x^(p-1)-1)/((x-1)(x-rho)), low coefficients first.
    std::vector<int> remainder(p, 0);
    remainder[0] = p - 1;
    remainder[p - 1] = 1;
    std::vector<int> quotient(p - 2, 0);
    const int den1 = mod(-(1 + rho), p);
    const int den0 = rho;
    for (int degree = p - 1; degree >= 2; --degree) {
        int q = remainder[degree];
        quotient[degree - 2] = q;
        remainder[degree] = 0;
        remainder[degree - 1] = mod(remainder[degree - 1] - 1LL * q * den1, p);
        remainder[degree - 2] = mod(remainder[degree - 2] - 1LL * q * den0, p);
    }
    if (remainder[0] != 0 || remainder[1] != 0) {
        throw std::runtime_error("nonzero polynomial-division remainder");
    }
    return quotient;
}

static std::vector<int> product_polynomial(
    const std::vector<int>& quotient, int alpha, int beta, int p
) {
    std::vector<int> out(p, 0);
    const int linear = mod(-(alpha + beta), p);
    const int constant = mod(1LL * alpha * beta, p);
    for (int i = 0; i < static_cast<int>(quotient.size()); ++i) {
        out[i] = mod(out[i] + 1LL * quotient[i] * constant, p);
        out[i + 1] = mod(out[i + 1] + 1LL * quotient[i] * linear, p);
        out[i + 2] = mod(out[i + 2] + quotient[i], p);
    }
    return out;
}

static bool square_shift(
    const std::vector<int>& polynomial,
    int p,
    std::vector<int>& root,
    int& shift
) {
    const int m = (p - 1) / 2;
    const int inv2 = (p + 1) / 2;
    root.assign(m + 1, 0);
    root[m] = 1;
    for (int degree = 2 * m - 1; degree >= m; --degree) {
        const int unknown = degree - m;
        long long known = 0;
        for (int i = unknown + 1; i < m; ++i) {
            int j = degree - i;
            if (j > unknown && j < m) known += 1LL * root[i] * root[j];
        }
        root[unknown] = mod(1LL * (polynomial[degree] - mod(known, p)) * inv2, p);
    }
    for (int degree = 1; degree <= 2 * m; ++degree) {
        long long coefficient = 0;
        int lo = std::max(0, degree - m);
        int hi = std::min(m, degree);
        for (int i = lo; i <= hi; ++i) coefficient += 1LL * root[i] * root[degree - i];
        if (mod(coefficient, p) != polynomial[degree]) return false;
    }
    shift = mod(1LL * root[0] * root[0] - polynomial[0], p);
    return shift != 0 && mod_pow(shift, m, p) == 1;
}

static std::string vector_json(const std::vector<int>& values) {
    std::ostringstream out;
    out << "[";
    for (size_t i = 0; i < values.size(); ++i) {
        if (i) out << ",";
        out << values[i];
    }
    out << "]";
    return out.str();
}

static void scan_prime(int p, int stop_after) {
    std::uint64_t checked = 0;
    std::vector<Example> examples;
    auto started = std::chrono::steady_clock::now();
    for (int rho = 2; rho < p; ++rho) {
        auto quotient = quotient_for_rho(p, rho);
        std::vector<int> available;
        for (int x = 0; x < p; ++x) {
            if (x != 0 && x != 1 && x != rho) available.push_back(x);
        }
        for (size_t ia = 0; ia < available.size(); ++ia) {
            for (size_t ib = ia + 1; ib < available.size(); ++ib) {
                ++checked;
                int alpha = available[ia], beta = available[ib];
                auto polynomial = product_polynomial(quotient, alpha, beta, p);
                std::vector<int> T;
                int d = 0;
                if (!square_shift(polynomial, p, T, d)) continue;
                int root_d = 1;
                while (root_d < p && mod(1LL * root_d * root_d, p) != d) ++root_d;
                if (root_d == p) throw std::runtime_error("square shift has no root");
                std::vector<int> profile(p, 0);
                for (int x = 0; x < p; ++x) {
                    int value = eval_poly(T, x, p);
                    if (value == mod(-root_d, p)) {
                        profile[x] = (x == alpha || x == beta) ? 2 : 1;
                    } else if (value == root_d) {
                        profile[x] = (x == alpha || x == beta) ? -2 : -1;
                    }
                }
                examples.push_back({rho, alpha, beta, d, root_d, profile});
                if (static_cast<int>(examples.size()) >= stop_after) goto done;
            }
        }
    }
done:
    double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started
    ).count();
    std::cout << "{\"p\":" << p
              << ",\"checked\":" << checked
              << ",\"complete\":" << (examples.size() < static_cast<size_t>(stop_after) ? "true" : "false")
              << ",\"seconds\":" << seconds
              << ",\"examples\":[";
    for (size_t i = 0; i < examples.size(); ++i) {
        if (i) std::cout << ",";
        const auto& e = examples[i];
        std::cout << "{\"rho\":" << e.rho
                  << ",\"alpha\":" << e.alpha
                  << ",\"beta\":" << e.beta
                  << ",\"d\":" << e.d
                  << ",\"root_d\":" << e.root_d
                  << ",\"profile\":" << vector_json(e.profile) << "}";
    }
    std::cout << "]}" << std::endl;
}

int main(int argc, char** argv) {
    std::string primes = argc > 1 ? argv[1] : "11,13,17,19,23,29,31";
    int stop_after = argc > 2 ? std::stoi(argv[2]) : 10;
    std::stringstream parser(primes);
    std::string item;
    while (std::getline(parser, item, ',')) scan_prime(std::stoi(item), stop_after);
    return 0;
}
