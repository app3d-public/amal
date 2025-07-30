constexpr vec() = default;
constexpr vec(const vec &rhs) noexcept = default;
constexpr vec(vec &&) noexcept = default;
constexpr vec &operator=(vec &&) noexcept = default;

template <
    typename V,
    std::enable_if_t<std::is_same<V, typename simd_type::value_type>::value && details::is_simd_enabled_v<V>, int> = 0>
constexpr explicit vec(V simd) noexcept : s(simd)
{
}

template <length_t L, AMAL_CONSTRUCT_SIMD>
constexpr explicit vec(const AMAL_VEC(L, T, P) & v) : s(v.s)
{
}

[[nodiscard]] static constexpr length_t length() { return vec_dimension; }

[[nodiscard]] constexpr T &operator[](length_t i)
{
    assert(i < vec_dimension);
    return data[i];
}

[[nodiscard]] constexpr T const &operator[](length_t i) const
{
    assert(i < vec_dimension);
    return data[i];
}

constexpr vec &operator=(vec const &v) = default;

template <typename U, enum Pack Q>
constexpr vec &operator=(AMAL_VEC(vec_dimension, U, Q) const &v)
{
    x = static_cast<T>(v.x);
    y = static_cast<T>(v.y);
    if constexpr (vec_dimension > 2) data[2] = static_cast<T>(v.z);
    if constexpr (vec_dimension > 3) data[3] = static_cast<T>(v.w);
    return *this;
}