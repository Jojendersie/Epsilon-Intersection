namespace ei { namespace details {
    /// \brief Template construct to create integer types by a size value.
    /// \details The specializations define signed and unsigned types derived
    ///     from a standard type. Using Int<2>::stype... will chose a
    ///     specialization with the required size automatically.
    ///
    ///     The DUMMY parameter avoids the creation of two equal
    ///     specializations in case one or multiple types have the same size.
    template <int BYTES, int DUMMY = 0> struct Int
    {
    };

    template <> struct Int<sizeof(char)>
    {
        typedef signed char stype;
        typedef unsigned char utype;
    };

    template <> struct Int<sizeof(short), sizeof(short) == sizeof(char) ? 1 : 0>
    {
        typedef signed short stype;
        typedef unsigned short utype;
    };

    template <> struct Int<sizeof(int), sizeof(int) == sizeof(short) ? 2 : 0>
    {
        typedef signed int stype;
        typedef unsigned int utype;
    };

    template <> struct Int<sizeof(long), sizeof(long) == sizeof(int) ? 3 : 0>
    {
        typedef signed long stype;
        typedef unsigned long utype;
    };

    template <> struct Int<sizeof(long long), sizeof(long long) == sizeof(long) ? 4 : 0>
    {
        typedef signed long long stype;
        typedef unsigned long long utype;
    };

} // namespace details

    template<int BYTES>
    using Sint = typename details::Int<BYTES>::stype;
    template<int BYTES>
    using Uint = typename details::Int<BYTES>::utype;

} // namespace ei
