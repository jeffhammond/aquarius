#ifndef _MARRAY_MITERATOR_HPP_
#define _MARRAY_MITERATOR_HPP_

#include <cassert>
#include <algorithm>
#include <vector>
#include <array>
#include <cstddef>
#include <type_traits>
#include <utility>

namespace MArray
{

namespace detail
{
    template <typename...>
    struct is_container_helper {};

    template <typename T, typename=void>
    struct is_container : std::false_type {};

    template <typename T>
    struct is_container<T,
        typename std::conditional<false,
                                  is_container_helper<typename T::value_type,
                                                      decltype(std::declval<T>().size()),
                                                      decltype(std::declval<T>().begin()),
                                                      decltype(std::declval<T>().end())>,
                                  void>::type>
    : std::true_type {};

    template <typename T, typename C, typename=void>
    struct is_container_of : std::false_type {};

    template <typename T, typename C>
    struct is_container_of<T, C, typename std::enable_if<is_container<C>::value>::type>
    : std::is_convertible<typename C::value_type, T> {};

    template <typename T, typename... Ts>
    struct are_containers_helper;

    template <typename T>
    struct are_containers_helper<T> : is_container<T> {};

    template <typename T, typename... Ts>
    struct are_containers_helper
    : std::conditional<is_container<T>::value,
                       are_containers_helper<Ts...>,
                       std::false_type>::type {};

    template <typename... Ts>
    struct are_containers;

    template <>
    struct are_containers<> : std::true_type {};

    template <typename... Ts>
    struct are_containers : are_containers_helper<Ts...> {};

    template <typename T, typename C, typename... Cs>
    struct are_containers_of_helper;

    template <typename T, typename C>
    struct are_containers_of_helper<T, C> : is_container_of<T, C> {};

    template <typename T, typename C, typename... Cs>
    struct are_containers_of_helper
    : std::conditional<is_container_of<T, C>::value,
                       are_containers_of_helper<T, Cs...>,
                       std::false_type>::type {};

    template <typename T, typename... Cs>
    struct are_containers_of;

    template <typename T>
    struct are_containers_of<T> : std::true_type {};

    template <typename T, typename... Cs>
    struct are_containers_of : are_containers_of_helper<T, Cs...> {};

    template <typename T, size_t N, typename C>
    struct is_array_of : std::false_type {};

    template <typename T, size_t N>
    struct is_array_of<T, N, std::array<T,N>> : std::true_type {};

    template <typename T, size_t N, typename C, typename... Cs>
    struct are_arrays_of_helper;

    template <typename T, size_t N, typename C>
    struct are_arrays_of_helper<T, N, C> : is_array_of<T, N, C> {};

    template <typename T, size_t N, typename C, typename... Cs>
    struct are_arrays_of_helper
    : std::conditional<is_array_of<T, N, C>::value,
                       are_arrays_of_helper<T, N, Cs...>,
                       std::false_type>::type {};

    template <typename T, size_t N, typename... Cs>
    struct are_arrays_of;

    template <typename T, size_t N>
    struct are_arrays_of<T, N> : std::true_type {};

    template <typename T, size_t N, typename... Cs>
    struct are_arrays_of : are_arrays_of_helper<T, N, Cs...> {};

    template <typename stride_type, size_t ndim, size_t N, size_t I, typename Offset, typename... Offsets>
    struct inc_offsets_helper;

    template <typename stride_type, size_t ndim, size_t N, typename Offset>
    struct inc_offsets_helper<stride_type, ndim, N, N, Offset>
    {
        inc_offsets_helper(int i,
                           Offset& off0,
                           const std::array<std::array<stride_type,ndim>,N>& strides)
        {
            off0 += strides[N-1][i];
        }
    };

    template <typename stride_type, size_t ndim, size_t N, size_t I, typename Offset, typename... Offsets>
    struct inc_offsets_helper
    {
        inc_offsets_helper(unsigned i,
                           Offset& off0, Offsets&... off,
                           const std::array<std::array<stride_type,ndim>,N>& strides)
        {
            off0 += strides[I-1][i];
            inc_offsets_helper<stride_type, ndim, N, I+1, Offsets...>(i, off..., strides);
        }
    };

    template <typename stride_type, size_t ndim, size_t N, typename... Offsets>
    void inc_offsets(unsigned i,
                     const std::array<std::array<stride_type,ndim>,N>& strides,
                     Offsets&... off)
    {
        inc_offsets_helper<stride_type, ndim, N, 1, Offsets...>(i, off..., strides);
    }

    template <typename idx_type, typename stride_type, size_t ndim, size_t N, size_t I, typename Offset, typename... Offsets>
    struct dec_offsets_helper;

    template <typename idx_type, typename stride_type, size_t ndim, size_t N, typename Offset>
    struct dec_offsets_helper<idx_type, stride_type, ndim, N, N, Offset>
    {
        dec_offsets_helper(unsigned i,
                           Offset& off0,
                           const std::array<idx_type,ndim>& pos,
                           const std::array<std::array<stride_type,ndim>,N>& strides)
        {
            off0 -= pos[i]*strides[N-1][i];
        }
    };

    template <typename idx_type, typename stride_type, size_t ndim, size_t N, size_t I, typename Offset, typename... Offsets>
    struct dec_offsets_helper
    {
        dec_offsets_helper(unsigned i,
                           Offset& off0, Offsets&... off,
                           const std::array<idx_type,ndim>& pos,
                           const std::array<std::array<stride_type,ndim>,N>& strides)
        {
            off0 -= pos[i]*strides[I-1][i];
            dec_offsets_helper<idx_type, stride_type, ndim, N, I+1, Offsets...>(i, off..., pos, strides);
        }
    };

    template <typename idx_type, typename stride_type, size_t ndim, size_t N, typename... Offsets>
    void dec_offsets(unsigned i,
                     const std::array<idx_type,ndim>& pos,
                     const std::array<std::array<stride_type,ndim>,N>& strides,
                     Offsets&... off)
    {
        dec_offsets_helper<idx_type, stride_type, ndim, N, 1, Offsets...>(i, off..., pos, strides);
    }

    template <typename idx_type, typename stride_type, size_t ndim, size_t N, size_t I, typename Offset, typename... Offsets>
    struct set_offsets_helper;

    template <typename idx_type, typename stride_type, size_t ndim, size_t N, typename Offset>
    struct set_offsets_helper<idx_type, stride_type, ndim, N, N, Offset>
    {
        set_offsets_helper(Offset& off0,
                           const std::array<idx_type,ndim>& pos,
                           const std::array<std::array<stride_type,ndim>,N>& strides)
        {
            off0 = 0;
            for (unsigned i = 0;i < ndim;i++) off0 += pos[i]*strides[N-1][i];
        }
    };

    template <typename idx_type, typename stride_type, size_t ndim, size_t N, size_t I, typename Offset, typename... Offsets>
    struct set_offsets_helper
    {
        set_offsets_helper(Offset& off0, Offsets&... off,
                           const std::array<idx_type,ndim>& pos,
                           const std::array<std::array<stride_type,ndim>,N>& strides)
        {
            off0 = 0;
            for (unsigned i = 0;i < ndim;i++) off0 += pos[i]*strides[I-1][i];
            set_offsets_helper<idx_type, stride_type, ndim, N, I+1, Offsets...>(off..., pos, strides);
        }
    };

    template <typename idx_type, typename stride_type, size_t ndim, size_t N, typename... Offsets>
    void set_offsets(const std::array<idx_type,ndim>& pos,
                     const std::array<std::array<stride_type,ndim>,N>& strides,
                     Offsets&... off)
    {
        set_offsets_helper<idx_type, stride_type, ndim, N, 1, Offsets...>(off..., pos, strides);
    }

    template <typename stride_type, size_t ndim, size_t N, size_t I, typename Stride, typename... Strides>
    struct set_strides_helper;

    template <typename stride_type, size_t ndim, size_t N, typename Stride>
    struct set_strides_helper<stride_type, ndim, N, N, Stride>
    {
            set_strides_helper(const Stride& stride0,
                               std::array<std::array<stride_type,ndim>,N>& _strides)
        {
                assert(stride0.size() == ndim);
                std::copy_n(stride0.begin(), ndim, _strides[N-1].begin());
        }
    };

    template <typename stride_type, size_t ndim, size_t N, size_t I, typename Stride, typename... Strides>
    struct set_strides_helper
    {
            set_strides_helper(const Stride& stride0, const Strides&... strides,
                               std::array<std::array<stride_type,ndim>,N>& _strides)
        {
            assert(stride0.size() == ndim);
            std::copy_n(stride0.begin(), ndim, _strides[I-1].begin());
            set_strides_helper<stride_type, ndim, N, I+1, Strides...>(strides..., _strides);
        }
    };

    template <typename stride_type, size_t ndim, size_t N, typename... Strides>
    void set_strides(std::array<std::array<stride_type,ndim>,N>& _strides,
                     const Strides&... strides)
    {
        set_strides_helper<stride_type, ndim, N, 1, Strides...>(strides..., _strides);
    }
}

template <typename idx_type, typename stride_type, unsigned ndim, unsigned N=1>
class miterator
{
    public:
        miterator(const miterator&) = default;

        miterator(miterator&&) = default;

        template <typename Len, typename... Strides,
                  typename=typename std::enable_if<detail::is_container_of<idx_type, Len>::value &&
                                                   detail::are_containers_of<stride_type, Strides...>::value &&
                                                   sizeof...(Strides) == N>::type>
        miterator(const Len& len, const Strides&... strides)
        : _pos{}, _first(true), _empty(false)
        {
            assert(len.size() == ndim);
            for (unsigned i = 0;i < ndim;i++) if (len[i] == 0) _empty = true;
            std::copy_n(len.begin(), ndim, _len.begin());
            detail::set_strides(_strides, strides...);
        }

        miterator& operator=(const miterator&) = default;

        miterator& operator=(miterator&&) = default;

        void reset()
        {
            _pos.fill(0);
            _first = true;
        }

        template <typename... Offsets,
                  typename=typename std::enable_if<sizeof...(Offsets) == N>::type>
        bool next(Offsets&... off)
        {
            if (_empty) return false;

            if (_first)
            {
                _first = false;
                return true;
            }

            if (ndim == 0)
            {
                _first = true;
                return false;
            }

            for (unsigned i = 0;i < ndim;i++)
            {
                if (_pos[i] == _len[i]-1)
                {
                    detail::dec_offsets(i, _pos, _strides, off...);
                    _pos[i] = 0;

                    if (i == ndim-1)
                    {
                        _first = true;
                        return false;
                    }
                }
                else
                {
                    detail::inc_offsets(i, _strides, off...);
                    _pos[i]++;
                    return true;
                }
            }

            return true;
        }

        template <typename... Offsets,
                  typename=typename std::enable_if<sizeof...(Offsets) == N>::type>
        void position(stride_type pos, Offsets&... off)
        {
            for (size_t i = 0;i < ndim;i++)
            {
                _pos[i] = pos%_len[i];
                pos = pos/_len[i];
            }
            assert(pos == 0);

            position(_pos, off...);
        }

        template <typename Pos, typename... Offsets,
                  typename=typename std::enable_if<detail::is_container_of<idx_type, Pos>::value &&
                                                   sizeof...(Offsets) == N>::type>
        void position(const Pos& pos, Offsets&... off)
        {
            assert(pos.size() == ndim);

            _pos = pos;

            for (size_t i = 0;i < ndim;i++)
            {
                assert(_pos[i] >= 0 && _pos[i] < _len[i]);
            }

            detail::set_offsets(_pos, _strides, off...);

            _first = true;
        }

        unsigned dimension() const
        {
            return ndim;
        }

        idx_type position(unsigned dim) const
        {
            return _pos[dim];
        }

        const std::array<idx_type,ndim>& position() const
        {
            return _pos;
        }

        idx_type length(unsigned dim) const
        {
            return _len[dim];
        }

        const std::array<idx_type,ndim>& lengths() const
        {
            return _len;
        }

        stride_type stride(unsigned i, unsigned dim) const
        {
            return _strides[i][dim];
        }

        const std::array<stride_type,ndim>& strides(unsigned i) const
        {
            return _strides[i];
        }

        friend void swap(miterator& i1, miterator& i2)
        {
            using std::swap;
            swap(i1._pos, i2._pos);
            swap(i1._len, i2._len);
            swap(i1._strides, i2._strides);
            swap(i1._first, i2._first);
        }

    private:
        std::array<idx_type,ndim> _pos;
        std::array<idx_type,ndim> _len;
        std::array<std::array<stride_type,ndim>,N> _strides;
        bool _first;
        bool _empty;
};

template <typename idx_type, typename stride_type, size_t ndim, typename... Strides,
          typename=typename std::enable_if<detail::are_arrays_of<stride_type, ndim, Strides...>::value>::type>
miterator<idx_type, stride_type, ndim, 1+sizeof...(Strides)>
make_iterator(const std::array<idx_type, ndim>& len,
              const std::array<stride_type, ndim>& stride0,
              const Strides&... strides)
{
    return {len, stride0, strides...};
}

}

#endif
