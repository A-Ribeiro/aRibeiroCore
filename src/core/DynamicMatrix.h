#ifndef dynamic_matrix_h___
#define dynamic_matrix_h___

#include <aRibeiroCore/common.h>

#if !defined(NDEBUG)
#define MATRIX_THROW_OUT_OF_BOUND_EXCEPTION
#endif

namespace aRibeiro
{

#if defined(MATRIX_THROW_OUT_OF_BOUND_EXCEPTION)
    template <typename T>
    class DynamicMatrix;

    template <typename _T>
    class DynamicMatrix_Proxy
    {
        int width;
        _T *ptr;
        DynamicMatrix_Proxy(int width, _T *ptr)
        {
            this->width = width;
            this->ptr = ptr;
        }

    public:
        _T &operator[](int x)
        {
            if (x < 0 || x >= width)
                throw std::runtime_error("x out of bounds exception.");
            return this->ptr[x];
        }
        const _T &operator[](int x) const
        {
            if (x < 0 || x >= width)
                throw std::runtime_error("x out of bounds exception.");
            return this->ptr[x];
        }
        friend class DynamicMatrix<_T>;
    };
#endif

    template <typename T>
    class DynamicMatrix
    {

    public:
        T *array;
        int width;
        int height;

        DynamicMatrix() : array(NULL), width(0), height(0)
        {
            setSize(32, 32);
        }

        DynamicMatrix(int width, int height) : array(NULL), width(0), height(0)
        {
            setSize(width, height);
        }

        DynamicMatrix(const DynamicMatrix &m) : array(NULL), width(0), height(0)
        {
            setSize(m.width, m.height);
            memcpy(array, m.array, sizeof(T) * width * height);
        }
        void operator=(const DynamicMatrix &m)
        {
            setSize(m.width, m.height);
            memcpy(array, m.array, sizeof(T) * width * height);
        }

#if defined(MATRIX_THROW_OUT_OF_BOUND_EXCEPTION)
        DynamicMatrix_Proxy<T> operator[](int y)
        {
            if (y < 0 || y >= height)
                throw std::runtime_error("y out of bounds exception.");
            return DynamicMatrix_Proxy<T>(width, &this->array[y * width]);
        }

        const DynamicMatrix_Proxy<T> operator[](int y) const
        {
            if (y < 0 || y >= height)
                throw std::runtime_error("y out of bounds exception.");
            return DynamicMatrix_Proxy<T>(width, &this->array[y * width]);
        }
#else
    T *operator[](int y)
    {
        return &this->array[y * width];
    }

    const T *operator[](int y) const
    {
        return &this->array[y * width];
    }
#endif

        ~DynamicMatrix()
        {
            if (array != NULL)
            {
                delete[] array;
                array = NULL;
            }
            this->width = 0;
            this->height = 0;
        }

        void setSize(int _w, int _h)
        {
            if (width == _w && height == _h)
                return;

            ARIBEIRO_ABORT(_w <= 0 || _h <= 0, "invalid matrix size: %i %i", _w, _h);

            // if (_w <= 0 || _h <= 0)
            // {
            //     _w = 32;
            //     _h = 32;
            // }

            width = _w;
            height = _h;

            if (array != NULL)
            {
                delete[] array;
                array = NULL;
            }

            array = new T[width * height];
        }

        void clear(const T &v)
        {
            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    array[x + y * width] = v;
                }
            }
        }
    };
}

#endif