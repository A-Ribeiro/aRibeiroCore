#ifndef __algorithms__h__
#define __algorithms__h__

#include <aRibeiroCore/common.h>

namespace aRibeiro {

    namespace PatternMatch {

        int32_t edit_distance(const char* s, const char* t);

    }

    namespace Sorting {

        // for counting sort, 
        // if you need more than 4.294.967.295 indexes, 
        // you can change this typedef
        typedef uint32_t counter_type;

        struct IndexInt32 {
            int32_t toSort;
            uint32_t index;

            static bool comparator(const IndexInt32 &i1, const IndexInt32 &i2) {
                return (i1.toSort < i2.toSort);
            }

            static IndexInt32 Create(uint32_t index, int32_t toSort) {
                IndexInt32 result;
                result.toSort = toSort;
                result.index = index;
                return result;
            }
        };

        struct IndexUInt32 {
            uint32_t toSort;
            uint32_t index;

            static bool comparator(const IndexUInt32 &i1, const IndexUInt32 &i2) {
                return (i1.toSort < i2.toSort);
            }

            static IndexInt32 Create(uint32_t index, uint32_t toSort) {
                IndexInt32 result;
                result.toSort = toSort;
                result.index = index;
                return result;
            }
        };


        ARIBEIRO_INLINE
        static uint32_t sort_spread_uint32(uint32_t min, uint32_t max, uint32_t v){
            v = v - min;
            uint64_t aux = v;
            uint64_t delta = max - min;
            aux = (aux * (uint64_t)UINT32_MAX) / delta;
            return (uint32_t)aux;
        }

        ARIBEIRO_INLINE
        static int32_t sort_spread_int32(int32_t min, int32_t max, int32_t v){
            uint64_t aux = (uint64_t)v - (uint64_t)min;
            uint64_t delta = (uint64_t)max - (uint64_t)min;
            aux = (aux * (uint64_t)UINT32_MAX) / delta;
            aux += (uint64_t)INT32_MIN;
            return (int32_t)aux;
        }


        // Float radix sort trick from: http://stereopsis.com/radix.html
        // ================================================================================================
        // flip a float for sorting
        //  finds SIGN of fp number.
        //  if it's 1 (negative float), it flips all bits
        //  if it's 0 (positive float), it flips the sign only
        // ================================================================================================
        ARIBEIRO_INLINE
        static uint32_t sort_float_to_uint32(const float &_f) {
            uint32_t f = *(uint32_t*)&_f;
            uint32_t mask = (-int32_t(f >> 31)) | 0x80000000;
            return f ^ mask;
        }

        ARIBEIRO_INLINE
        static int32_t sort_float_to_int32(const float& _f) {
            int32_t f = (int32_t)sort_float_to_uint32(_f) + INT32_MIN;
            return (int32_t)f;
        }

        // ================================================================================================
        // flip a float back (invert FloatFlip)
        //  signed was flipped from above, so:
        //  if sign is 1 (negative), it flips the sign bit back
        //  if sign is 0 (positive), it flips all bits back
        // ================================================================================================
        ARIBEIRO_INLINE
        static float sort_uint32_to_float(const uint32_t &f) {
            uint32_t mask = ((f >> 31) - 1) | 0x80000000;
            mask = f ^ mask;
            return *(float*)&mask;
        }
        ARIBEIRO_INLINE
        static float sort_int32_to_float(const int32_t& _f) {
            uint32_t f = (uint32_t)((uint32_t)_f - (uint32_t)INT32_MIN);
            return sort_uint32_to_float(f);
        }

        /// \brief Radix sort using Counting sort inside it for 'int32_t' type
        ///
        /// Radix sort implementation
        /// _________________________
        /// 
        /// The radix sort separates the number by its digit.
        /// 
        /// The natural base is 10 (didatic).
        /// 
        /// This algorithm uses the 'base_bits' to setup a binary base over the type 'int32_t'.
        /// 
        /// This binary base allow us to use the shift and bitwise operators to perform better on CPUs.
        /// 
        /// Counting sort implementation
        /// _________________________
        /// 
        /// The counting sort is used to sort the individual digits resulting from the radix.
        /// 
        /// This implementation allow the pre-allocation of the temporary buffer outside of the function.
        /// 
        /// If the 'tmp_array' is NULL, then the array is allocated and freed according the 'arrSize' parameter.
        ///
        /// Example:
        ///
        /// \code
        /// \endcode
        ///
        /// \author Alessandro Ribeiro
        ///
        void radix_counting_sort_signed(int32_t* _arr, uint32_t arrSize, int32_t* tmp_array = NULL);
        void radix_counting_sort_signed_index(IndexInt32* _arr, uint32_t arrSize, IndexInt32* tmp_array = NULL);

        /// \brief Radix sort using Counting sort inside it for 'uint32_t' type
        ///
        /// Radix sort implementation
        /// _________________________
        /// 
        /// The radix sort separates the number by its digit.
        /// 
        /// The natural base is 10 (didatic).
        /// 
        /// This algorithm uses the 'base_bits' to setup a binary base over the type 'uint32_t'.
        /// 
        /// This binary base allow us to use the shift and bitwise operators to perform better on CPUs.
        /// 
        /// Counting sort implementation
        /// _________________________
        /// 
        /// The counting sort is used to sort the individual digits resulting from the radix.
        /// 
        /// This implementation allow the pre-allocation of the temporary buffer outside of the function.
        /// 
        /// If the 'tmp_array' is NULL, then the array is allocated and freed according the 'arrSize' parameter.
        ///
        /// Example:
        ///
        /// \code
        /// \endcode
        ///
        /// \author Alessandro Ribeiro
        ///
        void radix_counting_sort_unsigned(uint32_t* _arr, uint32_t arrSize, uint32_t* tmp_array = NULL);
        void radix_counting_sort_unsigned_index(IndexUInt32* _arr, uint32_t arrSize, IndexUInt32* tmp_array = NULL);

        /// \brief Radix sort using Bucket sort inside it for 'int32_t' type
        ///
        /// Radix sort implementation
        /// _________________________
        /// 
        /// The radix sort separates the number by its digit.
        /// 
        /// The natural base is 10 (didatic).
        /// 
        /// This algorithm uses the 'base_bits' to setup a binary base over the type 'int32_t'.
        /// 
        /// This binary base allow us to use the shift and bitwise operators to perform better on CPUs.
        /// 
        /// Bucket sort implementation
        /// _________________________
        /// 
        /// The bucket sort is used to sort the individual digits resulting from the radix.
        /// 
        /// The main drawback of using bucket sort is that it allocates a huge number of individual lists.
        /// 
        /// But for some cases, it performs better than the Counting sort implementation.
        ///
        /// Example:
        ///
        /// \code
        /// \endcode
        ///
        /// \author Alessandro Ribeiro
        ///
        void radix_bucket_sort_signed(int32_t* arr, uint32_t arrSize);
        void radix_bucket_sort_signed_index(IndexInt32* arr, uint32_t arrSize);


        /// \brief Radix sort using Bucket sort inside it for 'uint32_t' type
        ///
        /// Radix sort implementation
        /// _________________________
        /// 
        /// The radix sort separates the number by its digit.
        /// 
        /// The natural base is 10 (didatic).
        /// 
        /// This algorithm uses the 'base_bits' to setup a binary base over the type 'uint32_t'.
        /// 
        /// This binary base allow us to use the shift and bitwise operators to perform better on CPUs.
        /// 
        /// Bucket sort implementation
        /// _________________________
        /// 
        /// The bucket sort is used to sort the individual digits resulting from the radix.
        /// 
        /// The main drawback of using bucket sort is that it allocates a huge number of individual lists.
        /// 
        /// But for some cases, it performs better than the Counting sort implementation.
        ///
        /// Example:
        ///
        /// \code
        /// \endcode
        ///
        /// \author Alessandro Ribeiro
        ///
        void radix_bucket_sort_unsigned(uint32_t* arr, uint32_t arrSize);
        void radix_bucket_sort_unsigned_index(IndexUInt32* arr, uint32_t arrSize);

        /// \brief Inplace replace
        ///
        /// Do the output buffer sort internal replace.
        /// No need to create a double buffer approach.
        ///
        /// Example:
        ///
        /// \code
        /// std::vector<int32_t> original_buffer();
        ///
        /// std::vector<Sorting::IndexInt32> sort_index_buffer(original_buffer.size());
        /// for (int i = 0; i < (int)sort_index_buffer.size(); i++)
        /// {
        ///     sort_index_buffer[i].index = i;
        ///     sort_index_buffer[i].toSort = original_buffer[i];
        /// }
        /// Sorting::radix_counting_sort_signed_index(&sort_index_buffer[0], sort_index_buffer.size());
        ///
        /// // Replace the original buffer with the sorted version
        /// Sorting::inplace_replace(&sort_index_buffer[0], &original_buffer[0], original_buffer.size());
        ///
        /// \endcode
        ///
        /// \author Alessandro Ribeiro
        ///
        template <typename T, typename U>
        void inplace_replace(T *index_buffer, U *output, uint32_t size)
        {
            std::vector<uint32_t> inplace_redirect(size);
            for (uint32_t i = 0; i < size; i++)
            {
                // make a swap

                uint32_t src_index = index_buffer[i].index;
                uint32_t target_index = i;

                if (src_index == target_index)
                    continue;

                // redirect walk
                while (src_index < target_index)
                    src_index = inplace_redirect[src_index];

                // swap
                U save = output[target_index];
                output[target_index] = output[src_index];
                output[src_index] = save;

                // save redirection
                inplace_redirect[target_index] = src_index;
            }
        }
    }
}

#endif
