#include "Algorithms.h"

namespace aRibeiro {

	namespace PatternMatch {

		ARIBEIRO_INLINE
			static int32_t intMin(int32_t a, int32_t b, int32_t c) {
#ifdef ARIBEIRO_SSE2

			__m128i _a = _mm_set1_epi32(a);
			__m128i _b = _mm_set1_epi32(b);
			__m128i _c = _mm_set1_epi32(c);

			__m128i _min = _mm_min_epi32(_a, _b);
			_min = _mm_min_epi32(_min, _c);

			return _mm_i32_(_min, 0);
#else
			int32_t min = (a < b) ? a : b;
			min = (min < c) ? min : c;
			return min;
#endif
		}

		int32_t edit_distance(const char* s, const char* t) {

			uint32_t s_len = (uint32_t)strlen(s);
			uint32_t t_len = (uint32_t)strlen(t);

			if (s_len < t_len)
				return edit_distance(t, s);

			int32_t result = 0;
			if (strcmp(s, t) == 0)
				return 0;
			else if (s_len == 0)
				return t_len;
			else if (t_len == 0)
				return s_len;

			std::vector<int32_t> row_prev(t_len + 1);
			std::vector<int32_t> row_curr(t_len + 1);

			for (int32_t i = 0; i < row_prev.size(); i++)
				row_prev[i] = i;

			for (int32_t i = 0; i < s_len; i++) {
				row_curr[0] = i + 1;
				for (int32_t j = 0; j < t_len; j++) {
					int32_t insertions = row_prev[j + 1] + 1;// j+1 instead of j since previous_row and current_row are one character longer
					int32_t deletions = row_curr[j] + 1;     // than s2
					int32_t substitutions = row_prev[j] + (s[i] != t[j]);
					row_curr[j + 1] = intMin(insertions, deletions, substitutions);
				}
				memcpy(&row_prev[0], &row_curr[0], sizeof(int32_t) * row_prev.size());
			}

			return row_curr[t_len];
		}

	}

	namespace Sorting {

		// The base implementation for all algorithms started from the stackoverflow discussion here:
		// https://stackoverflow.com/questions/15306665/radix-sort-for-negative-integers

		void radix_counting_sort_signed(int32_t* _arr, uint32_t arrSize, int32_t* tmp_array) {

			////
			//// OLD ALGORITHM - very fast, but can be better!!!
			//// 
			//// 

			//if (arrSize == 0)
			//	return;

			//// best performance on 8 bit of symbols...
			//const int32_t base_bits = 8;//1 byte
			//const int32_t max_bucket_symbols = 1 << base_bits;//256
			//const int32_t base_mask = max_bucket_symbols - 1;//0xff

			//const int32_t int_total_bytes = sizeof(int32_t);
			//const uint32_t int_mask_bits = 0xffffffff;
			//const int32_t int_total_bits = int_total_bytes << 3;// *8;

			//const int32_t radix_num = (int_total_bits / base_bits) + ((int_total_bits % base_bits) ? 1 : 0);
			//const uint32_t last_radix_mask = int_mask_bits >> ((radix_num - 1) * base_bits);

			//int32_t shift = 0;

			////Counting Sort
			//counter_type counting[max_bucket_symbols];
			//int32_t* aux;
			//if (tmp_array == NULL)
			//	aux = (int32_t*)malloc_aligned(arrSize * sizeof(int32_t));
			//else
			//	aux = tmp_array;

			//int32_t* in = _arr;
			//int32_t* out = aux;

			//for (int32_t i = 0; i < radix_num; i++) {

			//	// Cleaning counters
			//	memset(counting, 0, sizeof(counter_type) * max_bucket_symbols);

			//	// count the elements
			//	for (uint32_t j = 0; j < arrSize; j++) {
			//		int32_t currItem = in[j];
			//		uint32_t bucket_index = (((uint32_t)currItem >> shift) & base_mask);//0xff
			//		counting[bucket_index]++;
			//	}

			//	//compute offsets

			//	// Treat the last radix with sign bit specially
			//	// Output signed groups (128..256 = -128..-1) first
			//	// Other groups afterwards. No performance penalty, as compared to flipping sign bit
			//	// via (($currItem ^ 0x8000000000000000) >> $shift) & 0xFF)
			//	if (i == radix_num - 1) {
			//		uint32_t max = last_radix_mask;
			//		uint32_t max_half = (max >> 1) + 1;

			//		counter_type acc = counting[max_half];
			//		counting[max_half] = 0;
			//		for (uint32_t j = max_half + 1; j <= max; j++) {
			//			counter_type tmp = counting[j];
			//			counting[j] = acc;
			//			acc += tmp;
			//		}

			//		for (uint32_t j = 0; j < max_half; j++) {
			//			counter_type tmp = counting[j];
			//			counting[j] = acc;
			//			acc += tmp;
			//		}
			//	}
			//	else {
			//		counter_type acc = counting[0];
			//		counting[0] = 0;
			//		for (uint32_t j = 1; j < max_bucket_symbols; j++) {
			//			counter_type tmp = counting[j];
			//			counting[j] = acc;
			//			acc += tmp;
			//		}
			//	}

			//	// place elements in the output array
			//	for (uint32_t j = 0; j < arrSize; j++) {
			//		int32_t currItem = in[j];
			//		uint32_t bucket_index = (((uint32_t)currItem >> shift) & base_mask);//0xff
			//		counter_type out_index = counting[bucket_index];
			//		counting[bucket_index]++;
			//		out[out_index] = currItem;
			//	}

			//	//swap out, in
			//	int32_t* tmp = in;
			//	in = out;
			//	out = tmp;

			//	// Change shift value for next iterations
			//	shift += base_bits;
			//}

			//if (in != _arr)
			//	memcpy(_arr, in, arrSize * sizeof(int32_t));

			//if (tmp_array == NULL)
			//	free_aligned(aux);

			//*************************
			//* New Uber Optimized!!! *
			//*************************


#if defined(ARIBEIRO_SSE2)


			if (arrSize == 0)
				return;

			//Counting Sort
			__m128i counting[256];
			int32_t* aux;

			if (tmp_array == NULL)
				aux = (int32_t*)malloc_aligned(arrSize * sizeof(int32_t));
			else
				aux = tmp_array;

			int32_t* in = _arr;
			int32_t* out = aux;

			//for (int32_t i = 0; i < radix_num; i++) 
			{

				// Cleaning counters
				memset(counting, 0, sizeof(__m128i) * 256);


				const __m128i _c_0 = _mm_setr_epi32(1, 0, 0, 0);
				const __m128i _c_1 = _mm_setr_epi32(0, 1, 0, 0);
				const __m128i _c_2 = _mm_setr_epi32(0, 0, 1, 0);
				const __m128i _c_3 = _mm_setr_epi32(0, 0, 0, 1);

				// count the elements
				const uint32_t _minus_128 = (uint32_t)(uint8_t)-128;
				for (uint32_t j = 0; j < arrSize; j++) {


					uint32_t currItem = (uint32_t)in[j];
					uint32_t index_0 = (currItem) & 0xff;
					uint32_t index_1 = (currItem >> 8) & 0xff;
					uint32_t index_2 = (currItem >> 16) & 0xff;
					uint32_t index_3 = ((currItem >> 24) + _minus_128) & 0xff;

					_mm_i32_(counting[index_0], 0)++;
					_mm_i32_(counting[index_1], 1)++;
					_mm_i32_(counting[index_2], 2)++;
					_mm_i32_(counting[index_3], 3)++;

					/*counting[index_0] = _mm_add_epi32(counting[index_0], _c_0);
					counting[index_1] = _mm_add_epi32(counting[index_1], _c_1);
					counting[index_2] = _mm_add_epi32(counting[index_2], _c_2);
					counting[index_3] = _mm_add_epi32(counting[index_3], _c_3);*/
				}

				//compute offsets
				__m128i acc = counting[0];

				counting[0] = _mm_set1_epi32(0);

				for (uint32_t j = 1; j < 256; j++) {
					__m128i tmp = counting[j];
					counting[j] = acc;
					acc = _mm_add_epi32(acc, tmp);
				}

				//shift = 0
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						int32_t currItem = in[j];
						uint32_t bucket_index = (((uint32_t)currItem) & 0xff);

						uint32_t out_index = _mm_i32_(counting[bucket_index], 0);
						counting[bucket_index] = _mm_add_epi32(counting[bucket_index], _c_0);

						out[out_index] = currItem;
					}

					//swap out, in
					int32_t* tmp = in;
					in = out;
					out = tmp;
				}

				//shift = 8
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						int32_t currItem = in[j];
						uint32_t bucket_index = (((uint32_t)currItem >> 8) & 0xff);

						uint32_t out_index = _mm_i32_(counting[bucket_index], 1);

						counting[bucket_index] = _mm_add_epi32(counting[bucket_index], _c_1);

						out[out_index] = currItem;
					}

					//swap out, in
					int32_t* tmp = in;
					in = out;
					out = tmp;
				}

				//shift = 16
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						int32_t currItem = in[j];
						uint32_t bucket_index = (((uint32_t)currItem >> 16) & 0xff);

						uint32_t out_index = _mm_i32_(counting[bucket_index], 2);

						counting[bucket_index] = _mm_add_epi32(counting[bucket_index], _c_2);

						out[out_index] = currItem;
					}

					//swap out, in
					int32_t* tmp = in;
					in = out;
					out = tmp;
				}

				//shift = 24
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						int32_t currItem = in[j];
						uint32_t bucket_index = ((((int32_t)currItem >> 24) + _minus_128) & 0xff);

						uint32_t out_index = _mm_i32_(counting[bucket_index], 3);

						counting[bucket_index] = _mm_add_epi32(counting[bucket_index], _c_3);

						out[out_index] = currItem;
					}

					//swap out, in
					int32_t* tmp = in;
					in = out;
					out = tmp;
				}

			}

			if (tmp_array == NULL)
				free_aligned(aux);

#else

			if (arrSize == 0)
				return;

			//Counting Sort
			uint32_t counting[256][4];
			int32_t* aux;

			if (tmp_array == NULL)
				aux = (int32_t*)malloc_aligned(arrSize * sizeof(int32_t));
			else
				aux = tmp_array;

			int32_t* in = _arr;
			int32_t* out = aux;

			//for (int32_t i = 0; i < radix_num; i++) 
			{

				// Cleaning counters
				memset(&counting[0][0], 0, sizeof(uint32_t) * 256 * 4);

				// count the elements
				const uint32_t _minus_128 = (uint32_t)(uint8_t)-128;
				for (uint32_t j = 0; j < arrSize; j++) {
					uint32_t currItem = (uint32_t)in[j];
					counting[(currItem) & 0xff][0]++;
					counting[(currItem >> 8) & 0xff][1]++;
					counting[(currItem >> 16) & 0xff][2]++;
					counting[((currItem >> 24) + _minus_128) & 0xff][3]++;
				}

				//compute offsets
				//uint32_t acc[4];
				//memcpy(acc, counting, sizeof(uint32_t) << 2);
				//memset(counting, 0, sizeof(uint32_t) << 2);

				uint32_t acc[4] = { counting[0][0], counting[0][1], counting[0][2], counting[0][3] };

				counting[0][0] = 0;
				counting[0][1] = 0;
				counting[0][2] = 0;
				counting[0][3] = 0;

				//uint32_t tmp[4];

				for (uint32_t j = 1; j < 256; j++) {

					//memcpy(tmp, counting[j], sizeof(uint32_t) << 2);
					//memcpy(counting[j], acc, sizeof(uint32_t) << 2);

					uint32_t tmp[4] = { counting[j][0], counting[j][1], counting[j][2], counting[j][3] };

					counting[j][0] = acc[0];
					counting[j][1] = acc[1];
					counting[j][2] = acc[2];
					counting[j][3] = acc[3];

					acc[0] += tmp[0];
					acc[1] += tmp[1];
					acc[2] += tmp[2];
					acc[3] += tmp[3];
				}

				//for (uint32_t shift = 0; shift < 32; shift += 8) 
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						int32_t currItem = in[j];
						uint32_t bucket_index = (((uint32_t)currItem >> 0) & 0xff);
						uint32_t out_index = counting[bucket_index][0];
						counting[bucket_index][0]++;
						out[out_index] = currItem;
					}

					//swap out, in
					int32_t* tmp = in;
					in = out;
					out = tmp;
				}
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						int32_t currItem = in[j];
						uint32_t bucket_index = (((uint32_t)currItem >> 8) & 0xff);
						uint32_t out_index = counting[bucket_index][1];
						counting[bucket_index][1]++;
						out[out_index] = currItem;
					}

					//swap out, in
					int32_t* tmp = in;
					in = out;
					out = tmp;
				}
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						int32_t currItem = in[j];
						uint32_t bucket_index = (((uint32_t)currItem >> 16) & 0xff);
						uint32_t out_index = counting[bucket_index][2];
						counting[bucket_index][2]++;
						out[out_index] = currItem;
					}

					//swap out, in
					int32_t* tmp = in;
					in = out;
					out = tmp;
				}
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						int32_t currItem = in[j];
						uint32_t bucket_index = ((((uint32_t)currItem >> 24) + _minus_128) & 0xff);
						uint32_t out_index = counting[bucket_index][3];
						counting[bucket_index][3]++;
						out[out_index] = currItem;
					}

					//swap out, in
					int32_t* tmp = in;
					in = out;
					out = tmp;
				}
			}

			if (tmp_array == NULL)
				free_aligned(aux);

#endif
		}

		void radix_counting_sort_signed_index(IndexInt32* _arr, uint32_t arrSize, IndexInt32* tmp_array) {

			////
			//// OLD ALGORITHM - very fast, but can be better!!!
			//// 
			////

			//if (arrSize == 0)
			//	return;

			//// best performance on 8 bit of symbols...
			//const int32_t base_bits = 8;//1 byte
			//const int32_t max_bucket_symbols = 1 << base_bits;//256
			//const int32_t base_mask = max_bucket_symbols - 1;//0xff

			//const int32_t int_total_bytes = sizeof(int32_t);
			//const uint32_t int_mask_bits = 0xffffffff;
			//const int32_t int_total_bits = int_total_bytes << 3;// *8;

			//const int32_t radix_num = (int_total_bits / base_bits) + ((int_total_bits % base_bits) ? 1 : 0);
			//const uint32_t last_radix_mask = int_mask_bits >> ((radix_num - 1) * base_bits);

			//int32_t shift = 0;

			////Counting Sort
			//counter_type counting[max_bucket_symbols];
			//IndexInt32* aux;
			//if (tmp_array == NULL)
			//	aux = (IndexInt32*)malloc_aligned(arrSize * sizeof(IndexInt32));
			//else
			//	aux = tmp_array;

			//IndexInt32* in = _arr;
			//IndexInt32* out = aux;

			//for (int32_t i = 0; i < radix_num; i++) {

			//	// Cleaning counters
			//	memset(counting, 0, sizeof(counter_type) * max_bucket_symbols);

			//	// count the elements
			//	for (uint32_t j = 0; j < arrSize; j++) {
			//		const IndexInt32& currItem = in[j];
			//		uint32_t bucket_index = (((uint32_t)currItem.toSort >> shift) & base_mask);//0xff
			//		counting[bucket_index]++;
			//	}

			//	//compute offsets

			//	// Treat the last radix with sign bit specially
			//	// Output signed groups (128..256 = -128..-1) first
			//	// Other groups afterwards. No performance penalty, as compared to flipping sign bit
			//	// via (($currItem ^ 0x8000000000000000) >> $shift) & 0xFF)
			//	if (i == radix_num - 1) {
			//		uint32_t max = last_radix_mask;
			//		uint32_t max_half = (max >> 1) + 1;

			//		counter_type acc = counting[max_half];
			//		counting[max_half] = 0;
			//		for (uint32_t j = max_half + 1; j <= max; j++) {
			//			counter_type tmp = counting[j];
			//			counting[j] = acc;
			//			acc += tmp;
			//		}

			//		for (uint32_t j = 0; j < max_half; j++) {
			//			counter_type tmp = counting[j];
			//			counting[j] = acc;
			//			acc += tmp;
			//		}
			//	}
			//	else {
			//		counter_type acc = counting[0];
			//		counting[0] = 0;
			//		for (uint32_t j = 1; j < max_bucket_symbols; j++) {
			//			counter_type tmp = counting[j];
			//			counting[j] = acc;
			//			acc += tmp;
			//		}
			//	}

			//	// place elements in the output array
			//	for (uint32_t j = 0; j < arrSize; j++) {
			//		const IndexInt32& currItem = in[j];
			//		uint32_t bucket_index = (((uint32_t)currItem.toSort >> shift) & base_mask);//0xff
			//		counter_type out_index = counting[bucket_index];
			//		counting[bucket_index]++;
			//		out[out_index] = currItem;
			//	}

			//	//swap out, in
			//	IndexInt32* tmp = in;
			//	in = out;
			//	out = tmp;

			//	// Change shift value for next iterations
			//	shift += base_bits;
			//}

			//if (in != _arr)
			//	memcpy(_arr, in, arrSize * sizeof(IndexInt32));

			//if (tmp_array == NULL)
			//	free_aligned(aux);



			//*************************
			//* New Uber Optimized!!! *
			//*************************


#if defined(ARIBEIRO_SSE2)


			if (arrSize == 0)
				return;

			//Counting Sort
			__m128i counting[256];
			IndexInt32* aux;

			if (tmp_array == NULL)
				aux = (IndexInt32*)malloc_aligned(arrSize * sizeof(IndexInt32));
			else
				aux = tmp_array;

			IndexInt32* in = _arr;
			IndexInt32* out = aux;

			//for (int32_t i = 0; i < radix_num; i++) 
			{

				// Cleaning counters
				memset(counting, 0, sizeof(__m128i) * 256);


				const __m128i _c_0 = _mm_setr_epi32(1, 0, 0, 0);
				const __m128i _c_1 = _mm_setr_epi32(0, 1, 0, 0);
				const __m128i _c_2 = _mm_setr_epi32(0, 0, 1, 0);
				const __m128i _c_3 = _mm_setr_epi32(0, 0, 0, 1);

				// count the elements
				const uint32_t _minus_128 = (uint32_t)(uint8_t)-128;
				for (uint32_t j = 0; j < arrSize; j++) {


					uint32_t currItem = (uint32_t)in[j].toSort;
					uint32_t index_0 = (currItem) & 0xff;
					uint32_t index_1 = (currItem >> 8) & 0xff;
					uint32_t index_2 = (currItem >> 16) & 0xff;
					uint32_t index_3 = ((currItem >> 24) + _minus_128) & 0xff;

					_mm_i32_(counting[index_0], 0)++;
					_mm_i32_(counting[index_1], 1)++;
					_mm_i32_(counting[index_2], 2)++;
					_mm_i32_(counting[index_3], 3)++;

					/*counting[index_0] = _mm_add_epi32(counting[index_0], _c_0);
					counting[index_1] = _mm_add_epi32(counting[index_1], _c_1);
					counting[index_2] = _mm_add_epi32(counting[index_2], _c_2);
					counting[index_3] = _mm_add_epi32(counting[index_3], _c_3);*/
				}

				//compute offsets
				__m128i acc = counting[0];

				counting[0] = _mm_set1_epi32(0);

				for (uint32_t j = 1; j < 256; j++) {
					__m128i tmp = counting[j];
					counting[j] = acc;
					acc = _mm_add_epi32(acc, tmp);
				}

				//shift = 0
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						const IndexInt32& currItem = in[j];
						uint32_t bucket_index = (((uint32_t)currItem.toSort) & 0xff);

						uint32_t out_index = _mm_i32_(counting[bucket_index], 0);
						counting[bucket_index] = _mm_add_epi32(counting[bucket_index], _c_0);

						out[out_index] = currItem;
					}

					//swap out, in
					IndexInt32* tmp = in;
					in = out;
					out = tmp;
				}

				//shift = 8
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						const IndexInt32& currItem = in[j];
						uint32_t bucket_index = (((uint32_t)currItem.toSort >> 8) & 0xff);

						uint32_t out_index = _mm_i32_(counting[bucket_index], 1);

						counting[bucket_index] = _mm_add_epi32(counting[bucket_index], _c_1);

						out[out_index] = currItem;
					}

					//swap out, in
					IndexInt32* tmp = in;
					in = out;
					out = tmp;
				}

				//shift = 16
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						const IndexInt32& currItem = in[j];
						uint32_t bucket_index = (((uint32_t)currItem.toSort >> 16) & 0xff);

						uint32_t out_index = _mm_i32_(counting[bucket_index], 2);

						counting[bucket_index] = _mm_add_epi32(counting[bucket_index], _c_2);

						out[out_index] = currItem;
					}

					//swap out, in
					IndexInt32* tmp = in;
					in = out;
					out = tmp;
				}

				//shift = 24
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						const IndexInt32& currItem = in[j];
						uint32_t bucket_index = ((((int32_t)currItem.toSort >> 24) + _minus_128) & 0xff);

						uint32_t out_index = _mm_i32_(counting[bucket_index], 3);

						counting[bucket_index] = _mm_add_epi32(counting[bucket_index], _c_3);

						out[out_index] = currItem;
					}

					//swap out, in
					IndexInt32* tmp = in;
					in = out;
					out = tmp;
				}

			}

			if (tmp_array == NULL)
				free_aligned(aux);

#else

			if (arrSize == 0)
				return;

			//Counting Sort
			uint32_t counting[256][4];
			IndexInt32* aux;

			if (tmp_array == NULL)
				aux = (IndexInt32*)malloc_aligned(arrSize * sizeof(IndexInt32));
			else
				aux = tmp_array;

			IndexInt32* in = _arr;
			IndexInt32* out = aux;

			//for (int32_t i = 0; i < radix_num; i++) 
			{

				// Cleaning counters
				memset(&counting[0][0], 0, sizeof(uint32_t) * 256 * 4);

				// count the elements
				const uint32_t _minus_128 = (uint32_t)(uint8_t)-128;
				for (uint32_t j = 0; j < arrSize; j++) {
					uint32_t currItem = (uint32_t)in[j].toSort;
					counting[(currItem) & 0xff][0]++;
					counting[(currItem >> 8) & 0xff][1]++;
					counting[(currItem >> 16) & 0xff][2]++;
					counting[((currItem >> 24) + _minus_128) & 0xff][3]++;
				}

				//compute offsets
				//uint32_t acc[4];
				//memcpy(acc, counting, sizeof(uint32_t) << 2);
				//memset(counting, 0, sizeof(uint32_t) << 2);

				uint32_t acc[4] = { counting[0][0], counting[0][1], counting[0][2], counting[0][3] };

				counting[0][0] = 0;
				counting[0][1] = 0;
				counting[0][2] = 0;
				counting[0][3] = 0;

				//uint32_t tmp[4];

				for (uint32_t j = 1; j < 256; j++) {

					//memcpy(tmp, counting[j], sizeof(uint32_t) << 2);
					//memcpy(counting[j], acc, sizeof(uint32_t) << 2);

					uint32_t tmp[4] = { counting[j][0], counting[j][1], counting[j][2], counting[j][3] };

					counting[j][0] = acc[0];
					counting[j][1] = acc[1];
					counting[j][2] = acc[2];
					counting[j][3] = acc[3];

					acc[0] += tmp[0];
					acc[1] += tmp[1];
					acc[2] += tmp[2];
					acc[3] += tmp[3];
				}

				//for (uint32_t shift = 0; shift < 32; shift += 8) 
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						const IndexInt32& currItem = in[j];
						uint32_t bucket_index = (((uint32_t)currItem.toSort >> 0) & 0xff);
						uint32_t out_index = counting[bucket_index][0];
						counting[bucket_index][0]++;
						out[out_index] = currItem;
					}

					//swap out, in
					IndexInt32* tmp = in;
					in = out;
					out = tmp;
				}
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						const IndexInt32& currItem = in[j];
						uint32_t bucket_index = (((uint32_t)currItem.toSort >> 8) & 0xff);
						uint32_t out_index = counting[bucket_index][1];
						counting[bucket_index][1]++;
						out[out_index] = currItem;
					}

					//swap out, in
					IndexInt32* tmp = in;
					in = out;
					out = tmp;
				}
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						const IndexInt32& currItem = in[j];
						uint32_t bucket_index = (((uint32_t)currItem.toSort >> 16) & 0xff);
						uint32_t out_index = counting[bucket_index][2];
						counting[bucket_index][2]++;
						out[out_index] = currItem;
					}

					//swap out, in
					IndexInt32* tmp = in;
					in = out;
					out = tmp;
				}
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						const IndexInt32& currItem = in[j];
						uint32_t bucket_index = ((((uint32_t)currItem.toSort >> 24) + _minus_128) & 0xff);
						uint32_t out_index = counting[bucket_index][3];
						counting[bucket_index][3]++;
						out[out_index] = currItem;
					}

					//swap out, in
					IndexInt32* tmp = in;
					in = out;
					out = tmp;
				}
			}

			if (tmp_array == NULL)
				free_aligned(aux);

#endif

		}



		void radix_counting_sort_unsigned(uint32_t* _arr, uint32_t arrSize, uint32_t* tmp_array) {
			//
			// OLD ALGORITHM - very fast, but can be better!!!
			// 
			// 

			//if (arrSize == 0)
			//    return;

			//// best performance on 8 bit of symbols...
			//const int32_t base_bits = 8;//1 byte
			//const int32_t max_bucket_symbols = 1 << base_bits;//256
			//const int32_t base_mask = max_bucket_symbols - 1;//0xff

			//const int32_t int_total_bytes = sizeof(uint32_t);
			//const uint32_t int_mask_bits = 0xffffffff;
			//const int32_t int_total_bits = int_total_bytes << 3;// *8;

			//const int32_t radix_num = (int_total_bits / base_bits) + ((int_total_bits % base_bits) ? 1 : 0);
			//const uint32_t last_radix_mask = int_mask_bits >> ((radix_num - 1) * base_bits);

			//int32_t shift = 0;

			////Counting Sort
			//counter_type counting[max_bucket_symbols];
			//uint32_t* aux;

			//if (tmp_array == NULL)
			//    aux = (uint32_t*)malloc_aligned(arrSize * sizeof(uint32_t));
			//else
			//    aux = tmp_array;

			//uint32_t* in = _arr;
			//uint32_t* out = aux;

			//for (int32_t i = 0; i < radix_num; i++) {

			//    // Cleaning counters
			//    memset(counting, 0, sizeof(counter_type) * max_bucket_symbols);

			//    // count the elements
			//    for (uint32_t j = 0; j < arrSize; j++) {
			//        uint32_t currItem = in[j];
			//        uint32_t bucket_index = (((uint32_t)currItem >> shift) & base_mask);//0xff
			//        counting[bucket_index]++;
			//    }

			//    //compute offsets
			//    counter_type acc = counting[0];
			//    counting[0] = 0;
			//    for (uint32_t j = 1; j < max_bucket_symbols; j++) {
			//        counter_type tmp = counting[j];
			//        counting[j] = acc;
			//        acc += tmp;
			//    }

			//    // place elements in the output array
			//    for (uint32_t j = 0; j < arrSize; j++) {
			//        uint32_t currItem = in[j];
			//        uint32_t bucket_index = (((uint32_t)currItem >> shift) & base_mask);//0xff
			//        counter_type out_index = counting[bucket_index];
			//        counting[bucket_index]++;
			//        out[out_index] = currItem;
			//    }

			//    //swap out, in
			//    uint32_t* tmp = in;
			//    in = out;
			//    out = tmp;

			//    // Change shift value for next iterations
			//    shift += base_bits;
			//}

			//if (in != _arr)
			//    memcpy(_arr, in, arrSize * sizeof(uint32_t));

			//if (tmp_array == NULL)
			//    free_aligned(aux);

			//*************************
			//* New Uber Optimized!!! *
			//*************************


#if defined(ARIBEIRO_SSE2)


			if (arrSize == 0)
				return;

			//Counting Sort
			__m128i counting[256];
			uint32_t* aux;

			if (tmp_array == NULL)
				aux = (uint32_t*)malloc_aligned(arrSize * sizeof(uint32_t));
			else
				aux = tmp_array;

			uint32_t* in = _arr;
			uint32_t* out = aux;

			//for (int32_t i = 0; i < radix_num; i++) 
			{

				// Cleaning counters
				memset(counting, 0, sizeof(__m128i) * 256);


				const __m128i _c_0 = _mm_setr_epi32(1, 0, 0, 0);
				const __m128i _c_1 = _mm_setr_epi32(0, 1, 0, 0);
				const __m128i _c_2 = _mm_setr_epi32(0, 0, 1, 0);
				const __m128i _c_3 = _mm_setr_epi32(0, 0, 0, 1);

				// count the elements
				for (uint32_t j = 0; j < arrSize; j++) {


					uint32_t currItem = in[j];
					uint32_t index_0 = (currItem) & 0xff;
					uint32_t index_1 = (currItem >> 8) & 0xff;
					uint32_t index_2 = (currItem >> 16) & 0xff;
					uint32_t index_3 = (currItem >> 24) & 0xff;

					_mm_i32_(counting[index_0], 0)++;
					_mm_i32_(counting[index_1], 1)++;
					_mm_i32_(counting[index_2], 2)++;
					_mm_i32_(counting[index_3], 3)++;

					/*counting[index_0] = _mm_add_epi32(counting[index_0], _c_0);
					counting[index_1] = _mm_add_epi32(counting[index_1], _c_1);
					counting[index_2] = _mm_add_epi32(counting[index_2], _c_2);
					counting[index_3] = _mm_add_epi32(counting[index_3], _c_3);*/
				}

				//compute offsets
				__m128i acc = counting[0];

				counting[0] = _mm_set1_epi32(0);

				for (uint32_t j = 1; j < 256; j++) {
					__m128i tmp = counting[j];
					counting[j] = acc;
					acc = _mm_add_epi32(acc, tmp);
				}

				//shift = 0
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						uint32_t currItem = in[j];
						uint32_t bucket_index = ((currItem) & 0xff);

						uint32_t out_index = _mm_i32_(counting[bucket_index], 0);
						counting[bucket_index] = _mm_add_epi32(counting[bucket_index], _c_0);

						out[out_index] = currItem;
					}

					//swap out, in
					uint32_t* tmp = in;
					in = out;
					out = tmp;
				}

				//shift = 8
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						uint32_t currItem = in[j];
						uint32_t bucket_index = ((currItem >> 8) & 0xff);

						uint32_t out_index = _mm_i32_(counting[bucket_index], 1);

						counting[bucket_index] = _mm_add_epi32(counting[bucket_index], _c_1);

						out[out_index] = currItem;
					}

					//swap out, in
					uint32_t* tmp = in;
					in = out;
					out = tmp;
				}

				//shift = 16
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						uint32_t currItem = in[j];
						uint32_t bucket_index = ((currItem >> 16) & 0xff);

						uint32_t out_index = _mm_i32_(counting[bucket_index], 2);

						counting[bucket_index] = _mm_add_epi32(counting[bucket_index], _c_2);

						out[out_index] = currItem;
					}

					//swap out, in
					uint32_t* tmp = in;
					in = out;
					out = tmp;
				}

				//shift = 24
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						uint32_t currItem = in[j];
						uint32_t bucket_index = ((currItem >> 24) & 0xff);

						uint32_t out_index = _mm_i32_(counting[bucket_index], 3);

						counting[bucket_index] = _mm_add_epi32(counting[bucket_index], _c_3);

						out[out_index] = currItem;
					}

					//swap out, in
					uint32_t* tmp = in;
					in = out;
					out = tmp;
				}

			}

			if (tmp_array == NULL)
				free_aligned(aux);

#else

			if (arrSize == 0)
				return;

			//Counting Sort
			uint32_t counting[256][4];
			uint32_t* aux;

			if (tmp_array == NULL)
				aux = (uint32_t*)malloc_aligned(arrSize * sizeof(uint32_t));
			else
				aux = tmp_array;

			uint32_t* in = _arr;
			uint32_t* out = aux;

			//for (int32_t i = 0; i < radix_num; i++) 
			{

				// Cleaning counters
				memset(&counting[0][0], 0, sizeof(uint32_t) * 256 * 4);

				// count the elements
				for (uint32_t j = 0; j < arrSize; j++) {
					uint32_t currItem = in[j];
					counting[(currItem) & 0xff][0]++;
					counting[(currItem >> 8) & 0xff][1]++;
					counting[(currItem >> 16) & 0xff][2]++;
					counting[(currItem >> 24) & 0xff][3]++;
				}

				//compute offsets
				//uint32_t acc[4];
				//memcpy(acc, counting, sizeof(uint32_t) << 2);
				//memset(counting, 0, sizeof(uint32_t) << 2);

				uint32_t acc[4] = { counting[0][0], counting[0][1], counting[0][2], counting[0][3] };

				counting[0][0] = 0;
				counting[0][1] = 0;
				counting[0][2] = 0;
				counting[0][3] = 0;

				//uint32_t tmp[4];

				for (uint32_t j = 1; j < 256; j++) {

					//memcpy(tmp, counting[j], sizeof(uint32_t) << 2);
					//memcpy(counting[j], acc, sizeof(uint32_t) << 2);

					uint32_t tmp[4] = { counting[j][0], counting[j][1], counting[j][2], counting[j][3] };

					counting[j][0] = acc[0];
					counting[j][1] = acc[1];
					counting[j][2] = acc[2];
					counting[j][3] = acc[3];

					acc[0] += tmp[0];
					acc[1] += tmp[1];
					acc[2] += tmp[2];
					acc[3] += tmp[3];
				}

				//for (uint32_t shift = 0; shift < 32; shift += 8) 
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						uint32_t currItem = in[j];
						uint32_t bucket_index = ((currItem >> 0) & 0xff);
						uint32_t out_index = counting[bucket_index][0];
						counting[bucket_index][0]++;
						out[out_index] = currItem;
					}

					//swap out, in
					uint32_t* tmp = in;
					in = out;
					out = tmp;
				}
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						uint32_t currItem = in[j];
						uint32_t bucket_index = ((currItem >> 8) & 0xff);
						uint32_t out_index = counting[bucket_index][1];
						counting[bucket_index][1]++;
						out[out_index] = currItem;
					}

					//swap out, in
					uint32_t* tmp = in;
					in = out;
					out = tmp;
				}
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						uint32_t currItem = in[j];
						uint32_t bucket_index = ((currItem >> 16) & 0xff);
						uint32_t out_index = counting[bucket_index][2];
						counting[bucket_index][2]++;
						out[out_index] = currItem;
					}

					//swap out, in
					uint32_t* tmp = in;
					in = out;
					out = tmp;
				}
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						uint32_t currItem = in[j];
						uint32_t bucket_index = ((currItem >> 24) & 0xff);
						uint32_t out_index = counting[bucket_index][3];
						counting[bucket_index][3]++;
						out[out_index] = currItem;
					}

					//swap out, in
					uint32_t* tmp = in;
					in = out;
					out = tmp;
				}
			}

			if (tmp_array == NULL)
				free_aligned(aux);

#endif

		}


		void radix_counting_sort_unsigned_index(IndexUInt32* _arr, uint32_t arrSize, IndexUInt32* tmp_array) {


			////
			//// OLD ALGORITHM - very fast, but can be better!!!
			//// 
			////

			//if (arrSize == 0)
			//	return;

			//// best performance on 8 bit of symbols...
			//const int32_t base_bits = 8;//1 byte
			//const int32_t max_bucket_symbols = 1 << base_bits;//256
			//const int32_t base_mask = max_bucket_symbols - 1;//0xff

			//const int32_t int_total_bytes = sizeof(uint32_t);
			//const uint32_t int_mask_bits = 0xffffffff;
			//const int32_t int_total_bits = int_total_bytes << 3;// *8;

			//const int32_t radix_num = (int_total_bits / base_bits) + ((int_total_bits % base_bits) ? 1 : 0);
			//const uint32_t last_radix_mask = int_mask_bits >> ((radix_num - 1) * base_bits);

			//int32_t shift = 0;

			////Counting Sort
			//counter_type counting[max_bucket_symbols];
			//IndexUInt32* aux;

			//if (tmp_array == NULL)
			//	aux = (IndexUInt32*)malloc_aligned(arrSize * sizeof(IndexUInt32));
			//else
			//	aux = tmp_array;

			//IndexUInt32* in = _arr;
			//IndexUInt32* out = aux;

			//for (int32_t i = 0; i < radix_num; i++) {

			//	// Cleaning counters
			//	memset(counting, 0, sizeof(counter_type) * max_bucket_symbols);

			//	// count the elements
			//	for (uint32_t j = 0; j < arrSize; j++) {
			//		const IndexUInt32& currItem = in[j];
			//		uint32_t bucket_index = (((uint32_t)currItem.toSort >> shift) & base_mask);//0xff
			//		counting[bucket_index]++;
			//	}

			//	//compute offsets
			//	counter_type acc = counting[0];
			//	counting[0] = 0;
			//	for (uint32_t j = 1; j < max_bucket_symbols; j++) {
			//		counter_type tmp = counting[j];
			//		counting[j] = acc;
			//		acc += tmp;
			//	}

			//	// place elements in the output array
			//	for (uint32_t j = 0; j < arrSize; j++) {
			//		const IndexUInt32& currItem = in[j];
			//		uint32_t bucket_index = (((uint32_t)currItem.toSort >> shift) & base_mask);//0xff
			//		counter_type out_index = counting[bucket_index];
			//		counting[bucket_index]++;
			//		out[out_index] = currItem;
			//	}

			//	//swap out, in
			//	IndexUInt32* tmp = in;
			//	in = out;
			//	out = tmp;

			//	// Change shift value for next iterations
			//	shift += base_bits;
			//}

			//if (in != _arr)
			//	memcpy(_arr, in, arrSize * sizeof(IndexUInt32));

			//if (tmp_array == NULL)
			//	free_aligned(aux);

			//*************************
			//* New Uber Optimized!!! *
			//*************************

#if defined(ARIBEIRO_SSE2)


			if (arrSize == 0)
				return;

			//Counting Sort
			__m128i counting[256];
			IndexUInt32* aux;

			if (tmp_array == NULL)
				aux = (IndexUInt32*)malloc_aligned(arrSize * sizeof(IndexUInt32));
			else
				aux = tmp_array;

			IndexUInt32* in = _arr;
			IndexUInt32* out = aux;

			//for (int32_t i = 0; i < radix_num; i++) 
			{

				// Cleaning counters
				memset(counting, 0, sizeof(__m128i) * 256);


				const __m128i _c_0 = _mm_setr_epi32(1, 0, 0, 0);
				const __m128i _c_1 = _mm_setr_epi32(0, 1, 0, 0);
				const __m128i _c_2 = _mm_setr_epi32(0, 0, 1, 0);
				const __m128i _c_3 = _mm_setr_epi32(0, 0, 0, 1);

				// count the elements
				for (uint32_t j = 0; j < arrSize; j++) {


					uint32_t currItem = in[j].toSort;
					uint32_t index_0 = (currItem) & 0xff;
					uint32_t index_1 = (currItem >> 8) & 0xff;
					uint32_t index_2 = (currItem >> 16) & 0xff;
					uint32_t index_3 = (currItem >> 24) & 0xff;

					_mm_i32_(counting[index_0], 0)++;
					_mm_i32_(counting[index_1], 1)++;
					_mm_i32_(counting[index_2], 2)++;
					_mm_i32_(counting[index_3], 3)++;

					/*counting[index_0] = _mm_add_epi32(counting[index_0], _c_0);
					counting[index_1] = _mm_add_epi32(counting[index_1], _c_1);
					counting[index_2] = _mm_add_epi32(counting[index_2], _c_2);
					counting[index_3] = _mm_add_epi32(counting[index_3], _c_3);*/
				}

				//compute offsets
				__m128i acc = counting[0];

				counting[0] = _mm_set1_epi32(0);

				for (uint32_t j = 1; j < 256; j++) {
					__m128i tmp = counting[j];
					counting[j] = acc;
					acc = _mm_add_epi32(acc, tmp);
				}

				//shift = 0
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						const IndexUInt32& currItem = in[j];
						uint32_t bucket_index = ((currItem.toSort) & 0xff);

						uint32_t out_index = _mm_i32_(counting[bucket_index], 0);
						counting[bucket_index] = _mm_add_epi32(counting[bucket_index], _c_0);

						out[out_index] = currItem;
					}

					//swap out, in
					IndexUInt32* tmp = in;
					in = out;
					out = tmp;
				}

				//shift = 8
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						const IndexUInt32& currItem = in[j];
						uint32_t bucket_index = ((currItem.toSort >> 8) & 0xff);

						uint32_t out_index = _mm_i32_(counting[bucket_index], 1);

						counting[bucket_index] = _mm_add_epi32(counting[bucket_index], _c_1);

						out[out_index] = currItem;
					}

					//swap out, in
					IndexUInt32* tmp = in;
					in = out;
					out = tmp;
				}

				//shift = 16
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						const IndexUInt32& currItem = in[j];
						uint32_t bucket_index = ((currItem.toSort >> 16) & 0xff);

						uint32_t out_index = _mm_i32_(counting[bucket_index], 2);

						counting[bucket_index] = _mm_add_epi32(counting[bucket_index], _c_2);

						out[out_index] = currItem;
					}

					//swap out, in
					IndexUInt32* tmp = in;
					in = out;
					out = tmp;
				}

				//shift = 24
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						const IndexUInt32& currItem = in[j];
						uint32_t bucket_index = ((currItem.toSort >> 24) & 0xff);

						uint32_t out_index = _mm_i32_(counting[bucket_index], 3);

						counting[bucket_index] = _mm_add_epi32(counting[bucket_index], _c_3);

						out[out_index] = currItem;
					}

					//swap out, in
					IndexUInt32* tmp = in;
					in = out;
					out = tmp;
				}

			}

			if (tmp_array == NULL)
				free_aligned(aux);

#else

			if (arrSize == 0)
				return;

			//Counting Sort
			uint32_t counting[256][4];
			IndexUInt32* aux;

			if (tmp_array == NULL)
				aux = (IndexUInt32*)malloc_aligned(arrSize * sizeof(IndexUInt32));
			else
				aux = tmp_array;

			IndexUInt32* in = _arr;
			IndexUInt32* out = aux;

			//for (int32_t i = 0; i < radix_num; i++) 
			{

				// Cleaning counters
				memset(&counting[0][0], 0, sizeof(uint32_t) * 256 * 4);

				// count the elements
				for (uint32_t j = 0; j < arrSize; j++) {
					uint32_t currItem = in[j].toSort;
					counting[(currItem) & 0xff][0]++;
					counting[(currItem >> 8) & 0xff][1]++;
					counting[(currItem >> 16) & 0xff][2]++;
					counting[(currItem >> 24) & 0xff][3]++;
				}

				//compute offsets
				//uint32_t acc[4];
				//memcpy(acc, counting, sizeof(uint32_t) << 2);
				//memset(counting, 0, sizeof(uint32_t) << 2);

				uint32_t acc[4] = { counting[0][0], counting[0][1], counting[0][2], counting[0][3] };

				counting[0][0] = 0;
				counting[0][1] = 0;
				counting[0][2] = 0;
				counting[0][3] = 0;

				//uint32_t tmp[4];

				for (uint32_t j = 1; j < 256; j++) {

					//memcpy(tmp, counting[j], sizeof(uint32_t) << 2);
					//memcpy(counting[j], acc, sizeof(uint32_t) << 2);

					uint32_t tmp[4] = { counting[j][0], counting[j][1], counting[j][2], counting[j][3] };

					counting[j][0] = acc[0];
					counting[j][1] = acc[1];
					counting[j][2] = acc[2];
					counting[j][3] = acc[3];

					acc[0] += tmp[0];
					acc[1] += tmp[1];
					acc[2] += tmp[2];
					acc[3] += tmp[3];
				}

				//for (uint32_t shift = 0; shift < 32; shift += 8) 
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						const IndexUInt32& currItem = in[j];
						uint32_t bucket_index = ((currItem.toSort >> 0) & 0xff);
						uint32_t out_index = counting[bucket_index][0];
						counting[bucket_index][0]++;
						out[out_index] = currItem;
					}

					//swap out, in
					IndexUInt32* tmp = in;
					in = out;
					out = tmp;
				}
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						const IndexUInt32& currItem = in[j];
						uint32_t bucket_index = ((currItem.toSort >> 8) & 0xff);
						uint32_t out_index = counting[bucket_index][1];
						counting[bucket_index][1]++;
						out[out_index] = currItem;
					}

					//swap out, in
					IndexUInt32* tmp = in;
					in = out;
					out = tmp;
				}
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						const IndexUInt32& currItem = in[j];
						uint32_t bucket_index = ((currItem.toSort >> 16) & 0xff);
						uint32_t out_index = counting[bucket_index][2];
						counting[bucket_index][2]++;
						out[out_index] = currItem;
					}

					//swap out, in
					IndexUInt32* tmp = in;
					in = out;
					out = tmp;
				}
				{
					// place elements in the output array
					for (uint32_t j = 0; j < arrSize; j++) {
						const IndexUInt32& currItem = in[j];
						uint32_t bucket_index = ((currItem.toSort >> 24) & 0xff);
						uint32_t out_index = counting[bucket_index][3];
						counting[bucket_index][3]++;
						out[out_index] = currItem;
					}

					//swap out, in
					IndexUInt32* tmp = in;
					in = out;
					out = tmp;
				}
			}

			if (tmp_array == NULL)
				free_aligned(aux);

#endif

		}


		void radix_bucket_sort_signed(int32_t* arr, uint32_t arrSize) {

			if (arrSize == 0)
				return;

			const int32_t base_bits = 11;//1 byte
			const int32_t max_bucket_symbols = 1 << base_bits;//256
			const int32_t base_mask = max_bucket_symbols - 1;//0xff

			const int32_t int_total_bytes = sizeof(int32_t);
			const uint32_t int_mask_bits = 0xffffffff;
			const int32_t int_total_bits = int_total_bytes << 3;// *8;

			const int32_t radix_num = (int_total_bits / base_bits) + ((int_total_bits % base_bits) ? 1 : 0);
			const uint32_t last_radix_mask = int_mask_bits >> ((radix_num - 1) * base_bits);

			int32_t shift = 0;

			//bucket sort groups
			std::vector<int32_t>* bucket = new std::vector<int32_t>[max_bucket_symbols];//groups[256]

			for (int32_t i = 0; i < radix_num; i++) {

				// Cleaning groups
				for (int32_t j = 0; j < max_bucket_symbols; j++) {
					bucket[j].clear();
				}

				// Splitting items into radix groups
				for (uint32_t j = 0; j < arrSize; j++) {
					int32_t currItem = arr[j];
					uint32_t bucket_index = (((uint32_t)currItem >> shift) & base_mask);//0xff
					bucket[bucket_index].push_back(currItem);
				}

				// Copying sorted by radix items back into original array
				uint32_t arrPos = 0;

				// Treat the last radix with sign bit specially
				// Output signed groups (128..256 = -128..-1) first
				// Other groups afterwards. No performance penalty, as compared to flipping sign bit
				// via (($currItem ^ 0x8000000000000000) >> $shift) & 0xFF)
				if (i == radix_num - 1) {
					uint32_t max = last_radix_mask;
					uint32_t max_half = (max >> 1) + 1;
					for (uint32_t j = max_half; j <= max; j++) {
						std::vector<int32_t>& bucket_j = bucket[j];
						for (size_t k = 0; k < bucket_j.size(); k++) {
							arr[arrPos++] = bucket_j[k];
						}
					}

					for (uint32_t j = 0; j < max_half; j++) {
						std::vector<int32_t>& bucket_j = bucket[j];
						for (size_t k = 0; k < bucket_j.size(); k++) {
							arr[arrPos++] = bucket_j[k];
						}
					}
				}
				else {
					for (uint32_t j = 0; j < max_bucket_symbols; j++) {
						std::vector<int32_t>& bucket_j = bucket[j];
						for (size_t k = 0; k < bucket_j.size(); k++) {
							arr[arrPos++] = bucket_j[k];
						}
					}
				}

				// Change shift value for next iterations
				shift += base_bits;
			}

			delete[]bucket;
		}

		void radix_bucket_sort_signed_index(IndexInt32* arr, uint32_t arrSize) {

			if (arrSize == 0)
				return;

			const int32_t base_bits = 11;//1 byte
			const int32_t max_bucket_symbols = 1 << base_bits;//256
			const int32_t base_mask = max_bucket_symbols - 1;//0xff

			const int32_t int_total_bytes = sizeof(int32_t);
			const uint32_t int_mask_bits = 0xffffffff;
			const int32_t int_total_bits = int_total_bytes << 3;// *8;

			const int32_t radix_num = (int_total_bits / base_bits) + ((int_total_bits % base_bits) ? 1 : 0);
			const uint32_t last_radix_mask = int_mask_bits >> ((radix_num - 1) * base_bits);

			int32_t shift = 0;

			//bucket sort groups
			std::vector<IndexInt32>* bucket = new std::vector<IndexInt32>[max_bucket_symbols];//groups[256]

			for (int32_t i = 0; i < radix_num; i++) {

				// Cleaning groups
				for (int32_t j = 0; j < max_bucket_symbols; j++) {
					bucket[j].clear();
				}

				// Splitting items into radix groups
				for (uint32_t j = 0; j < arrSize; j++) {
					const IndexInt32& currItem = arr[j];
					uint32_t bucket_index = (((uint32_t)currItem.toSort >> shift) & base_mask);//0xff
					bucket[bucket_index].push_back(currItem);
				}

				// Copying sorted by radix items back into original array
				uint32_t arrPos = 0;

				// Treat the last radix with sign bit specially
				// Output signed groups (128..256 = -128..-1) first
				// Other groups afterwards. No performance penalty, as compared to flipping sign bit
				// via (($currItem ^ 0x8000000000000000) >> $shift) & 0xFF)
				if (i == radix_num - 1) {
					uint32_t max = last_radix_mask;
					uint32_t max_half = (max >> 1) + 1;
					for (uint32_t j = max_half; j <= max; j++) {
						std::vector<IndexInt32>& bucket_j = bucket[j];
						for (size_t k = 0; k < bucket_j.size(); k++) {
							arr[arrPos++] = bucket_j[k];
						}
					}

					for (uint32_t j = 0; j < max_half; j++) {
						std::vector<IndexInt32>& bucket_j = bucket[j];
						for (size_t k = 0; k < bucket_j.size(); k++) {
							arr[arrPos++] = bucket_j[k];
						}
					}
				}
				else {
					for (uint32_t j = 0; j < max_bucket_symbols; j++) {
						std::vector<IndexInt32>& bucket_j = bucket[j];
						for (size_t k = 0; k < bucket_j.size(); k++) {
							arr[arrPos++] = bucket_j[k];
						}
					}
				}

				// Change shift value for next iterations
				shift += base_bits;
			}

			delete[]bucket;
		}


		void radix_bucket_sort_unsigned(uint32_t* arr, uint32_t arrSize) {

			if (arrSize == 0)
				return;

			const int32_t base_bits = 11;//1 byte
			const int32_t max_bucket_symbols = 1 << base_bits;//256
			const int32_t base_mask = max_bucket_symbols - 1;//0xff

			const int32_t int_total_bytes = sizeof(int32_t);
			const uint32_t int_mask_bits = 0xffffffff;
			const int32_t int_total_bits = int_total_bytes << 3;// *8;

			const int32_t radix_num = (int_total_bits / base_bits) + ((int_total_bits % base_bits) ? 1 : 0);
			const uint32_t last_radix_mask = int_mask_bits >> ((radix_num - 1) * base_bits);

			int32_t shift = 0;

			//bucket sort groups
			std::vector<uint32_t>* bucket = new std::vector<uint32_t>[max_bucket_symbols];//groups[256]

			for (int32_t i = 0; i < radix_num; i++) {

				// Cleaning groups
				for (int32_t j = 0; j < max_bucket_symbols; j++) {
					bucket[j].clear();
				}

				// Splitting items into radix groups
				for (uint32_t j = 0; j < arrSize; j++) {
					uint32_t currItem = arr[j];
					uint32_t bucket_index = (((uint32_t)currItem >> shift) & base_mask);//0xff
					bucket[bucket_index].push_back(currItem);
				}

				// Copying sorted by radix items back into original array
				uint32_t arrPos = 0;

				for (uint32_t j = 0; j < max_bucket_symbols; j++) {
					std::vector<uint32_t>& bucket_j = bucket[j];
					for (size_t k = 0; k < bucket_j.size(); k++) {
						arr[arrPos++] = bucket_j[k];
					}
				}

				// Change shift value for next iterations
				shift += base_bits;
			}

			delete[]bucket;
		}

		void radix_bucket_sort_unsigned_index(IndexUInt32* arr, uint32_t arrSize) {

			if (arrSize == 0)
				return;

			const int32_t base_bits = 11;//1 byte
			const int32_t max_bucket_symbols = 1 << base_bits;//256
			const int32_t base_mask = max_bucket_symbols - 1;//0xff

			const int32_t int_total_bytes = sizeof(int32_t);
			const uint32_t int_mask_bits = 0xffffffff;
			const int32_t int_total_bits = int_total_bytes << 3;// *8;

			const int32_t radix_num = (int_total_bits / base_bits) + ((int_total_bits % base_bits) ? 1 : 0);
			const uint32_t last_radix_mask = int_mask_bits >> ((radix_num - 1) * base_bits);

			int32_t shift = 0;

			//bucket sort groups
			std::vector<IndexUInt32>* bucket = new std::vector<IndexUInt32>[max_bucket_symbols];//groups[256]

			for (int32_t i = 0; i < radix_num; i++) {

				// Cleaning groups
				for (int32_t j = 0; j < max_bucket_symbols; j++) {
					bucket[j].clear();
				}

				// Splitting items into radix groups
				for (uint32_t j = 0; j < arrSize; j++) {
					const IndexUInt32& currItem = arr[j];
					uint32_t bucket_index = (((uint32_t)currItem.toSort >> shift) & base_mask);//0xff
					bucket[bucket_index].push_back(currItem);
				}

				// Copying sorted by radix items back into original array
				uint32_t arrPos = 0;

				for (uint32_t j = 0; j < max_bucket_symbols; j++) {
					std::vector<IndexUInt32>& bucket_j = bucket[j];
					for (size_t k = 0; k < bucket_j.size(); k++) {
						arr[arrPos++] = bucket_j[k];
					}
				}

				// Change shift value for next iterations
				shift += base_bits;
			}

			delete[]bucket;
		}

	}
}
