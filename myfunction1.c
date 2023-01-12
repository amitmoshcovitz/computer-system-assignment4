#include <stdlib.h>
#include <stdbool.h>
#include "myfunction1.h"
#include "showBMP.h"
#include <string.h>
#include <stdio.h>

#pragma GCC optimize("Ofast", "fast-math", "unroll-loops")
#pragma GCC target ("avx2")

//great way I found online to calculate minimux and maximun without any comparisons.
#define macro_max(x,y) (x ^ ((x ^ y) & -(x < y)))
#define macro_min(x,y) (y ^ ((x ^ y) & -(x < y)))


inline void rowBlurKernel_func(Image *image, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {

	//saving n and data localy in order to save time later.
	unsigned char *data = image->data;
	short local_n = n;

	//we will use (3 * n) and (3 * n * n) a lot so we will save their value. -- notice it stays the same.
	short n_3 = 3*local_n;
	int n_n_3 = n_3 * local_n;

	//creating a temporary array the size of image for the calculations.
	unsigned char image_copy[n_n_3];
	memcpy(image_copy, data, n_n_3);

	//saving the for conditions so we will not have to calculations them each iteration.
	short condition_outer = local_n - 1;
	short condition_inner = n_3 - 3;

	//Varibles for the for loop optimization.
	unsigned char *data_index_temp = data + 3 + n_3;
	unsigned char *copy_index_temp = image_copy + 3 + n_3;
	unsigned char *data_index;
	unsigned char *copy_index;

	register int i,j;
	for (i = 1; i < condition_outer; ++i) {
		data_index = data_index_temp;
		copy_index = copy_index_temp;
		for (j = 3; j < condition_inner; j += 3) {
			// assign kernel's result to pixel at [i,j]
			*data_index = (*(copy_index - 3)  + *(copy_index)* 2 + *(copy_index + 3)) / 4;
			*(data_index + 1) = (*(copy_index - 2) + *(copy_index + 1) * 2 + *(copy_index + 4)) / 4;
			*(data_index + 2) = (*(copy_index - 1) + *(copy_index + 2) * 2 + *(copy_index + 5)) / 4;

			data_index += 3;
			copy_index += 3;
		}
		data_index_temp += n_3;
		copy_index_temp += n_3;
	}
}

inline void rowSharpKernel_func(Image *image, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {
	//saving n and data localy in order to save time later.
	unsigned char *data = image->data;
	short local_n = n;

	//we will use (3 * n) and (3 * n * n) a lot so we will save their value. -- notice it stays the same.
	short n_3 = 3*local_n;
	int n_n_3 = n_3 * local_n;

	//creating a temporary array the size of image for the calculations.
	unsigned char image_copy[n_n_3];
	memcpy(image_copy, data, n_n_3);

	//saving the for conditions so we will not have to calculations them each iteration.
	short condition_outer = local_n - 1;
	short condition_inner = n_3 - 3;

	//Varibles for the for loop optimization.
	unsigned char *data_index_temp = data + 3 + n_3;
	unsigned char *copy_index_temp = image_copy + 3 + n_3;
	unsigned char *data_index;
	unsigned char *copy_index;

	register int i,j;
	for (i = 1; i < condition_outer; ++i) {
		data_index = data_index_temp;
		copy_index = copy_index_temp;
		for (j = 3; j < condition_inner; j += 3) {
			// assign kernel's result to pixel at [i,j]
			*data_index = macro_min(macro_max((*(copy_index) * 3 - *(copy_index - 3) - *(copy_index + 3)), 0), 255);
			*(data_index + 1) = macro_min(macro_max((*(copy_index + 1) * 3 - *(copy_index - 2) - *(copy_index + 4)), 0), 255);
			*(data_index + 2) = macro_min(macro_max((*(copy_index + 2) * 3 - *(copy_index - 1) - *(copy_index + 5)), 0), 255);

			data_index += 3;
			copy_index += 3;
		}
		data_index_temp += n_3;
		copy_index_temp += n_3;
	}	
}

inline void blurKernel_withou_filter_func(Image *image, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {

	//saving n and data localy in order to save time later.
	unsigned char *data = image->data;
	short local_n = n;

	//we will use (3 * n) and (3 * n * n) a lot so we will save their value. -- notice it stays the same.
	short n_3 = 3*local_n;
	int n_n_3 = n_3 * local_n;

	//creating a temporary array the size of image for the calculations.
	unsigned char image_copy[n_n_3];
	memcpy(image_copy, data, n_n_3);

	//saving the for conditions so we will not have to calculations them each iteration.
	short condition_outer = local_n - 1;
	short condition_inner = n_3 - 3;



	//Varibles for the for loop optimization.
	unsigned char *i_j_dim;
	unsigned char* start_3;
	unsigned char* middle_3; 
	unsigned char* end_3;
	unsigned char* index_value;
	unsigned char* index_value_part_one = data + n_3 + 3;
	unsigned char *i_j_dim_help = image_copy + n_3 + 3;

	register int i, j;
	for (i = 1 ; i < condition_outer; ++i) {
		index_value = index_value_part_one;
		i_j_dim = i_j_dim_help;
		for (j = 3 ; j < condition_inner; j += 3) {

			start_3 = i_j_dim - n_3;
			middle_3 = i_j_dim;
			end_3 = i_j_dim + n_3;

			// assign kernel's result to pixel at [i,j]
			*(index_value) = ((*(middle_3)) + (*(start_3 - 3)) + (*(start_3)) + (*(start_3 + 3)) + (*(middle_3 - 3)) + (*(middle_3 + 3)) + (*(end_3 - 3)) + (*(end_3)) + (*(end_3 + 3))) / 9;
			*(index_value + 1) = ((*(middle_3 + 1)) + (*(start_3 - 2)) + (*(start_3 + 1)) + (*(start_3 + 4)) + (*(middle_3 - 2)) + (*(middle_3 + 4))+ (*(end_3 - 2)) + (*(end_3 + 1)) + (*(end_3 + 4))) / 9;
			*(index_value + 2) = ((*(middle_3 + 2)) + (*(start_3 - 1)) + (*(start_3 + 2)) + (*(start_3 + 5)) + (*(middle_3 - 1)) + (*(middle_3 + 5))+ (*(end_3 - 1)) + (*(end_3 + 2)) + (*(end_3 + 5))) / 9;

			index_value += 3;
			i_j_dim += 3;
		}	
		index_value_part_one += n_3;
		i_j_dim_help += n_3;
	}
}

inline void blurKernel_with_filter_func(Image *image, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {

	//saving n and data localy in order to save time later.
	unsigned char *data = image->data;
	short local_n = n;

	//we will use (3 * n) and (3 * n * n) a lot so we will save their value. -- notice it stays the same.
	short n_3 = 3*local_n;
	int n_n_3 = n_3 * local_n;

	//creating a temporary array the size of image for the calculations.
	unsigned char image_copy[n_n_3];
	memcpy(image_copy, data, n_n_3);

	//saving the for conditions so we will not have to calculations them each iteration.
	short condition_outer = local_n - 1;
	short condition_inner = n_3 - 3;



	//Varibles for the for loop optimization.
	unsigned char* start_3;
	unsigned char* middle_3; 
	unsigned char* end_3;

	short condition = local_n - 1;
	int index_value_part_one = n_3 + 3;

	//Varibles




	int index_value;


	//Varibles for the minimux and maximum.
	short min_intensity = 766; // arbitrary value that is higher than maximum possible intensity, which is 255*3=765
	int max_intensity = -1; // arbitrary value that is lower than minimum possible intensity, which is 0
	short min_row, min_col, max_row, max_col;
	unsigned char *loop_pixel;

	int i_j_dim_help = local_n;
	unsigned char *m_i_3 = image_copy + n_3;
	register int i, j;
	for (i = 1 ; i < condition_outer; ++i) {
		index_value = index_value_part_one;
		for (j = 3 ; j < condition_inner; j += 3) {


			loop_pixel = m_i_3 - n_3 + j;
			unsigned char i_0_0_r = *(loop_pixel - 3); unsigned char i_0_1_r = *(loop_pixel); unsigned char i_0_2_r = *(loop_pixel + 3);
			unsigned char i_0_0_g = *(loop_pixel - 2); unsigned char i_0_1_g = *(loop_pixel + 1); unsigned char i_0_2_g = *(loop_pixel + 4);
			unsigned char i_0_0_b = *(loop_pixel - 1); unsigned char i_0_1_b = *(loop_pixel + 2); unsigned char i_0_2_b = *(loop_pixel + 5);
			loop_pixel = m_i_3 + j;
			unsigned char i_1_0_r = *(loop_pixel - 3); unsigned char i_1_1_r = *(loop_pixel); unsigned char i_1_2_r = *(loop_pixel + 3);
			unsigned char i_1_0_g = *(loop_pixel - 2); unsigned char i_1_1_g = *(loop_pixel + 1); unsigned char i_1_2_g = *(loop_pixel + 4);
			unsigned char i_1_0_b = *(loop_pixel - 1); unsigned char i_1_1_b = *(loop_pixel + 2); unsigned char i_1_2_b = *(loop_pixel + 5);
			loop_pixel = m_i_3 + n_3 + j;
			unsigned char i_2_0_r = *(loop_pixel - 3); unsigned char i_2_1_r = *(loop_pixel); unsigned char i_2_2_r = *(loop_pixel + 3);
			unsigned char i_2_0_g = *(loop_pixel - 2); unsigned char i_2_1_g = *(loop_pixel + 1); unsigned char i_2_2_g = *(loop_pixel + 4);
			unsigned char i_2_0_b = *(loop_pixel - 1); unsigned char i_2_1_b = *(loop_pixel + 2); unsigned char i_2_2_b = *(loop_pixel + 5);


			register int minVal_r;			register int minVal_g;			register int minVal_b;
			register int maxVal_r;			register int maxVal_g;			register int maxVal_b;
			short sum = i_0_0_r + i_0_0_g + i_0_0_b;
			min_intensity = sum;
			minVal_r = i_0_0_r;
			minVal_g = i_0_0_g;
			minVal_b = i_0_0_b;
			max_intensity = sum;
			maxVal_r = i_0_0_r;
			maxVal_g = i_0_0_g;
			maxVal_b = i_0_0_b;
			sum = i_0_1_r + i_0_1_g + i_0_1_b;
			if (sum <= min_intensity) {
				min_intensity = sum;
				minVal_r = i_0_1_r;
				minVal_g = i_0_1_g;
				minVal_b = i_0_1_b;
			}
			if (sum  > max_intensity) {
				max_intensity = sum;
				maxVal_r = i_0_1_r;
				maxVal_g = i_0_1_g;
				maxVal_b = i_0_1_b;
			}
			sum = i_0_2_r + i_0_2_g + i_0_2_b;
			if (sum <= min_intensity) {
				min_intensity = sum;
				minVal_r = i_0_2_r;
				minVal_g = i_0_2_g;
				minVal_b = i_0_2_b;
			}
			if (sum  > max_intensity) {
				max_intensity = sum;
				maxVal_r = i_0_2_r;
				maxVal_g = i_0_2_g;
				maxVal_b = i_0_2_b;
			}
			sum = i_1_0_r + i_1_0_g + i_1_0_b;
			if (sum <= min_intensity) {
				min_intensity = sum;
				minVal_r = i_1_0_r;
				minVal_g = i_1_0_g;
				minVal_b = i_1_0_b;
			}
			if (sum  > max_intensity) {
				max_intensity = sum;;
				maxVal_r = i_1_0_r;
				maxVal_g = i_1_0_g;
				maxVal_b = i_1_0_b;
			}
			sum = i_1_1_r + i_1_1_g + i_1_1_b;
			if (sum <= min_intensity) {
				min_intensity = sum;
				minVal_r = i_1_1_r;
				minVal_g = i_1_1_g;
				minVal_b = i_1_1_b;
			}
			if (sum  > max_intensity) {
				max_intensity = sum;
				maxVal_r = i_1_1_r;
				maxVal_g = i_1_1_g;
				maxVal_b = i_1_1_b;
			}
			sum = i_1_2_r + i_1_2_g + i_1_2_b;
			if (sum <= min_intensity) {
				min_intensity = sum;
				minVal_r = i_1_2_r;
				minVal_g = i_1_2_g;
				minVal_b = i_1_2_b;
			}
			if (sum  > max_intensity) {
				max_intensity = sum;
				maxVal_r = i_1_2_r;
				maxVal_g = i_1_2_g;
				maxVal_b = i_1_2_b;
			}
			sum = i_2_0_r + i_2_0_g + i_2_0_b;
			if (sum <= min_intensity) {
				min_intensity = sum;
				minVal_r = i_2_0_r;
				minVal_g = i_2_0_g;
				minVal_b = i_2_0_b;
			}
			if (sum  > max_intensity) {
				max_intensity = sum;
				maxVal_r = i_2_0_r;
				maxVal_g = i_2_0_g;
				maxVal_b = i_2_0_b;
			}
			sum = i_2_1_r + i_2_1_g + i_2_1_b;
			if (sum <= min_intensity) {
				min_intensity = sum;
				minVal_r = i_2_1_r;
				minVal_g = i_2_1_g;
				minVal_b = i_2_1_b;
			}
			if (sum  > max_intensity) {
				max_intensity = sum;
				maxVal_r = i_2_1_r;
				maxVal_g = i_2_1_g;
				maxVal_b = i_2_1_b;
			}
			sum = i_2_2_r + i_2_2_g + i_2_2_b;
			if (sum <= min_intensity) {
				min_intensity = sum;
				minVal_r = i_2_2_r;
				minVal_g = i_2_2_g;
				minVal_b = i_2_2_b;
			}
			if (sum  > max_intensity) {
				max_intensity = sum;
				maxVal_r = i_2_2_r;
				maxVal_g = i_2_2_g;
				maxVal_b = i_2_2_b;
			}

			// assign kernel's result to pixel at [i,j]
			*(data + index_value) = (i_0_0_r + i_0_1_r + i_0_2_r + i_1_0_r + i_1_1_r + i_1_2_r + i_2_0_r + i_2_1_r + i_2_2_r - (minVal_r + maxVal_r)) / 7;
			*(data + index_value + 1) = (i_0_0_g + i_0_1_g + i_0_2_g + i_1_0_g + i_1_1_g + i_1_2_g + i_2_0_g + i_2_1_g + i_2_2_g - (minVal_g + maxVal_g)) / 7;
			*(data + index_value + 2) = (i_0_0_b + i_0_1_b + i_0_2_b + i_1_0_b + i_1_1_b + i_1_2_b + i_2_0_b + i_2_1_b + i_2_2_b - (minVal_b + maxVal_b)) / 7;

			index_value += 3;
		}	
		index_value_part_one += n_3;
		m_i_3 += n_3;
	}
}

inline void sharpKernel_func(Image *image, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {
		
	//saving n and data localy in order to save time later.
	unsigned char *data = image->data;
	short local_n = n;

	//we will use (3 * n) and (3 * n * n) a lot so we will save their value. -- notice it stays the same.
	short n_3 = 3*local_n;
	int n_n_3 = n_3 * local_n;

	//creating a temporary array the size of image for the calculations.
	unsigned char image_copy[n_n_3];
	memcpy(image_copy, data, n_n_3);

	//saving the for conditions so we will not have to calculations them each iteration.
	short condition_outer = local_n - 1;
	short condition_inner = n_3 - 3;



	//Varibles for the for loop optimization.
	unsigned char *i_j_dim;
	unsigned char* start_3;
	unsigned char* middle_3; 
	unsigned char* end_3;
	unsigned char* index_value;
	unsigned char* index_value_part_one = data + n_3 + 3;
	unsigned char *i_j_dim_help = image_copy + n_3 + 3;

	register int i, j;
	for (i = 1 ; i < condition_outer; ++i) {
		index_value = index_value_part_one;
		i_j_dim = i_j_dim_help;
		for (j = 3 ; j < condition_inner; j += 3) {

			start_3 = i_j_dim - n_3;
			middle_3 = i_j_dim;
			end_3 = i_j_dim + n_3;

			// assign kernel's result to pixel at [i,j]
			*(index_value) = (macro_min(macro_max((*(middle_3)) * 9 - ((*(start_3 - 3)) + (*(start_3)) + (*(start_3 + 3)) + (*(middle_3 - 3)) + (*(middle_3 + 3)) + (*(end_3 - 3)) + (*(end_3)) + (*(end_3 + 3))), 0), 255));
			*(index_value + 1) = (macro_min(macro_max((*(middle_3 + 1)) * 9 - ((*(start_3 - 2)) + (*(start_3 + 1)) + (*(start_3 + 4)) + (*(middle_3 - 2)) + (*(middle_3 + 4)) + (*(end_3 - 2)) + (*(end_3 + 1)) + (*(end_3 + 4))), 0), 255));
			*(index_value + 2) = (macro_min(macro_max((*(middle_3 + 2)) * 9 - ((*(start_3 - 1)) + (*(start_3 + 2)) + (*(start_3 + 5)) + (*(middle_3 - 1)) + (*(middle_3 + 5)) + (*(end_3 - 1)) + (*(end_3 + 2)) + (*(end_3 + 5))), 0), 255));
			index_value += 3;
			i_j_dim += 3;
		}
		index_value_part_one += n_3;
		i_j_dim_help += n_3;
	}
}





//TODO: remove parameters.
//TODO: check if I can remove more min-max.
//TODO: Change the name of min-max functions!!!!