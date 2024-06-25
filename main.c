#include <stdio.h> // Includes input/output functions (files and streams)
#include <stdlib.h> // Includes memory allocation functions
#include <math.h> // Includes mathematical functions and Pi constant

/*
Define N here so that it can be used in h2_result struct
Other constants are defined here so that they are immutable
Using constants instead of variables to store these constants is useful for
readability, optimization and preventing unintentional modifications during
execution time
The type of a constant usually depends on the context in which it is used,
even if it can be forced to be of a specific type
*/
#define N 100
#define N2 200
#define N3 4
#define T 2 * M_PI

// Definition of structures for complex numbers handling
typedef struct CN
{
  // Stores the real and imaginary part of a complex number.
  double re;
  double im;
} CN;

typedef struct h2_result
{
  // Stores the result of a function that returns a complex array.
  CN result[N];
} h2_result;

// Definition of functions, the whole functions are after main() for readability
double magnitude(CN c_num);

int in_array(int arr[], int size_arr, int target, int cond_type);

// The ** instead of the * is used to be able to change the array
void DFT(CN *input, CN **output, int num);

void IDFT(CN *input, CN **output, int num, int els[], int size_arr, 
  int type_con);

// a. Implement the mathematical functions h1 and h2
void h1(double step, CN result[]);

h2_result h2(double step);

void write_to_file(const char* filename, double step, CN *h_values, int num);

void find_max_magnitudes_indices(double magnitudes[], int indices[],
  int num_elements, int num_max);

int main()
{
  /*
  b. Sample these functions N = 100 times each over a time period T = 2 * PI
  and write the results to text files for plotting
  Example of using arrays to handle function results
  Since arrays cannot be returned, a void function is used instead
  In this case, the size of the result is set from the start, so static
  arrays are a natural data structure to choose
  Static arrays use the stack, the most accesible memory, used for temporary
  variables
  */
  double step = T / N;
  CN h1_values_static[N];
  h1(step, h1_values_static);
  
  /*
  Static arrays do not need to be freed like dynamic arrays
  Arrays cannot in general be returned from a function because they are local
  variables
  An indirect way of returning a static array is using a struct object of
  array type, without directly using dynamic memory
  */
  h2_result h2_values_custom = h2(step);
  
  /*
  Dynamic memory uses the heap, a more common choice when data has to persist
  and manual memory management is preferred. The arrays containing results are
  especially important to keep
  Static and custom type arrays can be converted to dynamic arrays by iterating
  over the elements. This is in general useful when the same function has to be
  applied to arguments of different types. This can be done with a for or a
  while loop:
  */
  CN* h1_values = (CN*)malloc(N * sizeof(CN));
  
  if (h1_values == NULL)
  {
    printf("Error: Memory allocation failed for h1_values.\n");
    
    return -1;
  }
  
  // Copy values from the existing array to the dynamically allocated memory
  for (int i = 0; i < N; i++)
  {
    h1_values[i] = h1_values_static[i];
  }
  
  // Allocate memory for h2_values dynamically
  CN* h2_values = (CN*)malloc(N * sizeof(CN));
  
  if (h2_values == NULL)
  {
    printf("Error: Memory allocation failed for h2_values.\n");
    
    return -1;
  }
  
  // Copy values from the h2 struct to the dynamically allocated memory
  for (int i = 0; i < N; i++)
  {
    h2_values[i] = h2_values_custom.result[i];
  }
  
  /*
  A function is used to write the results to files. This function will also
  be useful later to write the results of the inverse Fourier transform
  */
  write_to_file("h1.txt", step, h1_values, N);
  write_to_file("h2.txt", step, h2_values, N);

  /*
  d. Use discrete Fourier transform to obtain H1(w) and H2(w) from h1(t) and
  h2(t), respectively
  Dynamic arrays can be returned from a function, but it is not good practice
  not to free the space so the memory is allocated in main()
  If the memory is allocated inside a function, the memory can only be freed
  inside, before the return statement, and then the result cannot be returned
  From this point to the end of the program, dynamic memory is used to show
  different cases and ways in which dynamic memory can be used
  */
  CN* H1 = (CN*)malloc(N * sizeof(CN));
  CN* H2 = (CN*)malloc(N * sizeof(CN));

  if (H1 == NULL || H2 == NULL)
  {
    printf("Error: Memory allocation failed for H1 or H2.\n");
    
    return -1;
  }

  DFT(h1_values, &H1, N);
  DFT(h2_values, &H2, N);
  
  free(h1_values);
  free(h2_values);
  
  /*
  e. Print the results for H1(w) and H2(w)
  A table of the Fourier transforms is printed with some logic for a nicer
  formatting of complex numbers
  */
  printf("-----------------------------------------------------------\n");
  printf("|             H1             |             H2             |\n");
  printf("-----------------------------------------------------------\n");
  
  for (int i = 0; i < N; i++)
  {
    if (H1[i].im < 0 && H2[i].im < 0)
    {
      printf("|%11f - %11f i |%11f - %11f i |\n", H1[i].re, -H1[i].im, H2[i].re,
        -H2[i].im);
    } else if (H1[i].im < 0 && H2[i].im >= 0)
    {
      printf("|%11f - %11f i |%11f + %11f i |\n", H1[i].re, -H1[i].im, H2[i].re,
        H2[i].im);
    } else if (H1[i].im >= 0 && H2[i].im < 0)
    {
      printf("|%11f + %11f i |%11f - %11f i |\n", H1[i].re, H1[i].im, H2[i].re,
        -H2[i].im);
    } else
    {
      printf("|%11f + %11f i |%11f + %11f i |\n", H1[i].re, H1[i].im, H2[i].re,
        H2[i].im);
    }
  }
  
  printf("-----------------------------------------------------------\n");
  
  /*
  f. Apply an Inverse Fourier transform to H1(w) and H2(w) to obtain h'1 and
  h'2 respectively. Skipping n = 1 for H1(w) and n = 0 for H2(w)
  */
  CN* h1_prime = (CN*)malloc(N * sizeof(CN));
  CN* h2_prime = (CN*)malloc(N * sizeof(CN));
  int skip_1[] = {1};
  int skip_2[] = {0};
  
  IDFT(H1, &h1_prime, N, skip_1, 1, 0);
  IDFT(H2, &h2_prime, N, skip_2, 1, 0);

  free(H1);
  free(H2);
  
  /*
  g. Write the outcomes of the inverse Fourier transforms to text files for
  plotting
  The file outputting procedure with h1 and h2 is repeated with h1' and h2',
  but the variable types are different because of the demonstration made using
  different result handling methods/variable types
  */
  write_to_file("h1_prime.txt", step, h1_prime, N);
  write_to_file("h2_prime.txt", step, h2_prime, N);
  
  free(h1_prime);
  free(h2_prime);
  
  // Get data from file:
  FILE *fp3;
  CN* h3_values = (CN*)malloc(N2 * sizeof(CN));
  
  if (h3_values == NULL)
  {
    printf("Error: Memory allocation failed for h3_values.\n");
    
    return -1;
  }

  // i. Load the sampling data from the file "h3.txt" into memory
  const char* filename_h3 = "h3.txt";
  fp3 = fopen(filename_h3, "r");
  double times[N2];
  
  if (fp3 == NULL)
  {
    printf("Error: Unable to open file %s for reading.\n", filename_h3);
    
    return -1;
  }
  
  for (int i = 0; i < N2; i++)
  {
    // * in format to discard first element of row
    fscanf(fp3,"%*i,%lf,%lf,%lf",&times[i], &h3_values[i].re, &h3_values[i].im);
  }
  
  fclose(fp3);
  
  // j. Apply a discrete Fourier transform to h3 to obtain H3
  CN* H3 = (CN*)malloc(N2 * sizeof(CN));
  
  DFT(h3_values, &H3, N2);
  
  free(h3_values);
  
  // Calculate magnitudes of H3 values
  double magnitudes_H3[N2];
  
  for (int i = 0; i < N2; i++)
  {
    magnitudes_H3[i] = magnitude(H3[i]);
  }
  
  // Initialise array of indices of elements with maximum magnitudes
  int max_indices[N3];
  
  // Get indices of elements with maximum magnitude
  find_max_magnitudes_indices(magnitudes_H3, max_indices, N2, N3);
  
  /*
  k. Apply an inverse discrete Fourier transform to H3 to obtain h'3,
  including only the four terms out of the N terms for H3 with the largest
  amplitude
  */
  CN* h3_prime = (CN*)malloc(N2 * sizeof(CN));
  
  IDFT(H3, &h3_prime, N2, max_indices, N3, 1);
  
  // l. Write the outcome for h'3 to a text file
  double step2 = T / N2;
  
  write_to_file("h3_prime.txt", step2, h3_prime, N2);

  free(H3);
  free(h3_prime);

  return 0;
}

double magnitude(CN c_num)
{
  /*
  Calculates the magnitude of a complex number.
  
  Parameter: c_num [CN] - Complex number.
  Returns: [double] - Magnitude of the input complex number.
  */
  return sqrt(c_num.re * c_num.re + c_num.im * c_num.im);
}

void DFT(CN *input, CN **output, int num)
{
  /*
  Calculates the Discrete Fourier Transform (DFT) of an array of complex
  numbers, transforming from the time domain to the frequency domain.
  
  Parameters:
    input [CN*] - The input signal in the time domain.
    output [CN**] - Where the calculated DFT result will be stored.
    num [int] - The length of the input array, representing the number of
      samples in the time domain.
  */
  if (input == NULL || num <= 0)
  {
    printf("Error: Invalid parameter values in DFT function.\n");
    
    return;
  }

  *output = (CN*)malloc(num * sizeof(CN));
  
  if (*output == NULL)
  {
    printf("Error: Memory allocation failed in DFT function.\n");
    
    return;
  }

  for (int i = 0; i < num; i++)
  {
    // The result element is set to zero to begin the summation
    (*output)[i].re = 0;
    (*output)[i].im = 0;
    
    /*
    For each index, the real and imaginary part of the frequency component
    is calculated by iteration
    */
    for (int n = 0; n < num; n++)
    {
      double angle = - 2 * M_PI * i * n / num; //in rads
      
      // Separate exponential into real and imaginary parts
      (*output)[i].re += input[n].re * cos(angle) - input[n].im * sin(angle);
      (*output)[i].im += input[n].re * sin(angle) + input[n].im * cos(angle);
    }
  }
}

int in_array(int arr[], int size_arr, int target, int cond_type)
{
  /*
  Checks if a value is in an array.
  
  Parameters:
    arr [int[]] - The array to search.
    size_arr [int] - The size of the array.
    target [int] - The value to find.
    cond_type [int] - The condition type. If cond_type is 1, the function
      returns 1 if the target is in the array. If cond_type is 0, the function
      returns 1 if the target is not in the array.
  
  Returns:
    [int] - The corresponding result of the search (1 or 0).
  */
  if (arr == NULL || size_arr <= 0 || !(cond_type == 0 || cond_type == 1))
  {
    printf("Error: Invalid parameter values in in_array function.\n");
    
    return -1;
  }
  
  // If condition type is "1" (True if in array)
  if (cond_type)
  {
    for (int i = 0; i < size_arr; i++)
    {
      if (arr[i] == target)
      {
        return 1;
      }
      
      return 0;
    }
    
  // If condition type is "0" (True if not in array)
  } else
  {
    for (int i = 0; i < size_arr; i++)
    {
      if (arr[i] != target)
      {
        return 1;
      }
      
      return 0;
    }
  }
}

void IDFT(CN *input, CN **output, int num, int els[], int size_arr,
  int type_con)
{
  /*
  Calculates the Inverse Discrete Fourier Transform: transforms from frequency
  domain to time domain.
  
  Parameters:
    input [CN*] - The input complex numbers in the frequency domain.
    output [CN**] - Where the output complex numbers in the time domain will be
      stored.
    num [int] - Length of the input and output arrays.
    els [int[]] - The indexes of the elements to be included/excluded in the
      summation.
    size_arr [int] - Size of the 'els' array.
    type_con [char] - Condition type for the 'els' array: 1 for inclusion, 0 for
      exclusion.
  */
  if (input == NULL || num <= 0)
  {
    printf("Error: Invalid parameter values in IDFT function.\n");
    
    return;
  }

  *output = (CN*)malloc(num * sizeof(CN));
  
  if (*output == NULL)
  {
    printf("Error: Memory allocation failed in IDFT function.\n");
    
    return;
  }

  for (int n = 0; n < num; n++)
  {
    // The result element is set to zero to begin the summation
    (*output)[n].re = 0;
    (*output)[n].im = 0;
    
    int counter = 0;
    
    /*
    For each index, the real and imaginary part of the time component is
    calculated by iteration
    */
    for (int i = 0; i < num; i++)
    {
      /*
      It is decided if this index element should be included or not,
      depending on the result of the in_array() function and the condition
      type (skip or include)
      */
      if (in_array(els, size_arr, i, type_con))
      {
        double angle = 2 * M_PI * i * n / num;
        
        (*output)[n].re += input[i].re * cos(angle) - input[i].im * sin(angle);
        (*output)[n].im += input[i].re * sin(angle) + input[i].im * cos(angle);
        
        counter++;
      }
    }
    
    // The final result is divided by the number of values
    (*output)[n].re /= num;
    (*output)[n].im /= num;
  }
}

void h1(double step, CN result[])
{
  /*
  Generates an array of complex numbers defined by the function h1.
  
  Parameters:
    step [double] - Step size of the time between iterations.
    result [CN[]] - Where the results will be stored (array of complex numbers).
  */
  for (int i = 0; i < N; i++)
  {
    double t = i * step;
    
    result[i].re = cos(t) + cos(5 * t);
    result[i].im = sin(t) + sin(5 * t);
  }
}

h2_result h2(double step)
{
  /*
  Generates an array of complex numbers defined by the function h2.
  
  Parameters:
    step [double] - Step size of time between iterations.

  Returns:
    [h2_result] - The results of the function.
  */
  h2_result result;
  
  for (int i = 0; i < N; i++)
  {
    double t = i * step;
    
    result.result[i].re = exp(pow((t - M_PI), 2) / 2);
    result.result[i].im = 0;
  }
  
  return result;
}

void write_to_file(const char* filename, double step, CN *h_values, int num)
{
  /*
  Writes complex values to a file in CSV format (Time, Real part, Imaginary
  part).
  
  Parameters:
    filename [const char*] - Filename to write the data to.
    step [double] - Step size of time.
    h_values [CN*] - Complex values to write to the file.
    num [int] - Length of the array h_values.
  */
  if (filename == NULL || h_values == NULL || num <= 0)
  {
    printf("Error: Invalid parameter values in write_to_file function.\n");
    
    return;
  }

  FILE *fp;
  fp = fopen(filename, "w");
  
  if (fp == NULL)
  {
    printf("Error: Unable to open file %s for writing.\n", filename);
    
    return;
  }

  fprintf(fp, "Time,real,imaginary\n");

  double t = 0;
  
  for (int i = 0; i < num; i++)
  {
    // Write values to file
    fprintf(fp, "%lf,%lf,%lf\n", t, h_values[i].re, h_values[i].im);
    
    t += step;
  }

  fclose(fp);
}

void find_max_magnitudes_indices(double magnitudes[], int indices[],
  int num_elements, int num_max)
{
  /*
  Finds indexes of the elements with maximum magnitude in an array.
  
  Parameters:
    magnitudes [double[]] - The magnitudes of elements.
    indices [int[]] - Where the indices of elements with the largest magnitudes
      will be stored.
    num_elements [int] - The number of elements in the magnitudes array.
    num_max [int] - The number of maximum magnitudes to find.
  */
  /*
  Initialize array of indices with -1, any magnitude will be larger than this
  because all magnitudes are positive
  */
  for (int i = 0; i < num_max; i++)
  {
    indices[i] = -1;
  }
  
  // Find indices of elements with largest magnitudes
  for (int i = 0; i < num_elements; i++)
  {
    for (int j = 0; j < num_max; j++)
    {
      if (indices[j] == -1 || magnitudes[i] > magnitudes[indices[j]])
      {
        // Shift indices to make space for new maximum
        for (int k = num_max - 1; k > j; k--)
        {
          indices[k] = indices[k - 1];
        }
        
        indices[j] = i;
        
        break;
      }
    }
  }
}
