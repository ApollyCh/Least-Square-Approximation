# Least-Square-Approximation
Write a computer program in C++ programming language to compose a least square approximation for a given data set.

The input contains:

the length m of data set
m lines with experimental data t_i *space* b_i
the degree of the polynomial n

The output contains:

the matrix A itself after the line "A:"
the matrix (AT A) after the line "A_T*A:"
the matrix (AT A)-1 after the line "(A_T*A)_-1:"
the matrix AT b after the line "A_T*b:"
the answer itself after the line "x:"

Example: 

Input:

3

-2 1

0 2

2 4

1

Output:

A:

1.0000 -2.0000

1.0000 0.0000

1.0000 2.0000

A_T*A:

3.0000 0.0000

0.0000 8.0000

(A_T*A)^-1:

0.3333 0.0000

0.0000 0.1250

A_T*b:

7.0000

6.0000

x~:

2.3333

0.7500
