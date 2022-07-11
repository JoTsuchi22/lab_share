#include <vector>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <fstream>
#include <sstream>

#include <unistd.h>

using namespace std;

// class vector_matrix
template<class T>
class vector_matrix {
public:
    T *array;
    int length;
    int row;
    int col;
    int vec_bool;

    vector_matrix(int a)
    {
        vec_bool = 1;
        row = a;
        col = 1;
        length = row * col;
        array = (T *)calloc(length, sizeof(T));
    }
    vector_matrix(int a, int b)
    {
        vec_bool = 0;
        row = a;
        col = b;
        length = row * col;
        array = (T *)calloc(length, sizeof(T));
    }
    void initialization();
    void substitute(T num);
    T val(int i_row);
    T val(int i_row, int i_col);
    void print();
    void print(int i_row);
    void print(int i_row, int i_col);
    ~vector_matrix() { free(array); }
};


template<class T>
void vector_matrix<T>::initialization()
{
    int i;
    for (i = 0; i < length; i++)
        array[i] = i;
}


template<class T>
void vector_matrix<T>::substitute(T num)
{
    int i;
    for (i = 0; i < length; i++)
        array[i] = num;
}


template<class T>
T vector_matrix<T>::val(int i_row)
{
    return array[i_row];
}


template<class T>
T vector_matrix<T>::val(int i_row, int i_col)
{
    return array[i_row * col + i_col];
}


template<class T>
void vector_matrix<T>::print()
{
    int i, j;
    cout << scientific;
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            cout << array[i * col + j] << " ";
            if (col != 1 && j == col - 1)
                cout << endl;
        }
        if (i == row - 1)
                cout << endl;
    }
    if (col == 1)
        cout << endl;
}

template<class T>
void vector_matrix<T>::print(int i_row)
{
    cout << scientific;
    cout << array[i_row] << endl << endl;
}


template<class T>
void vector_matrix<T>::print(int i_row, int i_col)
{
    cout << scientific;
    cout << array[i_row * col + i_col] << endl << endl;
}


// dot
template<class T>
T dot(vector_matrix<T> *A, vector_matrix<T> *B)
{
    int i;
    T val = 0.0;
    if (A->vec_bool && B->vec_bool)
        for (i = 0; i < A->length; i++)
            val += A->array[i] * B->array[i];
    else
    {
        cout << "Error in dot product" << endl;
    }

    return val;
}


// matrix multiplication
template<class T>
void matrix_multiplication(vector_matrix<T> *C, vector_matrix<T> *A, vector_matrix<T> *B)
{
    int i, j, k;
    for (i = 0; i < C->length; i++)
        C->array[i] = 0;
    if (!(A->vec_bool && B->vec_bool) && (A->col == B->row))
        for (i = 0; i < A->row; i++)
            for (k = 0; k < A->col; k++)
                for (j = 0; j < B->col; j++)
                    C->array[i * C->col + j] += A->val(i, k) * B->val(k, j);
    else
    {
        cout << "Error in matrix multiplication" << endl;
    }
}


int main()
{
    vector_matrix <int>x(10);
    vector_matrix <int>y(10);

    x.initialization();
    x.print();
    
    y.initialization();
    y.print();
    
    double value;
    value = dot(&x, &y);

    cout << value << endl << endl;

    vector_matrix <double>A(2, 3);
    vector_matrix <double>B(3, 4);

    A.initialization();
    A.print();

    B.initialization();
    B.print();

    vector_matrix <double>C(2, 4);

    matrix_multiplication(&C, &A, &B);

    C.print();

    long int count = 0;
    for (;;)
    {
        printf("\rProgress %d / 100", (int)count++);
        usleep(30000);
    }
}