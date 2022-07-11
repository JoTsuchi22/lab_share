#include <vector>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <fstream>
#include <sstream>

using namespace std;

template<class T>
class vector_matrix {
    T *array;
    int length;
    int row;
    int col;

    vector_matrix(int a)
    {
        row = a;
        col = 1;
        length = row * col;
        array = (T *)calloc(length, sizeof(T));
    }
    vector_matrix(int a, int b)
    {
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
    ~vector_matrix()
    {
        free(array);
    }
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
            cout << array[j * row + i] << " ";
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


int main()
{
    vector_matrix <int>x(10);
    vector_matrix <double>A(5, 7);

    x.initialization();
    x.print();

    A.initialization();
    A.print();

    x.print(2);
    A.print(2, 3);

    for (int i = 0; i < x.length; i++)
        x.array[i] = 1;

    for (int i = 0; i < A.length; i++)
        A.array[i] = 3;

    double temp;
    temp = x.val(2) + A.val(2, 3);

    cout << temp << endl;
}