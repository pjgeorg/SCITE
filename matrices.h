/*
 * matrices.h
 *
 *  Created on: Mar 27, 2015
 *      Author: jahnka
 */

#ifndef MATRICES_H
#define MATRICES_H

#include <vector>
#include <array>

template<class T, std::size_t tSize>
using StaticArray = std::array<T, tSize>;

// Temporary, should be removed after refactor, causes mem-leaks
template<class T, std::size_t tSize>
auto toStaticArray(T const*const data)
{
    StaticArray<T, tSize> arr;
    for(std::size_t i = 0; i<tSize; ++i)
    {
        arr[i] = data[i];
    }

//    delete [] data;

    return arr;
}

template<class T>
using DynamicArray = std::vector<T>;

// Temporary, should be removed after refactor, causes mem-leaks
template<class T>
auto toDynamicArray(T const*const data, std::size_t size)
{
    DynamicArray<T> vec(size);
    for(std::size_t i = 0; i<size; ++i)
    {
        vec[i] = data[i];
    }

//    delete [] data;

    return vec;
}

template<class T, std::size_t N>
class Vector
{
    public:
        auto &operator[](std::size_t const i)
        {
            return mData[i];
        }

        auto &operator[](std::size_t const i) const
        {
            return mData[i];
        }

    private:
        std::array<T, N> mData;
};

template<class T, std::size_t M, std::size_t N>
class Matrix
{
    public:
        auto &operator[](std::size_t const i)
        {
            return mData[i];
        }

        auto &operator[](std::size_t const i) const
        {
            return mData[i];
        }

    private:
        Vector<Vector<T, N>, M> mData;
};

template<class T, std::size_t M, std::size_t N>
static inline auto toPointer(Matrix<T, M, N> &a)
{
    T** data = new T*[M];
    for(std::size_t i = 0; i<M; ++i)
        data[i] = reinterpret_cast<T*>(&a) + i*N;

    return data;
}

int** sumMatrices(int** first, int** second, int n, int m);
double getMaxEntry(double* array, int n);
void addToMatrix(int** first, int** second, int n, int m);
double** allocate_doubleMatrix(int n, int m);
int** allocate_intMatrix(int n, int m);
bool** allocate_boolMatrix(int n, int m);
double** init_doubleMatrix(int n, int m, double value);
int* init_intArray(int n, int value);
bool* init_boolArray(int n, bool value);
int** init_intMatrix(int n, int m, int value);
bool** init_boolMatrix(int n, int m, bool value);
double* init_doubleArray(int n, double value);
void reset_intMatrix(int** matrix, int n, int m, int value);
void free_boolMatrix(bool** matrix);
void free_intMatrix(int** matrix);
void free_doubleMatrix(double** matrix);
bool** deepCopy_boolMatrix(bool** matrix, int n, int m);
int** deepCopy_intMatrix(int** matrix, int n, int m);
int* deepCopy_intArray(int* array, int n);
double* deepCopy_doubleArray(double* array, int n);
double** deepCopy_doubleMatrix(double** matrix, int n, int m);
void print_boolMatrix(bool** array, int n, int m);
void print_doubleMatrix(double** matrix, int n, int m);
void print_intMatrix(int** matrix, int n, int m, char del);
void print_intArray(int* array, int n);
int* ancMatrixToParVector(bool** anc, int n);
bool identical_boolMatrices(bool** first, bool** second, int n, int m);
void delete_3D_intMatrix(int*** matrix, int n);
#endif
