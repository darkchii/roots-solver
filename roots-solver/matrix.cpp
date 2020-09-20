#include "matrix.h"


template<typename DType>
Matrix<DType>::Matrix()
	: mat(), size()
{
}

template<typename DType>
Matrix<DType>::Matrix(Size size)
	: mat(size.height, std::vector<DType>(size.width, 0))
	, size(size)
{
}

template<typename DType>
Matrix<DType>::Matrix(Size size, DType val)
	: mat(size.height, std::vector<DType>(size.width, val))
	, size(size)
{
}

template<typename DType>
Matrix<DType>::Matrix(size_t size)
	: mat(size, std::vector<DType>(size, 0))
	, size(Size(size, size))
{
}

template<typename DType>
Matrix<DType>::Matrix(size_t rows, size_t cols)
	: mat(rows, std::vector<DType>(cols, 0))
	, size(Size(cols, rows))
{
}

template<typename DType>
Matrix<DType>::Matrix(size_t rows, size_t cols, DType val)
	: mat(rows, std::vector<DType>(cols, val))
	, size(Size(cols, rows))
{
}

template<typename DType>
Matrix<DType> Matrix<DType>::zeros(Size size)
{
	return Matrix(size);
}

template<typename DType>
Matrix<DType> Matrix<DType>::ones(Size size)
{
	return Matrix<DType>(size, 1);
}

template<typename DType>
Matrix<DType> Matrix<DType>::diag(size_t size, DType val)
{
	Matrix<DType> mat = Matrix<DType>(size);
	for (size_t i = 0; i < size; i++)
		mat.at(i, i) = val;
	return mat;
}

template<>
Matrix<double>::Matrix()
	: mat(), size(size)
{
}

template<>
Matrix<double>::Matrix(Size size)
	: mat(size.height, std::vector<double>(size.width, 0))
	, size(size)
{
}

template<>
Matrix<double> Matrix<double>::zeros(Size size)
{
	return Matrix(size);
}

template<>
Matrix<double> Matrix<double>::ones(Size size)
{
	return Matrix(size, 1);
}

template<typename DType>
Matrix<DType> operator+(DType a, const Matrix<DType>& b)
{
	return b + a;
}

template<typename DType>
Matrix<DType> operator-(DType a, const Matrix<DType>& b)
{
	return b - a;
}

template<typename DType>
Matrix<DType> operator*(DType a, const Matrix<DType>& b)
{
	return b * a;
}

template<typename DType>
Matrix<DType> operator/(DType a, const Matrix<DType>& b)
{
	return b / a;
}

template<typename DType>
Matrix<DType>& Matrix<DType>::inverse()
{
	Matrix _umat = diag(shape());
	return {};
}