#pragma once
#include <vector>
#include <cassert>
#include <iomanip>
#include <ostream>
#include <istream>

class Size {
public:
	Size(size_t w, size_t h)
		: width(w), height(h)
	{
	}
	
	size_t width, height;
};

std::ostream& operator<<(std::ostream& os, const Size& size)
{
	os << '(' << size.width << ',' << size.height << ')';
	return os;
}

bool operator==(const Size& a, const Size& b)
{
	return a.width == b.width && a.height == b.height;
}

bool operator!=(const Size& a, const Size& b)
{
	return a.width != b.width && a.height != b.height;
}

template<typename DType>
class Matrix {
public:
	Matrix();
	Matrix(Size size);
	Matrix(Size size, DType val);
	Matrix(size_t size);
	Matrix(size_t rows, size_t cols);
	Matrix(size_t rows, size_t cols, DType val);

	static Matrix zeros(Size size);
	static Matrix ones(Size size);
	static Matrix diag(size_t size, DType val=1);

	Matrix& inverse();

	Matrix& transpose()
	{
		for (int i = 0; i < mat.shape().height; i++)
			for (int j = i + 1; j < mat.shape().width; j++)
				std::swap(mat[i][j], mat[j][i]);
		return *this;
	}

	Size shape() const
	{
		return size;
	}

	void reshape(Size size)
	{
		this->size = size;
	}

	void reshape(size_t rows, size_t cols)
	{
		this->size = Size(cols, rows);
	}

	size_t at(size_t i, size_t j) const
	{
		return mat[i][j];
	}

	size_t& at(size_t i, size_t j)
	{
		return mat[i][j];
	}

	Matrix add(DType a) const
	{
		assert(a.shape() == shape());
		Matrix c(shape());
		for (int i = 0; i < shape().height; i++)
			for (int j = 0; j < shape().width; j++)
				c.at(i, j) = at(i, j) + a;
		return c;
	}

	Matrix& add(DType a)
	{
		assert(a.shape() == shape());
		Matrix c(shape());
		for (int i = 0; i < shape().height; i++)
			for (int j = 0; j < shape().width; j++)
				c.at(i, j) = at(i, j) + a;
		return c;
	}

	Matrix add(const Matrix& a) const
	{
		assert(a.shape() == shape());
		Matrix c(shape());
		for (int i = 0; i < shape().height; i++)
			for (int j = 0; j < shape().width; j++)
				c.at(i, j) = at(i, j) + a.at(i, j);
		return c;
	}

	Matrix& add(const Matrix& a)
	{
		assert(a.shape() == shape());
		Matrix c(shape());
		for (int i = 0; i < shape().height; i++)
			for (int j = 0; j < shape().width; j++)
				c.at(i, j) = at(i, j) + a.at(i, j);
		return c;
	}

	Matrix sub(DType a) const
	{
		assert(a.shape() == shape());
		Matrix c(shape());
		for (int i = 0; i < shape().height; i++)
			for (int j = 0; j < shape().width; j++)
				c.at(i, j) = at(i, j) - a;
		return c;
	}

	Matrix& sub(DType a)
	{
		assert(a.shape() == shape());
		Matrix c(shape());
		for (int i = 0; i < shape().height; i++)
			for (int j = 0; j < shape().width; j++)
				c.at(i, j) = at(i, j) - a;
		return c;
	}

	Matrix sub(const Matrix& a) const
	{
		assert(a.shape() == shape());
		Matrix c(shape());
		for (int i = 0; i < shape().height; i++)
			for (int j = 0; j < shape().width; j++)
				c.at(i, j) = at(i, j) - a.at(i, j);
		return c;
	}

	Matrix& sub(const Matrix& a)
	{
		assert(a.shape() == shape());
		Matrix c(shape());
		for (int i = 0; i < shape().height; i++)
			for (int j = 0; j < shape().width; j++)
				c.at(i, j) = at(i, j) - a.at(i, j);
		return c;
	}

	Matrix mul(DType a) const
	{
		assert(a.shape() == shape());
		Matrix c(shape());
		for (int i = 0; i < shape().height; i++)
			for (int j = 0; j < shape().width; j++)
				c.at(i, j) = at(i, j) * a;
		return c;
	}

	Matrix& mul(DType a)
	{
		assert(a.shape() == shape());
		Matrix c(shape());
		for (int i = 0; i < shape().height; i++)
			for (int j = 0; j < shape().width; j++)
				c.at(i, j) = at(i, j) * a;
		return c;
	}

	Matrix mul(const Matrix& a) const
	{
		assert(a.shape() == shape());
		Matrix c(shape());
		for (int i = 0; i < shape().height; i++)
			for (int j = 0; j < shape().width; j++)
				c.at(i, j) = at(i, j) * a.at(i, j);
		return c;
	}

	Matrix& mul(const Matrix& a)
	{
		assert(a.shape() == shape());
		Matrix c(shape());
		for (int i = 0; i < shape().height; i++)
			for (int j = 0; j < shape().width; j++)
				c.at(i, j) = at(i, j) * a.at(i, j);
		return c;
	}

	Matrix matmul(const Matrix& a) const
	{
		assert(shape().width == a.shape().height);
		Matrix c(Size(shape().height, a.shape().width));
		for (int i = 0; i < c.shape().height; i++)
			for (int j = 0; j < c.shape().width; j++)
				for (int k = 0; k < shape().width; k++)
					c.at(i, j) += at(i, k) * a.at(k, j);
		return c;
	}

	Matrix& matmul(const Matrix& a)
	{
		assert(shape().width == a.shape().height);
		Matrix c(Size(shape().height, a.shape().width));
		for (int i = 0; i < c.shape().height; i++)
			for (int j = 0; j < c.shape().width; j++)
				for (int k = 0; k < shape().width; k++)
					c.at(i, j) += at(i, k) * a.at(k, j);
		return c;
	}

	Matrix div(DType a) const
	{
		assert(a.shape() == shape());
		Matrix c(shape());
		for (int i = 0; i < shape().height; i++)
			for (int j = 0; j < shape().width; j++)
				c.at(i, j) = at(i, j) / a;
		return c;
	}

	Matrix& div(DType a)
	{
		assert(a.shape() == shape());
		Matrix c(shape());
		for (int i = 0; i < shape().height; i++)
			for (int j = 0; j < shape().width; j++)
				c.at(i, j) = at(i, j) / a;
		return c;
	}

	Matrix div(const Matrix& a) const
	{
		assert(a.shape() == shape());
		Matrix c(shape());
		for (int i = 0; i < shape().height; i++)
			for (int j = 0; j < shape().width; j++)
				c.at(i, j) = at(i, j) / a.at(i, j);
		return c;
	}

	Matrix& div(const Matrix& a)
	{
		assert(a.shape() == shape());
		Matrix c(shape());
		for (int i = 0; i < shape().height; i++)
			for (int j = 0; j < shape().width; j++)
				c.at(i, j) = at(i, j) / a.at(i, j);
		return c;
	}

	Matrix operator+(DType a) const
	{
		return add(a);
	}

	Matrix& operator+(DType a)
	{
		return add(a);
	}

	Matrix operator+(const Matrix& a) const
	{
		return add(a);
	}

	Matrix& operator+(const Matrix& a)
	{
		return add(a);
	}

	Matrix operator-(DType a) const
	{
		return sub(a);
	}

	Matrix& operator-(DType a)
	{
		return sub(a);
	}

	Matrix operator-(const Matrix& a) const
	{
		return sub(a);
	}

	Matrix& operator-(const Matrix& a)
	{
		return sub(a);
	}

	Matrix operator*(DType a) const
	{
		return mul(a);
	}

	Matrix& operator*(DType a)
	{
		return mul(a);
	}

	Matrix operator*(const Matrix& a) const
	{
		return matmul(a);
	}

	Matrix& operator*(const Matrix& a)
	{
		return matmul(a);
	}

	Matrix operator/(DType a) const
	{
		return div(a);
	}

	Matrix& operator/(DType a)
	{
		return div(a);
	}

	Matrix operator/(const Matrix& a) const
	{
		return div(a);
	}

	Matrix& operator/(const Matrix& a)
	{
		return div(a);
	}

	Matrix& operator+=(DType a)
	{
		*this = *this + a;
		return *this;
	}

	Matrix& operator+=(const Matrix& a)
	{
		*this = *this + a;
		return *this;
	}

	Matrix& operator-=(DType a)
	{
		*this = *this - a;
		return *this;
	}

	Matrix& operator-=(const Matrix& a)
	{
		*this = *this - a;;
		return *this;
	}

	Matrix& operator*=(DType a)
	{
		*this = *this * a;
		return *this;
	}

	Matrix& operator*=(const Matrix& a)
	{
		*this = matmul(a);
		return *this;
	}

	Matrix& operator/=(DType a)
	{
		*this = *this / a;
		return *this;
	}

	Matrix& operator/=(const Matrix& a)
	{
		*this = *this / a;
		return *this;
	}


private:

	std::vector<std::vector<DType>> mat;
	Size size;
};
