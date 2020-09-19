#include <iostream>
#include <sstream>
#include <ostream>
#include <vector>
#include <cassert>
#include <iomanip>

using namespace std;

bool lup_decomposition(vector<vector<double>> a, vector<vector<double>>& l, vector<vector<double>>& u, vector<double>& p)
{
	for (int i = 0; i < a.size(); i++)
		p[i] = (double)i;
	for (int i = 0; i < a.size(); i++)
	{
		l[i][i] = 1;
		for (int j = i + 1; j < a.size(); j++)
			l[i][j] = 0;
	}
	for (int i = 1; i < a.size(); i++)
		for (int j = 0; j < i; j++)
			u[i][j] = 0;
	for (int i = 0; i < a.size(); i++)
	{
		double n = 0; int i_;
		for (int j = i; j < a.size(); j++)
			if (abs(a[j][i]) > n)
				n = abs(a[j][i]), i_ = j;
		if (n == 0) return false;
		swap(p[i], p[i_]);
		for (int j = 0; j < a.size(); j++)
			swap(a[i][j], a[i_][j]);
		u[i][i] = a[i][i];
		for (int j = i + 1; j < a.size(); j++)
		{
			l[j][i] = a[j][i] / u[i][i];
			u[i][j] = a[i][j];
		}
		for (int j = i + 1; j < a.size(); j++)
			for (int k = i + 1; k < a.size(); k++)
				a[j][k] -= l[j][i] * u[i][k];
	}
	return true;
}

void lu_solve(vector<vector<double>> a, vector<double> b, vector<vector<double>> l, vector<vector<double>> u, vector<double>& x)
{
	vector<double> y(a.size(), 0);
	for (int i = 0; i < a.size(); i++)
	{
		double sum = 0;
		for (int j = 0; j < i; j++)
			sum += l[i][j] * y[j];
		y[i] = b[i] - sum;
	}
	for (int i = a.size() - 1; i >= 0; i--)
	{
		double sum = 0;
		for (int j = i + 1; j < a.size(); j++)
			sum += u[i][j] * x[j];
		x[i] = (y[i] - sum) / u[i][i];
	}
}

void lup_solve(vector<vector<double>> a, vector<double> b, vector<vector<double>> l, vector<vector<double>> u, vector<double> p, vector<double>& x)
{
	vector<double> y(a.size(), 0);
	for (int i = 0; i < a.size(); i++)
	{
		double sum = 0;
		for (int j = 0; j < i; j++)
			sum += l[i][j] * y[j];
		y[i] = b[p[i]] - sum;
	}
	for (int i = a.size() - 1; i >= 0; i--)
	{
		double sum = 0;
		for (int j = i + 1; j < a.size(); j++)
			sum += u[i][j] * x[j];
		x[i] = (y[i] - sum) / u[i][i];
	}
}

vector<vector<double>> inverse(vector<vector<double>> a)
{
	vector<vector<double>>
		unit_matrix(a.size(), vector<double>(a.size(), 0)),
		x(a.size(), vector<double>(a.size(), 0)),
		l(a.size(), vector<double>(a.size(), 0)),
		u(a.size(), vector<double>(a.size(), 0));
	vector<double> p(a.size());
	for (int i = 0; i < a.size(); i++)
		unit_matrix[i][i] = 1;
	lup_decomposition(a, l, u, p);
	for (int i = 0; i < a.size(); i++)
		lup_solve(a, unit_matrix[i], l, u, p, x[i]);
	return x;
}

vector<double> func(vector<double> e, vector<double> a)
{
	vector<double> f(e.size(), 0);
	f[0] = a[0] / e[0] + e[2] + (a[5] - e[3]) * e[3] - a[4];
	f[1] = (a[5] - e[3]) * e[2] + e[1] + (a[0] * e[3]) / e[0] - a[3];
	f[2] = e[0] + (a[0] * e[2]) / e[0] + (a[5] - e[3]) * e[1] - a[2];
	f[3] = (a[0] * e[1]) / e[0] + (a[5] - e[3]) * e[0] - a[1];
	return f;
}

vector<vector<double>> jacob_matrix(vector<double> e, vector<double> a)
{
	vector<vector<double>> m(e.size(), vector<double>(e.size(), 0));
	double e2 = e[0] * e[0];

	m[0][0] = -a[0] / e2;
	m[0][2] = 1;
	m[0][3] = a[5] - 2 * e[3];

	m[1][0] = -(a[0] * e[3]) / e2;
	m[1][1] = 1;
	m[1][2] = a[5] - e[3];
	m[1][3] = a[0] / e[0] - e[2];

	m[2][0] = 1 - (a[0] * e[2]) / e2;
	m[2][1] = a[5] - e[3];
	m[2][2] = a[0] / e[0];
	m[2][3] = -e[1];

	m[3][0] = -(a[0] * e[1]) / e2 + a[5] - e[3];
	m[3][1] = a[0] / e[0];
	m[3][3] = -e[0];

	return m;
}

vector<double> sum(vector<vector<double>> a, int axis = 1)
{
	if (axis == 1)
	{
		vector<double> e(a.size());
		for (int i = 0; i < a.size(); i++)
			for (int j = 0; j < a[0].size(); j++)
				e[i] += a[i][j];
		return e;
	}
	else if (axis == -1)
	{
		vector<double> e(a[0].size());
		for (int i = 0; i < a[0].size(); i++)
			for (int j = 0; j < a.size(); j++)
				e[i] += a[i][j];
		return e;
	}
	else {
		vector<double> e(1);
		for (int i = 0; i < a.size(); i++)
			for (int j = 0; j < a[0].size(); j++)
				e[0] += a[i][j];
		return e;
	}
}

double sum(vector<double> a)
{
	double e = 0;
	for (auto i : a)
		e += i;
	return e;
}

vector<double> add(vector<double> a, vector<double> b)
{
	assert(a.size() == b.size());
	vector<double> c(a.size());
	for (int i = 0; i < a.size(); i++)
		c[i] = a[i] + b[i];
	return c;
}

vector<double> sub(vector<double> a, vector<double> b)
{
	assert(a.size() == b.size());
	vector<double> c(a.size());
	for (int i = 0; i < a.size(); i++)
		c[i] = a[i] - b[i];
	return c;
}

vector<vector<double>> sub(vector<vector<double>> a, vector<vector<double>> b)
{
	assert(a.size() == b.size());
	assert(a[0].size() == b[0].size());
	vector<vector<double>> c(a.size(), vector<double>(a[0].size(), 0));
	for (int i = 0; i < a.size(); i++)
		for (int j = 0; j < a[0].size(); j++)
			c[i][j] = a[i][j] - b[i][j];
	return c;
}

vector<vector<double>> add(vector<vector<double>> a, vector<vector<double>> b)
{
	assert(a.size() == b.size());
	assert(a[0].size() == b[0].size());
	vector<vector<double>> c(a.size(), vector<double>(a[0].size(), 0));
	for (int i = 0; i < a.size(); i++)
		for (int j = 0; j < a[0].size(); j++)
			c[i][j] = a[i][j] + b[i][j];
	return c;
}

vector<double> mul(vector<vector<double>> a, vector<double> b)
{
	assert(a[0].size() == b.size());
	vector<double> c(b.size(), 0);
	for (int i = 0; i < a.size(); i++)
		for (int j = 0; j < a[0].size(); j++)
			c[i] += a[i][j] * b[j];
	return c;
}

vector<vector<double>> mul(vector<vector<double>> a, vector<vector<double>> b)
{
	assert(a[0].size() == b.size());
	vector<vector<double>> r(a.size(), vector<double>(b[0].size(), 0));
	for (int i = 0; i < r.size(); i++)
		for (int j = 0; j < r[0].size(); j++)
			for (int k = 0; k < b.size(); k++)
				r[i][j] += a[i][k] * b[k][j];
	return r;
}

vector<vector<double>> transpose(vector<vector<double>> a)
{
	for (int i = 0; i < a.size(); i++)
		for (int j = i + 1; j < a[0].size(); j++)
			swap(a[i][j], a[j][i]);
	return a;
}

ostream& operator<<(ostream& os, vector<vector<double>> a)
{
	os << "[\n";
	for (auto rows : a)
	{
		for (auto col : rows)
			os << right << setw(14) << col << ' ';
		os << '\n';
	}
	os << ']';
	return os;
}

ostream& operator<<(ostream& os, vector<double> a)
{
	for (auto i : a) os << i << '\n';
	return os;
}

// Euclidean division of polynomials.
// Where u, v are coefficient and intercept of Quadratic equation x^2 + ux + v
vector<double> edp(double u, double v, vector<double> coeff)
{
	int n = coeff.size();
	vector<double> b(n);
	double rf = 0,
		   rs = 1;
	for (int i = 0; i < b.size(); i++)
	{
		b[i] = (coeff[n - i - 1] - v * rf) - (u * rs);
		rf = rs;
		rs = b[i];
	}
	return b;
}

int main()
{
	vector<double> a{ -4., 0., -2.};
	vector<double> u{0, 0},
		r{1, 1};
	double epsilon = 1e-5;
	while (abs(sum(r)) > epsilon)
	{
		cout << "多项式分解：\n";
		vector<double> b = edp(u[0], u[1], a);
		cout << b << '\n';

		cout << "求逆：\n";
		vector<vector<double>> z{ {-b[0], -1}, {0, -b[0]} };
		//cout << z << '\n';
		z = transpose(inverse(z));
		cout << z << '\n';
		
		cout << "余项：\n";
		r = { b[b.size() - 2], b[b.size() - 1] };
		cout << r << '\n';

		cout << "牛顿迭代：\n";
		u = sub(u, mul(z, r));
		cout << u << '\n';
	}
	
	//vector<double> e{ -2000., -100., 20., 1. };
	//int iteration = 13;
	//while (iteration--)
	//{
	//	vector<double> f = func(e, a);
	//	vector<vector<double>> z = jacob_matrix(e, a);
	//	vector<vector<double>> z_T = transpose(inverse(z));
	//	e = sub(e, mul(z_T, f));
	//	cout << "F:\n";
	//	cout << f << '\n';

	//	cout << "Z:\n";
	//	cout << z << '\n';

	//	cout << "Z^T:\n";
	//	cout << z_T << '\n';

	//	cout << "e:\n";
	//	cout << e << '\n';
	//}
	//double d1 = a[5] - e[3],
	//	d0 = a[0] / e[0];
	//cout << "(x^2 + " << d1 << "x + " << d0 << ")(x^4 + " << e[3] << "x^3 + " << e[2] << "x^2 + " << e[1] << "x + " << e[0] << ") = 0\n";
	return 0;
}
