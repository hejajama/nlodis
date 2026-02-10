// Unit tests for Vec class

#include "unit_test_framework.hpp"

#include <cmath>
#include <sstream>

#include "../dipole/vector.hpp"

namespace {
constexpr double kEps = 1e-12;
}

TEST(VEC_DEFAULT_CONSTRUCTOR)
{
	Vec v;
	ASSERT_ALMOST_EQUAL(v.GetX(), 0.0, kEps);
	ASSERT_ALMOST_EQUAL(v.GetY(), 0.0, kEps);
	ASSERT_ALMOST_EQUAL(v.GetZ(), 0.0, kEps);
}

TEST(VEC_CONSTRUCTORS_AND_COPY)
{
	Vec v2(1.5, -2.0);
	ASSERT_ALMOST_EQUAL(v2.GetX(), 1.5, kEps);
	ASSERT_ALMOST_EQUAL(v2.GetY(), -2.0, kEps);
	ASSERT_ALMOST_EQUAL(v2.GetZ(), 0.0, kEps);

	Vec v3(1.0, 2.0, 3.0);
	ASSERT_ALMOST_EQUAL(v3.GetX(), 1.0, kEps);
	ASSERT_ALMOST_EQUAL(v3.GetY(), 2.0, kEps);
	ASSERT_ALMOST_EQUAL(v3.GetZ(), 3.0, kEps);

	Vec vcopy(v3);
	ASSERT_ALMOST_EQUAL(vcopy.GetX(), 1.0, kEps);
	ASSERT_ALMOST_EQUAL(vcopy.GetY(), 2.0, kEps);
	ASSERT_ALMOST_EQUAL(vcopy.GetZ(), 3.0, kEps);
}

TEST(VEC_SETTERS_GETTERS)
{
	Vec v;
	v.SetX(3.0);
	v.SetY(-4.0);
	v.SetZ(5.0);
	ASSERT_ALMOST_EQUAL(v.GetX(), 3.0, kEps);
	ASSERT_ALMOST_EQUAL(v.GetY(), -4.0, kEps);
	ASSERT_ALMOST_EQUAL(v.GetZ(), 5.0, kEps);
}

TEST(VEC_ARITHMETIC_OPERATORS)
{
	Vec a(1.0, 2.0, 3.0);
	Vec b(-4.0, 5.0, -6.0);

	Vec sum = a + b;
	ASSERT_ALMOST_EQUAL(sum.GetX(), -3.0, kEps);
	ASSERT_ALMOST_EQUAL(sum.GetY(), 7.0, kEps);
	ASSERT_ALMOST_EQUAL(sum.GetZ(), -3.0, kEps);

	Vec diff = a - b;
	ASSERT_ALMOST_EQUAL(diff.GetX(), 5.0, kEps);
	ASSERT_ALMOST_EQUAL(diff.GetY(), -3.0, kEps);
	ASSERT_ALMOST_EQUAL(diff.GetZ(), 9.0, kEps);

	Vec c(1.0, 1.0, 1.0);
	c += b;
	ASSERT_ALMOST_EQUAL(c.GetX(), -3.0, kEps);
	ASSERT_ALMOST_EQUAL(c.GetY(), 6.0, kEps);
	ASSERT_ALMOST_EQUAL(c.GetZ(), -5.0, kEps);

	Vec d(1.0, 1.0, 1.0);
	d -= b;
	ASSERT_ALMOST_EQUAL(d.GetX(), 5.0, kEps);
	ASSERT_ALMOST_EQUAL(d.GetY(), -4.0, kEps);
	ASSERT_ALMOST_EQUAL(d.GetZ(), 7.0, kEps);
}

TEST(VEC_ASSIGNMENT)
{
	Vec a(1.0, 2.0, 3.0);
	Vec b(0.0, 0.0, 0.0);
	b = a;
	ASSERT_ALMOST_EQUAL(b.GetX(), 1.0, kEps);
	ASSERT_ALMOST_EQUAL(b.GetY(), 2.0, kEps);
	ASSERT_ALMOST_EQUAL(b.GetZ(), 3.0, kEps);
}

TEST(VEC_SCALAR_OPERATORS)
{
	Vec a(1.0, -2.0, 3.0);

	Vec scaled = a * 2.0;
	ASSERT_ALMOST_EQUAL(scaled.GetX(), 2.0, kEps);
	ASSERT_ALMOST_EQUAL(scaled.GetY(), -4.0, kEps);
	ASSERT_ALMOST_EQUAL(scaled.GetZ(), 6.0, kEps);

	a *= -3.0;
	ASSERT_ALMOST_EQUAL(a.GetX(), -3.0, kEps);
	ASSERT_ALMOST_EQUAL(a.GetY(), 6.0, kEps);
	ASSERT_ALMOST_EQUAL(a.GetZ(), -9.0, kEps);
}

TEST(VEC_DOT_PRODUCT)
{
	Vec a(1.0, 2.0, 3.0);
	Vec b(4.0, -5.0, 6.0);
	double dot = a * b;
	ASSERT_ALMOST_EQUAL(dot, 12.0, kEps); // 1*4 + 2*(-5) + 3*6 = 12
}

TEST(VEC_LENGTHS)
{
	Vec a(3.0, 4.0, 12.0);
	ASSERT_ALMOST_EQUAL(a.LenSqr(), 169.0, kEps);
	ASSERT_ALMOST_EQUAL(a.Len(), 13.0, kEps);
}

TEST(VEC_ROTATE_2D)
{
	Vec a(1.0, 0.0, 0.0);
	a.Rotate2D(M_PI / 2.0);
	ASSERT_ALMOST_EQUAL(a.GetX(), 0.0, 1e-10);
	ASSERT_ALMOST_EQUAL(a.GetY(), -1.0, 1e-10);
	ASSERT_ALMOST_EQUAL(a.GetZ(), 0.0, kEps);
}

TEST(VEC_STREAM_OUTPUT)
{
	Vec a(1.0, 2.0, 3.0);
	std::ostringstream os;
	os << a;
	const std::string out = os.str();
	ASSERT_TRUE(out.find("Vector (1") != std::string::npos);
	ASSERT_TRUE(out.find("2") != std::string::npos);
	ASSERT_TRUE(out.find("3") != std::string::npos);
}
