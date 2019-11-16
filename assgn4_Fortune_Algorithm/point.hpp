// NAME:Akash Tiwari
// ROLL:17CS10003

#include <iostream>
#include <time.h>
#include <iomanip>
#include <random>
#include <cmath>
#include <memory>
#include <vector>
#include <cassert>
#include <limits>

#ifndef point_hpp
#define point_hpp

#define POINT_EPSILON 1.0e-6

using namespace std;

class Pointxy 
{
	struct Pointxy_Compare 
    {
		bool operator()(const Pointxy &p1, const Pointxy &p2)
        {
			return (p1.x<p2.x||(p1.x==p2.x&&p1.y< p2.y));
		}
	};

public:
    double x, y;
    int id;
    const static double Inf;
	static Pointxy_Compare xy_compare;
    Pointxy &operator-=(const Pointxy &p);
    friend double scalar_prod(const Pointxy &p1, const Pointxy &p2);
    friend double vector_prod(const Pointxy &p1, const Pointxy &p2);
    friend Pointxy operator*(const Pointxy &p, double value);
    Pointxy(double x = 0.0, double y = 0.0);
    Pointxy(const Pointxy &point);
    Pointxy &operator*=(double value);
    friend Pointxy operator+(const Pointxy &p1, const Pointxy &p2);
    friend Pointxy operator/(const Pointxy &p1, const Pointxy &p2);
    friend Pointxy operator*(double value, const Pointxy &p);
    Pointxy &operator+=(const Pointxy &p);
    bool checkVertical();
    bool checkHorizontal();
    bool isValid();
    friend Pointxy operator/(const Pointxy &p, double value);
    friend vector<Pointxy> &operator<<(vector<Pointxy> &v, const Pointxy &p);
    friend Pointxy operator-(const Pointxy &p1, const Pointxy &p2);
    Pointxy &operator/=(double value);
    Pointxy Rotate90CW();
    Pointxy Rotate90CCW();
    friend ostream &operator<<(ostream &stream, const Pointxy &p);
    friend Pointxy operator-(const Pointxy &p);
    Pointxy normalized();
    void normalize();
    double norm();
    double norm2();
	static bool checkLeftTurn(const Pointxy &p1, const Pointxy &p2, const Pointxy &p3);
	static bool checkRightTurn(const Pointxy &p1, const Pointxy &p2, const Pointxy &p3);
    double operator[](int i);
    void setX(double x);
    void setY(double y);
    void setId(int id);
};

double scalar_prod(const Pointxy &p1, const Pointxy &p2);
double vector_prod(const Pointxy &p1, const Pointxy &p2);
bool equal(const Pointxy &p1, const Pointxy &p2, double EPSILON = POINT_EPSILON);
bool equal(double v1, double v2, double EPSILON = POINT_EPSILON);
int intersectionPoints(const Pointxy &f1, const Pointxy &f2, double directrix);
vector<Pointxy> findIntersectionPoints(const Pointxy &f1, const Pointxy &f2, double directrix);
#endif