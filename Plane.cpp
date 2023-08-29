#include <iomanip>
#include <iostream>
#include "Plane.hpp"



Plane::Plane(void) {

    a = 0;
    b = 0;
    c = 0;
    d = 0;
}

Plane::Plane(Vector pt1, Vector pt2, Vector pt3) {

    mpq_class& x1 = pt1.getX();
    mpq_class& y1 = pt1.getY();
    mpq_class& z1 = pt1.getZ();
    mpq_class& x2 = pt2.getX();
    mpq_class& y2 = pt2.getY();
    mpq_class& z2 = pt2.getZ();
    mpq_class& x3 = pt3.getX();
    mpq_class& y3 = pt3.getY();
    mpq_class& z3 = pt3.getZ();

    mpq_class a1 = x2 - x1;
    mpq_class b1 = y2 - y1;
    mpq_class c1 = z2 - z1;
    mpq_class a2 = x3 - x1;
    mpq_class b2 = y3 - y1;
    mpq_class c2 = z3 - z1;
    this->a = b1 * c2 - b2 * c1;
    this->b = a2 * c1 - a1 * c2;
    this->c = a1 * b2 - b1 * a2;
    this->d = (- a * x1 - b * y1 - c * z1);
}

Plane::Plane(Vector& pt1, Vector& pt2, Vector& pt3) {

    mpq_class& x1 = pt1.getX();
    mpq_class& y1 = pt1.getY();
    mpq_class& z1 = pt1.getZ();
    mpq_class& x2 = pt2.getX();
    mpq_class& y2 = pt2.getY();
    mpq_class& z2 = pt2.getZ();
    mpq_class& x3 = pt3.getX();
    mpq_class& y3 = pt3.getY();
    mpq_class& z3 = pt3.getZ();

    mpq_class a1 = x2 - x1;
    mpq_class b1 = y2 - y1;
    mpq_class c1 = z2 - z1;
    mpq_class a2 = x3 - x1;
    mpq_class b2 = y3 - y1;
    mpq_class c2 = z3 - z1;
    this->a = b1 * c2 - b2 * c1;
    this->b = a2 * c1 - a1 * c2;
    this->c = a1 * b2 - b1 * a2;
    this->d = (- a * x1 - b * y1 - c * z1);
}

void Plane::normal(Vector& n) {

    n.x = a;
    n.y = b;
    n.z = c;
}
