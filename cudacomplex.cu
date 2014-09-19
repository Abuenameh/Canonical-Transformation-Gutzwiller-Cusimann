/*
 * cudacomplex.cu
 *
 *  Created on: Sep 14, 2014
 *      Author: Abuenameh
 */

#include "cudacomplex.cuh"


//
// Non-member overloads for single complex
//

// subtract single complex from scalar
HOSTDEVICE _cudacomplex<float2, float> operator-(const float REF(a), const _cudacomplex<float2, float> REF(b)) {
    _cudacomplex<float2, float> result = {
        { a - b.value.x, -b.value.y}};
    return result;
}

// add single complex to scalar
HOSTDEVICE _cudacomplex<float2, float> operator+(const float REF(a), const _cudacomplex<float2, float> REF(b)) {
    _cudacomplex<float2, float> result = {
        { a + b.value.x, b.value.y}};
    return result;
}

// multiply scalar with single complex
HOSTDEVICE _cudacomplex<float2, float> operator*(const float REF(a), const _cudacomplex<float2, float> REF(b)) {
    _cudacomplex<float2, float> result = {
        { a * b.value.x, a * b.value.y}};
    return result;
}

// divide scalar by single complex
HOSTDEVICE _cudacomplex<float2, float> operator/(const float REF(a), const _cudacomplex<float2, float> REF(b)) {
    float tmp = (b.value.x * b.value.x + b.value.y * b.value.y);
    _cudacomplex<float2, float> result = {
        { (a * b.value.x) / tmp, (-a * b.value.y) / tmp}};
    return result;
}

//
// Non-member overloads for double complex
//

// subtract double complex from scalar
HOSTDEVICE _cudacomplex<double2, double> operator-(const double REF(a), const _cudacomplex<double2, double> REF(b)) {
    _cudacomplex<double2, double> result = {
        { a - b.value.x, -b.value.y}};
    return result;
}

// add double complex to scalar
HOSTDEVICE _cudacomplex<double2, double> operator+(const double REF(a), const _cudacomplex<double2, double> REF(b)) {
    _cudacomplex<double2, double> result = {
        { a + b.value.x, b.value.y}};
    return result;
}

// multiply scalar with double complex
HOSTDEVICE _cudacomplex<double2, double> operator*(const double REF(a), const _cudacomplex<double2, double> REF(b)) {
    _cudacomplex<double2, double> result = {
        { a * b.value.x, a * b.value.y}};
    return result;
}

// divide scalar by double complex
HOSTDEVICE _cudacomplex<double2, double> operator/(const double REF(a), const _cudacomplex<double2, double> REF(b)) {
    double tmp = (b.value.x * b.value.x + b.value.y * b.value.y);
    _cudacomplex<double2, double> result = {
        { (a * b.value.x) / tmp, (-a * b.value.y) / tmp}};
    return result;
}

// a possible alternative to a single complex constructor
HOSTDEVICE singlecomplex make_singlecomplex(float a, float b)
 {
    singlecomplex res;
    res.real() = a;
    res.imag() = b;
    return res;
}

// a possible alternative to a double complex constructor
HOSTDEVICE doublecomplex make_doublecomplex(double a, double b)
 {
    doublecomplex res;
    res.real() = a;
    res.imag() = b;
    return res;
}



