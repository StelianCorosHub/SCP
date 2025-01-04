#pragma once

template<unsigned b_stride>
double calculateLargeVectorDot (const double *a, const double *b, unsigned n) {
    double sum = 0;
    const double *a_end = a + (n & (int)(~3));
    for (; a != a_end; b += 4 * b_stride, a += 4) {
        double p0 = a[0], p1 = a[1], p2 = a[2], p3 = a[3];
        double q0 = b[0 * b_stride], q1 = b[1 * b_stride], q2 = b[2 * b_stride], q3 = b[3 * b_stride];
        double m0 = p0 * q0;
        double m1 = p1 * q1;
        double m2 = p2 * q2;
        double m3 = p3 * q3;
        sum += m0 + m1 + m2 + m3;
    }
    a_end += (n & 3);
    for (; a != a_end; b += b_stride, ++a) {
        sum += (*a) * (*b);
    }
    return sum;
}

