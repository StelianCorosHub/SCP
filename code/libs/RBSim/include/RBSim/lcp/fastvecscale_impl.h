/*
 * Vector scaling function implementation
 * Improvements and cooperative implementation copyright (c) 2017-2020 Oleh Derevenko, odar@eleks.com (change all "a" to "e")  
 */

#pragma once

#include <assert.h>

template<unsigned int a_stride, unsigned int d_stride>
void scaleLargeVector(double *aStart, const double *dStart, unsigned elementCount)
{
    assert (aStart && dStart && elementCount >= 0);
    
    const unsigned step = 4;

    double *ptrA = aStart;
    const double *ptrD = dStart;
    const double *const dStepsEnd = dStart + (size_t)(elementCount & ~(step - 1)) * d_stride;
    for (; ptrD != dStepsEnd; ptrA += step * a_stride, ptrD += step * d_stride) 
    {
        double a0 = ptrA[0], a1 = ptrA[1 * a_stride], a2 = ptrA[2 * a_stride], a3 = ptrA[3 * a_stride];
        double d0 = ptrD[0], d1 = ptrD[1 * d_stride], d2 = ptrD[2 * d_stride], d3 = ptrD[3 * d_stride];
        a0 *= d0;
        a1 *= d1;
        a2 *= d2;
        a3 *= d3;
        ptrA[0] = a0; ptrA[1 * a_stride] = a1; ptrA[2 * a_stride] = a2; ptrA[3 * a_stride] = a3;
        assert(step == 4);
    }

    switch (elementCount & (step - 1))
    {
        case 3:
        {
            double a2 = ptrA[2 * a_stride];
            double d2 = ptrD[2 * d_stride];
            ptrA[2 * a_stride] = a2 * d2;
            // break; -- proceed to case 2
        }

        case 2:
        {
            double a1 = ptrA[1 * a_stride];
            double d1 = ptrD[1 * d_stride];
            ptrA[1 * a_stride] = a1 * d1;
            // break; -- proceed to case 1
        }

        case 1:
        {
            double a0 = ptrA[0];
            double d0 = ptrD[0];
            ptrA[0] = a0 * d0;
            break;
        }
    }
    assert(step == 4);
}

