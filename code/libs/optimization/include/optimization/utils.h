#pragma once

#include <utils/utils.h>

// if the element at (row, col) is above the diagonal, it is skipped
inline void addMTripletToList_ignoreUpperElements(
        std::vector<SMTriplet> &triplets, int row, int col, double val) {
    if (row >= col) triplets.push_back(SMTriplet(row, col, val));
}

// if the element at (row, col) is above the diagonal, it is reflected on the
// lower diagonal
inline void addMTripletToList_reflectUpperElements(
        Array<SMTriplet> &triplets, int row, int col, double val) {
    if (IS_ZERO(val) && row != col)
        return;
    if (row >= col)
        triplets.push_back(SMTriplet(row, col, val));
    else
        triplets.push_back(SMTriplet(col, row, val));
}

// the element at (row, col) is mirrored, so it will be written symmetrically
// above and below the diagonal
inline void addMTripletToList_mirror(Array<SMTriplet> &triplets, int row,
                                     int col, double val) {
    if (row == col) {
        triplets.push_back(SMTriplet(row, col, val));
    } else {
        triplets.push_back(SMTriplet(row, col, val));
        triplets.push_back(SMTriplet(col, row, val));
    }
}

// write out the element as it comes
inline void addMTripletToList(Array<SMTriplet> &triplets, int row, int col,
                              double val) {
    triplets.push_back(SMTriplet(row, col, val));
}

#define ADD_HES_ELEMENT(list, i, j, v) \
    addMTripletToList_reflectUpperElements(list, i, j, (v))

