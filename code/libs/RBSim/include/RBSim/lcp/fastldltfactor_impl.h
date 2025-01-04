#pragma once

#include <assert.h>

static void solveL1Stripe_2 (const double *L, double *B, unsigned int rowCount, unsigned int rowSkip);

template<unsigned int d_stride>
void scaleAndFactorizeL1Stripe_2(double *ARow, double *d, unsigned int rowIndex, unsigned int rowSkip);

template<unsigned int d_stride>
inline void scaleAndFactorizeL1FirstRowStripe_2(double *ARow, double *d, unsigned int rowSkip);

static void solveStripeL1_1 (const double *L, double *B, unsigned rowCount, unsigned int rowSkip);

template<unsigned int d_stride>
void scaleAndFactorizeL1Stripe_1(double *ARow, double *d, unsigned rowIndex);

template<unsigned int d_stride>
inline void scaleAndFactorizeL1FirstRowStripe_1(double *ARow, double *d);


template<unsigned int d_stride>
void factorMatrixAsLDLT(double *A, double *d, unsigned rowCount, unsigned rowSkip)
{
    if (rowCount < 1) return;

    double *ARow = A;
    unsigned blockStartRow = 0;

    const unsigned blockStep = 2;
    const unsigned lastRowIndex = rowCount >= blockStep ? rowCount - blockStep + 1 : 0;

    /* compute blocks of 2 rows */
    bool subsequentPass = false;
    for (; blockStartRow < lastRowIndex; subsequentPass = true, ARow += blockStep * rowSkip, blockStartRow += blockStep) 
    {
        if (subsequentPass)
        {
            /* solve L*(D*l)=a, l is scaled elements in 2 x i block at A(i,0) */
            solveL1Stripe_2(A, ARow, blockStartRow, rowSkip);
            scaleAndFactorizeL1Stripe_2<d_stride>(ARow, d, blockStartRow, rowSkip);
        }
        else
        {
            scaleAndFactorizeL1FirstRowStripe_2<d_stride>(ARow, d, rowSkip);
        }
        assert(blockStep == 2);
        /* done factorizing 2 x 2 block */
    }

    /* compute the (less than 2) rows at the bottom */
    if (!subsequentPass || blockStartRow == lastRowIndex)
    {
        assert(blockStep == 2); // for the blockStartRow == lastRowIndex comparison above

        if (subsequentPass)
        {
            solveStripeL1_1(A, ARow, blockStartRow, rowSkip);
            scaleAndFactorizeL1Stripe_1<d_stride>(ARow, d, blockStartRow);
        }
        else
        {
            scaleAndFactorizeL1FirstRowStripe_1<d_stride>(ARow, d);
        }
        assert(blockStep == 2);
        /* done factorizing 1 x 1 block */
    }
}

/* solve L*X=B, with B containing 2 right hand sides.
 * L is an n*n lower triangular matrix with ones on the diagonal.
 * L is stored by rows and its leading dimension is rowSkip.
 * B is an n*2 matrix that contains the right hand sides.
 * B is stored by columns and its leading dimension is also rowSkip.
 * B is overwritten with X.
 * this processes blocks of 2*2.
 * if this is in the factorizer source file, n must be a multiple of 2.
 */
static 
void solveL1Stripe_2(const double *L, double *B, unsigned rowCount, unsigned rowSkip)
{
    assert(rowCount != 0);
    assert(rowCount % 2 == 0);

    /* compute all 2 x 2 blocks of X */
    unsigned blockStartRow = 0;
    for (bool exitLoop = false, subsequentPass = false; !exitLoop; subsequentPass = true, exitLoop = (blockStartRow += 2) == rowCount) 
    {
        const double *ptrLElement;
        double *ptrBElement;

        /* declare variables - Z matrix */
        double Z11, Z12, Z21, Z22;

        /* compute all 2 x 2 block of X, from rows i..i+2-1 */
        if (subsequentPass)
        {
            ptrLElement = L + blockStartRow * rowSkip;
            ptrBElement = B;

            /* set Z matrix to 0 */
            Z11 = 0; Z12 = 0; Z21 = 0; Z22 = 0;

            /* the inner loop that computes outer products and adds them to Z */
            // The iteration starts with even number and decreases it by 2. So, it must end in zero
            for (unsigned columnCounter = blockStartRow; ;) 
            {
                /* declare p and q vectors, etc */
                double p1, q1, p2, q2;

                /* compute outer product and add it to the Z matrix */
                p1 = ptrLElement[0];
                q1 = ptrBElement[0];
                Z11 += p1 * q1;
                q2 = ptrBElement[rowSkip];
                Z12 += p1 * q2;
                p2 = ptrLElement[rowSkip];
                Z21 += p2 * q1;
                Z22 += p2 * q2;

                /* compute outer product and add it to the Z matrix */
                p1 = ptrLElement[1];
                q1 = ptrBElement[1];
                Z11 += p1 * q1;
                q2 = ptrBElement[1 + rowSkip];
                Z12 += p1 * q2;
                p2 = ptrLElement[1 + rowSkip];
                Z21 += p2 * q1;
                Z22 += p2 * q2;

                if (columnCounter > 6)
                {
                    columnCounter -= 6;

                    /* advance pointers */
                    ptrLElement += 6;
                    ptrBElement += 6;

                    /* compute outer product and add it to the Z matrix */
                    p1 = ptrLElement[-4];
                    q1 = ptrBElement[-4];
                    Z11 += p1 * q1;
                    q2 = ptrBElement[-4 + rowSkip];
                    Z12 += p1 * q2;
                    p2 = ptrLElement[-4 + rowSkip];
                    Z21 += p2 * q1;
                    Z22 += p2 * q2;

                    /* compute outer product and add it to the Z matrix */
                    p1 = ptrLElement[-3];
                    q1 = ptrBElement[-3];
                    Z11 += p1 * q1;
                    q2 = ptrBElement[-3 + rowSkip];
                    Z12 += p1 * q2;
                    p2 = ptrLElement[-3 + rowSkip];
                    Z21 += p2 * q1;
                    Z22 += p2 * q2;

                    /* compute outer product and add it to the Z matrix */
                    p1 = ptrLElement[-2];
                    q1 = ptrBElement[-2];
                    Z11 += p1 * q1;
                    q2 = ptrBElement[-2 + rowSkip];
                    Z12 += p1 * q2;
                    p2 = ptrLElement[-2 + rowSkip];
                    Z21 += p2 * q1;
                    Z22 += p2 * q2;

                    /* compute outer product and add it to the Z matrix */
                    p1 = ptrLElement[-1];
                    q1 = ptrBElement[-1];
                    Z11 += p1 * q1;
                    q2 = ptrBElement[-1 + rowSkip];
                    Z12 += p1 * q2;
                    p2 = ptrLElement[-1 + rowSkip];
                    Z21 += p2 * q1;
                    Z22 += p2 * q2;
                }
                else
                {
                    /* advance pointers */
                    ptrLElement += 2;
                    ptrBElement += 2;

                    if ((columnCounter -= 2) == 0)
                    {
                        break;
                    }
                }
                /* end of inner loop */
            }
        }
        else
        {
            ptrLElement = L/* + blockStartRow * rowSkip*/; assert(blockStartRow == 0);
            ptrBElement = B;

            /* set Z matrix to 0 */
            Z11 = 0; Z12 = 0; Z21 = 0; Z22 = 0;
        }

        /* finish computing the X(i) block */
        
        double Y11 = ptrBElement[0] - Z11;
        double Y12 = ptrBElement[rowSkip] - Z12;

        double p2 = ptrLElement[rowSkip];

        ptrBElement[0] = Y11;
        ptrBElement[rowSkip] = Y12;

        double Y21 = ptrBElement[1] - Z21 - p2 * Y11;
        double Y22 = ptrBElement[1 + rowSkip] - Z22 - p2 * Y12;

        ptrBElement[1] = Y21;
        ptrBElement[1 + rowSkip] = Y22;
        /* end of outer loop */
    }
}

template<unsigned int d_stride>
void scaleAndFactorizeL1Stripe_2(double *ARow, double *d, unsigned factorizationRow, unsigned rowSkip)
{
    assert(factorizationRow != 0);
    assert(factorizationRow % 2 == 0);

    double *ptrAElement = ARow;
    double *ptrDElement = d;

    /* scale the elements in a 2 x i block at A(i,0), and also */
    /* compute Z = the outer product matrix that we'll need. */
    double Z11 = 0, Z21 = 0, Z22 = 0;

    for (unsigned columnCounter = factorizationRow; ; ) 
    {
        double p1, q1, p2, q2, dd;

        p1 = ptrAElement[0];
        p2 = ptrAElement[rowSkip];
        dd = ptrDElement[0 * d_stride];
        q1 = p1 * dd;
        q2 = p2 * dd;
        ptrAElement[0] = q1;
        ptrAElement[rowSkip] = q2;
        Z11 += p1 * q1;
        Z21 += p2 * q1;
        Z22 += p2 * q2;

        p1 = ptrAElement[1];
        p2 = ptrAElement[1 + rowSkip];
        dd = ptrDElement[1 * d_stride];
        q1 = p1 * dd;
        q2 = p2 * dd;
        ptrAElement[1] = q1;
        ptrAElement[1 + rowSkip] = q2;
        Z11 += p1 * q1;
        Z21 += p2 * q1;
        Z22 += p2 * q2;

        if (columnCounter > 6)
        {
            columnCounter -= 6;

            ptrAElement += 6;
            ptrDElement += 6 * d_stride;

            p1 = ptrAElement[-4];
            p2 = ptrAElement[-4 + rowSkip];
            dd = ptrDElement[-4 * (int)d_stride];
            q1 = p1 * dd;
            q2 = p2 * dd;
            ptrAElement[-4] = q1;
            ptrAElement[-4 + rowSkip] = q2;
            Z11 += p1 * q1;
            Z21 += p2 * q1;
            Z22 += p2 * q2;

            p1 = ptrAElement[-3];
            p2 = ptrAElement[-3 + rowSkip];
            dd = ptrDElement[-3 * (int)d_stride];
            q1 = p1 * dd;
            q2 = p2 * dd;
            ptrAElement[-3] = q1;
            ptrAElement[-3 + rowSkip] = q2;
            Z11 += p1 * q1;
            Z21 += p2 * q1;
            Z22 += p2 * q2;

            p1 = ptrAElement[-2];
            p2 = ptrAElement[-2 + rowSkip];
            dd = ptrDElement[-2 * (int)d_stride];
            q1 = p1 * dd;
            q2 = p2 * dd;
            ptrAElement[-2] = q1;
            ptrAElement[-2 + rowSkip] = q2;
            Z11 += p1 * q1;
            Z21 += p2 * q1;
            Z22 += p2 * q2;

            p1 = ptrAElement[-1];
            p2 = ptrAElement[-1 + rowSkip];
            dd = ptrDElement[-1 * (int)d_stride];
            q1 = p1 * dd;
            q2 = p2 * dd;
            ptrAElement[-1] = q1;
            ptrAElement[-1 + rowSkip] = q2;
            Z11 += p1 * q1;
            Z21 += p2 * q1;
            Z22 += p2 * q2;
        }
        else
        {
            ptrAElement += 2;
            ptrDElement += 2 * d_stride;

            if ((columnCounter -= 2) == 0)
            {
                break;
            }
        }
    }

    /* solve for diagonal 2 x 2 block at A(i,i) */
    double Y11 = ptrAElement[0] - Z11;
    double Y21 = ptrAElement[rowSkip] - Z21;
    double Y22 = ptrAElement[1 + rowSkip] - Z22;

    /* factorize 2 x 2 block Y, ptrDElement */
    /* factorize row 1 */
    double dd = 1.0/Y11;

    ptrDElement[0 * d_stride] = dd;
    assert(ptrDElement == d + (size_t)factorizationRow * d_stride);

    /* factorize row 2 */
    double q2 = Y21 * dd;
    ptrAElement[rowSkip] = q2;

    double sum = Y21 * q2;
    ptrDElement[1 * d_stride] = 1.0/(Y22 - sum);
}

template<unsigned int d_stride>
void scaleAndFactorizeL1FirstRowStripe_2(double *ARow, double *d, unsigned rowSkip)
{
    double *ptrAElement = ARow;
    double *ptrDElement = d;

    /* solve for diagonal 2 x 2 block at A(0,0) */
    double Y11 = ptrAElement[0]/* - Z11*/;
    double Y21 = ptrAElement[rowSkip]/* - Z21*/;
    double Y22 = ptrAElement[1 + rowSkip]/* - Z22*/;

    /* factorize 2 x 2 block Y, ptrDElement */
    /* factorize row 1 */
    double dd = 1.0/(Y11);

    ptrDElement[0 * d_stride] = dd;
    assert(ptrDElement == d/* + (size_t)factorizationRow * d_stride*/);

    /* factorize row 2 */
    double q2 = Y21 * dd;
    ptrAElement[rowSkip] = q2;

    double sum = Y21 * q2;
    ptrDElement[1 * d_stride] = 1.0/(Y22 - sum);
}


/* solve L*X=B, with B containing 1 right hand sides.
 * L is an n*n lower triangular matrix with ones on the diagonal.
 * L is stored by rows and its leading dimension is lskip.
 * B is an n*1 matrix that contains the right hand sides.
 * B is stored by columns and its leading dimension is also lskip.
 * B is overwritten with X.
 * this processes blocks of 2*2.
 * if this is in the factorizer source file, n must be a multiple of 2.
 */
static 
void solveStripeL1_1(const double *L, double *B, unsigned rowCount, unsigned rowSkip)
{
    assert(rowCount != 0);
    assert(rowCount % 2 == 0);

    /* compute all 2 x 1 blocks of X */
    unsigned blockStartRow = 0;
    for (bool exitLoop = false, subsequentPass = false; !exitLoop; subsequentPass = true, exitLoop = (blockStartRow += 2) == rowCount) 
    {
        const double *ptrLElement;
        double *ptrBElement;

        /* declare variables - Z matrix */
        double Z11, Z21;

        if (subsequentPass)
        {
            ptrLElement = L + (size_t)blockStartRow * rowSkip;
            ptrBElement = B;

            /* set the Z matrix to 0 */
            Z11 = 0; Z21 = 0;

            /* compute all 2 x 1 block of X, from rows i..i+2-1 */
            
            /* the inner loop that computes outer products and adds them to Z */
            // The iteration starts with even number and decreases it by 2. So, it must end in zero
            for (unsigned columnCounter = blockStartRow; ; ) 
            {
                /* declare p and q vectors, etc */
                double p1, q1, p2;

                /* compute outer product and add it to the Z matrix */
                p1 = ptrLElement[0];
                q1 = ptrBElement[0];
                Z11 += p1 * q1;
                p2 = ptrLElement[rowSkip];
                Z21 += p2 * q1;
                
                /* compute outer product and add it to the Z matrix */
                p1 = ptrLElement[1];
                q1 = ptrBElement[1];
                Z11 += p1 * q1;
                p2 = ptrLElement[1 + rowSkip];
                Z21 += p2 * q1;

                if (columnCounter > 6)
                {
                    columnCounter -= 6;

                    /* advance pointers */
                    ptrLElement += 6;
                    ptrBElement += 6;

                    /* compute outer product and add it to the Z matrix */
                    p1 = ptrLElement[-4];
                    q1 = ptrBElement[-4];
                    Z11 += p1 * q1;
                    p2 = ptrLElement[-4 + rowSkip];
                    Z21 += p2 * q1;

                    /* compute outer product and add it to the Z matrix */
                    p1 = ptrLElement[-3];
                    q1 = ptrBElement[-3];
                    Z11 += p1 * q1;
                    p2 = ptrLElement[-3 + rowSkip];
                    Z21 += p2 * q1;

                    /* compute outer product and add it to the Z matrix */
                    p1 = ptrLElement[-2];
                    q1 = ptrBElement[-2];
                    Z11 += p1 * q1;
                    p2 = ptrLElement[-2 + rowSkip];
                    Z21 += p2 * q1;

                    /* compute outer product and add it to the Z matrix */
                    p1 = ptrLElement[-1];
                    q1 = ptrBElement[-1];
                    Z11 += p1 * q1;
                    p2 = ptrLElement[-1 + rowSkip];
                    Z21 += p2 * q1;
                }
                else
                {
                    /* advance pointers */
                    ptrLElement += 2;
                    ptrBElement += 2;

                    if ((columnCounter -= 2) == 0)
                    {
                        break;
                    }
                }
                /* end of inner loop */
            }
        }
        else
        {
            ptrLElement = L/* + (size_t)blockStartRow * rowSkip*/; assert(blockStartRow == 0);
            ptrBElement = B;

            /* set the Z matrix to 0 */
            Z11 = 0; Z21 = 0;
        }
        
        /* finish computing the X(i) block */
        double p2 = ptrLElement[rowSkip];

        double Y11 = ptrBElement[0] - Z11;
        double Y21 = ptrBElement[1] - Z21 - p2 * Y11;

        ptrBElement[0] = Y11;
        ptrBElement[1] = Y21;
        /* end of outer loop */
    }
}

template<unsigned int d_stride>
void scaleAndFactorizeL1Stripe_1(double *ARow, double *d, unsigned factorizationRow)
{
    double *ptrAElement = ARow;
    double *ptrDElement = d;

    /* scale the elements in a 1 x i block at A(i,0), and also */
    /* compute Z = the outer product matrix that we'll need. */
    double Z11 = 0, Z22 = 0;

    for (unsigned columnCounter = factorizationRow; ; ) 
    {
        double p1, p2, q1, q2, dd1, dd2;

        p1 = ptrAElement[0];
        p2 = ptrAElement[1];
        dd1 = ptrDElement[0 * d_stride];
        dd2 = ptrDElement[1 * d_stride];
        q1 = p1 * dd1;
        q2 = p2 * dd2;
        ptrAElement[0] = q1;
        ptrAElement[1] = q2;
        Z11 += p1 * q1;
        Z22 += p2 * q2;

        if (columnCounter > 6)
        {
            columnCounter -= 6;

            ptrAElement += 6;
            ptrDElement += 6 * d_stride;

            p1 = ptrAElement[-4];
            p2 = ptrAElement[-3];
            dd1 = ptrDElement[-4 * (int)d_stride];
            dd2 = ptrDElement[-3 * (int)d_stride];
            q1 = p1 * dd1;
            q2 = p2 * dd2;
            ptrAElement[-4] = q1;
            ptrAElement[-3] = q2;
            Z11 += p1 * q1;
            Z22 += p2 * q2;

            p1 = ptrAElement[-2];
            p2 = ptrAElement[-1];
            dd1 = ptrDElement[-2 * (int)d_stride];
            dd2 = ptrDElement[-1 * (int)d_stride];
            q1 = p1 * dd1;
            q2 = p2 * dd2;
            ptrAElement[-2] = q1;
            ptrAElement[-1] = q2;
            Z11 += p1 * q1;
            Z22 += p2 * q2;
        }
        else
        {
            ptrAElement += 2;
            ptrDElement += 2 * d_stride;

            if ((columnCounter -= 2) == 0)
            {
                break;
            }
        }
    }

    double Y11 = ptrAElement[0] - (Z11 + Z22);

    /* solve for diagonal 1 x 1 block at A(i,i) */
    assert(ptrDElement == d + (size_t)factorizationRow * d_stride);
    /* factorize 1 x 1 block Y, ptrDElement */
    /* factorize row 1 */
    ptrDElement[0 * d_stride] = 1.0/(Y11);
}

template<unsigned int d_stride>
void scaleAndFactorizeL1FirstRowStripe_1(double *ARow, double *d)
{
    double *ptrAElement = ARow;
    double *ptrDElement = d;

    double Y11 = ptrAElement[0];

    /* solve for diagonal 1 x 1 block at A(0,0) */
    /* factorize 1 x 1 block Y, ptrDElement */
    /* factorize row 1 */
    ptrDElement[0 * d_stride] = 1.0/(Y11);
}
