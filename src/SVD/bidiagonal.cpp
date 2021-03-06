/*************************************************************************
Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.

Contributors:
    * Sergey Bochkanov (ALGLIB project). Translation from FORTRAN to
      pseudocode.

See subroutines comments for additional copyrights.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer listed
  in this license in the documentation and/or other materials
  provided with the distribution.

- Neither the name of the copyright holders nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*************************************************************************/

//#include "stdafx.h"
#include "bidiagonal.h"

/*************************************************************************
Reduction of a rectangular matrix to  bidiagonal form

The algorithm reduces the rectangular matrix A to  bidiagonal form by
orthogonal transformations P and Q: A = Q*B*P.

Input parameters:
    A       -   source matrix. Array whose indexes range within [1..N, 1..N].
    M       -   number of rows in matrix A.
    N       -   number of columns in matrix A.

Output parameters:
    A       -   matrices Q, B, P in compact form (see below).
    TauQ    -   scalar factors which are used to form matrix Q.
    TauP    -   scalar factors which are used to form matrix P.

The main diagonal and one of the secondary diagonals of matrix A are
replaced with bidiagonal matrix B. Other elements contain elementary
reflections which form MxM matrix Q and NxN matrix P, respectively.

If M>=N, B is the upper bidiagonal MxN matrix and is stored in the
corresponding  elements  of matrix A. Matrix Q is represented as a product
of elementary reflections Q = H(1)*H(2)*...*H(n), where H(i) = 1 - tau*v*v'.
Here tau is a scalar which is stored in TauQ[i], and vector v has the
following structure: v(1:i-1)=0, v(i)=1, v(i+1:m) is stored in elements
A(i+1:m,i). Matrix P is as follows: P = G(1)*G(2)*...*G(n-1), where
G(i) = 1 - tau*u*u'. Tau is stored in TauP[i], u(1:i)=0, u(i+1)=1,  u(i+2:n)
is stored in elements A(i,i+2:n).

If M<N, B is the lower bidiagonal MxN matrix and is stored in the
corresponding elements of matrix A. Q = H(1)*H(2)*...*H(m-1), where
H(i) = 1 - tau*v*v',  tau  is stored in TauQ, v(1:i)=0, v(i+1)=1, v(i+2:m)
is stored in elements A(i+1:m,i).  P = G(1)*G(2)*...*G(m),  G(i) = 1 - tau*u*u',
tau is stored in TauP,  u(1:i-1)=0, u(i)=1, u(i+1:n) is stored A(i,i+1:n).

EXAMPLE:

m=6, n=5 (m > n):               m=5, n=6 (m < n):

(  d   e   u1  u1  u1 )         (  d   u1  u1  u1  u1  u1 )
(  v1  d   e   u2  u2 )         (  e   d   u2  u2  u2  u2 )
(  v1  v2  d   e   u3 )         (  v1  e   d   u3  u3  u3 )
(  v1  v2  v3  d   e  )         (  v1  v2  e   d   u4  u4 )
(  v1  v2  v3  v4  d  )         (  v1  v2  v3  e   d   u5 )
(  v1  v2  v3  v4  v5 )

Here vi and ui are vectors which form H(i) and G(i), and d and e -
are the diagonal and off-diagonal elements of matrix B.
*************************************************************************/
void tobidiagonal(ap::real_2d_array& a,
     int m,
     int n,
     ap::real_1d_array& tauq,
     ap::real_1d_array& taup)
{
    ap::real_1d_array work;
    ap::real_1d_array t;
    int minmn;
    int maxmn;
    int i;
    double ltau;
    int mmip1;
    int nmi;
    int ip1;
    int nmip1;
    int mmi;

    minmn = ap::minint(m, n);
    maxmn = ap::maxint(m, n);
    work.setbounds(1, maxmn);
    t.setbounds(1, maxmn);
    taup.setbounds(1, minmn);
    tauq.setbounds(1, minmn);
    if( m>=n )
    {
        
        //
        // Reduce to upper bidiagonal form
        //
        for(i = 1; i <= n; i++)
        {
            
            //
            // Generate elementary reflector H(i) to annihilate A(i+1:m,i)
            //
            mmip1 = m-i+1;
            ap::vmove(t.getvector(1, mmip1), a.getcolumn(i, i, m));
            generatereflection(t, mmip1, ltau);
            tauq(i) = ltau;
            ap::vmove(a.getcolumn(i, i, m), t.getvector(1, mmip1));
            t(1) = 1;
            
            //
            // Apply H(i) to A(i:m,i+1:n) from the left
            //
            applyreflectionfromtheleft(a, ltau, t, i, m, i+1, n, work);
            if( i<n )
            {
                
                //
                // Generate elementary reflector G(i) to annihilate
                // A(i,i+2:n)
                //
                nmi = n-i;
                ip1 = i+1;
                ap::vmove(t.getvector(1, nmi), a.getrow(i, ip1, n));
                generatereflection(t, nmi, ltau);
                taup(i) = ltau;
                ap::vmove(a.getrow(i, ip1, n), t.getvector(1, nmi));
                t(1) = 1;
                
                //
                // Apply G(i) to A(i+1:m,i+1:n) from the right
                //
                applyreflectionfromtheright(a, ltau, t, i+1, m, i+1, n, work);
            }
            else
            {
                taup(i) = 0;
            }
        }
    }
    else
    {
        
        //
        // Reduce to lower bidiagonal form
        //
        for(i = 1; i <= m; i++)
        {
            
            //
            // Generate elementary reflector G(i) to annihilate A(i,i+1:n)
            //
            nmip1 = n-i+1;
            ap::vmove(t.getvector(1, nmip1), a.getrow(i, i, n));
            generatereflection(t, nmip1, ltau);
            taup(i) = ltau;
            ap::vmove(a.getrow(i, i, n), t.getvector(1, nmip1));
            t(1) = 1;
            
            //
            // Apply G(i) to A(i+1:m,i:n) from the right
            //
            applyreflectionfromtheright(a, ltau, t, i+1, m, i, n, work);
            if( i<m )
            {
                
                //
                // Generate elementary reflector H(i) to annihilate
                // A(i+2:m,i)
                //
                mmi = m-i;
                ip1 = i+1;
                ap::vmove(t.getvector(1, mmi), a.getcolumn(i, ip1, m));
                generatereflection(t, mmi, ltau);
                tauq(i) = ltau;
                ap::vmove(a.getcolumn(i, ip1, m), t.getvector(1, mmi));
                t(1) = 1;
                
                //
                // Apply H(i) to A(i+1:m,i+1:n) from the left
                //
                applyreflectionfromtheleft(a, ltau, t, i+1, m, i+1, n, work);
            }
            else
            {
                tauq(i) = 0;
            }
        }
    }
}


/*************************************************************************
Unpacking matrix Q which reduces a matrix to bidiagonal form.

Input parameters:
    QP          -   matrices Q and P in compact form.
                    Output of ToBidiagonal subroutine.
    M           -   number of rows in matrix A.
    N           -   number of columns in matrix A.
    TAUQ        -   scalar factors which are used to form Q.
                    Output of ToBidiagonal subroutine.
    QColumns    -   required number of columns in matrix Q.
                    M>=QColumns>=0.

Output parameters:
    Q           -   first QColumns columns of matrix Q.
                    Array whose indexes range within [1..M, 1..QColumns].
                    If QColumns=0, the array is not modified.

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************/
void unpackqfrombidiagonal(const ap::real_2d_array& qp,
     int m,
     int n,
     const ap::real_1d_array& tauq,
     int qcolumns,
     ap::real_2d_array& q)
{
    int i;
    int j;
    int ip1;
    ap::real_1d_array v;
    ap::real_1d_array work;
    int vm;

    ap::ap_error::make_assertion(qcolumns<=m);
    if( m==0||n==0||qcolumns==0 )
    {
        return;
    }
    
    //
    // init
    //
    q.setbounds(1, m, 1, qcolumns);
    v.setbounds(1, m);
    work.setbounds(1, qcolumns);
    
    //
    // prepare Q
    //
    for(i = 1; i <= m; i++)
    {
        for(j = 1; j <= qcolumns; j++)
        {
            if( i==j )
            {
                q(i,j) = 1;
            }
            else
            {
                q(i,j) = 0;
            }
        }
    }
    if( m>=n )
    {
        for(i = ap::minint(n, qcolumns); i >= 1; i--)
        {
            vm = m-i+1;
            ap::vmove(v.getvector(1, vm), qp.getcolumn(i, i, m));
            v(1) = 1;
            applyreflectionfromtheleft(q, tauq(i), v, i, m, 1, qcolumns, work);
        }
    }
    else
    {
        for(i = ap::minint(m-1, qcolumns-1); i >= 1; i--)
        {
            vm = m-i;
            ip1 = i+1;
            ap::vmove(v.getvector(1, vm), qp.getcolumn(i, ip1, m));
            v(1) = 1;
            applyreflectionfromtheleft(q, tauq(i), v, i+1, m, 1, qcolumns, work);
        }
    }
}


/*************************************************************************
Multiplication by matrix Q which reduces matrix A to  bidiagonal form.

The algorithm allows pre- or post-multiply by Q or Q'.

Input parameters:
    QP          -   matrices Q and P in compact form.
                    Output of ToBidiagonal subroutine.
    M           -   number of rows in matrix A.
    N           -   number of columns in matrix A.
    TAUQ        -   scalar factors which are used to form Q.
                    Output of ToBidiagonal subroutine.
    Z           -   multiplied matrix.
                    array whose indexes range within [1..ZRows, 1..ZColumns].
    ZRows       -   number of rows in matrix Z. If FromTheRight=False,
                    ZRows=M, otherwise ZRows can be arbitrary.
    ZColumns    -   number of columns in matrix Z. If FromTheRight=True,
                    ZColumns=M, otherwise ZColumns can be arbitrary.
    FromTheRight -  pre- or post-multiply.
    DoTranspose -   multiply by Q or Q'.

Output parameters:
    Z           -   product of Z and Q.
                    Array whose indexes range within [1..ZRows, 1..ZColumns].
                    If ZRows=0 or ZColumns=0, the array is not modified.

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************/
void multiplybyqfrombidiagonal(const ap::real_2d_array& qp,
     int m,
     int n,
     const ap::real_1d_array& tauq,
     ap::real_2d_array& z,
     int zrows,
     int zcolumns,
     bool fromtheright,
     bool dotranspose)
{
    int i;
    int ip1;
    int i1;
    int i2;
    int istep;
    ap::real_1d_array v;
    ap::real_1d_array work;
    int vm;
    int mx;

    if( m<=0||n<=0||zrows<=0||zcolumns<=0 )
    {
        return;
    }
    ap::ap_error::make_assertion(fromtheright&&zcolumns==m||!fromtheright&&zrows==m);
    
    //
    // init
    //
    mx = ap::maxint(m, n);
    mx = ap::maxint(mx, zrows);
    mx = ap::maxint(mx, zcolumns);
    v.setbounds(1, mx);
    work.setbounds(1, mx);
    if( m>=n )
    {
        
        //
        // setup
        //
        if( fromtheright )
        {
            i1 = 1;
            i2 = n;
            istep = +1;
        }
        else
        {
            i1 = n;
            i2 = 1;
            istep = -1;
        }
        if( dotranspose )
        {
            i = i1;
            i1 = i2;
            i2 = i;
            istep = -istep;
        }
        
        //
        // Process
        //
        i = i1;
        do
        {
            vm = m-i+1;
            ap::vmove(v.getvector(1, vm), qp.getcolumn(i, i, m));
            v(1) = 1;
            if( fromtheright )
            {
                applyreflectionfromtheright(z, tauq(i), v, 1, zrows, i, m, work);
            }
            else
            {
                applyreflectionfromtheleft(z, tauq(i), v, i, m, 1, zcolumns, work);
            }
            i = i+istep;
        }
        while(i!=i2+istep);
    }
    else
    {
        
        //
        // setup
        //
        if( fromtheright )
        {
            i1 = 1;
            i2 = m-1;
            istep = +1;
        }
        else
        {
            i1 = m-1;
            i2 = 1;
            istep = -1;
        }
        if( dotranspose )
        {
            i = i1;
            i1 = i2;
            i2 = i;
            istep = -istep;
        }
        
        //
        // Process
        //
        if( m-1>0 )
        {
            i = i1;
            do
            {
                vm = m-i;
                ip1 = i+1;
                ap::vmove(v.getvector(1, vm), qp.getcolumn(i, ip1, m));
                v(1) = 1;
                if( fromtheright )
                {
                    applyreflectionfromtheright(z, tauq(i), v, 1, zrows, i+1, m, work);
                }
                else
                {
                    applyreflectionfromtheleft(z, tauq(i), v, i+1, m, 1, zcolumns, work);
                }
                i = i+istep;
            }
            while(i!=i2+istep);
        }
    }
}


/*************************************************************************
Unpacking matrix P which reduces matrix A to bidiagonal form.
The subroutine returns transposed matrix P.

Input parameters:
    QP      -   matrices Q and P in compact form.
                Output of ToBidiagonal subroutine.
    M       -   number of rows in matrix A.
    N       -   number of columns in matrix A.
    TAUP    -   scalar factors which are used to form P.
                Output of ToBidiagonal subroutine.
    PTRows  -   required number of rows of matrix P^T. N >= PTRows >= 0.

Output parameters:
    PT      -   first PTRows columns of matrix P^T
                Array whose indexes range within [1..PTRows, 1..N].
                If PTRows=0, the array is not modified.

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************/
void unpackptfrombidiagonal(const ap::real_2d_array& qp,
     int m,
     int n,
     const ap::real_1d_array& taup,
     int ptrows,
     ap::real_2d_array& pt)
{
    int i;
    int j;
    int ip1;
    ap::real_1d_array v;
    ap::real_1d_array work;
    int vm;

    ap::ap_error::make_assertion(ptrows<=n);
    if( m==0||n==0||ptrows==0 )
    {
        return;
    }
    
    //
    // init
    //
    pt.setbounds(1, ptrows, 1, n);
    v.setbounds(1, n);
    work.setbounds(1, ptrows);
    
    //
    // prepare PT
    //
    for(i = 1; i <= ptrows; i++)
    {
        for(j = 1; j <= n; j++)
        {
            if( i==j )
            {
                pt(i,j) = 1;
            }
            else
            {
                pt(i,j) = 0;
            }
        }
    }
    if( m>=n )
    {
        for(i = ap::minint(n-1, ptrows-1); i >= 1; i--)
        {
            vm = n-i;
            ip1 = i+1;
            ap::vmove(v.getvector(1, vm), qp.getrow(i, ip1, n));
            v(1) = 1;
            applyreflectionfromtheright(pt, taup(i), v, 1, ptrows, i+1, n, work);
        }
    }
    else
    {
        for(i = ap::minint(m, ptrows); i >= 1; i--)
        {
            vm = n-i+1;
            ap::vmove(v.getvector(1, vm), qp.getrow(i, i, n));
            v(1) = 1;
            applyreflectionfromtheright(pt, taup(i), v, 1, ptrows, i, n, work);
        }
    }
}


/*************************************************************************
Multiplication by matrix P which reduces matrix A to  bidiagonal form.

The algorithm allows pre- or post-multiply by P or P'.

Input parameters:
    QP          -   matrices Q and P in compact form.
                    Output of ToBidiagonal subroutine.
    M           -   number of rows in matrix A.
    N           -   number of columns in matrix A.
    TAUP        -   scalar factors which are used to form P.
                    Output of ToBidiagonal subroutine.
    Z           -   multiplied matrix.
                    Array whose indexes range within [1..ZRows, 1..ZColumns].
    ZRows       -   number of rows in matrix Z. If FromTheRight=False,
                    ZRows=N, otherwise ZRows can be arbitrary.
    ZColumns    -   number of columns in matrix Z. If FromTheRight=True,
                    ZColumns=N, otherwise ZColumns can be arbitrary.
    FromTheRight -  pre- or post-multiply.
    DoTranspose -   multiply by P or P'.

Output parameters:
    Z - product of Z and P.
                Array whose indexes range within [1..ZRows,1..ZColumns].
                If ZRows=0 or ZColumns=0, the array is not modified.

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************/
void multiplybypfrombidiagonal(const ap::real_2d_array& qp,
     int m,
     int n,
     const ap::real_1d_array& taup,
     ap::real_2d_array& z,
     int zrows,
     int zcolumns,
     bool fromtheright,
     bool dotranspose)
{
    int i;
    int ip1;
    ap::real_1d_array v;
    ap::real_1d_array work;
    int vm;
    int mx;
    int i1;
    int i2;
    int istep;

    if( m<=0||n<=0||zrows<=0||zcolumns<=0 )
    {
        return;
    }
    ap::ap_error::make_assertion(fromtheright&&zcolumns==n||!fromtheright&&zrows==n);
    
    //
    // init
    //
    mx = ap::maxint(m, n);
    mx = ap::maxint(mx, zrows);
    mx = ap::maxint(mx, zcolumns);
    v.setbounds(1, mx);
    work.setbounds(1, mx);
    v.setbounds(1, mx);
    work.setbounds(1, mx);
    if( m>=n )
    {
        
        //
        // setup
        //
        if( fromtheright )
        {
            i1 = n-1;
            i2 = 1;
            istep = -1;
        }
        else
        {
            i1 = 1;
            i2 = n-1;
            istep = +1;
        }
        if( !dotranspose )
        {
            i = i1;
            i1 = i2;
            i2 = i;
            istep = -istep;
        }
        
        //
        // Process
        //
        if( n-1>0 )
        {
            i = i1;
            do
            {
                vm = n-i;
                ip1 = i+1;
                ap::vmove(v.getvector(1, vm), qp.getrow(i, ip1, n));
                v(1) = 1;
                if( fromtheright )
                {
                    applyreflectionfromtheright(z, taup(i), v, 1, zrows, i+1, n, work);
                }
                else
                {
                    applyreflectionfromtheleft(z, taup(i), v, i+1, n, 1, zcolumns, work);
                }
                i = i+istep;
            }
            while(i!=i2+istep);
        }
    }
    else
    {
        
        //
        // setup
        //
        if( fromtheright )
        {
            i1 = m;
            i2 = 1;
            istep = -1;
        }
        else
        {
            i1 = 1;
            i2 = m;
            istep = +1;
        }
        if( !dotranspose )
        {
            i = i1;
            i1 = i2;
            i2 = i;
            istep = -istep;
        }
        
        //
        // Process
        //
        i = i1;
        do
        {
            vm = n-i+1;
            ap::vmove(v.getvector(1, vm), qp.getrow(i, i, n));
            v(1) = 1;
            if( fromtheright )
            {
                applyreflectionfromtheright(z, taup(i), v, 1, zrows, i, n, work);
            }
            else
            {
                applyreflectionfromtheleft(z, taup(i), v, i, n, 1, zcolumns, work);
            }
            i = i+istep;
        }
        while(i!=i2+istep);
    }
}


/*************************************************************************
Unpacking of the main and secondary diagonals of bidiagonal decomposition
of matrix A.

Input parameters:
    B   -   output of ToBidiagonal subroutine.
    M   -   number of rows in matrix B.
    N   -   number of columns in matrix B.

Output parameters:
    IsUpper -   True, if the matrix is upper bidiagonal.
                otherwise IsUpper is False.
    D       -   the main diagonal.
                Array whose index ranges within [1..Min(M,N)].
    E       -   the secondary diagonal (upper or lower, depending on
                the value of IsUpper).
                Array index ranges within [1..Min(M,N)], the last
                element is not used.

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************/
void unpackdiagonalsfrombidiagonal(const ap::real_2d_array& b,
     int m,
     int n,
     bool& isupper,
     ap::real_1d_array& d,
     ap::real_1d_array& e)
{
    int i;

    isupper = m>=n;
    if( m==0||n==0 )
    {
        return;
    }
    if( isupper )
    {
        d.setbounds(1, n);
        e.setbounds(1, n);
        for(i = 1; i <= n-1; i++)
        {
            d(i) = b(i,i);
            e(i) = b(i,i+1);
        }
        d(n) = b(n,n);
    }
    else
    {
        d.setbounds(1, m);
        e.setbounds(1, m);
        for(i = 1; i <= m-1; i++)
        {
            d(i) = b(i,i);
            e(i) = b(i+1,i);
        }
        d(m) = b(m,m);
    }
}



