//hierarchy of linear algebraic operations (hlao) [matrix computations (linear algebra or matrix algebra)]

//ECMAScript module

//todo:
//   - matrix_multiplication() return a vector not a matrix from [[]] x [[],[],[]] (1 x m) x (m x n) will return [x0, x1, x2, ...] see I:\code\spatial_v2\js\regression\ex1.html (at theta update)
//   - not the correct implementation of matrix rank 'matrix_rank()'

//matrix, rectangular array of numbers containing:
//   - m rows
//   - n columns
//      - matrix has order m x n
//   - if m = n the matrix is a square matrix

//column vector, m x 1
//row vector, 1 x n

function matrix_arithmetic(A,B,operator){
    var m = A.length;
    var n = A[0].length; //0 - first row
    var r = B.length;
    var s = B[0].length; //0 - first row
    
    assert((m == r)&&(n == s),'Assertion failed: matrix A and B are not of the same order (matrix arithmetic).');
    
    var C = zeros_matrix(m,n);
    
    //Element-wise operation.
    switch(operator){
        case '+':
            for(var i=0;i<m;i=i+1){
                for(var j=0;j<n;j=j+1){
                    C[i][j] = A[i][j] + B[i][j];
                }
            }
            break;
        
        case '-':
            for(var i=0;i<m;i=i+1){
                for(var j=0;j<n;j=j+1){
                    C[i][j] = A[i][j] - B[i][j];
                }
            }
            break;
            
        case '*':
            for(var i=0;i<m;i=i+1){
                for(var j=0;j<n;j=j+1){
                    C[i][j] = A[i][j] * B[i][j];
                }
            }
            break;
            
        default:
            assert(false,'Assertion failed: operator not found (matrix arithmetic).');
    }

    return C;
}

function matrix_multiplication_scalar(A,k){
    var m = A.length;
    var n = A[0].length; //0 - first row
    //if(!(typeof k[0].length == 'undefined')) k = k[0][0];
    
    var C = zeros_matrix(m,n);
    
    for(var i=0;i<m;i=i+1){
        for(var j=0;j<n;j=j+1){
            C[i][j] = k * A[i][j];
        }
    }
    
    return C;
}

function matrix_multiplication(A,B){
    var m = A.length;
    var n = A[0].length; //0 - first row
    var r = B.length;
    var s = B[0].length; //0 - first row
    
    //if A is m x n
    //if B is r x s
    assert((n == r),'Assertion failed: the number of columns in A not equal to the number of rows in B (matrix multiplication).');
    
    var C = zeros_matrix(m,s);
    
    for(var i=0;i<m;i=i+1){
        for(var j=0;j<s;j=j+1){
            for(var k=0;k<n;k=k+1){
                C[i][j] = C[i][j] + A[i][k]*B[k][j];
            }
        }
    }
    
    if(m == 1){ //send back a row vector
        var D = [];
        for(var i=0;i<s;i=i+1){
            D.push(C[0][i]);
        }
        var C = D;
    }
    
    return C;
}

/*
var a = [
    [ 1, 2, 3, 4, 5],
    [ 6, 7, 8, 9,10],
    [11,12,13,14,15]
];
matrix_mean(a);
//returns [[3],[8],[13]]
*/
//Average of each row (DIM = 1).
function matrix_mean(A,DIM){
    var dim = size(A);
    var m = dim[0];
    var n = dim[1];
    var avg = [];
    
    switch(DIM){
        case 1: //rows
            for(var i=0;i<m;i=i+1){ //Row.
                if(typeof avg[i] == 'undefined'){
                    avg[i] = [];
                    avg[i][0] = 0.0;
                }
                for(var j=0;j<n;j=j+1){//Column.
                    avg[i][0] = avg[i][0] + A[i][j];
                }
            }
            
            avg = matrix_multiplication_scalar(avg,1.0/n);
            break;
            
        case 2: //columns
            for(var i=0;i<n;i=i+1){ //Column.
                avg[i] = 0.0;
                for(var j=0;j<m;j=j+1){//Row.
                    avg[i] = avg[i] + A[j][i];
                }
                avg[i] = avg[i]/n;
            }
            
            break;
            
        default:
            //assert() see matrixAlgebra.js
            assert(false,'Assertion failed: unknown dimension specified in function matrix_mean().');
    }
    
    return avg;
}

//Basically a copy of 'matrix_mean()'
function matrix_summation(A,DIM){
    var dim = size(A);
    var m = dim[0];
    var n = dim[1];
    var tot = [];
    
    switch(DIM){
        case 1: //rows
            for(var i=0;i<m;i=i+1){ //Row.
                if(typeof tot[i] == 'undefined'){
                    tot[i] = [];
                    tot[i][0] = 0.0;
                }
                for(var j=0;j<n;j=j+1){//Column.
                    tot[i][0] = tot[i][0] + A[i][j];
                }
            }
            
            //tot = matrix_multiplication_scalar(tot,1.0/n);
            break;
            
        case 2: //columns
            for(var i=0;i<n;i=i+1){ //Column.
                tot[i] = 0.0;
                for(var j=0;j<m;j=j+1){//Row.
                    tot[i] = tot[i] + A[j][i];
                }
                //tot[i] = tot[i]/n;
            }
            
            break;
            
        default:
            //assert() see matrixAlgebra.js
            assert(false,'Assertion failed: unknown dimension specified in function matrix_summation().');
    }
    
    return tot;
}

function vector_dot(a,b){
    //convert 'a' and 'b' to row vectors if required
    
    //'a' vector - row or column vector?
    if(!(typeof a[0].length == 'undefined')){ //if a[0].length 'undefined' - row vector DON'T transpose
        a = vector_transpose(a); //convert column vector to row vector
    }
    var a_length = a.length;
    
    //'b' vector - row or column vector?
    if(!(typeof b[0].length == 'undefined')){ //if b[0].length 'undefined' - row vector DON'T transpose
        b = vector_transpose(b); //convert column vector to row vector
    }
    var b_length = b.length;
    
    //both vectors should be the same length
    assert((a_length == b_length),'Assertion failed: vector A and B are not of the same order (dot product).');
    
    var ab = 0.0;
    
    for(var i=0;i<a_length;i=i+1){
        ab = ab + a[i] * b[i];
    }
    
    return ab;
}

function vector_cross(a,b){ //vector cross product
    //convert 'a' and 'b' to row vectors if required
    
    //'a' vector - row or column vector?
    if(!(typeof a[0].length == 'undefined')){ //if a[0].length 'undefined' - row vector DON'T transpose
        a = vector_transpose(a); //convert column vector to row vector
    }
    var a_length = a.length;
    
    //'b' vector - row or column vector?
    if(!(typeof b[0].length == 'undefined')){ //if b[0].length 'undefined' - row vector DON'T transpose
        b = vector_transpose(b); //convert column vector to row vector
    }
    var b_length = b.length;
    
    //both vectors should be the same length
    assert((a_length == b_length),'Assertion failed: vector A and B are not of the same order (vector cross).');
    
    var axbi = a[1]*b[2] - a[2]*b[1];
    var axbj = a[2]*b[0] - a[0]*b[2];
    var axbk = a[0]*b[1] - a[1]*b[0];
    
    return [[axbi],[axbj],[axbk]];
}

function vector_multiplication_scalar(a,r){
    var a_length = a.length;
    
    //'a' vector - row or column vector?
    if(typeof a[0].length == 'undefined'){ //if a[0].length 'undefined'
        //row vector
        var m = 1;
        var n = a_length;
        var b = zeros_matrix(m,n);
        for(var i=0;i<a_length;i=i+1){
            b[i] = a[i] * r;
        }
    } else {
        //column vector
        var m = a_length;
        var n = 1;
        var b = zeros_matrix(m,n);
        for(var i=0;i<a_length;i=i+1){
            b[i][0] = a[i][0] * r;
        }
    }
    
    return b;
}

function matrix_transpose(A){
    var m = A.length;    //rows
    var n = A[0].length; //columns
    
    var B = undefined_matrix(n,m);
    
    for(var i=0;i<m;i=i+1){
        for(var j=0;j<n;j=j+1){
            B[j][i] = A[i][j];
        }
    }
    
    return B;
}

function vector_transpose(a){
    //is 'a' a row or column vector?
    if(typeof a[0].length == 'undefined'){ //a[0].length is undefined if a row vector
        //row to column vector
        var m = 1;           //rows
        var n = a.length;    //columns
        var b = undefined_matrix(n,m);
        for(var i=0;i<n;i=i+1){
            b[i][0] = a[i];
        }
    } else {
        //column to row vector
        var m = a.length;    //rows
        var n = a[0].length; //columns (should be 1)
        var b = undefined_matrix(n,m);
        for(var i=0;i<m;i=i+1){
            b[i] = a[i][0];
        }
    }
    
    return b;
}

function zeros_matrix(m,n){
    if(arguments.length == 1){ //m = [m,n] (sent via an array)
        n = m[1];
        m = m[0];
    }
    
    var A = undefined_matrix(m,n);
    
    for(var i=0;i<m;i=i+1){
        for(var j=0;j<n;j=j+1){
            A[i][j] = 0.0;
        }
    }
    
    return A;
}

function ones_matrix(m,n){
    var A = undefined_matrix(m,n);
    
    for(var i=0;i<m;i=i+1){
        for(var j=0;j<n;j=j+1){
            A[i][j] = 1.0;
        }
    }
    
    return A;
}

function random_matrix(m,n){
    var A = undefined_matrix(m,n);
    
    for(var i=0;i<m;i=i+1){
        for(var j=0;j<n;j=j+1){
            A[i][j] = Math.random();
        }
    }
    
    return A;
}

function identity_matrix(N){
    var A = zeros_matrix(N,N);
    
    for(var i=0;i<N;i=i+1){
        for(var j=0;j<N;j=j+1){
            if(i == j) A[i][j] = 1.0;
        }
    }
    
    return A;
}

function undefined_matrix(m,n){
    var A = new Array(m);
    for(var i=0;i<A.length;i=i+1){
        A[i] = new Array(n);
    }
    return A;
}

function zeros_vector(r,type){
    var u = [];
    
    switch(type){
        case 'col':
            for(var i=0;i<r;i=i+1){
                u.push([0.0]);
            }
            break;
        case 'row':
            for(var i=0;i<r;i=i+1){
                u.push(0.0);
            }
            break;
        default:
            assert(false,'Assertion failed: vector type not one of column (col) or row (row) in function zeros_vector().');
    }
    
    return u;
}

function ones_vector(r,type){
    var u = [];
    
    switch(type){
        case 'col':
            for(var i=0;i<r;i=i+1){
                u.push([1.0]);
            }
            break;
        case 'row':
            for(var i=0;i<r;i=i+1){
                u.push(1.0);
            }
            break;
        default:
            assert(false,'Assertion failed: vector type not one of column (col) or row (row) in function ones_vector().');
    }
    
    return u;
}

function matrix_rank(W){
    var rank = 0;
    var zero_TOL = 1e-6;
    
    for(var i=0;i<W.length;i=i+1){
        var diag = Math.abs(W[i]); //W[i] is the diagonal.
        if (diag > zero_TOL){
            rank = rank + 1;
        }
    }
    
    return rank;
}

//For a row vector the second dimension will be 'undefined'.
//   - Example:
/*
var dummyArray = [1,2,3]; //row vector
//var dummyArray = [[1],[2],[3]];
var dim = getDimensions(dummyArray);
console.log(dim);
Array [ 3, undefined ]
if(typeof(dim[1]) === 'undefined'){
   console.log('Array is a row vector.');
} else {
   console.log('Array is a column vector or matrix.');
}
*/
function getDimensions(A){
    return[A.length, A[0].length];
}

function skew(a){
    var x = a[0][0]; var y = a[1][0]; var z = a[2][0];
    return([
        [0.0,  -z,   y],
        [  z, 0.0,  -x],
        [ -y,   x, 0.0]
    ]);
}

function vex(S){

    //S = [
    //        [    0.0, -1.0*wz,      wy],
    //        [     wz,     0.0, -1.0*wx],
    //        [-1.0*wy,      wx,     0.0]
    //    ];
    
    var wx = 0.5*(S[2][1] - S[1][2]);
    var wy = 0.5*(S[0][2] - S[2][0]);
    var wz = 0.5*(S[1][0] - S[0][1]);
    
    return [[wx],[wy],[wz]];
}

//REF: http://stackoverflow.com/questions/15313418/javascript-assert
function assert(condition, message){
    if(!condition){
        throw message || "Assertion failed";
    }
}

export {
    matrix_arithmetic,
    matrix_multiplication_scalar,
    matrix_multiplication,
    matrix_mean,
    matrix_summation,
    vector_dot,
    vector_cross,
    vector_multiplication_scalar,
    matrix_transpose,
    vector_transpose,
    zeros_matrix,
    ones_matrix,
    random_matrix,
    identity_matrix,
    undefined_matrix,
    zeros_vector,
    ones_vector,
    matrix_rank,
    skew,
    vex
};