function swap(arr, idx1, idx2){
    const temp = arr[idx1];
    arr[idx1] = arr[idx2];
    arr[idx2] = temp;
}

/**
 * Performs an LUP decomposition on A, in-place.
 * After execution `A` will contain the matrices L,U.
 * @param {Number[][]} A An nxn matrix. Note that `A` is updated in-place.
 * @returns {[Number, Number[]]} `[pi, sign]` The permutation `pi` represented as an array, and `sign`, the parity of the permutation
 * 
 * @example
 * // Logs [1, [3,0,2,1],[[5,-4,-3,1],[0.4,5.6,2.2,-3.4],[0.8,0.9285714285714287,-2.642857142857143,7.357142857142858],[-0.2,-0.5,-0.9459459459459459,9.45945945945946]]]
 * let test1 = [[2,4,1,-3],[-1,-2,2,4],[4,2,-3,5],[5,-4,-3,1]];
 * const [pi, sign] = LUPDecomposition(test1);
 * console.log([sign, pi, test1]);
 */
function LUPDecomposition(A){
    const n = A.length;
    const pi = Array(n).fill().map((_,i) => i);
    let sign = 1;
    for(let i=0; i<n; i++){
        let pivot_val = 0;
        let pivot_idx = i;
        for(let j=i; j<n; j++){
            if(Math.abs(A[j][i]) > pivot_val){
                pivot_val = Math.abs(A[j][i]);
                pivot_idx = j;
            }
        }
        if(pivot_val === 0){ throw("Error: singular matrix"); }
        swap(A, i, pivot_idx);
        if(i != pivot_idx){
            swap(pi, i, pivot_idx);
            sign *= -1;
        }
        for(let j=i+1; j<n; j++){
            A[j][i] = A[j][i] / A[i][i];
            for(let k=i+1; k<n; k++){
                A[j][k] = A[j][k] - (A[j][i] * A[i][k]);
            }
        }
    }
    return [pi, sign];
}

function vectorSolverHelper(A, pi, b){
    const n = A.length;
    // Solve Ly=Pb for the unknown y=UX using forward substitution
    let y = Array(n).fill();
    for(let i=0; i<n; i++){
        const sum = y.slice(0, i)
            .reduce((a, y_j, j) => a + (y_j * A[i][j]), 0);
        y[i] = b[pi[i]] - sum;
    }
    // Solve Ux=y for the unknown x using back substitution
    let x = Array(n).fill();
    for(let i=n-1; i>=0; i--){
        const sum = x.slice(i+1)
            .reduce((a, x_j, j) => a + (x_j * A[i][j+i+1]), 0);
        x[i] = (y[i] - sum) / A[i][i];
    }
    return x;
}

/**
 * Solves Ax=b for the unknown vector x.
 * @param {Number[][]} A An nxn matrix.
 * @param {Number[]} b A vector of length n.
 * 
 * @example
 * // returns something close to [21, 5, -1]
 * LUPSolveVector([[-2,8,-5],[3,-11,7],[9,-34,21]], [3,1,-2])
 */
function LUPSolveVector(A, b){
    // modifies A in-place to contain matrices L,U
    const pi = LUPDecomposition(A)[0];
    return vectorSolverHelper(A, pi, b);
}

/**
 * Solves Ax=b for the unknown matrix x.
 * @param {Number[][]} A An nxn matrix.
 * @param {Number[]} b An nxm matrix.
 */
function LUPSolveMatrix(A, b){
    // modifies A in-place to contain matrices L,U
    const pi = LUPDecomposition(A)[0];
    const n = A.length;
    const m = b[0].length;
    let inv = Array(n).fill().map(() => Array(m).fill());
    for(let j=0; j<m; j++){
        const b_col = b.map(r => r[j]);
        const col = vectorSolverHelper(A, pi, b_col);
        for(let i=0; i<n; i++){
            inv[i][j] = col[i];
        }
    }
    return inv;
}

/**
 * Computes the matrix inverse of `A`.
 * @param {Number[][]} A An nxn matrix.
 * 
 * @example
 * // the result should produce the identity when multiplied with the input
 * matrixInverse([[7,2,1],[0,3,-1],[-3,4,2]])
 */
function matrixInverse(A){
    const n = A.length;
    let b = Array(n).fill().map(() => Array(n).fill(0));
    for(let i=0; i<n; i++){ b[i][i] = 1; }
    return LUPSolveMatrix(A, b);
}

/**
 * Computes the determinant of `A`.
 * @param {Number[][]} A An nxn matrix.
 * 
 * @example
 * // returns something close to -85
 * determinant([[7,-2,1],[0,-3,-1],[-3,-4,2]])
 */
function determinant(A){
    try{
        let det = LUPDecomposition(A)[1];
        for(let i=0; i<A.length; i++){
            det *= A[i][i];
        }
        return det;
    } catch(e){ return 0; }
}


