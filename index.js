function swap(arr, idx1, idx2){
    const temp = arr[idx1];
    arr[idx1] = arr[idx2];
    arr[idx2] = temp
}

function LUPDecomposition(A){
    const n = A.length;
    const pi = Array(n).fill().map((_,i) => i);
    for(let i=0; i<n; i++){
        let p = 0;
        let i_pivot = i;
        for(let j=i; j<n; j++){
            if(Math.abs(A[j][i]) > p){
                p = Math.abs(A[j][i]);
                i_pivot = j;
            }
        }
        if(p === 0){ throw("Error: singular matrix"); }
        swap(pi, i, i_pivot);
        swap(A, i, i_pivot);
        for(let j=i+1; j<n; j++){
            A[j][i] = A[j][i] / A[i][i];
            for(let k=i+1; k<n; k++){
                A[j][k] = A[j][k] - (A[j][i] * A[i][k]);
            }
        }
    }
    return [A,pi];
}

/*
LUPDecomposition([[2,4,1,-3],[-1,-2,2,4],[4,2,-3,5],[5,-4,-3,1]])
*/
