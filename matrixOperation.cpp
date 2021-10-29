#include "matrixOperation.h"
// A-B, A is a matrix, B is a vector
std::vector<double> diff_vector(std::vector< double >  &in_A, std::vector<double> &in_B) {
    std::vector<double> out_c; out_c.resize(in_A.size());
    for (int i=0; i<in_A.size(); i++) {
        out_c[i]=in_A[i] - in_B[i];
    }
    return out_c;
}

// A+B, A is a vector, B is a vector
std::vector<double> add_vector(std::vector< double >  &in_A, std::vector<double> &in_B) {
    std::vector<double> out_c; out_c.resize(in_A.size());
    for (int i=0; i<in_A.size(); i++) {
        out_c[i]=in_A[i] + in_B[i];
    }
    return out_c;
}

// A is a matrix
std::vector< std::vector<double> > transpose_matrix(std::vector< std::vector<double> > &in_A ) {
    std::vector< std::vector<double> > out_c; out_c.resize(in_A[0].size());
    for (int i=0; i<in_A[0].size(); i++) {
        out_c[i].resize(in_A.size());
        for (int k=0; k<in_A.size(); k++) {
            out_c[i][k] = in_A[k][i];
        }
    }
    return out_c;
}

// A is a matrix, B is a matrix
std::vector< std::vector<double> > add_matrix(std::vector< std::vector<double> >  &in_A, std::vector< std::vector<double> > &in_B) {
    std::vector< std::vector<double> > out_c; out_c.resize(in_A.size());
    for (int i=0; i<in_A.size(); i++) {
        out_c[i].resize(in_A[i].size());
        for (int k=0; k<in_A[i].size(); k++) {
            out_c[i][k] = in_A[i][k] + in_B[i][k];
        }
    }
    return out_c;
}

// A is a matrix, B is a vector
std::vector<double> left_Mult_AB(std::vector< std::vector<double> >  &in_A, std::vector<double> &in_B) {

    std::vector<double> out_c; out_c.resize(in_A.size());
    for (int i=0; i<in_A.size(); i++) {
        double sum=0.0;
        for (int k=0; k<in_B.size(); k++) {
            sum +=  in_A[i][k] * in_B[k];
        }
        out_c[i]=sum;
    }

    return out_c;
}

// A is a vector, B is a matrix
std::vector<double> right_Mult_AB ( std::vector<double> &in_A, std::vector< std::vector<double> > &in_B) {
    if(in_A.size() != in_B.size())
        std::cout << "Dimension mismatch.."<< std::endl;
    std::vector<double> out_c; out_c.resize(in_B[0].size());
    for (int i=0; i<in_B[0].size(); i++) {
        double sum=0.0;
        for (int k=0; k<in_A.size(); k++) {
            sum +=  in_A[k] * in_B[k][i];
        }
        out_c[i]=sum;
    }
    return out_c;
}

std::vector< std::vector<double> > mult_ABC ( std::vector< std::vector<double> > &in_A, std::vector< std::vector<double> > &in_B, std::vector< std::vector<double> > &in_C) {

    std::vector< std::vector<double> > out_c, temp_c;
    temp_c = mult_AB ( in_B, in_C);
    out_c = mult_AB ( in_A, temp_c);

    return out_c;
}

std::vector< std::vector<double> > mult_AB ( std::vector< std::vector<double> > &in_A, std::vector< std::vector<double> > &in_B) {

    if(in_A[0].size() != in_B.size())
        std::cout << "Matrix dimension mismatch in mult_AB: " << in_A[0].size()  << " " << in_B.size() << std::endl;
    std::vector< std::vector<double> > out_c; out_c.resize(in_A.size());
    int D1 = in_A.size();
    int D2 = in_B.size();
    int D3 = in_B[0].size();
    for (int i=0; i<D1; i++) {
        out_c[i].resize(D3);
        for (int j=0; j<D3; j++) {
            double sum=0.0;
            for (int k=0; k<D2; k++) {
                sum +=  in_A[i][k] * in_B[k][j];
            }
            out_c[i][j]=sum;
        }
    }
    return out_c;
}



double dot_product ( std::vector<double>  &in_A, std::vector<double> &in_B) {

    double out_c = 0;
    for (int i=0; i<in_A.size(); i++) {
        out_c += in_A[i]*in_B[i];
    }
    return out_c;
}

double norm(std::vector<double> &in_A){
    double res = 0;
    for (int i = 0; i < in_A.size(); ++i) {
        res+=in_A[i]*in_A[i];
    }
    res = sqrt(res);
    return res;
}

double norm(std::vector<double> &in_A, std::vector<double> &in_B){
    double res = 0;
    for (int i = 0; i < in_A.size(); ++i) {
        res+=(in_A[i] - in_B[i])*(in_A[i] - in_B[i]);
    }

    res = sqrt(res);
    return res;
}

//Frobenius norm
double norm(std::vector<std::vector<double>> &in_A){
    int row = in_A.size();
    int col = in_A[0].size();
    double res = 0;
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            res += in_A[i][j] * in_A[i][j];
        }
    }

    res = sqrt(res);
    return res;
}
// Function to get cofactor of A[p][q] in temp[][]. n is current
// dimension of A[][]
void getCofactor( std::vector< std::vector<double> > &in_A, std::vector< std::vector<double> > &out_, int p, int q, int n)
{
    int i = 0, j = 0;

    // Looping for each element of the matrix
    for (int row = 0; row < n; row++)
    {
        for (int col = 0; col < n; col++)
        {
            // Copying into temporary matrix only those element
            // which are not in given row and column
            if (row != p && col != q)
            {
                out_[i][j++] = in_A[row][col];

                // Row is filled, so increase row index and
                // reset col index
                if (j == n - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

/* Recursive function for finding determinant of matrix.
n is current dimension of A[][]. */
double determinant(std::vector< std::vector<double> > &in_A, int n)
{
    double D = 0; // Initialize result

    int N = in_A.size();
    // Base case : if matrix contains single element
    if (n == 1)
        return in_A[0][0];

    std::vector< std::vector<double> > temp; // To store cofactors
    temp.resize(N);
    for (int i = 0; i < N; ++i) {
        temp[i].resize(N);
    }

    int sign = 1; // To store sign multiplier

    // Iterate for each element of first row
    for (int f = 0; f < n; f++)
    {
        // Getting Cofactor of A[0][f]
        getCofactor(in_A, temp, 0, f, n);
        D += sign * in_A[0][f] * determinant(temp, n - 1);

        // terms are to be added with alternate sign
        sign = -sign;
    }

    return D;
}

// Function to get adjoint of A[N][N] in adj[N][N].
void adjoint(std::vector< std::vector<double> > &in_A, std::vector< std::vector<double> > &adj)
{
    int N = in_A.size();
    if (N == 1)
    {
        adj[0][0] = 1;
        return;
    }

    // temp is used to store cofactors of A[][]
    int sign = 1;
    std::vector< std::vector<double> > temp; // To store cofactors
    temp.resize(N);
    for (int i = 0; i < N; ++i) {
        temp[i].resize(N);
    }


    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            // Get cofactor of A[i][j]
            getCofactor(in_A, temp, i, j, N);

            // sign of adj[j][i] positive if sum of row
            // and column indexes is even.
            sign = ((i+j)%2==0)? 1: -1;

            // Interchanging rows and columns to get the
            // transpose of the cofactor matrix
            adj[j][i] = (sign)*(determinant(temp, N-1));
        }
    }
}

// Function to calculate and store inverse, returns false if
// matrix is singular
bool inverse(std::vector< std::vector<double> > &in_A, std::vector< std::vector<double> > &out_inverse)
{

    // Find determinant of A[][]
    int N = in_A.size();
    out_inverse.resize(N);
    double det = determinant(in_A, N);

    if (fabs(det) < 1e-10)
    {
        std::cout << "Singular matrix, can't find its inverse" << std::endl;
        std::cout << "det: " << det << std::endl;
        return false;
    }

    // Find adjoint
    std::vector< std::vector<double> > adj; // To store cofactors
    adj.resize(N);
    for (int i = 0; i < N; ++i) {
        adj[i].resize(N);
    }
    adjoint(in_A, adj);

    // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
    for (int k=0; k<N; k++){
        out_inverse[k].resize(N);
        for (int j=0; j<N; j++)
            out_inverse[k][j] = adj[k][j]/det;
    }

    return true;
}
