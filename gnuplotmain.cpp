#include <iostream>
#include <utility>
#include <vector>
#include <iomanip>
#include <math.h>

using namespace std;


/** class for storing matrix **/
class Matrix {
public:
    int n;
    int m;

    vector<vector<double>> matx;
    vector<double> rowElements;

    Matrix(int n, int m) {
        this->n = n;
        this->m = m;
        createMatrix();

    }

    Matrix(int n, int m, vector<vector<double>> a) {
        this->n = n;
        this->m = m;
        this->matx = a;
    }

    void createMatrix() {
        vector<double> buf;

        for (int i = 0; i < this->m; i++)
            buf.push_back(0);

        for (int i = 0; i < this->n; i++)
            this->matx.push_back(buf);
    }

    void addToMatrix(double value) {
        if (rowElements.size() != m) {
            rowElements.push_back(value);
            if (rowElements.size() == m) {
                auto *buf = new vector<double>;
                for (auto el: rowElements) {
                    buf->push_back(el);
                }
                matx.push_back(*buf);
                rowElements.clear();
            }
        }
    }

    Matrix transpose() {
        Matrix gMatrix(this->m, this->n);
        for (int i = 0; i < gMatrix.n; i++) {
            for (int j = 0; j < gMatrix.m; j++) {
                gMatrix.addToMatrix(0);
            }
        }

        for (int i = 0; i < gMatrix.n; i++) {
            for (int j = 0; j < gMatrix.m; j++) {
                gMatrix.setValue(this->returnValue(j, i), i, j);
            }
        }
        return gMatrix;
    }


    virtual Matrix &operator=(Matrix diffMatrix) {
        if (this->n != diffMatrix.n || this->m != diffMatrix.m) {
            throw invalid_argument("Error: the dimensional problem occurred");
        }
        for (int i = 0; i < this->n; i++) {
            for (int j = 0; j < this->m; j++) {
                this->matx[i][j] = diffMatrix.returnValue(i, j);
            }
        }
        return *this;
    }

    virtual Matrix &operator*(Matrix aMatrix) {
        if (this->m != aMatrix.n) {
            throw invalid_argument("Error: the dimensional problem occurred");
        }

        auto *fMatrix = new Matrix(this->n, aMatrix.m);
        for (int i = 0; i < fMatrix->n; i++) {
            for (int j = 0; j < fMatrix->m; j++) {
                fMatrix->addToMatrix(0);

            }
        }
        for (int i = 0; i < fMatrix->n; i++) {
            for (int j = 0; j < fMatrix->m; j++) {
                fMatrix->addToMatrix(0);
                for (int k = 0; k < this->m; k++) {
                    double s = fMatrix->returnValue(i, j);
                    fMatrix->setValue(s + this->returnValue(i, k) * aMatrix.returnValue(k, j),
                                      i, j);
                }
            }
        }

        return *fMatrix;
    }

    Matrix &operator-(Matrix aMatrix) {
        if (aMatrix.m != this->m || aMatrix.n != this->n) {
            throw invalid_argument("Error: the dimensional problem occurred");
        }

        auto *eMatrix = new Matrix(aMatrix.n, aMatrix.m);
        for (int i = 0; i < eMatrix->n; i++) {
            for (int j = 0; j < eMatrix->m; j++) {
                eMatrix->addToMatrix(this->returnValue(i, j) - aMatrix.returnValue(i, j));
            }
        }

        return *eMatrix;
    }

    Matrix &operator+(Matrix bMatrix) {
        if (this->m != bMatrix.m || this->n != bMatrix.n) {
            throw invalid_argument("Error: the dimensional problem occurred");
        }

        auto *dMatrix = new Matrix(this->n, this->m);
        for (int i = 0; i < dMatrix->n; i++) {
            for (int j = 0; j < dMatrix->m; j++) {
                dMatrix->addToMatrix(this->returnValue(i, j) + bMatrix.returnValue(i, j));
            }
        }
        return *dMatrix;
    }

    double returnValue(int row, int col) {
        return matx[row][col];
    }

    void setValue(double value, int row, int col) {
        matx[row][col] = value;
    }

};


istream &operator>>(istream &in, Matrix &matrix) {
    double buf = 0;
    for (int i = 0; i < matrix.n; i++) {
        for (int j = 0; j < matrix.m; j++) {
            in >> buf;
            matrix.addToMatrix(buf);
        }
    }

    return in;
}

ostream &operator<<(ostream &out, Matrix &matrix) {
    for (int i = 0; i < matrix.n; i++) {
        for (int j = 0; j < matrix.m; j++) {
            double value = matrix.returnValue(i, j);
            if (abs(value) < 0.0000000001) {
                cout << 0.00;
            } else
                cout << value;

            if (j + 1 == matrix.m) {
                cout << endl;
            } else {
                cout << " ";
            }
        }
    }
    return out;
}



/** Class for storing square matrix NxN. Include all functions of Matrix class **/

class SquareMatrix : public Matrix {
public:
    explicit SquareMatrix(int n) : Matrix(n, n) {}

    SquareMatrix &operator=(SquareMatrix diffMatrix) {
        for (int i = 0; i < this->n; i++) {
            for (int j = 0; j < this->m; j++) {
                this->matx[i][j] = diffMatrix.returnValue(i, j);
            }
        }
        return *this;
    }

    virtual SquareMatrix &operator*(SquareMatrix aMatrix) {

        auto *fMatrix = new SquareMatrix(this->n);
        for (int i = 0; i < fMatrix->n; i++) {
            for (int j = 0; j < fMatrix->m; j++) {
                fMatrix->addToMatrix(0);
            }
        }
        for (int i = 0; i < fMatrix->n; i++) {
            for (int j = 0; j < fMatrix->m; j++) {
                fMatrix->addToMatrix(0);
                for (int k = 0; k < this->m; k++) {
                    double s = fMatrix->returnValue(i, j);
                    fMatrix->setValue(s + this->returnValue(i, k) * aMatrix.returnValue(k, j),
                                      i, j);
                }
            }
        }

        return *fMatrix;
    }


    SquareMatrix &operator-(SquareMatrix aMatrix) {
        auto *eMatrix = new SquareMatrix(aMatrix.n);
        for (int i = 0; i < aMatrix.n; i++) {
            for (int j = 0; j < eMatrix->m; j++) {
                eMatrix->addToMatrix(this->returnValue(i, j) - aMatrix.returnValue(i, j));
            }
        }

        return *eMatrix;
    }

    SquareMatrix &operator+(SquareMatrix bMatrix) {
        auto *dMatrix = new SquareMatrix(bMatrix.n);
        for (int i = 0; i < bMatrix.n; i++) {
            for (int j = 0; j < dMatrix->m; j++) {
                dMatrix->addToMatrix(this->returnValue(i, j) + bMatrix.returnValue(i, j));
            }
        }
        return *dMatrix;
    }

};

/** Class that create square identity matrix. Include the same functions as SquareMatrix class **/
class IdentityMatrix : public SquareMatrix {
public:
    explicit IdentityMatrix(int n) : SquareMatrix(n) {
        for (int i = 0; i < n; i++) {
            this->matx[i][i] = 1;
        }
    }

    Matrix &operator*(Matrix aMatrix) {

        Matrix *fMatrix = new SquareMatrix(this->n);
        for (int i = 0; i < fMatrix->n; i++) {
            for (int j = 0; j < fMatrix->m; j++) {
                fMatrix->addToMatrix(0);
            }
        }
        for (int i = 0; i < fMatrix->n; i++) {
            for (int j = 0; j < fMatrix->m; j++) {
                fMatrix->addToMatrix(0);
                for (int k = 0; k < this->m; k++) {
                    double s = fMatrix->returnValue(i, j);
                    fMatrix->setValue(s + this->returnValue(i, k) * aMatrix.returnValue(k, j),
                                      i, j);
                }
            }
        }

        return *fMatrix;
    }
};


/** Class that makes elimination operation of square matrix **/
class EliminationMatrix : public Matrix {
public:
    EliminationMatrix(int n, int e1, int e2, double num) : Matrix(n, n) {
        for (int i = 0; i < n; i++) {
            this->matx[i][i] = 1;
        }
        this->matx[e1][e2] -= num;
    }
};

/** Class that makes permutation operation of square matrix **/

class PermutationMatrix : public IdentityMatrix {
public:
    explicit PermutationMatrix(int n, int i1, int i2) : IdentityMatrix(n) {
        swap(this->matx[i1], this->matx[i2]);
    }
};

/** Class for storing column vector. Include all functions of Matrix class**/
class ColumnVector : public Matrix {
private:
    int size;
    int cols = 1;
    vector<double> column;
public:
    explicit ColumnVector(int n) : Matrix(n, 1) {}
    explicit ColumnVector(int n, vector<double> a) : Matrix(n, 1) {
        for (int i = 0; i < n; i++) {
            this->matx[i][0] = a[i];
        }
    }

    ColumnVector &operator=(Matrix diffMatrix) override {

        for (int i = 0; i < this->n; i++) {
            for (int j = 0; j < this->m; j++) {
                this->matx[i][j] = diffMatrix.returnValue(i, j);
            }
        }
        return *this;
    }

    ColumnVector &operator*(ColumnVector aMatrix) {
        if (this->m != aMatrix.n) {
            throw invalid_argument("Error: the dimensional problem occurred");
        }

        auto *fMatrix = new ColumnVector(this->n);
        for (int i = 0; i < fMatrix->n; i++) {
            for (int j = 0; j < fMatrix->m; j++) {
                fMatrix->addToMatrix(0);

            }
        }
        for (int i = 0; i < fMatrix->n; i++) {
            for (int j = 0; j < fMatrix->m; j++) {
                fMatrix->addToMatrix(0);
                for (int k = 0; k < this->m; k++) {
                    int s = fMatrix->returnValue(i, j);
                    fMatrix->setValue(s + this->returnValue(i, k) * aMatrix.returnValue(k, j),
                                      i, j);
                }
            }
        }

        return *fMatrix;
    }

    ColumnVector &operator-(ColumnVector aMatrix) {
        if (aMatrix.m != this->m || aMatrix.n != aMatrix.n) {
            throw invalid_argument("Error: the dimensional problem occurred");
        }

        auto *eMatrix = new ColumnVector(aMatrix.n);
        for (int i = 0; i < aMatrix.n; i++) {
            for (int j = 0; j < eMatrix->m; j++) {
                eMatrix->addToMatrix(this->returnValue(i, j) - aMatrix.returnValue(i, j));
            }
        }
        return *eMatrix;
    }

    ColumnVector &operator+(ColumnVector bMatrix) {
        auto *dMatrix = new ColumnVector(bMatrix.n);
        for (int i = 0; i < bMatrix.n; i++) {
            for (int j = 0; j < dMatrix->m; j++) {
                dMatrix->addToMatrix(this->returnValue(i, j) + bMatrix.returnValue(i, j));
            }
        }
        return *dMatrix;
    }

    /** count the norm **/
    double norm() const {
        double result = 0;
        for (int i = 0; i < this->column.size(); i++) {
            result += this->matx[i][0] * this->matx[i][0];
        }
        return sqrt(result);
    }

};


void twoMatricesPrint(Matrix matrix1, Matrix matrix2) {
    for (int i = 0; i < matrix1.n; i++) {
        for (int j = 0; j < matrix1.m; j++) {
            double value = matrix1.returnValue(i, j);
            if (abs(value) < 0.0000000001) {
                cout << 0.0000;
            } else
                cout << value;

            if (j + 1 == matrix1.m) {
                cout << " ";
                for (int q = 0; q < matrix2.m; q++) {
                    value = matrix2.returnValue(i, q);
                    if (abs(value) < 0.0000000001) {
                        cout << 0.0000;
                    } else {
                        cout << value;
                    }
                    if (q + 1 == matrix2.m) cout << endl;
                    else cout << " ";
                }

            } else {
                cout << " ";
            }
        }
    }

}

void diagonalization(Matrix *A, Matrix *b) {
    for (int i = 0; i < A->n; i++) {
        for (int j = 0; j < A->n; j++)
            b->matx[i][j] /= A->matx[i][i];
        A->matx[i][i] = 1;
    }
}

/** The function that find inverse matrix **/
Matrix& InverseMatrix(Matrix &matrix) {
    Matrix *A = new Matrix(matrix);
    Matrix *I = new IdentityMatrix(A->n);

    for (int i = 0; i < A->n; i++) {
        int pivot_row = i;
        for (int j = i + 1; j < A->n; j++) {
            double abs_value = abs(A->matx[j][i]);
            if (abs_value > abs(A->matx[pivot_row][i]))
                pivot_row = j;
        }

        // permutation part
        if (pivot_row != i) {
            Matrix* permutationMatrix = new PermutationMatrix(A->n, i, pivot_row);
            *A = *permutationMatrix * *A;
            *I = *permutationMatrix * *I;
        }

        for (int k = i + 1; k < A->n; k++) {
            if (A->matx[k][i] == 0) {
                continue;
            }

            double factor = A->matx[k][i] / A->matx[i][i];
            Matrix *eliminationMatrix = new EliminationMatrix(A->n, k, i, factor);

            *A = *eliminationMatrix * *A;
            *I = *eliminationMatrix * *I;
        }
    }
    // back elimination
    for (int i = A->n - 1; i >= 0; i--) {
        for (int j = i - 1; j >= 0; j--) {
            if (A->matx[j][i] == 0) {
                continue;
            }

            double factor = A->matx[j][i] / A->matx[i][i];
            Matrix* eliminationMatrix = new EliminationMatrix(A->n, j, i, factor);

            *A = *eliminationMatrix * *A;
            *I = *eliminationMatrix * *I;
        }
    }
    // diagonalization part
    diagonalization(A, I);
    return *I;
}


ColumnVector leastSquares(Matrix A, Matrix b, int degree) {
    cout << "A:" << endl;
    cout << A;

    cout << "A_T*A:" << endl;
    Matrix A_T(A.m, A.n);
    A_T = A.transpose();
    Matrix A_TT(A_T.n, A.m);

    A_TT = A_T * A;
    cout << A_TT;

    cout << "(A_T*A)^-1:" << endl;
    A_TT = InverseMatrix(A_TT);
    cout << A_TT;

    cout << "A_T*b:" << endl;
    ColumnVector A_Tb(A_T.n);
    A_Tb = A_T * b;
    cout << A_Tb;

    cout << "x~:" << endl;
    ColumnVector x(degree + 1);
    x = A_TT * A_Tb;
    cout << x;

    return x;
}

#ifdef WIN64
#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"
#else
#define GNUPLOT_NAME "gnuplot -persist"
#endif

int main() {

    cout << fixed << setprecision(4);

    int m, degree;
    cin >> m;

    vector<double> T;
    vector<double> b;

    for (int i = 0; i < m; i++) {
        double t, q;
        cin >> t >> q;
        T.push_back(t);
        b.push_back(q);
    }
    cin >> degree;

    vector<vector<double>> mtrx;
    for (int i = 0; i < m; i++) {
        vector<double> buf;
        for (int j = 0; j < degree + 1; j++) {
            buf.push_back(pow(T[i], j));
        }
        mtrx.push_back(buf);
    }

    Matrix A(m, degree + 1, mtrx);
    ColumnVector B(m, b);

    ColumnVector x = leastSquares(A, B, degree);
    cout << x;

#ifdef WIN64
    FILE* plotter = _popen(GNUPLOT_NAME, "w");
#else
    FILE* plotter = popen(GNUPLOT_NAME, "w");
#endif

    fprintf(plotter, "%s\n", "set border linewidth 1.5\n"
                             "set style line 1 linecolor rgb '#0060ad' linetype 1 linewidth 2\n"
                             "set style line 2 linecolor rgb '#dd181f' linetype 1 linewidth 2");

    fprintf(plotter, "%s\n", "set key at -50,70\n"
                             "set xlabel 'x'\n"
                             "set ylabel 'y'\n"
                             "set xrange [-100:100]\n"
                             "set yrange [-75:75]\n"
                             "set xtics 25\n"
                             "set ytics 25\n"
                             "set tics scale 0.75");

    fprintf(plotter, "%s", "f(x) = ");
    fprintf(plotter, "%f", x.matx[0][0]);
    for (int i = 1; i < degree + 1; i++) {
        fprintf(plotter, "%s", " + ");
        fprintf(plotter, "%f", x.matx[i][0]);
        fprintf(plotter, "%s", " * x**");
        fprintf(plotter, "%d", i);
    }
    fprintf(plotter, "%s\n", "");

    fprintf(plotter, "%s\n", "plot '-' title 'rawData' with lines linestyle 1, [x=-100:100] f(x) title 'fittedData' with "
                             "lines linestyle 2");
    for (int i = 0; i < m; i++) {
        fprintf(plotter, "%f%f\n", T[i], b[i]);
    }
    fprintf(plotter, "%c\n", 'e');


#ifdef WIN64
    _pclose(plotter);
#else
    pclose(plotter);
#endif

    return 0;
}
