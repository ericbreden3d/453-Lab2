class Matrix {
    int** matrix;
    int* arr;  // also maintain 1-D array for sending
    int size;
    Matrix get_subm(int split);
    int detrm_helper(Matrix m);
 public:
    Matrix(int n);
    ~Matrix();
    void fill_rand();
    Matrix operator*(Matrix& other);
    int determinant();
    int& operator()(int x, int y);
    void print();
};