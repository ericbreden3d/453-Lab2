class Matrix {
    int** matrix;
    int* arr;  // also maintain 1-D array for sending
    int size;
    Matrix get_detrm_subm(int split);
    int detrm_helper(Matrix& m);
 public:
    Matrix();
    Matrix(int n);
    Matrix(int* buff, int size);
    Matrix(const Matrix& other);
    ~Matrix();
    void fill_rand();
    Matrix operator*(Matrix& other);
    Matrix& operator=(const Matrix& other);
    Matrix get_subm(int len, int x, int y);
    int* get_1d();
    int determinant();
    int& operator()(int x, int y);
    void print();
};