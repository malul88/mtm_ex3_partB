
#include <functional>
#include <string>
#include <iostream>
#include <sstream>   
#include <ostream>
#include <map>
#include <iterator>
#include <string>
#include <fstream>

#include "Matrix.h"

//#define NDEPENDENCY

using namespace mtm;
using std::cout;
using std::endl;
using std::string;
using std::map;
using std::pair;
using std::ostream;

#define ASSERT_TEST(x) if(!(x)){ \
 cout<<"Failed assertion at line "<<__LINE__<<" in "<<__func__<<endl;\
  return false; }

#define ADD_TEST(x) tests[#x]=x;

class DependencyFinder{

private:

    int data;

    DependencyFinder(int data) : data(data){}

public:

    static map<string, bool> used;

    static void reset(){
        used["ctor"] = false;
        used["cctor"] = false;
        used["="] = false;
        used["dtor"] = false;
        used["+"] = false;
        used["+="] = false;
        used["uni -"] = false;
        used["bin -"] = false;
        used["<"] = false;
        used["<="] = false;
        used[">"] = false;
        used[">="] = false;
        used["=="] = false;
        used["!="] = false;
        used["bool"] = false;
        used["<<"] = false;
    }

    DependencyFinder(){data = 0; used["ctor"] = true;}
    DependencyFinder(const DependencyFinder& other){used["cctor"] = true; data = other.data;}
    DependencyFinder& operator= (const DependencyFinder& other){used["="] = true; data = other.data; return *this;}
    ~DependencyFinder(){used["dtor"] = true;}
    DependencyFinder operator+(const DependencyFinder& other) const {used["+"] = true; return DependencyFinder(data + other.data);}
    DependencyFinder& operator+=(const DependencyFinder& other){used["+="] = true; data += other.data; return *this;}
    DependencyFinder operator-() const{used["uni -"] = true; return DependencyFinder(-data);}
    DependencyFinder operator-(const DependencyFinder& other) const{used["bin -"] = true; return DependencyFinder(data - other.data);}
    bool operator<(const DependencyFinder& other) const{used["<"] = true; return data < other.data;}
    bool operator<=(const DependencyFinder& other) const{used["<="] = true; return data <= other.data;}
    bool operator>(const DependencyFinder& other) const{used[">"] = true; return data > other.data;}
    bool operator>=(const DependencyFinder& other) const{used[">="] = true; return data >= other.data;}
    bool operator==(const DependencyFinder& other) const{used["=="] = true; return data == other.data;}
    bool operator!=(const DependencyFinder& other) const{used["!="] = true; return data != other.data;}
    operator bool() const{used["bool"] = true; return data;}
    friend ostream& operator<<(ostream& os, DependencyFinder const d1){used["<<"] = true; return os << d1.data;}

};
map<string, bool> DependencyFinder::used;

DependencyFinder fnc (DependencyFinder val)
{
    return val;
}

void printCurrentDependency(ostream& os, string fname){

    os << fname << ",";
    for (pair<string, bool> element : DependencyFinder::used)
    {
        os << (element.second? "1":" ");
        os << ",";
    }
    os << endl;

}

bool findDependencies(ostream& os){

    os << "Searching for dependencies:"<<endl;

    DependencyFinder::reset();
    os << "Function" << ",";
    for (pair<string, bool> element : DependencyFinder::used)
    {
        os << element.first << ",";
    }
    os << endl;

    DependencyFinder::reset();
    Matrix<DependencyFinder> matrix1(Dimensions(3,3));
    Matrix<DependencyFinder> matrix2(Dimensions(3,3), DependencyFinder());
    printCurrentDependency(os, "Constructor");

    DependencyFinder::reset();
    Matrix<DependencyFinder> matrix3(matrix2);
    printCurrentDependency(os, "Copy Constructor");

    DependencyFinder::reset();
    matrix3 = matrix1;
    printCurrentDependency(os, "Assignment Operator");

    DependencyFinder::reset();
    Matrix<DependencyFinder>::Diagonal(10, DependencyFinder());
    printCurrentDependency(os, "Diagonal");

    DependencyFinder::reset();
    matrix1.transpose();
    printCurrentDependency(os, "Transpose");

    DependencyFinder::reset();
    -matrix1;
    printCurrentDependency(os, "Unary -");

    DependencyFinder::reset();
    matrix2-matrix1;
    printCurrentDependency(os, "Binary -");

    DependencyFinder::reset();
    matrix2+matrix1;
    printCurrentDependency(os, "Binary +");

    DependencyFinder::reset();
    std::stringstream buffer;
    buffer << matrix1;
    printCurrentDependency(os, "Operator<<");

    DependencyFinder::reset();
    matrix1(0,0);
    printCurrentDependency(os, "Operator()");

    DependencyFinder::reset();
    matrix1<DependencyFinder();
    printCurrentDependency(os, "Operator<");

    DependencyFinder::reset();
    matrix1<=DependencyFinder();
    printCurrentDependency(os, "Operator<=");

    DependencyFinder::reset();
    matrix1>DependencyFinder();
    printCurrentDependency(os, "Operator>");

    DependencyFinder::reset();
    matrix1>=DependencyFinder();
    printCurrentDependency(os, "Operator>=");

    DependencyFinder::reset();
    matrix1==DependencyFinder();
    printCurrentDependency(os, "Operator==");

    DependencyFinder::reset();
    matrix1!=DependencyFinder();
    printCurrentDependency(os, "Operator!=");

    DependencyFinder::reset();
    any(matrix1);
    printCurrentDependency(os, "Any");

    DependencyFinder::reset();
    all(matrix1);
    printCurrentDependency(os, "All");

    DependencyFinder::reset();
    matrix1.apply(fnc);
    printCurrentDependency(os, "Apply");

    return true;

}

template <class T>
bool checkAreEqual(const Matrix<T>& first, const Matrix<T>& second){

    if (first.height() != second.height() || first.width() != second.width()){
        return false;
    }

    for (int row = 0; row < first.height(); row++){
        for (int col = 0; col < first.width(); col++){
            if (first(row,col) != second(row,col)){
                return false;
            }
        }
    }

    return true;

}

const int N = 1000;
int sampleData[N];
void createSampleData(){

    // Sample Fibonacci sequence
    sampleData[0] = 1;
    sampleData[1] = 1;

    for (int i = 2; i < N; i++){
        sampleData[i] = (sampleData[i-1] + sampleData[i-2]) % 1000;
    }

}

class T1 : public string{
    public:
        T1(string name = "def"):string(name){}
};

bool testConstructor(){
    
    int rows = 57;
    int cols = 27;
    int scalar = 7;
    Matrix<int> matrix(Dimensions(rows, cols), scalar);

    ASSERT_TEST(matrix.height() == rows);
    ASSERT_TEST(matrix.width() == cols);
    ASSERT_TEST(matrix.size() == (rows * cols));

    for (int& element : matrix){
        ASSERT_TEST(element == scalar);
    }

    Matrix<int> matrix2(Dimensions(rows, cols));
    for (int& element : matrix2){
        ASSERT_TEST(element == 0);
    }

    ASSERT_TEST(!checkAreEqual(matrix, matrix2));

    Matrix<T1> matrix3(Dimensions(rows, cols));

    ASSERT_TEST(matrix3.height() == rows);
    ASSERT_TEST(matrix3.width() == cols);
    ASSERT_TEST(matrix3.size() == (rows * cols));

    for (T1& element : matrix3){
        ASSERT_TEST(element == "def");
    }

    Matrix<T1> matrix4(Dimensions(rows, cols), T1("other"));

    ASSERT_TEST(matrix4.height() == rows);
    ASSERT_TEST(matrix4.width() == cols);
    ASSERT_TEST(matrix4.size() == (rows * cols));

    for (T1& element : matrix4){
        ASSERT_TEST(element == "other");
    }

    try{
        Matrix<T1> mat(Dimensions(-3,-3));
        ASSERT_TEST(false);
    }
    catch(const Matrix<T1>::IllegalInitialization& e){
        // Tests that IllegalInitialization extends std::exception
        std::exception e_test = e;
        ASSERT_TEST(string(e.what()) == "Mtm matrix error: Illegal initialization values");
    }

    try{
        Matrix<T1> mat(Dimensions(0,0));
        ASSERT_TEST(false);
    }
    catch(const Matrix<T1>::IllegalInitialization& e){
        ASSERT_TEST(string(e.what()) == "Mtm matrix error: Illegal initialization values");
    }

    try{
        Matrix<T1> mat(Dimensions(1,0));
        ASSERT_TEST(false);
    }
    catch(const Matrix<T1>::IllegalInitialization& e){
        ASSERT_TEST(string(e.what()) == "Mtm matrix error: Illegal initialization values");
    }

    try{
        Matrix<T1> mat(Dimensions(0,1));
        ASSERT_TEST(false);
    }
    catch(const Matrix<T1>::IllegalInitialization& e){
        ASSERT_TEST(string(e.what()) == "Mtm matrix error: Illegal initialization values");
    }

    try{
        Matrix<T1> mat(Dimensions(0,1), T1());
        ASSERT_TEST(false);
    }
    catch(const Matrix<T1>::IllegalInitialization& e){
        ASSERT_TEST(string(e.what()) == "Mtm matrix error: Illegal initialization values");
    }

    return true;

}

bool testCopyConstructor(){

    int rows = 57;
    int cols = 27;
    int scalar = -5;
    Matrix<int> matrix(Dimensions(rows, cols), scalar);

    ASSERT_TEST(matrix.height() == rows);
    ASSERT_TEST(matrix.width() == cols);
    ASSERT_TEST(matrix.size() == (rows * cols));

    for (int& element : matrix){
        ASSERT_TEST(element == scalar);
    }

    Matrix<int> matrix2(matrix);
    ASSERT_TEST(matrix2.height() == rows);
    ASSERT_TEST(matrix2.width() == cols);
    ASSERT_TEST(matrix2.size() == (rows * cols));
    for (int& element : matrix2){
        ASSERT_TEST(element == scalar);
    }

    ASSERT_TEST(checkAreEqual(matrix, matrix2));

    Matrix<int> matrix3 = matrix;
    ASSERT_TEST(matrix3.height() == rows);
    ASSERT_TEST(matrix3.width() == cols);
    ASSERT_TEST(matrix3.size() == (rows * cols));
    for (int& element : matrix3){
        ASSERT_TEST(element == scalar);
    }

    ASSERT_TEST(checkAreEqual(matrix2, matrix3));
    ASSERT_TEST(checkAreEqual(matrix, matrix3));

    Matrix<T1> matrix4(Dimensions(rows,cols), T1("deff"));
    ASSERT_TEST(matrix4.height() == rows);
    ASSERT_TEST(matrix4.width() == cols);
    ASSERT_TEST(matrix4.size() == (rows * cols));

    ASSERT_TEST(checkAreEqual(matrix4, Matrix<T1>(Dimensions(rows,cols), T1("deff"))));

    Matrix<T1> matrix5(matrix4);
    ASSERT_TEST(matrix5.height() == rows);
    ASSERT_TEST(matrix5.width() == cols);
    ASSERT_TEST(matrix5.size() == (rows * cols));

    ASSERT_TEST(checkAreEqual(matrix5, Matrix<T1>(Dimensions(rows,cols), T1("deff"))));
    ASSERT_TEST(checkAreEqual(matrix4, matrix5));

    return true;

}

bool testAssignmentOperator(){

    int rows = 57;
    int cols = 27;
    int scalar = -5;
    Matrix<int> matrix(Dimensions(rows, cols), scalar);

    ASSERT_TEST(matrix.height() == rows);
    ASSERT_TEST(matrix.width() == cols);
    ASSERT_TEST(matrix.size() == (rows * cols));
    for (int& element : matrix){
        ASSERT_TEST(element == scalar);
    }

    Matrix<int> matrix2(Dimensions(1,2), -1);
    ASSERT_TEST(matrix2.height() == 1);
    ASSERT_TEST(matrix2.width() == 2);
    ASSERT_TEST(matrix2.size() == (1 * 2));
    for (int& element : matrix2){
        ASSERT_TEST(element == -1);
    }

    ASSERT_TEST(!checkAreEqual(matrix, matrix2));

    matrix2 = matrix;
    ASSERT_TEST(matrix2.height() == rows);
    ASSERT_TEST(matrix2.width() == cols);
    ASSERT_TEST(matrix2.size() == (rows * cols));
    for (int& element : matrix2){
        ASSERT_TEST(element == scalar);
    }

    ASSERT_TEST(checkAreEqual(matrix, matrix2));

    matrix = matrix2;
    ASSERT_TEST(matrix.height() == rows);
    ASSERT_TEST(matrix.width() == cols);
    ASSERT_TEST(matrix.size() == (rows * cols));
    for (int& element : matrix){
        ASSERT_TEST(element == scalar);
    }

    ASSERT_TEST(checkAreEqual(matrix, matrix2));

    Matrix<T1> matrix4(Dimensions(rows,cols), T1("deff"));
    ASSERT_TEST(matrix4.height() == rows);
    ASSERT_TEST(matrix4.width() == cols);
    ASSERT_TEST(matrix4.size() == (rows * cols));

    ASSERT_TEST(checkAreEqual(matrix4, Matrix<T1>(Dimensions(rows,cols), T1("deff"))));

    Matrix<T1> matrix5(Dimensions(rows,cols));
    ASSERT_TEST(!checkAreEqual(matrix5, Matrix<T1>(Dimensions(rows,cols), T1("deff"))));
    ASSERT_TEST(!checkAreEqual(matrix4, matrix5));

    matrix5 = matrix4;

    ASSERT_TEST(checkAreEqual(matrix5, Matrix<T1>(Dimensions(rows,cols), T1("deff"))));
    ASSERT_TEST(checkAreEqual(matrix4, matrix5));

    return true;

}

bool testDiagonal(){

    int n=107;
    Matrix<int> mat = Matrix<int>::Diagonal(n, 1);

    ASSERT_TEST(mat.height() == n);
    ASSERT_TEST(mat.width() == n);
    ASSERT_TEST(mat.size() == n*n);

    for (int row = 0; row < n; row++){
        for (int col = 0; col < n; col++){
            ASSERT_TEST(mat(row,col) == (row==col));
        }
    }

    Matrix<T1> mat2 = Matrix<T1>::Diagonal(n, T1("default"));

    ASSERT_TEST(mat2.height() == n);
    ASSERT_TEST(mat2.width() == n);
    ASSERT_TEST(mat2.size() == n*n);

    for (int row = 0; row < n; row++){
        for (int col = 0; col < n; col++){
            ASSERT_TEST(mat2(row,col) == ((row==col)?T1("default"):T1()));
        }
    }

    try{
        Matrix<T1> mat2 = Matrix<T1>::Diagonal(0, T1("default"));
        ASSERT_TEST(false);
    }
    catch(const Matrix<T1>::IllegalInitialization& e){
        ASSERT_TEST(string(e.what()) == "Mtm matrix error: Illegal initialization values");
    }

    try{
        Matrix<T1> mat2 = Matrix<T1>::Diagonal(-10, T1("default"));
        ASSERT_TEST(false);
    }
    catch(const Matrix<T1>::IllegalInitialization& e){
        ASSERT_TEST(string(e.what()) == "Mtm matrix error: Illegal initialization values");
    }

    return true;

}

bool testTranspose(){

    int n=107;
    Matrix<int> mat = mtm::Matrix<int>::Diagonal(n, 1);
    ASSERT_TEST(checkAreEqual(mat, mat.transpose()));

    int rows = 17;
    int cols = 5;
    Matrix<int> mat2 = Matrix<int>(Dimensions(rows,cols));
    ASSERT_TEST(mat2.height() == rows);
    ASSERT_TEST(mat2.width() == cols);
    ASSERT_TEST(mat2.size() == rows*cols);
    int i = 0;
    for (int& element : mat2){
        element = sampleData[i++];
    }
    
    Matrix<int> mat3 = mat2.transpose();
    ASSERT_TEST(mat3.height() == cols);
    ASSERT_TEST(mat3.width() == rows);
    ASSERT_TEST(mat3.size() == rows*cols);

    for (int row = 0; row < rows; row++){
        for (int col = 0; col < cols; col++){
            ASSERT_TEST(mat2(row,col)==mat3(col,row));
        }
    }

    ASSERT_TEST(checkAreEqual(mat2, (mat2.transpose()).transpose()));

    return true;

}

bool testOperatorNegative(){

    int rows = 17;
    int cols = 5;
    Matrix<int> mat = Matrix<int>(Dimensions(rows,cols));
    ASSERT_TEST(mat.height() == rows);
    ASSERT_TEST(mat.width() == cols);
    ASSERT_TEST(mat.size() == rows*cols);
    int i = 0;
    for (int& element : mat){
        element = sampleData[i++];
    }
    
    Matrix<int> mat2 = -mat;
    ASSERT_TEST(mat2.height() == rows);
    ASSERT_TEST(mat2.width() == cols);
    ASSERT_TEST(mat2.size() == rows*cols);

    for (int row = 0; row < rows; row++){
        for (int col = 0; col < cols; col++){
            ASSERT_TEST(mat2(row,col)==-mat(row,col));
        }
    }

    ASSERT_TEST(checkAreEqual(mat, -(-mat)));
    
    Matrix<int> mat3 = Matrix<int>(Dimensions(rows+1,cols));

    try{
        Matrix<int> mat4 = mat3-mat;
        ASSERT_TEST(false);
    }
    catch(const Matrix<int>::DimensionMismatch& e){
        // Tests that DimensionMismatch extends std::exception
        std::exception e_test = e;
        ASSERT_TEST(string(e.what()) == "Mtm matrix error: Dimension mismatch: (18,5) (17,5)");
    }

    return true;

}

bool testOperatorMatrixAddition(){

    int rows = 17;
    int cols = 5;
    Matrix<int> mat = Matrix<int>(Dimensions(rows,cols));
    ASSERT_TEST(mat.height() == rows);
    ASSERT_TEST(mat.width() == cols);
    ASSERT_TEST(mat.size() == rows*cols);
    int i = 0;
    for (int& element : mat){
        element = sampleData[i++];
    }
    
    Matrix<int> mat2 = -mat;
    ASSERT_TEST(mat2.height() == rows);
    ASSERT_TEST(mat2.width() == cols);
    ASSERT_TEST(mat2.size() == rows*cols);

    ASSERT_TEST(checkAreEqual((mat=mat+mat2), Matrix<int>(Dimensions(rows,cols))));
    ASSERT_TEST(!checkAreEqual(mat+mat2, Matrix<int>(Dimensions(rows,cols))));
    ASSERT_TEST(checkAreEqual((-mat2) + mat2, Matrix<int>(Dimensions(rows,cols))));

    Matrix<int> mat3 = mat2 + mat2;
    mat3 = mat3 + mat2;
    mat3 = mat3 - mat2;
    i = 0;
    for (int& element : mat3){
        ASSERT_TEST(element == -2*sampleData[i++]);
    }

    mat3=((mat3=(mat3 - mat2)) - mat2);

    ASSERT_TEST(checkAreEqual(mat3, Matrix<int>(Dimensions(rows,cols))));

    Matrix<int> mat5 = Matrix<int>(Dimensions(rows+1,cols));

    try{
        Matrix<int> mat4 = mat5+mat;
        ASSERT_TEST(false);
    }
    catch(const Matrix<int>::DimensionMismatch& e){
        ASSERT_TEST(string(e.what()) == "Mtm matrix error: Dimension mismatch: (18,5) (17,5)");
    }

    return true;

}

bool testOperatorScalarAddition(){

    int rows = 17;
    int cols = 5;
    Matrix<int> mat = Matrix<int>(Dimensions(rows,cols));
    ASSERT_TEST(mat.height() == rows);
    ASSERT_TEST(mat.width() == cols);
    ASSERT_TEST(mat.size() == rows*cols);
    int i = 0;
    for (int& element : mat){
        element = sampleData[i++];
    }

    ASSERT_TEST(checkAreEqual(mat, 1 + mat + (-1)));
    ASSERT_TEST(!checkAreEqual(mat, mat + 2 + (-1)));

    (mat += 1) += (-3);
    i = 0;
    for (int& element : mat){
        ASSERT_TEST(element == (sampleData[i++] - 2));
    }

    Matrix<int> mat2 = Matrix<int>(Dimensions(rows,cols), 5);
    mat2 += -5;
    ASSERT_TEST(checkAreEqual(1 + mat2 + (-1), Matrix<int>(Dimensions(rows, cols))));

    Matrix<string> mat5(Matrix<string>(Dimensions(1,2)));
    mat5(0,0) = "Hello world";
    mat5(0,1) = "Bye bye world";

    mat5 += string("!");
    ASSERT_TEST(mat5(0,0) == "Hello world!");
    ASSERT_TEST(mat5(0,1) == "Bye bye world!");

    Matrix<string> mat6 = mat5 + string("*");
    ASSERT_TEST(mat6(0,0) == "Hello world!*");
    ASSERT_TEST(mat6(0,1) == "Bye bye world!*");

    Matrix<string> mat7 = string("*") + mat6;
    ASSERT_TEST(mat7(0,0) == "*Hello world!*");
    ASSERT_TEST(mat7(0,1) == "*Bye bye world!*");

    return true;

}

bool testOperatorOutput(){

    int rows = 2;
    int cols = 5;
    Matrix<int> mat = Matrix<int>(Dimensions(rows,cols));
    ASSERT_TEST(mat.height() == rows);
    ASSERT_TEST(mat.width() == cols);
    ASSERT_TEST(mat.size() == rows*cols);
    int i = 0;
    for (int& element : mat){
        element = sampleData[i++];
    }

    std::string expected = "1 1 2 3 5 \n8 13 21 34 55 \n";
    std::stringstream buffer;
    buffer << mat;
    ASSERT_TEST(expected == buffer.str());
    return true;

}

bool testOperatorParenthesis(){

    Matrix<int> matrix(Dimensions(1,1));

    try{
        matrix(1,0);
        ASSERT_TEST(false);
    }
    catch(const Matrix<int>::AccessIllegalElement& e){
        // Tests that AccessIllegalElement extends std::exception
        std::exception e_test = e;
        ASSERT_TEST(string(e.what()) == "Mtm matrix error: An attempt to access an illegal element");
    }

    try{
        matrix(1,1);
        ASSERT_TEST(false);
    }
    catch(const Matrix<int>::AccessIllegalElement& e){
        ASSERT_TEST(string(e.what()) == "Mtm matrix error: An attempt to access an illegal element");
    }

    try{
        matrix(0,1);
        ASSERT_TEST(false);
    }
    catch(const Matrix<int>::AccessIllegalElement& e){
        ASSERT_TEST(string(e.what()) == "Mtm matrix error: An attempt to access an illegal element");
    }

    try{
        matrix(-1,0);
        ASSERT_TEST(false);
    }
    catch(const Matrix<int>::AccessIllegalElement& e){
        ASSERT_TEST(string(e.what()) == "Mtm matrix error: An attempt to access an illegal element");
    }

    try{
        matrix(0,-2);
        ASSERT_TEST(false);
    }
    catch(const Matrix<int>::AccessIllegalElement& e){
        ASSERT_TEST(string(e.what()) == "Mtm matrix error: An attempt to access an illegal element");
    }

    try{
        matrix(-1,-1);
        ASSERT_TEST(false);
    }
    catch(const Matrix<int>::AccessIllegalElement& e){
        ASSERT_TEST(string(e.what()) == "Mtm matrix error: An attempt to access an illegal element");
    }

    const Matrix<int> matrix2(Dimensions(1,1));

    try{
        matrix2(1,0);
        ASSERT_TEST(false);
    }
    catch(const Matrix<int>::AccessIllegalElement& e){
        ASSERT_TEST(string(e.what()) == "Mtm matrix error: An attempt to access an illegal element");
    }

    try{
        matrix2(1,1);
        ASSERT_TEST(false);
    }
    catch(const Matrix<int>::AccessIllegalElement& e){
        ASSERT_TEST(string(e.what()) == "Mtm matrix error: An attempt to access an illegal element");
    }

    try{
        matrix2(0,1);
        ASSERT_TEST(false);
    }
    catch(const Matrix<int>::AccessIllegalElement& e){
        ASSERT_TEST(string(e.what()) == "Mtm matrix error: An attempt to access an illegal element");
    }

    return true;

}

bool testIterator(){

    int rows = 2;
    int cols = 5;
    Matrix<int> mat = Matrix<int>(Dimensions(rows,cols));
    ASSERT_TEST(mat.height() == rows);
    ASSERT_TEST(mat.width() == cols);
    ASSERT_TEST(mat.size() == rows*cols);
    int i = 0;
    for (int& element : mat){
        element = sampleData[i++];
    }
    ASSERT_TEST(i==mat.size());

    i = 0;
    int last = 0;
    for (int& element : mat){
        ASSERT_TEST((last = element) == sampleData[i++]);
        element = 0;
    }
    ASSERT_TEST(i==mat.size());
    ASSERT_TEST(last == 55);
    ASSERT_TEST(checkAreEqual(mat, Matrix<int>(Dimensions(rows, cols))));

    mat = Matrix<int>(Dimensions(rows,cols));
    ASSERT_TEST(mat.height() == rows);
    ASSERT_TEST(mat.width() == cols);
    ASSERT_TEST(mat.size() == rows*cols);
    i = 0;
    for (int& element : mat){
        element = sampleData[i++];
    }
    ASSERT_TEST(i==mat.size());

    i = 0;
    last = 0;
    const Matrix<int> mat2 = mat;
    for (const int& element : mat2){
        ASSERT_TEST((last = element) == sampleData[i++]);
    }
    ASSERT_TEST(i==mat2.size());
    ASSERT_TEST(last == 55);

    Matrix<int>::iterator it = mat.begin();
    it++;
    ASSERT_TEST(it == it);
    ASSERT_TEST(it != mat.begin());
    ASSERT_TEST(it != mat.end());
    ASSERT_TEST(it == ++mat.begin());
    ASSERT_TEST(it != mat.begin()++);
    ASSERT_TEST((mat.begin())++ != ++(mat.begin()));

    try{
        Matrix<int>::iterator it = mat.end();
        *it;
        ASSERT_TEST(false);
    }
    catch(const Matrix<int>::AccessIllegalElement& e){
        ASSERT_TEST(string(e.what()) == "Mtm matrix error: An attempt to access an illegal element");
    }

    try{
        Matrix<int>::iterator it = mat.end();
        it++;
    }
    catch(const Matrix<int>::AccessIllegalElement& e){
        ASSERT_TEST(false);
    }

    try{
        Matrix<int>::iterator it = mat.end();
        ++it;
    }
    catch(const Matrix<int>::AccessIllegalElement& e){
        ASSERT_TEST(false);
    }

    try{
        Matrix<int>::const_iterator it = mat2.end();
        *it;
        ASSERT_TEST(false);
    }
    catch(const Matrix<int>::AccessIllegalElement& e){
        ASSERT_TEST(string(e.what()) == "Mtm matrix error: An attempt to access an illegal element");
    }

    try{
        Matrix<int>::const_iterator it = mat2.end();
        it++;
    }
    catch(const Matrix<int>::AccessIllegalElement& e){
        ASSERT_TEST(false);
    }

    try{
        Matrix<int>::const_iterator it = mat2.end();
        ++it;
    }
    catch(const Matrix<int>::AccessIllegalElement& e){
        ASSERT_TEST(false);
    }

    return true;

}

bool testLogical(){

    int rows = 4;
    int cols = 5;
    Matrix<int> mat = Matrix<int>(Dimensions(rows,cols));
    ASSERT_TEST(mat.height() == rows);
    ASSERT_TEST(mat.width() == cols);
    ASSERT_TEST(mat.size() == rows*cols);
    int i = 0;
    for (int& element : mat){
        element = sampleData[i++];
    }
    ASSERT_TEST(i==mat.size());

    int scalar = 8;

    Matrix<bool> lessThan = mat < scalar;
    Matrix<bool> lessThanOrEqualTo = mat <= scalar;
    Matrix<bool> greaterThan = mat > scalar;
    Matrix<bool> greaterThanOrEqualTo = mat >= scalar;
    Matrix<bool> equalTo = mat == scalar;
    Matrix<bool> different = mat != scalar;

    for (int row = 0; row < rows; row++){
        for (int col = 0; col < cols; col++){
            ASSERT_TEST(lessThan(row,col)==(mat(row,col) < scalar));
            ASSERT_TEST(lessThanOrEqualTo(row,col)==(mat(row,col) <= scalar));
            ASSERT_TEST(greaterThan(row,col)==(mat(row,col) > scalar));
            ASSERT_TEST(greaterThanOrEqualTo(row,col)==(mat(row,col) >= scalar));
            ASSERT_TEST(equalTo(row,col)==(mat(row,col) == scalar));
            ASSERT_TEST(different(row,col)==(mat(row,col) != scalar));
        }
    }

    ASSERT_TEST(any(lessThan));
    ASSERT_TEST(any(lessThanOrEqualTo));
    ASSERT_TEST(any(greaterThan));
    ASSERT_TEST(any(greaterThanOrEqualTo));
    ASSERT_TEST(any(equalTo));
    ASSERT_TEST(any(different));

    return true;

}

class T2 : public string{
    public:
        T2(string name = ""):string(name){}
        operator bool() const{return *this == "def";};
};

bool testLogicalAnyAll(){

    int rows = 51;
    int cols = 17;
    Dimensions dim(rows, cols);

    ASSERT_TEST(!any(Matrix<int>(dim)));
    ASSERT_TEST(!all(Matrix<int>(dim)));

    ASSERT_TEST(any(Matrix<int>(dim, 1)));
    ASSERT_TEST(any(Matrix<int>(dim, 1020292)));
    ASSERT_TEST(any(Matrix<int>(dim, -1)));
    ASSERT_TEST(all(Matrix<int>(dim, 1)));
    ASSERT_TEST(all(Matrix<int>(dim, 1020292)));
    ASSERT_TEST(all(Matrix<int>(dim, -1)));

    int n = 10;
    ASSERT_TEST(any(Matrix<int>::Diagonal(n,1)));
    ASSERT_TEST(!all(Matrix<int>::Diagonal(n,1)));
    ASSERT_TEST(any(Matrix<int>::Diagonal(n,1) + 1));
    ASSERT_TEST(all(Matrix<int>::Diagonal(n,1) + 1));

    ASSERT_TEST(any(Matrix<T2>::Diagonal(n,T2("def"))));
    ASSERT_TEST(!all(Matrix<T2>::Diagonal(n,T2("def"))));
    ASSERT_TEST(any(Matrix<T2>::Diagonal(n,T2("def")) + T2("def")));
    ASSERT_TEST(!any(Matrix<T2>::Diagonal(n,T2("def")) + T2("*")));
    ASSERT_TEST(!all(Matrix<T2>::Diagonal(n,T2("def")) + T2("def")));

    return true;

}

class SignedSquareRoot{
    public:
        int operator()(int val){
            return val>=0 ? sqrt(val) : -sqrt(-val);
        }
};

double square(double val){
    return val*val;
}

bool testApply(){

    Matrix<int> mat(Dimensions(4,4), 5);
    mat(0,0) = 0;
    mat(1,1) = -9;

    Matrix<int> mat2 = mat.apply(SignedSquareRoot());
    ASSERT_TEST(mat2(0,0) == 0);
    ASSERT_TEST(mat2(1,1) == -3);

    mat2(0,0) = sqrt(5);
    mat2(1,1) = sqrt(5);
    ASSERT_TEST(checkAreEqual(mat2, Matrix<int>(Dimensions(4,4), sqrt(5))));

    Matrix<double> mat3(Dimensions(4,4),9);
    mat3 = mat3.apply(square);
    ASSERT_TEST(checkAreEqual(mat3, Matrix<double>(Dimensions(4,4), 81)));

    return true;

}

bool run_test(std::function<bool()> test, std::string test_name){
    if(!test()){
        cout<<test_name<<" - FAILED."<<endl;
        return false;
    }
    cout<<test_name<<" - PASSED."<<endl;
    return true;
}
int main(){

    createSampleData();

    std::map<std::string, std::function<bool()>> tests;

    ADD_TEST(testConstructor);
    ADD_TEST(testCopyConstructor);
    ADD_TEST(testAssignmentOperator);
    ADD_TEST(testDiagonal);
    ADD_TEST(testTranspose);
    ADD_TEST(testOperatorNegative);
    ADD_TEST(testOperatorMatrixAddition);
    ADD_TEST(testOperatorScalarAddition);
    ADD_TEST(testOperatorOutput);
    ADD_TEST(testIterator);
    ADD_TEST(testLogical);
    ADD_TEST(testLogicalAnyAll);
    ADD_TEST(testApply);
    ADD_TEST(testOperatorParenthesis);

    int passed = 0;
    for (std::pair<std::string, std::function<bool()>> element : tests)
    {
        passed += run_test(element.second,"Running "+element.first);
    }

    cout<<"Passed "<<passed<<" out of "<<tests.size()<<" tests."<<endl;

    #ifndef NDEPENDENCY

        std::fstream fout; 
        fout.open("dependencies.csv", std::ios::out); 
        findDependencies(fout);
        fout.close();

    #endif

    return 0;
}
