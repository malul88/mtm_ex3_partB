#include <functional>
#include <string>
#include <iostream>
#include <sstream>   
#include <map>
#include <iterator>
#include <string>

#include "IntMatrix.h"
#include "../Auxiliaries.h"

using namespace mtm;
using std::cout;
using std::endl;
using std::string;
using std::to_string;

#define ASSERT_TEST(x) if(!(x)){ \
 cout<<("Failed assertion at line " + to_string(__LINE__) + " in " + __func__)<<endl;\
  return false; }

#define ADD_TEST(x) tests[#x]=x;

bool checkAreEqual(const IntMatrix& first, const IntMatrix& second){

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

bool testConstructor(){
    
    int rows = 57;
    int cols = 27;
    int scalar = 7;
    IntMatrix matrix(Dimensions(rows, cols), scalar);

    ASSERT_TEST(matrix.height() == rows);
    ASSERT_TEST(matrix.width() == cols);
    ASSERT_TEST(matrix.size() == (rows * cols));

    for (int& element : matrix){
        ASSERT_TEST(element == scalar);
    }

    IntMatrix matrix2(Dimensions(rows, cols));
    for (int& element : matrix2){
        ASSERT_TEST(element == 0);
    }

    ASSERT_TEST(!checkAreEqual(matrix, matrix2));

    return true;

}

bool testCopyConstructor(){

    int rows = 57;
    int cols = 27;
    int scalar = -5;
    IntMatrix matrix(Dimensions(rows, cols), scalar);

    ASSERT_TEST(matrix.height() == rows);
    ASSERT_TEST(matrix.width() == cols);
    ASSERT_TEST(matrix.size() == (rows * cols));

    for (int& element : matrix){
        ASSERT_TEST(element == scalar);
    }

    IntMatrix matrix2(matrix);
    ASSERT_TEST(matrix2.height() == rows);
    ASSERT_TEST(matrix2.width() == cols);
    ASSERT_TEST(matrix2.size() == (rows * cols));
    for (int& element : matrix2){
        ASSERT_TEST(element == scalar);
    }

    ASSERT_TEST(checkAreEqual(matrix, matrix2));

    IntMatrix matrix3 = matrix;
    ASSERT_TEST(matrix3.height() == rows);
    ASSERT_TEST(matrix3.width() == cols);
    ASSERT_TEST(matrix3.size() == (rows * cols));
    for (int& element : matrix3){
        ASSERT_TEST(element == scalar);
    }

    ASSERT_TEST(checkAreEqual(matrix2, matrix3));
    ASSERT_TEST(checkAreEqual(matrix, matrix3));

    return true;

}

bool testAssignmentOperator(){

    int rows = 57;
    int cols = 27;
    int scalar = -5;
    IntMatrix matrix(Dimensions(rows, cols), scalar);

    ASSERT_TEST(matrix.height() == rows);
    ASSERT_TEST(matrix.width() == cols);
    ASSERT_TEST(matrix.size() == (rows * cols));
    for (int& element : matrix){
        ASSERT_TEST(element == scalar);
    }

    IntMatrix matrix2(Dimensions(1,2), -1);
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

    return true;

}

bool testIdentity(){

    int n=107;
    IntMatrix mat = IntMatrix::Identity(n);

    ASSERT_TEST(mat.height() == n);
    ASSERT_TEST(mat.width() == n);
    ASSERT_TEST(mat.size() == n*n);

    for (int row = 0; row < n; row++){
        for (int col = 0; col < n; col++){
            ASSERT_TEST(mat(row,col) == (row==col));
        }
    }

    return true;

}

bool testTranspose(){

    int n=107;
    IntMatrix mat = mtm::IntMatrix::Identity(n);
    ASSERT_TEST(checkAreEqual(mat, mat.transpose()));

    int rows = 17;
    int cols = 5;
    IntMatrix mat2 = IntMatrix(Dimensions(rows,cols));
    ASSERT_TEST(mat2.height() == rows);
    ASSERT_TEST(mat2.width() == cols);
    ASSERT_TEST(mat2.size() == rows*cols);
    int i = 0;
    for (int& element : mat2){
        element = sampleData[i++];
    }
    
    IntMatrix mat3 = mat2.transpose();
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
    IntMatrix mat = IntMatrix(Dimensions(rows,cols));
    ASSERT_TEST(mat.height() == rows);
    ASSERT_TEST(mat.width() == cols);
    ASSERT_TEST(mat.size() == rows*cols);
    int i = 0;
    for (int& element : mat){
        element = sampleData[i++];
    }
    
    IntMatrix mat2 = -mat;
    ASSERT_TEST(mat2.height() == rows);
    ASSERT_TEST(mat2.width() == cols);
    ASSERT_TEST(mat2.size() == rows*cols);

    for (int row = 0; row < rows; row++){
        for (int col = 0; col < cols; col++){
            ASSERT_TEST(mat2(row,col)==-mat(row,col));
        }
    }

    ASSERT_TEST(checkAreEqual(mat, -(-mat)));

    return true;

}

bool testOperatorMatrixAddition(){

    int rows = 17;
    int cols = 5;
    IntMatrix mat = IntMatrix(Dimensions(rows,cols));
    ASSERT_TEST(mat.height() == rows);
    ASSERT_TEST(mat.width() == cols);
    ASSERT_TEST(mat.size() == rows*cols);
    int i = 0;
    for (int& element : mat){
        element = sampleData[i++];
    }
    
    IntMatrix mat2 = -mat;
    ASSERT_TEST(mat2.height() == rows);
    ASSERT_TEST(mat2.width() == cols);
    ASSERT_TEST(mat2.size() == rows*cols);

    ASSERT_TEST(checkAreEqual((mat=mat+mat2), IntMatrix(Dimensions(rows,cols))));
    ASSERT_TEST(!checkAreEqual(mat+mat2, IntMatrix(Dimensions(rows,cols))));
    ASSERT_TEST(checkAreEqual((-mat2) + mat2, IntMatrix(Dimensions(rows,cols))));

    IntMatrix mat3 = mat2 + mat2;
    mat3 = mat3 + mat2;
    mat3 = mat3 - mat2;
    i = 0;
    for (int& element : mat3){
        ASSERT_TEST(element == -2*sampleData[i++]);
    }

    mat3=((mat3=(mat3 - mat2)) - mat2);

    ASSERT_TEST(checkAreEqual(mat3, IntMatrix(Dimensions(rows,cols))));

    return true;

}

bool testOperatorScalarAddition(){

    int rows = 17;
    int cols = 5;
    IntMatrix mat = IntMatrix(Dimensions(rows,cols));
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

    IntMatrix mat2 = IntMatrix(Dimensions(rows,cols), 5);
    mat2 += -5;
    ASSERT_TEST(checkAreEqual(1 + mat2 + (-1), IntMatrix(Dimensions(rows, cols))));

    return true;

}

bool testOperatorOutput(){

    int rows = 2;
    int cols = 5;
    IntMatrix mat = IntMatrix(Dimensions(rows,cols));
    ASSERT_TEST(mat.height() == rows);
    ASSERT_TEST(mat.width() == cols);
    ASSERT_TEST(mat.size() == rows*cols);
    int i = 0;
    for (int& element : mat){
        element = sampleData[i++];
    }

    std::string expected = "1 1 2 3 5 \n8 13 21 34 55 \n\n";
    std::stringstream buffer;
    buffer << mat;
    ASSERT_TEST(expected == buffer.str());
    return true;

}

bool testIterator(){

    int rows = 2;
    int cols = 5;
    IntMatrix mat = IntMatrix(Dimensions(rows,cols));
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
    ASSERT_TEST(checkAreEqual(mat, IntMatrix(Dimensions(rows, cols))));

    mat = IntMatrix(Dimensions(rows,cols));
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
    const IntMatrix mat2 = mat;
    for (const int& element : mat2){
        ASSERT_TEST((last = element) == sampleData[i++]);
    }
    ASSERT_TEST(i==mat2.size());
    ASSERT_TEST(last == 55);

    IntMatrix::iterator it = mat.begin();
    it++;
    ASSERT_TEST(it == it);
    ASSERT_TEST(it != mat.begin());
    ASSERT_TEST(it != mat.end());
    ASSERT_TEST(it == ++mat.begin());
    ASSERT_TEST(it != mat.begin()++);
    ASSERT_TEST((mat.begin())++ != ++(mat.begin()));

    return true;

}

bool testLogical(){

    int rows = 4;
    int cols = 5;
    IntMatrix mat = IntMatrix(Dimensions(rows,cols));
    ASSERT_TEST(mat.height() == rows);
    ASSERT_TEST(mat.width() == cols);
    ASSERT_TEST(mat.size() == rows*cols);
    int i = 0;
    for (int& element : mat){
        element = sampleData[i++];
    }
    ASSERT_TEST(i==mat.size());

    int scalar = 8;

    IntMatrix lessThan = mat < scalar;
    IntMatrix lessThanOrEqualTo = mat <= scalar;
    IntMatrix greaterThan = mat > scalar;
    IntMatrix greaterThanOrEqualTo = mat >= scalar;
    IntMatrix equalTo = mat == scalar;
    IntMatrix different = mat != scalar;

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

bool testLogicalAnyAll(){

    int rows = 51;
    int cols = 17;
    Dimensions dim(rows, cols);

    ASSERT_TEST(!any(IntMatrix(dim)));
    ASSERT_TEST(!all(IntMatrix(dim)));

    ASSERT_TEST(any(IntMatrix(dim, 1)));
    ASSERT_TEST(any(IntMatrix(dim, 1020292)));
    ASSERT_TEST(any(IntMatrix(dim, -1)));
    ASSERT_TEST(all(IntMatrix(dim, 1)));
    ASSERT_TEST(all(IntMatrix(dim, 1020292)));
    ASSERT_TEST(all(IntMatrix(dim, -1)));

    int n = 10;
    ASSERT_TEST(any(IntMatrix::Identity(n)));
    ASSERT_TEST(!all(IntMatrix::Identity(n)));
    ASSERT_TEST(any(IntMatrix::Identity(n) + 1));
    ASSERT_TEST(all(IntMatrix::Identity(n) + 1));

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
    ADD_TEST(testIdentity);
    ADD_TEST(testTranspose);
    ADD_TEST(testOperatorNegative);
    ADD_TEST(testOperatorMatrixAddition);
    ADD_TEST(testOperatorScalarAddition);
    ADD_TEST(testOperatorOutput);
    ADD_TEST(testIterator);
    ADD_TEST(testLogical);
    ADD_TEST(testLogicalAnyAll);

    int passed = 0;
    for (std::pair<std::string, std::function<bool()>> element : tests)
    {
        passed += run_test(element.second,"Running "+element.first);
    }

    cout<<"Passed "<<passed<<" out of "<<tests.size()<<" tests."<<endl;

    return 0;
}